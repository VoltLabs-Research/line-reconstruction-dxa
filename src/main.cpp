#include <volt/cli/common.h>
#include <volt/core/progress.h>
#include <volt/lrdxa_service.h>
#include <volt/plugin/plugin_main.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/structures/crystal_topology_registry.h>
#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>
#include <algorithm>
#include <fstream>
#include <set>
#include <string>

using namespace Volt;
using namespace Volt::CLI;
using namespace Volt::Plugin;

namespace {

LatticeStructureType parseCrystalStructure(const std::string& str) {
    if(str == "FCC") return LATTICE_FCC;
    if(str == "BCC") return LATTICE_BCC;
    if(str == "HCP") return LATTICE_HCP;
    if(str == "SC") return LATTICE_SC;
    if(str == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if(str == "HEX_DIAMOND") return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", str);
    return LATTICE_FCC;
}

PluginDescriptor buildDescriptor() {
    return {
        "volt-line-reconstruction-dxa",
        "Line Reconstruction DXA",
        {
            {"--clusters-table", "path", "Path to *_clusters.table exported upstream.", "", {}, ""},
            {"--clusters-transitions", "path", "Path to *_cluster_transitions.table exported upstream.", "", {}, ""},
            {"--neighbor_lattice", "path", "Path to *_neighbor_lattice.parquet exported upstream.", "", {}, ""},
            {"--crystalStructure", "enum", "Reference crystal structure.", "FCC",
             {"BCC", "FCC", "HCP", "CUBIC_DIAMOND", "HEX_DIAMOND", "SC"}, ""},
            {"--lattice-dir", "path", "Directory containing lattice topology YAMLs.", "", {}, ""},
            {"--crystalPathSteps", "int", "Maximum crystal-path steps used for edge vectors.", "4", {}, ""},
            {"--tessellationGhostLayerScale", "float", "Ghost-layer scale relative to neighbor distance.", "3.5", {}, ""},
            {"--alphaScale", "float", "Alpha threshold scale relative to neighbor distance.", "3.5", {}, ""},
            {"--smoothingIterations", "int", "Taubin smoothing iterations for reconstructed lines.", "0", {}, ""},
            {"--linePointInterval", "float", "Line coarsening interval.", "1.2", {}, ""},
        }
    };
}

}

int main(int argc, char* argv[]) {
    const PluginDescriptor descriptor = buildDescriptor();

    if(argc < 2) {
        showPluginUsage(argv[0], descriptor);
        return 1;
    }

    std::string filename;
    std::string outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(auto exitCode = handleIntrospection(argv[0], descriptor, opts, filename)) {
        return *exitCode;
    }

    if(!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if(cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while(std::getline(cpuinfo, line)) {
                if(line.empty()) {
                    if(physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                int parsedValue = 0;
                const auto separator = line.find(':');
                if(separator != std::string::npos) {
                    try {
                        parsedValue = std::stoi(line.substr(separator + 1));
                    } catch(const std::exception&) {
                        parsedValue = 0;
                    }
                }
                if(line.rfind("physical id", 0) == 0) {
                    physicalId = parsedValue;
                } else if(line.rfind("core id", 0) == 0) {
                    coreId = parsedValue;
                } else if(line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, parsedValue);
                }
            }
            if(physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = physicalCoreIds.empty() ? fallbackCpuCores : static_cast<int>(physicalCoreIds.size());
        }
        opts["--threads"] = std::to_string(std::max(1, std::min(maxAvailableThreads, physicalCores > 0 ? physicalCores : maxAvailableThreads)));
    }

    const int requestedThreads = getInt(opts, "--threads");
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(std::max(1, requestedThreads))
    );
    initLogging("line-reconstruction-dxa");
    const Stopwatch runTimer;
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    const std::string crystalStructure = getString(opts, "--crystalStructure", "FCC");
    const int crystalPathSteps = getInt(opts, "--crystalPathSteps", 4);
    const double tessellationGhostLayerScale = getDouble(opts, "--tessellationGhostLayerScale", 3.5);
    const double alphaScale = getDouble(opts, "--alphaScale", 3.5);
    const int smoothingIterations = getInt(opts, "--smoothingIterations", 0);
    const double linePointInterval = getDouble(opts, "--linePointInterval", 1.2);
    logParameters(descriptor, opts);

    const std::string latticeDirectory = getString(opts, "--lattice-dir", "");
    if(!latticeDirectory.empty()){
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    LineReconstructionDXAService analyzer;
    analyzer.setClustersTablePath(getString(opts, "--clusters-table"));
    analyzer.setClusterTransitionsPath(getString(opts, "--clusters-transitions"));
    analyzer.setNeighborLatticePath(getString(opts, "--neighbor_lattice"));
    analyzer.setInputCrystalStructure(parseCrystalStructure(crystalStructure));
    analyzer.setCrystalPathSteps(crystalPathSteps);
    analyzer.setTessellationGhostLayerScale(tessellationGhostLayerScale);
    analyzer.setAlphaScale(alphaScale);
    analyzer.setSmoothingIterations(smoothingIterations);
    analyzer.setLinePointInterval(linePointInterval);

    spdlog::info("Starting line reconstruction DXA analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)) {
        spdlog::error("Analysis failed after {}: {}", runTimer.elapsed(), result.value("error", "Unknown error"));
        return 1;
    }
    spdlog::info("Line Reconstruction DXA completed in {}.", runTimer.elapsed());
    return 0;
}
