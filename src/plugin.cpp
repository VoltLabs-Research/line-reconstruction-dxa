#include <volt/lrdxa_plugin.h>
#include <volt/lrdxa_service.h>
#include <volt/plugin/plugin_entry.h>
#include <volt/structures/crystal_topology_registry.h>

#include <string>
#include <utility>
#include <vector>

using namespace Volt;
using namespace Volt::Plugin;
using S = LineReconstructionDXA;

static const std::vector<OptionBinding<S>> bindings = {
    opt<S>(CliOption{"--clusters-table", "path", "Path to *_clusters.table exported upstream.", ""},
        [](S& s, const OptsMap& opts) { s.setClustersTablePath(CLI::getString(opts, "--clusters-table")); }),
    opt<S>(CliOption{"--clusters-transitions", "path", "Path to *_cluster_transitions.table exported upstream.", ""},
        [](S& s, const OptsMap& opts) { s.setClusterTransitionsPath(CLI::getString(opts, "--clusters-transitions")); }),
    opt<S>(CliOption{"--neighbor_lattice", "path", "Path to *_neighbor_lattice.parquet exported upstream.", ""},
        [](S& s, const OptsMap& opts) { s.setNeighborLatticePath(CLI::getString(opts, "--neighbor_lattice")); }),
    opt<S>(CliOption{"--crystalStructure", "enum", "Reference crystal structure.", "FCC",
            {"BCC", "FCC", "HCP", "CUBIC_DIAMOND", "HEX_DIAMOND", "SC"}},
        [](S& s, const OptsMap& opts) {
            s.setInputCrystalStructure(
                LrdxaPlugin::parseCrystalStructure(CLI::getString(opts, "--crystalStructure", "FCC")));
        }),
    opt("--crystalPathSteps", "Maximum crystal-path steps used for edge vectors.", 4, &S::setCrystalPathSteps),
    opt("--tessellationGhostLayerScale", "Ghost-layer scale relative to neighbor distance.", 3.5,
        &S::setTessellationGhostLayerScale),
    opt("--alphaScale", "Alpha threshold scale relative to neighbor distance.", 3.5, &S::setAlphaScale),
    opt("--smoothingIterations", "Taubin smoothing iterations for reconstructed lines.", 0, &S::setSmoothingIterations),
    opt("--linePointInterval", "Line coarsening interval.", 1.2, &S::setLinePointInterval),
};

namespace Volt::LrdxaPlugin {

LatticeStructureType parseCrystalStructure(const std::string& name) {
    if(name == "FCC") return LATTICE_FCC;
    if(name == "BCC") return LATTICE_BCC;
    if(name == "HCP") return LATTICE_HCP;
    if(name == "SC") return LATTICE_SC;
    if(name == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if(name == "HEX_DIAMOND") return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", name);
    return LATTICE_FCC;
}

PluginDescriptor buildDescriptor() {
    auto meta = optionsMeta(bindings);
    meta.insert(meta.begin() + 4, {"--lattice-dir", "path", "Directory containing lattice topology YAMLs.", "", {}, ""});
    PluginDescriptor descriptor{"volt-line-reconstruction-dxa", "Line Reconstruction DXA", std::move(meta)};
    descriptor.defaultThreads = DefaultThreads::PhysicalCores;
    return descriptor;
}

}

VOLT_PLUGIN_MAIN(LrdxaPlugin::buildDescriptor(),
    [](const OptsMap& opts, const LammpsParser::Frame& frame,
       const LammpsParser::Frame*, const std::string& outputBase) -> json {
    const std::string latticeDirectory = CLI::getString(opts, "--lattice-dir", "");
    if(!latticeDirectory.empty()) {
        setCrystalTopologySearchRoot(latticeDirectory);
        spdlog::info("Using lattice directory: {}", latticeDirectory);
    }

    S analyzer;
    applyAll(analyzer, bindings, opts);
    return analyzer.compute(frame, outputBase);
})
