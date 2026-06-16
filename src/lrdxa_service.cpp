#include <volt/lrdxa_service.h>
#include <volt/lrdxa_pipeline.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/analysis/structure_analysis.h>
#include <spdlog/spdlog.h>
#include <chrono>
#include <utility>
#include <string_view>

namespace Volt {

using namespace Volt::Particles;

LineReconstructionDXA::LineReconstructionDXA()
    : _inputCrystalStructure(LATTICE_FCC)
    , _crystalPathSteps(4)
    , _tessellationGhostLayerScale(3.5)
    , _alphaScale(3.5)
    , _smoothingIterations(0)
    , _linePointInterval(1.2) {
}

void LineReconstructionDXA::setInputCrystalStructure(LatticeStructureType structure) {
    _inputCrystalStructure = structure;
}

void LineReconstructionDXA::setClustersTablePath(std::string path) {
    _clustersTablePath = std::move(path);
}

void LineReconstructionDXA::setClusterTransitionsPath(std::string path) {
    _clusterTransitionsPath = std::move(path);
}

void LineReconstructionDXA::setNeighborLatticePath(std::string path) {
    _neighborLatticePath = std::move(path);
}

void LineReconstructionDXA::setCrystalPathSteps(int crystalPathSteps) {
    _crystalPathSteps = crystalPathSteps;
}

void LineReconstructionDXA::setTessellationGhostLayerScale(double tessellationGhostLayerScale) {
    _tessellationGhostLayerScale = tessellationGhostLayerScale;
}

void LineReconstructionDXA::setAlphaScale(double alphaScale) {
    _alphaScale = alphaScale;
}

void LineReconstructionDXA::setSmoothingIterations(int smoothingIterations) {
    _smoothingIterations = smoothingIterations;
}

void LineReconstructionDXA::setLinePointInterval(double linePointInterval) {
    _linePointInterval = linePointInterval;
}

LineReconstructionDXAOptions LineReconstructionDXA::buildOptions() const {
    return LineReconstructionDXAOptions{
        .crystalPathSteps = _crystalPathSteps,
        .tessellationGhostLayerScale = _tessellationGhostLayerScale,
        .alphaScale = _alphaScale,
        .smoothingIterations = _smoothingIterations,
        .linePointInterval = _linePointInterval,
    };
}

json LineReconstructionDXA::compute(const LammpsParser::Frame& frame, const std::string& outputFile) {
    const auto startTime = std::chrono::high_resolution_clock::now();
    auto stageStart = startTime;
    const LineReconstructionDXAOptions options = buildOptions();

    json result;
    json stageMetrics = json::array();
    auto markStage = [&](std::string_view name) {
        const auto now = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - stageStart).count();
        stageMetrics.push_back({{"stage", name}, {"elapsed_ms", elapsed}});
        stageStart = now;
    };

    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError)) {
        return AnalysisResult::failure(frameError);
    }
    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()) {
        return AnalysisResult::failure(
            "LineReconstructionDXA requires --clusters-table and --clusters-transitions"
        );
    }
    if(_neighborLatticePath.empty()) {
        return AnalysisResult::failure(
            "LineReconstructionDXA requires --neighbor_lattice (per-atom neighbor topology parquet)"
        );
    }

    auto positions = std::move(prepared.positions);
    markStage("create_position_property");

    ReconstructedStructureContext context(positions.get(), frame.simulationCell);
    context.inputCrystalType = _inputCrystalStructure;

    StructureAnalysis structureAnalysis(context);
    std::string reconstructionError;
    if(!ReconstructedStructureLoader::load(
        frame,
        _neighborLatticePath,
        {_clustersTablePath, _clusterTransitionsPath},
        structureAnalysis,
        context,
        &reconstructionError
    )) {
        return AnalysisResult::failure(reconstructionError);
    }
    markStage("load_reconstructed_structure");

    result = DXA::LineReconstructionDXAPipeline::run(
        frame,
        outputFile,
        options,
        structureAnalysis,
        context,
        _jsonExporter,
        markStage
    );

    result["total_time"] = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - startTime
    ).count();
    result["stage_metrics"] = std::move(stageMetrics);
    return result;
}

}
