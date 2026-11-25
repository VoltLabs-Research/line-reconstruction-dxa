#include <volt/lrdxa_service.h>
#include <volt/lrdxa_log.h>
#include <volt/lrdxa_pipeline.h>
#include <volt/core/progress.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/analysis/structure_analysis.h>
#include <utility>

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
    DXA::LineReconstructionDXALog log;
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

    ReconstructedStructureContext context(prepared.positions.get(), frame.simulationCell);
    context.inputCrystalType = _inputCrystalStructure;

    StructureAnalysis structureAnalysis(context);
    ProgressSteps steps({
        {"Loading the reconstructed structure", 1.0},
        {"Reconstructing dislocation lines", 20.0},
    });
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
    log.dump(structureAnalysis);

    steps.next();
    return DXA::LineReconstructionDXAPipeline::run(
        frame,
        outputFile,
        buildOptions(),
        structureAnalysis,
        context,
        _jsonExporter,
        log
    );
}

}
