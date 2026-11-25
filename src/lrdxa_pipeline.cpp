#include <volt/lrdxa_pipeline.h>

#include <volt/lrdxa_engine.h>
#include <volt/lrdxa_log.h>
#include <volt/lrdxa_outputs.h>
#include <volt/core/progress.h>

namespace Volt::DXA {

json LineReconstructionDXAPipeline::run(
    const LammpsParser::Frame& frame,
    const std::string& outputFile,
    const LineReconstructionDXAOptions& options,
    StructureAnalysis& reconstructionStructureAnalysis,
    StructureContext& reconstructionContext,
    LineReconstructionJsonExporter& jsonExporter,
    LineReconstructionDXALog& log
) {
    LineReconstructionDXAAlgorithm algorithm(reconstructionStructureAnalysis, reconstructionContext);
    // Weighed by a sample run.
    ProgressSteps steps({
        {"Reconstructing dislocation lines", 16.0},
        {"Writing outputs", LineReconstructionOutputs::weight(), !outputFile.empty()},
    });
    algorithm.run(options, log);
    log.dump(algorithm);

    steps.next();
    if(!outputFile.empty()) {
        LineReconstructionOutputs{
            outputFile, frame, reconstructionStructureAnalysis, algorithm, jsonExporter, log
        }.write();
    }

    json result;
    result["is_failed"] = false;
    result["algorithm"] = "line-reconstruction-dxa";
    result["segments"] = algorithm.dislocationSegments().size();
    result["lines"] = algorithm.dislocationLineCount();
    result["unassigned_edges"] = algorithm.unassignedEdges().size();
    result["delaunay_vertices"] = algorithm.delaunayVertexCount();
    return result;
}

}
