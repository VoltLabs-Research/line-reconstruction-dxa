#include <volt/lrdxa_pipeline.h>

#include <volt/lrdxa_engine.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/utilities/json_utils.h>
#include <volt/lrdxa_json_exporter.h>

namespace Volt::DXA {

json LineReconstructionDXAPipeline::run(
    const LammpsParser::Frame& frame,
    const std::string& outputFile,
    const LineReconstructionDXAOptions& options,
    StructureAnalysis& reconstructionStructureAnalysis,
    StructureContext& reconstructionContext,
    LineReconstructionJsonExporter& jsonExporter,
    const std::function<void(std::string_view)>& markStage
) {
    LineReconstructionDXAAlgorithm algorithm(reconstructionStructureAnalysis, reconstructionContext);
    algorithm.run(options, markStage);

    if(!outputFile.empty()) {
        jsonExporter.writeDislocationLinesParquetToFile(
            algorithm.dislocationLines(),
            frame.simulationCell,
            outputFile + "_dislocations.parquet"
        );
        if(markStage) markStage("stream_dislocations_parquet");
        jsonExporter.writeLineSegmentsParquetToFile(
            algorithm.dislocationSegments(),
            frame.simulationCell,
            outputFile + "_dislocation_segments.parquet"
        );
        if(markStage) markStage("stream_dislocation_segments_parquet");
        jsonExporter.writeUnassignedEdgesParquetToFile(
            algorithm.unassignedEdges(),
            outputFile + "_unassigned_edges.parquet"
        );
        if(markStage) markStage("stream_unassigned_edges_parquet");
        jsonExporter.writeInterfaceMeshParquetToFile(
            algorithm.interfaceMesh(),
            reconstructionStructureAnalysis,
            outputFile + "_interface_mesh.parquet"
        );
        if(markStage) markStage("stream_interface_mesh_parquet");
        const auto simCellInfo = jsonExporter.getExtendedSimulationCellInfo(frame.simulationCell);
        JsonUtils::writeJsonToParquet(simCellInfo, outputFile + "_simulation_cell.parquet", false);
        if(markStage) markStage("stream_simulation_cell_parquet");
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
