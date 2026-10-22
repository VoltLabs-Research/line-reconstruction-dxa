#include <volt/lrdxa_pipeline.h>

#include <volt/lrdxa_engine.h>
#include <volt/core/progress.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/utilities/json_utils.h>
#include <volt/lrdxa_json_exporter.h>
#include <spdlog/spdlog.h>

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
        Stopwatch exportTimer;
        const std::string dislocationsPath = outputFile + "_dislocations.parquet";
        jsonExporter.writeDislocationLinesParquetToFile(
            algorithm.dislocationLines(),
            frame.simulationCell,
            dislocationsPath
        );
        spdlog::info("Exported {} dislocation lines to: {} ({})",
            algorithm.dislocationLineCount(), dislocationsPath, exportTimer.elapsed());
        if(markStage) markStage("stream_dislocations_parquet");
        exportTimer.restart();
        const std::string segmentsPath = outputFile + "_dislocation_segments.parquet";
        jsonExporter.writeLineSegmentsParquetToFile(
            algorithm.dislocationSegments(),
            frame.simulationCell,
            segmentsPath
        );
        spdlog::info("Exported {} snapshot segments to: {} ({})",
            algorithm.dislocationSegments().size(), segmentsPath, exportTimer.elapsed());
        if(markStage) markStage("stream_dislocation_segments_parquet");
        exportTimer.restart();
        const std::string unassignedEdgesPath = outputFile + "_unassigned_edges.parquet";
        jsonExporter.writeUnassignedEdgesParquetToFile(
            algorithm.unassignedEdges(),
            unassignedEdgesPath
        );
        spdlog::info("Exported {} unresolved edge records to: {} ({})",
            algorithm.unassignedEdges().size(), unassignedEdgesPath, exportTimer.elapsed());
        if(markStage) markStage("stream_unassigned_edges_parquet");
        spdlog::info("Exporting interface mesh ({} facets, {} vertices)...",
            algorithm.interfaceMesh().faceCount(), algorithm.interfaceMesh().vertexCount());
        exportTimer.restart();
        const std::string interfaceMeshPath = outputFile + "_interface_mesh.parquet";
        jsonExporter.writeInterfaceMeshParquetToFile(
            algorithm.interfaceMesh(),
            reconstructionStructureAnalysis,
            interfaceMeshPath
        );
        spdlog::info("Exported interface mesh to: {} ({})", interfaceMeshPath, exportTimer.elapsed());
        if(markStage) markStage("stream_interface_mesh_parquet");
        exportTimer.restart();
        const std::string simulationCellPath = outputFile + "_simulation_cell.parquet";
        const auto simCellInfo = jsonExporter.getExtendedSimulationCellInfo(frame.simulationCell);
        JsonUtils::writeJsonToParquet(simCellInfo, simulationCellPath, false);
        spdlog::info("Exported simulation cell to: {} ({})", simulationCellPath, exportTimer.elapsed());
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
