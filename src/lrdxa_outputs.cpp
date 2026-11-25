#include <volt/lrdxa_outputs.h>

#include <volt/lrdxa_engine.h>
#include <volt/lrdxa_json_exporter.h>
#include <volt/lrdxa_log.h>
#include <volt/utilities/json_utils.h>

namespace Volt::DXA {

const std::array<OutputFile<LineReconstructionExports, LineReconstructionOutputs>, 5> LineReconstructionOutputs::files{{
    {"Writing dislocations", 0.1, nullptr, [](const LineReconstructionOutputs& out) {
        const std::string path = out.path("dislocations");
        out.exporter.writeDislocationLinesParquetToFile(
            out.reconstruction.dislocationLines(), out.frame.simulationCell, path);
        out.log.dump(ExportedTable{
            ExportedTable::Contents::DislocationLines, path, out.reconstruction.dislocationLineCount()});
    }},
    {"Writing dislocation segments", 0.8, nullptr, [](const LineReconstructionOutputs& out) {
        const std::string path = out.path("dislocation_segments");
        out.exporter.writeLineSegmentsParquetToFile(
            out.reconstruction.dislocationSegments(), out.frame.simulationCell, path);
        out.log.dump(ExportedTable{
            ExportedTable::Contents::SnapshotSegments, path, out.reconstruction.dislocationSegments().size()});
    }},
    {"Writing unassigned edges", 0.1, nullptr, [](const LineReconstructionOutputs& out) {
        const std::string path = out.path("unassigned_edges");
        out.exporter.writeUnassignedEdgesParquetToFile(out.reconstruction.unassignedEdges(), path);
        out.log.dump(ExportedTable{
            ExportedTable::Contents::UnresolvedEdges, path, out.reconstruction.unassignedEdges().size()});
    }},
    {"Writing the interface mesh", 6.0, nullptr, [](const LineReconstructionOutputs& out) {
        const std::string path = out.path("interface_mesh");
        out.exporter.writeInterfaceMeshParquetToFile(out.reconstruction.interfaceMesh(), out.structure, path);
        out.log.dump(ExportedTable{ExportedTable::Contents::InterfaceMesh, path});
    }},
    {"Writing the simulation cell", 0.1, nullptr, [](const LineReconstructionOutputs& out) {
        const std::string path = out.path("simulation_cell");
        JsonUtils::writeJsonToParquet(out.exporter.getExtendedSimulationCellInfo(out.frame.simulationCell), path, false);
        out.log.dump(ExportedTable{ExportedTable::Contents::SimulationCell, path});
    }},
}};

double LineReconstructionOutputs::weight() {
    return totalWeight(outputSteps(files, LineReconstructionExports{}));
}

void LineReconstructionOutputs::write() const {
    writeOutputFiles(files, LineReconstructionExports{}, *this);
}

}
