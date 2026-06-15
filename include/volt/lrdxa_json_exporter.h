#pragma once

#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <volt/analysis/structure_analysis.h>
#include <volt/lrdxa_types.h>
#include <volt/math/lin_alg.h>

namespace Volt {

using json = nlohmann::json;

class LineReconstructionJsonExporter {
public:
    explicit LineReconstructionJsonExporter() = default;

    json getExtendedSimulationCellInfo(const SimulationCell& cell);

    void writeLineSegmentsParquetToFile(
        const std::vector<DXA::LineReconstructionSegment>& segments,
        const SimulationCell& simulationCell,
        const std::string& filePath
    );

    void writeDislocationLinesParquetToFile(
        const std::vector<DXA::LineReconstructionDislocationLine>& lines,
        const SimulationCell& simulationCell,
        const std::string& filePath
    );

    void writeUnassignedEdgesParquetToFile(
        const std::vector<DXA::LineReconstructionUnassignedEdge>& edges,
        const std::string& filePath
    );

    void writeInterfaceMeshParquetToFile(
        const DXA::LineReconstructionInterfaceMesh& mesh,
        const StructureAnalysis& structureAnalysis,
        const std::string& filePath
    );

private:
    json exportUnassignedEdgesToJson(const std::vector<DXA::LineReconstructionUnassignedEdge>& edges);

    template <typename MeshType>
    json getMeshData(const MeshType& mesh, const StructureAnalysis& structureAnalysis);

    json vectorToJson(const Vector3& vector);
    json affineTransformationToJson(const AffineTransformation& transform);
    json simulationCellToJson(const SimulationCell& cell);
    double calculateAngle(const Vector3& a, const Vector3& b);
};

}
