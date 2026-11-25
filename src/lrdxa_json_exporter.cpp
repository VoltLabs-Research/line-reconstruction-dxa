#include <volt/lrdxa_json_exporter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <volt/core/progress.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/dxa_serialization_detail.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_line_writer.h>

namespace Volt {

using DxaSerialization::Detail::getGlobalBurgersVector;
using DxaSerialization::Detail::polylineLength;
using DxaSerialization::Detail::UnwrappedMesh;
using DxaSerialization::Detail::unwrapFaces;

template <typename MeshType>
json LineReconstructionJsonExporter::getMeshData(const MeshType& mesh, const StructureAnalysis& structureAnalysis) {
    json meshData;
    const auto& originalFaces = mesh.faces();
    const UnwrappedMesh unwrapped = unwrapFaces(mesh, structureAnalysis.context().simCell);

    meshData["main_listing"] = {
        {"total_nodes", static_cast<int>(unwrapped.points.size())},
        {"total_facets", static_cast<int>(unwrapped.faces.size())}
    };

    // One progress span over both JSON loops: they dominate the export of large meshes.
    ProgressReporter progress("  Serializing interface mesh", unwrapped.points.size() + unwrapped.faces.size());
    json points = json::array();
    for(size_t i = 0; i < unwrapped.points.size(); ++i) {
        progress.update(i + 1);
        const auto& pos = unwrapped.points[i];
        points.push_back({
            {"index", static_cast<int>(i)},
            {"position", {pos.x(), pos.y(), pos.z()}}
        });
    }

    json facets = json::array();
    for(size_t faceIndex = 0; faceIndex < unwrapped.faces.size(); ++faceIndex) {
        progress.update(unwrapped.points.size() + faceIndex + 1);
        facets.push_back({
            {"vertices", unwrapped.faces[faceIndex]},
            {"region", originalFaces[faceIndex] ? originalFaces[faceIndex]->region : 0}
        });
    }

    meshData["sub_listings"] = {
        {"points", points},
        {"facets", facets}
    };
    meshData["export"]["MeshExporter"]["vertices"] = points;
    meshData["export"]["MeshExporter"]["facets"] = facets;
    return meshData;
}

void LineReconstructionJsonExporter::writeLineSegmentsParquetToFile(
    const std::vector<DXA::LineReconstructionSegment>& segments,
    const SimulationCell& simulationCell,
    const std::string& filePath
) {
    struct Row {
        std::vector<Point3> points;
        double length;
        const DXA::LineReconstructionSegment* segment;
    };
    std::vector<Row> rows;
    rows.reserve(segments.size());
    for(const auto& segment : segments) {
        for(auto& points : DxaSerialization::splitAtPeriodicBoundaries({segment.position1, segment.position2}, simulationCell)) {
            const double length = polylineLength(points);
            rows.push_back({std::move(points), length, &segment});
        }
    }

    streamLinesToParquet(
        filePath,
        rows.size(),
        [&](std::size_t i, std::vector<Point3>& outPoints) { outPoints = rows[i].points; },
        [&](ColumnarLineWriter& writer, std::size_t i) {
            const auto& row = rows[i];
            const auto& segment = *row.segment;
            const Vector3 burgersLocal = segment.burgersVector.localVec();
            const Vector3 burgersGlobal = getGlobalBurgersVector(segment.burgersVector);
            const std::string structureName = structureTypeName(segment.structureType);
            const auto family = DxaSerialization::classifyBurgersFamily(burgersLocal, structureName);
            writer.field("length", row.length);
            writer.field("num_points", row.points.size());
            writer.field("magnitude", burgersLocal.length());
            writer.field("burgers_vector_local", std::vector<double>{burgersLocal.x(), burgersLocal.y(), burgersLocal.z()});
            writer.field("burgers_vector_global", std::vector<double>{burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()});
            writer.field("crystal_structure", structureName);
            writer.field("burgers_family", family.name);
            writer.field("burgers_family_label", family.label);
            writer.field("cluster_id", segment.clusterId);
            writer.field("stage", segment.stage);
        }
    );
}

void LineReconstructionJsonExporter::writeDislocationLinesParquetToFile(
    const std::vector<DXA::LineReconstructionDislocationLine>& lines,
    const SimulationCell& simulationCell,
    const std::string& filePath
) {
    struct Row {
        std::vector<Point3> points;
        double length;
        const DXA::LineReconstructionDislocationLine* line;
    };
    std::vector<Row> rows;
    rows.reserve(lines.size());
    for(const auto& line : lines) {
        for(auto& points : DxaSerialization::splitAtPeriodicBoundaries(line.points, simulationCell)) {
            const double length = polylineLength(points);
            rows.push_back({std::move(points), length, &line});
        }
    }

    streamLinesToParquet(
        filePath,
        rows.size(),
        [&](std::size_t i, std::vector<Point3>& outPoints) { outPoints = rows[i].points; },
        [&](ColumnarLineWriter& writer, std::size_t i) {
            const auto& row = rows[i];
            const auto& line = *row.line;
            const Vector3 burgersLocal = line.burgersVector.localVec();
            const Vector3 burgersGlobal = getGlobalBurgersVector(line.burgersVector);
            const std::string structureName = structureTypeName(line.structureType);
            const auto family = DxaSerialization::classifyBurgersFamily(burgersLocal, structureName);
            writer.field("length", row.length);
            writer.field("num_points", row.points.size());
            writer.field("magnitude", burgersLocal.length());
            writer.field("burgers_vector_local", std::vector<double>{burgersLocal.x(), burgersLocal.y(), burgersLocal.z()});
            writer.field("burgers_vector_global", std::vector<double>{burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()});
            writer.field("crystal_structure", structureName);
            writer.field("burgers_family", family.name);
            writer.field("burgers_family_label", family.label);
            writer.field("cluster_id", line.clusterId);
            writer.field("is_closed", line.isClosed);
            writer.field("dislocation_type_id", line.dislocationTypeId);
        }
    );
}

json LineReconstructionJsonExporter::exportUnassignedEdgesToJson(
    const std::vector<DXA::LineReconstructionUnassignedEdge>& edges
) {
    json edgeArray = json::array();
    for(size_t index = 0; index < edges.size(); ++index) {
        const auto& edge = edges[index];
        edgeArray.push_back({
            {"edge_id", static_cast<int>(index)},
            {"position1", {edge.position1.x(), edge.position1.y(), edge.position1.z()}},
            {"position2", {edge.position2.x(), edge.position2.y(), edge.position2.z()}},
            {"atoms", {edge.atom1, edge.atom2}},
            {"stage", edge.stage}
        });
    }

    json data;
    data["main_listing"] = {{"unassigned_edges", static_cast<int>(edges.size())}};
    data["sub_listings"] = {{"unassigned_edges", edgeArray}};
    return data;
}

void LineReconstructionJsonExporter::writeUnassignedEdgesParquetToFile(
    const std::vector<DXA::LineReconstructionUnassignedEdge>& edges,
    const std::string& filePath
) {
    JsonUtils::writeJsonToParquet(exportUnassignedEdgesToJson(edges), filePath, false);
}

void LineReconstructionJsonExporter::writeInterfaceMeshParquetToFile(
    const DXA::LineReconstructionInterfaceMesh& mesh,
    const StructureAnalysis& structureAnalysis,
    const std::string& filePath
) {
    // The two take about equally long on a sample run.
    ProgressSteps steps({
        {"Serializing the interface mesh", 1.0},
        {"Writing the interface mesh to Parquet", 1.0},
    });
    const json meshData = getMeshData(mesh, structureAnalysis);
    steps.next();
    JsonUtils::writeJsonToParquet(meshData, filePath, false);
}

json LineReconstructionJsonExporter::vectorToJson(const Vector3& vector) {
    return json{{"x", vector.x()}, {"y", vector.y()}, {"z", vector.z()}};
}

json LineReconstructionJsonExporter::affineTransformationToJson(const AffineTransformation& transform) {
    json transformJson = json::array();
    for(size_t i = 0; i < 3; ++i) {
        json row = json::array();
        for(size_t j = 0; j < 3; ++j) {
            row.push_back(transform(i, j));
        }
        transformJson.push_back(row);
    }
    return transformJson;
}

json LineReconstructionJsonExporter::simulationCellToJson(const SimulationCell& cell) {
    json cellJson;
    cellJson["matrix"] = affineTransformationToJson(cell.matrix());
    cellJson["volume"] = cell.volume3D();
    cellJson["is_2d"] = cell.is2D();
    const Vector3 a = cell.matrix().column(0);
    const Vector3 b = cell.matrix().column(1);
    const Vector3 c = cell.matrix().column(2);
    cellJson["lattice_vectors"] = {
        {"a", vectorToJson(a)},
        {"b", vectorToJson(b)},
        {"c", vectorToJson(c)}
    };
    cellJson["lattice_parameters"] = {
        {"a_length", a.length()},
        {"b_length", b.length()},
        {"c_length", c.length()}
    };
    return cellJson;
}

json LineReconstructionJsonExporter::getExtendedSimulationCellInfo(const SimulationCell& cell) {
    json cellJson = simulationCellToJson(cell);
    const auto& pbcFlags = cell.pbcFlags();
    cellJson["periodic_boundary_conditions"] = {
        {"x", pbcFlags[0]},
        {"y", pbcFlags[1]},
        {"z", pbcFlags[2]}
    };
    const Vector3 a = cell.matrix().column(0);
    const Vector3 b = cell.matrix().column(1);
    const Vector3 c = cell.matrix().column(2);
    cellJson["angles"] = {
        {"alpha", calculateAngle(b, c)},
        {"beta", calculateAngle(a, c)},
        {"gamma", calculateAngle(a, b)}
    };
    cellJson["reciprocal_lattice"] = {
        {"matrix", affineTransformationToJson(cell.inverseMatrix())},
        {"volume", 1.0 / cell.volume3D()}
    };
    cellJson["dimensionality"] = {
        {"is_2d", cell.is2D()},
        {"effective_dimensions", cell.is2D() ? 2 : 3}
    };
    return wrapSimulationCellInfoJson(cell, std::move(cellJson));
}

json LineReconstructionJsonExporter::wrapSimulationCellInfoJson(const SimulationCell& cell, json&& cellJson) {
    const auto& pbcFlags = cell.pbcFlags();
    const Vector3 a = cell.matrix().column(0);
    const Vector3 b = cell.matrix().column(1);
    const Vector3 c = cell.matrix().column(2);
    json cellListing = json::array({
        json{
            {"volume", cell.volume3D()},
            {"is_2d", cell.is2D()},
            {"effective_dimensions", cell.is2D() ? 2 : 3},
            {"a_x", a.x()},
            {"a_y", a.y()},
            {"a_z", a.z()},
            {"b_x", b.x()},
            {"b_y", b.y()},
            {"b_z", b.z()},
            {"c_x", c.x()},
            {"c_y", c.y()},
            {"c_z", c.z()},
            {"a_length", a.length()},
            {"b_length", b.length()},
            {"c_length", c.length()},
            {"pbc_x", pbcFlags[0]},
            {"pbc_y", pbcFlags[1]},
            {"pbc_z", pbcFlags[2]}
        }
    });
    return {
        {"export", {{"SimulationCellExporter", {{"simulation_cell", cellJson}}}}},
        {"main_listing", {
            {"simulation_cells", 1},
            {"volume", cell.volume3D()},
            {"is_2d", cell.is2D()},
            {"effective_dimensions", cell.is2D() ? 2 : 3},
            {"periodic_dimensions", static_cast<int>(pbcFlags[0]) + static_cast<int>(pbcFlags[1]) + static_cast<int>(pbcFlags[2])}
        }},
        {"sub_listings", {{"simulation_cell", std::move(cellListing)}}}
    };
}

double LineReconstructionJsonExporter::calculateAngle(const Vector3& a, const Vector3& b) {
    const double magnitudes = a.length() * b.length();
    if(magnitudes == 0.0) return 0.0;
    double cosAngle = a.dot(b) / magnitudes;
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    return std::acos(cosAngle) * 180.0 / PI;
}

}
