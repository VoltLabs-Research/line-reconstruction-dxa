#include <volt/lrdxa_json_exporter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <volt/utilities/json_utils.h>

namespace Volt {

namespace {

Vector3 getGlobalBurgersVector(const ClusterVector& burgersVector) {
    if(burgersVector.cluster() == nullptr) {
        return burgersVector.localVec();
    }
    return burgersVector.toSpatialVector();
}

json wrapSimulationCellInfoJson(const SimulationCell& cell, json&& cellJson) {
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

void clipDislocationLine(
    const std::vector<Point3>& line,
    const SimulationCell& simulationCell,
    const std::function<void(const Point3&, const Point3&, bool)>& segmentCallback
) {
    if(line.size() < 2) return;
    bool isInitialSegment = true;

    auto v1Iter = line.cbegin();
    Point3 rp1 = simulationCell.absoluteToReduced(*v1Iter);
    Vector3 shiftVector = Vector3::Zero();
    for(size_t dimension = 0; dimension < 3; ++dimension) {
        if(simulationCell.pbcFlags()[dimension]) {
            const double shift = -std::floor(rp1[dimension]);
            rp1[dimension] += shift;
            shiftVector[dimension] += shift;
        }
    }

    for(auto v2Iter = v1Iter + 1; v2Iter != line.cend(); v1Iter = v2Iter, ++v2Iter) {
        Point3 rp2 = simulationCell.absoluteToReduced(*v2Iter) + shiftVector;
        int iterationCount = 0;
        while(true) {
            if(++iterationCount > 10) {
                segmentCallback(
                    simulationCell.reducedToAbsolute(rp1),
                    simulationCell.reducedToAbsolute(rp2),
                    isInitialSegment
                );
                break;
            }

            size_t crossDim = static_cast<size_t>(-1);
            double crossDir = 0.0;
            double smallestT = std::numeric_limits<double>::max();
            for(size_t dimension = 0; dimension < 3; ++dimension) {
                if(!simulationCell.pbcFlags()[dimension]) continue;
                const int d = static_cast<int>(std::floor(rp2[dimension])) - static_cast<int>(std::floor(rp1[dimension]));
                if(d == 0) continue;
                const double dr = rp2[dimension] - rp1[dimension];
                if(std::abs(dr) < 1e-9) continue;
                const double t = (d > 0)
                    ? (std::ceil(rp1[dimension]) - rp1[dimension]) / dr
                    : (std::floor(rp1[dimension]) - rp1[dimension]) / dr;
                if(t > 1e-9 && t < smallestT) {
                    smallestT = t;
                    crossDim = dimension;
                    crossDir = (d > 0) ? 1.0 : -1.0;
                }
            }

            if(smallestT < (1.0 - 1e-9)) {
                Point3 intersection = rp1 + smallestT * (rp2 - rp1);
                intersection[crossDim] = std::round(intersection[crossDim]);
                segmentCallback(
                    simulationCell.reducedToAbsolute(rp1),
                    simulationCell.reducedToAbsolute(intersection),
                    isInitialSegment
                );
                shiftVector[crossDim] -= crossDir;
                rp1 = intersection;
                rp1[crossDim] -= crossDir;
                rp2[crossDim] -= crossDir;
                isInitialSegment = true;
            } else {
                segmentCallback(
                    simulationCell.reducedToAbsolute(rp1),
                    simulationCell.reducedToAbsolute(rp2),
                    isInitialSegment
                );
                isInitialSegment = false;
                break;
            }
        }
        rp1 = rp2;
    }
}

}

template <typename MeshType>
json LineReconstructionJsonExporter::getMeshData(const MeshType& mesh, const StructureAnalysis& structureAnalysis) {
    json meshData;
    const auto& originalVertices = mesh.vertices();
    const auto& originalFaces = mesh.faces();
    const auto& cell = structureAnalysis.context().simCell;

    std::vector<Point3> exportPoints;
    exportPoints.reserve(originalVertices.size());
    std::vector<int> originalToExportVertexMap(originalVertices.size());
    for(size_t i = 0; i < originalVertices.size(); ++i) {
        exportPoints.push_back(originalVertices[i]->pos());
        originalToExportVertexMap[i] = static_cast<int>(i);
    }

    std::vector<std::vector<int>> exportFaces;
    exportFaces.reserve(originalFaces.size());
    for(const auto* face : originalFaces) {
        if(!face || !face->edges()) continue;
        std::vector<int> faceVertexIndices;
        std::vector<Point3> faceVertexPositions;
        auto* startEdge = face->edges();
        auto* currentEdge = startEdge;
        do {
            faceVertexIndices.push_back(currentEdge->vertex1()->index());
            faceVertexPositions.push_back(currentEdge->vertex1()->pos());
            currentEdge = currentEdge->nextFaceEdge();
        } while(currentEdge != startEdge);

        cell.unwrapPositions(faceVertexPositions.data(), faceVertexPositions.size());

        std::vector<int> newFaceIndices;
        for(size_t i = 0; i < faceVertexIndices.size(); ++i) {
            const int originalIndex = faceVertexIndices[i];
            const Point3& originalPos = originalVertices[originalIndex]->pos();
            const Point3& unwrappedPos = faceVertexPositions[i];
            if(!originalPos.equals(unwrappedPos, 1e-6)) {
                newFaceIndices.push_back(static_cast<int>(exportPoints.size()));
                exportPoints.push_back(unwrappedPos);
            } else {
                newFaceIndices.push_back(originalToExportVertexMap[originalIndex]);
            }
        }
        exportFaces.push_back(std::move(newFaceIndices));
    }

    meshData["main_listing"] = {
        {"total_nodes", static_cast<int>(exportPoints.size())},
        {"total_facets", static_cast<int>(exportFaces.size())}
    };

    json points = json::array();
    for(size_t i = 0; i < exportPoints.size(); ++i) {
        const auto& pos = exportPoints[i];
        points.push_back({
            {"index", static_cast<int>(i)},
            {"position", {pos.x(), pos.y(), pos.z()}}
        });
    }

    json facets = json::array();
    for(size_t faceIndex = 0; faceIndex < exportFaces.size(); ++faceIndex) {
        facets.push_back({
            {"vertices", exportFaces[faceIndex]},
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

json LineReconstructionJsonExporter::exportLineSegmentsToJson(
    const std::vector<DXA::LineReconstructionSegment>& segments,
    const SimulationCell* simulationCell
) {
    json segmentArray = json::array();
    int globalChunkId = 0;

    for(const auto& segment : segments) {
        const Vector3 burgersLocal = segment.burgersVector.localVec();
        const Vector3 burgersGlobal = getGlobalBurgersVector(segment.burgersVector);
        if(simulationCell) {
            std::vector<Point3> line = {segment.position1, segment.position2};
            std::vector<Point3> currentChunk;
            clipDislocationLine(line, *simulationCell, [&](const Point3& p1, const Point3& p2, bool isInitialSegment) {
                if(isInitialSegment && !currentChunk.empty()) {
                    const Point3 cp1 = currentChunk.front();
                    const Point3 cp2 = currentChunk.back();
                    double chunkLength = 0.0;
                    json points = json::array();
                    for(size_t i = 0; i < currentChunk.size(); ++i) {
                        points.push_back({currentChunk[i].x(), currentChunk[i].y(), currentChunk[i].z()});
                        if(i > 0) chunkLength += (currentChunk[i] - currentChunk[i - 1]).length();
                    }
                    segmentArray.push_back({
                        {"segment_id", globalChunkId++},
                        {"points", points},
                        {"num_points", static_cast<int>(currentChunk.size())},
                        {"length", chunkLength},
                        {"position1", {cp1.x(), cp1.y(), cp1.z()}},
                        {"position2", {cp2.x(), cp2.y(), cp2.z()}},
                        {"burgers_vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                        {"burgers_vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                        {"burgers_vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                        {"burgers", {
                            {"vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                            {"vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                            {"vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                            {"magnitude", burgersLocal.length()}
                        }},
                        {"magnitude", burgersLocal.length()},
                        {"cluster_id", segment.clusterId},
                        {"structure_type", segment.structureType},
                        {"stage", segment.stage}
                    });
                    currentChunk.clear();
                }
                if(currentChunk.empty()) {
                    currentChunk.push_back(p1);
                }
                currentChunk.push_back(p2);
            });
            if(!currentChunk.empty()) {
                const Point3 cp1 = currentChunk.front();
                const Point3 cp2 = currentChunk.back();
                double chunkLength = 0.0;
                json points = json::array();
                for(size_t i = 0; i < currentChunk.size(); ++i) {
                    points.push_back({currentChunk[i].x(), currentChunk[i].y(), currentChunk[i].z()});
                    if(i > 0) chunkLength += (currentChunk[i] - currentChunk[i - 1]).length();
                }
                segmentArray.push_back({
                    {"segment_id", globalChunkId++},
                    {"points", points},
                    {"num_points", static_cast<int>(currentChunk.size())},
                    {"length", chunkLength},
                    {"position1", {cp1.x(), cp1.y(), cp1.z()}},
                    {"position2", {cp2.x(), cp2.y(), cp2.z()}},
                    {"burgers_vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"burgers_vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"burgers_vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                    {"burgers", {
                        {"vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                        {"vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                        {"vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                        {"magnitude", burgersLocal.length()}
                    }},
                    {"magnitude", burgersLocal.length()},
                    {"cluster_id", segment.clusterId},
                    {"structure_type", segment.structureType},
                    {"stage", segment.stage}
                });
            }
        } else {
            const double chunkLength = (segment.position2 - segment.position1).length();
            segmentArray.push_back({
                {"segment_id", globalChunkId++},
                {"points", {{segment.position1.x(), segment.position1.y(), segment.position1.z()}, {segment.position2.x(), segment.position2.y(), segment.position2.z()}}},
                {"num_points", 2},
                {"length", chunkLength},
                {"position1", {segment.position1.x(), segment.position1.y(), segment.position1.z()}},
                {"position2", {segment.position2.x(), segment.position2.y(), segment.position2.z()}},
                {"burgers_vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                {"burgers_vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                {"burgers_vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                {"burgers", {
                    {"vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                    {"magnitude", burgersLocal.length()}
                }},
                {"magnitude", burgersLocal.length()},
                {"cluster_id", segment.clusterId},
                {"structure_type", segment.structureType},
                {"stage", segment.stage}
            });
        }
    }

    json data;
    data["main_listing"] = {
        {"dislocation_segments", static_cast<int>(segmentArray.size())},
        {"dislocations", static_cast<int>(segmentArray.size())}
    };
    data["sub_listings"] = {{"dislocation_segments", segmentArray}};
    data["export"]["DislocationExporter"]["segments"] = segmentArray;
    data["export"]["DislocationExporter"]["dislocation_segments"] = segmentArray;
    return data;
}

void LineReconstructionJsonExporter::writeLineSegmentsParquetToFile(
    const std::vector<DXA::LineReconstructionSegment>& segments,
    const SimulationCell& simulationCell,
    const std::string& filePath
) {
    JsonUtils::writeJsonToParquet(exportLineSegmentsToJson(segments, &simulationCell), filePath, false);
}

json LineReconstructionJsonExporter::exportDislocationLinesToJson(
    const std::vector<DXA::LineReconstructionDislocationLine>& lines,
    const SimulationCell* simulationCell
) {
    json lineArray = json::array();
    int globalChunkId = 0;
    int totalPoints = 0;
    double totalLength = 0.0;

    for(const auto& line : lines) {
        const Vector3 burgersLocal = line.burgersVector.localVec();
        const Vector3 burgersGlobal = getGlobalBurgersVector(line.burgersVector);
        auto appendLine = [&](const std::vector<Point3>& pointsVec) {
            if(pointsVec.size() < 2) {
                return;
            }

            double lineLength = 0.0;
            json points = json::array();
            for(size_t i = 0; i < pointsVec.size(); ++i) {
                points.push_back({pointsVec[i].x(), pointsVec[i].y(), pointsVec[i].z()});
                if(i > 0) {
                    lineLength += (pointsVec[i] - pointsVec[i - 1]).length();
                }
            }
            if(line.isClosed && pointsVec.size() > 1) {
                lineLength += (pointsVec.front() - pointsVec.back()).length();
            }

            lineArray.push_back({
                {"line_id", globalChunkId++},
                {"points", points},
                {"num_points", static_cast<int>(pointsVec.size())},
                {"length", lineLength},
                {"burgers_vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                {"burgers_vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                {"burgers_vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                {"burgers", {
                    {"vector", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"vector_local", {burgersLocal.x(), burgersLocal.y(), burgersLocal.z()}},
                    {"vector_global", {burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z()}},
                    {"magnitude", burgersLocal.length()}
                }},
                {"magnitude", burgersLocal.length()},
                {"cluster_id", line.clusterId},
                {"structure_type", line.structureType},
                {"is_closed", line.isClosed},
                {"dislocation_type_id", line.dislocationTypeId}
            });

            totalPoints += static_cast<int>(pointsVec.size());
            totalLength += lineLength;
        };

        if(simulationCell) {
            std::vector<Point3> currentChunk;
            clipDislocationLine(line.points, *simulationCell, [&](const Point3& p1, const Point3& p2, bool isInitialSegment) {
                if(isInitialSegment && !currentChunk.empty()) {
                    appendLine(currentChunk);
                    currentChunk.clear();
                }
                if(currentChunk.empty()) {
                    currentChunk.push_back(p1);
                }
                currentChunk.push_back(p2);
            });
            if(!currentChunk.empty()) {
                appendLine(currentChunk);
            }
        } else {
            appendLine(line.points);
        }
    }

    json data;
    data["main_listing"] = {
        {"dislocations", static_cast<int>(lineArray.size())},
        {"total_points", totalPoints},
        {"total_length", totalLength}
    };
    data["sub_listings"] = {{"dislocations", lineArray}};
    data["export"]["DislocationExporter"]["segments"] = lineArray;
    data["export"]["DislocationExporter"]["dislocations"] = lineArray;
    return data;
}

void LineReconstructionJsonExporter::writeDislocationLinesParquetToFile(
    const std::vector<DXA::LineReconstructionDislocationLine>& lines,
    const SimulationCell& simulationCell,
    const std::string& filePath
) {
    JsonUtils::writeJsonToParquet(exportDislocationLinesToJson(lines, &simulationCell), filePath, false);
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
    data["export"]["DislocationExporter"]["unassigned_edges"] = edgeArray;
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
    JsonUtils::writeJsonToParquet(getMeshData(mesh, structureAnalysis), filePath, false);
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

double LineReconstructionJsonExporter::calculateAngle(const Vector3& a, const Vector3& b) {
    const double magnitudes = a.length() * b.length();
    if(magnitudes == 0.0) return 0.0;
    double cosAngle = a.dot(b) / magnitudes;
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    return std::acos(cosAngle) * 180.0 / PI;
}

}
