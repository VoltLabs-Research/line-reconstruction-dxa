#include <volt/lrdxa_engine.h>

#include <volt/helpers/crystal_path_finder.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/structures/crystal_structure_types.h>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <array>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Volt::DXA {

namespace {

bool isFiniteCell(const DelaunayTessellation& tessellation, DelaunayTessellation::CellHandle cell) {
    return tessellation.isValidCell(cell);
}

std::array<DelaunayTessellation::VertexHandle, 4> getCellVertices(
    const DelaunayTessellation& tessellation,
    DelaunayTessellation::CellHandle cell) {
    return {
        tessellation.cellVertex(cell, 0),
        tessellation.cellVertex(cell, 1),
        tessellation.cellVertex(cell, 2),
        tessellation.cellVertex(cell, 3),
    };
}

int getInputPointIndex(const DelaunayTessellation& tessellation, DelaunayTessellation::VertexHandle vertex) {
    return tessellation.vertexIndex(vertex);
}

int findInputPointInCell(
    const DelaunayTessellation& tessellation,
    DelaunayTessellation::CellHandle cell,
    int atomIndex) {
    for(int localVertex = 0; localVertex < 4; ++localVertex) {
        if(getInputPointIndex(tessellation, tessellation.cellVertex(cell, localVertex)) == atomIndex) {
            return localVertex;
        }
    }
    return -1;
}

int computedSubTypeForCrystal(int latticeStructureType) {
    switch(latticeStructureType) {
        case LATTICE_FCC:
            return LATTICE_HCP;
        case LATTICE_HCP:
            return LATTICE_FCC;
        case LATTICE_CUBIC_DIAMOND:
            return LATTICE_HEX_DIAMOND;
        case LATTICE_HEX_DIAMOND:
            return LATTICE_CUBIC_DIAMOND;
        default:
            return LATTICE_OTHER;
    }
}

bool crystalUsesDirectTessellation(int latticeStructureType) {
    return computedSubTypeForCrystal(latticeStructureType) == LATTICE_OTHER;
}

}

LineReconstructionDXAAlgorithm::LineReconstructionDXAAlgorithm(StructureAnalysis& structureAnalysis, StructureContext& context)
    : _structureAnalysis(structureAnalysis)
    , _context(context)
    , _atomClusters(context.atomCount(), nullptr)
    , _atomOutboundEdges(context.atomCount(), nullptr)
    , _atomInboundEdges(context.atomCount(), nullptr) {
}

ClusterVector LineReconstructionDXAAlgorithm::OrientedEdge::vector() const {
    if(!_edge || !_edge->hasEdgeVector) {
        return ClusterVector{};
    }
    if(!_flipped) {
        return _edge->vector;
    }
    ClusterTransition* transition = _edge->transition;
    if(!transition) {
        return ClusterVector{};
    }
    return ClusterVector(-transition->transform(_edge->vector.localVec()), transition->cluster2);
}

size_t LineReconstructionDXAAlgorithm::AtomIndexTripletHash::operator()(const std::array<AtomIndex, 3>& arr) const noexcept {
    size_t seed = std::hash<AtomIndex>{}(arr[0]);
    seed ^= std::hash<AtomIndex>{}(arr[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<AtomIndex>{}(arr[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

void LineReconstructionDXAAlgorithm::run(
    const LineReconstructionDXAOptions& options,
    const std::function<void(std::string_view)>& markStage
) {
    const auto emitStage = [&](std::string_view stage) {
        if(markStage) {
            markStage(stage);
        }
    };

    markPerfectCrystallineRegions();
    emitStage("mark_perfect_crystalline_regions");

    const double ghostLayerSize = options.tessellationGhostLayerScale * _structureAnalysis.maximumNeighborDistance();
    buildTessellation(ghostLayerSize);
    emitStage("delaunay_tessellation");

    const double alpha = options.alphaScale * _structureAnalysis.maximumNeighborDistance();
    classifyTetrahedra(alpha);
    emitStage("classify_tetrahedra");

    generateTessellationEdges();
    emitStage("generate_tessellation_edges");

    assignAtomsToClusters();
    emitStage("assign_atoms_to_clusters");

    assignIdealVectorsToEdges(options.crystalPathSteps);
    emitStage("assign_ideal_vectors_to_edges");

    buildFacetLookupMap();
    emitStage("build_facet_lookup_map");

    _dislocationSegments.clear();
    _unassignedEdgesOutput.clear();
    _interfaceMesh.clear();
    _currentStage = 0;

    dumpDislocationSnapshot(alpha);
    emitStage("snapshot_initial_mapping");

    eliminateSpuriousDislocationLoops();
    emitStage("eliminate_spurious_dislocation_loops");

    dumpDislocationSnapshot(alpha);
    emitStage("snapshot_after_loop_cleanup");

    const bool skipComplementPasses = crystalUsesDirectTessellation(_context.inputCrystalType)
        && _structureAnalysis.clusterGraph().clusterTransitions().empty();

    if(skipComplementPasses) {
        emitStage("skip_complement_edge_vectors");
    } else {
        while(true) {
            bool changed = false;

            changed |= complementEdgeVectors(true);
            emitStage("complement_edge_vectors_forward");

            dumpDislocationSnapshot(alpha);
            emitStage("snapshot_after_forward_complement");

            changed |= complementEdgeVectors(false);
            emitStage("complement_edge_vectors_backward");

            dumpDislocationSnapshot(alpha);
            emitStage("snapshot_after_backward_complement");

            if(!changed) break;
        }
    }

    dumpDislocationSnapshot(alpha);
    emitStage("snapshot_final_mapping");

    stitchDislocationLines();
    emitStage("stitch_dislocation_lines");

    finishDislocationLines(options.smoothingIterations, options.linePointInterval);
    emitStage("finish_dislocation_lines");

    spdlog::info(
        "line-reconstruction-dxa produced {} snapshot segments, {} stitched lines, {} unresolved edges",
        _dislocationSegments.size(),
        _dislocationLines.size(),
        _unassignedEdgesOutput.size()
    );
}

void LineReconstructionDXAAlgorithm::markPerfectCrystallineRegions() {
    const size_t atomCount = _context.atomCount();
    const int computedSubType = computedSubTypeForCrystal(_context.inputCrystalType);
    const bool directCrystalTessellation = crystalUsesDirectTessellation(_context.inputCrystalType);
    _delaunayAtomMask = std::make_shared<ParticleProperty>(atomCount, DataType::Int, 1, 0, true);
    _delaunayVertexCount = 0;
    for(size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
        Cluster* cluster = _structureAnalysis.atomCluster(static_cast<int>(atomIndex));
        if(!cluster || cluster->id == 0) {
            _delaunayAtomMask->setInt(atomIndex, 1);
            _delaunayVertexCount++;
            continue;
        }

        const int structureType = cluster->structure;
        if(directCrystalTessellation && structureType == _context.inputCrystalType) {
            _delaunayAtomMask->setInt(atomIndex, 1);
            _delaunayVertexCount++;
            continue;
        }
        if(structureType == LATTICE_OTHER || (computedSubType != LATTICE_OTHER && structureType == computedSubType)) {
            _delaunayAtomMask->setInt(atomIndex, 1);
            _delaunayVertexCount++;
            continue;
        }

        int includeAtom = 0;
        const int numNeighbors = _structureAnalysis.numberOfNeighbors(static_cast<int>(atomIndex));
        for(int neighborIndex = 0; neighborIndex < numNeighbors; ++neighborIndex) {
            const int neighborAtomIndex = _structureAnalysis.getNeighbor(static_cast<int>(atomIndex), neighborIndex);
            Cluster* neighborCluster = _structureAnalysis.atomCluster(neighborAtomIndex);
            if(!neighborCluster || neighborCluster->id == 0) {
                includeAtom = 1;
                break;
            }
            if(computedSubType != LATTICE_OTHER && neighborCluster->structure == computedSubType) {
                includeAtom = 1;
                break;
            }
        }
        _delaunayAtomMask->setInt(atomIndex, includeAtom);
        _delaunayVertexCount += includeAtom;
    }
}

void LineReconstructionDXAAlgorithm::buildTessellation(double ghostLayerSize) {
    _tessellation.generateTessellation(
        _context.simCell,
        _context.positions->constDataPoint3(),
        _context.atomCount(),
        ghostLayerSize,
        false,
        _delaunayAtomMask ? _delaunayAtomMask->constDataInt() : nullptr
    );
}

void LineReconstructionDXAAlgorithm::classifyTetrahedra(double alpha) {
    _filledCells.assign(_tessellation.numberOfTetrahedra(), false);
    auto isInteriorTetrahedron = [&](DelaunayTessellation::CellHandle cell) -> bool {
        if(_tessellation.isGhostCell(cell)) {
            return false;
        }
        if(auto alphaTestResult = _tessellation.alphaTest(cell, alpha)) {
            return *alphaTestResult;
        }
        for(int facet = 0; facet < 4; ++facet) {
            DelaunayTessellation::CellHandle adjacentCell = _tessellation.mirrorFacet(cell, facet).first;
            if(!isFiniteCell(_tessellation, adjacentCell)) {
                return false;
            }
            auto adjacentAlphaTestResult = _tessellation.alphaTest(adjacentCell, alpha);
            if(adjacentAlphaTestResult.has_value() && !adjacentAlphaTestResult.value()) {
                return false;
            }
        }
        return true;
    };

    // Per-cell, read-only alpha test (alphaTest/mirrorFacet/isGhostCell are
    // const on the tessellation and called concurrently in opendxa too); each
    // cell writes a distinct _filledCells[cell], so this parallelizes cleanly.
    const size_t numCells = _tessellation.numberOfTetrahedra();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numCells, 4096),
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t cell = r.begin(); cell < r.end(); ++cell){
                _filledCells[cell] = isInteriorTetrahedron(static_cast<DelaunayTessellation::CellHandle>(cell)) ? 1 : 0;
            }
        });
}

void LineReconstructionDXAAlgorithm::generateTessellationEdges() {
    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isFilledCell(cell)) {
            continue;
        }
        const auto cellVertices = getCellVertices(_tessellation, cell);
        std::array<AtomIndex, 4> cellAtoms{};
        std::array<Point3, 4> cellPositions{};
        for(int vertex = 0; vertex < 4; ++vertex) {
            cellAtoms[vertex] = getInputPointIndex(_tessellation, cellVertices[vertex]);
            cellPositions[vertex] = _tessellation.vertexPosition(cellVertices[vertex]);
        }
        for(int edgeIndex = 0; edgeIndex < 6; ++edgeIndex) {
            const int localVertex1 = kCellEdgeVertices[edgeIndex][0];
            const int localVertex2 = kCellEdgeVertices[edgeIndex][1];
            AtomIndex atom1 = cellAtoms[localVertex1];
            AtomIndex atom2 = cellAtoms[localVertex2];
            if(atom1 == atom2) {
                continue;
            }
            if(_context.simCell.isWrappedVector(cellPositions[localVertex2] - cellPositions[localVertex1])) {
                continue;
            }
            if(findEdge(atom1, atom2)) {
                continue;
            }
            TessellationEdge* edge = _edgePool.construct(atom1, atom2, cell);
            edge->nextOutboundEdge = std::exchange(_atomOutboundEdges[atom1], edge);
            edge->nextInboundEdge = std::exchange(_atomInboundEdges[atom2], edge);
        }
    }
}

void LineReconstructionDXAAlgorithm::assignAtomsToClusters() {
    std::vector<AtomIndex> unassignedAtoms;
    for(AtomIndex atomIndex = 0; atomIndex < static_cast<AtomIndex>(_context.atomCount()); ++atomIndex) {
        Cluster* cluster = _structureAnalysis.atomCluster(atomIndex);
        _atomClusters[atomIndex] = cluster;
        if(!cluster || cluster->structure == LATTICE_OTHER) {
            unassignedAtoms.push_back(atomIndex);
        }
    }

    while(!unassignedAtoms.empty()) {
        size_t remaining = 0;
        for(AtomIndex atomIndex : unassignedAtoms) {
            Cluster*& assignedCluster = _atomClusters[atomIndex];
            for(TessellationEdge* edge = _atomOutboundEdges[atomIndex]; edge != nullptr; edge = edge->nextOutboundEdge) {
                Cluster* cluster = clusterOfAtom(edge->atom2);
                if(cluster && cluster->structure != LATTICE_OTHER) {
                    assignedCluster = cluster;
                    break;
                }
            }
            if(assignedCluster && assignedCluster->structure != LATTICE_OTHER) {
                continue;
            }
            for(TessellationEdge* edge = _atomInboundEdges[atomIndex]; edge != nullptr; edge = edge->nextInboundEdge) {
                Cluster* cluster = clusterOfAtom(edge->atom1);
                if(cluster && cluster->structure != LATTICE_OTHER) {
                    assignedCluster = cluster;
                    break;
                }
            }
            if(assignedCluster && assignedCluster->structure != LATTICE_OTHER) {
                continue;
            }
            unassignedAtoms[remaining++] = atomIndex;
        }
        if(remaining == unassignedAtoms.size()) {
            break;
        }
        unassignedAtoms.resize(remaining);
    }
}

void LineReconstructionDXAAlgorithm::eliminateSpuriousDislocationLoops() {
    std::unordered_map<TessellationEdge*, int> compatibleCount;
    std::unordered_map<TessellationEdge*, int> incompatibleCount;

    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isFilledCell(cell)) {
            continue;
        }
        auto cellEdges = getOrientedEdges(cell);

        bool allAssigned = true;
        for(int i = 0; i < 6; ++i) {
            if(!cellEdges[i] || !cellEdges[i].hasEdgeVector()) {
                allAssigned = false;
                break;
            }
        }
        if(!allAssigned) {
            continue;
        }

        bool compatible = isElasticMappingCompatible(cellEdges);

        for(int i = 0; i < 6; ++i) {
            TessellationEdge* edge = cellEdges[i].undirectedEdge();
            if(!edge) continue;
            if(compatible) {
                compatibleCount[edge]++;
            } else {
                incompatibleCount[edge]++;
            }
        }
    }

    int invalidated = 0;
    for(auto& [edge, incompat] : incompatibleCount) {
        if(compatibleCount.count(edge) > 0) {
            continue;
        }
        if(!edge->hasEdgeVector) {
            continue;
        }
        edge->hasEdgeVector = false;
        _unassignedEdges.push_back(edge);
        invalidated++;
    }

    if(invalidated > 0) {
        spdlog::info("eliminateSpuriousDislocationLoops: invalidated {} edges", invalidated);
    }
}

void LineReconstructionDXAAlgorithm::assignIdealVectorsToEdges(int crystalPathSteps) {
    CrystalPathFinder pathFinder(_structureAnalysis, crystalPathSteps);

    for(size_t atomIndex = 0; atomIndex < _atomOutboundEdges.size(); ++atomIndex) {
        for(TessellationEdge* edge = _atomOutboundEdges[atomIndex]; edge != nullptr; edge = edge->nextOutboundEdge) {
            Cluster* cluster1 = clusterOfAtom(edge->atom1);
            Cluster* cluster2 = clusterOfAtom(edge->atom2);
            if(!cluster1 || !cluster2) {
                continue;
            }
            if(cluster1->structure == LATTICE_OTHER || cluster2->structure == LATTICE_OTHER) {
                continue;
            }
            edge->transition = _structureAnalysis.clusterGraph().determineClusterTransition(cluster1, cluster2);
            if(!edge->transition) {
                continue;
            }
            std::optional<ClusterVector> idealVector = pathFinder.findPath(edge->atom1, edge->atom2);
            if(idealVector.has_value() && idealVector->transformToCluster(cluster1, _structureAnalysis.clusterGraph())) {
                edge->vector = *idealVector;
                edge->hasEdgeVector = true;
            } else {
                _unassignedEdges.push_back(edge);
            }
        }
    }
}

void LineReconstructionDXAAlgorithm::buildFacetLookupMap() {
    _primaryFacetLookupMap.clear();
    _primaryFacetLookupMap.reserve(_unassignedEdges.size() * 6);
    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isFilledCell(cell)) {
            continue;
        }
        const auto cellVertices = getCellVertices(_tessellation, cell);
        std::array<AtomIndex, 4> cellAtoms{};
        for(int vertex = 0; vertex < 4; ++vertex) {
            cellAtoms[vertex] = getInputPointIndex(_tessellation, cellVertices[vertex]);
        }
        for(int facet = 0; facet < 4; ++facet) {
            std::array<AtomIndex, 3> facetVertices = {
                cellAtoms[DelaunayTessellation::cellFacetVertexIndex(facet, 0)],
                cellAtoms[DelaunayTessellation::cellFacetVertexIndex(facet, 1)],
                cellAtoms[DelaunayTessellation::cellFacetVertexIndex(facet, 2)]
            };
            TessellationEdge* e1 = findEdge(facetVertices[0], facetVertices[1]);
            TessellationEdge* e2 = findEdge(facetVertices[1], facetVertices[2]);
            TessellationEdge* e3 = findEdge(facetVertices[2], facetVertices[0]);
            if(!e1 || !e2 || !e3) {
                continue;
            }
            if(e1->hasEdgeVector && e2->hasEdgeVector && e3->hasEdgeVector) {
                continue;
            }
            reorderFacetVertices(facetVertices);
            _primaryFacetLookupMap.emplace(facetVertices, std::make_pair(cell, facet));
        }
    }
}

bool LineReconstructionDXAAlgorithm::complementEdgeVectors(bool forward) {
    size_t numAssignedEdges = 0;
    constexpr int tab_next_around_edge[4][4] = {
        {5, 2, 3, 1},
        {3, 5, 0, 2},
        {1, 3, 5, 0},
        {2, 0, 1, 5}
    };

    auto visitEdgeIncidentFacets = [&](TessellationEdge* edge, auto&& visitor) {
        AtomIndex atom1 = edge->atom1;
        AtomIndex atom2 = edge->atom2;
        DelaunayTessellation::CellHandle startCell = edge->adjacentCell;
        DelaunayTessellation::CellHandle cell = startCell;
        int iterLimit = 1000;
        do {
            if(--iterLimit <= 0) break;
            int localVertex1 = findInputPointInCell(_tessellation, cell, atom1);
            int localVertex2 = findInputPointInCell(_tessellation, cell, atom2);
            if(localVertex1 < 0 || localVertex2 < 0) {
                break;
            }
            int nextFacet = tab_next_around_edge[localVertex1][localVertex2];
            if(nextFacet < 0 || nextFacet > 3) {
                break;
            }
            int facet;
            std::tie(cell, facet) = _tessellation.mirrorFacet(cell, nextFacet);
            if(!isFiniteCell(_tessellation, cell)) {
                break;
            }
            if(_tessellation.isGhostCell(cell)) {
                std::array<AtomIndex, 3> vertices{};
                for(int i = 0; i < 3; ++i) {
                    vertices[i] = getInputPointIndex(_tessellation, _tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(facet, i)));
                }
                reorderFacetVertices(vertices);
                auto iter = _primaryFacetLookupMap.find(vertices);
                if(iter == _primaryFacetLookupMap.end()) {
                    continue;
                }
                std::tie(cell, facet) = iter->second;
            }
            if(!isFilledCell(cell)) {
                continue;
            }
            visitor(cell, facet);
        } while(cell != startCell);
    };

    class UniqueVectorsList {
    public:
        void add(const ClusterVector& item) {
            for(int i = 0; i < _numItems; ++i) {
                if(_items[i].localVec().equals(item.localVec(), CA_LATTICE_VECTOR_EPSILON)) {
                    _counts[i]++;
                    return;
                }
            }
            if(_numItems < 8) {
                _items[_numItems] = item;
                _counts[_numItems] = 1;
                _numItems++;
            }
        }
        bool empty() const { return _numItems == 0; }
        const ClusterVector& mostFrequent() const {
            int maxIndex = 0;
            for(int i = 1; i < _numItems; ++i) {
                if(_counts[i] > _counts[maxIndex]) {
                    maxIndex = i;
                }
            }
            return _items[maxIndex];
        }

    private:
        ClusterVector _items[8];
        int _counts[8] = {};
        int _numItems = 0;
    };

    std::deque<TessellationEdge*> queue;
    for(TessellationEdge* unresolvedEdge : _unassignedEdges) {
        if(unresolvedEdge->hasEdgeVector || !unresolvedEdge->transition) {
            continue;
        }
        queue.clear();
        queue.push_back(unresolvedEdge);
        while(!queue.empty()) {
            TessellationEdge* edge = queue.front();
            queue.pop_front();
            if(edge->hasEdgeVector) {
                continue;
            }
            UniqueVectorsList candidates;
            OrientedEdge secondaryCandidate;
            visitEdgeIncidentFacets(edge, [&](DelaunayTessellation::CellHandle cell, int facet) {
                std::array<OrientedEdge, 3> facetEdges = getFacetCircuitEdges(cell, facet);
                int rotations = 0;
                while(facetEdges[0].atom2() != edge->atom1 && rotations < 3) {
                    std::rotate(facetEdges.begin(), facetEdges.begin() + 1, facetEdges.end());
                    rotations++;
                }
                if(rotations >= 3) return;
                if(facetEdges[1].hasEdgeVector() && facetEdges[2].hasEdgeVector()) {
                    if(!facetEdges[0].transition() || !facetEdges[2].transition()) return;
                    Vector3 inferred = -facetEdges[0].transition()->reverseTransform(facetEdges[1].vector().localVec())
                        - facetEdges[2].transition()->transform(facetEdges[2].vector().localVec());
                    candidates.add(ClusterVector(inferred, facetEdges[0].transition()->cluster1));
                } else if(forward && facetEdges[1].hasEdgeVector()) {
                    secondaryCandidate = facetEdges[1];
                    if(facetEdges[2].undirectedEdge() && facetEdges[2].undirectedEdge()->transition) {
                        queue.push_back(facetEdges[2].undirectedEdge());
                    }
                } else if(forward && facetEdges[2].hasEdgeVector()) {
                    secondaryCandidate = -facetEdges[2];
                    if(facetEdges[1].undirectedEdge() && facetEdges[1].undirectedEdge()->transition) {
                        queue.push_back(facetEdges[1].undirectedEdge());
                    }
                }
            });

            if(!candidates.empty()) {
                edge->vector = ClusterVector(-edge->transition->reverseTransform(candidates.mostFrequent().localVec()), edge->transition->cluster1);
            } else if(forward && secondaryCandidate) {
                AtomIndex oppositeAtom = secondaryCandidate.atom1() == edge->atom1 ? edge->atom2 : edge->atom1;
                if(isVertexFullyIsolated(oppositeAtom)) {
                    if(secondaryCandidate.atom1() == edge->atom1) {
                        edge->vector = secondaryCandidate.vector();
                    } else {
                        edge->vector = ClusterVector(edge->transition->reverseTransform(secondaryCandidate.vector().localVec()), edge->transition->cluster1);
                    }
                } else {
                    queue.clear();
                    continue;
                }
            } else {
                queue.clear();
                break;
            }
            edge->hasEdgeVector = true;
            numAssignedEdges++;
        }
    }
    return numAssignedEdges != 0;
}

std::array<bool, 6> LineReconstructionDXAAlgorithm::determineConservativeGlideEdgesFCC(
    const std::array<OrientedEdge, 6>& cellEdges,
    const ClusterVector facetBurgers[4],
    const bool activeFacet[4])
{
    std::array<bool, 6> result = {false, false, false, false, false, false};

    static const Vector3 FCC_SLIP_PLANE_NORMALS[4] = {
        Vector3(1.0, 1.0, 1.0),
        Vector3(1.0, 1.0, -1.0),
        Vector3(1.0, -1.0, 1.0),
        Vector3(-1.0, 1.0, 1.0),
    };

    constexpr double EXPECTED_SQ_MAG = 1.0 / 6.0;
    constexpr double TOLERANCE = 0.001;

    bool isFCC = false;
    for(int f = 0; f < 4; ++f) {
        if(activeFacet[f] && facetBurgers[f].cluster()) {
            isFCC = (facetBurgers[f].cluster()->structure == LATTICE_FCC);
            break;
        }
    }
    if(!isFCC) return result;

    int consensusPlane = 0;

    for(int f = 0; f < 4; ++f) {
        if(!activeFacet[f]) continue;

        const Vector3& b = facetBurgers[f].localVec();
        double sqMag = b.x() * b.x() + b.y() * b.y() + b.z() * b.z();

        if(std::abs(sqMag - EXPECTED_SQ_MAG) > TOLERANCE) continue;

        int plane = 0;
        if(std::abs(b.x() + b.y() + b.z()) < TOLERANCE)
            plane = 1;
        else if(std::abs(b.x() + b.y() - b.z()) < TOLERANCE)
            plane = 2;
        else if(std::abs(b.x() - b.y() + b.z()) < TOLERANCE)
            plane = 3;
        else if(std::abs(-b.x() + b.y() + b.z()) < TOLERANCE)
            plane = 4;

        if(plane == 0) continue;

        if(consensusPlane == 0) {
            consensusPlane = plane;
        } else if(plane != consensusPlane) {
            return result;
        }
    }

    if(consensusPlane == 0) return result;

    const Vector3& normal = FCC_SLIP_PLANE_NORMALS[consensusPlane - 1];
    for(int e = 0; e < 6; ++e) {
        if(!cellEdges[e] || !cellEdges[e].hasEdgeVector()) continue;
        const Vector3& ev = cellEdges[e].vector().localVec();
        double dot = ev.x() * normal.x() + ev.y() * normal.y() + ev.z() * normal.z();
        if(std::abs(dot) < TOLERANCE) {
            result[e] = true;
        }
    }

    return result;
}

void LineReconstructionDXAAlgorithm::dumpDislocationSnapshot(double alpha) {
    appendDislocationSegmentSnapshot();
    appendUnassignedEdgeSnapshot();
    appendInterfaceMeshSnapshot(alpha);
    ++_currentStage;
}

void LineReconstructionDXAAlgorithm::appendDislocationSegmentSnapshot() {
    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isFilledCell(cell)) {
            continue;
        }
        const auto cellEdges = getOrientedEdges(cell);
        for(int facet = 0; facet < 4; ++facet) {
            auto facetEdges = getFacetCircuitEdges(cell, facet);
            if(!facetEdges[0] || !facetEdges[1] || !facetEdges[2]) {
                continue;
            }
            if(facetEdges[0].isBlocked() || facetEdges[1].isBlocked() || facetEdges[2].isBlocked()) {
                continue;
            }
            if(!facetEdges[0].hasEdgeVector() || !facetEdges[1].hasEdgeVector() || !facetEdges[2].hasEdgeVector()) {
                continue;
            }

            ClusterVector burgersVector(Vector3::Zero(), clusterOfAtom(facetEdges[0].atom1()));
            burgersVector.localVec() += facetEdges[0].vector().localVec();
            burgersVector.localVec() += facetEdges[0].transition()->reverseTransform(facetEdges[1].vector().localVec());
            burgersVector.localVec() += facetEdges[2].transition()->transform(facetEdges[2].vector().localVec());
            const Vector3 localBurgers = burgersVector.localVec();
            const bool isActiveFacet = (facet == 0)
                ? (std::abs(localBurgers.x()) > CA_LATTICE_VECTOR_EPSILON || std::abs(localBurgers.y()) > CA_LATTICE_VECTOR_EPSILON)
                : !localBurgers.isZero(CA_LATTICE_VECTOR_EPSILON);
            if(!isActiveFacet) {
                continue;
            }

            if(facetEdges[0].hasEdgeVector() && facetEdges[1].hasEdgeVector() && facetEdges[2].hasEdgeVector()) {
                ClusterTransition* t0 = facetEdges[0].transition();
                ClusterTransition* t1 = facetEdges[1].transition();
                ClusterTransition* t2 = facetEdges[2].transition();
                if(!t0->isSelfTransition() || !t1->isSelfTransition() || !t2->isSelfTransition()) {
                    Matrix3 frankRotation = t2->tm * t1->tm * t0->tm;
                    if(!frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON)) {
                        continue;
                    }
                }
            }

            const auto cellVertices = getCellVertices(_tessellation, cell);
            std::array<Point3, 4> tetPoints = {
                _tessellation.vertexPosition(cellVertices[0]),
                _tessellation.vertexPosition(cellVertices[1]),
                _tessellation.vertexPosition(cellVertices[2]),
                _tessellation.vertexPosition(cellVertices[3])
            };
            _context.simCell.unwrapPositions(tetPoints.data(), tetPoints.size());

            Vector3 facetSum = Vector3::Zero();
            Vector3 cellSum = Vector3::Zero();
            for(int vertex = 0; vertex < 4; ++vertex) {
                cellSum += tetPoints[vertex] - Point3::Origin();
            }
            for(int vertex = 0; vertex < 3; ++vertex) {
                facetSum += tetPoints[DelaunayTessellation::cellFacetVertexIndex(facet, vertex)] - Point3::Origin();
            }
            Point3 faceCenter = Point3::Origin() + facetSum / 3.0;
            Point3 cellCenter = Point3::Origin() + cellSum / 4.0;
            Cluster* cluster = clusterOfAtom(facetEdges[0].atom1());
            ClusterVector normalizedBurgers = normalizeBurgersVector(burgersVector, cluster);
            Cluster* effectiveCluster = normalizedBurgers.cluster();
            int structureType = effectiveCluster ? effectiveCluster->structure : LATTICE_OTHER;
            _dislocationSegments.push_back({
                .position1 = faceCenter,
                .position2 = cellCenter,
                .burgersVector = normalizedBurgers,
                .clusterId = effectiveCluster ? effectiveCluster->id : 0,
                .structureType = structureType,
                .stage = _currentStage,
            });
        }
    }
}

void LineReconstructionDXAAlgorithm::stitchDislocationLines() {
    _dislocationLines.clear();

    struct TmpSegment {
        ClusterVector burgers;
        int node0;
        int node1;
        int clusterId;
        int structureType;
        bool visited;
    };

    std::vector<TmpSegment> tmpSegments;
    std::vector<Point3> nodePositions;
    std::vector<bool> nodeIsCenter;
    int glideFilterCalls = 0;
    int glideFilterHits = 0;
    int glideAdjustedFacets = 0;

    std::unordered_map<std::array<AtomIndex, 3>, int, AtomIndexTripletHash> facetNodeMap;

    auto getOrCreateFacetNode = [&](const std::array<AtomIndex, 3>& key, const Point3& pos) -> std::pair<int, Point3> {
        auto [it, inserted] = facetNodeMap.emplace(key, static_cast<int>(nodePositions.size()));
        if(inserted) {
            nodePositions.push_back(pos);
            nodeIsCenter.push_back(false);
            return {it->second, pos};
        }

        const Point3& stored = nodePositions[it->second];
        Point3 aligned = stored + _context.simCell.unwrapVector(pos - stored);
        return {it->second, aligned};
    };

    auto createCellNode = [&](const Point3& pos) -> int {
        const int nodeId = static_cast<int>(nodePositions.size());
        nodePositions.push_back(pos);
        nodeIsCenter.push_back(true);
        return nodeId;
    };

    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isFilledCell(cell)) {
            continue;
        }
        const auto cellVertices = getCellVertices(_tessellation, cell);
        std::array<Point3, 4> tetPoints = {
            _tessellation.vertexPosition(cellVertices[0]),
            _tessellation.vertexPosition(cellVertices[1]),
            _tessellation.vertexPosition(cellVertices[2]),
            _tessellation.vertexPosition(cellVertices[3])
        };
        _context.simCell.unwrapPositions(tetPoints.data(), tetPoints.size());
        const auto cellEdges = getOrientedEdges(cell);

        bool cellFullyMapped = true;
        for(const auto& edge : cellEdges) {
            if(!edge || edge.isBlocked()) {
                cellFullyMapped = false;
                break;
            }
        }
        if(!cellFullyMapped) {
            continue;
        }

        std::array<AtomIndex, 4> cellAtoms = {
            getInputPointIndex(_tessellation, cellVertices[0]),
            getInputPointIndex(_tessellation, cellVertices[1]),
            getInputPointIndex(_tessellation, cellVertices[2]),
            getInputPointIndex(_tessellation, cellVertices[3])
        };
        Cluster* referenceCluster = nullptr;
        for(AtomIndex atom : cellAtoms) {
            Cluster* cluster = clusterOfAtom(atom);
            if(cluster && cluster->structure != LATTICE_OTHER) {
                referenceCluster = cluster;
                break;
            }
        }
        if(!referenceCluster) {
            continue;
        }

        auto edgeVectorInReferenceCluster = [&](const OrientedEdge& edge) -> Vector3 {
            if(!edge || !edge.hasEdgeVector()) {
                return Vector3::Zero();
            }
            ClusterVector orientedVector = edge.vector();
            Cluster* sourceCluster = orientedVector.cluster();
            if(!sourceCluster) {
                return Vector3::Zero();
            }
            if(sourceCluster == referenceCluster) {
                return orientedVector.localVec();
            }
            ClusterTransition* transform = _structureAnalysis.clusterGraph().determineClusterTransition(sourceCluster, referenceCluster);
            if(!transform) {
                return Vector3::Zero();
            }
            return transform->transform(orientedVector.localVec());
        };

        const Vector3 localFacetVectors[4] = {
            edgeVectorInReferenceCluster(cellEdges[3]) + edgeVectorInReferenceCluster(cellEdges[5]) - edgeVectorInReferenceCluster(cellEdges[4]),
            edgeVectorInReferenceCluster(cellEdges[1]) + edgeVectorInReferenceCluster(cellEdges[5]) - edgeVectorInReferenceCluster(cellEdges[2]),
            edgeVectorInReferenceCluster(cellEdges[0]) + edgeVectorInReferenceCluster(cellEdges[4]) - edgeVectorInReferenceCluster(cellEdges[2]),
            edgeVectorInReferenceCluster(cellEdges[0]) + edgeVectorInReferenceCluster(cellEdges[3]) - edgeVectorInReferenceCluster(cellEdges[1])
        };

        bool activeFacet[4] = {};
        ClusterVector rawFacetBurgers[4];
        ClusterVector facetBurgers[4];
        int facetStructureType[4] = {};
        int facetClusterId[4] = {};
        std::array<AtomIndex, 3> facetKeys[4];

        for(int facet = 0; facet < 4; ++facet) {
            auto facetEdges = getFacetCircuitEdges(cell, facet);
            if(!facetEdges[0] || !facetEdges[1] || !facetEdges[2]) {
                continue;
            }
            if(facetEdges[0].isBlocked() || facetEdges[1].isBlocked() || facetEdges[2].isBlocked()) {
                continue;
            }
            if(!facetEdges[0].hasEdgeVector() || !facetEdges[1].hasEdgeVector() || !facetEdges[2].hasEdgeVector()) {
                continue;
            }

            ClusterVector burgersVector(Vector3::Zero(), clusterOfAtom(facetEdges[0].atom1()));
            if(facetEdges[0].hasEdgeVector()) {
                burgersVector.localVec() += facetEdges[0].vector().localVec();
            }
            if(facetEdges[1].hasEdgeVector()) {
                if(facetEdges[0].hasEdgeVector()) {
                    burgersVector.localVec() += facetEdges[0].transition()->reverseTransform(facetEdges[1].vector().localVec());
                } else if(facetEdges[2].hasEdgeVector()) {
                    burgersVector.localVec() += facetEdges[2].transition()->transform(
                        facetEdges[1].transition()->transform(facetEdges[1].vector().localVec())
                    );
                } else {
                    burgersVector = facetEdges[1].vector();
                }
            }
            if(facetEdges[2].hasEdgeVector()) {
                burgersVector.localVec() += facetEdges[2].transition()->transform(facetEdges[2].vector().localVec());
            }
            const Vector3& localFacetVector = localFacetVectors[facet];
            const bool isActiveFacet = (std::abs(localFacetVector.x()) > CA_LATTICE_VECTOR_EPSILON
                || std::abs(localFacetVector.y()) > CA_LATTICE_VECTOR_EPSILON
                || std::abs(localFacetVector.z()) > CA_LATTICE_VECTOR_EPSILON);
            if(!isActiveFacet) {
                continue;
            }

            ClusterTransition* t0 = facetEdges[0].transition();
            ClusterTransition* t1 = facetEdges[1].transition();
            ClusterTransition* t2 = facetEdges[2].transition();
            if(!t0 || !t1 || !t2) {
                continue;
            }
            if(!t0->isSelfTransition() || !t1->isSelfTransition() || !t2->isSelfTransition()) {
                Matrix3 frankRotation = t2->tm * t1->tm * t0->tm;
                if(!frankRotation.equals(Matrix3::Identity(), 1e-4f)) {
                    continue;
                }
            }

            std::array<AtomIndex, 3> key = {
                static_cast<AtomIndex>(cellVertices[DelaunayTessellation::cellFacetVertexIndex(facet, 0)]),
                static_cast<AtomIndex>(cellVertices[DelaunayTessellation::cellFacetVertexIndex(facet, 1)]),
                static_cast<AtomIndex>(cellVertices[DelaunayTessellation::cellFacetVertexIndex(facet, 2)])
            };
            std::sort(key.begin(), key.end());

            Cluster* rawCluster = clusterOfAtom(facetEdges[0].atom1());
            ClusterVector normalizedBurgers = normalizeBurgersVector(burgersVector, rawCluster);
            Cluster* effectiveCluster = normalizedBurgers.cluster();

            activeFacet[facet] = true;
            rawFacetBurgers[facet] = ClusterVector(localFacetVector, referenceCluster);
            facetBurgers[facet] = normalizedBurgers;
            facetStructureType[facet] = effectiveCluster ? effectiveCluster->structure : LATTICE_OTHER;
            facetClusterId[facet] = effectiveCluster ? effectiveCluster->id : 0;
            facetKeys[facet] = key;
        }

        int activeCount = 0;
        for(int facet = 0; facet < 4; ++facet) {
            if(activeFacet[facet]) activeCount++;
        }

        std::array<bool, 6> glideEdge = {false, false, false, false, false, false};
        int glideEdgeCount = 0;
        if(activeCount >= 2) {
            glideEdge = determineConservativeGlideEdgesFCC(cellEdges, rawFacetBurgers, activeFacet);
            for(int e = 0; e < 6; ++e) if(glideEdge[e]) ++glideEdgeCount;
            ++glideFilterCalls;
            if(glideEdgeCount > 0) ++glideFilterHits;

        }

        int facetNodeIds[4] = {-1, -1, -1, -1};
        for(int facet = 0; facet < 4; ++facet) {
            if(!activeFacet[facet]) continue;

            Point3 fv0 = tetPoints[DelaunayTessellation::cellFacetVertexIndex(facet, 0)];
            Point3 fv1 = tetPoints[DelaunayTessellation::cellFacetVertexIndex(facet, 1)];
            Point3 fv2 = tetPoints[DelaunayTessellation::cellFacetVertexIndex(facet, 2)];

            Point3 faceCenter = Point3::Origin() + ((fv0 - Point3::Origin()) + (fv1 - Point3::Origin()) + (fv2 - Point3::Origin())) / 3.0;

            const int facetLocalVerts[3] = {
                DelaunayTessellation::cellFacetVertexIndex(facet, 0),
                DelaunayTessellation::cellFacetVertexIndex(facet, 1),
                DelaunayTessellation::cellFacetVertexIndex(facet, 2)
            };
            int oppositeVertex = -1;
            int glideCount = 0;
            for(int localOpposite = 0; localOpposite < 3; ++localOpposite) {
                const int localA = facetLocalVerts[(localOpposite + 1) % 3];
                const int localB = facetLocalVerts[(localOpposite + 2) % 3];

                int edgeIndex = -1;
                for(int e = 0; e < 6; ++e) {
                    const int ev0 = kCellEdgeVertices[e][0];
                    const int ev1 = kCellEdgeVertices[e][1];
                    if((ev0 == localA && ev1 == localB) || (ev0 == localB && ev1 == localA)) {
                        edgeIndex = e;
                        break;
                    }
                }
                if(edgeIndex >= 0 && glideEdge[edgeIndex]) {
                    oppositeVertex = localOpposite;
                    glideCount++;
                }
            }

            if(glideCount == 1 && oppositeVertex >= 0) {
                const Point3 weightedOpposite = tetPoints[facetLocalVerts[oppositeVertex]];
                const Point3 other1 = tetPoints[facetLocalVerts[(oppositeVertex + 1) % 3]];
                const Point3 other2 = tetPoints[facetLocalVerts[(oppositeVertex + 2) % 3]];
                faceCenter = Point3::Origin() +
                    ((weightedOpposite - Point3::Origin()) * 2.0 +
                     (other1 - Point3::Origin()) +
                     (other2 - Point3::Origin())) / 4.0;
                ++glideAdjustedFacets;
            }

            facetNodeIds[facet] = getOrCreateFacetNode(facetKeys[facet], faceCenter).first;
        }

        Point3 cellCenter;
        if(activeCount <= 1) {
            Vector3 sum = Vector3::Zero();
            for(int v = 0; v < 4; ++v) {
                sum += tetPoints[v] - Point3::Origin();
            }
            cellCenter = Point3::Origin() + sum / 4.0;
        } else {
            int firstActiveFacet = -1;
            for(int facet = 0; facet < 4; ++facet) {
                if(activeFacet[facet]) {
                    firstActiveFacet = facet;
                    break;
                }
            }
            if(firstActiveFacet < 0) {
                continue;
            }

            const Point3 baseNodePos = nodePositions[facetNodeIds[firstActiveFacet]];
            Vector3 offsetSum = Vector3::Zero();
            int contributingNodes = 1;
            for(int facet = 0; facet < 4; ++facet) {
                if(!activeFacet[facet] || facet == firstActiveFacet) {
                    continue;
                }
                offsetSum += _context.simCell.unwrapVector(nodePositions[facetNodeIds[facet]] - baseNodePos);
                contributingNodes++;
            }
            cellCenter = baseNodePos + offsetSum / static_cast<double>(contributingNodes);
        }

        int centerNode = createCellNode(cellCenter);

        if(activeCount == 0) {
            continue;
        }

        for(int facet = 0; facet < 4; ++facet) {
            if(!activeFacet[facet]) continue;
            tmpSegments.push_back({
                facetBurgers[facet],
                centerNode,
                facetNodeIds[facet],
                facetClusterId[facet],
                facetStructureType[facet],
                false
            });
        }
    }

    if(tmpSegments.empty()) {
        spdlog::info("stitchDislocationLines: no temporary segments to stitch");
        return;
    }

    const int nodeCount = static_cast<int>(nodePositions.size());
    const int segCount = static_cast<int>(tmpSegments.size());
    std::vector<int> nodeDegree(nodeCount, 0);
    for(const auto& segment : tmpSegments) {
        if(segment.node0 >= 0) nodeDegree[segment.node0]++;
        if(segment.node1 >= 0) nodeDegree[segment.node1]++;
    }
    std::vector<int> nodeHead(nodeCount, -1);
    std::vector<int> nodeTail(nodeCount, -1);
    std::vector<int> halfEdgeNext(2 * segCount, -1);

    for(int s = 0; s < segCount; ++s) {
        int he0 = 2 * s;
        int he1 = 2 * s + 1;
        int n0 = tmpSegments[s].node0;
        int n1 = tmpSegments[s].node1;

        if(nodeHead[n1] < 0) {
            nodeHead[n1] = nodeTail[n1] = he0;
        } else {
            halfEdgeNext[nodeTail[n1]] = he0;
            nodeTail[n1] = he0;
        }

        if(nodeHead[n0] < 0) {
            nodeHead[n0] = nodeTail[n0] = he1;
        } else {
            halfEdgeNext[nodeTail[n0]] = he1;
            nodeTail[n0] = he1;
        }
    }

    auto heSegIdx = [](int he) -> int { return he / 2; };
    auto heIsNode1Side = [](int he) -> bool { return (he % 2) == 0; };
    auto heOtherNode = [&](int he) -> int {
        int seg = heSegIdx(he);
        return heIsNode1Side(he) ? tmpSegments[seg].node0 : tmpSegments[seg].node1;
    };

    auto walkLine = [&](int currentNode, std::vector<Point3>& linePoints, bool pushBack) -> bool {
        while(true) {
            int foundHE = -1;
            bool sawAlternative = false;
            const bool enforceUniqueContinuation = nodeIsCenter[currentNode] || nodeDegree[currentNode] > 2;
            for(int he = nodeHead[currentNode]; he >= 0; he = halfEdgeNext[he]) {
                int seg = heSegIdx(he);
                if(tmpSegments[seg].visited) continue;
                if(foundHE >= 0 && enforceUniqueContinuation) {
                    sawAlternative = true;
                    break;
                }
                foundHE = he;
                if(!enforceUniqueContinuation) {
                    break;
                }
            }
            if(foundHE < 0 || sawAlternative) break;

            int seg = heSegIdx(foundHE);
            tmpSegments[seg].visited = true;

            const Point3 currentPos = pushBack ? linePoints.back() : linePoints.front();
            int nextNode = heOtherNode(foundHE);
            Point3 nextPos = currentPos + _context.simCell.unwrapVector(nodePositions[nextNode] - currentPos);

            if(pushBack) {
                linePoints.push_back(nextPos);
            } else {
                linePoints.insert(linePoints.begin(), nextPos);
            }

            currentNode = nextNode;
        }
        return false;
    };

    for(int s = 0; s < segCount; ++s) {
        if(tmpSegments[s].visited) continue;
        tmpSegments[s].visited = true;

        int n0 = tmpSegments[s].node0;
        int n1 = tmpSegments[s].node1;

        std::vector<Point3> linePoints;
        linePoints.push_back(nodePositions[n1]);
        walkLine(n1, linePoints, true);

        Point3 centerPos = linePoints.front() + _context.simCell.unwrapVector(nodePositions[n0] - linePoints.front());
        linePoints.insert(linePoints.begin(), centerPos);
        walkLine(n0, linePoints, false);

        bool isClosed = false;
        if(linePoints.size() >= 3) {
            Vector3 diff = _context.simCell.unwrapVector(linePoints.front() - linePoints.back());
            if(diff.length() < 1e-6) {
                isClosed = true;
                linePoints.pop_back();
            }
        }

        _dislocationLines.push_back({
            std::move(linePoints),
            tmpSegments[s].burgers,
            tmpSegments[s].clusterId,
            tmpSegments[s].structureType,
            isClosed
        });
    }

    int closedCount = 0;
    int openCount = 0;
    double totalLength = 0.0;
    for(const auto& line : _dislocationLines) {
        if(line.isClosed) closedCount++;
        else openCount++;
        for(size_t i = 1; i < line.points.size(); ++i) {
            totalLength += (line.points[i] - line.points[i - 1]).length();
        }
        if(line.isClosed && line.points.size() > 1) {
            totalLength += (line.points.front() - line.points.back()).length();
        }
    }

    spdlog::info("stitchDislocationLines: {} lines ({} closed, {} open), total length {:.2f}, from {} segments with {} nodes",
        _dislocationLines.size(), closedCount, openCount, totalLength, segCount, nodeCount);
    if(glideFilterCalls > 0) {
        spdlog::info("  FCC glide filter: {} cells checked, {} had glide edges, {} facets adjusted",
            glideFilterCalls, glideFilterHits, glideAdjustedFacets);
    }
}

void LineReconstructionDXAAlgorithm::finishDislocationLines(int smoothingLevel, double linePointInterval) {
    for(auto& line : _dislocationLines) {
        line.burgersVector = ClusterVector(-line.burgersVector.localVec(), line.burgersVector.cluster());
    }

    if(linePointInterval > 0.0) {
        constexpr int CORE_SIZE_DEFAULT = 6;

        for(auto& line : _dislocationLines) {
            const int N = static_cast<int>(line.points.size());
            if(N < (line.isClosed ? 3 : 4)) continue;

            std::vector<Point3> coarsened;

            if(!line.isClosed) {
                coarsened.push_back(line.points.front());
            }

            Vector3 sumPos = Vector3::Zero();
            int groupSize = 0;
            int cumulativeWeight = 0;
            int numEmitted = line.isClosed ? 0 : 1;

            int startIdx = line.isClosed ? 0 : 1;
            int endIdx = line.isClosed ? N : N - 1;

            for(int i = startIdx; i < endIdx; ++i) {
                sumPos += (line.points[i] - Point3::Origin());
                groupSize++;
                cumulativeWeight += CORE_SIZE_DEFAULT;

                double threshold = static_cast<double>(cumulativeWeight) * linePointInterval;
                bool shouldEmit = (static_cast<double>(groupSize * groupSize) >= threshold);

                int remaining = endIdx - i;
                int maxGroupSize = (numEmitted > 0) ? (N / (2 * numEmitted)) : N;
                if(groupSize >= maxGroupSize && remaining > 0) {
                    shouldEmit = true;
                }

                if(shouldEmit) {
                    Vector3 centroid = sumPos / static_cast<double>(groupSize);
                    coarsened.push_back(Point3::Origin() + centroid);
                    numEmitted++;
                    sumPos = Vector3::Zero();
                    groupSize = 0;
                    cumulativeWeight = 0;
                }
            }

            if(groupSize > 0) {
                Vector3 centroid = sumPos / static_cast<double>(groupSize);
                coarsened.push_back(Point3::Origin() + centroid);
            }

            if(!line.isClosed) {
                coarsened.push_back(line.points.back());
            }

            if(line.isClosed && coarsened.size() >= 3) {
                double closingDist = (coarsened.back() - coarsened.front()).length();
                if(closingDist < linePointInterval * 0.5) {
                    coarsened.pop_back();
                }
            }

            if(coarsened.size() >= 2) {
                line.points = std::move(coarsened);
            }
        }
    }

    if(smoothingLevel > 0) {
        constexpr double lambda = 0.5;
        constexpr double mu = -10.0 / 19.0;

        for(int iter = 0; iter < smoothingLevel; ++iter) {
            for(int pass = 0; pass < 2; ++pass) {
                double weight = (pass == 0) ? lambda : mu;

                for(auto& line : _dislocationLines) {
                    const int n = static_cast<int>(line.points.size());
                    if(n < 3) continue;

                    std::vector<Vector3> delta(n, Vector3::Zero());

                    if(line.isClosed) {
                        for(int i = 0; i < n; ++i) {
                            int prev = (i - 1 + n) % n;
                            int next = (i + 1) % n;
                            delta[i] = ((line.points[prev] - Point3::Origin()) +
                                       (line.points[next] - Point3::Origin())) * 0.5
                                       - (line.points[i] - Point3::Origin());
                        }
                        for(int i = 0; i < n; ++i) {
                            line.points[i] = Point3::Origin() +
                                (line.points[i] - Point3::Origin()) + delta[i] * weight;
                        }
                    } else {
                        for(int i = 1; i < n - 1; ++i) {
                            delta[i] = ((line.points[i - 1] - Point3::Origin()) +
                                       (line.points[i + 1] - Point3::Origin())) * 0.5
                                       - (line.points[i] - Point3::Origin());
                        }
                        for(int i = 1; i < n - 1; ++i) {
                            line.points[i] = Point3::Origin() +
                                (line.points[i] - Point3::Origin()) + delta[i] * weight;
                        }
                    }
                }
            }
        }
    }

    classifyDislocationTypes();

    double finalLength = 0.0;
    int closedCount = 0;
    for(const auto& line : _dislocationLines) {
        double lineLength = 0.0;
        if(line.isClosed) {
            closedCount++;
        }
        for(size_t i = 1; i < line.points.size(); ++i) {
            lineLength += (line.points[i] - line.points[i - 1]).length();
        }
        if(line.isClosed && line.points.size() > 1) {
            lineLength += (line.points.front() - line.points.back()).length();
        }
        finalLength += lineLength;
    }
    spdlog::info(
        "finishDislocationLines: {} final lines ({} closed, {} open), total length {:.2f}",
        _dislocationLines.size(),
        closedCount,
        static_cast<int>(_dislocationLines.size()) - closedCount,
        finalLength
    );
}

void LineReconstructionDXAAlgorithm::classifyDislocationTypes() {
    constexpr double TOL = 0.01;

    int perfect = 0, shockley = 0, frank = 0, stairrod = 0, hirth = 0, other = 0;

    for(auto& line : _dislocationLines) {
        if(line.structureType != LATTICE_FCC) {
            other++;
            continue;
        }
        const Vector3& b = line.burgersVector.localVec();
        double sqMag = b.x() * b.x() + b.y() * b.y() + b.z() * b.z();

        if(std::abs(sqMag - 0.5) < TOL) {
            line.dislocationTypeId = 1;
            perfect++;
        } else if(std::abs(sqMag - 1.0 / 6.0) < TOL) {
            line.dislocationTypeId = 2;
            shockley++;
        } else if(std::abs(sqMag - 1.0 / 3.0) < TOL) {
            line.dislocationTypeId = 3;
            frank++;
        } else if(std::abs(sqMag - 1.0 / 18.0) < TOL) {
            line.dislocationTypeId = 4;
            stairrod++;
        } else if(std::abs(sqMag - 1.0 / 9.0) < TOL) {
            line.dislocationTypeId = 5;
            hirth++;
        } else {
            line.dislocationTypeId = 0;
            other++;
        }
    }

    spdlog::info("FCC dislocation types: {} perfect, {} Shockley, {} Frank, {} stair-rod, {} Hirth, {} other",
        perfect, shockley, frank, stairrod, hirth, other);
}

void LineReconstructionDXAAlgorithm::appendUnassignedEdgeSnapshot() {
    std::unordered_set<TessellationEdge*> uniqueEdges;
    uniqueEdges.reserve(_unassignedEdges.size());

    for(TessellationEdge* edge : _unassignedEdges) {
        if(!uniqueEdges.insert(edge).second) {
            continue;
        }
        if(!edge->transition || edge->hasEdgeVector) {
            continue;
        }
        _unassignedEdgesOutput.push_back({
            .position1 = _context.positions->getPoint3(edge->atom1),
            .position2 = _context.positions->getPoint3(edge->atom2),
            .atom1 = edge->atom1,
            .atom2 = edge->atom2,
            .stage = _currentStage,
        });
    }
}

void LineReconstructionDXAAlgorithm::appendInterfaceMeshSnapshot(double alpha) {
    for(DelaunayTessellation::CellHandle cell : _tessellation.cells()) {
        if(!isInteriorInterfaceCell(cell, alpha)) {
            continue;
        }
        for(int facet = 0; facet < 4; ++facet) {
            DelaunayTessellation::Facet mirrorFacet = _tessellation.mirrorFacet(cell, facet);
            if(isInteriorInterfaceCell(mirrorFacet.first, alpha)) {
                continue;
            }

            std::array<LineReconstructionInterfaceMesh::Vertex*, 3> facetVertices{};
            for(int vertex = 0; vertex < 3; ++vertex) {
                const auto handle = _tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(facet, 2 - vertex));
                facetVertices[vertex] = _interfaceMesh.createVertex(_tessellation.vertexPosition(handle));
            }

            auto* face = _interfaceMesh.createFace(facetVertices.begin(), facetVertices.end());
            face->region = _currentStage;
        }
    }
}

bool LineReconstructionDXAAlgorithm::isInteriorInterfaceCell(DelaunayTessellation::CellHandle cell, double alpha) const {
    if(!isFiniteCell(_tessellation, cell)) {
        return false;
    }

    bool allCrystalline = true;
    for(int vertex = 0; vertex < 4; ++vertex) {
        const AtomIndex atomIndex = getInputPointIndex(_tessellation, _tessellation.cellVertex(cell, vertex));
        Cluster* cluster = _structureAnalysis.atomCluster(atomIndex);
        if(!cluster || cluster->id == 0) {
            allCrystalline = false;
            break;
        }
    }
    if(allCrystalline) {
        return true;
    }

    if(auto alphaTestResult = _tessellation.alphaTest(cell, alpha)) {
        return *alphaTestResult;
    }

    for(int facet = 0; facet < 4; ++facet) {
        DelaunayTessellation::CellHandle adjacentCell = _tessellation.mirrorFacet(cell, facet).first;
        if(!isFiniteCell(_tessellation, adjacentCell)) {
            return false;
        }
        auto adjacentAlphaTestResult = _tessellation.alphaTest(adjacentCell, alpha);
        if(adjacentAlphaTestResult.has_value() && !adjacentAlphaTestResult.value()) {
            return false;
        }
    }
    return true;
}

ClusterVector LineReconstructionDXAAlgorithm::normalizeBurgersVector(const ClusterVector& burgersVector, Cluster* cluster) const {
    ClusterVector normalized = burgersVector;
    if(!cluster) {
        return normalized;
    }
    if(cluster->structure == _context.inputCrystalType) {
        return normalized;
    }
    for(ClusterTransition* transition = cluster->transitions; transition != nullptr && transition->distance <= 1; transition = transition->next) {
        if(transition->cluster2 && transition->cluster2->structure == _context.inputCrystalType) {
            return ClusterVector(transition->transform(burgersVector.localVec()), transition->cluster2);
        }
    }
    return normalized;
}

bool LineReconstructionDXAAlgorithm::isFilledCell(DelaunayTessellation::CellHandle cell) const {
    return cell < _filledCells.size() ? _filledCells[cell] : false;
}

bool LineReconstructionDXAAlgorithm::isElasticMappingCompatible(const std::array<OrientedEdge, 6>& cellEdges) const {
    for(int facet = 0; facet < 4; ++facet) {
        auto facetEdges = getFacetCircuitEdges(cellEdges, facet);
        if(!facetEdges[0] || !facetEdges[1] || !facetEdges[2]) {
            return false;
        }
        ClusterTransition* t0 = facetEdges[0].transition();
        ClusterTransition* t1 = facetEdges[1].transition();
        ClusterTransition* t2 = facetEdges[2].transition();
        if(!t0 || !t1 || !t2) {
            return false;
        }
        Vector3 burgersVector = facetEdges[0].vector().localVec()
            + t0->reverseTransform(facetEdges[1].vector().localVec())
            + t2->transform(facetEdges[2].vector().localVec());
        if(!burgersVector.isZero(CA_LATTICE_VECTOR_EPSILON)) {
            return false;
        }
        if(!t0->isSelfTransition() || !t1->isSelfTransition() || !t2->isSelfTransition()) {
            Matrix3 frankRotation = t2->tm * t1->tm * t0->tm;
            if(!frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON)) {
                return false;
            }
        }
    }
    return true;
}

std::array<LineReconstructionDXAAlgorithm::OrientedEdge, 6> LineReconstructionDXAAlgorithm::getOrientedEdges(DelaunayTessellation::CellHandle cell) const {
    const auto cellVertices = getCellVertices(_tessellation, cell);
    std::array<OrientedEdge, 6> edges;
    for(int index = 0; index < 6; ++index) {
        edges[index] = getOrientedEdge(cellVertices, index);
    }
    return edges;
}

LineReconstructionDXAAlgorithm::OrientedEdge LineReconstructionDXAAlgorithm::getOrientedEdge(AtomIndex atomIndex1, AtomIndex atomIndex2) const {
    TessellationEdge* edge = findEdge(atomIndex1, atomIndex2);
    if(!edge) {
        return {};
    }
    return OrientedEdge(edge, edge->atom1 != atomIndex1);
}

LineReconstructionDXAAlgorithm::OrientedEdge LineReconstructionDXAAlgorithm::getOrientedEdge(const std::array<DelaunayTessellation::VertexHandle, 4>& cellVertices, int localEdgeIndex) const {
    AtomIndex atom1 = getInputPointIndex(_tessellation, cellVertices[kCellEdgeVertices[localEdgeIndex][0]]);
    AtomIndex atom2 = getInputPointIndex(_tessellation, cellVertices[kCellEdgeVertices[localEdgeIndex][1]]);
    return getOrientedEdge(atom1, atom2);
}

std::array<LineReconstructionDXAAlgorithm::OrientedEdge, 3> LineReconstructionDXAAlgorithm::getFacetCircuitEdges(DelaunayTessellation::CellHandle cell, int facetIndex) const {
    const auto cellVertices = getCellVertices(_tessellation, cell);
    switch(facetIndex) {
        case 0: return { getOrientedEdge(cellVertices, 3), getOrientedEdge(cellVertices, 5), -getOrientedEdge(cellVertices, 4) };
        case 1: return { getOrientedEdge(cellVertices, 2), -getOrientedEdge(cellVertices, 5), -getOrientedEdge(cellVertices, 1) };
        case 2: return { getOrientedEdge(cellVertices, 0), getOrientedEdge(cellVertices, 4), -getOrientedEdge(cellVertices, 2) };
        case 3: return { getOrientedEdge(cellVertices, 1), -getOrientedEdge(cellVertices, 3), -getOrientedEdge(cellVertices, 0) };
        default: return {};
    }
}

std::array<LineReconstructionDXAAlgorithm::OrientedEdge, 3> LineReconstructionDXAAlgorithm::getFacetCircuitEdges(const std::array<OrientedEdge, 6>& cellEdges, int facetIndex) {
    switch(facetIndex) {
        case 0: return { cellEdges[3], cellEdges[5], -cellEdges[4] };
        case 1: return { cellEdges[2], -cellEdges[5], -cellEdges[1] };
        case 2: return { cellEdges[0], cellEdges[4], -cellEdges[2] };
        case 3: return { cellEdges[1], -cellEdges[3], -cellEdges[0] };
        default: return {};
    }
}

LineReconstructionDXAAlgorithm::TessellationEdge* LineReconstructionDXAAlgorithm::findEdge(AtomIndex atomIndex1, AtomIndex atomIndex2) const {
    for(TessellationEdge* edge = _atomOutboundEdges[atomIndex1]; edge != nullptr; edge = edge->nextOutboundEdge) {
        if(edge->atom2 == atomIndex2) {
            return edge;
        }
    }
    for(TessellationEdge* edge = _atomInboundEdges[atomIndex1]; edge != nullptr; edge = edge->nextInboundEdge) {
        if(edge->atom1 == atomIndex2) {
            return edge;
        }
    }
    return nullptr;
}

Cluster* LineReconstructionDXAAlgorithm::clusterOfAtom(AtomIndex atomIndex) const {
    return atomIndex < static_cast<AtomIndex>(_atomClusters.size()) ? _atomClusters[atomIndex] : nullptr;
}

void LineReconstructionDXAAlgorithm::reorderFacetVertices(std::array<AtomIndex, 3>& atomIndices) {
    std::sort(atomIndices.begin(), atomIndices.end());
}

bool LineReconstructionDXAAlgorithm::isVertexFullyIsolated(AtomIndex atom) const {
    for(TessellationEdge* edge = _atomOutboundEdges[atom]; edge != nullptr; edge = edge->nextOutboundEdge) {
        if(edge->hasEdgeVector) {
            return false;
        }
    }
    for(TessellationEdge* edge = _atomInboundEdges[atom]; edge != nullptr; edge = edge->nextInboundEdge) {
        if(edge->hasEdgeVector) {
            return false;
        }
    }
    return true;
}

}
