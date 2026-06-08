#pragma once

#include <volt/core/volt.h>
#include <volt/core/particle_property.h>
#include <volt/lrdxa_options.h>
#include <volt/lrdxa_types.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/helpers/cluster_vector.h>
#include <volt/utilities/memory_pool.h>
#include <functional>
#include <string_view>

namespace Volt {

class StructureContext;
class StructureAnalysis;

namespace DXA {

class LineReconstructionDXAAlgorithm {
public:
    using AtomIndex = int;

    explicit LineReconstructionDXAAlgorithm(StructureAnalysis& structureAnalysis, StructureContext& context);

    void run(const LineReconstructionDXAOptions& options, const std::function<void(std::string_view)>& markStage = {});

    const std::vector<LineReconstructionSegment>& dislocationSegments() const noexcept { return _dislocationSegments; }
    const std::vector<LineReconstructionDislocationLine>& dislocationLines() const noexcept { return _dislocationLines; }
    size_t dislocationLineCount() const noexcept { return _dislocationLines.size(); }
    const std::vector<LineReconstructionUnassignedEdge>& unassignedEdges() const noexcept { return _unassignedEdgesOutput; }
    const LineReconstructionInterfaceMesh& interfaceMesh() const noexcept { return _interfaceMesh; }
    int delaunayVertexCount() const noexcept { return _delaunayVertexCount; }

private:
    struct TessellationEdge {
        AtomIndex atom1;
        AtomIndex atom2;
        ClusterVector vector;
        ClusterTransition* transition = nullptr;
        DelaunayTessellation::CellHandle adjacentCell{};
        TessellationEdge* nextOutboundEdge = nullptr;
        TessellationEdge* nextInboundEdge = nullptr;
        bool hasEdgeVector = false;

        TessellationEdge(AtomIndex atom1, AtomIndex atom2, DelaunayTessellation::CellHandle adjacentCell) noexcept
            : atom1(atom1), atom2(atom2), adjacentCell(adjacentCell) {}
    };

    class OrientedEdge {
    public:
        OrientedEdge() = default;
        OrientedEdge(TessellationEdge* edge, bool flipped) : _edge(edge), _flipped(flipped) {}

        explicit operator bool() const noexcept { return _edge != nullptr; }
        AtomIndex atom1() const { return _flipped ? _edge->atom2 : _edge->atom1; }
        AtomIndex atom2() const { return _flipped ? _edge->atom1 : _edge->atom2; }
        bool hasEdgeVector() const { return _edge && _edge->hasEdgeVector; }
        bool isBlocked() const { return !_edge || _edge->transition == nullptr; }
        ClusterTransition* transition() const { return !_flipped ? _edge->transition : (_edge->transition ? _edge->transition->reverse : nullptr); }
        TessellationEdge* undirectedEdge() const { return _edge; }
        ClusterVector vector() const;
        OrientedEdge operator-() const { return OrientedEdge(_edge, !_flipped); }

    private:
        TessellationEdge* _edge = nullptr;
        bool _flipped = false;
    };

    struct AtomIndexTripletHash {
        size_t operator()(const std::array<AtomIndex, 3>& arr) const noexcept;
    };

    static constexpr int kCellEdgeVertices[6][2] = {
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
    };

    void markPerfectCrystallineRegions();
    void buildTessellation(double ghostLayerSize);
    void classifyTetrahedra(double alpha);
    void generateTessellationEdges();
    void assignAtomsToClusters();
    void assignIdealVectorsToEdges(int crystalPathSteps);
    void eliminateSpuriousDislocationLoops();
    void buildFacetLookupMap();
    bool complementEdgeVectors(bool forward);
    void dumpDislocationSnapshot(double alpha);
    void appendDislocationSegmentSnapshot();
    void stitchDislocationLines();
    void finishDislocationLines(int smoothingLevel, double linePointInterval);
    void classifyDislocationTypes();
    void appendUnassignedEdgeSnapshot();
    void appendInterfaceMeshSnapshot(double alpha);
    bool isInteriorInterfaceCell(DelaunayTessellation::CellHandle cell, double alpha) const;
    ClusterVector normalizeBurgersVector(const ClusterVector& burgersVector, Cluster* cluster) const;

    static std::array<bool, 6> determineConservativeGlideEdgesFCC(
        const std::array<OrientedEdge, 6>& cellEdges,
        const ClusterVector facetBurgers[4],
        const bool activeFacet[4]);

    bool isFilledCell(DelaunayTessellation::CellHandle cell) const;
    bool isElasticMappingCompatible(const std::array<OrientedEdge, 6>& cellEdges) const;
    std::array<OrientedEdge, 6> getOrientedEdges(DelaunayTessellation::CellHandle cell) const;
    OrientedEdge getOrientedEdge(AtomIndex atomIndex1, AtomIndex atomIndex2) const;
    OrientedEdge getOrientedEdge(const std::array<DelaunayTessellation::VertexHandle, 4>& cellVertices, int localEdgeIndex) const;
    std::array<OrientedEdge, 3> getFacetCircuitEdges(DelaunayTessellation::CellHandle cell, int facetIndex) const;
    static std::array<OrientedEdge, 3> getFacetCircuitEdges(const std::array<OrientedEdge, 6>& cellEdges, int facetIndex);
    TessellationEdge* findEdge(AtomIndex atomIndex1, AtomIndex atomIndex2) const;
    Cluster* clusterOfAtom(AtomIndex atomIndex) const;
    static void reorderFacetVertices(std::array<AtomIndex, 3>& atomIndices);
    bool isVertexFullyIsolated(AtomIndex atom) const;

private:
    StructureAnalysis& _structureAnalysis;
    StructureContext& _context;
    DelaunayTessellation _tessellation;
    std::vector<Cluster*> _atomClusters;
    std::vector<TessellationEdge*> _atomOutboundEdges;
    std::vector<TessellationEdge*> _atomInboundEdges;
    std::vector<TessellationEdge*> _unassignedEdges;
    // vector<char> (not vector<bool>): classifyTetrahedra writes distinct
    // elements from multiple threads, which the bit-packed vector<bool> cannot
    // do safely.
    std::vector<char> _filledCells;
    std::unordered_map<std::array<AtomIndex, 3>, DelaunayTessellation::Facet, AtomIndexTripletHash> _primaryFacetLookupMap;
    MemoryPool<TessellationEdge> _edgePool;
    std::shared_ptr<ParticleProperty> _delaunayAtomMask;
    std::vector<LineReconstructionSegment> _dislocationSegments;
    std::vector<LineReconstructionDislocationLine> _dislocationLines;
    std::vector<LineReconstructionUnassignedEdge> _unassignedEdgesOutput;
    LineReconstructionInterfaceMesh _interfaceMesh;
    int _currentStage = 0;
    int _delaunayVertexCount = 0;
};

}

}

namespace Volt {
using LineReconstructionDXAEngine = DXA::LineReconstructionDXAAlgorithm;
}
