#pragma once

#include <volt/core/volt.h>
#include <volt/helpers/half_edge_mesh.h>
#include <volt/helpers/cluster_vector.h>
#include <volt/structures/crystal_structure_types.h>

#include <cstddef>
#include <string>
#include <vector>

namespace Volt::DXA {

struct LineReconstructionSegment {
    Point3 position1;
    Point3 position2;
    ClusterVector burgersVector;
    int clusterId = 0;
    int structureType = LATTICE_OTHER;
    int stage = 0;
};

struct LineReconstructionUnassignedEdge {
    Point3 position1;
    Point3 position2;
    int atom1 = -1;
    int atom2 = -1;
    int stage = 0;
};

struct LineReconstructionDislocationLine {
    std::vector<Point3> points;
    ClusterVector burgersVector;
    int clusterId = 0;
    int structureType = LATTICE_OTHER;
    bool isClosed = false;
    int dislocationTypeId = 0;
};

struct LineReconstructionInterfaceMeshVertex {};

struct LineReconstructionInterfaceMeshFace {
    int region = 0;
};

struct LineReconstructionInterfaceMeshEdge {};

using LineReconstructionInterfaceMesh = HalfEdgeMesh<LineReconstructionInterfaceMeshEdge, LineReconstructionInterfaceMeshFace, LineReconstructionInterfaceMeshVertex>;

// The atoms the tessellation is built on: those of the defect regions, or every atom of a crystal tessellated whole.
struct TessellationVertices {
    int selected = 0;
    std::size_t atoms = 0;
    bool wholeCrystal = false;
};

// The tetrahedra the alpha-shape test sorted into filled and empty ones.
struct ClassifiedTetrahedra {
    std::size_t tetrahedra = 0;
    double alpha = 0.0;
};

struct TessellationEdges {
    std::size_t edges = 0;
    std::size_t filledTetrahedra = 0;
    std::size_t primaryTetrahedra = 0;
};

// The non-crystalline atoms given the cluster of a crystalline neighbor, pass after pass until a pass assigns none.
struct ClusterAssignment {
    std::size_t nonCrystallineAtoms = 0;
    int passes = 0;
    std::size_t unassignedAtoms = 0;
};

// The edges given the ideal lattice vector of a crystal path; the unresolved ones found no path.
struct IdealVectorAssignment {
    std::size_t edges = 0;
    std::size_t assigned = 0;
    std::size_t unresolved = 0;
};

// The facets with a missing edge vector, indexed by their atoms.
struct FacetIndex {
    std::size_t facets = 0;
};

struct SpuriousLoopElimination {
    int invalidatedEdges = 0;
    std::size_t incompatibleTetrahedra = 0;
    std::size_t mappedTetrahedra = 0;
};

// What one snapshot of the dislocation network added. Each snapshot appends its stage's complete set of segments,
// unresolved edges and interface facets, so these are that stage's counts.
struct NetworkSnapshot {
    int stage = 0;
    std::size_t segments = 0;
    std::size_t unresolvedEdges = 0;
    std::size_t interfaceFacets = 0;
};

// The snapshot taken once one of the mapping passes is over.
struct MappedNetwork {
    enum class Pass { Initial, AfterLoopCleanup, Final };

    Pass pass;
    NetworkSnapshot snapshot;
};

// A crystal tessellated whole whose cluster graph has no transitions leaves no edge vector to complete.
struct SkippedEdgeVectorCompletion {
    int crystalType = LATTICE_OTHER;
};

// One forward and one backward completion pass, and the unresolved edges each of them left.
struct EdgeVectorCompletionPass {
    int iteration = 0;
    std::size_t unresolvedAfterForward = 0;
    std::size_t unresolvedAfterBackward = 0;
    bool changed = false;
};

struct EdgeVectorCompletion {
    int iterations = 0;
    std::size_t unresolvedBefore = 0;
    std::size_t unresolvedAfter = 0;
};

// The lines stitched from the segments and nodes of the fully mapped tetrahedra, and what the FCC glide filter did on
// the way.
struct StitchedLines {
    const std::vector<LineReconstructionDislocationLine>& lines;
    std::size_t segments = 0;
    std::size_t nodes = 0;
    int glideFilterCalls = 0;
    int glideFilterHits = 0;
    int glideAdjustedFacets = 0;
};

// A table the reconstruction wrote, and how many rows of its contents it holds when it lists them.
struct ExportedTable {
    enum class Contents { DislocationLines, SnapshotSegments, UnresolvedEdges, InterfaceMesh, SimulationCell };

    Contents contents;
    std::string path;
    std::size_t rows = 0;
};

}
