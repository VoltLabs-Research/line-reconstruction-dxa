#include <volt/lrdxa_log.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/lrdxa_engine.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <spdlog/spdlog.h>

#include <cstdint>
#include <string>

namespace Volt::DXA {

void LineReconstructionDXALog::addLength(const LineReconstructionDislocationLine& line, double& length) {
    for(size_t i = 1; i < line.points.size(); ++i) {
        length += (line.points[i] - line.points[i - 1]).length();
    }
    if(line.isClosed && line.points.size() > 1) {
        length += (line.points.front() - line.points.back()).length();
    }
}

void LineReconstructionDXALog::dump(const StructureAnalysis& structure) {
    // clusters() includes the graph's built-in null cluster (id 0).
    spdlog::info("Reconstructed structure loaded in {}: {} clusters, {} cluster transitions, max neighbor distance {:.3f}",
        phaseElapsed(),
        structure.clusterGraph().clusters().size() - 1,
        structure.clusterGraph().clusterTransitions().size(),
        structure.maximumNeighborDistance());
    _reconstruction.restart();
}

void LineReconstructionDXALog::dump(const TessellationVertices& vertices) {
    spdlog::info("Selected {} of {} atoms ({:.1f}%) as tessellation vertices ({}) in {}",
        vertices.selected, vertices.atoms,
        percentOf(static_cast<std::uint64_t>(vertices.selected), vertices.atoms),
        vertices.wholeCrystal ? "whole crystal" : "defect regions",
        phaseElapsed());
}

void LineReconstructionDXALog::dump(const DelaunayTessellation& tessellation) {
    spdlog::info("Delaunay tessellation ready in {}: {} primary tetrahedra",
        phaseElapsed(), tessellation.numberOfPrimaryTetrahedra());
}

void LineReconstructionDXALog::dump(const ClassifiedTetrahedra& tetrahedra) {
    spdlog::info("Classified {} tetrahedra by alpha-shape test (alpha {:.3f}) in {}",
        tetrahedra.tetrahedra, tetrahedra.alpha, phaseElapsed());
}

void LineReconstructionDXALog::dump(const TessellationEdges& edges) {
    spdlog::info("Generated {} tessellation edges from {} filled tetrahedra ({:.1f}% of primary) in {}",
        edges.edges, edges.filledTetrahedra, percentOf(edges.filledTetrahedra, edges.primaryTetrahedra), phaseElapsed());
}

void LineReconstructionDXALog::dump(const ClusterAssignment& assignment) {
    spdlog::info("Assigned {} of {} non-crystalline atoms to a neighboring crystal cluster in {} ({} passes, {} left unassigned)",
        assignment.nonCrystallineAtoms - assignment.unassignedAtoms, assignment.nonCrystallineAtoms, phaseElapsed(),
        assignment.passes, assignment.unassignedAtoms);
}

void LineReconstructionDXALog::dump(const IdealVectorAssignment& assignment) {
    spdlog::info("Assigned ideal vectors to {} of {} edges in {}: {} unresolved, {} skipped (non-crystalline atom or no cluster transition)",
        assignment.assigned, assignment.edges, phaseElapsed(), assignment.unresolved,
        assignment.edges - assignment.assigned - assignment.unresolved);
}

void LineReconstructionDXALog::dump(const FacetIndex& index) {
    spdlog::info("Indexed {} facets with missing edge vectors in {}", index.facets, phaseElapsed());
}

void LineReconstructionDXALog::dump(const MappedNetwork& network) {
    const char* label = "initial mapping";
    if(network.pass == MappedNetwork::Pass::AfterLoopCleanup) label = "after loop cleanup";
    if(network.pass == MappedNetwork::Pass::Final) label = "final mapping";
    const NetworkSnapshot& snapshot = network.snapshot;
    spdlog::info("Stage {} ({}): {} dislocation segments, {} unresolved edges, {} interface facets in {}",
        snapshot.stage, label, snapshot.segments, snapshot.unresolvedEdges, snapshot.interfaceFacets, phaseElapsed());
    _sinceSnapshot.restart();
}

void LineReconstructionDXALog::dump(const SpuriousLoopElimination& elimination) {
    spdlog::info("Spurious loop elimination: invalidated {} edges ({} of {} fully mapped tetrahedra incompatible) in {}",
        elimination.invalidatedEdges, elimination.incompatibleTetrahedra, elimination.mappedTetrahedra, phaseElapsed());
}

void LineReconstructionDXALog::dump(const SkippedEdgeVectorCompletion& skipped) {
    spdlog::info("Skipping edge vector completion: {} crystal without cluster transitions",
        latticeStructureTypeName(skipped.crystalType));
    phaseMilliseconds();
}

void LineReconstructionDXALog::dump(const EdgeVectorCompletionPass& pass) {
    // Read for every pass, logged or not, so that each line times its own pass.
    const std::string elapsed = phaseElapsed();
    if(pass.iteration <= 10 || pass.iteration % 10 == 0) {
        spdlog::info("  Iteration {}: {} unresolved edges after the forward pass, {} after the backward pass ({})",
            pass.iteration, pass.unresolvedAfterForward, pass.unresolvedAfterBackward, elapsed);
    }
}

void LineReconstructionDXALog::dump(const EdgeVectorCompletion& completion) {
    spdlog::info("Edge vector completion converged after {} iteration{} in {}: {} edges recovered, {} still unresolved",
        completion.iterations, completion.iterations == 1 ? "" : "s", _sinceSnapshot.elapsed(),
        completion.unresolvedBefore - completion.unresolvedAfter, completion.unresolvedAfter);
    phaseMilliseconds();
}

void LineReconstructionDXALog::dump(const StitchedLines& stitched) {
    if(stitched.segments == 0) {
        spdlog::info("No dislocation segments to stitch ({})", phaseElapsed());
        return;
    }
    int closedCount = 0;
    double totalLength = 0.0;
    for(const auto& line : stitched.lines) {
        if(line.isClosed) closedCount++;
        addLength(line, totalLength);
    }
    spdlog::info("Stitched {} dislocation lines ({} closed, {} open), total length {:.2f}, from {} segments and {} nodes in {}",
        stitched.lines.size(), closedCount, static_cast<int>(stitched.lines.size()) - closedCount, totalLength,
        stitched.segments, stitched.nodes, phaseElapsed());
    if(stitched.glideFilterCalls > 0) {
        spdlog::info("  FCC glide filter: {} cells checked, {} had glide edges, {} facets adjusted",
            stitched.glideFilterCalls, stitched.glideFilterHits, stitched.glideAdjustedFacets);
    }
}

void LineReconstructionDXALog::dump(const std::vector<LineReconstructionDislocationLine>& finishedLines) {
    // Indexed by dislocation type id: other, perfect, Shockley, Frank, stair-rod, Hirth.
    int types[6] = {};
    int fccLines = 0;
    for(const auto& line : finishedLines) {
        if(line.structureType != LATTICE_FCC) {
            types[0]++;
            continue;
        }
        fccLines++;
        types[line.dislocationTypeId]++;
    }
    // Lines of non-FCC crystals all land in "other", so the breakdown says nothing without FCC lines.
    if(fccLines > 0) {
        spdlog::info("FCC dislocation types: {} perfect, {} Shockley, {} Frank, {} stair-rod, {} Hirth, {} other",
            types[1], types[2], types[3], types[4], types[5], types[0]);
    }

    double finalLength = 0.0;
    int closedCount = 0;
    for(const auto& line : finishedLines) {
        if(line.isClosed) closedCount++;
        double lineLength = 0.0;
        addLength(line, lineLength);
        finalLength += lineLength;
    }
    spdlog::info("Finished {} dislocation lines ({} closed, {} open), final total length {:.2f}, in {}",
        finishedLines.size(), closedCount, static_cast<int>(finishedLines.size()) - closedCount, finalLength,
        phaseElapsed());
}

void LineReconstructionDXALog::dump(const LineReconstructionDXAAlgorithm& reconstruction) {
    spdlog::info(
        "Line reconstruction finished in {}: {} dislocation lines, {} snapshot segments and {} unresolved edge records over {} stages",
        _reconstruction.elapsed(),
        reconstruction.dislocationLineCount(),
        reconstruction.dislocationSegments().size(),
        reconstruction.unassignedEdges().size(),
        reconstruction.stageCount()
    );
    phaseMilliseconds();
}

void LineReconstructionDXALog::dump(const ExportedTable& table) {
    switch(table.contents) {
        case ExportedTable::Contents::DislocationLines:
            spdlog::info("Exported {} dislocation lines to: {} ({})", table.rows, table.path, phaseElapsed());
            break;
        case ExportedTable::Contents::SnapshotSegments:
            spdlog::info("Exported {} snapshot segments to: {} ({})", table.rows, table.path, phaseElapsed());
            break;
        case ExportedTable::Contents::UnresolvedEdges:
            spdlog::info("Exported {} unresolved edge records to: {} ({})", table.rows, table.path, phaseElapsed());
            break;
        case ExportedTable::Contents::InterfaceMesh:
            spdlog::info("Exported interface mesh to: {} ({})", table.path, phaseElapsed());
            break;
        case ExportedTable::Contents::SimulationCell:
            spdlog::info("Exported simulation cell to: {} ({})", table.path, phaseElapsed());
            break;
    }
}

}
