#pragma once

#include <volt/core/run_log.h>
#include <volt/core/progress.h>
#include <volt/lrdxa_types.h>

#include <vector>

namespace Volt {

class DelaunayTessellation;
class StructureAnalysis;

namespace DXA {

class LineReconstructionDXAAlgorithm;

// The execution log of one line reconstruction. Each phase hands it what it built, and the log writes what that holds
// and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class LineReconstructionDXALog : public RunLog {
public:
    using RunLog::dump;

    // The structure of the upstream contract; the reconstruction starts once it is loaded.
    void dump(const StructureAnalysis& structure);
    void dump(const TessellationVertices& vertices);
    void dump(const DelaunayTessellation& tessellation);
    void dump(const ClassifiedTetrahedra& tetrahedra);
    void dump(const TessellationEdges& edges);
    void dump(const ClusterAssignment& assignment);
    void dump(const IdealVectorAssignment& assignment);
    void dump(const FacetIndex& index);
    void dump(const MappedNetwork& network);
    void dump(const SpuriousLoopElimination& elimination);
    void dump(const SkippedEdgeVectorCompletion& skipped);
    void dump(const EdgeVectorCompletionPass& pass);
    void dump(const EdgeVectorCompletion& completion);
    void dump(const StitchedLines& stitched);
    // The lines once they are smoothed and their FCC dislocation types are known.
    void dump(const std::vector<LineReconstructionDislocationLine>& finishedLines);
    // The reconstruction once it is over, timed from the loaded structure.
    void dump(const LineReconstructionDXAAlgorithm& reconstruction);
    void dump(const ExportedTable& table);

private:
    // Adds the length of the line to `length` segment by segment, the segment that closes a closed line included.
    static void addLength(const LineReconstructionDislocationLine& line, double& length);

    Stopwatch _reconstruction;
    // Edge vector completion runs from the snapshot after the loop cleanup, the last one dumped before it.
    Stopwatch _sinceSnapshot;
};

}

}
