#pragma once

#include <volt/core/volt.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/lrdxa_options.h>
#include <volt/lrdxa_json_exporter.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt {

using json = nlohmann::json;

class LineReconstructionDXA {
public:
    LineReconstructionDXA();

    void setInputCrystalStructure(LatticeStructureType structure);
    void setClustersTablePath(std::string path);
    void setClusterTransitionsPath(std::string path);
    void setNeighborLatticePath(std::string path);
    void setCrystalPathSteps(int crystalPathSteps);
    void setTessellationGhostLayerScale(double tessellationGhostLayerScale);
    void setAlphaScale(double alphaScale);
    void setSmoothingIterations(int smoothingIterations);
    void setLinePointInterval(double linePointInterval);

    json compute(const LammpsParser::Frame& frame, const std::string& outputFile = "");

private:
    LineReconstructionDXAOptions buildOptions() const;

private:
    LatticeStructureType _inputCrystalStructure;
    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
    std::string _neighborLatticePath;
    int _crystalPathSteps;
    double _tessellationGhostLayerScale;
    double _alphaScale;
    int _smoothingIterations;
    double _linePointInterval;

    LineReconstructionJsonExporter _jsonExporter;
};

}

namespace Volt {
using LineReconstructionDXAService = LineReconstructionDXA;
}
