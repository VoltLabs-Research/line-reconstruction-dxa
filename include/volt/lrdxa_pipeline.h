#pragma once

#include <volt/core/lammps_parser.h>
#include <volt/lrdxa_options.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt {

using json = nlohmann::json;

class StructureContext;
class StructureAnalysis;
class LineReconstructionJsonExporter;

namespace DXA {

class LineReconstructionDXALog;

class LineReconstructionDXAPipeline {
public:
    static json run(
        const LammpsParser::Frame& frame,
        const std::string& outputFile,
        const LineReconstructionDXAOptions& options,
        StructureAnalysis& reconstructionStructureAnalysis,
        StructureContext& reconstructionContext,
        LineReconstructionJsonExporter& jsonExporter,
        LineReconstructionDXALog& log
    );
};

}

}
