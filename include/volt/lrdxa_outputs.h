#pragma once

#include <volt/core/lammps_parser.h>
#include <volt/plugin/output_files.h>

#include <array>
#include <string>
#include <string_view>

namespace Volt {

class StructureAnalysis;
class LineReconstructionJsonExporter;

namespace DXA {

class LineReconstructionDXAAlgorithm;
class LineReconstructionDXALog;

// The files a reconstruction writes, which are all of them whenever it is given an output base.
struct LineReconstructionExports {};

// What the files are written from: where they go and what the reconstruction built.
struct LineReconstructionOutputs {
    const std::string& base;
    const LammpsParser::Frame& frame;
    const StructureAnalysis& structure;
    const LineReconstructionDXAAlgorithm& reconstruction;
    LineReconstructionJsonExporter& exporter;
    LineReconstructionDXALog& log;

    // What writing the files weighs next to the reconstruction.
    static double weight();

    // Writes every file, one progress step per file.
    void write() const;

    std::string path(std::string_view name) const { return outputPath(base, name); }

private:
    static const std::array<OutputFile<LineReconstructionExports, LineReconstructionOutputs>, 5> files;
};

}

}
