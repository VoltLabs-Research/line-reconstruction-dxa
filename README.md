# LineReconstructionDXA

`LineReconstructionDXA` reconstructs dislocation lines from an upstream cluster package and DXA-style geometric parameters.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- LineReconstructionDXA
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `opendxa/1.0.0` (from the `OpenDXA` repository)
- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)
- `structure-identification/1.0.0` (from the `StructureIdentification` repository)
- `common-neighbor-analysis/1.0.0` (from the `CommonNeighborAnalysis` repository)
- `polyhedral-template-matching/1.0.0` (from the `PolyhedralTemplateMatching` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/line-reconstruction-dxa --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

Usage:

```bash
line-reconstruction-dxa <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--clusters-table <path>` | Yes | Path to `*_clusters.table` exported upstream. | |
| `--clusters-transitions <path>` | Yes | Path to `*_cluster_transitions.table` exported upstream. | |
| `--crystalStructure <type>` | No | Reference crystal structure: `BCC`, `FCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`, `SC`. | `FCC` |
| `--crystalPathSteps <int>` | No | Maximum crystal-path steps used for edge vectors. | `4` |
| `--tessellationGhostLayerScale <float>` | No | Ghost-layer scale relative to neighbor distance. | `3.5` |
| `--alphaScale <float>` | No | Alpha threshold scale relative to neighbor distance. | `3.5` |
| `--smoothingIterations <int>` | No | Taubin smoothing iterations for reconstructed lines. | `0` |
| `--linePointInterval <float>` | No | Line coarsening interval. | `1.2` |
| `--threads <int>` | No | Maximum worker threads. | auto capped to physical cores |
| `--help` | No | Print CLI help. | |
