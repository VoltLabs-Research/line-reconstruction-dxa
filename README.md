# Line Reconstruction DXA

Reconstructs dislocation lines via DXA from upstream cluster tables and neighbor topology.

## Install

```bash
vpm install @voltlabs/line-reconstruction-dxa
```

## CLI

```bash
line-reconstruction-dxa <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--clusters-table <path>` | yes | — | `*_clusters.table` exported upstream. |
| `--clusters-transitions <path>` | yes | — | `*_cluster_transitions.table` exported upstream. |
| `--neighbor_lattice <path>` | yes | — | `*_neighbor_lattice.parquet` exported upstream. |
| `--crystalStructure <type>` | no | `FCC` | Reference crystal structure: `BCC`, `FCC`, `HCP`, `CUBIC_DIAMOND`, `HEX_DIAMOND`, `SC`. |
| `--crystalPathSteps <int>` | no | `4` | Maximum crystal-path steps used for edge vectors. |
| `--tessellationGhostLayerScale <float>` | no | `3.5` | Ghost-layer scale relative to neighbor distance. |
| `--alphaScale <float>` | no | `3.5` | Alpha threshold scale relative to neighbor distance. |
| `--smoothingIterations <int>` | no | `0` | Taubin smoothing iterations for reconstructed lines. |
| `--linePointInterval <float>` | no | `1.2` | Line coarsening interval. |
| `--threads <int>` | no | auto capped to physical cores | Max worker threads (TBB/OMP). |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_dislocations.parquet` | Dislocation Lines | LineExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/line-reconstruction-dxa
