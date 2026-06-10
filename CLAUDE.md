# CLAUDE.md — line-reconstruction-dxa

> Read the root `../../CLAUDE.md` and `../../ECOSYSTEM-DECISIONS.md` first. Plugin-specific notes only.

## ⚠️ Licensing — BINDING (D-012)
This plugin contains **`LICENSE.OVITO`** and incorporates **OVITO-derived code**. The OVITO license terms
are binding. **Do NOT copy code from here into differently-licensed parts of the ecosystem** without
checking those terms.

## What this is
A VOLT analysis plugin: **Line Reconstruction DXA**. Reconstructs dislocation lines via DXA using cluster
tables. C++23. Package + dir `line-reconstruction-dxa`. Sits at the **top of the dislocation dependency
chain** — it pulls in opendxa + the structure providers.

## Build
```bash
# Linux shorthand (drop the boost flag on macOS/Windows — see D-004):
conan install . --build=missing -s compiler.cppstd=17 -o boost/*:without_stacktrace=True
cmake --preset conan-release && cmake --build --preset conan-release
```
`volt_add_plugin(...)` macro; C++23 from CMake, deps `cppstd=17` (D-004). No test suite.

## Dependencies (conanfile.py) — the deepest chain in the ecosystem
- Intra-ecosystem (D-005): `opendxa/[>=2.0]`, `coretoolkit/[>=2.0]`, `structure-identification/[>=2.0]`,
  `common-neighbor-analysis/[>=2.0]`, `polyhedral-template-matching/[>=2.0]`.
- Third-party (exact): `boost/1.88.0`, `onetbb/2021.12.0`, `spdlog/1.14.1`, `nlohmann_json/3.11.3`.
- The workflow `dependency_repos` must list the **full transitive closure** as PascalCase repo names —
  remember `common-neighbor-analysis` lives in repo **`AdaptiveCommonNeighborAnalysis`** (D-005).

## plugin.json (workflow DAG)
- entrypoint: `type: executable`, binary `line-reconstruction-dxa`.
- exposure: `dislocations.parquet` → `DislocationExporter`→glb (no `atoms.parquet` here).

## scripts/
- `lrdxa_plugin_wrapper.py` — locates binary, filters VOLT flags, runs, emits parquet.

## Gotchas
- Because it depends on opendxa (OVITO-derived), the licensing note applies transitively. Parquet via
  duckdb (D-003).
