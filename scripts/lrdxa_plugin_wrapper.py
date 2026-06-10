#!/usr/bin/env python3
"""Volt plugin wrapper for LineReconstructionDXA. STDOUT is reserved for the
stub's binary IPC protocol — never write there from Python or subprocesses.

Pipeline integration: when the daemon supplies upstream cluster tables,
the wrapper passes them through. Standalone runs fall back to invoking the
bundled PTM binary inline."""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BINARY_NAME = "line-reconstruction-dxa"
PTM_BINARY_NAME = "polyhedral-template-matching"
PLUGIN_REPO_DIRNAME = "LineReconstructionDXA"
ENV_BINARY_OVERRIDE = "VOLT_LRDXA_BINARY"
ENV_PTM_BINARY_OVERRIDE = "VOLT_PTM_BINARY"
REQUIRED_OUTPUTS = ["_dislocations.parquet"]
LOG_TAG = "lrdxa-plugin"

SCRIPT_DIR = Path(__file__).resolve().parent
PLUGIN_ROOT = Path(os.environ.get("PLUGIN_PROJECT_DIR", SCRIPT_DIR.parent)).resolve()
PLUGINS_ROOT = PLUGIN_ROOT.parent if PLUGIN_ROOT.parent.name == "plugins" else None
EMBEDDED_LOADER = (PLUGIN_ROOT / "lib/ld-linux-x86-64.so.2").resolve()
EMBEDDED_LIBRARY_DIR = (PLUGIN_ROOT / "lib").resolve()
LATTICE_DIR_CANDIDATES = [
    (PLUGIN_ROOT / "share" / "volt" / "lattices").resolve(),
    (PLUGIN_ROOT / "lattices").resolve(),
]


class WrapperError(RuntimeError):
    pass


def _resolve_binary(binary_name: str, env_var: str, repo_subdir: str) -> Path:
    env_value = os.environ.get(env_var, "").strip()
    candidates: list[Path] = []
    if env_value:
        candidates.append(Path(env_value))
    candidates.append(PLUGIN_ROOT / "bin" / binary_name)
    if PLUGINS_ROOT is not None:
        candidates.extend([
            PLUGINS_ROOT / repo_subdir / "build/build/Release" / binary_name,
            PLUGINS_ROOT / repo_subdir / "build-local/build/Release" / binary_name,
            PLUGINS_ROOT / repo_subdir / "build-manual/build/Release" / binary_name,
        ])
    which = shutil.which(binary_name)
    if which:
        candidates.append(Path(which))
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()
    probed = "\n".join(f"  - {c}" for c in candidates)
    raise WrapperError(f"No pude resolver binario {binary_name}. Paths probados:\n{probed}")


def _resolve_embedded_runtime_command(command: list[str]) -> list[str]:
    if not command:
        return command
    if not EMBEDDED_LOADER.exists() or not EMBEDDED_LIBRARY_DIR.exists():
        return command
    binary_path = Path(command[0]).resolve()
    try:
        binary_path.relative_to(PLUGIN_ROOT)
    except ValueError:
        return command
    if binary_path == EMBEDDED_LOADER:
        return command
    return [
        str(EMBEDDED_LOADER),
        "--library-path", str(EMBEDDED_LIBRARY_DIR),
        str(binary_path),
        *command[1:],
    ]


def _run(command: list[str]) -> None:
    command = _resolve_embedded_runtime_command(command)
    sys.stderr.write(f"[{LOG_TAG}] {' '.join(shlex.quote(part) for part in command)}\n")
    sys.stderr.flush()
    completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=None)
    if completed.returncode != 0:
        raise WrapperError(f"El comando fallo con exit code {completed.returncode}: {command[0]}")


def _require_outputs(output_base: str) -> None:
    for suffix in REQUIRED_OUTPUTS:
        expected = Path(f"{output_base}{suffix}")
        if not expected.exists():
            raise WrapperError(f"Falta el archivo requerido: {expected}")


_VOLT_RUNTIME_FLAGS_WITH_VALUE = {
    "--selectedTimesteps",
    "--selected-timesteps",
}


def _filter_runtime_flags(args: list[str]) -> list[str]:
    filtered: list[str] = []
    i = 0
    while i < len(args):
        token = str(args[i])
        if token in _VOLT_RUNTIME_FLAGS_WITH_VALUE:
            i += 2
            continue
        filtered.append(token)
        i += 1
    return filtered


def _ensure_executable(path: Path) -> None:
    try:
        mode = path.stat().st_mode
        if not (mode & 0o111):
            path.chmod(mode | 0o755)
    except OSError:
        pass


def _has_flag(args: list[str], flag: str) -> bool:
    return any(a == flag for a in args)


def _flag_value(args: list[str], flag: str, default: str | None = None) -> str | None:
    for index, token in enumerate(args):
        if token == flag and index + 1 < len(args):
            return str(args[index + 1])
    return default


_STRUCTURE_TOKEN_MAP = {
    "FCC": "FCC",
    "BCC": "BCC",
    "HCP": "HCP",
    "SC": "SC",
    "DIAMOND": "CUBIC_DIAMOND",
    "CUBIC_DIAMOND": "CUBIC_DIAMOND",
    "HEX_DIAMOND": "HEX_DIAMOND",
    "HEXAGONAL_DIAMOND": "HEX_DIAMOND",
}


def _normalize_structure_token(value: str | None, default: str = "FCC") -> str:
    token = str(value or default).strip().replace("-", "_").replace(" ", "_").upper()
    return _STRUCTURE_TOKEN_MAP.get(token, default)


def _crystal_structure_from_args(args: list[str]) -> str:
    return _normalize_structure_token(_flag_value(args, "--crystalStructure", "FCC"))


def _artifact_path(artifacts: dict, name: str) -> str | None:
    artifact = artifacts.get(name)
    if not isinstance(artifact, dict):
        return None
    path = artifact.get("path")
    return str(path) if path else None


def _append_flag_if_missing(args: list[str], flag: str, value: str | None) -> list[str]:
    if not value or _has_flag(args, flag):
        return args
    return [*args, flag, value]


def _with_default_lattice_dir(args: list[str]) -> list[str]:
    if _has_flag(args, "--lattice-dir"):
        return args
    for lattice_dir in LATTICE_DIR_CANDIDATES:
        if lattice_dir.is_dir():
            return [*args, "--lattice-dir", str(lattice_dir)]
    return args


def _with_structure_id_artifacts(
    args: list[str],
    annotated_dump: str,
    clusters_table: str,
    cluster_transitions: str,
) -> list[str]:
    new_args = [annotated_dump, *args[1:]] if args else [annotated_dump]
    new_args = _append_flag_if_missing(new_args, "--clusters-table", clusters_table)
    new_args = _append_flag_if_missing(new_args, "--clusters-transitions", cluster_transitions)
    return new_args


def _find_parent_output_artifacts(args: list[str]) -> tuple[str, str, str] | None:
    if len(args) < 2:
        return None
    output_base = args[1]
    annotated_dump = f"{output_base}_annotated.dump"
    clusters_table = f"{output_base}_clusters.table"
    cluster_transitions = f"{output_base}_cluster_transitions.table"
    if all(Path(path).exists() for path in (annotated_dump, clusters_table, cluster_transitions)):
        return annotated_dump, clusters_table, cluster_transitions
    return None


def _run_inline_ptm(input_dump: str, scratch_dir: Path, crystal_structure: str) -> tuple[str, str, str]:
    ptm_binary = _resolve_binary(PTM_BINARY_NAME, ENV_PTM_BINARY_OVERRIDE, "PolyhedralTemplateMatching")
    _ensure_executable(ptm_binary)
    if EMBEDDED_LOADER.exists():
        _ensure_executable(EMBEDDED_LOADER)
    intermediate_base = str(scratch_dir / "structure_id")
    command = [str(ptm_binary), input_dump, intermediate_base, "--crystalStructure", crystal_structure, "--rmsd", "0.1"]
    _run(command)
    annotated_dump = f"{intermediate_base}_annotated.dump"
    clusters_table = f"{intermediate_base}_clusters.table"
    cluster_transitions = f"{intermediate_base}_cluster_transitions.table"
    for path in (annotated_dump, clusters_table, cluster_transitions):
        if not Path(path).exists():
            raise WrapperError(f"Inline PTM did not produce {path}")
    return annotated_dump, clusters_table, cluster_transitions


def _ensure_structure_id_artifacts(args: list[str], scratch_dir: Path) -> list[str]:
    if _has_flag(args, "--clusters-table") and _has_flag(args, "--clusters-transitions"):
        return args
    if len(args) < 2:
        raise WrapperError("Se esperaban al menos 2 argumentos: <input_dump> <output_base>")
    parent_artifacts = _find_parent_output_artifacts(args)
    if parent_artifacts is not None:
        return _with_structure_id_artifacts(args, *parent_artifacts)
    crystal_structure = _crystal_structure_from_args(args)
    annotated_dump, clusters_table, cluster_transitions = _run_inline_ptm(args[0], scratch_dir, crystal_structure)
    return _with_structure_id_artifacts(args, annotated_dump, clusters_table, cluster_transitions)


def _apply_pipeline_input_overrides(args: list, config: dict) -> list:
    artifacts = config.get("pluginInputArtifacts") or config.get("plugin_input_artifacts")
    if not isinstance(artifacts, dict):
        return args
    annotated_path = _artifact_path(artifacts, "annotatedDump")
    if not annotated_path:
        return args
    next_args = [annotated_path, *args[1:]] if args else [annotated_path]
    next_args = _append_flag_if_missing(next_args, "--clusters-table", _artifact_path(artifacts, "clustersTable"))
    next_args = _append_flag_if_missing(
        next_args,
        "--clusters-transitions",
        _artifact_path(artifacts, "clustersTransitions"),
    )
    return next_args


def _run_binary_with_args(args: list[str], scratch_dir: Path) -> dict:
    args = _filter_runtime_flags([str(a) for a in args])
    if len(args) < 2:
        raise WrapperError("Se esperaban al menos 2 argumentos: <input_dump> <output_base>")
    args = _ensure_structure_id_artifacts(args, scratch_dir)
    args = _with_default_lattice_dir(args)

    output_base = args[1]
    Path(output_base).parent.mkdir(parents=True, exist_ok=True)
    binary = _resolve_binary(BINARY_NAME, ENV_BINARY_OVERRIDE, PLUGIN_REPO_DIRNAME)
    _ensure_executable(binary)
    if EMBEDDED_LOADER.exists():
        _ensure_executable(EMBEDDED_LOADER)
    command = [str(binary), *args]
    _run(command)
    _require_outputs(output_base)
    return {
        "ok": True,
        "outputBase": output_base,
        "binary": str(binary),
        "outputs": [f"{output_base}{suffix}" for suffix in REQUIRED_OUTPUTS],
    }


def process(frame, config):
    del frame
    if not isinstance(config, dict):
        raise WrapperError("config debe ser un dict")
    args = config.get("args")
    if not isinstance(args, list):
        raise WrapperError("config['args'] debe ser una lista de strings")
    args = _apply_pipeline_input_overrides(args, config)
    with tempfile.TemporaryDirectory(prefix="lrdxa-inline-") as scratch:
        return _run_binary_with_args(args, Path(scratch))


def _main_cli() -> int:
    argv = sys.argv[1:]
    with tempfile.TemporaryDirectory(prefix="lrdxa-inline-") as scratch:
        _run_binary_with_args(argv, Path(scratch))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(_main_cli())
    except WrapperError as error:
        sys.stderr.write(f"[{LOG_TAG}] error: {error}\n")
        raise SystemExit(1)
