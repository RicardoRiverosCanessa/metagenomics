#!/usr/bin/env python3
"""inspect_versions.py — MetaNexus version, runtime, and database inspector.

Inspects pipeline tools, container runtimes (Apptainer/Singularity), host/mamba
environments, and reference databases in a strictly read-only manner.
Supports --versions-all, --versions-deep, --versions-format text|json, and --strict-versions.
Uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Canonical identity defaults
DEFAULT_PIPELINE_NAME = "MetaNexus"
DEFAULT_PIPELINE_VERSION = "0.3.2"
DEFAULT_PIPELINE_RELEASE_STATUS = "development"

DEFAULT_TIMEOUT = 20.0
MAX_OUTPUT_BYTES = 65536

# Exit codes
EXIT_OK = 0
EXIT_STRICT_FAILED = 1
EXIT_USAGE_ERROR = 2


@dataclass
class PipelineIdentity:
    name: str = DEFAULT_PIPELINE_NAME
    version: str = DEFAULT_PIPELINE_VERSION
    release_status: str = DEFAULT_PIPELINE_RELEASE_STATUS

    def to_dict(self) -> Dict[str, str]:
        return {
            "name": self.name,
            "version": self.version,
            "release_status": self.release_status,
        }


@dataclass
class ToolResult:
    stage: str
    tool: str
    role: str
    version: str = "unknown"
    version_source: str = "unknown"
    status: str = "unknown"  # available, unavailable, unknown, error
    command: str = ""
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "tool": self.tool,
            "role": self.role,
            "version": self.version,
            "version_source": self.version_source,
            "status": self.status,
            "command": self.command,
            "error_message": self.error_message,
        }


@dataclass
class RuntimeResult:
    stage: str
    tool: str
    mode: str
    name: str
    uri: str
    path: str = "unknown"
    size_bytes: Optional[int] = None
    mtime_ns: Optional[int] = None
    fingerprint_method: str = "unknown"
    fingerprint: str = "unknown"
    status: str = "unknown"  # available, unavailable, unknown, error
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "tool": self.tool,
            "mode": self.mode,
            "name": self.name,
            "uri": self.uri,
            "path": self.path,
            "size_bytes": self.size_bytes,
            "mtime_ns": self.mtime_ns,
            "fingerprint_method": self.fingerprint_method,
            "fingerprint": self.fingerprint,
            "status": self.status,
            "error_message": self.error_message,
        }


@dataclass
class DatabaseResult:
    stage: str
    name: str
    path: str
    resolved_path: Optional[str] = None
    release: str = "unknown"
    version_source: str = "unknown"
    size_bytes: Optional[int] = None
    fingerprint_method: str = "unknown"
    fingerprint: str = "unknown"
    status: str = "unknown"  # available, unavailable, not-configured, unknown, error
    details: str = ""
    error_message: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "stage": self.stage,
            "name": self.name,
            "path": self.path,
            "resolved_path": self.resolved_path,
            "release": self.release,
            "version_source": self.version_source,
            "size_bytes": self.size_bytes,
            "fingerprint_method": self.fingerprint_method,
            "fingerprint": self.fingerprint,
            "status": self.status,
            "details": self.details,
            "error_message": self.error_message,
        }


@dataclass
class InspectionReport:
    pipeline: PipelineIdentity
    inspected_at_utc: str
    execution_mode: str
    profile_scope: str
    deep: bool
    all_stages: bool
    strict: bool
    tools: List[ToolResult] = field(default_factory=list)
    runtimes: List[RuntimeResult] = field(default_factory=list)
    databases: List[DatabaseResult] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pipeline": self.pipeline.to_dict(),
            "inspected_at_utc": self.inspected_at_utc,
            "execution_mode": self.execution_mode,
            "profile_scope": self.profile_scope,
            "deep": self.deep,
            "all_stages": self.all_stages,
            "strict": self.strict,
            "summary": self.summary,
            "tools": [t.to_dict() for t in self.tools],
            "runtimes": [r.to_dict() for r in self.runtimes],
            "databases": [d.to_dict() for d in self.databases],
        }


# Canonical tool specifications
DEFAULT_TOOL_SPECS: List[Dict[str, str]] = [
    {
        "stage": "01_fastp",
        "tool": "fastp",
        "role": "primary",
        "sif_name": "fastp.sif",
        "uri": "docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5",
        "executable": "fastp",
        "query": "--version",
        "env_name": "",
    },
    {
        "stage": "02_megahit",
        "tool": "megahit",
        "role": "primary",
        "sif_name": "megahit.sif",
        "uri": "docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_4",
        "executable": "megahit",
        "query": "--version",
        "env_name": "",
    },
    {
        "stage": "03_kraken2",
        "tool": "kraken2",
        "role": "primary",
        "sif_name": "kraken2_bracken.sif",
        "uri": "docker://staphb/kraken2:2.1.3",
        "executable": "kraken2",
        "query": "--version",
        "env_name": "",
    },
    {
        "stage": "03_kraken2",
        "tool": "bracken",
        "role": "auxiliary",
        "sif_name": "kraken2_bracken.sif",
        "uri": "docker://staphb/kraken2:2.1.3",
        "executable": "bracken",
        "query": "-v",
        "env_name": "",
    },
    {
        "stage": "04_mapping",
        "tool": "bowtie2",
        "role": "primary",
        "sif_name": "mapping_depth.sif",
        "uri": "docker://staphb/bowtie2:2.5.4",
        "executable": "bowtie2",
        "query": "--version",
        "env_name": "",
    },
    {
        "stage": "04_mapping",
        "tool": "samtools",
        "role": "auxiliary",
        "sif_name": "mapping_depth.sif",
        "uri": "docker://staphb/bowtie2:2.5.4",
        "executable": "samtools",
        "query": "--version",
        "env_name": "",
    },
    {
        "stage": "05a_vamb",
        "tool": "vamb",
        "role": "primary",
        "sif_name": "vamb.sif",
        "uri": "docker://ghcr.io/RasmussenLab/vamb:latest",
        "executable": "vamb",
        "query": "--version",
        "env_name": "vamb_env",
    },
    {
        "stage": "05b_metabinner",
        "tool": "metabinner",
        "role": "primary",
        "sif_name": "metabinner.sif",
        "uri": "docker://quay.io/biocontainers/metabinner:1.4.4--pyhdfd78af_0",
        "executable": "metabinner",
        "query": "--version",
        "env_name": "metabinner_env",
    },
    {
        "stage": "05c_comebin",
        "tool": "comebin",
        "role": "primary",
        "sif_name": "comebin.sif",
        "uri": "docker://quay.io/biocontainers/comebin:1.0.4--pyhdfd78af_0",
        "executable": "run_comebin.sh",
        "query": "-V",
        "env_name": "comebin_env",
    },
    {
        "stage": "06_refinement",
        "tool": "metawrap",
        "role": "primary",
        "sif_name": "metawrap.sif",
        "uri": "docker://quay.io/biocontainers/metawrap-mg:1.3.2--hdfd78af_1",
        "executable": "metawrap",
        "query": "--version",
        "env_name": "metawrap_env",
    },
    {
        "stage": "07_reassembly",
        "tool": "metawrap",
        "role": "primary",
        "sif_name": "metawrap.sif",
        "uri": "docker://quay.io/biocontainers/metawrap-mg:1.3.2--hdfd78af_1",
        "executable": "metawrap",
        "query": "--version",
        "env_name": "metawrap_env",
    },
    {
        "stage": "08_checkm2",
        "tool": "checkm2",
        "role": "primary",
        "sif_name": "checkm2.sif",
        "uri": "docker://quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0",
        "executable": "checkm2",
        "query": "--version",
        "env_name": "checkm2_env",
    },
    {
        "stage": "10_gtdbtk",
        "tool": "gtdbtk",
        "role": "primary",
        "sif_name": "gtdbtk.sif",
        "uri": "docker://quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_0",
        "executable": "gtdbtk",
        "query": "--version",
        "env_name": "gtdbtk_env",
    },
    {
        "stage": "11_coverm",
        "tool": "coverm",
        "role": "primary",
        "sif_name": "coverm.sif",
        "uri": "docker://quay.io/biocontainers/coverm:0.6.1--h9f5acd7_1",
        "executable": "coverm",
        "query": "--version",
        "env_name": "",
    },
]

# Canonical database specifications
DEFAULT_DB_SPECS: List[Dict[str, Any]] = [
    {
        "stage": "03_kraken2",
        "name": "kraken2-bracken",
        "env_var": "KRAKEN2_DB",
        "cli_flag": "--kraken2-db",
        "default_managed_subdir": "installed/kraken2/standard-8/20260626",
        "required_files": ["hash.k2d", "opts.k2d", "taxo.k2d"],
    },
    {
        "stage": "08_checkm2",
        "name": "checkm2",
        "env_var": "CHECKM2_DB",
        "cli_flag": "--checkm2-db",
        "is_file": True,
    },
    {
        "stage": "10_gtdbtk",
        "name": "gtdbtk",
        "env_var": "GTDBTK_DATA_PATH",
        "cli_flag": "--gtdbtk-data-path",
        "version_files": ["VERSION", "release.txt", "gtdbtk.package_version"],
    },
]


def extract_reported_version(output_text: str) -> str:
    """Extract semantic or revision version string from raw command output."""
    if not output_text:
        return "unknown"
    # Match standard semantic versions: e.g. 1.2.3, v0.23.4, 2.5.4-beta, 1.19.2
    m = re.search(r"(?<![A-Za-z0-9])v?([0-9]+(?:\.[0-9A-Za-z][0-9A-Za-z+._-]*)+)", output_text)
    if m:
        return m.group(1)
    # Match commit hashes (7 to 40 hex digits)
    m_commit = re.search(r"\b([0-9a-fA-F]{7,40})\b", output_text)
    if m_commit:
        return m_commit.group(1)
    return "unknown"


def run_command_with_limits(
    cmd: Sequence[str],
    timeout: float = DEFAULT_TIMEOUT,
    max_bytes: int = MAX_OUTPUT_BYTES,
) -> Tuple[int, str, bool]:
    """Execute a command with timeout and output capture limit.

    Returns:
        (exit_code, combined_output, timed_out)
    """
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            close_fds=True,
        )
    except FileNotFoundError:
        return (127, "command not found", False)
    except Exception as ex:
        return (1, f"failed to execute command: {ex}", False)

    timed_out = False
    output_chunks: List[bytes] = []
    bytes_read = 0

    try:
        raw_out, _ = proc.communicate(timeout=timeout)
        if raw_out:
            if len(raw_out) > max_bytes:
                raw_out = raw_out[:max_bytes]
            output_chunks.append(raw_out)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            proc.kill()
        except OSError:
            pass
        try:
            raw_out, _ = proc.communicate(timeout=1.0)
            if raw_out:
                if len(raw_out) > max_bytes:
                    raw_out = raw_out[:max_bytes]
                output_chunks.append(raw_out)
        except Exception:
            pass

    returncode = proc.returncode if proc.returncode is not None else (124 if timed_out else 1)
    combined_bytes = b"".join(output_chunks)
    text = combined_bytes.decode("utf-8", errors="replace").strip()
    return (returncode, text, timed_out)


def compute_file_sha256(path: str, max_chunk_size: int = 4 * 1024 * 1024) -> str:
    """Compute sha256 checksum of a regular file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(max_chunk_size), b""):
            h.update(chunk)
    return h.hexdigest().lower()


def compute_conda_meta_fingerprint(prefix_path: str) -> str:
    """Compute deterministic fingerprint of a conda-meta directory."""
    meta_dir = os.path.join(prefix_path, "conda-meta")
    if not os.path.isdir(meta_dir):
        return "unknown"
    rows: List[str] = []
    try:
        for name in os.listdir(meta_dir):
            if not name.endswith(".json"):
                continue
            fpath = os.path.join(meta_dir, name)
            try:
                with open(fpath, "r", encoding="utf-8") as f:
                    data = json.load(f)
                p_name = str(data.get("name") or "")
                p_ver = str(data.get("version") or "")
                p_build = str(data.get("build") or "")
                channel = str(data.get("channel") or "")
                rows.append(f"{p_name}\t{p_ver}\t{p_build}\t{channel}")
            except (OSError, ValueError, TypeError):
                continue
    except OSError:
        return "unknown"

    if not rows:
        return "unknown"
    payload = ("\n".join(sorted(rows)) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


class InspectionCache:
    """In-memory cache for runtime identities, tool versions, and database states."""

    def __init__(self) -> None:
        self.file_hashes: Dict[str, str] = {}
        self.runtime_cache: Dict[str, RuntimeResult] = {}
        self.tool_cache: Dict[str, ToolResult] = {}
        self.database_cache: Dict[str, DatabaseResult] = {}


class VersionInspector:
    """Read-only inspector for MetaNexus tools, runtimes, and databases."""

    def __init__(
        self,
        pipeline: Optional[PipelineIdentity] = None,
        engine: str = "apptainer",
        sif_dir: str = "./containers",
        deep: bool = False,
        all_stages: bool = False,
        strict: bool = False,
        timeout: float = DEFAULT_TIMEOUT,
        cache: Optional[InspectionCache] = None,
    ) -> None:
        self.pipeline = pipeline or PipelineIdentity()
        self.engine = engine.strip().lower()
        self.sif_dir = Path(sif_dir).expanduser().resolve()
        self.deep = deep
        self.all_stages = all_stages
        self.strict = strict
        self.timeout = timeout
        self.cache = cache or InspectionCache()

        # Locate apptainer / singularity on host
        self.apptainer_bin = shutil.which("apptainer") or shutil.which("singularity")
        # Locate mamba on host
        self.mamba_bin = shutil.which("mamba")

    def inspect_runtime(
        self,
        stage: str,
        tool: str,
        sif_name: str,
        uri: str,
        env_name: str,
        executable: str,
    ) -> RuntimeResult:
        cache_key = f"{self.engine}|{stage}|{tool}|{sif_name}|{env_name}"
        if cache_key in self.cache.runtime_cache:
            return self.cache.runtime_cache[cache_key]

        if self.engine == "apptainer":
            mode = "apptainer"
            name = tool
            sif_path = self.sif_dir / sif_name
            runtime_path_str = str(sif_path)
            res = RuntimeResult(
                stage=stage,
                tool=tool,
                mode=mode,
                name=name,
                uri=uri,
                path=runtime_path_str,
            )

            if not self.apptainer_bin:
                res.status = "unavailable"
                res.error_message = "Apptainer or Singularity binary not found on host"
            elif not sif_path.is_file():
                res.status = "unavailable"
                res.error_message = f"Container image '{sif_name}' not found at {sif_path}"
            elif not os.access(sif_path, os.R_OK):
                res.status = "unavailable"
                res.error_message = f"Container image '{sif_name}' is not readable"
            else:
                try:
                    st = sif_path.stat()
                    res.size_bytes = st.st_size
                    res.mtime_ns = getattr(st, "st_mtime_ns", int(st.st_mtime * 1e9))
                    res.status = "available"

                    if self.deep:
                        if runtime_path_str not in self.cache.file_hashes:
                            self.cache.file_hashes[runtime_path_str] = compute_file_sha256(runtime_path_str)
                        res.fingerprint_method = "sha256-file"
                        res.fingerprint = self.cache.file_hashes[runtime_path_str]
                    else:
                        res.fingerprint_method = "stat-quick"
                        res.fingerprint = f"size:{res.size_bytes};mtime:{res.mtime_ns}"
                except OSError as ex:
                    res.status = "error"
                    res.error_message = f"Cannot stat SIF image: {ex}"

        else:
            # Mamba mode
            mode = "mamba"
            name = env_name if env_name else "host"
            res = RuntimeResult(
                stage=stage,
                tool=tool,
                mode=mode,
                name=name,
                uri="not-applicable",
            )

            if env_name:
                if not self.mamba_bin:
                    res.status = "unavailable"
                    res.error_message = "mamba binary not found on host"
                else:
                    # Query env prefix
                    rc, out, _ = run_command_with_limits(
                        [self.mamba_bin, "run", "-n", env_name, "python", "-c", "import sys; print(sys.prefix)"],
                        timeout=self.timeout,
                    )
                    prefix = out.splitlines()[0].strip() if out and rc == 0 else ""
                    if prefix and os.path.isdir(prefix):
                        res.path = prefix
                        res.status = "available"
                        if self.deep:
                            res.fingerprint_method = "conda-meta-v1"
                            res.fingerprint = compute_conda_meta_fingerprint(prefix)
                        else:
                            res.fingerprint_method = "prefix-present"
                            res.fingerprint = "available"
                    else:
                        res.status = "unavailable"
                        res.error_message = f"Mamba environment '{env_name}' not found or invalid"
            else:
                # Host tool
                host_path = shutil.which(executable)
                if host_path:
                    res.path = host_path
                    res.status = "available"
                    try:
                        st = os.stat(host_path)
                        res.size_bytes = st.st_size
                        res.mtime_ns = getattr(st, "st_mtime_ns", int(st.st_mtime * 1e9))
                        if self.deep:
                            if host_path not in self.cache.file_hashes:
                                self.cache.file_hashes[host_path] = compute_file_sha256(host_path)
                            res.fingerprint_method = "sha256-file"
                            res.fingerprint = self.cache.file_hashes[host_path]
                        else:
                            res.fingerprint_method = "stat-quick"
                            res.fingerprint = f"size:{res.size_bytes};mtime:{res.mtime_ns}"
                    except OSError:
                        res.fingerprint_method = "which"
                        res.fingerprint = "available"
                else:
                    res.status = "unavailable"
                    res.error_message = f"Host executable '{executable}' not found in PATH"

        self.cache.runtime_cache[cache_key] = res
        return res

    def inspect_tool(
        self,
        spec: Dict[str, str],
        runtime_res: RuntimeResult,
    ) -> ToolResult:
        stage = spec["stage"]
        tool = spec["tool"]
        role = spec.get("role", "primary")
        executable = spec["executable"]
        query = spec["query"]
        env_name = spec.get("env_name", "")

        cache_key = f"{self.engine}|{stage}|{tool}|{executable}|{query}"
        if cache_key in self.cache.tool_cache:
            return self.cache.tool_cache[cache_key]

        res = ToolResult(
            stage=stage,
            tool=tool,
            role=role,
            command=f"{executable} {query}".strip(),
        )

        if runtime_res.status != "available":
            res.status = "unavailable"
            res.version = "unavailable"
            res.version_source = "unknown"
            res.error_message = runtime_res.error_message or "Runtime unavailable"
            self.cache.tool_cache[cache_key] = res
            return res

        # Run query inside runtime
        output_text = ""
        timed_out = False
        rc = 0

        if self.engine == "apptainer":
            cmd = [
                self.apptainer_bin,  # type: ignore[list-item]
                "exec",
                "--cleanenv",
                runtime_res.path,
                executable,
            ]
            if query:
                cmd.append(query)
            rc, output_text, timed_out = run_command_with_limits(cmd, timeout=self.timeout)

            # Fallback for Python-based tools like metabinner or comebin
            if (timed_out or extract_reported_version(output_text) == "unknown") and tool in ("metabinner", "comebin"):
                py_pkg = "metabinner" if tool == "metabinner" else "comebin"
                py_cmd = [
                    self.apptainer_bin,  # type: ignore[list-item]
                    "exec",
                    "--cleanenv",
                    runtime_res.path,
                    "python3",
                    "-c",
                    f"import {py_pkg}; print(getattr({py_pkg}, '__version__', getattr({py_pkg}, 'VERSION', 'unknown')))",
                ]
                py_rc, py_out, _ = run_command_with_limits(py_cmd, timeout=self.timeout)
                if py_rc == 0 and extract_reported_version(py_out) != "unknown":
                    output_text = py_out
                    res.command = f"python3 -c import {py_pkg}"

            # Fallback to SIF inspect labels if version is still unknown
            if extract_reported_version(output_text) == "unknown":
                label_cmd = [self.apptainer_bin, "inspect", "--json", runtime_res.path]  # type: ignore[list-item]
                l_rc, l_out, _ = run_command_with_limits(label_cmd, timeout=self.timeout)
                if l_rc == 0:
                    try:
                        labels = json.loads(l_out).get("data", {}).get("attributes", {}).get("labels", {})
                        for k in ("org.opencontainers.image.version", "VERSION", "version"):
                            if labels.get(k):
                                output_text = str(labels[k])
                                res.command = "apptainer inspect --json (labels)"
                                break
                    except Exception:
                        pass
        else:
            # Mamba mode
            if env_name:
                cmd = [self.mamba_bin, "run", "-n", env_name, executable]  # type: ignore[list-item]
                if query:
                    cmd.append(query)
            else:
                cmd = [executable]
                if query:
                    cmd.append(query)
            rc, output_text, timed_out = run_command_with_limits(cmd, timeout=self.timeout)

        if timed_out:
            res.status = "error"
            res.version = "timeout"
            res.version_source = "timeout"
            res.error_message = f"Query timed out after {self.timeout}s"
        else:
            version = extract_reported_version(output_text)
            if version != "unknown":
                res.status = "available"
                res.version = version
                res.version_source = "runtime-command"
            else:
                res.status = "unknown"
                res.version = "unknown"
                res.version_source = "unknown"
                if output_text:
                    res.error_message = f"Version pattern not found in output: {output_text[:120]}..."
                else:
                    res.error_message = "No output produced by version query"

        self.cache.tool_cache[cache_key] = res
        return res

    def inspect_database(
        self,
        spec: Dict[str, Any],
        user_path: Optional[str] = None,
        project_root: Optional[Path] = None,
    ) -> DatabaseResult:
        stage = spec["stage"]
        name = spec["name"]
        cache_key = f"{stage}|{name}|{user_path}"
        if cache_key in self.cache.database_cache:
            return self.cache.database_cache[cache_key]

        raw_path = user_path or os.environ.get(spec.get("env_var", ""), "")
        resolved_path: Optional[Path] = None

        # Check default managed location if user did not provide an explicit path
        if not raw_path and spec.get("default_managed_subdir") and project_root:
            default_path = project_root / "db" / spec["default_managed_subdir"]
            if default_path.exists():
                raw_path = str(default_path)

        res = DatabaseResult(
            stage=stage,
            name=name,
            path=raw_path if raw_path else "<unset>",
        )

        if not raw_path:
            res.status = "not-configured"
            res.details = f"Optional database path not configured ({spec.get('cli_flag', '')} or {spec.get('env_var', '')})"
            self.cache.database_cache[cache_key] = res
            return res

        p = Path(raw_path).expanduser().resolve()
        res.resolved_path = str(p)

        if spec.get("is_file"):
            # Single file database (e.g. CheckM2 diamond file)
            if not p.is_file():
                res.status = "unavailable"
                res.error_message = f"Database file not found: {p}"
            elif not os.access(p, os.R_OK):
                res.status = "unavailable"
                res.error_message = f"Database file is not readable: {p}"
            else:
                try:
                    st = p.stat()
                    res.size_bytes = st.st_size
                    res.status = "available"
                    res.release = p.name
                    res.version_source = "filename"
                    res.details = f"Valid diamond database file ({p.name})"
                    if self.deep:
                        p_str = str(p)
                        if p_str not in self.cache.file_hashes:
                            self.cache.file_hashes[p_str] = compute_file_sha256(p_str)
                        res.fingerprint_method = "sha256-file"
                        res.fingerprint = self.cache.file_hashes[p_str]
                    else:
                        res.fingerprint_method = "stat-quick"
                        res.fingerprint = f"size:{res.size_bytes}"
                except OSError as ex:
                    res.status = "error"
                    res.error_message = f"Cannot inspect file: {ex}"
        else:
            # Directory database (e.g. Kraken2/Bracken, GTDB-Tk)
            if not p.is_dir():
                res.status = "unavailable"
                res.error_message = f"Database directory not found: {p}"
            elif not os.access(p, os.R_OK):
                res.status = "unavailable"
                res.error_message = f"Database directory is not readable: {p}"
            else:
                try:
                    if name == "kraken2-bracken":
                        req_files = spec.get("required_files", [])
                        missing = [f for f in req_files if not (p / f).is_file()]
                        if missing:
                            res.status = "unavailable"
                            res.details = f"Missing required Kraken2 files: {', '.join(missing)}"
                        else:
                            res.status = "available"
                            # Check for installation.json
                            inst_json = p / "installation.json"
                            if inst_json.is_file():
                                try:
                                    inst_data = json.loads(inst_json.read_text(encoding="utf-8"))
                                    res.release = str(inst_data.get("release", "unknown"))
                                    res.version_source = "installation-manifest"
                                except Exception:
                                    res.version_source = "unmanaged"
                            else:
                                res.release = "unmanaged"
                                res.version_source = "directory"

                            res.details = "All core Kraken2 index files present."
                            if self.deep:
                                # Deep fingerprint: hash required and kmer_distrib files
                                key_rows = []
                                for kf in sorted(req_files + [f.name for f in p.glob("*.kmer_distrib")]):
                                    fp = p / kf
                                    if fp.is_file():
                                        fst = fp.stat()
                                        key_rows.append(f"{kf}\t{fst.st_size}\t{getattr(fst, 'st_mtime_ns', int(fst.st_mtime * 1e9))}")
                                payload = ("\n".join(key_rows) + "\n").encode("utf-8")
                                res.fingerprint_method = "metadata-keyfiles-v1"
                                res.fingerprint = hashlib.sha256(payload).hexdigest()
                            else:
                                res.fingerprint_method = "keyfiles-present"
                                res.fingerprint = "available"

                    elif name == "gtdbtk":
                        # GTDB-Tk inspection
                        version_files = spec.get("version_files", [])
                        found_rel = "unknown"
                        for vf in version_files:
                            vfp = p / vf
                            if vfp.is_file():
                                try:
                                    line = vfp.read_text(encoding="utf-8").strip().splitlines()[0]
                                    if line:
                                        found_rel = line
                                        res.version_source = f"file:{vf}"
                                        break
                                except Exception:
                                    pass
                        res.release = found_rel
                        res.status = "available"
                        res.details = "GTDB-Tk reference directory located."
                        # Prompt requirement: "No implementes todavía: fingerprint profundo de GTDB-Tk;"
                        res.fingerprint_method = "not-applicable"
                        res.fingerprint = "omitted-gtdbtk"

                except OSError as ex:
                    res.status = "error"
                    res.error_message = f"Cannot inspect database directory: {ex}"

        self.cache.database_cache[cache_key] = res
        return res

    def inspect_all(
        self,
        active_stages: Sequence[str],
        tool_specs: Sequence[Dict[str, str]] = DEFAULT_TOOL_SPECS,
        db_specs: Sequence[Dict[str, Any]] = DEFAULT_DB_SPECS,
        database_paths: Optional[Dict[str, str]] = None,
        project_root: Optional[Path] = None,
    ) -> InspectionReport:
        database_paths = database_paths or {}
        now_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        scope_name = "all" if self.all_stages else (",".join(active_stages) if active_stages else "none")

        report = InspectionReport(
            pipeline=self.pipeline,
            inspected_at_utc=now_utc,
            execution_mode=self.engine,
            profile_scope=scope_name,
            deep=self.deep,
            all_stages=self.all_stages,
            strict=self.strict,
        )

        # Filter active tool specs
        target_tools = [
            t for t in tool_specs
            if self.all_stages or t["stage"] in active_stages
        ]

        # Filter active db specs
        target_dbs = [
            d for d in db_specs
            if self.all_stages or d["stage"] in active_stages
        ]

        # Inspect runtimes and tools
        for tspec in target_tools:
            r_res = self.inspect_runtime(
                stage=tspec["stage"],
                tool=tspec["tool"],
                sif_name=tspec["sif_name"],
                uri=tspec["uri"],
                env_name=tspec.get("env_name", ""),
                executable=tspec["executable"],
            )
            report.runtimes.append(r_res)

            t_res = self.inspect_tool(spec=tspec, runtime_res=r_res)
            report.tools.append(t_res)

        # Inspect databases
        for dspec in target_dbs:
            user_p = database_paths.get(dspec["name"])
            db_res = self.inspect_database(spec=dspec, user_path=user_p, project_root=project_root)
            report.databases.append(db_res)

        # Compute summary metrics
        tool_available = sum(1 for t in report.tools if t.status == "available")
        tool_unavailable = sum(1 for t in report.tools if t.status == "unavailable")
        tool_unknown = sum(1 for t in report.tools if t.status in ("unknown", "error"))

        db_available = sum(1 for d in report.databases if d.status == "available")
        db_unavailable = sum(1 for d in report.databases if d.status in ("unavailable", "error"))
        db_not_configured = sum(1 for d in report.databases if d.status == "not-configured")

        # Status: ok if everything inspected is available; warning if some are unavailable/unknown
        has_issues = (tool_unavailable > 0 or tool_unknown > 0 or db_unavailable > 0 or (self.strict and db_not_configured > 0))
        overall_status = "ok" if not has_issues else ("failed" if self.strict else "warning")

        report.summary = {
            "total_tools": len(report.tools),
            "available_tools": tool_available,
            "unavailable_tools": tool_unavailable,
            "unknown_tools": tool_unknown,
            "total_databases": len(report.databases),
            "available_databases": db_available,
            "unavailable_databases": db_unavailable,
            "not_configured_databases": db_not_configured,
            "status": overall_status,
        }

        return report


def format_report_json(report: InspectionReport) -> str:
    """Format report as a JSON string."""
    return json.dumps(report.to_dict(), indent=2, sort_keys=True, ensure_ascii=False)


def format_report_text(report: InspectionReport) -> str:
    """Format report as clean, human-readable text."""
    lines: List[str] = []
    bar = "=" * 76
    sub_bar = "-" * 76

    lines.append(bar)
    lines.append(f"{report.pipeline.name} {report.pipeline.version} ({report.pipeline.release_status}) — Version & Runtime Inspection")
    lines.append(bar)
    lines.append(f"Execution engine  : {report.execution_mode}")
    lines.append(f"Stages in scope   : {report.profile_scope}")
    lines.append(f"Deep inspection   : {'true' if report.deep else 'false'}")
    lines.append(f"Strict mode       : {'true' if report.strict else 'false'}")
    lines.append(f"Timestamp (UTC)   : {report.inspected_at_utc}")
    lines.append("")

    # Tools section
    lines.append("[TOOLS]")
    lines.append(f"{'Stage':<14} {'Tool':<12} {'Version':<14} {'Status':<12} {'Source':<18}")
    lines.append(sub_bar)
    if not report.tools:
        lines.append("  (No tools in active scope)")
    else:
        for t in report.tools:
            lines.append(
                f"{t.stage:<14} {t.tool:<12} {t.version:<14} {t.status:<12} {t.version_source:<18}"
            )
            if t.error_message and t.status != "available":
                lines.append(f"  └─ note: {t.error_message}")
    lines.append("")

    # Runtimes section
    lines.append("[RUNTIMES]")
    lines.append(f"{'Tool':<12} {'Mode':<10} {'Status':<12} {'Path / Image':<40}")
    lines.append(sub_bar)
    if not report.runtimes:
        lines.append("  (No runtimes in active scope)")
    else:
        seen_runtimes = set()
        for r in report.runtimes:
            rt_key = (r.tool, r.mode, r.path)
            if rt_key in seen_runtimes:
                continue
            seen_runtimes.add(rt_key)
            display_path = r.path if r.path != "unknown" else r.uri
            lines.append(f"{r.tool:<12} {r.mode:<10} {r.status:<12} {display_path:<40}")
            if r.fingerprint_method not in ("unknown", "not-applicable", "stat-quick", "which", "prefix-present"):
                lines.append(f"  └─ fingerprint [{r.fingerprint_method}]: {r.fingerprint}")
            if r.error_message and r.status != "available":
                lines.append(f"  └─ error: {r.error_message}")
    lines.append("")

    # Databases section
    lines.append("[DATABASES]")
    lines.append(f"{'Stage':<14} {'Database':<16} {'Release':<14} {'Status':<14} {'Path':<20}")
    lines.append(sub_bar)
    if not report.databases:
        lines.append("  (No databases in active scope)")
    else:
        for d in report.databases:
            lines.append(f"{d.stage:<14} {d.name:<16} {d.release:<14} {d.status:<14} {d.path:<20}")
            if d.fingerprint_method not in ("unknown", "not-applicable", "stat-quick", "keyfiles-present"):
                lines.append(f"  └─ fingerprint [{d.fingerprint_method}]: {d.fingerprint}")
            if d.details:
                lines.append(f"  └─ details: {d.details}")
            if d.error_message:
                lines.append(f"  └─ error: {d.error_message}")
    lines.append("")

    # Summary
    s = report.summary
    lines.append("[SUMMARY]")
    lines.append(
        f"Tools     : {s.get('total_tools', 0)} total "
        f"({s.get('available_tools', 0)} available, "
        f"{s.get('unavailable_tools', 0)} unavailable, "
        f"{s.get('unknown_tools', 0)} unknown)"
    )
    lines.append(
        f"Databases : {s.get('total_databases', 0)} total "
        f"({s.get('available_databases', 0)} available, "
        f"{s.get('unavailable_databases', 0)} unavailable, "
        f"{s.get('not_configured_databases', 0)} not configured)"
    )
    status_label = "OK" if s.get("status") == "ok" else ("FAILED" if report.strict else "ISSUES DETECTED")
    lines.append(f"Status    : {status_label}")
    lines.append(bar)

    return "\n".join(lines)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="inspect_versions.py",
        description="Inspect MetaNexus tools, runtimes, and databases without side-effects.",
    )
    parser.add_argument("--pipeline-name", default=DEFAULT_PIPELINE_NAME, help="Pipeline canonical name.")
    parser.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION, help="Pipeline semantic version.")
    parser.add_argument("--pipeline-status", default=DEFAULT_PIPELINE_RELEASE_STATUS, help="Pipeline release status.")
    parser.add_argument("-e", "--engine", default="apptainer", choices=["apptainer", "mamba"], help="Execution engine.")
    parser.add_argument("--sif-dir", default="./containers", help="Directory containing SIF container images.")
    parser.add_argument("--stages", default="all", help="Comma-separated active stage IDs or 'all'.")
    parser.add_argument("--kraken2-db", default="", help="Path to Kraken2/Bracken database directory.")
    parser.add_argument("--checkm2-db", default="", help="Path to CheckM2 diamond database file.")
    parser.add_argument("--gtdbtk-data-path", default="", help="Path to GTDB-Tk reference data directory.")
    parser.add_argument("--all", "--versions-all", dest="versions_all", action="store_true", help="Inspect all stages.")
    parser.add_argument("--deep", "--versions-deep", dest="versions_deep", action="store_true", help="Perform deep fingerprinting.")
    parser.add_argument(
        "--format",
        "--versions-format",
        dest="versions_format",
        default="text",
        choices=["text", "json"],
        help="Output format: text or json.",
    )
    parser.add_argument(
        "--strict",
        "--strict-versions",
        dest="strict_versions",
        action="store_true",
        help="Exit with code 1 if any tool or database is unavailable or unknown.",
    )
    parser.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, help="Per-command execution timeout in seconds.")
    parser.add_argument("--spec-file", default=None, help="Optional JSON file defining tool/database specifications.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    identity = PipelineIdentity(
        name=args.pipeline_name,
        version=args.pipeline_version,
        release_status=args.pipeline_status,
    )

    # Determine stages
    if args.versions_all or args.stages.strip() == "all":
        all_stages = True
        active_stages: List[str] = []
    else:
        all_stages = False
        active_stages = [s.strip() for s in args.stages.split(",") if s.strip()]

    # Load custom specs if provided
    tool_specs = DEFAULT_TOOL_SPECS
    db_specs = DEFAULT_DB_SPECS
    if args.spec_file:
        try:
            with open(args.spec_file, "r", encoding="utf-8") as f:
                spec_data = json.load(f)
            if "tools" in spec_data:
                tool_specs = spec_data["tools"]
            if "databases" in spec_data:
                db_specs = spec_data["databases"]
        except Exception as ex:
            print(f"[ERROR] Failed to load spec file '{args.spec_file}': {ex}", file=sys.stderr)
            return EXIT_USAGE_ERROR

    database_paths = {
        "kraken2-bracken": args.kraken2_db,
        "checkm2": args.checkm2_db,
        "gtdbtk": args.gtdbtk_data_path,
    }

    project_root = Path(__file__).resolve().parent.parent

    inspector = VersionInspector(
        pipeline=identity,
        engine=args.engine,
        sif_dir=args.sif_dir,
        deep=args.versions_deep,
        all_stages=all_stages,
        strict=args.strict_versions,
        timeout=args.timeout,
    )

    report = inspector.inspect_all(
        active_stages=active_stages,
        tool_specs=tool_specs,
        db_specs=db_specs,
        database_paths=database_paths,
        project_root=project_root,
    )

    if args.versions_format == "json":
        print(format_report_json(report))
    else:
        print(format_report_text(report))

    if args.strict_versions and report.summary.get("status") == "failed":
        return EXIT_STRICT_FAILED

    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
