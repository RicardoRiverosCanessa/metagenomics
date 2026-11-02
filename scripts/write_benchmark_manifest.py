#!/usr/bin/env python3
"""write_benchmark_manifest.py

Centralized, structured benchmark and execution manifest manager for the
shotgun metagenomics pipeline.

Manages:
  1. OUTDIR/logs/run_manifest.tsv (1 row per RUN_ID)
  2. OUTDIR/logs/stage_status.tsv (1 row per stage and RUN_ID)
  3. OUTDIR/run_status.tsv (compatibility view derived from stage_status.tsv)

Guarantees atomic file updates, strict schema column ordering, and monotonic
timing isolation.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional

RUN_MANIFEST_COLUMNS: List[str] = [
    "schema_version",
    "run_id",
    "pipeline_version",
    "profile",
    "execution_mode",
    "dry_run",
    "start_time_utc",
    "end_time_utc",
    "wall_time_seconds",
    "final_status",
    "exit_code",
    "outdir",
    "command_line",
    "hostname",
    "cpu_requested",
    "threads_effective",
    "memory_limit_bytes",
    "active_sample_count",
    "active_reads_source",
    "sample_manifest_path",
    "sample_manifest_fingerprint",
    "assembly_source",
    "assembly_path",
    "assembly_fingerprint",
    "bins_source_requested",
    "bins_source_resolved",
    "bins_fingerprint",
    "min_completeness",
    "max_contamination",
    "tmpdir_source",
    "tmpdir_base",
    "keep_temporary",
    "resume_plan_path",
    "stage_status_path",
    "software_inventory_path",
    "database_inventory_path",
    "error_history_path",
    "last_error_path",
    "taxonomy_summary_path",
    "mag_summary_path",
    "run_overview_path",
    "warning_count",
    "error_count",
]

STAGE_STATUS_COLUMNS: List[str] = [
    "schema_version",
    "run_id",
    "stage_id",
    "stage_name",
    "applicability",
    "planned_action",
    "checkpoint_action",
    "final_status",
    "start_time_utc",
    "end_time_utc",
    "wall_time_seconds",
    "user_cpu_seconds",
    "system_cpu_seconds",
    "cpu_utilization_percent",
    "max_rss_kib",
    "io_read_bytes",
    "io_write_bytes",
    "exit_code",
    "checkpoint_path",
    "checkpoint_status",
    "checkpoint_reason",
    "input_state",
    "validation_status",
    "primary_log",
    "primary_output",
    "output_fingerprint",
    "runtime",
    "tool_version",
    "runtime_fingerprint",
    "database_fingerprint",
    "error_id",
    "message",
]

COMPAT_RUN_STATUS_COLUMNS: List[str] = [
    "stage_internal_id",
    "stage_public_name",
    "planned_action",
    "final_status",
    "start_time",
    "end_time",
    "exit_code",
    "checkpoint_written",
]

DEFAULT_KNOWN_STAGES: List[tuple[str, str]] = [
    ("01_fastp", "Read quality control with fastp"),
    ("02_megahit", "Metagenome assembly with MEGAHIT"),
    ("03_kraken2", "Read taxonomic profiling with Kraken2/Bracken"),
    ("04_mapping", "Read mapping and coverage estimation"),
    ("05a_vamb", "Genome binning with Vamb"),
    ("05b_metabinner", "Genome binning with MetaBinner"),
    ("05c_comebin", "Genome binning with COMEBin"),
    ("06_refinement", "Consensus bin refinement with MetaWRAP"),
    ("07_reassembly", "MAG reassembly with MetaWRAP"),
    ("08_checkm2", "MAG quality assessment with CheckM2"),
    ("09_filter_mags", "MAG quality filtering"),
    ("10_gtdbtk", "Taxonomic classification with GTDB-Tk"),
    ("11_coverm", "MAG abundance estimation with CoverM"),
    ("12a_taxonomy_summary", "Read-based taxonomic summary"),
    ("12b_mag_summary", "MAG quality, taxonomy, and abundance summary"),
    ("12c_run_overview", "Pipeline run overview and execution metrics"),
]


def sanitize_tsv(val: Any) -> str:
    """Sanitize string for TSV cell safety: replace newlines/tabs with space."""
    if val is None:
        return "NA"
    s = str(val).strip()
    if not s:
        return "NA"
    return re.sub(r"[\r\n\t]+", " ", s).strip()


def atomic_write_tsv(path: str, columns: List[str], rows: List[Dict[str, Any]]) -> None:
    """Atomically write a TSV file via a temporary file in the same directory."""
    parent_dir = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=parent_dir, prefix=".tmp_tsv_")
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for r in rows:
                writer.writerow([r.get(c, "NA") for c in columns])
        os.replace(tmp_path, path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def load_tsv(path: str) -> List[Dict[str, str]]:
    """Load an existing TSV file if present, returning a list of dicts."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return []
    rows: List[Dict[str, str]] = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            return []
        for row in reader:
            if not row or not any(row):
                continue
            d: Dict[str, str] = {}
            for i, col in enumerate(header):
                d[col] = row[i] if i < len(row) else "NA"
            rows.append(d)
    return rows


def parse_gnu_time_file(path: str) -> Dict[str, Any]:
    """Parse output of GNU time -v."""
    metrics: Dict[str, Any] = {
        "user_cpu_seconds": 0.0,
        "system_cpu_seconds": 0.0,
        "cpu_utilization_percent": 0.0,
        "max_rss_kib": 0,
        "found": False,
    }
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return metrics

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if "User time (seconds):" in line:
                    m = re.search(r"User time \(seconds\):\s*([0-9.]+)", line)
                    if m:
                        metrics["user_cpu_seconds"] = float(m.group(1))
                        metrics["found"] = True
                elif "System time (seconds):" in line:
                    m = re.search(r"System time \(seconds\):\s*([0-9.]+)", line)
                    if m:
                        metrics["system_cpu_seconds"] = float(m.group(1))
                        metrics["found"] = True
                elif "Percent of CPU this job got:" in line:
                    m = re.search(r"Percent of CPU this job got:\s*([0-9]+)%", line)
                    if m:
                        metrics["cpu_utilization_percent"] = float(m.group(1))
                        metrics["found"] = True
                elif "Maximum resident set size (kbytes):" in line:
                    m = re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", line)
                    if m:
                        metrics["max_rss_kib"] = int(m.group(1))
                        metrics["found"] = True
    except Exception:
        pass
    return metrics


def aggregate_metrics_from_dir(metrics_dir: str) -> Dict[str, Any]:
    """Aggregate metrics across multiple GNU time log files in a directory."""
    if not os.path.isdir(metrics_dir):
        return {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
                "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
                "io_read_bytes": "NA", "io_write_bytes": "NA"}

    time_files = [
        os.path.join(metrics_dir, f)
        for f in os.listdir(metrics_dir)
        if f.startswith("time_") and f.endswith(".txt")
    ]
    if not time_files:
        return {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
                "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
                "io_read_bytes": "NA", "io_write_bytes": "NA"}

    total_user = 0.0
    total_sys = 0.0
    max_rss = 0
    max_cpu_pct = 0.0
    measured = False

    for tf in time_files:
        parsed = parse_gnu_time_file(tf)
        if parsed["found"]:
            measured = True
            total_user += parsed["user_cpu_seconds"]
            total_sys += parsed["system_cpu_seconds"]
            if parsed["max_rss_kib"] > max_rss:
                max_rss = parsed["max_rss_kib"]
            if parsed["cpu_utilization_percent"] > max_cpu_pct:
                max_cpu_pct = parsed["cpu_utilization_percent"]

    if not measured:
        return {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
                "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
                "io_read_bytes": "NA", "io_write_bytes": "NA"}

    return {
        "user_cpu_seconds": f"{total_user:.3f}",
        "system_cpu_seconds": f"{total_sys:.3f}",
        "cpu_utilization_percent": f"{max_cpu_pct:.1f}",
        "max_rss_kib": str(max_rss),
        "io_read_bytes": "NA",
        "io_write_bytes": "NA",
    }


def write_compat_run_status(compat_path: str, stage_rows: List[Dict[str, Any]], run_id: str) -> None:
    """Generate run_status.tsv compatibility view from stage_status rows for run_id."""
    if not compat_path:
        return
    rows_for_run = [r for r in stage_rows if r.get("run_id") == run_id]
    if not rows_for_run:
        rows_for_run = stage_rows

    compat_rows: List[Dict[str, Any]] = []
    for sr in rows_for_run:
        st_utc = sr.get("start_time_utc", "NA")
        end_utc = sr.get("end_time_utc", "NA")
        exit_c = sr.get("exit_code", "NA")
        chk_stat = sr.get("checkpoint_status", "missing")
        chk_written = "true" if chk_stat == "written" else "false"

        compat_rows.append({
            "stage_internal_id": sr.get("stage_id", "NA"),
            "stage_public_name": sr.get("stage_name", "NA"),
            "planned_action": sr.get("planned_action", "NA"),
            "final_status": sr.get("final_status", "NA"),
            "start_time": "" if st_utc == "NA" else st_utc,
            "end_time": "" if end_utc == "NA" else end_utc,
            "exit_code": "" if exit_c == "NA" else exit_c,
            "checkpoint_written": chk_written,
        })

    parent_dir = os.path.dirname(os.path.abspath(compat_path))
    os.makedirs(parent_dir, exist_ok=True)
    atomic_write_tsv(compat_path, COMPAT_RUN_STATUS_COLUMNS, compat_rows)


def parse_resume_plan(resume_plan_path: str) -> Dict[str, Dict[str, str]]:
    """Parse resume_plan.tsv if available to extract planned actions and reasons."""
    plan: Dict[str, Dict[str, str]] = {}
    if not resume_plan_path or not os.path.isfile(resume_plan_path):
        return plan
    try:
        with open(resume_plan_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if not header:
                return plan
            stg_idx = 0
            act_idx = 4
            reas_idx = 5
            for i, col in enumerate(header):
                if col in ("stage_internal_id", "stage_id", "stage"):
                    stg_idx = i
                elif col in ("planned_action", "action"):
                    act_idx = i
                elif col == "reason":
                    reas_idx = i
            for row in reader:
                if not row or not any(row):
                    continue
                if len(row) > stg_idx:
                    stg = row[stg_idx].strip()
                    action = row[act_idx].strip() if len(row) > act_idx else "run"
                    reason = row[reas_idx].strip() if len(row) > reas_idx else "none"
                    plan[stg] = {"action": action, "reason": reason}
    except Exception:
        pass
    return plan


# ── Subcommand Implementations ───────────────────────────────────────────────

def cmd_init(args: argparse.Namespace) -> int:
    """Initialize run_manifest.tsv and stage_status.tsv for the execution."""
    run_id = args.run_id
    if not run_id:
        print("[ERROR] --run-id is required for init.", file=sys.stderr)
        return 1

    # 1. Update/Add run_manifest row
    run_rows = load_tsv(args.run_manifest_path)
    # Check if run_id already exists in run_manifest
    found_run = False
    new_run_row: Dict[str, Any] = {
        "schema_version": "1",
        "run_id": run_id,
        "pipeline_version": sanitize_tsv(args.pipeline_version or "unknown"),
        "profile": sanitize_tsv(args.profile or "custom"),
        "execution_mode": sanitize_tsv(args.execution_mode or "apptainer"),
        "dry_run": "true" if args.dry_run else "false",
        "start_time_utc": sanitize_tsv(args.start_time_utc or "NA"),
        "end_time_utc": "NA",
        "wall_time_seconds": "NA",
        "final_status": "running",
        "exit_code": "NA",
        "outdir": sanitize_tsv(args.outdir or "NA"),
        "command_line": sanitize_tsv(args.command_line or "NA"),
        "hostname": sanitize_tsv(args.hostname or "unknown"),
        "cpu_requested": sanitize_tsv(args.cpu_requested or "NA"),
        "threads_effective": sanitize_tsv(args.threads_effective or "1"),
        "memory_limit_bytes": sanitize_tsv(args.memory_limit_bytes or "NA"),
        "active_sample_count": sanitize_tsv(args.active_sample_count or "0"),
        "active_reads_source": sanitize_tsv(args.active_reads_source or "none"),
        "sample_manifest_path": sanitize_tsv(args.sample_manifest_path or "NA"),
        "sample_manifest_fingerprint": sanitize_tsv(args.sample_manifest_fingerprint or "NA"),
        "assembly_source": sanitize_tsv(args.assembly_source or "none"),
        "assembly_path": sanitize_tsv(args.assembly_path or "NA"),
        "assembly_fingerprint": sanitize_tsv(args.assembly_fingerprint or "NA"),
        "bins_source_requested": sanitize_tsv(args.bins_source_requested or "NA"),
        "bins_source_resolved": sanitize_tsv(args.bins_source_resolved or "NA"),
        "bins_fingerprint": sanitize_tsv(args.bins_fingerprint or "NA"),
        "min_completeness": sanitize_tsv(args.min_completeness or "NA"),
        "max_contamination": sanitize_tsv(args.max_contamination or "NA"),
        "tmpdir_source": sanitize_tsv(args.tmpdir_source or "NA"),
        "tmpdir_base": sanitize_tsv(args.tmpdir_base or "NA"),
        "keep_temporary": "true" if args.keep_temporary else "false",
        "resume_plan_path": sanitize_tsv(args.resume_plan_path or "NA"),
        "stage_status_path": sanitize_tsv(args.stage_status_path or "NA"),
        "software_inventory_path": sanitize_tsv(args.software_inventory_path or "NA"),
        "database_inventory_path": sanitize_tsv(args.database_inventory_path or "NA"),
        "error_history_path": sanitize_tsv(args.error_history_path or "NA"),
        "last_error_path": sanitize_tsv(args.last_error_path or "NA"),
        "taxonomy_summary_path": sanitize_tsv(args.taxonomy_summary_path or "NA"),
        "mag_summary_path": sanitize_tsv(args.mag_summary_path or "NA"),
        "run_overview_path": sanitize_tsv(args.run_overview_path or "NA"),
        "warning_count": sanitize_tsv(args.warning_count or "0"),
        "error_count": sanitize_tsv(args.error_count or "0"),
    }

    for i, r in enumerate(run_rows):
        if r.get("run_id") == run_id:
            run_rows[i] = new_run_row
            found_run = True
            break
    if not found_run:
        run_rows.append(new_run_row)

    atomic_write_tsv(args.run_manifest_path, RUN_MANIFEST_COLUMNS, run_rows)

    # 2. Update/Add stage_status rows
    plan = parse_resume_plan(args.resume_plan_path)
    stage_rows = load_tsv(args.stage_status_path)

    # Filter out any prior entries for this run_id to ensure clean idempotence
    stage_rows = [r for r in stage_rows if r.get("run_id") != run_id]

    # Build stage rows in deterministic order
    known_stages = DEFAULT_KNOWN_STAGES
    for stg_id, stg_name in known_stages:
        stg_plan = plan.get(stg_id, {})
        action = stg_plan.get("action", "run")
        reason = stg_plan.get("reason", "none")

        if action == "disabled":
            applicability = "disabled"
            planned_action = "disabled"
            chk_action = "skipped-disabled"
            final_st = "disabled"
        elif action == "not-applicable":
            applicability = "not-applicable"
            planned_action = "not-applicable"
            chk_action = "skipped-not-applicable"
            final_st = "not-applicable"
        else:
            applicability = "applicable"
            planned_action = action
            chk_action = "not-reached"
            final_st = "pending"

        chk_path = os.path.join(args.checkpoint_dir or "", f"{stg_id}.done") if args.checkpoint_dir else "NA"

        stage_rows.append({
            "schema_version": "1",
            "run_id": run_id,
            "stage_id": stg_id,
            "stage_name": stg_name,
            "applicability": applicability,
            "planned_action": planned_action,
            "checkpoint_action": chk_action,
            "final_status": final_st,
            "start_time_utc": "NA",
            "end_time_utc": "NA",
            "wall_time_seconds": "NA",
            "user_cpu_seconds": "NA",
            "system_cpu_seconds": "NA",
            "cpu_utilization_percent": "NA",
            "max_rss_kib": "NA",
            "io_read_bytes": "NA",
            "io_write_bytes": "NA",
            "exit_code": "NA",
            "checkpoint_path": chk_path,
            "checkpoint_status": "missing",
            "checkpoint_reason": sanitize_tsv(reason),
            "input_state": "pending",
            "validation_status": "pending" if applicability == "applicable" else "not-applicable",
            "primary_log": "NA",
            "primary_output": "NA",
            "output_fingerprint": "NA",
            "runtime": sanitize_tsv(args.execution_mode or "apptainer"),
            "tool_version": "NA",
            "runtime_fingerprint": "NA",
            "database_fingerprint": "NA",
            "error_id": "none",
            "message": "none",
        })

    atomic_write_tsv(args.stage_status_path, STAGE_STATUS_COLUMNS, stage_rows)
    if args.compat_run_status_path:
        write_compat_run_status(args.compat_run_status_path, stage_rows, run_id)

    return 0


def cmd_stage_start(args: argparse.Namespace) -> int:
    """Mark an applicable stage as running."""
    run_id = args.run_id
    stage_id = args.stage_id
    stage_rows = load_tsv(args.stage_status_path)
    updated = False

    for r in stage_rows:
        if r.get("run_id") == run_id and r.get("stage_id") == stage_id:
            r["final_status"] = "running"
            r["start_time_utc"] = sanitize_tsv(args.start_time_utc or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            updated = True
            break

    if updated:
        atomic_write_tsv(args.stage_status_path, STAGE_STATUS_COLUMNS, stage_rows)
        if args.compat_run_status_path:
            write_compat_run_status(args.compat_run_status_path, stage_rows, run_id)
    return 0


def cmd_stage_finish(args: argparse.Namespace) -> int:
    """Mark a stage as successfully completed or reused."""
    run_id = args.run_id
    stage_id = args.stage_id
    action = args.planned_action or "run"
    stage_rows = load_tsv(args.stage_status_path)
    updated = False

    # Collect resource metrics
    if action == "reuse":
        metrics = {
            "user_cpu_seconds": "NA",
            "system_cpu_seconds": "NA",
            "cpu_utilization_percent": "NA",
            "max_rss_kib": "NA",
            "io_read_bytes": "NA",
            "io_write_bytes": "NA",
        }
        start_utc = "NA"
        end_utc = "NA"
        wall_time = "NA"
        chk_action = "reused"
        final_st = "reused"
        chk_stat = "reused"
        exit_code = "0"
    else:
        # Executed stage
        chk_action = args.checkpoint_action or ("forced" if action == "force" else ("restarted" if action == "restart" else "executed"))
        final_st = args.final_status or "success"
        start_utc = sanitize_tsv(args.start_time_utc)
        end_utc = sanitize_tsv(args.end_time_utc)
        wall_time = sanitize_tsv(args.wall_time_seconds)
        chk_stat = sanitize_tsv(args.checkpoint_status or "written")
        exit_code = str(args.exit_code if args.exit_code is not None else 0)

        if args.metrics_dir and os.path.isdir(args.metrics_dir):
            metrics = aggregate_metrics_from_dir(args.metrics_dir)
        elif args.time_file and os.path.isfile(args.time_file):
            parsed = parse_gnu_time_file(args.time_file)
            if parsed["found"]:
                metrics = {
                    "user_cpu_seconds": f"{parsed['user_cpu_seconds']:.3f}",
                    "system_cpu_seconds": f"{parsed['system_cpu_seconds']:.3f}",
                    "cpu_utilization_percent": f"{parsed['cpu_utilization_percent']:.1f}",
                    "max_rss_kib": str(parsed["max_rss_kib"]),
                    "io_read_bytes": "NA",
                    "io_write_bytes": "NA",
                }
            else:
                metrics = {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
                           "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
                           "io_read_bytes": "NA", "io_write_bytes": "NA"}
        else:
            metrics = {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
                       "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
                       "io_read_bytes": "NA", "io_write_bytes": "NA"}

    for r in stage_rows:
        if r.get("run_id") == run_id and r.get("stage_id") == stage_id:
            r["checkpoint_action"] = chk_action
            r["final_status"] = final_st
            r["start_time_utc"] = start_utc
            r["end_time_utc"] = end_utc
            r["wall_time_seconds"] = wall_time
            r["user_cpu_seconds"] = metrics["user_cpu_seconds"]
            r["system_cpu_seconds"] = metrics["system_cpu_seconds"]
            r["cpu_utilization_percent"] = metrics["cpu_utilization_percent"]
            r["max_rss_kib"] = metrics["max_rss_kib"]
            r["io_read_bytes"] = metrics["io_read_bytes"]
            r["io_write_bytes"] = metrics["io_write_bytes"]
            r["exit_code"] = exit_code
            r["checkpoint_status"] = chk_stat
            if args.checkpoint_path:
                r["checkpoint_path"] = sanitize_tsv(args.checkpoint_path)
            if args.checkpoint_reason:
                r["checkpoint_reason"] = sanitize_tsv(args.checkpoint_reason)
            if args.input_state:
                r["input_state"] = sanitize_tsv(args.input_state)
            if args.validation_status:
                r["validation_status"] = sanitize_tsv(args.validation_status)
            if args.primary_log:
                r["primary_log"] = sanitize_tsv(args.primary_log)
            if args.primary_output:
                r["primary_output"] = sanitize_tsv(args.primary_output)
            if args.output_fingerprint:
                r["output_fingerprint"] = sanitize_tsv(args.output_fingerprint)
            if args.runtime:
                r["runtime"] = sanitize_tsv(args.runtime)
            if args.tool_version:
                r["tool_version"] = sanitize_tsv(args.tool_version)
            if args.runtime_fingerprint:
                r["runtime_fingerprint"] = sanitize_tsv(args.runtime_fingerprint)
            if args.database_fingerprint:
                r["database_fingerprint"] = sanitize_tsv(args.database_fingerprint)
            if args.message:
                r["message"] = sanitize_tsv(args.message)
            updated = True
            break

    if updated:
        atomic_write_tsv(args.stage_status_path, STAGE_STATUS_COLUMNS, stage_rows)
        if args.compat_run_status_path:
            write_compat_run_status(args.compat_run_status_path, stage_rows, run_id)
    return 0


def cmd_stage_fail(args: argparse.Namespace) -> int:
    """Record stage failure, marking downstream pending stages as not-reached."""
    run_id = args.run_id
    stage_id = args.stage_id
    stage_rows = load_tsv(args.stage_status_path)

    metrics = {"user_cpu_seconds": "NA", "system_cpu_seconds": "NA",
               "cpu_utilization_percent": "NA", "max_rss_kib": "NA",
               "io_read_bytes": "NA", "io_write_bytes": "NA"}
    if args.metrics_dir and os.path.isdir(args.metrics_dir):
        metrics = aggregate_metrics_from_dir(args.metrics_dir)

    for r in stage_rows:
        if r.get("run_id") != run_id:
            continue
        stg = r.get("stage_id")
        if stg == stage_id:
            r["checkpoint_action"] = "failed-before-checkpoint"
            r["final_status"] = "failed"
            r["start_time_utc"] = sanitize_tsv(args.start_time_utc)
            r["end_time_utc"] = sanitize_tsv(args.end_time_utc or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            r["wall_time_seconds"] = sanitize_tsv(args.wall_time_seconds)
            r["user_cpu_seconds"] = metrics["user_cpu_seconds"]
            r["system_cpu_seconds"] = metrics["system_cpu_seconds"]
            r["cpu_utilization_percent"] = metrics["cpu_utilization_percent"]
            r["max_rss_kib"] = metrics["max_rss_kib"]
            r["io_read_bytes"] = "NA"
            r["io_write_bytes"] = "NA"
            r["exit_code"] = str(args.exit_code if args.exit_code is not None else 1)
            r["validation_status"] = "failed"
            r["error_id"] = sanitize_tsv(args.error_id or "none")
            r["message"] = sanitize_tsv(args.message or "stage execution failed")
            if args.primary_log:
                r["primary_log"] = sanitize_tsv(args.primary_log)
        elif r.get("final_status") == "pending" or r.get("checkpoint_action") == "not-reached":
            # Downstream unexecuted stages
            r["final_status"] = "not-reached"
            r["checkpoint_action"] = "not-reached"
            r["start_time_utc"] = "NA"
            r["end_time_utc"] = "NA"
            r["wall_time_seconds"] = "NA"
            r["user_cpu_seconds"] = "NA"
            r["system_cpu_seconds"] = "NA"
            r["cpu_utilization_percent"] = "NA"
            r["max_rss_kib"] = "NA"
            r["io_read_bytes"] = "NA"
            r["io_write_bytes"] = "NA"
            r["exit_code"] = "NA"
            r["validation_status"] = "not-applicable"

    atomic_write_tsv(args.stage_status_path, STAGE_STATUS_COLUMNS, stage_rows)
    if args.compat_run_status_path:
        write_compat_run_status(args.compat_run_status_path, stage_rows, run_id)
    return 0


def cmd_stage_signal(args: argparse.Namespace) -> int:
    """Record interruption by signal, marking downstream pending stages as not-reached."""
    run_id = args.run_id
    stage_id = args.stage_id
    sig = (args.signal or "INT").upper()
    status_label = "interrupted" if sig == "INT" else "terminated"
    stage_rows = load_tsv(args.stage_status_path)

    for r in stage_rows:
        if r.get("run_id") != run_id:
            continue
        stg = r.get("stage_id")
        if stg == stage_id:
            r["checkpoint_action"] = "failed-before-checkpoint"
            r["final_status"] = status_label
            r["end_time_utc"] = sanitize_tsv(args.end_time_utc or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            r["exit_code"] = str(args.exit_code if args.exit_code is not None else 130)
            r["validation_status"] = "failed"
            r["error_id"] = sanitize_tsv(args.error_id or "none")
            r["message"] = sanitize_tsv(args.message or f"pipeline interrupted by signal {sig}")
        elif r.get("final_status") == "pending":
            r["final_status"] = "not-reached"
            r["checkpoint_action"] = "not-reached"
            r["start_time_utc"] = "NA"
            r["end_time_utc"] = "NA"
            r["wall_time_seconds"] = "NA"
            r["user_cpu_seconds"] = "NA"
            r["system_cpu_seconds"] = "NA"
            r["cpu_utilization_percent"] = "NA"
            r["max_rss_kib"] = "NA"
            r["io_read_bytes"] = "NA"
            r["io_write_bytes"] = "NA"
            r["exit_code"] = "NA"
            r["validation_status"] = "not-applicable"

    atomic_write_tsv(args.stage_status_path, STAGE_STATUS_COLUMNS, stage_rows)
    if args.compat_run_status_path:
        write_compat_run_status(args.compat_run_status_path, stage_rows, run_id)
    return 0


def cmd_finalize_run(args: argparse.Namespace) -> int:
    """Finalize run_manifest.tsv with final status, monotonic duration, and exit code."""
    run_id = args.run_id
    run_rows = load_tsv(args.run_manifest_path)
    updated = False

    for r in run_rows:
        if r.get("run_id") == run_id:
            r["final_status"] = sanitize_tsv(args.final_status or "success")
            r["exit_code"] = str(args.exit_code if args.exit_code is not None else 0)
            r["end_time_utc"] = sanitize_tsv(args.end_time_utc or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            r["wall_time_seconds"] = sanitize_tsv(args.wall_time_seconds or "NA")
            if args.warning_count is not None:
                r["warning_count"] = sanitize_tsv(args.warning_count)
            if args.error_count is not None:
                r["error_count"] = sanitize_tsv(args.error_count)
            if args.taxonomy_summary_path:
                r["taxonomy_summary_path"] = sanitize_tsv(args.taxonomy_summary_path)
            if args.mag_summary_path:
                r["mag_summary_path"] = sanitize_tsv(args.mag_summary_path)
            if args.run_overview_path:
                r["run_overview_path"] = sanitize_tsv(args.run_overview_path)
            updated = True
            break

    if updated:
        atomic_write_tsv(args.run_manifest_path, RUN_MANIFEST_COLUMNS, run_rows)
    return 0


def cmd_sync_compat(args: argparse.Namespace) -> int:
    """Synchronize run_status.tsv compatibility view from stage_status.tsv."""
    stage_rows = load_tsv(args.stage_status_path)
    write_compat_run_status(args.compat_run_status_path, stage_rows, args.run_id)
    return 0


# ── Argument Parsing ─────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pipeline benchmark and manifest manager")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # init
    p_init = subparsers.add_parser("init")
    p_init.add_argument("--run-manifest-path", required=True)
    p_init.add_argument("--stage-status-path", required=True)
    p_init.add_argument("--compat-run-status-path", default="")
    p_init.add_argument("--resume-plan-path", default="")
    p_init.add_argument("--checkpoint-dir", default="")
    p_init.add_argument("--run-id", required=True)
    p_init.add_argument("--pipeline-version", default="unknown")
    p_init.add_argument("--profile", default="custom")
    p_init.add_argument("--execution-mode", default="apptainer")
    p_init.add_argument("--dry-run", action="store_true")
    p_init.add_argument("--start-time-utc", default="")
    p_init.add_argument("--outdir", default="")
    p_init.add_argument("--command-line", default="")
    p_init.add_argument("--hostname", default="")
    p_init.add_argument("--cpu-requested", default="NA")
    p_init.add_argument("--threads-effective", default="1")
    p_init.add_argument("--memory-limit-bytes", default="NA")
    p_init.add_argument("--active-sample-count", default="0")
    p_init.add_argument("--active-reads-source", default="none")
    p_init.add_argument("--sample-manifest-path", default="NA")
    p_init.add_argument("--sample-manifest-fingerprint", default="NA")
    p_init.add_argument("--assembly-source", default="none")
    p_init.add_argument("--assembly-path", default="NA")
    p_init.add_argument("--assembly-fingerprint", default="NA")
    p_init.add_argument("--bins-source-requested", default="NA")
    p_init.add_argument("--bins-source-resolved", default="NA")
    p_init.add_argument("--bins-fingerprint", default="NA")
    p_init.add_argument("--min-completeness", default="NA")
    p_init.add_argument("--max-contamination", default="NA")
    p_init.add_argument("--tmpdir-source", default="NA")
    p_init.add_argument("--tmpdir-base", default="NA")
    p_init.add_argument("--keep-temporary", action="store_true")
    p_init.add_argument("--software-inventory-path", default="NA")
    p_init.add_argument("--database-inventory-path", default="NA")
    p_init.add_argument("--error-history-path", default="NA")
    p_init.add_argument("--last-error-path", default="NA")
    p_init.add_argument("--taxonomy-summary-path", default="NA")
    p_init.add_argument("--mag-summary-path", "--mag_summary-path", dest="mag_summary_path", default="NA")
    p_init.add_argument("--run-overview-path", default="NA")
    p_init.add_argument("--warning-count", default="0")
    p_init.add_argument("--error-count", default="0")

    # stage-start
    p_stg_start = subparsers.add_parser("stage-start")
    p_stg_start.add_argument("--stage-status-path", required=True)
    p_stg_start.add_argument("--compat-run-status-path", default="")
    p_stg_start.add_argument("--run-id", required=True)
    p_stg_start.add_argument("--stage-id", required=True)
    p_stg_start.add_argument("--start-time-utc", default="")

    # stage-finish
    p_stg_fin = subparsers.add_parser("stage-finish")
    p_stg_fin.add_argument("--stage-status-path", required=True)
    p_stg_fin.add_argument("--compat-run-status-path", default="")
    p_stg_fin.add_argument("--run-id", required=True)
    p_stg_fin.add_argument("--stage-id", required=True)
    p_stg_fin.add_argument("--planned-action", default="run")
    p_stg_fin.add_argument("--checkpoint-action", default="")
    p_stg_fin.add_argument("--final-status", default="success")
    p_stg_fin.add_argument("--start-time-utc", default="")
    p_stg_fin.add_argument("--end-time-utc", default="")
    p_stg_fin.add_argument("--wall-time-seconds", default="NA")
    p_stg_fin.add_argument("--metrics-dir", default="")
    p_stg_fin.add_argument("--time-file", default="")
    p_stg_fin.add_argument("--exit-code", type=int, default=0)
    p_stg_fin.add_argument("--checkpoint-path", default="")
    p_stg_fin.add_argument("--checkpoint-status", default="")
    p_stg_fin.add_argument("--checkpoint-reason", default="")
    p_stg_fin.add_argument("--input-state", default="")
    p_stg_fin.add_argument("--validation-status", default="passed")
    p_stg_fin.add_argument("--primary-log", default="")
    p_stg_fin.add_argument("--primary-output", default="")
    p_stg_fin.add_argument("--output-fingerprint", default="")
    p_stg_fin.add_argument("--runtime", default="")
    p_stg_fin.add_argument("--tool-version", default="")
    p_stg_fin.add_argument("--runtime-fingerprint", default="")
    p_stg_fin.add_argument("--database-fingerprint", default="")
    p_stg_fin.add_argument("--message", default="")

    # stage-fail
    p_stg_fail = subparsers.add_parser("stage-fail")
    p_stg_fail.add_argument("--stage-status-path", required=True)
    p_stg_fail.add_argument("--compat-run-status-path", default="")
    p_stg_fail.add_argument("--run-id", required=True)
    p_stg_fail.add_argument("--stage-id", required=True)
    p_stg_fail.add_argument("--planned-action", default="run")
    p_stg_fail.add_argument("--exit-code", type=int, default=1)
    p_stg_fail.add_argument("--start-time-utc", default="")
    p_stg_fail.add_argument("--end-time-utc", default="")
    p_stg_fail.add_argument("--wall-time-seconds", default="NA")
    p_stg_fail.add_argument("--metrics-dir", default="")
    p_stg_fail.add_argument("--primary-log", default="")
    p_stg_fail.add_argument("--error-id", default="")
    p_stg_fail.add_argument("--message", default="")

    # stage-signal
    p_stg_sig = subparsers.add_parser("stage-signal")
    p_stg_sig.add_argument("--stage-status-path", required=True)
    p_stg_sig.add_argument("--compat-run-status-path", default="")
    p_stg_sig.add_argument("--run-id", required=True)
    p_stg_sig.add_argument("--stage-id", required=True)
    p_stg_sig.add_argument("--signal", default="INT")
    p_stg_sig.add_argument("--exit-code", type=int, default=130)
    p_stg_sig.add_argument("--end-time-utc", default="")
    p_stg_sig.add_argument("--error-id", default="")
    p_stg_sig.add_argument("--message", default="")

    # finalize-run
    p_fin = subparsers.add_parser("finalize-run")
    p_fin.add_argument("--run-manifest-path", required=True)
    p_fin.add_argument("--run-id", required=True)
    p_fin.add_argument("--final-status", default="success")
    p_fin.add_argument("--exit-code", type=int, default=0)
    p_fin.add_argument("--end-time-utc", default="")
    p_fin.add_argument("--wall-time-seconds", default="NA")
    p_fin.add_argument("--warning-count", default=None)
    p_fin.add_argument("--error-count", default=None)
    p_fin.add_argument("--taxonomy-summary-path", default="")
    p_fin.add_argument("--mag-summary-path", "--mag_summary-path", dest="mag_summary_path", default="")
    p_fin.add_argument("--run-overview-path", default="")

    # sync-compat
    p_sync = subparsers.add_parser("sync-compat")
    p_sync.add_argument("--stage-status-path", required=True)
    p_sync.add_argument("--compat-run-status-path", required=True)
    p_sync.add_argument("--run-id", required=True)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.subcommand == "init":
        return cmd_init(args)
    elif args.subcommand == "stage-start":
        return cmd_stage_start(args)
    elif args.subcommand == "stage-finish":
        return cmd_stage_finish(args)
    elif args.subcommand == "stage-fail":
        return cmd_stage_fail(args)
    elif args.subcommand == "stage-signal":
        return cmd_stage_signal(args)
    elif args.subcommand == "finalize-run":
        return cmd_finalize_run(args)
    elif args.subcommand == "sync-compat":
        return cmd_sync_compat(args)
    else:
        parser.print_help(sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
