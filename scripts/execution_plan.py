#!/usr/bin/env python3
"""scripts/execution_plan.py - Canonical Execution Plan Validator, Summarizer, and Renderer.

MetaNexus 0.3.2 (development)
Adheres strictly to the execution plan contract:
- Validates schema, closed enums, types, and lists.
- Validates graph consistency and detects cycles (exits with code 3 on invariant violation).
- Calculates summary counts and executability.
- Renders deterministic JSON and plain-text outputs.
- DOES NOT decide stage applicability, actions, force, restart, compatibility, or invalidation.
"""

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Set

SCHEMA_VERSION = "1"
PIPELINE_NAME = "MetaNexus"
PIPELINE_VERSION = "0.3.2"
PIPELINE_RELEASE_STATUS = "development"

STAGE_IDS: List[str] = [
    "01_fastp",
    "03_kraken2",
    "02_megahit",
    "04_mapping",
    "05a_vamb",
    "05b_metabinner",
    "05c_comebin",
    "06_refinement",
    "07_reassembly",
    "08_checkm2",
    "09_filter_mags",
    "10_gtdbtk",
    "11_coverm",
    "12a_taxonomy_summary",
    "12b_mag_summary",
    "12c_run_overview",
]

# Closed enumerations
CHECKPOINT_STATUSES: Set[str] = {
    "missing",
    "valid",
    "invalid",
    "stale",
    "legacy",
    "not-applicable",
    "not-inspected",
}

RESOURCE_STATUSES: Set[str] = {
    "available",
    "unknown",
    "unavailable",
    "invalid",
    "not-configured",
    "not-applicable",
    "indeterminate",
}

OUTPUT_STATUSES: Set[str] = {
    "available",
    "missing",
    "invalid",
    "unknown",
    "not-applicable",
    "not-inspected",
}

VALIDATION_STATUSES: Set[str] = {
    "valid",
    "valid-with-warnings",
    "invalid",
    "failed",
    "not-performed",
    "not-applicable",
    "not-inspected",
}

PLANNED_ACTIONS: Set[str] = {
    "run",
    "reuse",
    "force",
    "restart",
    "disabled",
    "not-applicable",
    "blocked",
}

REQUESTED_ACTIONS: Set[str] = {
    "normal",
    "force",
    "restart",
}

APPLICABILITY_VALUES: Set[str] = {
    "applicable",
    "not-applicable",
    "indeterminate",
}

KNOWN_REASON_CODES: Set[str] = {
    "checkpoint-valid",
    "checkpoint-missing",
    "checkpoint-invalid",
    "checkpoint-stale",
    "checkpoint-legacy",
    "checkpoint-not-inspected",
    "configuration-changed",
    "input-fingerprint-changed",
    "runtime-fingerprint-changed",
    "database-fingerprint-changed",
    "product-missing",
    "product-validation-failed",
    "upstream-invalidated",
    "forced-by-user",
    "upstream-forced",
    "restart-boundary",
    "downstream-of-restart",
    "stage-disabled",
    "stage-not-applicable",
    "external-input-missing",
    "external-input-invalid",
    "runtime-unavailable",
    "database-not-configured",
    "dependency-blocked",
    "applicability-indeterminate",
    "output-parent-not-writable",
    "output-writability-unknown",
}


# Exit codes
EXIT_OK = 0
EXIT_BLOCKED_OR_STRICT = 1
EXIT_USAGE_ERROR = 2
EXIT_INTERNAL_ERROR = 3


class ExecutionPlanError(Exception):
    """Base exception for execution plan errors."""
    pass


class ExecutionPlanValidationError(ExecutionPlanError):
    """Validation or configuration error (exit 2)."""
    pass


class ExecutionPlanInternalError(ExecutionPlanError):
    """Internal implementation or invariant failure (exit 3)."""
    pass


def fail_validation(message: str) -> None:
    """Raise validation error (exit code 2)."""
    raise ExecutionPlanValidationError(message)


def fail_invariant(message: str) -> None:
    """Raise invariant failure (exit code 3)."""
    raise ExecutionPlanInternalError(message)


ALL_STAGE_FIELDS: List[str] = [
    "stage_id",
    "stage_name",
    "configured_enabled",
    "required",
    "requested_action",
    "applicability",
    "planned_action",
    "checkpoint_path",
    "checkpoint_status",
    "checkpoint_reason",
    "input_status",
    "runtime_status",
    "database_status",
    "output_status",
    "validation_status",
    "invalidated_by",
    "blocking_stage",
    "reason_code",
    "reason_detail",
    "diagnostic_codes",
    "dependencies",
    "downstream",
    "primary_output",
    "primary_log",
]

CANONICAL_STAGE_IDS = STAGE_IDS
VALID_PLANNED_ACTIONS = PLANNED_ACTIONS
VALID_REQUESTED_ACTIONS = REQUESTED_ACTIONS
VALID_APPLICABILITY = APPLICABILITY_VALUES
VALID_CHECKPOINT_STATUS = CHECKPOINT_STATUSES
VALID_INPUT_STATUS = RESOURCE_STATUSES
VALID_RUNTIME_STATUS = RESOURCE_STATUSES
VALID_DATABASE_STATUS = RESOURCE_STATUSES
VALID_OUTPUT_STATUS = OUTPUT_STATUSES
VALID_VALIDATION_STATUS = VALIDATION_STATUSES


def is_ordered_and_unique(items: List[str]) -> bool:
    """Check that a list contains no duplicates."""
    return len(items) == len(set(items))


def detect_cycles(plan: List[Dict[str, Any]]) -> bool:
    """Detect cycles in dependency graph using DFS. Returns True if cycle found."""
    graph: Dict[str, List[str]] = {}
    for stg in plan:
        sid = stg.get("stage_id", "")
        deps = stg.get("dependencies", [])
        graph[sid] = list(deps)

    # 0 = unvisited, 1 = visiting (in recursion stack), 2 = visited
    state: Dict[str, int] = {k: 0 for k in graph}

    def dfs(node: str) -> bool:
        state[node] = 1
        for parent in graph.get(node, []):
            if parent not in state:
                continue
            if state[parent] == 1:
                return True
            if state[parent] == 0:
                if dfs(parent):
                    return True
        state[node] = 2
        return False

    for node in graph:
        if state[node] == 0:
            if dfs(node):
                return True
    return False


def validate_stage_record(stg: Dict[str, Any], seen_ids: Set[str]) -> None:
    """Validate a single canonical stage record."""
    for field in ALL_STAGE_FIELDS:
        if field not in stg:
            fail_validation(f"Missing required field '{field}' in stage record")

    sid = stg["stage_id"]
    if not isinstance(sid, str) or not sid:
        fail_invariant(f"stage_id must be non-empty string, got {sid!r}")
    if sid in seen_ids:
        fail_invariant(f"Duplicated stage_id '{sid}'")
    seen_ids.add(sid)

    if not isinstance(stg["stage_name"], str) or not stg["stage_name"]:
        fail_invariant(f"stage_name must be non-empty string for stage {sid}")

    if not isinstance(stg["configured_enabled"], bool):
        fail_invariant(f"configured_enabled must be boolean for stage {sid}")

    if not isinstance(stg["required"], bool):
        fail_invariant(f"required must be boolean for stage {sid}")

    req_act = stg["requested_action"]
    if req_act not in REQUESTED_ACTIONS:
        fail_validation(f"Invalid requested_action '{req_act}' for stage {sid}")

    appl = stg["applicability"]
    if appl not in APPLICABILITY_VALUES:
        fail_validation(f"Invalid applicability '{appl}' for stage {sid}")

    plan_act = stg["planned_action"]
    if plan_act not in PLANNED_ACTIONS:
        fail_validation(f"Invalid planned_action '{plan_act}' for stage {sid}")

    ckpt_st = stg["checkpoint_status"]
    if ckpt_st not in CHECKPOINT_STATUSES:
        fail_validation(f"Invalid checkpoint_status '{ckpt_st}' for stage {sid}")

    for res_name in ["input_status", "runtime_status", "database_status"]:
        val = stg[res_name]
        if val not in RESOURCE_STATUSES:
            fail_validation(f"Invalid {res_name} '{val}' for stage {sid}")

    out_st = stg["output_status"]
    if out_st not in OUTPUT_STATUSES:
        fail_validation(f"Invalid output_status '{out_st}' for stage {sid}")

    val_st = stg["validation_status"]
    if val_st not in VALIDATION_STATUSES:
        fail_validation(f"Invalid validation_status '{val_st}' for stage {sid}")

    for list_field in ["invalidated_by", "diagnostic_codes", "dependencies", "downstream"]:
        val = stg[list_field]
        if not isinstance(val, list):
            fail_invariant(f"{list_field} must be a list in stage {sid}")
        if not is_ordered_and_unique(val):
            fail_invariant(f"{list_field} contains duplicates in stage {sid}")

    # Specific invariants
    if plan_act == "disabled" and stg["configured_enabled"] is not False:
        fail_invariant(f"Stage {sid} has planned_action=disabled but configured_enabled is not false")

    if appl == "not-applicable" and stg["required"] is True:
        fail_invariant(f"Stage {sid} has applicability=not-applicable but required is true")

    if plan_act == "force" and req_act != "force":
        fail_invariant(f"Stage {sid} has planned_action=force but requested_action is '{req_act}'")

    if plan_act == "restart" and req_act != "restart":
        fail_invariant(f"Stage {sid} has planned_action=restart but requested_action is '{req_act}'")

    if plan_act == "reuse" and ckpt_st != "valid":
        fail_invariant(f"Stage {sid} has planned_action=reuse but checkpoint_status is '{ckpt_st}'")


def validate_plan_invariants(plan: List[Dict[str, Any]]) -> None:
    """Validate cross-stage graph, dependency, and invalidation invariants."""
    stage_by_id: Dict[str, Dict[str, Any]] = {s["stage_id"]: s for s in plan}

    # Verify dependency existence
    for stg in plan:
        sid = stg["stage_id"]
        for dep in stg["dependencies"]:
            if dep not in stage_by_id:
                fail_invariant(f"Stage {sid} depends on unknown stage '{dep}'")
        for down in stg["downstream"]:
            if down not in stage_by_id:
                fail_invariant(f"Stage {sid} references unknown downstream stage '{down}'")

    # Verify acyclicity
    if detect_cycles(plan):
        fail_invariant("Cycle detected in static dependency graph")

    # Invariants on reuse vs rebuilt dependencies
    for stg in plan:
        sid = stg["stage_id"]
        if stg["planned_action"] == "reuse":
            for dep_id in stg["dependencies"]:
                dep_act = stage_by_id[dep_id]["planned_action"]
                if dep_act in {"run", "force", "restart"}:
                    fail_invariant(
                        f"Stage {sid} has planned_action=reuse but dependency '{dep_id}' is {dep_act}"
                    )

    # Invariants on restart descendants
    for stg in plan:
        sid = stg["stage_id"]
        if "restart-boundary" in stg["diagnostic_codes"] or stg["planned_action"] == "restart":
            # Descendants must not be restart
            for down_id in stg["downstream"]:
                down_act = stage_by_id[down_id]["planned_action"]
                if down_act == "restart":
                    fail_invariant(
                        f"Restart descendant '{down_id}' has planned_action=restart (must be run)"
                    )


def compute_summary(plan: List[Dict[str, Any]]) -> Dict[str, int]:
    """Calculate summary counts for all planned actions."""
    counts = {
        "run": 0,
        "reuse": 0,
        "force": 0,
        "restart": 0,
        "disabled": 0,
        "not_applicable": 0,
        "blocked": 0,
    }
    for stg in plan:
        action = stg.get("planned_action", "")
        if action == "run":
            counts["run"] += 1
        elif action == "reuse":
            counts["reuse"] += 1
        elif action == "force":
            counts["force"] += 1
        elif action == "restart":
            counts["restart"] += 1
        elif action == "disabled":
            counts["disabled"] += 1
        elif action == "not-applicable":
            counts["not_applicable"] += 1
        elif action == "blocked":
            counts["blocked"] += 1
        else:
            fail_invariant(f"Unknown planned_action '{action}' in summary count")
    return counts


def is_executable(plan: List[Dict[str, Any]]) -> bool:
    """Plan is executable if and only if no stage has required=true and planned_action=blocked."""
    for stg in plan:
        if stg.get("required") is True and stg.get("planned_action") == "blocked":
            return False
    return True


compute_plan_summary = compute_summary
is_plan_executable = is_executable


def validate_execution_plan(plan_data: Dict[str, Any]) -> None:
    """Validate full plan document against schema and invariants."""
    if not isinstance(plan_data, dict):
        fail_validation("Top-level execution plan must be a JSON object")

    if "schema_version" not in plan_data:
        fail_validation("Missing required field 'schema_version'")
    ver = plan_data.get("schema_version")
    if ver != SCHEMA_VERSION:
        fail_invariant(f"Expected schema_version '{SCHEMA_VERSION}', got '{ver}'")

    if "pipeline" not in plan_data:
        fail_validation("Missing required field 'pipeline'")
    pipeline = plan_data.get("pipeline", {})
    if not isinstance(pipeline, dict):
        fail_validation("pipeline must be an object")

    pname = pipeline.get("name")
    if pname != PIPELINE_NAME:
        fail_invariant(f"Expected pipeline name '{PIPELINE_NAME}', got '{pname}'")

    pver = pipeline.get("version")
    if pver != PIPELINE_VERSION:
        fail_invariant(f"Expected pipeline version '{PIPELINE_VERSION}', got '{pver}'")

    pstatus = pipeline.get("release_status")
    if pstatus != PIPELINE_RELEASE_STATUS:
        fail_invariant(f"Expected pipeline release_status '{PIPELINE_RELEASE_STATUS}', got '{pstatus}'")

    if "plan" not in plan_data:
        fail_validation("Missing required field 'plan'")
    plan = plan_data.get("plan")
    if not isinstance(plan, list):
        fail_validation("plan must be a list of stage records")

    seen_ids: Set[str] = set()
    for stg in plan:
        if not isinstance(stg, dict):
            fail_invariant("Each item in plan must be an object")
        validate_stage_record(stg, seen_ids)

    validate_plan_invariants(plan)

    # Validate summary consistency if present
    calc_summary = compute_summary(plan)
    if "summary" in plan_data:
        spec_summary = plan_data["summary"]
        for k, v in calc_summary.items():
            if spec_summary.get(k) != v:
                fail_invariant(f"Summary mismatch for '{k}': recorded {spec_summary.get(k)}, computed {v}")

    # Validate executable consistency if present
    calc_exec = is_executable(plan)
    if "executable" in plan_data:
        spec_exec = plan_data["executable"]
        if spec_exec is not calc_exec:
            fail_invariant(f"Executable mismatch: recorded {spec_exec}, computed {calc_exec}")


def canonicalize_plan_data(plan_data: Dict[str, Any]) -> Dict[str, Any]:
    """Ensure list fields are ordered and unique, and summary and executable fields are present."""
    plan = plan_data.get("plan", [])
    for stg in plan:
        for list_fld in ["invalidated_by", "diagnostic_codes", "dependencies", "downstream"]:
            if list_fld in stg and isinstance(stg[list_fld], list):
                seen = set()
                dedup = []
                for x in stg[list_fld]:
                    if x not in seen:
                        seen.add(x)
                        dedup.append(x)
                stg[list_fld] = dedup

    summary = compute_summary(plan)
    executable = is_executable(plan)
    plan_data["schema_version"] = SCHEMA_VERSION
    plan_data["summary"] = summary
    plan_data["executable"] = executable
    return plan_data
    return plan_data


def assemble_plan_from_tsv(content: str, profile: str) -> Dict[str, Any]:
    """Assemble execution plan JSON structure from TSV rows produced by collect_execution_plan."""
    plan = []
    for line in content.splitlines():
        line = line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 24:
            continue
        (
            sid, sname, conf_en, req, req_act, appl, plan_act,
            chk_path, chk_st, chk_rs, in_st, run_st, db_st,
            out_st, val_st, inval_by, blk_stg, r_code, r_detail,
            diag_codes, deps, down, prim_out, prim_log
        ) = fields[:24]

        def parse_list(val: str) -> List[str]:
            if not val or val == "NA":
                return []
            parts = [p.strip() for p in val.replace(",", " ").split() if p.strip()]
            seen = set()
            out = []
            for p in parts:
                if p not in seen:
                    seen.add(p)
                    out.append(p)
            return out

        def parse_str(val: str) -> Optional[str]:
            if not val or val == "NA":
                return None
            return val

        record = {
            "stage_id": sid,
            "stage_name": sname,
            "configured_enabled": conf_en.lower() == "true",
            "required": req.lower() == "true",
            "requested_action": req_act,
            "applicability": appl,
            "planned_action": plan_act,
            "checkpoint_path": chk_path,
            "checkpoint_status": chk_st,
            "checkpoint_reason": parse_str(chk_rs),
            "input_status": in_st,
            "runtime_status": run_st,
            "database_status": db_st,
            "output_status": out_st,
            "validation_status": val_st,
            "invalidated_by": parse_list(inval_by),
            "blocking_stage": parse_str(blk_stg),
            "reason_code": r_code,
            "reason_detail": parse_str(r_detail),
            "diagnostic_codes": parse_list(diag_codes),
            "dependencies": parse_list(deps),
            "downstream": parse_list(down),
            "primary_output": prim_out,
            "primary_log": prim_log,
        }
        plan.append(record)

    doc = {
        "schema_version": "1",
        "pipeline": {
            "name": "MetaNexus",
            "version": "0.3.2",
            "release_status": "development",
            "profile": profile,
        },
        "plan": plan,
    }
    return canonicalize_plan_data(doc)


def render_json(plan_data: Dict[str, Any]) -> str:
    """Render canonical JSON execution plan."""
    validate_execution_plan(plan_data)
    canonical = canonicalize_plan_data(plan_data)
    return json.dumps(canonical, indent=2)


def render_text(plan_data: Dict[str, Any], details: str = "full") -> str:
    """Render plain-text execution plan (summary or full)."""
    validate_execution_plan(plan_data)
    canonical = canonicalize_plan_data(plan_data)
    plan = canonical.get("plan", [])
    summary = canonical.get("summary", {})
    exec_val = canonical.get("executable", True)

    lines: List[str] = []
    lines.append(f"{PIPELINE_NAME} {PIPELINE_VERSION} execution plan")
    lines.append("")

    headers = [
        "STAGE",
        "REQUESTED",
        "ACTION",
        "REQUIRED",
        "APPLICABILITY",
        "CHECKPOINT",
        "REASON",
        "INVALIDATED_BY",
    ]

    rows: List[List[str]] = []
    for stg in plan:
        inval = stg.get("invalidated_by", [])
        inval_str = ",".join(inval) if inval else "NA"
        rows.append([
            str(stg.get("stage_id", "NA")),
            str(stg.get("requested_action", "NA")),
            str(stg.get("planned_action", "NA")),
            "true" if stg.get("required") is True else "false",
            str(stg.get("applicability", "NA")),
            str(stg.get("checkpoint_status", "NA")),
            str(stg.get("reason_code", "NA")),
            inval_str,
        ])

    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    header_line = "  ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    lines.append(header_line)
    lines.append("  ".join("-" * col_widths[i] for i in range(len(headers))))

    for row in rows:
        lines.append("  ".join(val.ljust(col_widths[i]) for i, val in enumerate(row)))

    lines.append("")
    lines.append("Summary:")
    lines.append(f"  Run:            {summary.get('run', 0)}")
    lines.append(f"  Reuse:          {summary.get('reuse', 0)}")
    lines.append(f"  Force:          {summary.get('force', 0)}")
    lines.append(f"  Restart:        {summary.get('restart', 0)}")
    lines.append(f"  Disabled:       {summary.get('disabled', 0)}")
    lines.append(f"  Not applicable: {summary.get('not_applicable', 0)}")
    lines.append(f"  Blocked:        {summary.get('blocked', 0)}")
    lines.append(f"  Executable:     {'true' if exec_val else 'false'}")

    if details == "full":
        lines.append("")
        lines.append("Stage details:")
        for stg in plan:
            sid = stg.get("stage_id", "NA")
            lines.append(f"  {sid}:")
            deps = stg.get("dependencies", [])
            lines.append(f"    dependencies:       {json.dumps(deps)}")
            down = stg.get("downstream", [])
            lines.append(f"    downstream:         {json.dumps(down)}")
            lines.append(f"    input_status:       {stg.get('input_status', 'NA')}")
            lines.append(f"    runtime_status:     {stg.get('runtime_status', 'NA')}")
            lines.append(f"    database_status:    {stg.get('database_status', 'NA')}")
            lines.append(f"    output_status:      {stg.get('output_status', 'NA')}")
            lines.append(f"    validation_status:  {stg.get('validation_status', 'NA')}")
            lines.append(f"    checkpoint_path:    {stg.get('checkpoint_path') or 'NA'}")
            lines.append(f"    primary_output:     {stg.get('primary_output') or 'NA'}")
            lines.append(f"    primary_log:        {stg.get('primary_log') or 'NA'}")
            r_detail = stg.get("reason_detail")
            lines.append(f"    reason_detail:      {r_detail if r_detail is not None else 'NA'}")
            diag = stg.get("diagnostic_codes", [])
            lines.append(f"    diagnostic_codes:   {json.dumps(diag)}")
            blk = stg.get("blocking_stage")
            lines.append(f"    blocking_stage:     {blk if blk is not None else 'NA'}")

    return "\n".join(lines)


def load_input_data(input_path: Optional[str]) -> Dict[str, Any]:
    """Read plan JSON data from file or stdin."""
    try:
        if not input_path or input_path == "-":
            content = sys.stdin.read()
        else:
            with open(input_path, "r", encoding="utf-8") as f:
                content = f.read()
        return json.loads(content)
    except json.JSONDecodeError as exc:
        sys.stderr.write(f"[ERROR] Failed to parse input plan JSON: {exc}\n")
        sys.exit(3)
    except Exception as exc:
        sys.stderr.write(f"[ERROR] Failed to read input plan: {exc}\n")
        sys.exit(3)


def main() -> None:
    """CLI entry point for scripts/execution_plan.py."""
    parser = argparse.ArgumentParser(
        description="MetaNexus Execution Plan Validator, Summarizer, and Renderer"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # validate
    p_val = subparsers.add_parser("validate", help="Validate plan invariants")
    p_val.add_argument("--input", default="-", help="Input JSON file or '-' for stdin")

    # render
    p_ren = subparsers.add_parser("render", help="Render plan in text or json")
    p_ren.add_argument("--format", choices=["text", "json"], default="text", help="Output format")
    p_ren.add_argument("--details", choices=["summary", "full"], default="full", help="Detail level")
    p_ren.add_argument("--input", default="-", help="Input JSON file or '-' for stdin")

    # is-executable
    p_exe = subparsers.add_parser("is-executable", help="Exit 0 if executable, 1 if not")
    p_exe.add_argument("--input", default="-", help="Input JSON file or '-' for stdin")

    # summary
    p_sum = subparsers.add_parser("summary", help="Print summary JSON")
    p_sum.add_argument("--input", default="-", help="Input JSON file or '-' for stdin")

    # assemble
    p_asm = subparsers.add_parser("assemble", help="Assemble JSON document from TSV stream")
    p_asm.add_argument("--profile", default="full", help="Pipeline profile")
    p_asm.add_argument("--input", default="-", help="Input TSV file or '-' for stdin")

    try:
        args = parser.parse_args()
    except SystemExit as exc:
        # User CLI error
        sys.exit(2 if exc.code != 0 else 0)

    try:
        if args.command == "assemble":
            if not args.input or args.input == "-":
                tsv_content = sys.stdin.read()
            else:
                with open(args.input, "r", encoding="utf-8") as f:
                    tsv_content = f.read()
            doc = assemble_plan_from_tsv(tsv_content, args.profile)
            sys.stdout.write(json.dumps(doc, indent=2) + "\n")
            sys.exit(0)

        data = load_input_data(args.input)
        validate_execution_plan(data)

        if args.command == "validate":
            sys.exit(0)
        elif args.command == "render":
            if args.format == "json":
                sys.stdout.write(render_json(data) + "\n")
            else:
                sys.stdout.write(render_text(data, details=args.details) + "\n")
            sys.exit(0)
        elif args.command == "is-executable":
            executable = is_executable(data.get("plan", []))
            sys.exit(0 if executable else 1)
        elif args.command == "summary":
            summary = compute_summary(data.get("plan", []))
            sys.stdout.write(json.dumps(summary, indent=2) + "\n")
            sys.exit(0)
    except ExecutionPlanValidationError as exc:
        sys.stderr.write(f"[ERROR] Plan validation error: {exc}\n")
        sys.exit(2)
    except ExecutionPlanInternalError as exc:
        sys.stderr.write(f"[ERROR] Plan invariant failure: {exc}\n")
        sys.exit(3)


if __name__ == "__main__":
    main()
