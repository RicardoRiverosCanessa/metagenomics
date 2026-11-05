#!/usr/bin/env python3
"""write_disk_space_report.py

Structured disk space report manager for the metagenomics pipeline.
Maintains OUTDIR/logs/disk_space.current.tsv and OUTDIR/logs/disk_space.tsv.

Guarantees:
- Strict 22-column schema with invariant header and column ordering
- Strict validation of enums, integer byte fields, and sequence numbers
- Monotonic, strictly increasing measurement sequence per RUN_ID starting at 1
- Single RUN_ID validation across the entire file
- Sanitization of message strings (no tabs, carriage returns, or line feeds)
- Safe, atomic file replacement (temp file + fsync + os.replace)
- Summary status and warning count calculation based on latest check status
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

SCHEMA_COLUMNS: List[str] = [
    "run_id",
    "measurement_sequence",
    "check_id",
    "timestamp_utc",
    "stage_id",
    "storage_role",
    "purpose",
    "path",
    "filesystem_id",
    "measurement_type",
    "estimate_confidence",
    "input_bytes_read",
    "estimated_required_bytes",
    "filesystem_available_bytes",
    "filesystem_used_bytes",
    "filesystem_total_bytes",
    "required_free_inodes",
    "available_free_inodes",
    "safety_factor",
    "shared_capacity_group",
    "status",
    "message",
]

ALLOWED_STORAGE_ROLES: Set[str] = {"outdir", "tmpdir", "other"}

ALLOWED_MEASUREMENT_TYPES: Set[str] = {
    "preflight-estimate",
    "stage-precheck",
    "post-stage-observation",
    "final-observation",
}

ALLOWED_CONFIDENCE: Set[str] = {
    "hard-minimum",
    "conservative",
    "informational",
    "unknown",
    "not-applicable",
}

ALLOWED_STATUS: Set[str] = {
    "sufficient",
    "warning",
    "insufficient",
    "unknown",
    "not-applicable",
}


def sanitize_cell(val: Any) -> str:
    """Sanitize string for TSV cell safety: replace newlines/tabs with spaces/semicolons."""
    if val is None:
        return "NA"
    s = str(val).replace("\r\n", " ").replace("\n", " ").replace("\r", " ").replace("\t", " ")
    s = s.strip()
    return s if s else "NA"


def sanitize_message(val: Any) -> str:
    """Sanitize message field: remove tabs, newlines, trim, and collapse spaces."""
    if val is None:
        return ""
    s = str(val).replace("\r\n", "; ").replace("\n", "; ").replace("\r", "; ").replace("\t", " ")
    s = re.sub(r"[ \t]+", " ", s).strip()
    return s


def validate_integer_or_na(val: str, col_name: str) -> str:
    """Validate non-negative integer string or NA."""
    val = val.strip()
    if val == "NA":
        return "NA"
    if not re.match(r"^[0-9]+$", val):
        raise ValueError(f"Column '{col_name}' must be a non-negative integer or 'NA', got: '{val}'")
    return str(int(val))


def validate_float_or_na(val: str, col_name: str) -> str:
    """Validate positive float decimal or NA."""
    val = val.strip()
    if val == "NA":
        return "NA"
    try:
        f = float(val)
        if f <= 0.0:
            raise ValueError()
        return f"{f:.2f}"
    except Exception:
        raise ValueError(f"Column '{col_name}' must be a positive decimal or 'NA', got: '{val}'")


def validate_row_dict(row: Dict[str, str], expected_run_id: Optional[str] = None) -> Dict[str, str]:
    """Validate a single dictionary against the strict schema contract."""
    validated: Dict[str, str] = {}

    # 1. run_id
    run_id = sanitize_cell(row.get("run_id", ""))
    if not run_id or run_id == "NA":
        raise ValueError("run_id cannot be empty or 'NA'")
    if expected_run_id and run_id != expected_run_id:
        raise ValueError(f"run_id mismatch: expected '{expected_run_id}', found '{run_id}'")
    validated["run_id"] = run_id

    # 2. measurement_sequence
    seq_str = str(row.get("measurement_sequence", "")).strip()
    if not seq_str.isdigit() or int(seq_str) < 1:
        raise ValueError(f"measurement_sequence must be an integer >= 1, got '{seq_str}'")
    validated["measurement_sequence"] = str(int(seq_str))

    # 3. check_id
    check_id = sanitize_cell(row.get("check_id", ""))
    if not check_id or check_id == "NA":
        raise ValueError("check_id cannot be empty or 'NA'")
    validated["check_id"] = check_id

    # 4. timestamp_utc
    ts = sanitize_cell(row.get("timestamp_utc", ""))
    if not ts or ts == "NA":
        raise ValueError("timestamp_utc cannot be empty")
    validated["timestamp_utc"] = ts

    # 5. stage_id
    stage_id = sanitize_cell(row.get("stage_id", ""))
    if not stage_id or stage_id == "NA":
        raise ValueError("stage_id cannot be empty")
    validated["stage_id"] = stage_id

    # 6. storage_role
    role = sanitize_cell(row.get("storage_role", "")).lower()
    if role not in ALLOWED_STORAGE_ROLES:
        raise ValueError(f"storage_role '{role}' not in allowed: {sorted(ALLOWED_STORAGE_ROLES)}")
    validated["storage_role"] = role

    # 7. purpose
    purpose = sanitize_cell(row.get("purpose", ""))
    if not purpose or purpose == "NA":
        raise ValueError("purpose cannot be empty")
    validated["purpose"] = purpose

    # 8. path
    path = sanitize_cell(row.get("path", ""))
    if not path or path == "NA":
        raise ValueError("path cannot be empty")
    validated["path"] = path

    # 9. filesystem_id
    fs_id = sanitize_cell(row.get("filesystem_id", ""))
    if not fs_id:
        fs_id = "unknown"
    validated["filesystem_id"] = fs_id

    # 10. measurement_type
    mtype = sanitize_cell(row.get("measurement_type", "")).lower()
    if mtype not in ALLOWED_MEASUREMENT_TYPES:
        raise ValueError(f"measurement_type '{mtype}' not in allowed: {sorted(ALLOWED_MEASUREMENT_TYPES)}")
    validated["measurement_type"] = mtype

    # 11. estimate_confidence
    conf = sanitize_cell(row.get("estimate_confidence", "")).lower()
    if conf not in ALLOWED_CONFIDENCE:
        raise ValueError(f"estimate_confidence '{conf}' not in allowed: {sorted(ALLOWED_CONFIDENCE)}")
    validated["estimate_confidence"] = conf

    # 12-16. Byte sizes
    validated["input_bytes_read"] = validate_integer_or_na(str(row.get("input_bytes_read", "NA")), "input_bytes_read")
    validated["estimated_required_bytes"] = validate_integer_or_na(str(row.get("estimated_required_bytes", "NA")), "estimated_required_bytes")
    validated["filesystem_available_bytes"] = validate_integer_or_na(str(row.get("filesystem_available_bytes", "NA")), "filesystem_available_bytes")
    validated["filesystem_used_bytes"] = validate_integer_or_na(str(row.get("filesystem_used_bytes", "NA")), "filesystem_used_bytes")
    validated["filesystem_total_bytes"] = validate_integer_or_na(str(row.get("filesystem_total_bytes", "NA")), "filesystem_total_bytes")

    # 17-18. Inodes
    validated["required_free_inodes"] = validate_integer_or_na(str(row.get("required_free_inodes", "NA")), "required_free_inodes")
    validated["available_free_inodes"] = validate_integer_or_na(str(row.get("available_free_inodes", "NA")), "available_free_inodes")

    # 19. safety_factor
    validated["safety_factor"] = validate_float_or_na(str(row.get("safety_factor", "NA")), "safety_factor")

    # 20. shared_capacity_group
    scg = sanitize_cell(row.get("shared_capacity_group", "none"))
    validated["shared_capacity_group"] = scg if scg else "none"

    # 21. status
    st = sanitize_cell(row.get("status", "")).lower()
    if st not in ALLOWED_STATUS:
        raise ValueError(f"status '{st}' not in allowed: {sorted(ALLOWED_STATUS)}")
    validated["status"] = st

    # 22. message
    msg = sanitize_message(row.get("message", ""))
    validated["message"] = msg

    return validated


def load_report_tsv(path: str, expected_run_id: Optional[str] = None) -> List[Dict[str, str]]:
    """Load existing disk space TSV report. Returns empty list if missing or empty.
    Raises ValueError if corrupted or if rows do not match schema/expected_run_id."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []

    rows: List[Dict[str, str]] = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return []

        if header != SCHEMA_COLUMNS:
            raise ValueError(f"Header in '{path}' does not match required schema. Found: {header}")

        for line_num, r in enumerate(reader, start=2):
            if not r or not any(r):
                continue
            if len(r) != len(SCHEMA_COLUMNS):
                raise ValueError(f"Malformed row at line {line_num} in '{path}': expected {len(SCHEMA_COLUMNS)} columns, got {len(r)}")
            row_dict = dict(zip(SCHEMA_COLUMNS, r))
            validated = validate_row_dict(row_dict, expected_run_id=expected_run_id)
            rows.append(validated)

    return rows


def atomic_write_tsv(dest_path: str, rows: List[Dict[str, str]]) -> None:
    """Atomically write rows to dest_path using a temporary file in the same directory."""
    parent_dir = os.path.dirname(os.path.abspath(dest_path))
    if parent_dir and not os.path.exists(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)

    tmp_path = f"{dest_path}.tmp.{os.getpid()}"
    try:
        with open(tmp_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=SCHEMA_COLUMNS, delimiter="\t")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
            f.flush()
            os.fsync(f.fileno())

        os.replace(tmp_path, dest_path)
    finally:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass


def validate_and_sort_run_rows(rows: List[Dict[str, str]], run_id: str) -> List[Dict[str, str]]:
    """Validate unique sequences, sort by sequence strictly, and check uniqueness keys."""
    if not rows:
        return []

    # Sort strictly by integer measurement_sequence
    sorted_rows = sorted(rows, key=lambda x: int(x["measurement_sequence"]))

    seen_seqs: Set[int] = set()
    seen_unique_keys: Set[Tuple[str, str, str, str]] = set()

    for r in sorted_rows:
        cur_run_id = r["run_id"]
        if cur_run_id != run_id:
            raise ValueError(f"Multiple RUN_IDs detected in single report: '{cur_run_id}' vs '{run_id}'")

        seq = int(r["measurement_sequence"])
        if seq in seen_seqs:
            raise ValueError(f"Duplicate measurement_sequence '{seq}' found in report for run_id '{run_id}'")
        seen_seqs.add(seq)

        ukey = (cur_run_id, r["check_id"], r["measurement_type"], r["measurement_sequence"])
        if ukey in seen_unique_keys:
            raise ValueError(f"Duplicate unique key {ukey} found in report")
        seen_unique_keys.add(ukey)

    # Check that sequence is strictly increasing
    last_seq = 0
    for r in sorted_rows:
        seq = int(r["measurement_sequence"])
        if seq <= last_seq:
            raise ValueError(f"measurement_sequence is not strictly increasing: {seq} after {last_seq}")
        last_seq = seq

    return sorted_rows


def cmd_record(args: argparse.Namespace) -> int:
    """Record a disk measurement into disk_space.current.tsv and optionally publish."""
    report_current = os.path.abspath(args.report_current)
    run_id = sanitize_cell(args.run_id)
    if not run_id or run_id == "NA":
        print("[ERROR] --run-id is required and cannot be empty", file=sys.stderr)
        return 1

    # Load existing rows
    existing_rows: List[Dict[str, str]] = []
    if os.path.exists(report_current) and os.path.getsize(report_current) > 0:
        try:
            # If the file belongs to another run_id, we start a fresh report for the new run_id
            first_rows = load_report_tsv(report_current, expected_run_id=None)
            if first_rows and first_rows[0]["run_id"] == run_id:
                existing_rows = first_rows
            else:
                # Prior run_id report — do not keep previous run's rows
                existing_rows = []
        except Exception as e:
            print(f"[WARN] Failed to load existing report '{report_current}': {e}. Starting fresh for {run_id}.", file=sys.stderr)
            existing_rows = []

    # Determine sequence number
    if args.measurement_sequence == "auto" or not args.measurement_sequence:
        next_seq = len(existing_rows) + 1
    else:
        try:
            next_seq = int(args.measurement_sequence)
        except ValueError:
            print(f"[ERROR] Invalid --measurement-sequence '{args.measurement_sequence}'", file=sys.stderr)
            return 1

    row_data: Dict[str, str] = {
        "run_id": run_id,
        "measurement_sequence": str(next_seq),
        "check_id": args.check_id,
        "timestamp_utc": args.timestamp_utc,
        "stage_id": args.stage_id,
        "storage_role": args.storage_role,
        "purpose": args.purpose,
        "path": args.path,
        "filesystem_id": args.filesystem_id,
        "measurement_type": args.measurement_type,
        "estimate_confidence": args.estimate_confidence,
        "input_bytes_read": args.input_bytes_read,
        "estimated_required_bytes": args.estimated_required_bytes,
        "filesystem_available_bytes": args.filesystem_available_bytes,
        "filesystem_used_bytes": args.filesystem_used_bytes,
        "filesystem_total_bytes": args.filesystem_total_bytes,
        "required_free_inodes": args.required_free_inodes,
        "available_free_inodes": args.available_free_inodes,
        "safety_factor": args.safety_factor,
        "shared_capacity_group": args.shared_capacity_group,
        "status": args.status,
        "message": args.message,
    }

    try:
        validated_row = validate_row_dict(row_data, expected_run_id=run_id)
    except Exception as e:
        print(f"[ERROR] Validation failed for new row: {e}", file=sys.stderr)
        return 1

    existing_rows.append(validated_row)

    try:
        final_rows = validate_and_sort_run_rows(existing_rows, run_id=run_id)
    except Exception as e:
        print(f"[ERROR] Validation failed for cumulative report: {e}", file=sys.stderr)
        return 1

    # Atomically write to current.tsv
    try:
        atomic_write_tsv(report_current, final_rows)
    except Exception as e:
        print(f"[ERROR] Failed to write report to '{report_current}': {e}", file=sys.stderr)
        return 1

    # If publish path specified, atomically sync
    if args.publish_path:
        publish_path = os.path.abspath(args.publish_path)
        try:
            atomic_write_tsv(publish_path, final_rows)
        except Exception as e:
            print(f"[ERROR] Failed to publish report to '{publish_path}': {e}", file=sys.stderr)
            return 1

    print(f"SEQUENCE={next_seq}")
    return 0


def cmd_publish(args: argparse.Namespace) -> int:
    """Publish current report to public path atomically after validating."""
    report_current = os.path.abspath(args.report_current)
    publish_path = os.path.abspath(args.publish_path)
    run_id = sanitize_cell(args.run_id)

    if not os.path.exists(report_current) or os.path.getsize(report_current) == 0:
        print(f"[ERROR] Current report does not exist or is empty: {report_current}", file=sys.stderr)
        return 1

    try:
        rows = load_report_tsv(report_current, expected_run_id=run_id if run_id and run_id != "NA" else None)
        if not rows:
            print(f"[ERROR] No valid rows found in {report_current}", file=sys.stderr)
            return 1
        cur_run_id = rows[0]["run_id"]
        sorted_rows = validate_and_sort_run_rows(rows, run_id=cur_run_id)
        atomic_write_tsv(publish_path, sorted_rows)
    except Exception as e:
        print(f"[ERROR] Publish failed: {e}", file=sys.stderr)
        return 1

    print(f"[Disk Space] Published {len(sorted_rows)} measurements to {publish_path}")
    return 0


def cmd_summary(args: argparse.Namespace) -> int:
    """Calculate warning count and overall status from latest measurement per check_id."""
    report_current = os.path.abspath(args.report_current)
    run_id = sanitize_cell(args.run_id) if args.run_id else None

    if not os.path.exists(report_current) or os.path.getsize(report_current) == 0:
        print("TOTAL_MEASUREMENTS=0")
        print("DISTINCT_CHECKS=0")
        print("WARNING_COUNT=0")
        print("INSUFFICIENT_COUNT=0")
        print("OVERALL_STATUS=not-applicable")
        return 0

    try:
        rows = load_report_tsv(report_current, expected_run_id=run_id)
        if not rows:
            print("TOTAL_MEASUREMENTS=0")
            print("DISTINCT_CHECKS=0")
            print("WARNING_COUNT=0")
            print("INSUFFICIENT_COUNT=0")
            print("OVERALL_STATUS=not-applicable")
            return 0
    except Exception as e:
        print(f"[ERROR] Failed to read report: {e}", file=sys.stderr)
        return 1

    # Map each check_id to its latest status (rows are ordered by sequence)
    latest_check_status: Dict[str, str] = {}
    for r in rows:
        cid = r["check_id"]
        st = r["status"]
        latest_check_status[cid] = st

    warn_count = sum(1 for st in latest_check_status.values() if st == "warning")
    insufficient_count = sum(1 for st in latest_check_status.values() if st == "insufficient")
    unknown_count = sum(1 for st in latest_check_status.values() if st == "unknown")
    sufficient_count = sum(1 for st in latest_check_status.values() if st == "sufficient")

    if insufficient_count > 0:
        overall_status = "insufficient"
    elif warn_count > 0:
        overall_status = "warning"
    elif sufficient_count > 0:
        overall_status = "sufficient"
    elif unknown_count > 0:
        overall_status = "unknown"
    else:
        overall_status = "not-applicable"

    if args.json:
        data = {
            "total_measurements": len(rows),
            "distinct_checks": len(latest_check_status),
            "warning_count": warn_count,
            "insufficient_count": insufficient_count,
            "overall_status": overall_status,
            "latest_check_status": latest_check_status,
        }
        print(json.dumps(data, indent=2))
    else:
        print(f"TOTAL_MEASUREMENTS={len(rows)}")
        print(f"DISTINCT_CHECKS={len(latest_check_status)}")
        print(f"WARNING_COUNT={warn_count}")
        print(f"INSUFFICIENT_COUNT={insufficient_count}")
        print(f"OVERALL_STATUS={overall_status}")

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Structured disk space report manager.")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # Subcommand: record
    p_rec = subparsers.add_parser("record", help="Record a measurement")
    p_rec.add_argument("--report-current", required=True, help="Path to disk_space.current.tsv")
    p_rec.add_argument("--publish-path", default="", help="Path to disk_space.tsv")
    p_rec.add_argument("--run-id", required=True)
    p_rec.add_argument("--measurement-sequence", default="auto")
    p_rec.add_argument("--check-id", required=True)
    p_rec.add_argument("--timestamp-utc", required=True)
    p_rec.add_argument("--stage-id", required=True)
    p_rec.add_argument("--storage-role", required=True)
    p_rec.add_argument("--purpose", required=True)
    p_rec.add_argument("--path", required=True)
    p_rec.add_argument("--filesystem-id", default="unknown")
    p_rec.add_argument("--measurement-type", required=True)
    p_rec.add_argument("--estimate-confidence", required=True)
    p_rec.add_argument("--input-bytes-read", default="NA")
    p_rec.add_argument("--estimated-required-bytes", default="NA")
    p_rec.add_argument("--filesystem-available-bytes", default="NA")
    p_rec.add_argument("--filesystem-used-bytes", default="NA")
    p_rec.add_argument("--filesystem-total-bytes", default="NA")
    p_rec.add_argument("--required-free-inodes", default="NA")
    p_rec.add_argument("--available-free-inodes", default="NA")
    p_rec.add_argument("--safety-factor", default="NA")
    p_rec.add_argument("--shared-capacity-group", default="none")
    p_rec.add_argument("--status", required=True)
    p_rec.add_argument("--message", default="")

    # Subcommand: publish
    p_pub = subparsers.add_parser("publish", help="Publish current report to public path")
    p_pub.add_argument("--report-current", required=True)
    p_pub.add_argument("--publish-path", required=True)
    p_pub.add_argument("--run-id", default="")

    # Subcommand: summary
    p_sum = subparsers.add_parser("summary", help="Summary status and warning counts")
    p_sum.add_argument("--report-current", required=True)
    p_sum.add_argument("--run-id", default="")
    p_sum.add_argument("--json", action="store_true", default=False)

    args = parser.parse_args()

    if args.subcommand == "record":
        return cmd_record(args)
    elif args.subcommand == "publish":
        return cmd_publish(args)
    elif args.subcommand == "summary":
        return cmd_summary(args)
    else:
        parser.print_help()
        return 2


if __name__ == "__main__":
    sys.exit(main())
