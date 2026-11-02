#!/usr/bin/env python3
"""write_error_record.py - Structured, safe, non-recursive fatal termination record manager.

Manages persistent error tracking for the metagenomics pipeline:
- Atomically tracks error state via pending/committed markers.
- Deduplicates errors between subshells, functions, and orchestrator process.
- Implements exclusive flock locking on errors.tsv.
- Atomically publishes last_error.json materialized view.
- Redacts secrets and sanitizes TSV fields.
- Extracts bounded log tails (text only, <= 50 lines, <= 8 KiB).
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import time

try:
    import fcntl
    HAS_FCNTL = True
except ImportError:
    HAS_FCNTL = False

SCHEMA_VERSION = 1

TSV_HEADER = [
    "run_id",
    "error_id",
    "timestamp_utc",
    "stage_id",
    "stage_name",
    "stage_action",
    "exit_code",
    "signal",
    "trap_source_file",
    "trap_source_line",
    "call_source_file",
    "call_source_line",
    "function_name",
    "command",
    "log_path",
    "runtime",
    "last_valid_checkpoint",
    "last_valid_checkpoint_path",
    "last_valid_product",
    "last_valid_product_fingerprint",
    "message",
    "status",
    "run_status_update_status",
]

ALLOWED_STATUSES = {
    "failed",
    "interrupted",
    "terminated",
    "preflight-failed",
    "integration-error",
}

ALLOWED_RUN_STATUS_UPDATE = {
    "updated",
    "failed",
    "not-available",
}

ALLOWED_ACTIONS = {
    "run",
    "reuse",
    "force",
    "restart",
    "preflight",
    "summary",
    "unknown",
}

ALLOWED_SIGNALS = {
    "none",
    "INT",
    "TERM",
    "HUP",
}

RECOVERY_HINTS = {
    "01_fastp": "Fix the fastp input or runtime problem and rerun with --force-stage fastp.",
    "02_megahit": "Fix MEGAHIT or input reads and rerun with --force-stage megahit.",
    "03_kraken2": "Fix Kraken2/Bracken or its database and rerun with --force-stage kraken2.",
    "04_mapping": "Fix Bowtie2/Samtools mapping and rerun with --force-stage mapping.",
    "05a_vamb": "Fix Vamb and rerun with --force-stage vamb.",
    "05b_metabinner": "Fix MetaBinner and rerun with --force-stage metabinner.",
    "05c_comebin": "Fix the COMEBin error and rerun with --force-stage comebin.",
    "06_refinement": "Fix MetaWRAP refinement and rerun with --force-stage refinement.",
    "07_reassembly": "Fix reassembly and rerun with --force-stage reassembly.",
    "08_checkm2": "Fix CheckM2 or its database and rerun with --force-stage checkm2.",
    "09_filter_mags": "Fix MAG filtering and rerun with --force-stage mag-filtering.",
    "10_gtdbtk": "Fix GTDB-Tk or its database and rerun with --force-stage gtdbtk.",
    "11_coverm": "Fix CoverM and rerun with --force-stage coverm.",
    "12a_taxonomy_summary": "Fix the summary inputs and rerun with --force-stage taxonomy-summary.",
    "12b_mag_summary": "Fix the summary inputs and rerun with --force-stage mag-summary.",
    "12c_run_overview": "Fix the summary inputs and rerun with --force-stage run-overview.",
    "pipeline": "Fix the preflight error and rerun the original command.",
}

SECRET_PATTERNS = [
    (re.compile(r'(--(?:token|password|secret|credential|api-key|key|access-token))(?:[=\s]+)(\S+)', re.IGNORECASE), r'\1=<redacted>'),
    (re.compile(r'\b(TOKEN|PASSWORD|SECRET|CREDENTIAL|API_KEY|ACCESS_TOKEN)=(\S+)', re.IGNORECASE), r'\1=<redacted>'),
]


def redact_secrets(text: str) -> str:
    """Redacts common credentials, tokens, and passwords."""
    if not text:
        return ""
    result = str(text)
    for pattern, repl in SECRET_PATTERNS:
        result = pattern.sub(repl, result)
    return result


def sanitize_tsv_cell(val: object, default: str = "unknown") -> str:
    """Sanitizes TSV cell text: strips tabs, newlines, CRs and redacts secrets."""
    if val is None:
        return default
    s = redact_secrets(str(val))
    s = s.replace("\r", " ").replace("\n", " ").replace("\t", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s if s else default


def compute_error_key(
    run_id: str,
    stage_id: str,
    exit_code: int | str,
    trap_file: str,
    trap_line: int | str,
    call_file: str,
    call_line: int | str,
    command: str,
) -> str:
    """Computes a deterministic SHA-256 error key over failure identity."""
    sanitized_cmd = sanitize_tsv_cell(command)
    payload = (
        f"{run_id}|{stage_id}|{exit_code}|"
        f"{trap_file}|{trap_line}|{call_file}|{call_line}|{sanitized_cmd}"
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def get_recovery_hint(stage_id: str, custom_hint: str | None = None) -> str:
    """Returns canonical recovery hint for stage."""
    if custom_hint and custom_hint not in ("unknown", "none"):
        return custom_hint
    if stage_id in RECOVERY_HINTS:
        return RECOVERY_HINTS[stage_id]
    if stage_id.startswith("0") or stage_id.startswith("1"):
        # Strip numeric prefix for force-stage name
        parts = stage_id.split("_", 1)
        sub_name = parts[1] if len(parts) > 1 else stage_id
        return f"Fix stage {stage_id} failure and rerun with --force-stage {sub_name}."
    return f"Fix stage {stage_id} failure and rerun the pipeline."


def extract_log_tail(log_path: str, max_lines: int = 20, max_bytes: int = 8192) -> list[str]:
    """Extracts bounded text tail from log file (max 50 lines, max 8 KiB)."""
    if not log_path or log_path in ("unknown", "none"):
        return []
    path = Path(log_path)
    if not path.is_file():
        return []
    try:
        # Check for binary file content
        with path.open("rb") as f:
            chunk = f.read(1024)
            if b"\x00" in chunk:
                return []
            f.seek(0, os.SEEK_END)
            size = f.tell()
            # Read at most 64 KiB from end to parse lines efficiently
            read_size = min(size, 65536)
            f.seek(max(0, size - read_size))
            tail_bytes = f.read()

        text = tail_bytes.decode("utf-8", errors="replace")
        lines = text.splitlines()
        capped_lines = min(max(1, max_lines), 50)
        selected = lines[-capped_lines:]

        # Cap total bytes <= max_bytes
        total_bytes = sum(len(line.encode("utf-8", errors="replace")) for line in selected)
        while selected and total_bytes > max_bytes:
            removed = selected.pop(0)
            total_bytes -= len(removed.encode("utf-8", errors="replace"))

        return [redact_secrets(line) for line in selected]
    except Exception:
        return []


def is_pid_alive(pid: int) -> bool:
    """Checks if a given process PID is currently running."""
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        # Process exists, cannot signal
        return True


def parse_marker(marker_path: Path) -> dict[str, str] | None:
    """Parses marker file into a dictionary."""
    if not marker_path.is_file():
        return None
    try:
        res = {}
        for line in marker_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            res[k.strip()] = v.strip()
        return res if res else None
    except Exception:
        return None


def write_marker_atomic(marker_path: Path, data: dict[str, str]) -> None:
    """Atomically writes marker file via tempfile in the same directory."""
    dir_path = marker_path.parent
    dir_path.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=dir_path, delete=False, prefix=".marker_tmp_", encoding="utf-8") as f:
        tmp_name = f.name
        for k, v in data.items():
            f.write(f"{k}={v}\n")
        f.flush()
        try:
            os.fsync(f.fileno())
        except OSError:
            pass
    os.replace(tmp_name, marker_path)


class ErrorRecordManager:
    """Manages the full lifecycle of persistent error recording."""

    def __init__(self, outdir: Path | str, run_id: str):
        self.outdir = Path(outdir).resolve()
        self.logs_dir = self.outdir / "logs"
        self.run_id = run_id
        self.marker_path = self.logs_dir / f".error_recorded.{self.run_id}"
        self.tsv_path = self.logs_dir / "errors.tsv"
        self.json_path = self.logs_dir / "last_error.json"

    def reserve(
        self,
        writer_pid: int,
        error_key: str,
        stage_id: str,
    ) -> tuple[bool, str, str, str]:
        """Attempts to reserve marker as pending.

        Returns: (is_primary_writer, error_id, current_state, error_key)
        """
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        now_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

        # 1. Attempt exclusive creation
        try:
            fd = os.open(
                str(self.marker_path),
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o644,
            )
            try:
                # Assign initial ERROR_ID under lock
                error_id = self._assign_next_error_id_under_lock(stage_id)
                content = (
                    f"state=pending\n"
                    f"writer_pid={writer_pid}\n"
                    f"created_at={now_utc}\n"
                    f"error_key={error_key}\n"
                    f"error_id={error_id}\n"
                )
                os.write(fd, content.encode("utf-8"))
                try:
                    os.fsync(fd)
                except OSError:
                    pass
            finally:
                os.close(fd)
            return True, error_id, "pending", error_key
        except FileExistsError:
            pass

        # 2. Marker already exists: Inspect state
        marker = parse_marker(self.marker_path)
        if marker and marker.get("state") == "committed":
            return False, marker.get("error_id", f"{self.run_id}:E0001"), "committed", marker.get("error_key", error_key)

        # 3. Marker pending: Check if PID is alive
        existing_pid = int(marker.get("writer_pid", "0")) if marker else 0
        if marker and marker.get("state") == "pending" and is_pid_alive(existing_pid):
            # Wait up to 2 seconds for writer to commit
            deadline = time.time() + 2.0
            while time.time() < deadline:
                time.sleep(0.05)
                cur = parse_marker(self.marker_path)
                if cur and cur.get("state") == "committed":
                    return False, cur.get("error_id", f"{self.run_id}:E0001"), "committed", cur.get("error_key", error_key)

            # Still pending with living PID
            cur = parse_marker(self.marker_path)
            err_id = cur.get("error_id", f"{self.run_id}:E0001") if cur else f"{self.run_id}:E0001"
            return False, err_id, "pending", error_key

        # 4. Marker pending but writer PID is dead: Recover abandoned pending
        error_id = marker.get("error_id") if marker else None
        if not error_id:
            error_id = self._assign_next_error_id_under_lock(stage_id)

        rec_data = {
            "state": "pending",
            "writer_pid": str(writer_pid),
            "created_at": now_utc,
            "error_key": error_key,
            "error_id": error_id,
        }
        write_marker_atomic(self.marker_path, rec_data)
        return True, error_id, "pending", error_key

    def _assign_next_error_id_under_lock(self, stage_id: str) -> str:
        """Determines next ERROR_ID for this RUN_ID under lock."""
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        with open(self.tsv_path, "a+", encoding="utf-8") as f:
            if HAS_FCNTL:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                f.seek(0)
                lines = f.readlines()
                max_seq = 0
                for line in lines:
                    cols = line.split("\t")
                    if len(cols) > 1 and cols[0] == self.run_id:
                        eid = cols[1]
                        if ":E" in eid:
                            try:
                                seq = int(eid.split(":E")[-1])
                                if seq > max_seq:
                                    max_seq = seq
                            except ValueError:
                                pass
                next_seq = max_seq + 1
                return f"{self.run_id}:E{next_seq:04d}"
            finally:
                if HAS_FCNTL:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

    def commit(
        self,
        record: dict[str, object],
    ) -> tuple[bool, str]:
        """Appends to errors.tsv and atomically updates last_error.json and marker.

        Returns (success: bool, error_id: str)
        """
        error_id = str(record["error_id"])
        error_key = str(record.get("error_key", ""))
        writer_pid = str(record.get("pid", os.getpid()))
        created_at = str(record.get("timestamp_utc", ""))

        self.logs_dir.mkdir(parents=True, exist_ok=True)

        # 1. Lock and append to errors.tsv
        tsv_row = [sanitize_tsv_cell(record.get(col, "unknown")) for col in TSV_HEADER]
        with open(self.tsv_path, "a+", encoding="utf-8") as f:
            if HAS_FCNTL:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                if size == 0:
                    # Write exact header
                    f.write("\t".join(TSV_HEADER) + "\n")

                # Check if this exact error_key was already persisted for this RUN_ID
                f.seek(0)
                reader = csv.reader(f, delimiter="\t")
                rows = list(reader)
                already_persisted = False
                for r in rows:
                    if r and r[0] == self.run_id and len(r) > 1 and r[1] == error_id:
                        already_persisted = True
                        break

                if not already_persisted:
                    f.seek(0, os.SEEK_END)
                    f.write("\t".join(tsv_row) + "\n")
                    f.flush()
                    try:
                        os.fsync(f.fileno())
                    except OSError:
                        pass
            finally:
                if HAS_FCNTL:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

        # 2. Materialize last_error.json atomically
        json_obj = {col: sanitize_tsv_cell(record.get(col, "unknown")) for col in TSV_HEADER}
        json_obj.update(
            schema_version=SCHEMA_VERSION,
            pipeline_version=str(record.get("pipeline_version", "2.0.0")),
            hostname=str(record.get("hostname", "localhost")),
            pid=int(record.get("pid", os.getpid())),
            error_key=error_key,
            checkpoint_recovery_hint=str(record.get("checkpoint_recovery_hint", "")),
            log_tail=record.get("log_tail", []),
            temporary_directory=str(record.get("temporary_directory", "none")),
            temporary_files_preserved=bool(record.get("temporary_files_preserved", False)),
            temporary_cleanup_status=str(record.get("temporary_cleanup_status", "not-applicable")),
        )

        try:
            with tempfile.NamedTemporaryFile("w", dir=self.logs_dir, delete=False, prefix=".last_error_tmp_", encoding="utf-8") as jf:
                tmp_json = jf.name
                json.dump(json_obj, jf, indent=2, ensure_ascii=False)
                jf.write("\n")
                jf.flush()
                try:
                    os.fsync(jf.fileno())
                except OSError:
                    pass
            os.replace(tmp_json, self.json_path)
        except Exception as e:
            # Secondary view failure must not revert primary TSV
            print(f"[WARN] Failed to materialize last_error.json: {e}", file=sys.stderr)

        # 3. Mark committed
        committed_at = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        marker_data = {
            "state": "committed",
            "writer_pid": writer_pid,
            "created_at": created_at or committed_at,
            "error_key": error_key,
            "error_id": error_id,
            "committed_at": committed_at,
        }
        write_marker_atomic(self.marker_path, marker_data)
        return True, error_id


def build_parser() -> argparse.ArgumentParser:
    """Builds argument parser for write_error_record.py."""
    parser = argparse.ArgumentParser(description="Structured pipeline error recorder")
    subparsers = parser.add_subparsers(dest="action", required=True)

    # Subcommand: reserve
    p_res = subparsers.add_parser("reserve", help="Atomically reserve error marker")
    p_res.add_argument("--outdir", required=True, help="Pipeline outdir")
    p_res.add_argument("--run-id", required=True, help="Current RUN_ID")
    p_res.add_argument("--stage-id", default="pipeline", help="Stage internal ID")
    p_res.add_argument("--stage-name", default="unknown")
    p_res.add_argument("--stage-action", default="unknown")
    p_res.add_argument("--exit-code", type=int, default=1)
    p_res.add_argument("--trap-source-file", default="unknown")
    p_res.add_argument("--trap-source-line", default="0")
    p_res.add_argument("--call-source-file", default="unknown")
    p_res.add_argument("--call-source-line", default="0")
    p_res.add_argument("--command", default="unknown")
    p_res.add_argument("--pid", type=int, default=os.getpid())

    # Subcommand: commit
    p_com = subparsers.add_parser("commit", help="Commit error record to TSV, JSON, and marker")
    p_com.add_argument("--outdir", required=True)
    p_com.add_argument("--run-id", required=True)
    p_com.add_argument("--error-id", required=True)
    p_com.add_argument("--error-key", default="")
    p_com.add_argument("--stage-id", default="pipeline")
    p_com.add_argument("--stage-name", default="unknown")
    p_com.add_argument("--stage-action", default="unknown")
    p_com.add_argument("--exit-code", type=int, default=1)
    p_com.add_argument("--signal", default="none")
    p_com.add_argument("--trap-source-file", default="unknown")
    p_com.add_argument("--trap-source-line", default="0")
    p_com.add_argument("--call-source-file", default="unknown")
    p_com.add_argument("--call-source-line", default="0")
    p_com.add_argument("--function-name", default="unknown")
    p_com.add_argument("--command", default="unknown")
    p_com.add_argument("--log-path", default="unknown")
    p_com.add_argument("--runtime", default="unknown")
    p_com.add_argument("--last-valid-checkpoint", default="none")
    p_com.add_argument("--last-valid-checkpoint-path", default="none")
    p_com.add_argument("--last-valid-product", default="none")
    p_com.add_argument("--last-valid-product-fingerprint", default="none")
    p_com.add_argument("--message", default="unknown")
    p_com.add_argument("--status", default="failed")
    p_com.add_argument("--run-status-update-status", default="not-available")
    p_com.add_argument("--pipeline-version", default="2.0.0")
    p_com.add_argument("--hostname", default="localhost")
    p_com.add_argument("--pid", type=int, default=os.getpid())
    p_com.add_argument("--timestamp-utc", default="")
    p_com.add_argument("--recovery-hint", default="")
    p_com.add_argument("--log-tail-lines", type=int, default=20)
    p_com.add_argument("--temporary-directory", default="none")
    p_com.add_argument("--temporary-files-preserved", default="false")
    p_com.add_argument("--temporary-cleanup-status", default="not-applicable")

    # Subcommand: record (all-in-one)
    p_rec = subparsers.add_parser("record", help="Atomically reserve and commit error in one invocation")
    p_rec.add_argument("--outdir", required=True)
    p_rec.add_argument("--run-id", required=True)
    p_rec.add_argument("--stage-id", default="pipeline")
    p_rec.add_argument("--stage-name", default="unknown")
    p_rec.add_argument("--stage-action", default="unknown")
    p_rec.add_argument("--exit-code", type=int, default=1)
    p_rec.add_argument("--signal", default="none")
    p_rec.add_argument("--trap-source-file", default="unknown")
    p_rec.add_argument("--trap-source-line", default="0")
    p_rec.add_argument("--call-source-file", default="unknown")
    p_rec.add_argument("--call-source-line", default="0")
    p_rec.add_argument("--function-name", default="unknown")
    p_rec.add_argument("--command", default="unknown")
    p_rec.add_argument("--log-path", default="unknown")
    p_rec.add_argument("--runtime", default="unknown")
    p_rec.add_argument("--last-valid-checkpoint", default="none")
    p_rec.add_argument("--last-valid-checkpoint-path", default="none")
    p_rec.add_argument("--last-valid-product", default="none")
    p_rec.add_argument("--last-valid-product-fingerprint", default="none")
    p_rec.add_argument("--message", default="unknown")
    p_rec.add_argument("--status", default="failed")
    p_rec.add_argument("--run-status-update-status", default="not-available")
    p_rec.add_argument("--pipeline-version", default="2.0.0")
    p_rec.add_argument("--hostname", default="localhost")
    p_rec.add_argument("--pid", type=int, default=os.getpid())
    p_rec.add_argument("--timestamp-utc", default="")
    p_rec.add_argument("--recovery-hint", default="")
    p_rec.add_argument("--log-tail-lines", type=int, default=20)
    p_rec.add_argument("--temporary-directory", default="none")
    p_rec.add_argument("--temporary-files-preserved", default="false")
    p_rec.add_argument("--temporary-cleanup-status", default="not-applicable")

    # Subcommand: check
    p_chk = subparsers.add_parser("check", help="Check current marker status")
    p_chk.add_argument("--outdir", required=True)
    p_chk.add_argument("--run-id", required=True)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    mgr = ErrorRecordManager(args.outdir, args.run_id)

    if args.action == "check":
        marker = parse_marker(mgr.marker_path)
        if not marker:
            print("STATE=none")
            return 0
        for k, v in marker.items():
            print(f"{k.upper()}={v}")
        return 0

    if args.action == "reserve":
        err_key = compute_error_key(
            run_id=args.run_id,
            stage_id=args.stage_id,
            exit_code=args.exit_code,
            trap_file=args.trap_source_file,
            trap_line=args.trap_source_line,
            call_file=args.call_source_file,
            call_line=args.call_source_line,
            command=args.command,
        )
        is_primary, error_id, state, key = mgr.reserve(
            writer_pid=args.pid,
            error_key=err_key,
            stage_id=args.stage_id,
        )
        print(f"IS_PRIMARY={'true' if is_primary else 'false'}")
        print(f"ERROR_ID={error_id}")
        print(f"ERROR_KEY={key}")
        print(f"STATE={state}")
        return 0

    if args.action == "commit":
        ts = args.timestamp_utc or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        err_key = args.error_key or compute_error_key(
            run_id=args.run_id,
            stage_id=args.stage_id,
            exit_code=args.exit_code,
            trap_file=args.trap_source_file,
            trap_line=args.trap_source_line,
            call_file=args.call_source_file,
            call_line=args.call_source_line,
            command=args.command,
        )
        rec_hint = get_recovery_hint(args.stage_id, args.recovery_hint)
        log_tail = extract_log_tail(args.log_path, max_lines=args.log_tail_lines)

        rec = {
            "run_id": args.run_id,
            "error_id": args.error_id,
            "timestamp_utc": ts,
            "stage_id": args.stage_id,
            "stage_name": args.stage_name,
            "stage_action": args.stage_action,
            "exit_code": args.exit_code,
            "signal": args.signal if args.signal in ALLOWED_SIGNALS else "none",
            "trap_source_file": args.trap_source_file,
            "trap_source_line": args.trap_source_line,
            "call_source_file": args.call_source_file,
            "call_source_line": args.call_source_line,
            "function_name": args.function_name,
            "command": args.command,
            "log_path": args.log_path,
            "runtime": args.runtime,
            "last_valid_checkpoint": args.last_valid_checkpoint,
            "last_valid_checkpoint_path": args.last_valid_checkpoint_path,
            "last_valid_product": args.last_valid_product,
            "last_valid_product_fingerprint": args.last_valid_product_fingerprint,
            "message": args.message,
            "status": args.status if args.status in ALLOWED_STATUSES else "failed",
            "run_status_update_status": args.run_status_update_status if args.run_status_update_status in ALLOWED_RUN_STATUS_UPDATE else "not-available",
            "pipeline_version": args.pipeline_version,
            "hostname": args.hostname,
            "pid": args.pid,
            "error_key": err_key,
            "checkpoint_recovery_hint": rec_hint,
            "log_tail": log_tail,
            "temporary_directory": getattr(args, "temporary_directory", "none"),
            "temporary_files_preserved": getattr(args, "temporary_files_preserved", "false").lower() in ("true", "1", "yes"),
            "temporary_cleanup_status": getattr(args, "temporary_cleanup_status", "not-applicable"),
        }

        mgr.commit(rec)
        print("IS_PRIMARY=true")
        print(f"ERROR_ID={args.error_id}")
        print(f"ERROR_KEY={err_key}")
        print("STATE=committed")
        return 0

    if args.action == "record":
        err_key = compute_error_key(
            run_id=args.run_id,
            stage_id=args.stage_id,
            exit_code=args.exit_code,
            trap_file=args.trap_source_file,
            trap_line=args.trap_source_line,
            call_file=args.call_source_file,
            call_line=args.call_source_line,
            command=args.command,
        )
        is_primary, error_id, state, key = mgr.reserve(
            writer_pid=args.pid,
            error_key=err_key,
            stage_id=args.stage_id,
        )
        if not is_primary:
            print("IS_PRIMARY=false")
            print(f"ERROR_ID={error_id}")
            print(f"ERROR_KEY={key}")
            print(f"STATE={state}")
            return 0

        ts = args.timestamp_utc or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        rec_hint = get_recovery_hint(args.stage_id, args.recovery_hint)
        log_tail = extract_log_tail(args.log_path, max_lines=args.log_tail_lines)

        rec = {
            "run_id": args.run_id,
            "error_id": error_id,
            "timestamp_utc": ts,
            "stage_id": args.stage_id,
            "stage_name": args.stage_name,
            "stage_action": args.stage_action,
            "exit_code": args.exit_code,
            "signal": args.signal if args.signal in ALLOWED_SIGNALS else "none",
            "trap_source_file": args.trap_source_file,
            "trap_source_line": args.trap_source_line,
            "call_source_file": args.call_source_file,
            "call_source_line": args.call_source_line,
            "function_name": args.function_name,
            "command": args.command,
            "log_path": args.log_path,
            "runtime": args.runtime,
            "last_valid_checkpoint": args.last_valid_checkpoint,
            "last_valid_checkpoint_path": args.last_valid_checkpoint_path,
            "last_valid_product": args.last_valid_product,
            "last_valid_product_fingerprint": args.last_valid_product_fingerprint,
            "message": args.message,
            "status": args.status if args.status in ALLOWED_STATUSES else "failed",
            "run_status_update_status": args.run_status_update_status if args.run_status_update_status in ALLOWED_RUN_STATUS_UPDATE else "not-available",
            "pipeline_version": args.pipeline_version,
            "hostname": args.hostname,
            "pid": args.pid,
            "error_key": err_key,
            "checkpoint_recovery_hint": rec_hint,
            "log_tail": log_tail,
            "temporary_directory": getattr(args, "temporary_directory", "none"),
            "temporary_files_preserved": getattr(args, "temporary_files_preserved", "false").lower() in ("true", "1", "yes"),
            "temporary_cleanup_status": getattr(args, "temporary_cleanup_status", "not-applicable"),
        }
        mgr.commit(rec)
        print("IS_PRIMARY=true")
        print(f"ERROR_ID={error_id}")
        print(f"ERROR_KEY={err_key}")
        print("STATE=committed")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
