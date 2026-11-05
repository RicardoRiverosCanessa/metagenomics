#!/usr/bin/env python3
"""Unit tests for scripts/write_disk_space_report.py"""

import os
import shutil
import tempfile
import unittest

from scripts.write_disk_space_report import (
    ALLOWED_CONFIDENCE,
    ALLOWED_MEASUREMENT_TYPES,
    ALLOWED_STATUS,
    ALLOWED_STORAGE_ROLES,
    SCHEMA_COLUMNS,
    cmd_publish,
    cmd_record,
    cmd_summary,
    load_report_tsv,
    sanitize_message,
    validate_row_dict,
)


class TestWriteDiskSpaceReport(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_disk_space_")
        self.current_tsv = os.path.join(self.test_dir, "disk_space.current.tsv")
        self.public_tsv = os.path.join(self.test_dir, "disk_space.tsv")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_schema_columns_count_and_order(self):
        self.assertEqual(len(SCHEMA_COLUMNS), 22)
        expected = [
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
        self.assertEqual(SCHEMA_COLUMNS, expected)

    def test_sanitize_message(self):
        raw = "msg\twith\ttabs\r\nand\nnewlines   and   extra spaces"
        sanitized = sanitize_message(raw)
        self.assertNotIn("\t", sanitized)
        self.assertNotIn("\r", sanitized)
        self.assertNotIn("\n", sanitized)
        self.assertEqual(sanitized, "msg with tabs; and; newlines and extra spaces")

    def test_validate_row_dict_enums(self):
        base_row = {
            "run_id": "RUN_001",
            "measurement_sequence": "1",
            "check_id": "pipeline.outdir.outputs",
            "timestamp_utc": "2026-11-03T18:00:00Z",
            "stage_id": "pipeline",
            "storage_role": "outdir",
            "purpose": "pipeline-outputs",
            "path": "/data/out",
            "filesystem_id": "2049|/data",
            "measurement_type": "preflight-estimate",
            "estimate_confidence": "hard-minimum",
            "input_bytes_read": "1000",
            "estimated_required_bytes": "2000",
            "filesystem_available_bytes": "50000",
            "filesystem_used_bytes": "10000",
            "filesystem_total_bytes": "100000",
            "required_free_inodes": "16",
            "available_free_inodes": "5000",
            "safety_factor": "1.25",
            "shared_capacity_group": "none",
            "status": "sufficient",
            "message": "filesystem_source=real",
        }
        val = validate_row_dict(base_row, expected_run_id="RUN_001")
        self.assertEqual(val["run_id"], "RUN_001")
        self.assertEqual(val["storage_role"], "outdir")
        self.assertEqual(val["status"], "sufficient")

        # Invalid storage_role
        bad_row = dict(base_row)
        bad_row["storage_role"] = "input"
        with self.assertRaises(ValueError):
            validate_row_dict(bad_row)

        # Invalid measurement_type
        bad_row = dict(base_row)
        bad_row["measurement_type"] = "mid-stage"
        with self.assertRaises(ValueError):
            validate_row_dict(bad_row)

        # Invalid estimate_confidence
        bad_row = dict(base_row)
        bad_row["estimate_confidence"] = "guess"
        with self.assertRaises(ValueError):
            validate_row_dict(bad_row)

        # Invalid status
        bad_row = dict(base_row)
        bad_row["status"] = "ok"
        with self.assertRaises(ValueError):
            validate_row_dict(bad_row)

    def test_record_single_and_cumulative(self):
        class Args:
            pass

        args1 = Args()
        args1.report_current = self.current_tsv
        args1.publish_path = self.public_tsv
        args1.run_id = "RUN_001"
        args1.measurement_sequence = "auto"
        args1.check_id = "01_fastp.outdir.cleaned-reads"
        args1.timestamp_utc = "2026-11-03T18:00:00Z"
        args1.stage_id = "01_fastp"
        args1.storage_role = "outdir"
        args1.purpose = "cleaned-reads"
        args1.path = "/data/out/01_fastp"
        args1.filesystem_id = "2049|/data"
        args1.measurement_type = "stage-precheck"
        args1.estimate_confidence = "informational"
        args1.input_bytes_read = "100000"
        args1.estimated_required_bytes = "110000"
        args1.filesystem_available_bytes = "500000"
        args1.filesystem_used_bytes = "200000"
        args1.filesystem_total_bytes = "1000000"
        args1.required_free_inodes = "NA"
        args1.available_free_inodes = "10000"
        args1.safety_factor = "1.10"
        args1.shared_capacity_group = "none"
        args1.status = "sufficient"
        args1.message = "input_reads=100000"

        rc1 = cmd_record(args1)
        self.assertEqual(rc1, 0)
        self.assertTrue(os.path.exists(self.current_tsv))
        self.assertTrue(os.path.exists(self.public_tsv))

        rows1 = load_report_tsv(self.current_tsv, expected_run_id="RUN_001")
        self.assertEqual(len(rows1), 1)
        self.assertEqual(rows1[0]["measurement_sequence"], "1")

        # Second record
        args2 = Args()
        args2.report_current = self.current_tsv
        args2.publish_path = self.public_tsv
        args2.run_id = "RUN_001"
        args2.measurement_sequence = "auto"
        args2.check_id = "07_reassembly.tmpdir.concatenated-reads"
        args2.timestamp_utc = "2026-11-03T18:05:00Z"
        args2.stage_id = "07_reassembly"
        args2.storage_role = "tmpdir"
        args2.purpose = "read-concatenation"
        args2.path = "/scratch/tmp"
        args2.filesystem_id = "2050|/scratch"
        args2.measurement_type = "stage-precheck"
        args2.estimate_confidence = "hard-minimum"
        args2.input_bytes_read = "200000"
        args2.estimated_required_bytes = "220000"
        args2.filesystem_available_bytes = "250000"
        args2.filesystem_used_bytes = "500000"
        args2.filesystem_total_bytes = "1000000"
        args2.required_free_inodes = "16"
        args2.available_free_inodes = "8000"
        args2.safety_factor = "1.10"
        args2.shared_capacity_group = "07_reassembly"
        args2.status = "warning"
        args2.message = "headroom_low"

        rc2 = cmd_record(args2)
        self.assertEqual(rc2, 0)

        rows2 = load_report_tsv(self.current_tsv, expected_run_id="RUN_001")
        self.assertEqual(len(rows2), 2)
        self.assertEqual(rows2[0]["measurement_sequence"], "1")
        self.assertEqual(rows2[1]["measurement_sequence"], "2")

    def test_single_run_id_enforcement_and_new_run_reset(self):
        class Args:
            pass

        args1 = Args()
        args1.report_current = self.current_tsv
        args1.publish_path = ""
        args1.run_id = "OLD_RUN"
        args1.measurement_sequence = "1"
        args1.check_id = "01_fastp.outdir.cleaned-reads"
        args1.timestamp_utc = "2026-11-03T18:00:00Z"
        args1.stage_id = "01_fastp"
        args1.storage_role = "outdir"
        args1.purpose = "cleaned-reads"
        args1.path = "/data/out/01_fastp"
        args1.filesystem_id = "2049|/data"
        args1.measurement_type = "stage-precheck"
        args1.estimate_confidence = "informational"
        args1.input_bytes_read = "100000"
        args1.estimated_required_bytes = "110000"
        args1.filesystem_available_bytes = "500000"
        args1.filesystem_used_bytes = "200000"
        args1.filesystem_total_bytes = "1000000"
        args1.required_free_inodes = "NA"
        args1.available_free_inodes = "10000"
        args1.safety_factor = "1.10"
        args1.shared_capacity_group = "none"
        args1.status = "sufficient"
        args1.message = ""

        cmd_record(args1)
        self.assertEqual(len(load_report_tsv(self.current_tsv)), 1)

        # Now record with NEW_RUN
        args2 = Args()
        args2.report_current = self.current_tsv
        args2.publish_path = self.public_tsv
        args2.run_id = "NEW_RUN"
        args2.measurement_sequence = "auto"
        args2.check_id = "pipeline.outdir.outputs"
        args2.timestamp_utc = "2026-11-03T19:00:00Z"
        args2.stage_id = "pipeline"
        args2.storage_role = "outdir"
        args2.purpose = "pipeline-outputs"
        args2.path = "/data/out"
        args2.filesystem_id = "2049|/data"
        args2.measurement_type = "preflight-estimate"
        args2.estimate_confidence = "hard-minimum"
        args2.input_bytes_read = "NA"
        args2.estimated_required_bytes = "1000"
        args2.filesystem_available_bytes = "500000"
        args2.filesystem_used_bytes = "200000"
        args2.filesystem_total_bytes = "1000000"
        args2.required_free_inodes = "NA"
        args2.available_free_inodes = "10000"
        args2.safety_factor = "NA"
        args2.shared_capacity_group = "none"
        args2.status = "sufficient"
        args2.message = "new run fresh start"

        cmd_record(args2)
        new_rows = load_report_tsv(self.current_tsv, expected_run_id="NEW_RUN")
        # Should only contain NEW_RUN rows, not OLD_RUN
        self.assertEqual(len(new_rows), 1)
        self.assertEqual(new_rows[0]["run_id"], "NEW_RUN")
        self.assertEqual(new_rows[0]["measurement_sequence"], "1")

    def test_summary_warning_deduplication(self):
        class Args:
            pass

        # Record preflight warning for check A
        a1 = Args()
        a1.report_current = self.current_tsv
        a1.publish_path = ""
        a1.run_id = "RUN_001"
        a1.measurement_sequence = "1"
        a1.check_id = "check_a"
        a1.timestamp_utc = "2026-11-03T18:00:00Z"
        a1.stage_id = "01_fastp"
        a1.storage_role = "outdir"
        a1.purpose = "cleaned-reads"
        a1.path = "/data"
        a1.filesystem_id = "2049|/data"
        a1.measurement_type = "preflight-estimate"
        a1.estimate_confidence = "informational"
        a1.input_bytes_read = "100"
        a1.estimated_required_bytes = "110"
        a1.filesystem_available_bytes = "120"
        a1.filesystem_used_bytes = "10"
        a1.filesystem_total_bytes = "1000"
        a1.required_free_inodes = "NA"
        a1.available_free_inodes = "NA"
        a1.safety_factor = "1.10"
        a1.shared_capacity_group = "none"
        a1.status = "warning"
        a1.message = "headroom low"
        cmd_record(a1)

        # Record stage precheck warning for check A again (should NOT count as 2 warnings!)
        a2 = Args()
        a2.report_current = self.current_tsv
        a2.publish_path = ""
        a2.run_id = "RUN_001"
        a2.measurement_sequence = "2"
        a2.check_id = "check_a"
        a2.timestamp_utc = "2026-11-03T18:01:00Z"
        a2.stage_id = "01_fastp"
        a2.storage_role = "outdir"
        a2.purpose = "cleaned-reads"
        a2.path = "/data"
        a2.filesystem_id = "2049|/data"
        a2.measurement_type = "stage-precheck"
        a2.estimate_confidence = "informational"
        a2.input_bytes_read = "100"
        a2.estimated_required_bytes = "110"
        a2.filesystem_available_bytes = "120"
        a2.filesystem_used_bytes = "10"
        a2.filesystem_total_bytes = "1000"
        a2.required_free_inodes = "NA"
        a2.available_free_inodes = "NA"
        a2.safety_factor = "1.10"
        a2.shared_capacity_group = "none"
        a2.status = "warning"
        a2.message = "headroom low"
        cmd_record(a2)

        s_args = Args()
        s_args.report_current = self.current_tsv
        s_args.run_id = "RUN_001"
        s_args.json = True
        # Summary should see 1 warning, not 2
        import io
        from contextlib import redirect_stdout
        buf = io.StringIO()
        with redirect_stdout(buf):
            cmd_summary(s_args)
        import json
        res = json.loads(buf.getvalue())
        self.assertEqual(res["total_measurements"], 2)
        self.assertEqual(res["distinct_checks"], 1)
        self.assertEqual(res["warning_count"], 1)
        self.assertEqual(res["overall_status"], "warning")

        # Now check A resolves to sufficient
        a3 = Args()
        a3.report_current = self.current_tsv
        a3.publish_path = ""
        a3.run_id = "RUN_001"
        a3.measurement_sequence = "3"
        a3.check_id = "check_a"
        a3.timestamp_utc = "2026-11-03T18:02:00Z"
        a3.stage_id = "01_fastp"
        a3.storage_role = "outdir"
        a3.purpose = "cleaned-reads"
        a3.path = "/data"
        a3.filesystem_id = "2049|/data"
        a3.measurement_type = "stage-precheck"
        a3.estimate_confidence = "informational"
        a3.input_bytes_read = "100"
        a3.estimated_required_bytes = "110"
        a3.filesystem_available_bytes = "5000"
        a3.filesystem_used_bytes = "10"
        a3.filesystem_total_bytes = "1000"
        a3.required_free_inodes = "NA"
        a3.available_free_inodes = "NA"
        a3.safety_factor = "1.10"
        a3.shared_capacity_group = "none"
        a3.status = "sufficient"
        a3.message = "space recovered"
        cmd_record(a3)

        buf2 = io.StringIO()
        with redirect_stdout(buf2):
            cmd_summary(s_args)
        res2 = json.loads(buf2.getvalue())
        self.assertEqual(res2["total_measurements"], 3)
        self.assertEqual(res2["warning_count"], 0)
        self.assertEqual(res2["overall_status"], "sufficient")


if __name__ == "__main__":
    unittest.main()
