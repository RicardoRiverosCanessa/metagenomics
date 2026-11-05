#!/usr/bin/env python3
"""Unit tests for scripts/inspect_versions.py."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from scripts.inspect_versions import (
    DEFAULT_PIPELINE_NAME,
    DEFAULT_PIPELINE_RELEASE_STATUS,
    DEFAULT_PIPELINE_VERSION,
    EXIT_OK,
    EXIT_STRICT_FAILED,
    EXIT_USAGE_ERROR,
    DatabaseResult,
    InspectionCache,
    InspectionReport,
    PipelineIdentity,
    RuntimeResult,
    ToolResult,
    VersionInspector,
    compute_file_sha256,
    extract_reported_version,
    format_report_json,
    format_report_text,
    main,
    run_command_with_limits,
)


class TestPipelineIdentity(unittest.TestCase):
    def test_default_identity(self):
        ident = PipelineIdentity()
        self.assertEqual(ident.name, "MetaNexus")
        self.assertEqual(ident.version, "0.3.2")
        self.assertEqual(ident.release_status, "development")

    def test_identity_dict(self):
        ident = PipelineIdentity(name="MetaNexus", version="0.3.2", release_status="development")
        d = ident.to_dict()
        self.assertEqual(d["name"], "MetaNexus")
        self.assertEqual(d["version"], "0.3.2")
        self.assertEqual(d["release_status"], "development")


class TestVersionExtraction(unittest.TestCase):
    def test_extract_standard_semantic_versions(self):
        cases = [
            ("fastp 0.23.4", "0.23.4"),
            ("fastp version 0.23.4\n", "0.23.4"),
            ("MEGAHIT v1.2.9", "1.2.9"),
            ("Kraken version 2.1.3", "2.1.3"),
            ("Bracken version 2.9", "2.9"),
            ("bowtie2-align-s version 2.5.4\n64-bit", "2.5.4"),
            ("samtools 1.19.2\nUsing htslib 1.19.2", "1.19.2"),
            ("Vamb 4.1.3", "4.1.3"),
            ("MetaBinner v1.4.4", "1.4.4"),
            ("COMEBin 1.0.4", "1.0.4"),
            ("metawrap v1.3.2", "1.3.2"),
            ("CheckM2 version 1.0.2", "1.0.2"),
            ("gtdbtk 2.4.0", "2.4.0"),
            ("CoverM 0.6.1", "0.6.1"),
            ("R scripting front-end version 4.3.1 (2023-06-16)", "4.3.1"),
        ]
        for text, expected in cases:
            with self.subTest(text=text):
                self.assertEqual(extract_reported_version(text), expected)

    def test_extract_commit_hash(self):
        text = "commit a1b2c3d4e5f678901234567890abcdef12345678"
        self.assertEqual(extract_reported_version(text), "a1b2c3d4e5f678901234567890abcdef12345678")

    def test_extract_short_commit_hash(self):
        text = "revision 7a8b9c0"
        self.assertEqual(extract_reported_version(text), "7a8b9c0")

    def test_extract_unknown_and_empty(self):
        self.assertEqual(extract_reported_version(""), "unknown")
        self.assertEqual(extract_reported_version("No numeric version here"), "unknown")


class TestCommandExecutionWithLimits(unittest.TestCase):
    def test_run_command_success(self):
        rc, out, timed_out = run_command_with_limits([sys.executable, "-c", "print('hello 1.0.0')"], timeout=5.0)
        self.assertEqual(rc, 0)
        self.assertFalse(timed_out)
        self.assertIn("hello 1.0.0", out)

    def test_run_command_timeout(self):
        rc, out, timed_out = run_command_with_limits(
            [sys.executable, "-c", "import time; time.sleep(3)"],
            timeout=0.2,
        )
        self.assertTrue(timed_out)

    def test_run_command_truncation(self):
        # Generate 100 KB of text and verify it is truncated at 1000 bytes
        rc, out, timed_out = run_command_with_limits(
            [sys.executable, "-c", "print('A' * 100000)"],
            timeout=5.0,
            max_bytes=1000,
        )
        self.assertEqual(rc, 0)
        self.assertFalse(timed_out)
        self.assertLessEqual(len(out.encode("utf-8")), 1000)

    def test_run_command_not_found(self):
        rc, out, timed_out = run_command_with_limits(["__non_existent_command_12345__"], timeout=2.0)
        self.assertEqual(rc, 127)
        self.assertFalse(timed_out)


class TestInspectionCache(unittest.TestCase):
    def test_cache_reuse(self):
        cache = InspectionCache()
        cache.file_hashes["/fake/path.sif"] = "dummy_sha256"

        ident = PipelineIdentity()
        inspector = VersionInspector(pipeline=ident, cache=cache)
        self.assertEqual(inspector.cache.file_hashes.get("/fake/path.sif"), "dummy_sha256")


class TestDatabaseInspection(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_db_inspect_")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_not_configured_database(self):
        inspector = VersionInspector()
        spec = {"stage": "03_kraken2", "name": "kraken2-bracken", "cli_flag": "--kraken2-db"}
        res = inspector.inspect_database(spec, user_path=None)
        self.assertEqual(res.status, "not-configured")
        self.assertEqual(res.release, "unknown")

    def test_checkm2_diamond_file_exists(self):
        dmnd_path = os.path.join(self.test_dir, "uniref100.KO.1.dmnd")
        with open(dmnd_path, "wb") as f:
            f.write(b"mock diamond database content")

        inspector = VersionInspector(deep=True)
        spec = {"stage": "08_checkm2", "name": "checkm2", "is_file": True}
        res = inspector.inspect_database(spec, user_path=dmnd_path)
        self.assertEqual(res.status, "available")
        self.assertEqual(res.release, "uniref100.KO.1.dmnd")
        self.assertEqual(res.fingerprint_method, "sha256-file")
        self.assertTrue(len(res.fingerprint) == 64)

    def test_checkm2_diamond_file_missing(self):
        missing_path = os.path.join(self.test_dir, "missing.dmnd")
        inspector = VersionInspector()
        spec = {"stage": "08_checkm2", "name": "checkm2", "is_file": True}
        res = inspector.inspect_database(spec, user_path=missing_path)
        self.assertEqual(res.status, "unavailable")

    def test_kraken2_database_with_missing_keyfiles(self):
        k2_dir = os.path.join(self.test_dir, "k2_db")
        os.makedirs(k2_dir, exist_ok=True)
        # Create only hash.k2d, omitting opts.k2d and taxo.k2d
        Path(k2_dir, "hash.k2d").write_text("fake")

        inspector = VersionInspector()
        spec = {
            "stage": "03_kraken2",
            "name": "kraken2-bracken",
            "required_files": ["hash.k2d", "opts.k2d", "taxo.k2d"],
        }
        res = inspector.inspect_database(spec, user_path=k2_dir)
        self.assertEqual(res.status, "unavailable")
        self.assertIn("opts.k2d", res.details)

    def test_kraken2_database_valid_with_manifest(self):
        k2_dir = os.path.join(self.test_dir, "k2_db")
        os.makedirs(k2_dir, exist_ok=True)
        for req in ["hash.k2d", "opts.k2d", "taxo.k2d"]:
            Path(k2_dir, req).write_text("fake-content")
        # Write installation.json
        manifest = {"catalog_id": "kraken2-standard-8", "release": "20260626", "status": "installed"}
        Path(k2_dir, "installation.json").write_text(json.dumps(manifest))

        inspector = VersionInspector(deep=True)
        spec = {
            "stage": "03_kraken2",
            "name": "kraken2-bracken",
            "required_files": ["hash.k2d", "opts.k2d", "taxo.k2d"],
        }
        res = inspector.inspect_database(spec, user_path=k2_dir)
        self.assertEqual(res.status, "available")
        self.assertEqual(res.release, "20260626")
        self.assertEqual(res.version_source, "installation-manifest")
        self.assertEqual(res.fingerprint_method, "metadata-keyfiles-v1")

    def test_gtdbtk_database_inspection(self):
        gtdb_dir = os.path.join(self.test_dir, "gtdb")
        os.makedirs(gtdb_dir, exist_ok=True)
        Path(gtdb_dir, "VERSION").write_text("r220\n")

        inspector = VersionInspector(deep=True)
        spec = {
            "stage": "10_gtdbtk",
            "name": "gtdbtk",
            "version_files": ["VERSION", "release.txt"],
        }
        res = inspector.inspect_database(spec, user_path=gtdb_dir)
        self.assertEqual(res.status, "available")
        self.assertEqual(res.release, "r220")
        # Prompt requirement: deep fingerprint for GTDB-Tk must be omitted
        self.assertEqual(res.fingerprint_method, "not-applicable")
        self.assertEqual(res.fingerprint, "omitted-gtdbtk")


class TestRuntimeAndToolInspection(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_runtime_")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_apptainer_runtime_missing_sif(self):
        inspector = VersionInspector(engine="apptainer", sif_dir=self.test_dir)
        r_res = inspector.inspect_runtime(
            stage="01_fastp",
            tool="fastp",
            sif_name="fastp.sif",
            uri="docker://fastp",
            env_name="",
            executable="fastp",
        )
        self.assertEqual(r_res.status, "unavailable")
        self.assertIn("not found", r_res.error_message or "")

        t_res = inspector.inspect_tool(
            {"stage": "01_fastp", "tool": "fastp", "executable": "fastp", "query": "--version"},
            r_res,
        )
        self.assertEqual(t_res.status, "unavailable")

    def test_mamba_runtime_missing_env(self):
        inspector = VersionInspector(engine="mamba")
        r_res = inspector.inspect_runtime(
            stage="05a_vamb",
            tool="vamb",
            sif_name="vamb.sif",
            uri="",
            env_name="__non_existent_env_12345__",
            executable="vamb",
        )
        self.assertEqual(r_res.status, "unavailable")

    def test_host_runtime_existing_tool(self):
        inspector = VersionInspector(engine="mamba")
        r_res = inspector.inspect_runtime(
            stage="01_test",
            tool="python3",
            sif_name="",
            uri="",
            env_name="",
            executable=sys.executable,
        )
        self.assertEqual(r_res.status, "available")
        t_res = inspector.inspect_tool(
            {"stage": "01_test", "tool": "python3", "executable": sys.executable, "query": "--version"},
            r_res,
        )
        self.assertEqual(t_res.status, "available")
        self.assertNotEqual(t_res.version, "unknown")


class TestReportFormatting(unittest.TestCase):
    def test_report_json_format(self):
        ident = PipelineIdentity(name="MetaNexus", version="0.3.2", release_status="development")
        report = InspectionReport(
            pipeline=ident,
            inspected_at_utc="2026-11-04T18:00:00Z",
            execution_mode="apptainer",
            profile_scope="full",
            deep=False,
            all_stages=False,
            strict=False,
            tools=[
                ToolResult(stage="01_fastp", tool="fastp", role="primary", version="0.23.4", status="available")
            ],
            runtimes=[
                RuntimeResult(stage="01_fastp", tool="fastp", mode="apptainer", name="fastp", uri="docker://fastp", status="available")
            ],
            databases=[
                DatabaseResult(stage="03_kraken2", name="kraken2-bracken", path="<unset>", status="not-configured")
            ],
            summary={
                "total_tools": 1,
                "available_tools": 1,
                "unavailable_tools": 0,
                "unknown_tools": 0,
                "total_databases": 1,
                "available_databases": 0,
                "unavailable_databases": 0,
                "not_configured_databases": 1,
                "status": "ok",
            },
        )
        json_str = format_report_json(report)
        data = json.loads(json_str)
        self.assertEqual(data["pipeline"]["name"], "MetaNexus")
        self.assertEqual(data["pipeline"]["version"], "0.3.2")
        self.assertEqual(data["summary"]["status"], "ok")
        self.assertEqual(len(data["tools"]), 1)
        self.assertEqual(data["tools"][0]["version"], "0.23.4")

    def test_report_text_format(self):
        ident = PipelineIdentity(name="MetaNexus", version="0.3.2", release_status="development")
        report = InspectionReport(
            pipeline=ident,
            inspected_at_utc="2026-11-04T18:00:00Z",
            execution_mode="apptainer",
            profile_scope="full",
            deep=False,
            all_stages=False,
            strict=False,
            tools=[
                ToolResult(stage="01_fastp", tool="fastp", role="primary", version="0.23.4", status="available")
            ],
            runtimes=[
                RuntimeResult(stage="01_fastp", tool="fastp", mode="apptainer", name="fastp", uri="docker://fastp", status="available")
            ],
            databases=[],
            summary={"status": "ok", "total_tools": 1, "available_tools": 1},
        )
        text_str = format_report_text(report)
        self.assertIn("MetaNexus 0.3.2 (development)", text_str)
        self.assertIn("[TOOLS]", text_str)
        self.assertIn("01_fastp", text_str)
        self.assertIn("0.23.4", text_str)
        self.assertIn("[SUMMARY]", text_str)


class TestCliExecution(unittest.TestCase):
    def test_cli_format_json(self):
        cmd = [
            sys.executable,
            "-m",
            "scripts.inspect_versions",
            "--stages",
            "none",
            "--format",
            "json",
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        data = json.loads(res.stdout)
        self.assertEqual(data["pipeline"]["name"], "MetaNexus")
        self.assertEqual(data["pipeline"]["version"], "0.3.2")

    def test_cli_strict_mode_failure(self):
        # Inspect a stage with a non-existent SIF in a dummy directory under strict mode
        temp_sif_dir = tempfile.mkdtemp(prefix="test_sif_strict_")
        try:
            cmd = [
                sys.executable,
                "-m",
                "scripts.inspect_versions",
                "--stages",
                "01_fastp",
                "--sif-dir",
                temp_sif_dir,
                "--strict",
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, EXIT_STRICT_FAILED)
        finally:
            shutil.rmtree(temp_sif_dir, ignore_errors=True)

    def test_cli_non_strict_mode_success_with_warnings(self):
        # In non-strict mode, missing components do not exit with error
        temp_sif_dir = tempfile.mkdtemp(prefix="test_sif_non_strict_")
        try:
            cmd = [
                sys.executable,
                "-m",
                "scripts.inspect_versions",
                "--stages",
                "01_fastp",
                "--sif-dir",
                temp_sif_dir,
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, EXIT_OK)
        finally:
            shutil.rmtree(temp_sif_dir, ignore_errors=True)

    def test_cli_invalid_format_exits_2(self):
        cmd = [
            sys.executable,
            "-m",
            "scripts.inspect_versions",
            "--format",
            "invalid_format_name",
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, EXIT_USAGE_ERROR)


class TestBashPipelineVersionsCli(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script_path = str(Path(__file__).resolve().parent.parent / "run_metagenomics_pipeline.sh")

    def test_bash_version_flag(self):
        res = subprocess.run([self.script_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("0.3.2", res.stdout)
        self.assertIn("development", res.stdout)

    def test_bash_short_v_flag(self):
        res = subprocess.run([self.script_path, "-V"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("0.3.2", res.stdout)

    def test_bash_profile_and_version(self):
        res = subprocess.run([self.script_path, "--profile", "full", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)

    def test_bash_version_and_versions_conflict_exits_2(self):
        res1 = subprocess.run([self.script_path, "--version", "--versions"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res1.returncode, 2)
        self.assertIn("Cannot combine", res1.stderr)

        res2 = subprocess.run([self.script_path, "-V", "--versions"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res2.returncode, 2)
        self.assertIn("Cannot combine", res2.stderr)

    def test_bash_versions_help_exits_0(self):
        res = subprocess.run([self.script_path, "--versions", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("Version & Runtime Inspection", res.stdout)
        self.assertIn("--versions", res.stdout)

    def test_bash_general_help_exits_0(self):
        res = subprocess.run([self.script_path, "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("--versions", res.stdout)

    def test_bash_versions_format_invalid_exits_2(self):
        res = subprocess.run([self.script_path, "--versions", "--versions-format", "unsupported_fmt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 2)

    def test_bash_versions_json_does_not_create_outdir(self):
        test_outdir = os.path.join(tempfile.gettempdir(), f"test_outdir_check_{os.getpid()}")
        try:
            res = subprocess.run(
                [self.script_path, "--versions", "--profile", "taxonomy", "--versions-format", "json", "-o", test_outdir],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(res.returncode, 0)
            data = json.loads(res.stdout)
            self.assertEqual(data["pipeline"]["name"], "MetaNexus")
            self.assertFalse(os.path.exists(test_outdir))
        finally:
            shutil.rmtree(test_outdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
