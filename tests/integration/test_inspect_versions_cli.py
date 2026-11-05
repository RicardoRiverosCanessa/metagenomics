# -*- coding: utf-8 -*-
"""Integration tests for scripts.inspect_versions (CLI layer).
These tests invoke the module as a script, exercise argument parsing,
conflict handling, and ensure no side‑effects such as output directory
creation or external tool invocation.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest
try:
    import pytest
    pytestmark = pytest.mark.integration
except ImportError:
    pass

from scripts.inspect_versions import (
    EXIT_OK,
    EXIT_STRICT_FAILED,
    EXIT_USAGE_ERROR,
)


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

    def test_cli_bash_version_flag(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res = subprocess.run([script_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("0.3.2", res.stdout)
        self.assertIn("development", res.stdout)

    def test_cli_bash_short_v_flag(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res = subprocess.run([script_path, "-V"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("0.3.2", res.stdout)

    def test_cli_bash_version_and_versions_conflict_exits_2(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res1 = subprocess.run([script_path, "--version", "--versions"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res1.returncode, 2)
        self.assertIn("Cannot combine", res1.stderr)
        res2 = subprocess.run([script_path, "-V", "--versions"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res2.returncode, 2)
        self.assertIn("Cannot combine", res2.stderr)

    def test_cli_bash_versions_help_exits_0(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res = subprocess.run([script_path, "--versions", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("Version & Runtime Inspection", res.stdout)
        self.assertIn("--versions", res.stdout)

    def test_cli_bash_general_help_exits_0(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res = subprocess.run([script_path, "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("MetaNexus", res.stdout)
        self.assertIn("--versions", res.stdout)

    def test_cli_bash_versions_format_invalid_exits_2(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        res = subprocess.run([script_path, "--versions", "--versions-format", "unsupported_fmt"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 2)

    def test_cli_bash_versions_json_does_not_create_outdir(self):
        script_path = str(Path(__file__).resolve().parents[2] / "run_metagenomics_pipeline.sh")
        test_outdir = os.path.join(tempfile.gettempdir(), f"test_outdir_check_{os.getpid()}")
        try:
            cmd = [
                script_path,
                "--versions",
                "--profile",
                "taxonomy",
                "--versions-format",
                "json",
                "-o",
                test_outdir,
            ]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.assertEqual(res.returncode, 0)
            data = json.loads(res.stdout)
            self.assertEqual(data["pipeline"]["name"], "MetaNexus")
            self.assertFalse(os.path.exists(test_outdir))
        finally:
            shutil.rmtree(test_outdir, ignore_errors=True)
