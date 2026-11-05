#!/usr/bin/env python3
"""Smoke test for MetaNexus execution plan and dry-run orchestration.

Sets up a minimal environment with stub executables and smoke reads fixtures
to verify pipeline orchestration without heavy computation or side effects.
"""

import os
import pathlib
import shutil
import subprocess
import tempfile
import unittest

try:
    import pytest
    pytestmark = pytest.mark.smoke
except ImportError:
    pass

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURES_DIR = pathlib.Path(__file__).resolve().parents[1] / "fixtures" / "smoke"
STUB_BIN_DIR = FIXTURES_DIR / "bin"
READS_DIR = FIXTURES_DIR / "reads"


class TestSmokePipeline(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="smoke_out_")
        self.orig_path = os.environ.get("PATH", "")
        # Prepend stub binaries
        self.env = dict(os.environ)
        self.env["PATH"] = f"{STUB_BIN_DIR}:{self.orig_path}"

    def tearDown(self):
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_smoke_execution_plan(self):
        cmd = [
            "bash",
            "./run_metagenomics_pipeline.sh",
            "--plan",
            "--reads-dir", str(READS_DIR),
            "--outdir", self.temp_dir,
        ]
        result = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=self.env,
        )
        self.assertEqual(
            result.returncode, 0,
            f"Smoke plan exited with {result.returncode}\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        self.assertIn("execution plan", result.stdout.lower())
        self.assertIn("01_fastp", result.stdout)

    def test_smoke_dry_run(self):
        cmd = [
            "bash",
            "./run_metagenomics_pipeline.sh",
            "--dry-run",
            "--reads-dir", str(READS_DIR),
            "--outdir", self.temp_dir,
        ]
        result = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=self.env,
        )
        self.assertEqual(
            result.returncode, 0,
            f"Smoke dry-run exited with {result.returncode}\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )
        self.assertIn("dry-run", result.stdout.lower())


if __name__ == "__main__":
    unittest.main()
