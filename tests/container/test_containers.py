#!/usr/bin/env python3
"""Container validation tests for MetaNexus Apptainer (.sif) images."""

import os
import shutil
import subprocess
import unittest
from pathlib import Path

try:
    import pytest
    pytestmark = pytest.mark.container
except ImportError:
    pass

REPO_ROOT = Path(__file__).resolve().parents[2]
SIF_DIR = REPO_ROOT / "containers"

CONTAINER_CHECKS = [
    ("fastp.sif", "fastp --version"),
    ("megahit.sif", "megahit --version"),
    ("kraken2_bracken.sif", "kraken2 --version && which bracken"),
    ("mapping_depth.sif", "bowtie2 --version && samtools --version && which jgi_summarize_bam_contig_depths"),
    ("vamb.sif", "vamb --version"),
    ("metabinner.sif", "which run_metabinner.sh && test -f /opt/metabinner/scripts/Filter_tooshort.py"),
    ("comebin.sif", 'python -c "import torch" && which run_comebin.sh'),
    ("metawrap.sif", "metawrap --version || which metawrap"),
    ("checkm2.sif", "checkm2 --version"),
    ("gtdbtk.sif", "gtdbtk --version"),
    ("coverm.sif", "coverm --version && samtools --version"),
    ("r_base.sif", "Rscript --version"),
]


class TestContainerImages(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = shutil.which("apptainer") or shutil.which("singularity")

    def test_container_engine_available(self):
        self.assertIsNotNone(
            self.engine,
            "Neither apptainer nor singularity found in PATH"
        )

    def test_all_12_containers_present_and_executable(self):
        if not self.engine:
            self.skipTest("Container engine not installed on host")

        failures = []
        for sif_name, check_cmd in CONTAINER_CHECKS:
            sif_path = SIF_DIR / sif_name
            if not sif_path.is_file():
                failures.append(f"MISSING: {sif_name} not found in {SIF_DIR}")
                continue

            res = subprocess.run(
                [self.engine, "exec", "--cleanenv", str(sif_path), "bash", "-c", check_cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if res.returncode != 0:
                failures.append(f"FAILED: {sif_name} running '{check_cmd}' (exit {res.returncode}): {res.stderr.strip()}")

        self.assertEqual(
            len(failures),
            0,
            f"Container verification failed:\n" + "\n".join(failures)
        )


if __name__ == "__main__":
    unittest.main()
