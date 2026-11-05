#!/usr/bin/env bash
# run_container_tests.sh – wrapper for container verification with explicit PASS/FAIL/SKIP outcomes
# It is invoked by `make test-container`.

set -euo pipefail

# Locate apptainer (or singularity)
APPTAINER_CMD=$(command -v apptainer || command -v singularity || true)
if [[ -z "$APPTAINER_CMD" ]]; then
  echo "[SKIP] Apptainer/Singularity not installed – container tests skipped."
  exit 0
fi

SIF_DIR="$(dirname "${0}")/../containers"
if [[ ! -d "$SIF_DIR" ]]; then
  echo "[SKIP] Containers directory not found: $SIF_DIR"
  exit 0
fi

FAIL=0
SKIP=0

check_sif() {
  local sif_name="$1"
  local cmd="$2"
  local sif_path="$SIF_DIR/$sif_name"
  if [[ ! -f "$sif_path" ]]; then
    echo "[SKIP] $sif_name not found"
    SKIP=$((SKIP+1))
    return
  fi
  if "$APPTAINER_CMD" exec --cleanenv "$sif_path" bash -c "$cmd" >/dev/null 2>&1; then
    echo "[PASS] $sif_name ($cmd)"
  else
    echo "[FAIL] $sif_name ($cmd)"
    FAIL=$((FAIL+1))
  fi
}

# List of containers to verify – order does not matter
check_sif fastp.sif            "fastp --version"
check_sif megahit.sif          "megahit --version"
check_sif kraken2_bracken.sif  "kraken2 --version && which bracken"
check_sif mapping_depth.sif    "bowtie2 --version && samtools --version && which jgi_summarize_bam_contig_depths"
check_sif vamb.sif             "vamb --version"
check_sif metabinner.sif       "which run_metabinner.sh && test -f /opt/metabinner/scripts/Filter_tooshort.py"
check_sif comebin.sif          "python -c \"import torch\" && which run_comebin.sh"
check_sif metawrap.sif         "metawrap --version || which metawrap"
check_sif checkm2.sif          "checkm2 --version"
check_sif gtdbtk.sif           "gtdbtk --version"
check_sif coverm.sif           "coverm --version && samtools --version"
check_sif r_base.sif           "Rscript --version"

if (( FAIL > 0 )); then
  echo "[ERROR] $FAIL container(s) failed verification."
  exit 1
else
  echo "[SUCCESS] All existing containers passed verification ($SKIP skipped)."
  exit 0
fi
