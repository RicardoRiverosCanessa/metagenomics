#!/usr/bin/env bash
###############################################################################
# run_metagenomics_pipeline.sh — MetaNexus Shotgun Metagenomics Pipeline
#
# Internal checkpoint IDs/directories are stable for existing results.
# Visible execution order begins with fastp, Kraken2/Bracken, then MEGAHIT.
#   01_fastp / 01_qc
#   03_kraken2 / 03_kraken2 (executed before assembly)
#   02_megahit / 02_assembly
#   04_mapping, 05a_vamb, 05b_metabinner, 05c_comebin
#   06_refinement, 07_reassembly, 08_checkm2, 09_filter_mags
#   10_gtdbtk, 11_coverm, 12_summary
#
# Usage:
#   1. Edit all variables in the CONFIGURATION block below.
#   2. bash run_metagenomics_pipeline.sh [--profile PROFILE] [--dry-run] [--force]
#
# Execution model:
#   Container steps: fastp, megahit, kraken2, bracken, bowtie2, samtools,
#                     jgi_summarize_bam_contig_depths, vamb, metabinner,
#                     comebin, metawrap, checkm2, gtdbtk, coverm, Rscript
#   mamba mode     : all tools retain their original host/mamba execution paths
###############################################################################

set -Eeuo pipefail

early_error_handler() {
    local rc="$1"
    local line="$2"

    trap - ERR
    set +e

    local ts pid
    ts="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "unknown-time")"
    pid="$$"

    printf '[FATAL] %s [%s] Early pipeline failure: exit=%s line=%s command=<redacted-early-command>\n' \
        "${ts}" "${pid}" "${rc}" "${line}" >&2

    exit "${rc}"
}

early_signal_handler() {
    local sig="$1"
    local rc="$2"

    trap - INT TERM HUP ERR
    set +e

    local ts pid
    ts="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "unknown-time")"
    pid="$$"

    printf '[FATAL] %s [%s] Early pipeline interrupted by %s (exit=%s)\n' \
        "${ts}" "${pid}" "${sig}" "${rc}" >&2

    exit "${rc}"
}

trap 'early_error_handler "$?" "$LINENO"' ERR
trap 'early_signal_handler INT 130' INT
trap 'early_signal_handler TERM 143' TERM
trap 'early_signal_handler HUP 129' HUP

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
SUMMARY_GENERATOR="${PIPELINE_SCRIPTS_DIR}/generate_pipeline_summaries.py"
WRITE_ERROR_RECORD="${PIPELINE_SCRIPTS_DIR}/write_error_record.py"

readonly PIPELINE_NAME="MetaNexus"
readonly PIPELINE_VERSION="0.3.2"
readonly PIPELINE_RELEASE_STATUS="development"
PIPELINE_START_TIME_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "unknown-time")"
export PIPELINE_START_TIME_UTC
PIPELINE_PID="$$"
PIPELINE_HOSTNAME="$(hostname 2>/dev/null || uname -n 2>/dev/null || echo "localhost")"
PIPELINE_INVOCATION_RAW=("$0" "$@")
PIPELINE_INVOCATION_COMMAND=""
PIPELINE_START_MONOTONIC=""

_rnd="$(head -c 2 /dev/urandom 2>/dev/null | od -An -tx1 2>/dev/null | tr -d ' \n' || true)"
[[ -z "${_rnd}" ]] && _rnd="$(printf '%04x' "$(( RANDOM % 65536 ))")"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ 2>/dev/null || echo "19700101T000000Z")_${PIPELINE_PID}_${_rnd}"

###############################################################################
#  CONFIGURATION  —  Edit everything in this block before running             #
###############################################################################

# ── Input / Output ────────────────────────────────────────────────────────────
READS_DIR="${READS_DIR:-${RAW_READS_DIR:-}}" # Raw reads directory (compatibility alias)
RAW_READS_DIR="${READS_DIR}"
CLEAN_READS_DIR="${CLEAN_READS_DIR:-}"
ACTIVE_READS_DIR=""
ACTIVE_READS_SOURCE=""
ACTIVE_R1_FILES=()
ACTIVE_R2_FILES=()
ACTIVE_SAMPLE_IDS=()
OUTDIR="${OUTDIR:-./metagenomics_out}"     # All pipeline output goes here
ASSEMBLY_IN="${ASSEMBLY_IN:-}"   # Pre-assembled FASTA when MEGAHIT is disabled.
BAM_DIR_INPUT="${BAM_DIR_INPUT:-}"
BAM_MANIFEST_INPUT="${BAM_MANIFEST_INPUT:-}"
JGI_DEPTH_FILE_INPUT="${JGI_DEPTH_FILE_INPUT:-}"
BAM_DIR_ACTIVE=""
BAM_MANIFEST_ACTIVE=""
JGI_DEPTH_FILE_ACTIVE=""
COMEBIN_BAM_DIR=""
COVERAGE_SOURCE=""
RESUME_CONTROL=resume
FORCE_STAGE_NAMES=()
RESTART_FROM_NAME=""
RESUME_REQUESTED=false
FORCE_REQUESTED=false
RESTART_REQUESTED=false
RESTART_REQUEST_COUNT=0
RUN_MANIFEST_FILE=""
STAGE_STATUS_FILE=""
RUN_STATUS_FILE=""
RESUME_PLAN_FILE=""
BENCHMARK_WRITER="${PIPELINE_SCRIPTS_DIR}/write_benchmark_manifest.py"
WRITE_DISK_SPACE_REPORT="${PIPELINE_SCRIPTS_DIR}/write_disk_space_report.py"
CURRENT_STAGE_METRICS_DIR=""
GNU_TIME_BIN=""
GNU_TIME_AVAILABLE=false
CURRENT_STAGE=""
CURRENT_STAGE_NAME=""
CURRENT_STAGE_ACTION=""
CURRENT_STAGE_LOG=""
CURRENT_STAGE_RUNTIME="unknown"
CURRENT_STAGE_COMMAND=""
CURRENT_COMMAND_SOURCE_FILE=""
CURRENT_COMMAND_SOURCE_LINE=""

LAST_VALID_CHECKPOINT="none"
LAST_VALID_CHECKPOINT_PATH="none"
LAST_VALID_PRODUCT="none"
LAST_VALID_PRODUCT_FINGERPRINT="none"

ERROR_HANDLER_ACTIVE=0
ERROR_ALREADY_RECORDED=0
RUN_STATUS_INITIALIZED=false
PIPELINE_TERMINATION_RECORDED=false
declare -A STAGE_ACTION=()
declare -A STAGE_REASON=()
declare -A STAGE_STATE=()
declare -A GRAPH_CHILDREN=()
declare -A GRAPH_PARENTS=()
declare -A STAGE_FINAL_STATUS=()
declare -A STAGE_STAGE_NAME=()
declare -A STAGE_CHECKPOINT_PATH=()
declare -A STAGE_REQUESTED_ACTION=()
declare -A STAGE_APPLICABILITY=()
declare -A STAGE_REQUIRED=()
declare -A STAGE_CONFIGURED_ENABLED=()
declare -A STAGE_CHECKPOINT_STATUS=()
declare -A STAGE_CHECKPOINT_REASON=()
declare -A STAGE_INPUT_STATUS=()
declare -A STAGE_RUNTIME_STATUS=()
declare -A STAGE_DATABASE_STATUS=()
declare -A STAGE_OUTPUT_STATUS=()
declare -A STAGE_VALIDATION_STATUS=()
declare -A STAGE_INVALIDATED_BY=()
declare -A STAGE_BLOCKING_STAGE=()
declare -A STAGE_REASON_CODE=()
declare -A STAGE_REASON_DETAIL=()
declare -A STAGE_DIAGNOSTIC_CODES=()
declare -A STAGE_DEPENDENCIES=()
declare -A STAGE_DOWNSTREAM=()
declare -A STAGE_PRIMARY_OUTPUT=()
declare -A STAGE_PRIMARY_LOG=()
declare -A STAGE_PRECONDITIONS=()
# shellcheck disable=SC2034
declare -A STAGE_START_TIME=()
# shellcheck disable=SC2034
declare -A STAGE_END_TIME=()
# shellcheck disable=SC2034
declare -A STAGE_EXIT_CODE=()
declare -A STAGE_CHECKPOINT_WRITTEN=()
declare -A STAGE_FROZEN_EXTERNAL=()
declare -A PREPARED_STAGES=()
declare -A ID_VERSION=() ID_VERSION_SOURCE=() ID_EXECUTION_MODE=()
declare -A ID_RUNTIME_NAME=() ID_RUNTIME_URI=() ID_RUNTIME_PATH=()
declare -A ID_RUNTIME_SIZE=() ID_RUNTIME_MTIME=() ID_RUNTIME_METHOD=()
declare -A ID_RUNTIME_FINGERPRINT=() ID_COMMAND=() ID_STATUS=()
declare -A ID_DATABASE_RELEASE=() ID_DATABASE_VERSION_SOURCE=()
declare -A ID_DATABASE_PATH=() ID_DATABASE_SIZE=() ID_DATABASE_MTIME=()
declare -A ID_DATABASE_METHOD=() ID_DATABASE_FINGERPRINT=() ID_DATABASE_STATUS=()
declare -A SIF_IDENTITY_CACHE=() MAMBA_IDENTITY_CACHE=() VERSION_CACHE=() DB_IDENTITY_CACHE=() DB_SNAPSHOT=()
declare -A SIF_OBJECT_CACHE=() MAMBA_OBJECT_CACHE=()
SOFTWARE_INVENTORY_ROWS=()
DATABASE_INVENTORY_ROWS=()
INVENTORY_UNKNOWN_COUNT=0
INVENTORY_UNAVAILABLE_COUNT=0

# ── Read File Patterns ────────────────────────────────────────────────────────
READ_SUFFIX_R1="_R1.fastq.gz"
READ_SUFFIX_R2="_R2.fastq.gz"

# ── Compute Resources ─────────────────────────────────────────────────────────
# Use "auto" to detect from hardware, a positive integer to set manually, or
# export THREADS / MEMORY as environment variables before running.
THREADS="${THREADS:-auto}"   # logical CPU cores to use (or "auto")
MEMORY="${MEMORY:-auto}"     # GB RAM cap for MEGAHIT / MetaWRAP (or "auto")

# ── Resource Safety Margins ───────────────────────────────────────────────────
CPU_RESERVE=2    # logical cores always kept free for the OS
MEM_FRACTION=85  # maximum % of available RAM the pipeline may use

# ── Tool Binaries Directory ───────────────────────────────────────────────────
# Compiled binaries accessed directly. Leave empty ("") if already on $PATH.
BIN_DIR=""

# ── Container / Environment Execution Engine ──────────────────────────────────
# EXEC_MODE: "apptainer" — run isolated steps via Apptainer/Singularity .sif
#            "mamba"     — fall back to  mamba run -n <env>  (original behaviour)
EXEC_MODE="${EXEC_MODE:-apptainer}"

# ── Apptainer/Singularity Container Images (.sif) ────────────────────────────
# SIF_DIR  : local directory where .sif files are stored / downloaded.
# AUTO_PULL_SIF : if true, missing images are pulled from their URI automatically.
SIF_DIR="${SIF_DIR:-./containers}"
AUTO_PULL_SIF=true

# Per-tool image paths and their BioContainers / Docker Hub source URIs.
# Override any path to point to a locally-built or pre-cached .sif.
SIF_FASTP="${SIF_DIR}/fastp.sif"
SIF_FASTP_URI="docker://quay.io/biocontainers/fastp:0.23.4--h125f33a_5"

SIF_MEGAHIT="${SIF_DIR}/megahit.sif"
SIF_MEGAHIT_URI="docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_4"

SIF_KRAKEN2="${SIF_DIR}/kraken2_bracken.sif"
SIF_KRAKEN2_URI="docker://staphb/kraken2:2.1.3"

SIF_MAPPING="${SIF_DIR}/mapping_depth.sif"
SIF_MAPPING_URI="docker://staphb/bowtie2:2.5.4"

SIF_VAMB="${SIF_DIR}/vamb.sif"
SIF_VAMB_URI="docker://ghcr.io/RasmussenLab/vamb:latest"
VAMB_BINSPLIT_SEPARATOR="${VAMB_BINSPLIT_SEPARATOR:-}"
VAMB_BINSPLIT_SEPARATOR_WAS_SET=false
# Vamb --minfasta is a minimum output-bin size, independent of MIN_CONTIG_LEN.
VAMB_MINFASTA="${VAMB_MINFASTA:-1000}"
VAMB_VERSION=""
VAMB_INTERNAL_MIN_CONTIG_LEN=""
VAMB_DEFAULT_THREADS=""

SIF_METABINNER="${SIF_DIR}/metabinner.sif"
SIF_METABINNER_URI="docker://quay.io/biocontainers/metabinner:1.4.4--pyhdfd78af_0"

SIF_COMEBIN="${SIF_DIR}/comebin.sif"
SIF_COMEBIN_URI="docker://quay.io/biocontainers/comebin:1.0.4--pyhdfd78af_0"

SIF_METAWRAP="${SIF_DIR}/metawrap.sif"
SIF_METAWRAP_URI="docker://quay.io/biocontainers/metawrap-mg:1.3.2--hdfd78af_1"

SIF_CHECKM2="${SIF_DIR}/checkm2.sif"
SIF_CHECKM2_URI="docker://quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0"

SIF_GTDBTK="${SIF_DIR}/gtdbtk.sif"
SIF_GTDBTK_URI="docker://quay.io/biocontainers/gtdbtk:2.4.0--pyhdfd78af_0"

SIF_COVERM="${SIF_DIR}/coverm.sif"
SIF_COVERM_URI="docker://quay.io/biocontainers/coverm:0.6.1--h9f5acd7_1"

# ── Mamba Environment Names (used when EXEC_MODE="mamba") ────────────────────
ENV_VAMB="vamb"
ENV_METABINNER="metabinner_env"
ENV_COMEBIN="comebin_env"
ENV_METAWRAP="metawrap_env"
ENV_CHECKM2="checkm2_env"
ENV_GTDBTK="gtdbtk_env"

# ── Database Paths ────────────────────────────────────────────────────────────
KRAKEN2_DB="${KRAKEN2_DB:-}"
CHECKM2_DB="${CHECKM2_DB:-}"
GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-}"

# ── MetaBinner Script Paths ───────────────────────────────────────────────────
# When EXEC_MODE="apptainer" the scripts live inside the container at
# /opt/metabinner/scripts; they do NOT need to be present on the host.
# When EXEC_MODE="mamba" set this to the actual host path.
#
# Smart default: if the variable has not been set (or still holds the
# placeholder), the value is resolved later in main() after EXEC_MODE is known.
# You can always override it by exporting before running:
#   export METABINNER_SCRIPTS_DIR=/my/host/path
_METABINNER_SCRIPTS_DIR_ENV="${METABINNER_SCRIPTS_DIR:-}"
METABINNER_SCRIPTS_DIR=""
METABINNER_RUNNER="${METABINNER_RUNNER:-run_metabinner.sh}"

# ── COMEBin Script Path ───────────────────────────────────────────────────────
# Same logic: run_comebin.sh lives inside the container when EXEC_MODE="apptainer".
COMEBIN_SCRIPT="${COMEBIN_SCRIPT:-run_comebin.sh}"
COMEBIN_VERSION="1.1.0"
COMEBIN_DEVICE="${COMEBIN_DEVICE:-cpu}"
COMEBIN_SEED="${COMEBIN_SEED:-42}"
COMEBIN_DEVICE_CLI_SET=false
COMEBIN_PYTORCH_VERSION="uninspected"
COMEBIN_PYTORCH_CUDA_BUILD="uninspected"
COMEBIN_PYTORCH_CUDA_AVAILABLE="uninspected"
COMEBIN_CUDA_DEVICE_COUNT="uninspected"
COMEBIN_APPTAINER_NV=false

# ── Assembly Parameters ───────────────────────────────────────────────────────
MIN_CONTIG_LEN=1000    # bp — contigs shorter than this are discarded

# ── Kraken2 + Bracken Parameters ─────────────────────────────────────────────
BRACKEN_READ_LEN=150
BRACKEN_TAX_LEVEL="S"   # S=Species G=Genus F=Family O=Order C=Class P=Phylum
BRACKEN_THRESHOLD=10

# ── MAG Quality Filtering (CheckM2) ──────────────────────────────────────────
MIN_COMPLETENESS="${MIN_COMPLETENESS:-70}"
MAX_CONTAMINATION="${MAX_CONTAMINATION:-10}"
MAG_FILTERING_APPLICABLE=false

# ── MetaWrap Refinement Thresholds ────────────────────────────────────────────
METAWRAP_COMPLETENESS=50    # passed to bin_refinement -c
METAWRAP_CONTAMINATION=10   # passed to bin_refinement -x

# ── CoverM Parameters ─────────────────────────────────────────────────────────
COVERM_METHODS="relative_abundance"   # leave empty for CoverM built-in default
COVERM_MIN_COVERED_FRACTION="0.0"

# ── Pipeline Step Flags ───────────────────────────────────────────────────────
RUN_FASTP="${RUN_FASTP:-true}"                # Step 01
RUN_MEGAHIT="${RUN_MEGAHIT:-true}"            # Internal checkpoint 02_megahit
RUN_KRAKEN2="${RUN_KRAKEN2:-true}"            # Internal checkpoint 03_kraken2
RUN_BRACKEN="${RUN_BRACKEN:-true}"
RUN_REASSEMBLY="${RUN_REASSEMBLY:-false}"     # Step 07
RUN_GTDBTK="${RUN_GTDBTK:-true}"              # Step 10
RUN_MAPPING="${RUN_MAPPING:-true}"
RUN_BINNING="${RUN_BINNING:-true}"
_REFINEMENT_MODE_WAS_SET="${REFINEMENT_MODE+x}"
if [[ -z "${_REFINEMENT_MODE_WAS_SET}" && -n "${RUN_REFINEMENT+x}" ]]; then
    case "${RUN_REFINEMENT}" in
        true) REFINEMENT_MODE=auto ;;
        false) REFINEMENT_MODE=never ;;
        *) echo "[ERROR] RUN_REFINEMENT must be true or false." >&2; exit 1 ;;
    esac
else
    REFINEMENT_MODE="${REFINEMENT_MODE:-auto}"
fi
unset RUN_REFINEMENT
SELECTED_BINNERS="${SELECTED_BINNERS:-vamb,metabinner,comebin}"
BINS_DIR_INPUT="${BINS_DIR_INPUT:-}"
BINS_SOURCE="${BINS_SOURCE:-auto}"
RUN_REFINEMENT_EFFECTIVE=false
SELECTED_BINNER_ARRAY=()
BINS_SOURCE_RESOLVED=""
BINS_FOR_REASSEMBLY=""
BINS_FOR_QC=""
RUN_CHECKM2="${RUN_CHECKM2:-true}"
RUN_MAG_FILTERING="${RUN_MAG_FILTERING:-true}"
RUN_COVERM="${RUN_COVERM:-true}"
RUN_SUMMARY="${RUN_SUMMARY:-true}"
RUN_TAXONOMY_SUMMARY="${RUN_TAXONOMY_SUMMARY:-true}"
RUN_MAG_SUMMARY="${RUN_MAG_SUMMARY:-true}"
RUN_RUN_OVERVIEW="${RUN_RUN_OVERVIEW:-true}"
STOP_AFTER="${STOP_AFTER:-}"
REQUESTED_PROFILE="${PROFILE:-full}"
RESULTING_PROFILE="${REQUESTED_PROFILE}"
NEEDS_READS=false
NEEDS_ASSEMBLY=false
NEEDS_MAGS=false

# ── Execution Flags ───────────────────────────────────────────────────────────
DRY_RUN="${DRY_RUN:-false}"                   # Print commands without running them
FORCE_RERUN="${FORCE_RERUN:-false}"           # Ignore checkpoints and rerun every step
PLAN_MODE="${PLAN_MODE:-false}"
PLAN_FORMAT="text"
PLAN_DETAILS="full"
STRICT_PLAN=false
PLAN_IS_FROZEN=false
COLLECTING_PLAN=false

# ── Temporary Workspace ────────────────────────────────────────────────────────
# Captured before any export; used as the environment-tier candidate.
TMPDIR_ENV_ORIGINAL="${TMPDIR:-}"
# CLI option value for --tmpdir (empty means not provided).
PIPELINE_TMPDIR_CLI=""
# Resolved tier: cli | environment | fallback
PIPELINE_TMPDIR_SOURCE=""
# Canonical base directory chosen after precedence resolution.
PIPELINE_TMPDIR_BASE=""
# ${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline
PIPELINE_TMP_PARENT=""
# ${PIPELINE_TMP_PARENT}/${RUN_ID}  — set after plan is frozen
PIPELINE_RUN_TMPDIR=""
# ${PIPELINE_RUN_TMPDIR}/reassembly  — set after plan is frozen
REASSEMBLY_TMPDIR=""
# ${PIPELINE_RUN_TMPDIR}/.pipeline_tmp_owner
PIPELINE_TMP_MARKER=""
# false = remove on controlled success/failure; true = always preserve
KEEP_TEMPORARY=false
# Concatenation strategy name and version (frozen at preflight, not changed at runtime)
REASSEMBLY_TEMP_READS_STRATEGY="gzip-member-concat"
REASSEMBLY_TEMP_READS_STRATEGY_VERSION="1"
# Validation name and version
REASSEMBLY_TEMP_READS_VALIDATION="gzip-integrity-and-source-manifest"
REASSEMBLY_TEMP_READS_VALIDATION_VERSION="1"
# Safety margin multiplier applied to estimated space requirement (≥1.0)
REASSEMBLY_CONCAT_SPACE_FACTOR="1.10"
# Minimum inodes required in the pipeline workspace (pipeline files only)
REQUIRED_FREE_INODES=16
# ── Disk Space Assessment Configuration ───────────────────────────────────────
DISK_WARNING_HEADROOM_FACTOR="1.25"
FASTP_OUTPUT_FACTOR="1.10"
MEGAHIT_WORK_FACTOR="2.50"
KRAKEN_OUTPUT_FACTOR="1.00"
MAPPING_OUTPUT_FACTOR="1.50"
DISK_PRECHECK_STAGES=(
  01_fastp
  02_megahit
  03_kraken2
  04_mapping
  05a_vamb
  05b_metabinner
  05c_comebin
  06_refinement
  07_reassembly
  08_checkm2
  10_gtdbtk
  11_coverm
)

###############################################################################
#  END OF CONFIGURATION                                                       #
###############################################################################

# ── CLI Parsing ───────────────────────────────────────────────────────────────
_show_help() {
    cat <<HELP
${PIPELINE_NAME} ${PIPELINE_VERSION} (${PIPELINE_RELEASE_STATUS}) — Shotgun Metagenomics Pipeline

Usage:
  bash run_metagenomics_pipeline.sh [OPTIONS]

Identity & Inspection:
  -V, --version                Show pipeline identity and exit
      --versions               Inspect tools, runtimes, and databases (read-only)
      --versions-all           Inspect all pipeline components across all stages
      --versions-deep          Compute cryptographic fingerprints during inspection
      --versions-format <text|json>
                               Format inspection report as text or json (default: text)
      --strict-versions        Exit with code 1 if any component is unavailable

Options:
  -r, --reads-dir <directory>  Directory containing raw paired-end reads
      --clean-reads-dir <dir>  Directory containing externally cleaned reads
      --kraken2-db <directory> Kraken2 database directory
      --checkm2-db <file>      CheckM2 database file
      --gtdbtk-data-path <dir> GTDB-Tk reference data directory
  -a, --assembly-in <file>     External assembly FASTA; requires --no-megahit
      --bam-dir <directory>    External <sample>.sorted.bam files (use --no-mapping)
      --bam-manifest <file>    Tab-separated sample_id, bam, index manifest
      --jgi-depth-file <file>  External JGI depth matrix (use --no-mapping)
  -o, --outdir <directory>     Output directory (default: ./metagenomics_out)
      --metabinner-scripts-dir <directory>
                               MetaBinner scripts path (host path for mamba;
                               defaults to /opt/metabinner/scripts for Apptainer)
  -t, --threads <int>          Override THREADS (overrides "auto" and env var)
  -m, --memory  <int>          Override MEMORY in GB (overrides "auto" and env var)
  -e, --engine  <apptainer|mamba>
                               Select execution engine for pipeline tools.
                               Default: apptainer (uses .sif images from SIF_DIR).
                               Fallback: mamba (uses  mamba run -n <env>).
  --fastp / --no-fastp         Enable/disable fastp (default: enabled)
  --profile <full|taxonomy|assembly|custom>
                               Select both branches or individual stages
  --taxonomy-only              Alias for --profile taxonomy
  --assembly-only              Alias for --profile assembly
  --megahit / --no-megahit     Enable/disable MEGAHIT (default: enabled)
  --kraken2 / --no-kraken2     Enable/disable Kraken2 (default: enabled)
  --bracken / --no-bracken     Enable/disable Bracken abundance estimation (default: enabled)
  --reassembly / --no-reassembly
                               Enable/disable MetaWrap reassembly (default: disabled)
  --tmpdir DIR                 Base directory for temporary pipeline workspaces.
                               The pipeline creates an isolated RUN_ID directory
                               below this location. Priority: CLI > TMPDIR env >
                               fallback (OUTDIR/.tmp). Must be absolute or
                               resolvable to an absolute path.
  --keep-temporary             Preserve the current execution's temporary workspace
                               after success or failure.
  --no-keep-temporary          Remove the current execution's temporary workspace
                               during controlled cleanup (default).
  --gtdbtk / --no-gtdbtk       Enable/disable GTDB-Tk (default: enabled)
  --mapping / --no-mapping     Enable/disable mapping (default: enabled)
  --binning / --no-binning     Enable/disable binning (default: enabled)
  --binners <list>             Comma-separated: vamb,metabinner,comebin
  --comebin-device DEVICE      Training device passed to COMEBin through
                               run_comebin.sh -d. Accepted values: cpu,
                               cuda, cuda:N. Default: cpu.
  --comebin-seed INTEGER       COMEBin random seed (default: 42)
  --min-completeness VALUE     Minimum CheckM2 completeness required to retain
                               a MAG. Inclusive threshold. Default: 70.
  --max-contamination VALUE    Maximum CheckM2 contamination allowed to retain
                               a MAG. Inclusive threshold. Default: 10.
  --vamb-binsplit-separator CHAR
                               Advanced Vamb binsplitting; disabled by default
  --refinement-mode <mode>     Consensus refinement: auto, always, or never
  --refinement                 Alias for --refinement-mode always
  --no-refinement              Alias for --refinement-mode never
  --bins-dir <directory>      Directory containing one external final bin collection
  --bins-source <source>       auto, vamb, metabinner, comebin, refined, or external
  --checkm2 / --no-checkm2     Enable/disable CheckM2 (default: enabled)
  --mag-filtering / --no-mag-filtering
                               Enable/disable MAG filtering (default: enabled)
  --coverm / --no-coverm       Enable/disable CoverM (default: enabled)
  --summary / --no-summary     Enable/disable all summaries (default: enabled)
  --taxonomy-summary / --no-taxonomy-summary
                               Enable/disable read-based taxonomy summary (default: enabled)
  --mag-summary / --no-mag-summary
                               Enable/disable MAG quality/taxonomy/abundance summary (default: enabled)
  --run-overview / --no-run-overview
                               Enable/disable technical run overview (default: enabled)
  --stop-after <stage>         Stop after: mapping, binning, refinement,
                               reassembly, checkm2, mag-filtering, gtdbtk,
                               coverm, taxonomy-summary, mag-summary, or run-overview
  --plan                       Display the execution plan and decisions without executing analysis
  --plan-format <text|json>    Execution plan format (default: text; requires --plan)
  --plan-details <summary|full>
                               Execution plan detail level (default: full; requires --plan)
  --strict-plan                Exit 1 if any required stage is blocked (requires --plan)
  --dry-run / --no-dry-run     Validate and print the frozen plan without writing results
  --force / --no-force         Force all active stages / cancel prior --force
  --resume                     Reuse compatible checkpoints (default)
  --force-stage <stage>        Force a stage and run its effective descendants (canonical names:
                               fastp, kraken2, megahit, mapping, vamb, metabinner, comebin,
                               refinement, reassembly, checkm2, mag-filtering, gtdbtk, coverm,
                               taxonomy-summary, mag-summary, run-overview)
  --restart-from <stage>       Destructively restart one effective chain
  -h, --help                   Show this help message and exit

Example:
  bash run_metagenomics_pipeline.sh \\
    --reads-dir /data/reads \\
    --kraken2-db /databases/kraken2 \\
    --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \\
    --gtdbtk-data-path /databases/gtdbtk \\
    --outdir ./metagenomics_out

External assembly (only when downstream stages need one):
  bash run_metagenomics_pipeline.sh --no-megahit \\
    --assembly-in /data/assembly/final.contigs.fa --reads-dir /data/reads

Resource auto-detection:
  When THREADS or MEMORY are set to "auto" (the default), the pipeline
  queries the host hardware and applies safety margins:
    CPU_RESERVE  (default 2)  — cores reserved for the OS
    MEM_FRACTION (default 85) — max % of available RAM to use
  Override any of these in the CONFIGURATION block or via environment variables.

Container engine (EXEC_MODE=apptainer):
  .sif images are resolved from SIF_DIR (default: ./containers/).
  If a .sif is missing and AUTO_PULL_SIF=true, it is pulled automatically
  via  apptainer pull  (or  singularity pull  if apptainer is unavailable).
  Set EXEC_MODE=mamba to bypass containers entirely and use conda environments.
HELP
    exit 0
}

_show_plan_help() {
    cat <<PLAN_HELP
${PIPELINE_NAME} ${PIPELINE_VERSION} (${PIPELINE_RELEASE_STATUS}) — Execution Plan Inspection

Usage:
  ./run_metagenomics_pipeline.sh --plan [OPTIONS]

Options:
  --plan                       Inspect and display the canonical execution plan in a
                               strictly read-only manner without executing analysis.
  --plan-format <text|json>    Output format: text (default) or json.
  --plan-details <summary|full>
                               Detail level: full (default) or summary.
  --strict-plan                Exit with code 1 if any required stage is blocked.
  --profile <profile>          Target pipeline profile (full, taxonomy, assembly, custom).
  --force-stage <stage>        Force a stage and its effective descendants.
  --restart-from <stage>       Restart from a boundary stage.
  -h, --help                   Show this help message and exit.

Examples:
  ./run_metagenomics_pipeline.sh --plan
  ./run_metagenomics_pipeline.sh --plan --plan-format json
  ./run_metagenomics_pipeline.sh --plan --plan-details summary
  ./run_metagenomics_pipeline.sh --plan --strict-plan
PLAN_HELP
    exit 0
}

_show_versions_help() {
    cat <<VERSIONS_HELP
${PIPELINE_NAME} ${PIPELINE_VERSION} (${PIPELINE_RELEASE_STATUS}) — Version & Runtime Inspection

Usage:
  ./run_metagenomics_pipeline.sh --versions [OPTIONS]

Options:
  --versions                   Inspect pipeline tools, runtimes, and databases in a
                               strictly read-only manner without executing analysis.
  --versions-all               Inspect all tools, runtimes, and reference databases
                               across all stages, regardless of profile or disabled stages.
  --versions-deep              Perform cryptographic fingerprinting of runtime images
                               and reference databases (except deep GTDB-Tk).
  --versions-format <text|json>
                               Output format: text (default) or json.
  --strict-versions            Exit with non-zero code (1) if any inspected tool, runtime,
                               or database is unavailable or unknown.
  -e, --engine <apptainer|mamba>
                               Execution engine to inspect (default: apptainer).
  --profile <full|taxonomy|assembly|custom>
                               Target profile whose stages should be inspected (default: full).
  --taxonomy-only              Alias for --profile taxonomy.
  --assembly-only              Alias for --profile assembly.
  --kraken2-db <dir>           Path to Kraken2/Bracken database directory.
  --checkm2-db <file>          Path to CheckM2 diamond database file.
  --gtdbtk-data-path <dir>     Path to GTDB-Tk reference data directory.
  --sif-dir <dir>              Directory containing .sif images (default: ./containers).
  --timeout <seconds>          Per-command execution timeout (default: 20s).
  -h, --help                   Show this help message and exit.

Examples:
  ./run_metagenomics_pipeline.sh --versions
  ./run_metagenomics_pipeline.sh --versions --versions-format json
  ./run_metagenomics_pipeline.sh --versions --profile taxonomy
  ./run_metagenomics_pipeline.sh --versions --versions-all --versions-deep
  ./run_metagenomics_pipeline.sh --versions --strict-versions
VERSIONS_HELP
    exit 0
}

# shellcheck disable=SC2034
_handle_versions_mode() {
    local -a raw_args=("$@")
    local req_profile="full"
    local run_megahit=true
    local run_fastp=true
    local run_kraken2=true
    local run_mapping=true
    local run_binning=true
    local run_reassembly=false
    local run_checkm2=true
    local run_mag_filtering=true
    local run_gtdbtk=true
    local run_coverm=true
    local run_summary=true
    local run_taxonomy_summary=true
    local run_mag_summary=true
    local run_run_overview=true
    local selected_binners="${SELECTED_BINNERS:-vamb,metabinner,comebin}"
    local refinement_mode="auto"
    local exec_mode="${EXEC_MODE:-apptainer}"
    local sif_dir="${SIF_DIR:-./containers}"
    local kraken2_db="${KRAKEN2_DB:-}"
    local checkm2_db="${CHECKM2_DB:-}"
    local gtdbtk_data_path="${GTDBTK_DATA_PATH:-}"
    local versions_all=false
    local versions_deep=false
    local versions_format="text"
    local strict_versions=false
    local timeout_val=""

    # First pass: check for profile
    local idx=0
    while (( idx < ${#raw_args[@]} )); do
        local arg="${raw_args[idx]}"
        case "${arg}" in
            --profile)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    req_profile="${raw_args[idx+1]}"
                    idx=$((idx + 2))
                else
                    echo "[ERROR] --profile requires a value (full, taxonomy, assembly, custom)." >&2
                    exit 1
                fi
                ;;
            --taxonomy-only) req_profile="taxonomy"; idx=$((idx + 1)) ;;
            --assembly-only) req_profile="assembly"; idx=$((idx + 1)) ;;
            *) idx=$((idx + 1)) ;;
        esac
    done

    # Apply profile defaults
    case "${req_profile}" in
        full)
            run_megahit=true; run_kraken2=true; run_mapping=true; run_binning=true
            refinement_mode="auto"; run_checkm2=true; run_mag_filtering=true
            run_coverm=true; run_summary=true; run_taxonomy_summary=true
            run_mag_summary=true; run_run_overview=true ;;
        taxonomy)
            run_kraken2=true; run_megahit=false; run_mapping=false; run_binning=false
            refinement_mode="never"; run_reassembly=false; run_checkm2=false
            run_mag_filtering=false; run_gtdbtk=false; run_coverm=false
            run_summary=true; run_taxonomy_summary=true; run_mag_summary=false
            run_run_overview=true ;;
        assembly)
            run_kraken2=false; run_megahit=true; run_mapping=true; run_binning=true
            refinement_mode="auto"; run_checkm2=true; run_mag_filtering=true
            run_coverm=true; run_summary=true; run_taxonomy_summary=false
            run_mag_summary=true; run_run_overview=true ;;
        custom) ;;
        *)
            echo "[ERROR] Invalid profile '${req_profile}' (expected full, taxonomy, assembly, or custom)." >&2
            exit 1 ;;
    esac

    # Second pass: parse overrides and versions flags
    idx=0
    while (( idx < ${#raw_args[@]} )); do
        local arg="${raw_args[idx]}"
        case "${arg}" in
            --versions) idx=$((idx + 1)) ;;
            --versions-all) versions_all=true; idx=$((idx + 1)) ;;
            --versions-deep) versions_deep=true; idx=$((idx + 1)) ;;
            --versions-format)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    local fmt_val="${raw_args[idx+1]}"
                    case "${fmt_val}" in
                        text|json) versions_format="${fmt_val}"; idx=$((idx + 2)) ;;
                        *) echo "[ERROR] Invalid format '${fmt_val}' for --versions-format (expected text or json)." >&2; exit 2 ;;
                    esac
                else
                    echo "[ERROR] --versions-format requires a value (text or json)." >&2
                    exit 2
                fi
                ;;
            --strict-versions) strict_versions=true; idx=$((idx + 1)) ;;
            --timeout)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    timeout_val="${raw_args[idx+1]}"
                    idx=$((idx + 2))
                else
                    echo "[ERROR] --timeout requires a value." >&2
                    exit 1
                fi
                ;;
            -e|--engine)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    local eng_val="${raw_args[idx+1]}"
                    case "${eng_val}" in
                        apptainer|mamba) exec_mode="${eng_val}"; idx=$((idx + 2)) ;;
                        *) echo "[ERROR] --engine must be 'apptainer' or 'mamba'." >&2; exit 1 ;;
                    esac
                else
                    echo "[ERROR] --engine requires a value (apptainer or mamba)." >&2
                    exit 1
                fi
                ;;
            --sif-dir)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    sif_dir="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    echo "[ERROR] --sif-dir requires a directory path." >&2
                    exit 1
                fi
                ;;
            --kraken2-db)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    kraken2_db="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    echo "[ERROR] --kraken2-db requires a directory value." >&2
                    exit 1
                fi
                ;;
            --checkm2-db)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    checkm2_db="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    echo "[ERROR] --checkm2-db requires a file value." >&2
                    exit 1
                fi
                ;;
            --gtdbtk-data-path)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    gtdbtk_data_path="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    echo "[ERROR] --gtdbtk-data-path requires a directory value." >&2
                    exit 1
                fi
                ;;
            --fastp) run_fastp=true; idx=$((idx + 1)) ;;
            --no-fastp) run_fastp=false; idx=$((idx + 1)) ;;
            --megahit) run_megahit=true; idx=$((idx + 1)) ;;
            --no-megahit) run_megahit=false; idx=$((idx + 1)) ;;
            --kraken2) run_kraken2=true; idx=$((idx + 1)) ;;
            --no-kraken2) run_kraken2=false; idx=$((idx + 1)) ;;
            --mapping) run_mapping=true; idx=$((idx + 1)) ;;
            --no-mapping) run_mapping=false; idx=$((idx + 1)) ;;
            --binning) run_binning=true; idx=$((idx + 1)) ;;
            --no-binning) run_binning=false; idx=$((idx + 1)) ;;
            --binners)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    selected_binners="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    idx=$((idx + 1))
                fi
                ;;
            --refinement-mode)
                if (( idx + 1 < ${#raw_args[@]} )); then
                    refinement_mode="${raw_args[idx+1]}"; idx=$((idx + 2))
                else
                    idx=$((idx + 1))
                fi
                ;;
            --refinement) refinement_mode="always"; idx=$((idx + 1)) ;;
            --no-refinement) refinement_mode="never"; idx=$((idx + 1)) ;;
            --reassembly) run_reassembly=true; idx=$((idx + 1)) ;;
            --no-reassembly) run_reassembly=false; idx=$((idx + 1)) ;;
            --checkm2) run_checkm2=true; idx=$((idx + 1)) ;;
            --no-checkm2) run_checkm2=false; idx=$((idx + 1)) ;;
            --mag-filtering) run_mag_filtering=true; idx=$((idx + 1)) ;;
            --no-mag-filtering) run_mag_filtering=false; idx=$((idx + 1)) ;;
            --gtdbtk) run_gtdbtk=true; idx=$((idx + 1)) ;;
            --no-gtdbtk) run_gtdbtk=false; idx=$((idx + 1)) ;;
            --coverm) run_coverm=true; idx=$((idx + 1)) ;;
            --no-coverm) run_coverm=false; idx=$((idx + 1)) ;;
            --summary) run_summary=true; idx=$((idx + 1)) ;;
            --no-summary) run_summary=false; idx=$((idx + 1)) ;;
            --profile) idx=$((idx + 2)) ;;
            --taxonomy-only|--assembly-only) idx=$((idx + 1)) ;;
            -r|--reads-dir|--clean-reads-dir|-a|--assembly-in|--bam-dir|--bam-manifest|\
            --jgi-depth-file|--force-stage|--restart-from|--bins-dir|--bins-source|\
            --vamb-binsplit-separator|--comebin-device|--comebin-seed|--min-completeness|\
            --max-contamination|--tmpdir|-o|--outdir|--metabinner-scripts-dir|-t|--threads|\
            -m|--memory|--stop-after)
                idx=$((idx + 2)) ;;
            --dry-run|--no-dry-run|--force|--no-force|--resume|--keep-temporary|--no-keep-temporary|\
            --bracken|--no-bracken|--taxonomy-summary|--no-taxonomy-summary|--mag-summary|\
            --no-mag-summary|--run-overview|--no-run-overview)
                idx=$((idx + 1)) ;;
            *)
                idx=$((idx + 1)) ;;
        esac
    done

    # Compute active stages
    local -a active_stages=()
    if [[ "${versions_all}" == true ]]; then
        active_stages=(01_fastp 03_kraken2 02_megahit 04_mapping 05a_vamb 05b_metabinner 05c_comebin 06_refinement 07_reassembly 08_checkm2 09_filter_mags 10_gtdbtk 11_coverm 12a_taxonomy_summary 12b_mag_summary 12c_run_overview)
    else
        [[ "${run_fastp}" == true ]] && active_stages+=(01_fastp)
        [[ "${run_kraken2}" == true ]] && active_stages+=(03_kraken2)
        [[ "${run_megahit}" == true ]] && active_stages+=(02_megahit)
        [[ "${run_mapping}" == true ]] && active_stages+=(04_mapping)
        if [[ "${run_binning}" == true ]]; then
            IFS=',' read -r -a b_list <<< "${selected_binners}"
            for b in "${b_list[@]}"; do
                case "${b}" in
                    vamb) active_stages+=(05a_vamb) ;;
                    metabinner) active_stages+=(05b_metabinner) ;;
                    comebin) active_stages+=(05c_comebin) ;;
                esac
            done
        fi
        if [[ "${refinement_mode}" == auto || "${refinement_mode}" == always ]]; then
            active_stages+=(06_refinement)
        fi
        [[ "${run_reassembly}" == true ]] && active_stages+=(07_reassembly)
        [[ "${run_checkm2}" == true ]] && active_stages+=(08_checkm2)
        [[ "${run_mag_filtering}" == true ]] && active_stages+=(09_filter_mags)
        [[ "${run_gtdbtk}" == true ]] && active_stages+=(10_gtdbtk)
        [[ "${run_coverm}" == true ]] && active_stages+=(11_coverm)
    fi

    local stages_joined
    stages_joined="$(IFS=,; echo "${active_stages[*]}")"
    [[ -z "${stages_joined}" ]] && stages_joined="none"

    local -a py_cmd=(
        python3 "${PIPELINE_SCRIPTS_DIR}/inspect_versions.py"
        --pipeline-name "${PIPELINE_NAME}"
        --pipeline-version "${PIPELINE_VERSION}"
        --pipeline-status "${PIPELINE_RELEASE_STATUS}"
        --engine "${exec_mode}"
        --sif-dir "${sif_dir}"
        --stages "${stages_joined}"
        --format "${versions_format}"
    )
    [[ "${versions_all}" == true ]] && py_cmd+=(--versions-all)
    [[ "${versions_deep}" == true ]] && py_cmd+=(--versions-deep)
    [[ "${strict_versions}" == true ]] && py_cmd+=(--strict-versions)
    [[ -n "${kraken2_db}" ]] && py_cmd+=(--kraken2-db "${kraken2_db}")
    [[ -n "${checkm2_db}" ]] && py_cmd+=(--checkm2-db "${checkm2_db}")
    [[ -n "${gtdbtk_data_path}" ]] && py_cmd+=(--gtdbtk-data-path "${gtdbtk_data_path}")
    [[ -n "${timeout_val}" ]] && py_cmd+=(--timeout "${timeout_val}")

    trap - ERR
    set +e
    "${py_cmd[@]}"
    local rc=$?
    exit "${rc}"
}

# ── Terminal Pre-parsing (--version, -V, --versions, --plan, --help) ─────────
_has_version=false
_has_versions=false
_has_plan=false
_has_plan_format=false
_has_plan_details=false
_has_strict_plan=false
_has_dry_run=false
_has_help=false

for _arg in "$@"; do
    case "${_arg}" in
        --version|-V)
            _has_version=true
            ;;
        --versions)
            _has_versions=true
            ;;
        --plan)
            _has_plan=true
            ;;
        --plan-format|--plan-format=*)
            _has_plan_format=true
            ;;
        --plan-details|--plan-details=*)
            _has_plan_details=true
            ;;
        --strict-plan)
            _has_strict_plan=true
            ;;
        --dry-run)
            _has_dry_run=true
            ;;
        --help|-h)
            _has_help=true
            ;;
    esac
done

if [[ "${_has_version}" == true && "${_has_versions}" == true ]]; then
    echo "[ERROR] Cannot combine --version/-V with --versions." >&2
    exit 2
fi

if [[ "${_has_plan}" == true && ( "${_has_version}" == true || "${_has_versions}" == true ) ]]; then
    echo "[ERROR] Cannot combine --plan with --version or --versions." >&2
    exit 2
fi

if [[ "${_has_plan}" == true && "${_has_dry_run}" == true ]]; then
    echo "[ERROR] Cannot combine --plan with --dry-run." >&2
    exit 2
fi

if [[ "${_has_plan}" == false && ( "${_has_plan_format}" == true || "${_has_plan_details}" == true || "${_has_strict_plan}" == true ) ]]; then
    echo "[ERROR] --plan-format, --plan-details, and --strict-plan require --plan." >&2
    exit 2
fi

if [[ "${_has_version}" == true ]]; then
    printf '%s %s (%s)\n' "${PIPELINE_NAME}" "${PIPELINE_VERSION}" "${PIPELINE_RELEASE_STATUS}"
    exit 0
fi

if [[ "${_has_plan}" == true && "${_has_help}" == true ]]; then
    _show_plan_help
fi

if [[ "${_has_versions}" == true && "${_has_help}" == true ]]; then
    _show_versions_help
fi

if [[ "${_has_help}" == true && "${_has_versions}" == false && "${_has_plan}" == false ]]; then
    _show_help
fi

if [[ "${_has_versions}" == true ]]; then
    _handle_versions_mode "$@"
fi

# The last profile is applied first; individual CLI options are then applied
# in original order and always take precedence.
_profile_seen=""
_profile_set() {
    local value="$1"
    case "${value}" in
        full|taxonomy|assembly|custom) ;;
        *) echo "[ERROR] Invalid profile '${value}' (expected full, taxonomy, assembly, or custom)." >&2; exit 2 ;;
    esac
    _profile_seen="${value}"
}
_profile_scan=("$@")
while [[ ${#_profile_scan[@]} -gt 0 ]]; do
    case "${_profile_scan[0]}" in
        --profile)
            [[ ${#_profile_scan[@]} -ge 2 && -n "${_profile_scan[1]}" ]] || {
                echo "[ERROR] --profile requires a value." >&2; exit 2; }
            _profile_set "${_profile_scan[1]}"
            _profile_scan=("${_profile_scan[@]:2}") ;;
        --taxonomy-only) _profile_set taxonomy; _profile_scan=("${_profile_scan[@]:1}") ;;
        --assembly-only) _profile_set assembly; _profile_scan=("${_profile_scan[@]:1}") ;;
        -r|--reads-dir|--clean-reads-dir|--kraken2-db|--checkm2-db|--gtdbtk-data-path|\
        -a|--assembly-in|\
        --bam-dir|--bam-manifest|--jgi-depth-file|--force-stage|--restart-from|\
        --binners|--refinement-mode|--bins-dir|--bins-source|--vamb-binsplit-separator|\
        --comebin-device|--comebin-seed|--min-completeness|--max-contamination|\
        --tmpdir|--plan-format|--plan-details|\
        -o|--outdir|--metabinner-scripts-dir|-t|--threads|-m|--memory|-e|--engine|--stop-after)
            _profile_scan=("${_profile_scan[@]:2}") ;;
        *) _profile_scan=("${_profile_scan[@]:1}") ;;
    esac
done
[[ -n "${_profile_seen}" ]] && REQUESTED_PROFILE="${_profile_seen}"
RESULTING_PROFILE="${REQUESTED_PROFILE}"
case "${REQUESTED_PROFILE}" in
    full|taxonomy|assembly|custom) ;;
    *) echo "[ERROR] Invalid profile '${REQUESTED_PROFILE}' (expected full, taxonomy, assembly, or custom)." >&2; exit 2 ;;
esac

case "${REQUESTED_PROFILE}" in
    full)
        RUN_MEGAHIT=true
        RUN_KRAKEN2=true
        RUN_MAPPING=true
        RUN_BINNING=true
        REFINEMENT_MODE=auto
        RUN_CHECKM2=true
        RUN_MAG_FILTERING=true
        RUN_COVERM=true
        RUN_SUMMARY=true
        RUN_TAXONOMY_SUMMARY=true
        RUN_MAG_SUMMARY=true
        RUN_RUN_OVERVIEW=true ;;
    taxonomy)
        RUN_KRAKEN2=true
        RUN_MEGAHIT=false
        RUN_MAPPING=false
        RUN_BINNING=false
        REFINEMENT_MODE=never
        RUN_REASSEMBLY=false
        RUN_CHECKM2=false
        RUN_MAG_FILTERING=false
        RUN_GTDBTK=false
        RUN_COVERM=false
        RUN_SUMMARY=true
        RUN_TAXONOMY_SUMMARY=true
        RUN_MAG_SUMMARY=false
        RUN_RUN_OVERVIEW=true ;;
    assembly)
        RUN_KRAKEN2=false
        RUN_MEGAHIT=true
        RUN_MAPPING=true
        RUN_BINNING=true
        REFINEMENT_MODE=auto
        RUN_CHECKM2=true
        RUN_MAG_FILTERING=true
        RUN_COVERM=true
        RUN_SUMMARY=true
        RUN_TAXONOMY_SUMMARY=false
        RUN_MAG_SUMMARY=true
        RUN_RUN_OVERVIEW=true ;;
    custom) ;;
esac

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile) shift 2 ;;
        --taxonomy-only|--assembly-only) shift ;;
        -r|--reads-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] $1 requires a directory value." >&2; exit 1; }
            READS_DIR="$2"; RAW_READS_DIR="$2"; shift 2 ;;
        --clean-reads-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --clean-reads-dir requires a directory value." >&2; exit 1; }
            CLEAN_READS_DIR="$2"; shift 2 ;;
        --kraken2-db)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --kraken2-db requires a directory value." >&2; exit 1; }
            KRAKEN2_DB="$2"; shift 2 ;;
        --checkm2-db)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --checkm2-db requires a file value." >&2; exit 1; }
            CHECKM2_DB="$2"; shift 2 ;;
        --gtdbtk-data-path)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --gtdbtk-data-path requires a directory value." >&2; exit 1; }
            GTDBTK_DATA_PATH="$2"; shift 2 ;;
        --binners)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --binners requires a comma-separated method list." >&2; exit 1; }
            SELECTED_BINNERS="$2"; shift 2 ;;
        --comebin-device)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --comebin-device requires exactly one value (cpu, cuda, or cuda:N)." >&2; exit 1; }
            COMEBIN_DEVICE="$2"
            COMEBIN_DEVICE_CLI_SET=true
            shift 2 ;;
        --comebin-seed)
            [[ -z "${2:-}" || ! "${2}" =~ ^[0-9]+$ ]] && {
                echo "[ERROR] --comebin-seed requires a non-negative integer." >&2; exit 1; }
            COMEBIN_SEED="$2"; shift 2 ;;
        --min-completeness)
            [[ "${2+x}" != x ]] && {
                echo "[ERROR] --min-completeness requires a value." >&2; exit 1; }
            MIN_COMPLETENESS="$2"
            shift 2 ;;
        --max-contamination)
            [[ "${2+x}" != x ]] && {
                echo "[ERROR] --max-contamination requires a value." >&2; exit 1; }
            MAX_CONTAMINATION="$2"
            shift 2 ;;
        --refinement-mode)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --refinement-mode requires auto, always, or never." >&2; exit 1; }
            REFINEMENT_MODE="$2"; shift 2 ;;
        --bins-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --bins-dir requires a directory value." >&2; exit 1; }
            BINS_DIR_INPUT="$2"; shift 2 ;;
        --bins-source)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --bins-source requires auto, vamb, metabinner, comebin, refined, or external." >&2; exit 1; }
            case "${2}" in
                auto|vamb|metabinner|comebin|refined|external) BINS_SOURCE="$2" ;;
                *)
                    echo "\`--bins-source\` selects a logical collection. Use \`--bins-dir DIR\` to provide the location of external bins." >&2
                    exit 1
                    ;;
            esac
            shift 2 ;;
        --vamb-binsplit-separator)
            [[ "${2+x}" != x ]] && {
                echo "[ERROR] --vamb-binsplit-separator requires exactly one non-empty character." >&2; exit 1; }
            VAMB_BINSPLIT_SEPARATOR="$2"
            VAMB_BINSPLIT_SEPARATOR_WAS_SET=true
            shift 2 ;;
        -a|--assembly-in)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] $1 requires a FASTA file value." >&2; exit 1; }
            ASSEMBLY_IN="$2"; shift 2 ;;
        --bam-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --bam-dir requires a directory value." >&2; exit 1; }
            BAM_DIR_INPUT="$2"; shift 2 ;;
        --bam-manifest)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --bam-manifest requires a file value." >&2; exit 1; }
            BAM_MANIFEST_INPUT="$2"; shift 2 ;;
        --jgi-depth-file)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --jgi-depth-file requires a file value." >&2; exit 1; }
            JGI_DEPTH_FILE_INPUT="$2"; shift 2 ;;
        -o|--outdir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] $1 requires a directory value." >&2; exit 1; }
            OUTDIR="$2"; shift 2 ;;
        --metabinner-scripts-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] --metabinner-scripts-dir requires a directory value." >&2; exit 1; }
            _METABINNER_SCRIPTS_DIR_CLI="$2"; shift 2 ;;
        -t|--threads)
            [[ -z "${2:-}" || ! "${2}" =~ ^[1-9][0-9]*$ ]] && {
                echo "[ERROR] --threads requires a positive integer." >&2; exit 1; }
            THREADS="$2"; shift 2 ;;
        -m|--memory)
            [[ -z "${2:-}" || ! "${2}" =~ ^[1-9][0-9]*$ ]] && {
                echo "[ERROR] --memory requires a positive integer (GB)." >&2; exit 1; }
            MEMORY="$2"; shift 2 ;;
        -e|--engine)
            [[ -z "${2:-}" ]] && { echo "[ERROR] --engine requires a value (apptainer|mamba)." >&2; exit 1; }
            case "${2}" in
                apptainer|mamba) EXEC_MODE="${2}"; shift 2 ;;
                *) echo "[ERROR] --engine must be 'apptainer' or 'mamba'." >&2; exit 1 ;;
            esac ;;
        --fastp)       RUN_FASTP=true;       shift ;;
        --no-fastp)    RUN_FASTP=false;      shift ;;
        --megahit)     RUN_MEGAHIT=true;     shift ;;
        --no-megahit)  RUN_MEGAHIT=false;    shift ;;
        --kraken2)     RUN_KRAKEN2=true;     shift ;;
        --no-kraken2)  RUN_KRAKEN2=false;    shift ;;
        --bracken)     RUN_BRACKEN=true;     shift ;;
        --no-bracken)  RUN_BRACKEN=false;    shift ;;
        --reassembly)     RUN_REASSEMBLY=true;  shift ;;
        --no-reassembly)  RUN_REASSEMBLY=false; shift ;;
        --tmpdir)
            # Reject absent or empty argument; also reject another recognized
            # CLI option appearing where the value should be, but allow paths
            # that begin with '-' (e.g. ./-temporary/work).
            if [[ -z "${2+set}" || -z "${2}" ]]; then
                echo "[ERROR] --tmpdir requires a directory value." >&2; exit 1
            fi
            # Check whether the next token is a recognized pipeline option.
            # We test a curated set of single-token flags and two-token keys;
            # any genuine path cannot match this list even if it starts with '-'.
            case "${2}" in
                --tmpdir|--keep-temporary|--no-keep-temporary|\
                --reassembly|--no-reassembly|--dry-run|--no-dry-run|\
                --force|--no-force|--resume|--force-stage|--restart-from|\
                --plan|--plan-format|--plan-details|--strict-plan|\
                --fastp|--no-fastp|--megahit|--no-megahit|--kraken2|--no-kraken2|\
                --bracken|--no-bracken|--gtdbtk|--no-gtdbtk|--mapping|--no-mapping|\
                --binning|--no-binning|--refinement|--no-refinement|--checkm2|--no-checkm2|\
                --mag-filtering|--no-mag-filtering|--coverm|--no-coverm|\
                --summary|--no-summary|--taxonomy-summary|--no-taxonomy-summary|\
                --mag-summary|--no-mag-summary|--run-overview|--no-run-overview|\
                --taxonomy-only|--assembly-only|\
                -r|--reads-dir|--clean-reads-dir|--kraken2-db|--checkm2-db|\
                --gtdbtk-data-path|-a|--assembly-in|--bam-dir|--bam-manifest|\
                --jgi-depth-file|--binners|\
                --refinement-mode|--bins-dir|--bins-source|--vamb-binsplit-separator|\
                --comebin-device|--comebin-seed|--min-completeness|--max-contamination|\
                -o|--outdir|--metabinner-scripts-dir|-t|--threads|-m|--memory|\
                -e|--engine|--stop-after|--profile|-h|--help)
                    echo "[ERROR] --tmpdir requires a directory value (got recognized option '${2}')." >&2
                    exit 1 ;;
            esac
            PIPELINE_TMPDIR_CLI="${2}"; shift 2 ;;
        --keep-temporary)    KEEP_TEMPORARY=true;  shift ;;
        --no-keep-temporary) KEEP_TEMPORARY=false; shift ;;
        --gtdbtk)     RUN_GTDBTK=true;    shift ;;
        --no-gtdbtk)  RUN_GTDBTK=false;   shift ;;
        --mapping) RUN_MAPPING=true; shift ;;
        --no-mapping) RUN_MAPPING=false; shift ;;
        --binning) RUN_BINNING=true; shift ;;
        --no-binning) RUN_BINNING=false; shift ;;
        --refinement) REFINEMENT_MODE=always; shift ;;
        --no-refinement) REFINEMENT_MODE=never; shift ;;
        --checkm2) RUN_CHECKM2=true; shift ;;
        --no-checkm2) RUN_CHECKM2=false; shift ;;
        --mag-filtering) RUN_MAG_FILTERING=true; shift ;;
        --no-mag-filtering) RUN_MAG_FILTERING=false; shift ;;
        --coverm) RUN_COVERM=true; shift ;;
        --no-coverm) RUN_COVERM=false; shift ;;
        --summary) RUN_SUMMARY=true; shift ;;
        --no-summary) RUN_SUMMARY=false; shift ;;
        --taxonomy-summary) RUN_TAXONOMY_SUMMARY=true; shift ;;
        --no-taxonomy-summary) RUN_TAXONOMY_SUMMARY=false; shift ;;
        --mag-summary) RUN_MAG_SUMMARY=true; shift ;;
        --no-mag-summary) RUN_MAG_SUMMARY=false; shift ;;
        --run-overview) RUN_RUN_OVERVIEW=true; shift ;;
        --no-run-overview) RUN_RUN_OVERVIEW=false; shift ;;
        --stop-after)
            [[ -z "${2:-}" || "${2:-}" == -* ]] && {
                echo "[ERROR] --stop-after requires a stage value." >&2; exit 2; }
            STOP_AFTER="$2"; shift 2 ;;
        --plan)
            PLAN_MODE=true; shift ;;
        --plan-format)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --plan-format requires a value (text or json)." >&2; exit 2; }
            case "${2}" in
                text|json) PLAN_FORMAT="${2}" ;;
                *) echo "[ERROR] Invalid --plan-format '${2}' (expected text or json)." >&2; exit 2 ;;
            esac
            shift 2 ;;
        --plan-details)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --plan-details requires a value (summary or full)." >&2; exit 2; }
            case "${2}" in
                summary|full) PLAN_DETAILS="${2}" ;;
                *) echo "[ERROR] Invalid --plan-details '${2}' (expected summary or full)." >&2; exit 2 ;;
            esac
            shift 2 ;;
        --strict-plan)
            STRICT_PLAN=true; shift ;;
        --dry-run)     DRY_RUN=true;      shift ;;
        --no-dry-run)  DRY_RUN=false;     shift ;;
        --force)     FORCE_RERUN=true; FORCE_REQUESTED=true; shift ;;
        --no-force)  FORCE_RERUN=false; FORCE_REQUESTED=false; shift ;;
        --resume) RESUME_REQUESTED=true; shift ;;
        --force-stage)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --force-stage requires a stage name." >&2; exit 2; }
            resolve_stage_id "${2}" >/dev/null 2>&1 || { echo "[ERROR] Unknown stage '${2}' for --force-stage." >&2; exit 2; }
            FORCE_STAGE_NAMES+=("$2"); shift 2 ;;
        --restart-from)
            [[ -z "${2:-}" || "${2}" == -* ]] && { echo "[ERROR] --restart-from requires a stage name." >&2; exit 2; }
            resolve_stage_id "${2}" >/dev/null 2>&1 || { echo "[ERROR] Unknown stage '${2}' for --restart-from." >&2; exit 2; }
            RESTART_REQUEST_COUNT=$((RESTART_REQUEST_COUNT + 1))
            (( RESTART_REQUEST_COUNT == 1 )) || {
                echo "[ERROR] --restart-from may be specified only once." >&2; exit 2; }
            RESTART_FROM_NAME="$2"; RESTART_REQUESTED=true; shift 2 ;;
        -h|--help)  _show_help ;;
        *)  echo "[ERROR] Unknown option: $1. Use -h for help." >&2; exit 2 ;;
    esac
done

if [[ "${VAMB_BINSPLIT_SEPARATOR_WAS_SET}" == true ]]; then
    python3 - "${VAMB_BINSPLIT_SEPARATOR}" <<'PY_VAMB_SEPARATOR_OPTION'
import sys
separator = sys.argv[1]
if len(separator) != 1 or separator.isspace() or separator in "/>":
    raise SystemExit(
        "[ERROR] --vamb-binsplit-separator requires exactly one non-whitespace "
        "character other than '/' or '>'."
    )
PY_VAMB_SEPARATOR_OPTION
fi

if [[ "${PLAN_MODE}" != true ]]; then
    if [[ "${_has_plan_format}" == true || "${_has_plan_details}" == true || "${STRICT_PLAN}" == true ]]; then
        echo "[ERROR] --plan-format, --plan-details, and --strict-plan require --plan." >&2
        exit 2
    fi
fi

if [[ "${PLAN_MODE}" == true && "${DRY_RUN}" == true ]]; then
    echo "[ERROR] --plan cannot be combined with --dry-run." >&2
    exit 2
fi

if [[ "${FORCE_RERUN}" == true ]]; then FORCE_REQUESTED=true; fi
if [[ "${RESUME_REQUESTED}" == true && ( "${FORCE_REQUESTED}" == true || \
      ${#FORCE_STAGE_NAMES[@]} -gt 0 || "${RESTART_REQUESTED}" == true ) ]]; then
    echo "[ERROR] --resume cannot be combined with --force, --force-stage, or --restart-from." >&2
    exit 2
fi
if [[ "${FORCE_REQUESTED}" == true && ( ${#FORCE_STAGE_NAMES[@]} -gt 0 || "${RESTART_REQUESTED}" == true ) ]]; then
    echo "[ERROR] --force cannot be combined with --force-stage or --restart-from." >&2
    exit 2
fi
if (( ${#FORCE_STAGE_NAMES[@]} > 0 )) && [[ "${RESTART_REQUESTED}" == true ]]; then
    echo "[ERROR] --force-stage cannot be combined with --restart-from." >&2
    exit 2
fi
if [[ "${RESTART_REQUESTED}" == true ]]; then
    RESUME_CONTROL=restart
elif (( ${#FORCE_STAGE_NAMES[@]} > 0 )); then
    RESUME_CONTROL=force-stage
elif [[ "${FORCE_REQUESTED}" == true ]]; then
    RESUME_CONTROL=force
else
    RESUME_CONTROL=resume
fi

# Resolve MetaBinner's path only after parsing, since --engine may override the environment.
if [[ -n "${_METABINNER_SCRIPTS_DIR_CLI:-}" ]]; then
    METABINNER_SCRIPTS_DIR="${_METABINNER_SCRIPTS_DIR_CLI}"
elif [[ -n "${_METABINNER_SCRIPTS_DIR_ENV}" ]]; then
    METABINNER_SCRIPTS_DIR="${_METABINNER_SCRIPTS_DIR_ENV}"
elif [[ "${EXEC_MODE}" == "apptainer" ]]; then
    METABINNER_SCRIPTS_DIR="/opt/metabinner/scripts"
fi

[[ -n "${BIN_DIR}" ]] && export PATH="${BIN_DIR}:${PATH}"

# ── Colour Helpers & Logging ──────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

visible_stage_number() {
    local requested="$1" internal id enabled number=0
    case "${requested}" in
        01) internal=01_fastp ;;
        02) internal=02_megahit ;;
        03) internal=03_kraken2 ;;
        04) internal=04_mapping ;;
        05a) internal=05a_vamb ;;
        05b) internal=05b_metabinner ;;
        05c) internal=05c_comebin ;;
        06) internal=06_refinement ;;
        07) internal=07_reassembly ;;
        08) internal=08_checkm2 ;;
        09) internal=09_filter_mags ;;
        10) internal=10_gtdbtk ;;
        11) internal=11_coverm ;;
        12a) internal=12a_taxonomy_summary ;;
        12b) internal=12b_mag_summary ;;
        12c) internal=12c_run_overview ;;
        12) internal=12c_run_overview ;;
        *) printf '%s' "${requested}"; return ;;
    esac
    local -a order=(01_fastp 03_kraken2 02_megahit 04_mapping 05a_vamb 05b_metabinner
                   05c_comebin 06_refinement 07_reassembly 08_checkm2 09_filter_mags
                   10_gtdbtk 11_coverm 12a_taxonomy_summary 12b_mag_summary 12c_run_overview)
    for id in "${order[@]}"; do
        case "${id}" in
            01_fastp) enabled="${RUN_FASTP}" ;;
            03_kraken2) enabled="${RUN_KRAKEN2}" ;;
            02_megahit) enabled="${RUN_MEGAHIT}" ;;
            04_mapping) enabled="${RUN_MAPPING}" ;;
            05a_vamb) enabled=false; [[ "${RUN_BINNING}" == true ]] && binner_selected vamb && enabled=true ;;
            05b_metabinner) enabled=false; [[ "${RUN_BINNING}" == true ]] && binner_selected metabinner && enabled=true ;;
            05c_comebin) enabled=false; [[ "${RUN_BINNING}" == true ]] && binner_selected comebin && enabled=true ;;
            06_refinement) enabled="${RUN_REFINEMENT_EFFECTIVE}" ;;
            07_reassembly) enabled="${RUN_REASSEMBLY}" ;;
            08_checkm2) enabled="${RUN_CHECKM2}" ;;
            09_filter_mags) enabled="${RUN_MAG_FILTERING}" ;;
            10_gtdbtk) enabled="${RUN_GTDBTK}" ;;
            11_coverm) enabled="${RUN_COVERM}" ;;
            12a_taxonomy_summary) enabled=false; [[ "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true && "${RUN_KRAKEN2}" == true ]] && enabled=true ;;
            12b_mag_summary) enabled=false; [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true && ( "${RUN_CHECKM2}" == true || "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || "${RUN_COVERM}" == true ) ]] && enabled=true ;;
            12c_run_overview) enabled=false; [[ "${RUN_SUMMARY}" == true && "${RUN_RUN_OVERVIEW}" == true ]] && enabled=true ;;
        esac
        [[ "${enabled}" == true ]] || continue
        number=$((number + 1))
        [[ "${id}" == "${internal}" ]] && { printf '%02d' "${number}"; return; }
    done
    printf '%s' "${requested}"
}

_format_visible_steps() {
    local text="$*" internal visible
    local -a replacements=()
    for internal in 01 02 03 04 05a 05b 05c 06 07 08 09 10 11 12 12a 12b 12c; do
        replacements+=(-e "s/\\[Step ${internal}\\]/__STEP_${internal}__/g")
    done
    for internal in 01 02 03 04 05a 05b 05c 06 07 08 09 10 11 12 12a 12b 12c; do
        visible="$(visible_stage_number "${internal}")"
        replacements+=(-e "s/__STEP_${internal}__/[Step ${visible}]/g")
    done
    text="$(printf '%s' "${text}" | sed "${replacements[@]}")"
    printf '%s' "${text}"
}

log_info()  { echo -e "${GREEN}[INFO]${NC}  $(date '+%Y-%m-%d %H:%M:%S')  $(_format_visible_steps "$*")"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $(date '+%Y-%m-%d %H:%M:%S')  $(_format_visible_steps "$*")"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S')  $(_format_visible_steps "$*")" >&2; }
log_step()  {
    local n="$1"; shift
    n="$(visible_stage_number "${n}")"
    echo -e "\n${CYAN}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}${BOLD}║  [Step ${n}]  $*${NC}"
    echo -e "${CYAN}${BOLD}╚══════════════════════════════════════════════════════════════╝${NC}\n"
}

validate_percentage_threshold() {
    local name="$1"
    local value="$2"
    if ! LC_ALL=C awk -v value="${value}" '
        BEGIN {
            if (value !~ /^[0-9]+([.][0-9]+)?$/ || value < 0 || value > 100) {
                exit 1
            }
        }
    '; then
        log_error "${name} must be a number between 0 and 100: ${value}"
        return 1
    fi
}

# ── Sanitization Helpers ──────────────────────────────────────────────────────
sanitize_command() {
    local -a cmd=("$@")
    local str=""
    if (( ${#cmd[@]} > 0 )); then
        printf -v str '%q ' "${cmd[@]}"
        str="${str% }"
    fi
    str="$(printf '%s' "${str}" | sed -E \
        -e 's/(--(token|password|secret|credential|api-key|key|access-token))([= ])([^ ]+)/\1\3<redacted>/gI' \
        -e 's/\b(TOKEN|PASSWORD|SECRET|CREDENTIAL|API_KEY|ACCESS_TOKEN)=([^ ]+)/\1=<redacted>/gI')"
    printf '%s' "${str}"
}

sanitize_shell_command() {
    local str="${1:-}"
    str="$(printf '%s' "${str}" | sed -E \
        -e 's/(--(token|password|secret|credential|api-key|key|access-token))([= ])([^ ]+)/\1\3<redacted>/gI' \
        -e 's/\b(TOKEN|PASSWORD|SECRET|CREDENTIAL|API_KEY|ACCESS_TOKEN)=([^ ]+)/\1=<redacted>/gI')"
    printf '%s' "${str}"
}

sanitize_tsv_field() {
    local val="${1:-unknown}"
    local default="${2:-unknown}"
    val="$(printf '%s' "${val}" | tr '\r\n\t' '   ' | sed -E 's/[[:space:]]+/ /g' | sed -E 's/^ +| +$//g')"
    val="$(sanitize_shell_command "${val}")"
    [[ -z "${val}" ]] && val="${default}"
    printf '%s' "${val}"
}

# ── Temporary Workspace Management ───────────────────────────────────────────
#
# Tracking variables for error records and status reporting
TEMPORARY_CLEANUP_STATUS="not-applicable"
TEMPORARY_FILES_PRESERVED=false

#
# resolve_pipeline_tmpdir_base
#   Applies precedence: CLI > TMPDIR env > fallback (OUTDIR/.tmp).
#   Sets PIPELINE_TMPDIR_SOURCE and PIPELINE_TMPDIR_BASE.
#   Must be called during preflight after OUTDIR is known.
#   Does NOT create any directories — that happens in create_pipeline_run_tmpdir.
#
resolve_pipeline_tmpdir_base() {
    local candidate="" source="" outdir_abs=""
    outdir_abs="$(realpath -m -- "${OUTDIR}" 2>/dev/null || printf '%s' "${OUTDIR}")"

    # ── Canonicalize a path (may or may not exist) ────────────────────────────
    _canonicalize_tmpdir_path() {
        local raw="$1"
        local abs=""

        if [[ -L "${raw}" ]]; then
            log_error "[Tmpdir] Base directory is a symlink — refusing: ${raw}"
            return 1
        fi

        if [[ -e "${raw}" ]]; then
            abs="$(realpath -- "${raw}" 2>/dev/null)" || {
                log_error "[Tmpdir] realpath failed for existing path: ${raw}"
                return 1
            }
        else
            # realpath -m accepts non-existent paths; fall back to manual
            if realpath --help 2>&1 | grep -q '\-m'; then
                abs="$(realpath -m -- "${raw}" 2>/dev/null)" || {
                    log_error "[Tmpdir] realpath -m failed for: ${raw}"
                    return 1
                }
            else
                # Manual: resolve the deepest existing ancestor, then append
                local dir basename rest
                dir="${raw}"
                rest=""
                while [[ -n "${dir}" && "${dir}" != "/" && ! -e "${dir}" ]]; do
                    basename="${dir##*/}"
                    dir="${dir%/*}"
                    rest="${basename}${rest:+/}${rest}"
                done
                [[ -z "${dir}" ]] && dir="/"
                local dir_abs
                dir_abs="$(realpath -- "${dir}" 2>/dev/null || printf '%s' "${dir}")"
                abs="${dir_abs}${rest:+/}${rest}"
            fi
        fi
        # Must be absolute
        if [[ "${abs:0:1}" != "/" ]]; then
            log_error "[Tmpdir] Resolved path is not absolute: ${abs} (from: ${raw})"
            return 1
        fi
        # Reject ambiguous ..
        if [[ "${abs}" == *"/../"* || "${abs}" == *"/.." ]]; then
            log_error "[Tmpdir] Resolved path still contains '..': ${abs}"
            return 1
        fi
        printf '%s' "${abs}"
    }

    # ── Location restriction checks ───────────────────────────────────────────
    _check_tmpdir_location() {
        local path="$1"
        local src="$2"

        if [[ -z "${path}" ]]; then
            log_error "[Tmpdir] Resolved tmpdir path is empty."
            return 1
        fi
        if [[ "${path}" == "/" ]]; then
            log_error "[Tmpdir] Refusing to use filesystem root as tmpdir: ${path}"
            return 1
        fi
        if [[ -L "${path}" ]]; then
            log_error "[Tmpdir] Refusing symlink as tmpdir base: ${path}"
            return 1
        fi
        local home_abs=""
        home_abs="$(realpath -- "${HOME:-/nonexistent}" 2>/dev/null || printf '%s' "${HOME:-}")"
        if [[ -n "${home_abs}" && "${path}" == "${home_abs}" ]]; then
            log_error "[Tmpdir] Refusing to use HOME as tmpdir: ${path}"
            return 1
        fi
        if [[ "${path}" == "${outdir_abs}" ]]; then
            log_error "[Tmpdir] Refusing to use OUTDIR itself as tmpdir: ${path}"
            return 1
        fi
        # Reject paths inside OUTDIR — unless this is the fallback tier
        if [[ "${src}" != "fallback" && ( "${path}" == "${outdir_abs}/"* || "${path}" == "${outdir_abs}" ) ]]; then
            log_error "[Tmpdir] Refusing temporary directory inside OUTDIR (source=${src}): ${path}"
            log_error "         Use --tmpdir with a path outside ${outdir_abs}."
            return 1
        fi
        # Even for fallback, only ${OUTDIR}/.tmp is allowed
        if [[ "${src}" == "fallback" && "${path}" != "${outdir_abs}/.tmp" ]]; then
            log_error "[Tmpdir] Internal: fallback path mismatch: expected ${outdir_abs}/.tmp, got ${path}"
            return 1
        fi
        # Repo root guard
        local repo_root=""
        repo_root="$(realpath -m -- "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null || true)"
        if [[ -n "${repo_root}" && "${path}" == "${repo_root}" ]]; then
            log_error "[Tmpdir] Refusing to use the repository root as tmpdir: ${path}"
            return 1
        fi
        return 0
    }

    # ── Tier: CLI ─────────────────────────────────────────────────────────────
    if [[ -n "${PIPELINE_TMPDIR_CLI}" ]]; then
        source="cli"
        candidate="$(_canonicalize_tmpdir_path "${PIPELINE_TMPDIR_CLI}")" || return 1
        _check_tmpdir_location "${candidate}" "${source}" || return 1
        PIPELINE_TMPDIR_SOURCE="${source}"
        PIPELINE_TMPDIR_BASE="${candidate}"
        PIPELINE_TMP_PARENT="${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline"
        log_info "[Tmpdir] Source: cli | Base: ${PIPELINE_TMPDIR_BASE}"
        return 0
    fi

    # ── Tier: environment ─────────────────────────────────────────────────────
    if [[ -n "${TMPDIR_ENV_ORIGINAL}" ]]; then
        source="environment"
        candidate="$(_canonicalize_tmpdir_path "${TMPDIR_ENV_ORIGINAL}")" || return 1
        _check_tmpdir_location "${candidate}" "${source}" || return 1
        PIPELINE_TMPDIR_SOURCE="${source}"
        PIPELINE_TMPDIR_BASE="${candidate}"
        PIPELINE_TMP_PARENT="${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline"
        log_info "[Tmpdir] Source: environment | Base: ${PIPELINE_TMPDIR_BASE}"
        return 0
    fi

    # ── Tier: fallback ────────────────────────────────────────────────────────
    source="fallback"
    candidate="${outdir_abs}/.tmp"
    PIPELINE_TMPDIR_SOURCE="${source}"
    PIPELINE_TMPDIR_BASE="${candidate}"
    PIPELINE_TMP_PARENT="${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline"
    log_info "[Tmpdir] Source: fallback | Base: ${PIPELINE_TMPDIR_BASE}"
    return 0
}

# validate_pipeline_tmpdir_base
#   Validates an existing base directory: must be a non-symlink, writable dir
#   whose canonical path matches what was resolved.
#   If the base does not yet exist, validates the nearest existing ancestor.
validate_pipeline_tmpdir_base() {
    local base="${PIPELINE_TMPDIR_BASE}"
    if [[ -z "${base}" ]]; then
        log_error "[Tmpdir] PIPELINE_TMPDIR_BASE is empty."
        return 1
    fi
    if [[ -e "${base}" ]]; then
        if [[ -L "${base}" ]]; then
            log_error "[Tmpdir] Base is a symlink — refusing: ${base}"
            return 1
        fi
        if [[ ! -d "${base}" ]]; then
            log_error "[Tmpdir] Base exists but is not a directory: ${base}"
            return 1
        fi
        if [[ ! -w "${base}" ]]; then
            log_error "[Tmpdir] Base directory is not writable: ${base}"
            return 1
        fi
        local real_base=""
        real_base="$(realpath -- "${base}" 2>/dev/null)" || {
            log_error "[Tmpdir] Cannot resolve canonical path for: ${base}"
            return 1
        }
        if [[ "${real_base}" != "${base}" ]]; then
            log_error "[Tmpdir] Base resolves to a different canonical path."
            log_error "         Expected: ${base}"
            log_error "         Resolved: ${real_base}"
            return 1
        fi
    else
        # Walk up to the nearest existing ancestor
        local ancestor="${base}"
        while [[ -n "${ancestor}" && "${ancestor}" != "/" && ! -e "${ancestor}" ]]; do
            ancestor="${ancestor%/*}"
        done
        [[ -z "${ancestor}" ]] && ancestor="/"
        if [[ -L "${ancestor}" ]]; then
            log_error "[Tmpdir] Nearest existing ancestor is a symlink: ${ancestor}"
            return 1
        fi
        if [[ ! -d "${ancestor}" ]]; then
            log_error "[Tmpdir] Nearest existing ancestor is not a directory: ${ancestor}"
            return 1
        fi
        if [[ ! -w "${ancestor}" ]]; then
            log_error "[Tmpdir] Nearest existing ancestor is not writable: ${ancestor}"
            return 1
        fi
    fi
    return 0
}

# ── Structured Disk Space Assessment ──────────────────────────────────────────

# find_existing_parent PATH
#   Finds the nearest existing ancestor directory for PATH.
#   If PATH exists and is a directory, returns PATH.
#   If PATH is a regular file, returns its dirname.
#   Walks up until an existing directory is found, defaulting to "/".
find_existing_parent() {
    local target="$1"
    [[ -z "${target}" ]] && target="/"
    local dir="${target}"
    if [[ -f "${dir}" ]]; then
        dir="${dir%/*}"
        [[ -z "${dir}" ]] && dir="/"
    fi
    while [[ -n "${dir}" && "${dir}" != "/" && ! -d "${dir}" ]]; do
        dir="${dir%/*}"
    done
    [[ -z "${dir}" ]] && dir="/"
    realpath -m -- "${dir}" 2>/dev/null || printf '%s' "${dir}"
}

# canonicalize_storage_path PATH
#   Resolves an absolute, canonical path for PATH without trailing slashes.
canonicalize_storage_path() {
    local raw="$1"
    [[ -z "${raw}" ]] && { printf '/'; return 0; }
    local abs=""
    if [[ -e "${raw}" ]]; then
        abs="$(realpath -- "${raw}" 2>/dev/null)" || abs="$(realpath -m -- "${raw}" 2>/dev/null)" || abs="${raw}"
    else
        abs="$(realpath -m -- "${raw}" 2>/dev/null)" || {
            local parent rest
            parent="$(find_existing_parent "${raw}")"
            rest="${raw#"${parent}"}"
            rest="${rest#/}"
            abs="${parent%/}/${rest}"
        }
    fi
    while [[ "${abs}" != "/" && "${abs}" == */ ]]; do
        abs="${abs%/}"
    done
    printf '%s' "${abs}"
}

# get_filesystem_device_id PATH
#   Queries the numeric device ID using stat -c '%d' on nearest existing ancestor.
get_filesystem_device_id() {
    local path="$1"
    local check_path
    check_path="$(find_existing_parent "${path}")"
    local dev_id
    dev_id="$(stat -c '%d' -- "${check_path}" 2>/dev/null || echo "")"
    if [[ -n "${dev_id}" && "${dev_id}" =~ ^[0-9]+$ ]]; then
        printf '%s' "${dev_id}"
    else
        printf 'unknown'
    fi
}

# get_filesystem_mount_point PATH
#   Queries the mount point using POSIX df -Pk on nearest existing ancestor.
get_filesystem_mount_point() {
    local path="$1"
    local check_path
    check_path="$(find_existing_parent "${path}")"
    local mount_pt
    mount_pt="$(LC_ALL=C df -Pk -- "${check_path}" 2>/dev/null | awk 'NR>1 { if (NF>=6) { print $6; exit } }')"
    if [[ -n "${mount_pt}" ]]; then
        printf '%s' "${mount_pt}"
    else
        printf 'unknown'
    fi
}

# get_filesystem_id PATH
#   Returns <device_id>|<mount_point> or unknown
get_filesystem_id() {
    local path="$1"
    local dev mnt
    dev="$(get_filesystem_device_id "${path}")"
    mnt="$(get_filesystem_mount_point "${path}")"
    if [[ "${dev}" == "unknown" && "${mnt}" == "unknown" ]]; then
        printf 'unknown'
    else
        printf '%s|%s' "${dev}" "${mnt}"
    fi
}

# compare_filesystem_identities FS_ID_1 FS_ID_2
#   Implements Section 9 logic to detect shared capacity between filesystems.
compare_filesystem_identities() {
    local fs1="$1" fs2="$2"
    if [[ "${fs1}" == "unknown" || "${fs2}" == "unknown" ]]; then
        printf 'unknown'
        return 0
    fi
    local dev1="${fs1%%|*}" mnt1="${fs1#*|}"
    local dev2="${fs2%%|*}" mnt2="${fs2#*|}"

    # 1. if device_id and mount_point match: shared=true
    if [[ "${dev1}" != "unknown" && "${dev2}" != "unknown" && "${dev1}" == "${dev2}" && \
          "${mnt1}" != "unknown" && "${mnt2}" != "unknown" && "${mnt1}" == "${mnt2}" ]]; then
        printf 'true'
        return 0
    fi

    # 2. if device_id missing, but mount_point matches: shared=true
    if [[ ( "${dev1}" == "unknown" || "${dev2}" == "unknown" ) && \
          "${mnt1}" != "unknown" && "${mnt2}" != "unknown" && "${mnt1}" == "${mnt2}" ]]; then
        printf 'true'
        return 0
    fi

    # 3. if device_id matches, but mount_point differs: shared=unknown
    if [[ "${dev1}" != "unknown" && "${dev2}" != "unknown" && "${dev1}" == "${dev2}" && \
          "${mnt1}" != "${mnt2}" ]]; then
        printf 'unknown'
        return 0
    fi

    # 4. if mount_point matches, but device_id differs: shared=unknown
    if [[ "${mnt1}" != "unknown" && "${mnt2}" != "unknown" && "${mnt1}" == "${mnt2}" && \
          "${dev1}" != "${dev2}" ]]; then
        printf 'unknown'
        return 0
    fi

    # 5. only declare false when consistent evidence of different filesystems
    if [[ "${dev1}" != "unknown" && "${dev2}" != "unknown" && "${dev1}" != "${dev2}" && \
          "${mnt1}" != "unknown" && "${mnt2}" != "unknown" && "${mnt1}" != "${mnt2}" ]]; then
        printf 'false'
        return 0
    fi

    # 6. on conflict / edge case: shared=unknown
    printf 'unknown'
}

# get_filesystem_capacity PATH
#   Queries capacity via LC_ALL=C df -Pk -- PATH.
#   Returns tab-separated: total_bytes used_bytes available_bytes mount_point
get_filesystem_capacity() {
    local path="$1"
    local check_path
    check_path="$(find_existing_parent "${path}")"
    local df_out
    df_out="$(LC_ALL=C df -Pk -- "${check_path}" 2>/dev/null | awk 'NR>1 { if (NF>=6 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/) { print $2, $3, $4, $6; exit } }')"
    if [[ -z "${df_out}" ]]; then
        printf 'NA\tNA\tNA\tunknown\n'
        return 0
    fi
    local tot_k used_k avail_k mnt
    tot_k="$(awk '{print $1}' <<< "${df_out}")"
    used_k="$(awk '{print $2}' <<< "${df_out}")"
    avail_k="$(awk '{print $3}' <<< "${df_out}")"
    mnt="$(awk '{print $4}' <<< "${df_out}")"

    local tot_bytes used_bytes avail_bytes
    tot_bytes="$(awk -v k="${tot_k}" 'BEGIN { printf "%.0f", k * 1024 }')"
    used_bytes="$(awk -v k="${used_k}" 'BEGIN { printf "%.0f", k * 1024 }')"
    avail_bytes="$(awk -v k="${avail_k}" 'BEGIN { printf "%.0f", k * 1024 }')"

    printf '%s\t%s\t%s\t%s\n' "${tot_bytes}" "${used_bytes}" "${avail_bytes}" "${mnt}"
}

# get_available_inodes PATH
#   Queries free inodes via LC_ALL=C df -Pi -- PATH.
#   Returns non-negative integer or NA.
get_available_inodes() {
    local path="$1"
    local check_path
    check_path="$(find_existing_parent "${path}")"
    local ifree
    ifree="$(LC_ALL=C df -Pi -- "${check_path}" 2>/dev/null | awk 'NR>1 { if (NF>=5) { print $4; exit } }')"
    if [[ -n "${ifree}" && "${ifree}" =~ ^[0-9]+$ ]]; then
        printf '%s' "${ifree}"
    else
        printf 'NA'
    fi
}

# sum_file_sizes FILE...
#   Sums file sizes in bytes.
sum_file_sizes() {
    local total=0 s f
    for f in "$@"; do
        if [[ -f "${f}" ]]; then
            s="$(stat -c '%s' -- "${f}" 2>/dev/null || echo 0)"
            if [[ "${s}" =~ ^[0-9]+$ ]]; then
                total=$(( total + s ))
            fi
        fi
    done
    printf '%s' "${total}"
}

# get_canonical_reads_bytes
#   Returns the total size in bytes of active reads from manifest or active arrays.
get_canonical_reads_bytes() {
    local manifest="${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}"
    if [[ -s "${manifest}" && -f "${manifest}" ]]; then
        local sum
        sum="$(awk -F'\t' 'NR>1 { if ($8 ~ /^[0-9]+$/) s += $8; if ($9 ~ /^[0-9]+$/) s += $9 } END { printf "%.0f", s }' "${manifest}" 2>/dev/null)"
        if [[ -n "${sum}" && "${sum}" =~ ^[0-9]+$ && "${sum}" -gt 0 ]]; then
            printf '%s' "${sum}"
            return 0
        fi
    fi
    local -a files=()
    if (( ${#ACTIVE_R1_FILES[@]} > 0 )); then
        files=("${ACTIVE_R1_FILES[@]}" "${ACTIVE_R2_FILES[@]}")
    elif (( ${#INPUT_R1_FILES[@]} > 0 )); then
        files=("${INPUT_R1_FILES[@]}" "${INPUT_R2_FILES[@]}")
    fi
    if (( ${#files[@]} > 0 )); then
        sum_file_sizes "${files[@]}"
        return 0
    fi
    printf '0'
}

# estimate_stage_disk_requirement STAGE_ID STORAGE_ROLE
#   Returns tab-separated:
#   estimated_required_bytes safety_factor estimate_confidence required_inodes blocking
estimate_stage_disk_requirement() {
    local stage="$1" role="$2"
    local reads_bytes
    reads_bytes="$(get_canonical_reads_bytes)"

    case "${stage}" in
        01_fastp)
            local factor="${FASTP_OUTPUT_FACTOR:-1.10}"
            local est
            est="$(awk -v b="${reads_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"
            printf '%s\t%s\t%s\t%s\t%s\n' "${est}" "${factor}" "informational" "NA" "false"
            ;;
        02_megahit)
            local factor="${MEGAHIT_WORK_FACTOR:-2.50}"
            local est
            est="$(awk -v b="${reads_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"
            printf '%s\t%s\t%s\t%s\t%s\n' "${est}" "${factor}" "informational" "NA" "false"
            ;;
        03_kraken2)
            local factor="${KRAKEN_OUTPUT_FACTOR:-1.00}"
            local est
            est="$(awk -v b="${reads_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"
            printf '%s\t%s\t%s\t%s\t%s\n' "${est}" "${factor}" "informational" "NA" "false"
            ;;
        04_mapping)
            local factor="${MAPPING_OUTPUT_FACTOR:-1.50}"
            local est
            est="$(awk -v b="${reads_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"
            printf '%s\t%s\t%s\t%s\t%s\n' "${est}" "${factor}" "informational" "NA" "false"
            ;;
        05a_vamb|05b_metabinner|05c_comebin)
            printf 'NA\tNA\tunknown\tNA\tfalse\n'
            ;;
        06_refinement)
            printf 'NA\tNA\tunknown\tNA\tfalse\n'
            ;;
        07_reassembly)
            if [[ "${role}" == "tmpdir" ]]; then
                local factor="${REASSEMBLY_CONCAT_SPACE_FACTOR:-1.10}"
                local est
                est="$(awk -v b="${reads_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"
                printf '%s\t%s\t%s\t%s\t%s\n' "${est}" "${factor}" "hard-minimum" "${REQUIRED_FREE_INODES:-16}" "true"
            else
                printf 'NA\tNA\tunknown\tNA\tfalse\n'
            fi
            ;;
        08_checkm2|10_gtdbtk|11_coverm)
            printf 'NA\tNA\tunknown\tNA\tfalse\n'
            ;;
        pipeline)
            printf 'NA\tNA\tunknown\tNA\tfalse\n'
            ;;
        *)
            printf 'NA\tNA\tnot-applicable\tNA\tfalse\n'
            ;;
    esac
}

# classify_disk_requirement AVAIL_BYTES REQ_BYTES CONFIDENCE BLOCKING AVAIL_INODES REQ_INODES
#   Determines sufficient | warning | insufficient | unknown | not-applicable.
classify_disk_requirement() {
    local avail_bytes="$1" req_bytes="$2" confidence="$3" blocking="$4" avail_inodes="$5" req_inodes="$6"
    local headroom_factor="${DISK_WARNING_HEADROOM_FACTOR:-1.25}"

    if [[ "${confidence}" == "not-applicable" ]]; then
        printf 'not-applicable'
        return 0
    fi
    if [[ "${req_bytes}" == "NA" ]]; then
        printf 'unknown'
        return 0
    fi

    local st="sufficient"

    if [[ "${avail_bytes}" == "NA" ]]; then
        st="warning"
    else
        local warn_threshold
        warn_threshold="$(awk -v r="${req_bytes}" -v h="${headroom_factor}" 'BEGIN { printf "%.0f", r * h }')"

        if [[ "${blocking}" == "true" ]]; then
            if (( avail_bytes < req_bytes )); then
                st="insufficient"
            elif (( avail_bytes < warn_threshold )); then
                st="warning"
            else
                st="sufficient"
            fi
        else
            # Informational never produces insufficient
            if (( avail_bytes < warn_threshold )); then
                st="warning"
            else
                st="sufficient"
            fi
        fi
    fi

    # Inode checks
    if [[ "${st}" != "insufficient" && "${req_inodes}" != "NA" && "${req_inodes}" =~ ^[0-9]+$ && "${req_inodes}" -gt 0 ]]; then
        if [[ "${avail_inodes}" == "NA" ]]; then
            st="warning"
        elif [[ "${avail_inodes}" =~ ^[0-9]+$ ]]; then
            if (( avail_inodes < req_inodes )); then
                if [[ "${blocking}" == "true" ]]; then
                    st="insufficient"
                else
                    st="warning"
                fi
            fi
        fi
    fi

    printf '%s' "${st}"
}

# record_disk_measurement ARGS...
#   Appends a measurement to disk_space.current.tsv (and disk_space.tsv if published).
record_disk_measurement() {
    if [[ "${DRY_RUN}" == true ]]; then
        return 0
    fi
    python3 "${WRITE_DISK_SPACE_REPORT}" record "$@" >/dev/null 2>&1
}

# publish_disk_space_report
#   Atomically publishes current.tsv to disk_space.tsv.
publish_disk_space_report() {
    if [[ "${DRY_RUN}" == true ]]; then
        return 0
    fi
    if [[ ! -s "${DISK_SPACE_CURRENT_TSV}" ]]; then
        return 0
    fi
    python3 "${WRITE_DISK_SPACE_REPORT}" publish \
        --report-current "${DISK_SPACE_CURRENT_TSV}" \
        --publish-path "${DISK_SPACE_TSV}" \
        --run-id "${RUN_ID}" >/dev/null 2>&1
}

# get_latest_disk_status
#   Returns summary of disk checks.
get_latest_disk_status() {
    if [[ ! -s "${DISK_SPACE_CURRENT_TSV}" ]]; then
        printf 'OVERALL_STATUS=not-applicable\nWARNING_COUNT=0\n'
        return 0
    fi
    python3 "${WRITE_DISK_SPACE_REPORT}" summary \
        --report-current "${DISK_SPACE_CURRENT_TSV}" \
        --run-id "${RUN_ID}" 2>/dev/null || printf 'OVERALL_STATUS=unknown\nWARNING_COUNT=0\n'
}

# check_disk_requirement STAGE_ID STORAGE_ROLE PURPOSE TARGET_PATH SHARED_GROUP MEASUREMENT_TYPE [OVERRIDE_REQ_BYTES] [OVERRIDE_IN_BYTES]
#   Executes full check: capacity inquiry, classification, logging, and recording.
check_disk_requirement() {
    local stage="$1" role="$2" purpose="$3" target_path="$4" group="$5" mtype="$6"
    local override_req="${7:-}" override_in="${8:-}"

    local check_id="${stage}.${role}.${purpose}"
    local existing_parent
    existing_parent="$(find_existing_parent "${target_path}")"
    local fs_id
    fs_id="$(get_filesystem_id "${existing_parent}")"

    local cap_info tot_bytes used_bytes avail_bytes
    cap_info="$(get_filesystem_capacity "${existing_parent}")"
    tot_bytes="$(awk '{print $1}' <<< "${cap_info}")"
    used_bytes="$(awk '{print $2}' <<< "${cap_info}")"
    avail_bytes="$(awk '{print $3}' <<< "${cap_info}")"

    local avail_inodes
    avail_inodes="$(get_available_inodes "${existing_parent}")"

    # Estimation
    local req_bytes safety_factor confidence req_inodes blocking
    if [[ -n "${override_req}" ]]; then
        req_bytes="${override_req}"
        safety_factor="${REASSEMBLY_CONCAT_SPACE_FACTOR:-1.10}"
        confidence="hard-minimum"
        req_inodes="${REQUIRED_FREE_INODES:-16}"
        blocking="true"
    else
        local est_line
        est_line="$(estimate_stage_disk_requirement "${stage}" "${role}")"
        req_bytes="$(awk '{print $1}' <<< "${est_line}")"
        safety_factor="$(awk '{print $2}' <<< "${est_line}")"
        confidence="$(awk '{print $3}' <<< "${est_line}")"
        req_inodes="$(awk '{print $4}' <<< "${est_line}")"
        blocking="$(awk '{print $5}' <<< "${est_line}")"
    fi

    local in_bytes="NA"
    if [[ -n "${override_in}" ]]; then
        in_bytes="${override_in}"
    elif [[ "${stage}" != "pipeline" && "${confidence}" != "not-applicable" ]]; then
        in_bytes="$(get_canonical_reads_bytes)"
    fi

    # Classify
    local status
    status="$(classify_disk_requirement "${avail_bytes}" "${req_bytes}" "${confidence}" "${blocking}" "${avail_inodes}" "${req_inodes}")"

    CURRENT_STAGE_DISK_STATUS="${status}"
    CURRENT_STAGE_DISK_REQ="${req_bytes}"
    CURRENT_STAGE_DISK_AVAIL="${avail_bytes}"

    # Build message
    local msg=""
    if [[ "${existing_parent}" != "${target_path}" ]]; then
        msg="filesystem_source=existing-parent;target_path=${target_path}"
    else
        msg="filesystem_source=real"
    fi
    if [[ "${confidence}" == "informational" && "${status}" == "warning" ]]; then
        msg="${msg};headroom_warning=available_below_headroom"
    fi

    local ts_utc
    ts_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "1970-01-01T00:00:00Z")"

    # Record measurement
    if [[ "${DRY_RUN}" != true ]]; then
        record_disk_measurement \
            --report-current "${DISK_SPACE_CURRENT_TSV}" \
            --publish-path "${DISK_SPACE_TSV}" \
            --run-id "${RUN_ID}" \
            --measurement-sequence "auto" \
            --check-id "${check_id}" \
            --timestamp-utc "${ts_utc}" \
            --stage-id "${stage}" \
            --storage-role "${role}" \
            --purpose "${purpose}" \
            --path "${target_path}" \
            --filesystem-id "${fs_id}" \
            --measurement-type "${mtype}" \
            --estimate-confidence "${confidence}" \
            --input-bytes-read "${in_bytes}" \
            --estimated-required-bytes "${req_bytes}" \
            --filesystem-available-bytes "${avail_bytes}" \
            --filesystem-used-bytes "${used_bytes}" \
            --filesystem-total-bytes "${tot_bytes}" \
            --required-free-inodes "${req_inodes}" \
            --available-free-inodes "${avail_inodes}" \
            --safety-factor "${safety_factor}" \
            --shared-capacity-group "${group}" \
            --status "${status}" \
            --message "${msg}"
    fi

    if [[ "${status}" == "insufficient" && "${blocking}" == "true" ]]; then
        log_error "Insufficient disk space for ${check_id}. Stage: ${stage}; filesystem: ${fs_id}; required: ${req_bytes} bytes; available: ${avail_bytes} bytes; path: ${target_path}; confidence: ${confidence}."
        return 1
    elif [[ "${status}" == "warning" ]]; then
        log_warn "Disk space warning for ${check_id}: available (${avail_bytes} B) is near or below estimated requirement (${req_bytes} B) on ${fs_id}."
    fi

    return 0
}

# check_tmpdir_space SOURCE_R1_ARRAY_NAME SOURCE_R2_ARRAY_NAME
#   Canonical implementation of reassembly temporary workspace disk check.
check_tmpdir_space() {
    local -n _r1_arr="$1"
    local -n _r2_arr="$2"
    local total_bytes
    total_bytes="$(sum_file_sizes "${_r1_arr[@]}" "${_r2_arr[@]}")"
    local factor="${REASSEMBLY_CONCAT_SPACE_FACTOR:-1.10}"
    local required_bytes
    required_bytes="$(awk -v b="${total_bytes}" -v f="${factor}" 'BEGIN { val = b * f; ival = int(val); if (val > ival) ival++; printf "%.0f", ival }')"

    local check_path="${PIPELINE_TMPDIR_BASE}"
    if ! check_disk_requirement "07_reassembly" "tmpdir" "read-concatenation" "${check_path}" "07_reassembly" "stage-precheck" "${required_bytes}" "${total_bytes}"; then
        return 1
    fi
    log_info "[Tmpdir] Space check for reassembly temporary concatenation passed (required: ${required_bytes} bytes, factor: ${factor}x)."
    return 0
}

# run_stage_precheck STAGE_ID
#   Remeasures filesystem capacity and executes precheck immediately before a stage.
run_stage_precheck() {
    local stage="$1"
    local r_group="none"
    [[ "${SHARED_CAPACITY_FS}" == "true" ]] && r_group="${stage}"

    case "${stage}" in
        01_fastp)
            check_disk_requirement "01_fastp" "outdir" "cleaned-reads" "${QC_DIR}" "none" "stage-precheck" || return 1
            ;;
        02_megahit)
            check_disk_requirement "02_megahit" "outdir" "assembly" "${ASSEMBLY_DIR}" "none" "stage-precheck" || return 1
            ;;
        03_kraken2)
            check_disk_requirement "03_kraken2" "outdir" "taxonomy-output" "${KRAKEN2_DIR}" "none" "stage-precheck" || return 1
            ;;
        04_mapping)
            check_disk_requirement "04_mapping" "outdir" "mapping-bam" "${MAPPING_DIR}" "none" "stage-precheck" || return 1
            ;;
        05a_vamb)
            check_disk_requirement "05a_vamb" "outdir" "binning" "${VAMB_DIR}" "none" "stage-precheck" || return 1
            ;;
        05b_metabinner)
            check_disk_requirement "05b_metabinner" "outdir" "binning" "${METABINNER_DIR}" "none" "stage-precheck" || return 1
            ;;
        05c_comebin)
            check_disk_requirement "05c_comebin" "outdir" "binning" "${COMEBIN_DIR}" "none" "stage-precheck" || return 1
            ;;
        06_refinement)
            check_disk_requirement "06_refinement" "outdir" "refinement" "${REFINEMENT_DIR}" "none" "stage-precheck" || return 1
            ;;
        07_reassembly)
            check_disk_requirement "07_reassembly" "tmpdir" "read-concatenation" "${PIPELINE_TMPDIR_BASE}" "${r_group}" "stage-precheck" || return 1
            check_disk_requirement "07_reassembly" "outdir" "reassembly-output" "${REASSEMBLY_DIR}" "${r_group}" "stage-precheck" || return 1
            if [[ "${SHARED_CAPACITY_FS}" == "true" ]]; then
                local req_concat
                req_concat="$(estimate_stage_disk_requirement 07_reassembly tmpdir | awk '{print $1}')"
                check_disk_requirement "07_reassembly" "other" "shared-capacity" "${OUTDIR}" "${r_group}" "stage-precheck" "${req_concat}" || return 1
            fi
            ;;
        08_checkm2)
            check_disk_requirement "08_checkm2" "outdir" "quality-assessment" "${CHECKM2_DIR}" "none" "stage-precheck" || return 1
            ;;
        10_gtdbtk)
            check_disk_requirement "10_gtdbtk" "outdir" "classification-output" "${GTDBTK_DIR}" "none" "stage-precheck" || return 1
            ;;
        11_coverm)
            check_disk_requirement "11_coverm" "outdir" "coverage-output" "${ABUNDANCE_DIR}" "none" "stage-precheck" || return 1
            ;;
        *)
            return 0
            ;;
    esac
    return 0
}

# record_post_stage_observation STAGE_ID
#   Records filesystem capacity observation following successful stage completion.
record_post_stage_observation() {
    local stage="$1"
    if [[ "${DRY_RUN}" == true ]]; then return 0; fi

    local check_path="${OUTDIR}"
    local existing_parent
    existing_parent="$(find_existing_parent "${check_path}")"
    local fs_id
    fs_id="$(get_filesystem_id "${existing_parent}")"

    local cap_info tot_bytes used_bytes avail_bytes
    cap_info="$(get_filesystem_capacity "${existing_parent}")"
    tot_bytes="$(awk '{print $1}' <<< "${cap_info}")"
    used_bytes="$(awk '{print $2}' <<< "${cap_info}")"
    avail_bytes="$(awk '{print $3}' <<< "${cap_info}")"

    local avail_inodes
    avail_inodes="$(get_available_inodes "${existing_parent}")"

    local ts_utc
    ts_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "1970-01-01T00:00:00Z")"

    record_disk_measurement \
        --report-current "${DISK_SPACE_CURRENT_TSV}" \
        --publish-path "${DISK_SPACE_TSV}" \
        --run-id "${RUN_ID}" \
        --measurement-sequence "auto" \
        --check-id "${stage}.outdir.observation" \
        --timestamp-utc "${ts_utc}" \
        --stage-id "${stage}" \
        --storage-role "outdir" \
        --purpose "filesystem-observation" \
        --path "${check_path}" \
        --filesystem-id "${fs_id}" \
        --measurement-type "post-stage-observation" \
        --estimate-confidence "not-applicable" \
        --input-bytes-read "NA" \
        --estimated-required-bytes "NA" \
        --filesystem-available-bytes "${avail_bytes}" \
        --filesystem-used-bytes "${used_bytes}" \
        --filesystem-total-bytes "${tot_bytes}" \
        --required-free-inodes "NA" \
        --available-free-inodes "${avail_inodes}" \
        --safety-factor "NA" \
        --shared-capacity-group "none" \
        --status "sufficient" \
        --message "post_stage_observation=stage_finished;stage_id=${stage}"

    if [[ "${stage}" == "07_reassembly" && -n "${PIPELINE_TMPDIR_BASE}" ]]; then
        local tmp_parent tmp_fs_id tmp_cap_info t_tot t_used t_avail t_inodes
        tmp_parent="$(find_existing_parent "${PIPELINE_TMPDIR_BASE}")"
        tmp_fs_id="$(get_filesystem_id "${tmp_parent}")"
        tmp_cap_info="$(get_filesystem_capacity "${tmp_parent}")"
        t_tot="$(awk '{print $1}' <<< "${tmp_cap_info}")"
        t_used="$(awk '{print $2}' <<< "${tmp_cap_info}")"
        t_avail="$(awk '{print $3}' <<< "${tmp_cap_info}")"
        t_inodes="$(get_available_inodes "${tmp_parent}")"

        record_disk_measurement \
            --report-current "${DISK_SPACE_CURRENT_TSV}" \
            --publish-path "${DISK_SPACE_TSV}" \
            --run-id "${RUN_ID}" \
            --measurement-sequence "auto" \
            --check-id "${stage}.tmpdir.observation" \
            --timestamp-utc "${ts_utc}" \
            --stage-id "${stage}" \
            --storage-role "tmpdir" \
            --purpose "filesystem-observation" \
            --path "${PIPELINE_TMPDIR_BASE}" \
            --filesystem-id "${tmp_fs_id}" \
            --measurement-type "post-stage-observation" \
            --estimate-confidence "not-applicable" \
            --input-bytes-read "NA" \
            --estimated-required-bytes "NA" \
            --filesystem-available-bytes "${t_avail}" \
            --filesystem-used-bytes "${t_used}" \
            --filesystem-total-bytes "${t_tot}" \
            --required-free-inodes "NA" \
            --available-free-inodes "${t_inodes}" \
            --safety-factor "NA" \
            --shared-capacity-group "none" \
            --status "sufficient" \
            --message "post_stage_observation=tmpdir_observation;stage_id=${stage}"
    fi
}

# record_final_disk_observations
#   Records final filesystem observation for OUTDIR and TMPDIR at pipeline completion.
record_final_disk_observations() {
    if [[ "${DRY_RUN}" == true ]]; then return 0; fi

    local check_path="${OUTDIR}"
    local existing_parent
    existing_parent="$(find_existing_parent "${check_path}")"
    local fs_id
    fs_id="$(get_filesystem_id "${existing_parent}")"

    local cap_info tot_bytes used_bytes avail_bytes
    cap_info="$(get_filesystem_capacity "${existing_parent}")"
    tot_bytes="$(awk '{print $1}' <<< "${cap_info}")"
    used_bytes="$(awk '{print $2}' <<< "${cap_info}")"
    avail_bytes="$(awk '{print $3}' <<< "${cap_info}")"

    local avail_inodes
    avail_inodes="$(get_available_inodes "${existing_parent}")"

    local ts_utc
    ts_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "1970-01-01T00:00:00Z")"

    record_disk_measurement \
        --report-current "${DISK_SPACE_CURRENT_TSV}" \
        --publish-path "${DISK_SPACE_TSV}" \
        --run-id "${RUN_ID}" \
        --measurement-sequence "auto" \
        --check-id "pipeline.outdir.final" \
        --timestamp-utc "${ts_utc}" \
        --stage-id "pipeline" \
        --storage-role "outdir" \
        --purpose "filesystem-observation" \
        --path "${check_path}" \
        --filesystem-id "${fs_id}" \
        --measurement-type "final-observation" \
        --estimate-confidence "not-applicable" \
        --input-bytes-read "NA" \
        --estimated-required-bytes "NA" \
        --filesystem-available-bytes "${avail_bytes}" \
        --filesystem-used-bytes "${used_bytes}" \
        --filesystem-total-bytes "${tot_bytes}" \
        --required-free-inodes "NA" \
        --available-free-inodes "${avail_inodes}" \
        --safety-factor "NA" \
        --shared-capacity-group "none" \
        --status "sufficient" \
        --message "final_observation=pipeline_completed"

    if [[ -n "${PIPELINE_TMPDIR_BASE}" ]]; then
        local tmp_parent tmp_fs_id tmp_cap_info t_tot t_used t_avail t_inodes
        tmp_parent="$(find_existing_parent "${PIPELINE_TMPDIR_BASE}")"
        tmp_fs_id="$(get_filesystem_id "${tmp_parent}")"
        tmp_cap_info="$(get_filesystem_capacity "${tmp_parent}")"
        t_tot="$(awk '{print $1}' <<< "${tmp_cap_info}")"
        t_used="$(awk '{print $2}' <<< "${tmp_cap_info}")"
        t_avail="$(awk '{print $3}' <<< "${tmp_cap_info}")"
        t_inodes="$(get_available_inodes "${tmp_parent}")"

        record_disk_measurement \
            --report-current "${DISK_SPACE_CURRENT_TSV}" \
            --publish-path "${DISK_SPACE_TSV}" \
            --run-id "${RUN_ID}" \
            --measurement-sequence "auto" \
            --check-id "pipeline.tmpdir.final" \
            --timestamp-utc "${ts_utc}" \
            --stage-id "pipeline" \
            --storage-role "tmpdir" \
            --purpose "filesystem-observation" \
            --path "${PIPELINE_TMPDIR_BASE}" \
            --filesystem-id "${tmp_fs_id}" \
            --measurement-type "final-observation" \
            --estimate-confidence "not-applicable" \
            --input-bytes-read "NA" \
            --estimated-required-bytes "NA" \
            --filesystem-available-bytes "${t_avail}" \
            --filesystem-used-bytes "${t_used}" \
            --filesystem-total-bytes "${t_tot}" \
            --required-free-inodes "NA" \
            --available-free-inodes "${t_inodes}" \
            --safety-factor "NA" \
            --shared-capacity-group "none" \
            --status "sufficient" \
            --message "final_observation=tmpdir_status"
    fi

    publish_disk_space_report
}

# run_disk_preflight
#   Constructs complete preflight assessment for OUTDIR, TMPDIR, and applicable stages.
run_disk_preflight() {
    local preflight_errors=0

    # 1. OUTDIR resolution & initial measurement
    OUTDIR="$(canonicalize_storage_path "${OUTDIR}")"
    local outdir_parent
    outdir_parent="$(find_existing_parent "${OUTDIR}")"
    local parent_fs_id
    parent_fs_id="$(get_filesystem_id "${outdir_parent}")"

    local parent_cap p_avail
    parent_cap="$(get_filesystem_capacity "${outdir_parent}")"
    p_avail="$(awk '{print $3}' <<< "${parent_cap}")"

    # Step 3: verify basic minimum space exists to create basic structure
    if [[ "${p_avail}" != "NA" && "${p_avail}" -le 0 ]]; then
        log_error "No disk space available in ${outdir_parent} (${p_avail} bytes) to create OUTDIR structure."
        preflight_errors=$(( preflight_errors + 1 ))
    fi

    # Step 4 & 5: If normal execution, create OUTDIR and measure real filesystem
    if [[ "${DRY_RUN}" != true ]]; then
        mkdir -p "${OUTDIR}" "${LOGS_DIR}" || {
            log_error "Cannot create OUTDIR or LOGS_DIR: ${OUTDIR}"
            return 1
        }
        OUTDIR_FS_ID="$(get_filesystem_id "${OUTDIR}")"
        if [[ "${OUTDIR_FS_ID}" != "${parent_fs_id}" ]]; then
            log_info "[Disk] OUTDIR filesystem resolved to ${OUTDIR_FS_ID} (inferred from parent was ${parent_fs_id})."
        fi
    else
        OUTDIR_FS_ID="${parent_fs_id}"
    fi

    # TMPDIR resolution
    local tmpdir_target="${PIPELINE_TMPDIR_BASE:-${OUTDIR}/.tmp}"
    local tmpdir_parent
    tmpdir_parent="$(find_existing_parent "${tmpdir_target}")"
    TMPDIR_FS_ID="$(get_filesystem_id "${tmpdir_parent}")"

    # Detect shared capacity
    SHARED_CAPACITY_FS="$(compare_filesystem_identities "${OUTDIR_FS_ID}" "${TMPDIR_FS_ID}")"

    # Pipeline general checks
    check_disk_requirement "pipeline" "outdir" "pipeline-outputs" "${OUTDIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
    check_disk_requirement "pipeline" "tmpdir" "workspace" "${tmpdir_target}" "none" "preflight-estimate" || (( preflight_errors += 1 ))

    local largest_blocking_bytes=0 largest_blocking_stg="none"
    local largest_info_bytes=0 largest_info_stg="none"

    # Stage prechecks
    local stg action
    for stg in "${DISK_PRECHECK_STAGES[@]}"; do
        action="${STAGE_ACTION[$stg]:-disabled}"
        if [[ "${action}" == "reuse" ]]; then
            if [[ "${DRY_RUN}" != true ]]; then
                record_disk_measurement \
                    --report-current "${DISK_SPACE_CURRENT_TSV}" \
                    --publish-path "${DISK_SPACE_TSV}" \
                    --run-id "${RUN_ID}" \
                    --measurement-sequence "auto" \
                    --check-id "${stg}.outdir.execution" \
                    --timestamp-utc "$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "1970-01-01T00:00:00Z")" \
                    --stage-id "${stg}" \
                    --storage-role "outdir" \
                    --purpose "reused-checkpoint" \
                    --path "${OUTDIR}" \
                    --filesystem-id "${OUTDIR_FS_ID}" \
                    --measurement-type "preflight-estimate" \
                    --estimate-confidence "not-applicable" \
                    --input-bytes-read "NA" \
                    --estimated-required-bytes "NA" \
                    --filesystem-available-bytes "NA" \
                    --filesystem-used-bytes "NA" \
                    --filesystem-total-bytes "NA" \
                    --required-free-inodes "NA" \
                    --available-free-inodes "NA" \
                    --safety-factor "NA" \
                    --shared-capacity-group "none" \
                    --status "not-applicable" \
                    --message "action=reuse;no_new_space_required"
            fi
            continue
        fi

        if ! stage_action_requires_execution "${stg}"; then
            continue
        fi

        # Track largest requirements
        local est_line est_req est_block
        est_line="$(estimate_stage_disk_requirement "${stg}" "outdir")"
        est_req="$(awk '{print $1}' <<< "${est_line}")"
        est_block="$(awk '{print $5}' <<< "${est_line}")"

        if [[ "${est_req}" != "NA" && "${est_req}" =~ ^[0-9]+$ ]]; then
            if [[ "${est_block}" == "true" ]]; then
                if (( est_req > largest_blocking_bytes )); then
                    largest_blocking_bytes="${est_req}"
                    largest_blocking_stg="${stg}"
                fi
            else
                if (( est_req > largest_info_bytes )); then
                    largest_info_bytes="${est_req}"
                    largest_info_stg="${stg}"
                fi
            fi
        fi

        local r_group="none"
        [[ "${SHARED_CAPACITY_FS}" == "true" ]] && r_group="${stg}"

        case "${stg}" in
            01_fastp)
                check_disk_requirement "01_fastp" "outdir" "cleaned-reads" "${QC_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            02_megahit)
                check_disk_requirement "02_megahit" "outdir" "assembly" "${ASSEMBLY_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            03_kraken2)
                check_disk_requirement "03_kraken2" "outdir" "taxonomy-output" "${KRAKEN2_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            04_mapping)
                check_disk_requirement "04_mapping" "outdir" "mapping-bam" "${MAPPING_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            05a_vamb)
                check_disk_requirement "05a_vamb" "outdir" "binning" "${VAMB_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            05b_metabinner)
                check_disk_requirement "05b_metabinner" "outdir" "binning" "${METABINNER_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            05c_comebin)
                check_disk_requirement "05c_comebin" "outdir" "binning" "${COMEBIN_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            06_refinement)
                check_disk_requirement "06_refinement" "outdir" "refinement" "${REFINEMENT_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            07_reassembly)
                local est_tmp_line est_tmp_req est_tmp_block
                est_tmp_line="$(estimate_stage_disk_requirement "07_reassembly" "tmpdir")"
                est_tmp_req="$(awk '{print $1}' <<< "${est_tmp_line}")"
                est_tmp_block="$(awk '{print $5}' <<< "${est_tmp_line}")"
                if [[ "${est_tmp_req}" != "NA" && "${est_tmp_req}" =~ ^[0-9]+$ && "${est_tmp_block}" == "true" ]]; then
                    if (( est_tmp_req > largest_blocking_bytes )); then
                        largest_blocking_bytes="${est_tmp_req}"
                        largest_blocking_stg="07_reassembly"
                    fi
                fi
                check_disk_requirement "07_reassembly" "tmpdir" "read-concatenation" "${PIPELINE_TMPDIR_BASE}" "${r_group}" "preflight-estimate" || (( preflight_errors += 1 ))
                check_disk_requirement "07_reassembly" "outdir" "reassembly-output" "${REASSEMBLY_DIR}" "${r_group}" "preflight-estimate" || (( preflight_errors += 1 ))
                if [[ "${SHARED_CAPACITY_FS}" == "true" ]]; then
                    local req_concat
                    req_concat="$(estimate_stage_disk_requirement 07_reassembly tmpdir | awk '{print $1}')"
                    check_disk_requirement "07_reassembly" "other" "shared-capacity" "${OUTDIR}" "${r_group}" "preflight-estimate" "${req_concat}" || (( preflight_errors += 1 ))
                fi
                ;;
            08_checkm2)
                check_disk_requirement "08_checkm2" "outdir" "quality-assessment" "${CHECKM2_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            10_gtdbtk)
                check_disk_requirement "10_gtdbtk" "outdir" "classification-output" "${GTDBTK_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
            11_coverm)
                check_disk_requirement "11_coverm" "outdir" "coverage-output" "${ABUNDANCE_DIR}" "none" "preflight-estimate" || (( preflight_errors += 1 ))
                ;;
        esac
    done

    # Publish report if normal execution
    if [[ "${DRY_RUN}" != true ]]; then
        publish_disk_space_report
    fi

    # Determine summary status
    local sum_out
    sum_out="$(get_latest_disk_status)"
    local overall_st warn_cnt
    overall_st="$(grep -E '^OVERALL_STATUS=' <<< "${sum_out}" | cut -d'=' -f2)"
    warn_cnt="$(grep -E '^WARNING_COUNT=' <<< "${sum_out}" | cut -d'=' -f2)"
    DISK_PREFLIGHT_STATUS="${overall_st:-unknown}"
    DISK_WARNING_COUNT="${warn_cnt:-0}"

    local outdir_avail tmpdir_avail
    outdir_avail="$(get_filesystem_capacity "${OUTDIR}" | awk '{print $3}')"
    tmpdir_avail="$(get_filesystem_capacity "${tmpdir_target}" | awk '{print $3}')"

    local largest_block_disp="none:NA"
    if (( largest_blocking_bytes > 0 )); then
        largest_block_disp="${largest_blocking_stg}:${largest_blocking_bytes}"
    fi
    local largest_info_disp="none:NA"
    if (( largest_info_bytes > 0 )); then
        largest_info_disp="${largest_info_stg}:${largest_info_bytes}"
    fi

    # Print visible block
    if [[ "${DRY_RUN}" == true ]]; then
        local in_bytes
        in_bytes="$(get_canonical_reads_bytes)"
        echo "Disk report                    : would write"
        echo "OUTDIR filesystem              : ${OUTDIR_FS_ID}"
        echo "OUTDIR available bytes         : ${outdir_avail}"
        echo "TMPDIR filesystem              : ${TMPDIR_FS_ID}"
        echo "TMPDIR available bytes         : ${tmpdir_avail}"
        echo "Shared filesystem              : ${SHARED_CAPACITY_FS}"
        echo "Estimated input bytes          : ${in_bytes}"
        echo "Largest blocking requirement   : ${largest_block_disp}"
        echo "Largest informational estimate : ${largest_info_disp}"
        echo "Disk status                    : ${DISK_PREFLIGHT_STATUS}"
    else
        echo "Disk report                    : ${DISK_SPACE_TSV}"
        echo "OUTDIR filesystem              : ${OUTDIR_FS_ID}"
        echo "OUTDIR available bytes         : ${outdir_avail}"
        echo "TMPDIR filesystem              : ${TMPDIR_FS_ID}"
        echo "TMPDIR available bytes         : ${tmpdir_avail}"
        echo "Shared filesystem              : ${SHARED_CAPACITY_FS}"
        echo "Largest blocking requirement   : ${largest_block_disp}"
        echo "Largest informational estimate : ${largest_info_disp}"
        echo "Disk preflight status          : ${DISK_PREFLIGHT_STATUS}"
    fi

    if (( preflight_errors > 0 )); then
        return 1
    fi
    return 0
}

# create_pipeline_run_tmpdir
#   Creates PIPELINE_RUN_TMPDIR and REASSEMBLY_TMPDIR exclusively (no -p).
#   Writes the ownership marker.
create_pipeline_run_tmpdir() {
    local parent="${PIPELINE_TMP_PARENT}"
    local run_dir="${PIPELINE_RUN_TMPDIR}"
    local reassembly_dir="${REASSEMBLY_TMPDIR}"
    local marker="${PIPELINE_TMP_MARKER}"

    if [[ -z "${run_dir}" || -z "${reassembly_dir}" || -z "${marker}" ]]; then
        log_error "[Tmpdir] PIPELINE_RUN_TMPDIR / REASSEMBLY_TMPDIR / PIPELINE_TMP_MARKER not set."
        return 1
    fi

    # Create the shared parent (idempotent)
    mkdir -p -- "${parent}" || {
        log_error "[Tmpdir] Cannot create parent directory: ${parent}"
        return 1
    }

    # Create the per-run directory exclusively (fail if it exists)
    local old_umask
    old_umask="$(umask)"
    umask 077
    if ! mkdir --mode=700 -- "${run_dir}" 2>/dev/null; then
        umask "${old_umask}"
        log_error "[Tmpdir] Per-run workspace already exists — refusing to reuse."
        log_error "         Path: ${run_dir}"
        log_error "         Inspect it manually. If it is orphaned, remove it and retry."
        return 1
    fi
    mkdir --mode=700 -- "${reassembly_dir}" || {
        umask "${old_umask}"
        log_error "[Tmpdir] Cannot create reassembly subdirectory: ${reassembly_dir}"
        return 1
    }
    umask "${old_umask}"

    # Write ownership marker (atomic via temporary file)
    local ts_utc
    ts_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "1970-01-01T00:00:00Z")"
    local marker_tmp="${marker}.tmp.$$"
    {
        printf 'schema_version=1\n'
        printf 'run_id=%s\n'          "${RUN_ID}"
        printf 'pipeline_pid=%s\n'    "$$"
        printf 'created_at=%s\n'      "${ts_utc}"
        printf 'tmpdir_source=%s\n'   "${PIPELINE_TMPDIR_SOURCE}"
        printf 'tmpdir_base=%s\n'     "${PIPELINE_TMPDIR_BASE}"
        printf 'run_tmpdir=%s\n'      "${run_dir}"
    } > "${marker_tmp}" || {
        rm -f -- "${marker_tmp}" 2>/dev/null || true
        log_error "[Tmpdir] Cannot write ownership marker (tmp): ${marker_tmp}"
        return 1
    }
    chmod 600 -- "${marker_tmp}" || true
    mv -- "${marker_tmp}" "${marker}" || {
        rm -f -- "${marker_tmp}" 2>/dev/null || true
        log_error "[Tmpdir] Cannot install ownership marker: ${marker}"
        return 1
    }

    log_info "[Tmpdir] Workspace created: ${run_dir}"
    log_info "[Tmpdir] Reassembly tmpdir : ${reassembly_dir}"
    log_info "[Tmpdir] Marker            : ${marker}"
    return 0
}

# validate_run_tmpdir_for_cleanup
#   Returns 0 iff PIPELINE_RUN_TMPDIR is safe to remove:
#   - non-empty, non-symlink directory
#   - absolute path, not /, not HOME, not OUTDIR, not base, not parent
#   - exact match: ${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline/${RUN_ID}
#   - reassembly subdir is not a symlink
#   - ownership marker exists, is regular file, not symlink
#   - ownership marker fields match RUN_ID, PID, base, and run_tmpdir
#   - canonical path matches
#   - concatenated reads/parts are not symlinks
validate_run_tmpdir_for_cleanup() {
    local run_dir="${PIPELINE_RUN_TMPDIR:-}"
    local marker="${PIPELINE_TMP_MARKER:-}"
    local parent="${PIPELINE_TMP_PARENT:-}"
    local base="${PIPELINE_TMPDIR_BASE:-}"
    local reassembly_dir="${REASSEMBLY_TMPDIR:-}"

    [[ -z "${run_dir}" ]] && return 1
    [[ -z "${marker}" ]]  && return 1
    [[ -z "${parent}" ]]  && return 1
    [[ -z "${base}" ]]    && return 1

    # Absolute path check
    [[ "${run_dir:0:1}" == "/" ]] || return 1
    # Forbidden roots
    [[ "${run_dir}" == "/" ]] && return 1
    [[ -n "${HOME:-}" && "${run_dir}" == "${HOME}" ]] && return 1

    local outdir_abs=""
    outdir_abs="$(realpath -m -- "${OUTDIR}" 2>/dev/null || printf '%s' "${OUTDIR}")"
    [[ -n "${outdir_abs}" && "${run_dir}" == "${outdir_abs}" ]] && return 1
    [[ "${run_dir}" == "${base}" ]] && return 1
    [[ "${run_dir}" == "${parent}" ]] && return 1

    # Exact pattern match
    local expected_run_dir="${base}/metagenomics_pipeline/${RUN_ID}"
    if [[ "${run_dir}" != "${expected_run_dir}" ]]; then
        log_warn "[Tmpdir] Workspace path does not match expected pattern — skipping cleanup."
        log_warn "         Expected: ${expected_run_dir}"
        log_warn "         Actual  : ${run_dir}"
        return 1
    fi

    # Symlink checks on directories
    [[ -L "${run_dir}" ]] && { log_warn "[Tmpdir] Workspace is a symlink — skipping cleanup."; return 1; }
    [[ -d "${run_dir}" ]] || return 1
    if [[ -n "${reassembly_dir}" && -e "${reassembly_dir}" ]]; then
        [[ -L "${reassembly_dir}" ]] && { log_warn "[Tmpdir] Reassembly subdir is a symlink — skipping cleanup."; return 1; }
    fi

    # Canonical path equality
    local canon_run=""
    canon_run="$(realpath -- "${run_dir}" 2>/dev/null || printf '%s' "${run_dir}")"
    if [[ "${canon_run}" != "${run_dir}" ]]; then
        log_warn "[Tmpdir] Workspace resolves to a different canonical path — skipping cleanup."
        return 1
    fi

    # Verify parent is direct child
    local actual_parent="${run_dir%/*}"
    local canon_parent=""
    canon_parent="$(realpath -- "${parent}" 2>/dev/null || printf '%s' "${parent}")"
    local canon_actual=""
    canon_actual="$(realpath -- "${actual_parent}" 2>/dev/null || printf '%s' "${actual_parent}")"
    if [[ "${canon_actual}" != "${canon_parent}" ]]; then
        log_warn "[Tmpdir] Workspace parent mismatch — skipping cleanup."
        return 1
    fi

    # Ownership marker checks
    [[ -f "${marker}" && ! -L "${marker}" ]] || { log_warn "[Tmpdir] Ownership marker missing or symlink."; return 1; }

    local m_run_id="" m_pid="" m_base="" m_run_tmpdir=""
    while IFS='=' read -r k v; do
        case "${k}" in
            run_id)       m_run_id="${v}"       ;;
            pipeline_pid) m_pid="${v}"          ;;
            tmpdir_base)  m_base="${v}"         ;;
            run_tmpdir)   m_run_tmpdir="${v}"   ;;
        esac
    done < "${marker}"

    if [[ "${m_run_id}" != "${RUN_ID}" ]]; then
        log_warn "[Tmpdir] Ownership marker RUN_ID mismatch — refusing cleanup (marker: ${m_run_id}, expected: ${RUN_ID})."
        return 1
    fi
    if [[ "${m_pid}" != "$$" ]]; then
        log_warn "[Tmpdir] Ownership marker PID mismatch — refusing cleanup (marker: ${m_pid}, expected: $$)."
        return 1
    fi
    if [[ -n "${m_base}" && "${m_base}" != "${base}" ]]; then
        log_warn "[Tmpdir] Ownership marker base mismatch — refusing cleanup."
        return 1
    fi
    if [[ -n "${m_run_tmpdir}" && "${m_run_tmpdir}" != "${run_dir}" ]]; then
        log_warn "[Tmpdir] Ownership marker run_tmpdir mismatch — refusing cleanup."
        return 1
    fi

    # Concatenated files symlink checks
    if [[ -n "${reassembly_dir}" && -d "${reassembly_dir}" ]]; then
        [[ -L "${reassembly_dir}/concatenated_R1.fastq.gz" ]] && return 1
        [[ -L "${reassembly_dir}/concatenated_R2.fastq.gz" ]] && return 1
        [[ -L "${reassembly_dir}/concatenated_R1.fastq.gz.part" ]] && return 1
        [[ -L "${reassembly_dir}/concatenated_R2.fastq.gz.part" ]] && return 1
    fi

    return 0
}

# cleanup_run_tmpdir
#   Idempotent wrapper that validates and removes PIPELINE_RUN_TMPDIR.
cleanup_run_tmpdir() {
    local run_dir="${PIPELINE_RUN_TMPDIR:-}"
    [[ -z "${run_dir}" || ! -d "${run_dir}" ]] && return 0
    if ! validate_run_tmpdir_for_cleanup; then
        return 1
    fi
    rm -rf -- "${run_dir}"
}

# cleanup_safe_temporaries CONTEXT
#   Removes PIPELINE_RUN_TMPDIR when it is safe to do so.
#   CONTEXT: success | failure | signal
#   Signal context always preserves (cannot guarantee child processes exited).
#   When KEEP_TEMPORARY=true always preserves.
cleanup_safe_temporaries() {
    local context="${1:-unknown}"
    local run_dir="${PIPELINE_RUN_TMPDIR:-}"

    [[ -z "${run_dir}" ]] && return 0
    [[ ! -d "${run_dir}" ]] && return 0

    if [[ "${context}" == "signal" ]]; then
        log_warn "[Tmpdir] Signal received — preserving workspace (cannot guarantee children exited)."
        log_warn "         Preserved at: ${run_dir}"
        log_warn "         Inspect and remove manually when safe."
        TEMPORARY_CLEANUP_STATUS="preserved-after-signal"
        TEMPORARY_FILES_PRESERVED=true
        return 0
    fi

    if [[ "${KEEP_TEMPORARY}" == true ]]; then
        log_info "[Tmpdir] KEEP_TEMPORARY=true — preserving workspace: ${run_dir}"
        TEMPORARY_CLEANUP_STATUS="preserved-by-flag"
        TEMPORARY_FILES_PRESERVED=true
        return 0
    fi

    if [[ "${context}" == "failure" ]]; then
        # Check if child processes of this PID are still active
        if pgrep -P "$$" >/dev/null 2>&1; then
            log_warn "[Tmpdir] Child processes still active — preserving workspace: ${run_dir}"
            log_warn "         Inspect and remove manually when safe."
            TEMPORARY_CLEANUP_STATUS="preserved-active-children"
            TEMPORARY_FILES_PRESERVED=true
            return 0
        fi
    fi

    if ! validate_run_tmpdir_for_cleanup; then
        log_warn "[Tmpdir] Cleanup validation failed — preserving workspace: ${run_dir}"
        log_warn "         Inspect and remove manually when safe."
        TEMPORARY_CLEANUP_STATUS="validation-failed"
        TEMPORARY_FILES_PRESERVED=true
        return 0
    fi

    log_info "[Tmpdir] Removing temporary workspace (context=${context}): ${run_dir}"
    if rm -rf -- "${run_dir}"; then
        log_info "[Tmpdir] Temporary workspace removed."
        TEMPORARY_CLEANUP_STATUS="removed"
        TEMPORARY_FILES_PRESERVED=false
    else
        log_warn "[Tmpdir] rm -rf failed for workspace: ${run_dir}"
        TEMPORARY_CLEANUP_STATUS="removal-failed"
        TEMPORARY_FILES_PRESERVED=true
    fi
}

# show_tmpdir_dry_run_info
#   Prints temporary workspace plan during dry-run (no directories created).
show_tmpdir_dry_run_info() {
    local total_in
    total_in="$(get_canonical_reads_bytes)"
    local factor="${REASSEMBLY_CONCAT_SPACE_FACTOR:-1.10}"
    local est_bytes
    est_bytes="$(awk -v b="${total_in}" -v f="${factor}" 'BEGIN {
        val = b * f
        ival = int(val)
        if (val > ival) ival++
        printf "%.0f", ival
    }')"

    local check_path="${PIPELINE_TMPDIR_BASE}"
    local parent
    parent="$(find_existing_parent "${check_path}")"
    local cap_info avail_bytes_str="unknown"
    cap_info="$(get_filesystem_capacity "${parent}")"
    local a_bytes
    a_bytes="$(awk '{print $3}' <<< "${cap_info}")"
    if [[ "${a_bytes}" != "NA" ]]; then
        avail_bytes_str="${a_bytes}"
    fi

    local avail_inodes_str="unknown"
    local a_inodes
    a_inodes="$(get_available_inodes "${parent}")"
    if [[ "${a_inodes}" != "NA" ]]; then
        avail_inodes_str="${a_inodes}"
    fi

    log_info "[DRY-RUN] Temporary workspace plan:"
    log_info "  Temporary base                   : ${PIPELINE_TMPDIR_BASE:-<unresolved>}"
    log_info "  Temporary source                 : ${PIPELINE_TMPDIR_SOURCE:-<unresolved>}"
    log_info "  Run temporary directory         : would create ${PIPELINE_RUN_TMPDIR:-<unresolved>}"
    log_info "  Reassembly temp strategy         : ${REASSEMBLY_TEMP_READS_STRATEGY}"
    log_info "  Strategy version                 : ${REASSEMBLY_TEMP_READS_STRATEGY_VERSION}"
    log_info "  Temporary validation             : ${REASSEMBLY_TEMP_READS_VALIDATION}"
    log_info "  Estimated concatenation bytes    : ${est_bytes}"
    log_info "  Available temporary bytes        : ${avail_bytes_str}"
    log_info "  Required free inodes             : ${REQUIRED_FREE_INODES:-16}"
    log_info "  Available free inodes            : ${avail_inodes_str}"
    log_info "  MetaWRAP workspace estimate      : not included"
    log_info "  Keep temporary                   : ${KEEP_TEMPORARY}"
    log_info "  Reassembly bind                  : would bind read-write ${PIPELINE_RUN_TMPDIR:-<unresolved>}"
}

# prepare_reassembly_tmpdir R1_RESULT_VAR R2_RESULT_VAR
#   Validates inputs, checks space, and concatenates paired-end reads into REASSEMBLY_TMPDIR.
#   Uses binary (gzip multi-member) concatenation with .part + gzip -t.
#   Sets the named result variables to the final paths.
prepare_reassembly_tmpdir() {
    local _r1_var="$1"
    local _r2_var="$2"

    local r1_out="${REASSEMBLY_TMPDIR}/concatenated_R1.fastq.gz"
    local r2_out="${REASSEMBLY_TMPDIR}/concatenated_R2.fastq.gz"
    local r1_part="${r1_out}.part"
    local r2_part="${r2_out}.part"

    # Step 1: Validate active sample arrays from manifest
    local n_ids="${#ACTIVE_SAMPLE_IDS[@]}"
    local n_r1="${#ACTIVE_R1_FILES[@]}"
    local n_r2="${#ACTIVE_R2_FILES[@]}"
    if (( n_ids == 0 )); then
        log_error "[Step 07] No active samples found in manifest for reassembly."
        return 1
    fi
    if ! (( n_ids == n_r1 && n_ids == n_r2 )); then
        log_error "[Step 07] Sample count mismatch in active manifest: ${n_ids} IDs, ${n_r1} R1s, ${n_r2} R2s."
        return 1
    fi

    local i r1_f r2_f
    for (( i=0; i<${#ACTIVE_SAMPLE_IDS[@]}; i++ )); do
        r1_f="${ACTIVE_R1_FILES[i]}"
        r2_f="${ACTIVE_R2_FILES[i]}"
        if [[ ! -f "${r1_f}" || ! -r "${r1_f}" ]]; then
            log_error "[Step 07] Active R1 file missing or unreadable: ${r1_f}"
            return 1
        fi
        if [[ ! -f "${r2_f}" || ! -r "${r2_f}" ]]; then
            log_error "[Step 07] Active R2 file missing or unreadable: ${r2_f}"
            return 1
        fi
        if [[ "${r1_f}" == "${r2_f}" ]]; then
            log_error "[Step 07] R1 and R2 are identical path for sample ${ACTIVE_SAMPLE_IDS[i]}: ${r1_f}"
            return 1
        fi
    done

    # Step 2: Check that final files do not already exist (no accidental reuse)
    if [[ -e "${r1_out}" || -e "${r2_out}" ]]; then
        log_error "[Step 07] Reassembly temporary file already exists — refusing to overwrite."
        log_error "         R1: ${r1_out}"
        log_error "         R2: ${r2_out}"
        return 1
    fi

    # Step 3: Space pre-flight
    check_tmpdir_space ACTIVE_R1_FILES ACTIVE_R2_FILES || return 1

    # Step 4: Concatenate R1 in manifest order to .part
    log_info "[Step 07] Concatenating ${#ACTIVE_R1_FILES[@]} R1 reads into temporary ${r1_part}..."
    rm -f -- "${r1_part}"
    if ! cat -- "${ACTIVE_R1_FILES[@]}" > "${r1_part}"; then
        log_error "[Step 07] Failed during R1 concatenation."
        rm -f -- "${r1_part}" 2>/dev/null || true
        return 1
    fi

    # Step 5: Concatenate R2 in manifest order to .part
    log_info "[Step 07] Concatenating ${#ACTIVE_R2_FILES[@]} R2 reads into temporary ${r2_part}..."
    rm -f -- "${r2_part}"
    if ! cat -- "${ACTIVE_R2_FILES[@]}" > "${r2_part}"; then
        log_error "[Step 07] Failed during R2 concatenation."
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    fi

    # Step 6: Verify both part files are non-empty
    if [[ ! -s "${r1_part}" ]]; then
        log_error "[Step 07] Concatenated R1 part file is empty: ${r1_part}"
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    fi
    if [[ ! -s "${r2_part}" ]]; then
        log_error "[Step 07] Concatenated R2 part file is empty: ${r2_part}"
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    fi

    # Step 7: Validate gzip integrity via gzip -t
    log_info "[Step 07] Validating R1 gzip integrity..."
    if ! gzip -t -- "${r1_part}" 2>/dev/null; then
        log_error "[Step 07] gzip integrity check failed for concatenated R1."
        log_error "          Note: gzip -t checks file integrity only; it does not"
        log_error "          validate that R1 and R2 reads are correctly paired."
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    fi
    log_info "[Step 07] Validating R2 gzip integrity..."
    if ! gzip -t -- "${r2_part}" 2>/dev/null; then
        log_error "[Step 07] gzip integrity check failed for concatenated R2."
        log_error "          Note: gzip -t checks file integrity only; it does not"
        log_error "          validate that R1 and R2 reads are correctly paired."
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    fi

    # Step 8: Atomic rename within same directory
    mv -- "${r1_part}" "${r1_out}" || {
        log_error "[Step 07] Cannot rename R1 part to final path: ${r1_out}"
        rm -f -- "${r1_part}" "${r2_part}" 2>/dev/null || true
        return 1
    }
    mv -- "${r2_part}" "${r2_out}" || {
        log_error "[Step 07] Cannot rename R2 part to final path: ${r2_out}"
        rm -f -- "${r1_out}" "${r2_part}" 2>/dev/null || true
        return 1
    }

    # Step 9: Final product check
    if [[ ! -f "${r1_out}" || ! -s "${r1_out}" || ! -f "${r2_out}" || ! -s "${r2_out}" ]]; then
        log_error "[Step 07] Concatenated final files failed final existence/non-empty check."
        return 1
    fi

    printf -v "${_r1_var}" '%s' "${r1_out}"
    printf -v "${_r2_var}" '%s' "${r2_out}"
    log_info "[Step 07] Concatenated reads ready in temporary workspace:"
    log_info "          R1: ${r1_out}"
    log_info "          R2: ${r2_out}"
    return 0
}

# ── Utility Functions & Wrappers ─────────────────────────────────────────────
run_cmd() {
    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_RUNTIME="${CURRENT_STAGE_RUNTIME:-host}"
    CURRENT_STAGE_COMMAND="$(sanitize_command "$@")"

    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] %s\n' "${CURRENT_STAGE_COMMAND}"
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi

    local cmd_display
    printf -v cmd_display '%q ' "$@"
    log_info "Running: ${cmd_display}"

    local had_errexit=0
    local rc=0
    [[ $- == *e* ]] && had_errexit=1

    set +e
    "$@"
    rc=$?
    (( had_errexit )) && set -e

    if (( rc == 0 )); then
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    return "${rc}"
}

run_cmd_redirect() {
    local output="$1"
    shift
    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_RUNTIME="${CURRENT_STAGE_RUNTIME:-host}"
    CURRENT_STAGE_COMMAND="$(sanitize_command "$@") > $(printf '%q' "${output}")"

    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] %s\n' "${CURRENT_STAGE_COMMAND}"
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi

    local cmd_display
    printf -v cmd_display '%q ' "$@"
    log_info "Running: ${cmd_display}> $(printf '%q' "${output}")"

    local had_errexit=0
    local rc=0
    [[ $- == *e* ]] && had_errexit=1

    set +e
    "$@" > "${output}"
    rc=$?
    (( had_errexit )) && set -e

    if (( rc == 0 )); then
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    return "${rc}"
}

run_cmd_with_tee() {
    local log="$1"
    shift
    [[ "${1:-}" == "--" ]] && shift

    CURRENT_STAGE_LOG="${log}"
    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_RUNTIME="${CURRENT_STAGE_RUNTIME:-host}"
    CURRENT_STAGE_COMMAND="$(sanitize_command "$@")"

    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] %s 2>&1 | tee %q\n' "${CURRENT_STAGE_COMMAND}" "${log}"
        "$@"
        local mock_rc=$?
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return "${mock_rc}"
    fi

    mkdir -p "$(dirname "${log}")"

    local had_errexit=0
    local pipe_rc=0
    [[ $- == *e* ]] && had_errexit=1

    set +e
    if [[ "${GNU_TIME_AVAILABLE:-false}" == true && -n "${CURRENT_STAGE_METRICS_DIR:-}" && -d "${CURRENT_STAGE_METRICS_DIR:-}" ]]; then
        local time_out="${CURRENT_STAGE_METRICS_DIR}/time_${RANDOM}_$$.txt"
        "${GNU_TIME_BIN}" -v -o "${time_out}" "$@" 2>&1 | tee "${log}"
        pipe_rc="${PIPESTATUS[0]}"
    else
        "$@" 2>&1 | tee "${log}"
        pipe_rc="${PIPESTATUS[0]}"
    fi
    (( had_errexit )) && set -e

    if (( pipe_rc == 0 )); then
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    return "${pipe_rc}"
}

# ── run_tool(step_label, sif_path, uri, command...) ───────────────────────────
# Runs a direct binary in its SIF image in Apptainer mode, preserving the
# original direct command path in mamba mode.
run_tool() {
    local step_label="$1"
    local sif_path="$2"
    local uri="$3"
    shift 3

    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_COMMAND="$(sanitize_command "$@")"

    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        CURRENT_STAGE_RUNTIME="host"
        run_cmd "$@"
        return $?
    fi

    ensure_sif_image "${sif_path}" "${uri}" || return 1
    local runtime="apptainer"
    command -v apptainer &>/dev/null || runtime="singularity"
    local -a bind_args=()
    build_bind_mounts bind_args

    local abs_sif
    abs_sif="$(realpath -m "${sif_path}" 2>/dev/null || echo "${sif_path}")"
    CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"

    log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] %q ' "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" "$@"
        printf '\n'
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi

    local had_errexit=0
    local rc=0
    [[ $- == *e* ]] && had_errexit=1

    set +e
    "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" "$@"
    rc=$?
    (( had_errexit )) && set -e

    if (( rc == 0 )); then
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    return "${rc}"
}

# ── run_tool_shell(step_label, sif_path, uri, command) ─────────────────────────
# Runs compound commands (such as the mapping pipeline) inside one SIF image.
run_tool_shell() {
    local step_label="$1"
    local sif_path="$2"
    local uri="$3"
    local command="$4"

    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_COMMAND="$(sanitize_shell_command "${command}")"

    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        CURRENT_STAGE_RUNTIME="host"
        if [[ "${DRY_RUN}" == true ]]; then
            printf '[DRY-RUN] bash -o pipefail -c %q\n' "${command}"
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
            return 0
        fi

        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        bash -o pipefail -c "${command}"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi

        return "${rc}"
    fi

    ensure_sif_image "${sif_path}" "${uri}" || return 1
    local runtime="apptainer"
    command -v apptainer &>/dev/null || runtime="singularity"
    local -a bind_args=()
    build_bind_mounts bind_args

    local abs_sif
    abs_sif="$(realpath -m "${sif_path}" 2>/dev/null || echo "${sif_path}")"
    CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"

    log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] %q ' "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" bash -o pipefail -c "${command}"
        printf '\n'
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi

    local had_errexit=0
    local rc=0
    [[ $- == *e* ]] && had_errexit=1

    set +e
    "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" bash -o pipefail -c "${command}"
    rc=$?
    (( had_errexit )) && set -e

    if (( rc == 0 )); then
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    return "${rc}"
}

# ── ensure_sif_image(sif_path, uri) ──────────────────────────────────────────
# Verifies that a .sif image exists in SIF_DIR. If missing and AUTO_PULL_SIF
# is true, pulls it automatically using apptainer pull or singularity pull.
ensure_sif_image() {
    local sif_path="$1"
    local uri="$2"

    if [[ -f "${sif_path}" ]]; then
        log_info "[Container] Image found: ${sif_path}"
        return 0
    fi

    if [[ "${AUTO_PULL_SIF}" != true ]]; then
        log_error "[Container] Image not found: ${sif_path}"
        log_error "[Container] Set AUTO_PULL_SIF=true or place the image manually."
        return 1
    fi

    local runtime=""
    if command -v apptainer &>/dev/null; then
        runtime="apptainer"
    elif command -v singularity &>/dev/null; then
        runtime="singularity"
    else
        log_error "[Container] Neither 'apptainer' nor 'singularity' found in PATH."
        return 1
    fi

    mkdir -p "$(dirname "${sif_path}")"
    log_info "[Container] Pulling image: ${uri} → ${sif_path}"

    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} ${runtime} pull ${sif_path} ${uri}"
        return 0
    fi

    "${runtime}" pull "${sif_path}" "${uri}"
}

run_in_sif_exec() {
    local sif_file="$1"
    shift
    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    local abs_sif
    abs_sif="$(realpath -m "${sif_file}" 2>/dev/null || echo "${sif_file}")"
    CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"
    CURRENT_STAGE_COMMAND="$(sanitize_command "$@")"

    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} [$(basename "${sif_file}")] $*"
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi
    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        local apptainer_bin
        apptainer_bin=$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null)
        local -a bind_args=()
        build_bind_mounts bind_args

        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        "${apptainer_bin}" exec --cleanenv "${bind_args[@]}" "${sif_file}" "$@"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi
        return "${rc}"
    else
        run_cmd "$@"
    fi
}

run_in_sif_shell() {
    local sif_file="$1"
    local command_str="$2"
    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    local abs_sif
    abs_sif="$(realpath -m "${sif_file}" 2>/dev/null || echo "${sif_file}")"
    CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"
    CURRENT_STAGE_COMMAND="$(sanitize_shell_command "${command_str}")"

    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} [$(basename "${sif_file}")] ${command_str}"
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
        return 0
    fi
    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        local apptainer_bin
        apptainer_bin=$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null)
        local -a bind_args=()
        build_bind_mounts bind_args

        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        "${apptainer_bin}" exec --cleanenv "${bind_args[@]}" "${sif_file}" bash -o pipefail -c "${command_str}"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi
        return "${rc}"
    else
        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        bash -o pipefail -c "${command_str}"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi
        return "${rc}"
    fi
}

run_in_sif() {
    run_in_sif_shell "$1" "${*:2}"
}


build_bind_mounts() {
    local output_name="$1"
    local -n output_args="${output_name}"
    local -a bind_paths=()
    local -A bind_modes=()

    _add_bind() {
        local raw_path="$1" mode="${2:-rw}"
        [[ -z "${raw_path}" ]] && return
        # realpath -m resolves even if the path does not yet exist
        local abs
        abs=$(realpath -m "${raw_path}" 2>/dev/null || echo "${raw_path}")
        # Avoid duplicates
        local path
        for path in "${bind_paths[@]+"${bind_paths[@]}"}"; do
            if [[ "${path}" == "${abs}" ]]; then
                [[ "${mode}" == ro ]] && bind_modes["${path}"]=ro
                return
            fi
        done
        bind_paths+=("${abs}")
        bind_modes["${abs}"]="${mode}"
    }

    _add_bind "${OUTDIR}"
    if [[ "${NEEDS_READS}" == true ]]; then
        local read_file
        if [[ "${RUN_FASTP}" == true && ${#INPUT_R1_FILES[@]} -gt 0 ]]; then
            for read_file in "${INPUT_R1_FILES[@]}" "${INPUT_R2_FILES[@]}"; do
                _add_bind "${read_file}" ro
            done
        elif [[ "${ACTIVE_READS_SOURCE}" == raw || "${ACTIVE_READS_SOURCE}" == external-clean ]]; then
            for read_file in "${ACTIVE_R1_FILES[@]}" "${ACTIVE_R2_FILES[@]}"; do
                _add_bind "${read_file}" ro
            done
        fi
    fi
    [[ "${RUN_CHECKM2}" == true && -n "${CHECKM2_DB}" ]] && stage_in_scope checkm2 && _add_bind "${CHECKM2_DB}" ro
    [[ "${RUN_KRAKEN2}" == true ]] && [[ -d "${KRAKEN2_DB}" ]] && _add_bind "${KRAKEN2_DB}" ro
    [[ "${RUN_GTDBTK}" == true ]] && stage_in_scope gtdbtk && \
        [[ -d "${GTDBTK_DATA_PATH}" ]] && _add_bind "${GTDBTK_DATA_PATH}" ro
    [[ "${RUN_BINNING}" == true ]] && stage_in_scope binning && \
        [[ -d "${METABINNER_SCRIPTS_DIR}" ]] && _add_bind "${METABINNER_SCRIPTS_DIR}" ro
    [[ "${NEEDS_ASSEMBLY}" == true && "${RUN_MEGAHIT}" != true && \
       -n "${ASSEMBLY_IN}" ]] && _add_bind "${ASSEMBLY_IN}" ro
    local binner external_bam
    if [[ "${RUN_BINNING}" == true && "${COVERAGE_SOURCE}" == external ]]; then
        for binner in "${SELECTED_BINNER_ARRAY[@]}"; do
            if binner_requires_input "${binner}" jgi-depth && [[ -n "${JGI_DEPTH_FILE_ACTIVE}" ]]; then
                _add_bind "${JGI_DEPTH_FILE_ACTIVE}" ro
            fi
            if binner_requires_input "${binner}" bam; then
                for external_bam in "${ACTIVE_BAMS[@]}" "${ACTIVE_BAM_INDEXES[@]}"; do
                    [[ -n "${external_bam}" ]] && _add_bind "${external_bam}" ro
                done
            fi
        done
    fi

    local path
    for path in "${bind_paths[@]+"${bind_paths[@]}"}"; do
        if [[ "${bind_modes[$path]}" == ro ]]; then
            output_args+=("-B" "${path}:${path}:ro")
        else
            output_args+=("-B" "${path}:${path}")
        fi
    done
    return 0
}

# ── _env_to_sif(env_name) → sif_path / uri ───────────────────────────────────
# Maps a mamba environment name to its SIF_* variable pair.
_env_to_sif() {
    local env_name="$1"
    case "${env_name}" in
        "${ENV_VAMB}") echo "${SIF_VAMB}" "${SIF_VAMB_URI}" ;;
        "${ENV_METABINNER}") echo "${SIF_METABINNER}" "${SIF_METABINNER_URI}" ;;
        "${ENV_COMEBIN}") echo "${SIF_COMEBIN}" "${SIF_COMEBIN_URI}" ;;
        "${ENV_METAWRAP}") echo "${SIF_METAWRAP}" "${SIF_METAWRAP_URI}" ;;
        "${ENV_CHECKM2}") echo "${SIF_CHECKM2}" "${SIF_CHECKM2_URI}" ;;
        "${ENV_GTDBTK}") echo "${SIF_GTDBTK}" "${SIF_GTDBTK_URI}" ;;
        *)
            log_error "[Container] No SIF mapping for environment: ${env_name}"
            return 1 ;;
    esac
}

# ── _has_gpu() ────────────────────────────────────────────────────────────────
# Returns 0 (true) if at least one NVIDIA GPU is visible on the host.
_has_gpu() {
    command -v nvidia-smi &>/dev/null && nvidia-smi --query-gpu=name \
        --format=csv,noheader &>/dev/null
}

validate_comebin_device() {
    if [[ "${COMEBIN_DEVICE}" == cpu || "${COMEBIN_DEVICE}" == cuda || \
          "${COMEBIN_DEVICE}" =~ ^cuda:[0-9]+$ ]]; then
        return 0
    fi

    log_error "Invalid COMEBin device: ${COMEBIN_DEVICE}"
    log_error "Allowed values: cpu, cuda, cuda:N"
    return 1
}

validate_comebin_preflight() {
    local runtime="apptainer" wrapper_path="" runtime_info="" python_check="" wrapper_check="" version_output=""
    local torch_cuda_requested=false device_index=""
    local -a bind_args=() runtime_args=() visible_inputs=()
    local input

    [[ "${RUN_BINNING}" == true ]] && binner_selected comebin || return 0
    validate_comebin_device || return 1
    stage_in_scope binning || return 0

    if [[ "${COMEBIN_DEVICE}" == cuda || "${COMEBIN_DEVICE}" == cuda:* ]]; then
        torch_cuda_requested=true
    fi
    COMEBIN_APPTAINER_NV=false
    if [[ "${EXEC_MODE}" == apptainer ]]; then
        if [[ "${torch_cuda_requested}" == true ]]; then
            COMEBIN_APPTAINER_NV=true
            runtime_args+=(--nv)
        fi
        command -v apptainer &>/dev/null || runtime="singularity"
        ensure_sif_image "${SIF_COMEBIN}" "${SIF_COMEBIN_URI}" || return 1
        if [[ ! -f "${SIF_COMEBIN}" ]]; then
            log_info "[COMEBin] Runtime inspection deferred; image will be pulled before execution."
            log_info "COMEBin selected       : true"
            log_info "COMEBin version        : ${COMEBIN_VERSION}"
            log_info "COMEBin device         : ${COMEBIN_DEVICE}"
            log_info "COMEBin seed           : ${COMEBIN_SEED}"
            log_info "PyTorch version        : uninspected"
            log_info "PyTorch CUDA build     : uninspected"
            log_info "CUDA available         : uninspected"
            log_info "CUDA device count      : uninspected"
            log_info "Apptainer --nv         : ${COMEBIN_APPTAINER_NV}"
            [[ "${DRY_RUN}" == true ]] && show_comebin_command
            return 0
        fi
        build_bind_mounts bind_args
        if [[ ! -d "${OUTDIR}" ]]; then
            local outdir_abs
            outdir_abs="$(realpath -m "${OUTDIR}")"
            local -a existing_binds=()
            local bind_index
            for (( bind_index=0; bind_index<${#bind_args[@]}; bind_index++ )); do
                if [[ "${bind_args[$bind_index]}" == -B && \
                      "${bind_args[$((bind_index + 1))]:-}" == "${outdir_abs}:${outdir_abs}" ]]; then
                    bind_index=$((bind_index + 1))
                    continue
                fi
                existing_binds+=("${bind_args[$bind_index]}")
            done
            bind_args=("${existing_binds[@]}")
        fi
    else
        runtime="mamba"
    fi
    wrapper_check="$(cat <<'SH_COMEBIN_WRAPPER'
runner=$(command -v "$1") || exit 1
test -x "$runner" || exit 1
help_output=$(bash "$runner" -h 2>&1) || exit 1
grep -Eq "^[[:space:]]*-d[[:space:]]" <<< "$help_output" || exit 1
printf "%s\n" "$runner"
SH_COMEBIN_WRAPPER
)"
    if [[ "${EXEC_MODE}" == apptainer ]]; then
        if ! wrapper_path="$("${runtime}" exec --cleanenv "${bind_args[@]}" \
                "${SIF_COMEBIN}" bash -c "${wrapper_check}" _ "${COMEBIN_SCRIPT}" 2>&1)"; then
            log_error "[COMEBin] run_comebin.sh is missing, not executable, or does not accept -d in ${SIF_COMEBIN}."
            log_error "${wrapper_path}"
            return 1
        fi
    else
        if ! wrapper_path="$(mamba run -n "${ENV_COMEBIN}" \
                bash -c "${wrapper_check}" _ "${COMEBIN_SCRIPT}" 2>&1)"; then
            log_error "[COMEBin] run_comebin.sh is missing, not executable, or does not accept -d in ${ENV_COMEBIN}."
            log_error "${wrapper_path}"
            return 1
        fi
    fi
    if [[ "${EXEC_MODE}" == apptainer ]]; then
        if ! version_output="$("${runtime}" exec --cleanenv "${bind_args[@]}" \
                "${SIF_COMEBIN}" "${wrapper_path}" -V 2>&1)"; then
            log_error "[COMEBin] Could not query the version from ${SIF_COMEBIN}."
            log_error "${version_output}"
            return 1
        fi
    else
        if ! version_output="$(mamba run -n "${ENV_COMEBIN}" "${wrapper_path}" -V 2>&1)"; then
            log_error "[COMEBin] Could not query the version from ${ENV_COMEBIN}."
            log_error "${version_output}"
            return 1
        fi
    fi
    if [[ "${version_output}" =~ COMEBin[[:space:]]version:[[:space:]]*([^[:space:]]+) ]]; then
        COMEBIN_VERSION="${BASH_REMATCH[1]}"
    else
        log_error "[COMEBin] Could not parse the version reported by run_comebin.sh."
        log_error "${version_output}"
        return 1
    fi

    [[ -s "${CONTIGS}" ]] && visible_inputs+=("${CONTIGS}")
    if [[ "${COVERAGE_SOURCE}" == external ]]; then
        visible_inputs+=("${ACTIVE_BAMS[@]}" "${ACTIVE_BAM_INDEXES[@]}")
    elif [[ -d "${MAPPING_DIR}" ]]; then
        while IFS= read -r input; do
            [[ -s "${input}" ]] && visible_inputs+=("${input}")
        done < <(find "${MAPPING_DIR}" -maxdepth 1 -type f -name '*_sorted.bam' -print | sort)
    fi

    python_check="$(cat <<'PY_COMEBIN_RUNTIME'
import os
import sys
import torch
import Bio
import click
import hnswlib
import igraph
import leidenalg
import numpy
import pandas
import sklearn
import tensorboard
import yaml
from torch.utils.tensorboard import SummaryWriter

runner = os.path.realpath(sys.argv[1])
sys.path.insert(0, os.path.dirname(runner))
import cluster
import simclr
import train_CLmodel
import utils

torch.tensor([1], device="cpu")
for path in sys.argv[2:]:
    if not os.path.exists(path):
        raise SystemExit("COMEBin input is not visible in its runtime: " + path)
print("\t".join((
    torch.__version__,
    torch.version.cuda or "none",
    "true" if torch.cuda.is_available() else "false",
    str(torch.cuda.device_count()),
)))
PY_COMEBIN_RUNTIME
)"
    if [[ "${EXEC_MODE}" == apptainer ]]; then
        if ! runtime_info="$("${runtime}" exec --cleanenv "${runtime_args[@]}" \
                "${bind_args[@]}" "${SIF_COMEBIN}" \
                /opt/env/bin/python -c "${python_check}" "${wrapper_path}" "${visible_inputs[@]}" 2>&1)"; then
            log_error "[COMEBin] PyTorch/dependency/input preflight failed in ${SIF_COMEBIN}."
            log_error "${runtime_info}"
            if [[ "${torch_cuda_requested}" == true ]]; then
                log_error "COMEBin device '${COMEBIN_DEVICE}' was requested, but the selected runtime does not provide a usable CUDA-enabled PyTorch environment. Use --comebin-device cpu or provide a CUDA-enabled COMEBin image."
            fi
            return 1
        fi
    else
        if ! runtime_info="$(mamba run -n "${ENV_COMEBIN}" \
                python -c "${python_check}" "${wrapper_path}" "${visible_inputs[@]}" 2>&1)"; then
            log_error "[COMEBin] PyTorch/dependency/input preflight failed in ${ENV_COMEBIN}."
            log_error "${runtime_info}"
            if [[ "${torch_cuda_requested}" == true ]]; then
                log_error "COMEBin device '${COMEBIN_DEVICE}' was requested, but the selected runtime does not provide a usable CUDA-enabled PyTorch environment. Use --comebin-device cpu or provide a CUDA-enabled COMEBin image."
            fi
            return 1
        fi
    fi

    runtime_info="${runtime_info##*$'\n'}"
    IFS=$'\t' read -r COMEBIN_PYTORCH_VERSION COMEBIN_PYTORCH_CUDA_BUILD \
        COMEBIN_PYTORCH_CUDA_AVAILABLE COMEBIN_CUDA_DEVICE_COUNT <<< "${runtime_info}"
    if [[ -z "${COMEBIN_PYTORCH_VERSION}" || \
          ( "${COMEBIN_PYTORCH_CUDA_AVAILABLE}" != true && \
            "${COMEBIN_PYTORCH_CUDA_AVAILABLE}" != false ) || \
          ! "${COMEBIN_CUDA_DEVICE_COUNT}" =~ ^[0-9]+$ ]]; then
        log_error "[COMEBin] Could not parse PyTorch runtime inspection output."
        log_error "${runtime_info}"
        return 1
    fi
    if [[ "${torch_cuda_requested}" == true ]]; then
        if [[ "${COMEBIN_PYTORCH_CUDA_BUILD}" == none || \
             "${COMEBIN_PYTORCH_CUDA_AVAILABLE}" != true || \
             ! "${COMEBIN_CUDA_DEVICE_COUNT}" =~ ^[1-9][0-9]*$ ]]; then
            log_error "COMEBin device '${COMEBIN_DEVICE}' was requested, but the selected runtime does not provide a usable CUDA-enabled PyTorch environment. Use --comebin-device cpu or provide a CUDA-enabled COMEBin image."
            return 1
        fi
        if [[ "${COMEBIN_DEVICE}" == cuda:* ]]; then
            device_index="${COMEBIN_DEVICE#cuda:}"
            if ! python3 -c 'import sys; raise SystemExit(int(sys.argv[1]) >= int(sys.argv[2]))' \
                    "${device_index}" "${COMEBIN_CUDA_DEVICE_COUNT}"; then
                log_error "COMEBin device index ${device_index} is out of range; CUDA device count is ${COMEBIN_CUDA_DEVICE_COUNT}."
                return 1
            fi
        fi
    fi

    log_info "COMEBin selected       : true"
    log_info "COMEBin version        : ${COMEBIN_VERSION}"
    log_info "COMEBin device         : ${COMEBIN_DEVICE}"
    log_info "COMEBin seed           : ${COMEBIN_SEED}"
    log_info "PyTorch version        : ${COMEBIN_PYTORCH_VERSION}"
    log_info "PyTorch CUDA build     : ${COMEBIN_PYTORCH_CUDA_BUILD}"
    log_info "CUDA available         : ${COMEBIN_PYTORCH_CUDA_AVAILABLE}"
    log_info "CUDA device count      : ${COMEBIN_CUDA_DEVICE_COUNT}"
    log_info "Apptainer --nv         : ${COMEBIN_APPTAINER_NV}"
    if [[ -s "${CONTIGS}" ]]; then
        log_info "[COMEBin] Assembly is present and readable in the selected runtime."
    else
        log_info "[COMEBin] Assembly visibility check deferred until upstream assembly completes."
    fi
    if (( ${#visible_inputs[@]} == 1 )) && [[ -s "${CONTIGS}" ]]; then
        log_info "[COMEBin] BAM visibility check deferred until mapping completes."
    elif (( ${#visible_inputs[@]} > 1 )); then
        log_info "[COMEBin] Existing BAM inputs are visible in the selected runtime."
    fi
    [[ "${DRY_RUN}" == true ]] && show_comebin_command
    return 0
}

# ── run_mamba_file(step_label, env_name, script_file) ────────────────────────
# Hybrid dispatcher: routes execution through Apptainer or mamba depending on
# EXEC_MODE. Signature is unchanged so all step functions work without edits.
run_mamba_file() {
    local step_label="$1"
    local env_name="$2"
    local script_file="$3"

    CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[1]:-unknown}"
    CURRENT_COMMAND_SOURCE_LINE="${BASH_LINENO[0]:-0}"
    CURRENT_STAGE_COMMAND="bash $(printf '%q' "${script_file}")"

    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        local sif_path uri
        read -r sif_path uri < <(_env_to_sif "${env_name}")
        local abs_sif
        abs_sif="$(realpath -m "${sif_path}" 2>/dev/null || echo "${sif_path}")"
        CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"

        log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
        log_info "[Step ${step_label}] Script: ${script_file}"

        ensure_sif_image "${sif_path}" "${uri}" || return 1
        local -a bind_args=()
        build_bind_mounts bind_args

        local -a gpu_args=()
        if [[ "${env_name}" == "${ENV_VAMB}" ]] && _has_gpu; then
            gpu_args=(--nv)
            log_info "[Step ${step_label}] GPU detected — adding --nv flag."
        elif [[ "${env_name}" == "${ENV_COMEBIN}" && \
                ( "${COMEBIN_DEVICE}" == cuda || "${COMEBIN_DEVICE}" == cuda:* ) ]]; then
            gpu_args=(--nv)
            log_info "[Step ${step_label}] COMEBin CUDA requested — adding --nv flag."
        fi

        local runtime="apptainer"
        command -v apptainer &>/dev/null || runtime="singularity"

        if [[ "${DRY_RUN}" == true ]]; then
            printf '[DRY-RUN] %q ' "${runtime}" exec --cleanenv "${gpu_args[@]}" "${bind_args[@]}" "${sif_path}" bash "${script_file}"
            printf '\n'
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
            return 0
        fi

        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        "${runtime}" exec --cleanenv "${gpu_args[@]}" "${bind_args[@]}" \
            "${sif_path}" bash "${script_file}"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi

        return "${rc}"

    else
        CURRENT_STAGE_RUNTIME="mamba:${env_name}"
        log_info "[Step ${step_label}] Engine: mamba | Env: ${BOLD}${env_name}${NC}"
        log_info "[Step ${step_label}] Script: ${script_file}"
        if [[ "${DRY_RUN}" == true ]]; then
            echo -e "${YELLOW}[DRY-RUN]${NC} mamba run -n ${env_name} bash ${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
            return 0
        fi

        local had_errexit=0
        local rc=0
        [[ $- == *e* ]] && had_errexit=1

        set +e
        mamba run -n "${env_name}" bash "${script_file}"
        rc=$?
        (( had_errexit )) && set -e

        if (( rc == 0 )); then
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
        fi

        return "${rc}"
    fi
}

check_tool() {
    if ! command -v "$1" &>/dev/null; then
        log_error "Required tool not found in PATH: $1"
        return 1
    fi
}

# ── Checkpoint Helpers ────────────────────────────────────────────────────────
CHECKPOINT_DIR=""  # initialized after OUTDIR is created

require_nonempty_file() {
    local path="$1"
    local description="$2"
    if [[ ! -s "${path}" ]]; then
        log_error "[Validation] ${description} is missing or empty: ${path}"
        return 1
    fi
}

require_directory() {
    local path="$1"
    local description="$2"
    if [[ ! -d "${path}" ]]; then
        log_error "[Validation] ${description} directory is missing: ${path}"
        return 1
    fi
}

collect_files_by_suffix() {
    local directory="$1"
    local suffix="$2"
    local output_name="$3"
    # shellcheck disable=SC2178,SC2034
    local -n result_ref="${output_name}"
    result_ref=()
    [[ -d "${directory}" ]] || return 0
    # shellcheck disable=SC2034
    mapfile -d '' -t result_ref < <(
        find "${directory}" -maxdepth 1 -type f -name "*${suffix}" -print0 | sort -z
    )
}

compute_file_fingerprint() {
    local path="$1"
    python3 - "$path" <<'PY'
import hashlib, os, sys
path = sys.argv[1]
st = os.stat(path)
print(f"{st.st_size}:{st.st_mtime_ns}")
PY
}

validate_sample_id() {
    local sample_id="$1"
    if [[ -z "${sample_id}" ]]; then
        log_error "Sample ID is empty."
        return 1
    fi
    if [[ ! "${sample_id}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]]; then
        log_error "Invalid sample ID '${sample_id}'. It must match ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ and contain no spaces or reserved characters."
        return 1
    fi
    local lower="${sample_id,,}"
    case "${lower}" in
        na|null|none|true|false|unknown)
            log_error "Sample ID '${sample_id}' is reserved and cannot be used."
            return 1
            ;;
    esac
    return 0
}

discover_paired_reads() {
    local source_dir="$1"
    local out_sample_name="$2"
    local out_r1_name="$3"
    local out_r2_name="$4"
    local -n out_sample_ref="${out_sample_name}"
    local -n out_r1_ref="${out_r1_name}"
    local -n out_r2_ref="${out_r2_name}"
    out_sample_ref=()
    out_r1_ref=()
    out_r2_ref=()
    local -a files=()
    local -A seen_exact=()
    local -A seen_case=()
    local file sample_id r1 r2 lower

    [[ -d "${source_dir}" ]] || {
        log_error "Read directory does not exist: ${source_dir}"
        return 1
    }
    mapfile -d '' -t files < <(find "${source_dir}" -maxdepth 1 -type f -name "*${READ_SUFFIX_R1}" -print0 | LC_ALL=C sort -z)

    for file in "${files[@]}"; do
        if [[ -L "${file}" ]]; then
            log_error "Symbolic links for read files are not supported in this version: ${file}"
            return 1
        fi
        sample_id="$(basename "${file}" "${READ_SUFFIX_R1}")"
        if ! validate_sample_id "${sample_id}"; then
            return 1
        fi
        lower="${sample_id,,}"
        if [[ -n "${seen_exact[${sample_id}]:-}" ]]; then
            log_error "Duplicate sample identifier: ${sample_id}"
            return 1
        fi
        if [[ -n "${seen_case[${lower}]:-}" ]]; then
            log_error "Sample IDs must also be unique when compared case-insensitively: '${sample_id}' conflicts with '${seen_case[${lower}]}'."
            return 1
        fi
        seen_exact["${sample_id}"]=1
        seen_case["${lower}"]="${sample_id}"

        r1="${file}"
        r2="${file/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}"
        if [[ "${r1}" == "${r2}" ]]; then
            log_error "Read pair resolution failed for sample ${sample_id}: ${r1}"
            return 1
        fi
        if [[ ! -f "${r1}" || ! -s "${r1}" ]]; then
            log_error "R1 read is missing or empty: ${r1}"
            return 1
        fi
        if [[ ! -f "${r2}" || ! -s "${r2}" ]]; then
            log_error "Matching R2 read is missing or empty: ${r2}"
            return 1
        fi
        if [[ -L "${r2}" ]]; then
            log_error "Symbolic links for read files are not supported in this version: ${r2}"
            return 1
        fi
        out_sample_ref+=("${sample_id}")
        out_r1_ref+=("${r1}")
        out_r2_ref+=("${r2}")
    done

    if (( ${#out_sample_ref[@]} == 0 )); then
        log_error "No reads matching *${READ_SUFFIX_R1} found in ${source_dir}."
        return 1
    fi
}

build_sample_manifest_candidate() {
    local source_dir="$1"
    local reads_source="$2"
    local tmp_path="$3"
    local -a sample_ids=() r1_files=() r2_files=()
    discover_paired_reads "${source_dir}" sample_ids r1_files r2_files || return 1
    local sample_id r1 r2 resolved_r1 resolved_r2 size1 size2 fp1 fp2
    : > "${tmp_path}"
    {
        printf '%s\n' "sample_id	reads_layout	reads_source	read1_path	read1_resolved_path	read2_path	read2_resolved_path	read1_size_bytes	read2_size_bytes	read1_fingerprint	read2_fingerprint	fingerprint_method	input_read1_name	input_read2_name"
        for i in "${!sample_ids[@]}"; do
            sample_id="${sample_ids[$i]}"
            r1="${r1_files[$i]}"
            r2="${r2_files[$i]}"
            resolved_r1="$(realpath "${r1}")"
            resolved_r2="$(realpath "${r2}")"
            [[ -n "${resolved_r1}" && -n "${resolved_r2}" ]] || {
                log_error "Resolved read paths are empty for ${sample_id}."
                return 1
            }
            size1="$(stat -c '%s' "${resolved_r1}")"
            size2="$(stat -c '%s' "${resolved_r2}")"
            fp1="$(compute_file_fingerprint "${resolved_r1}")"
            fp2="$(compute_file_fingerprint "${resolved_r2}")"
            printf '%s\tpaired\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tsize-mtime-ns\t%s\t%s\n' \
                "${sample_id}" "${reads_source}" "${r1}" "${resolved_r1}" "${r2}" "${resolved_r2}" "${size1}" "${size2}" "${fp1}" "${fp2}" "$(basename "${r1}")" "$(basename "${r2}")"
        done
    } >> "${tmp_path}"
}

validate_sample_manifest() {
    local manifest_path="$1"
    local expected_source="$2"
    local header expected_line
    expected_line="sample_id	reads_layout	reads_source	read1_path	read1_resolved_path	read2_path	read2_resolved_path	read1_size_bytes	read2_size_bytes	read1_fingerprint	read2_fingerprint	fingerprint_method	input_read1_name	input_read2_name"
    [[ -s "${manifest_path}" ]] || {
        log_error "Sample manifest is missing or empty: ${manifest_path}"
        return 1
    }
    header="$(head -n 1 "${manifest_path}")"
    [[ "${header}" == "${expected_line}" ]] || {
        log_error "Sample manifest header mismatch in ${manifest_path}."
        return 1
    }
    local row expected_rows=0 sample_id reads_layout reads_source read1_path read1_resolved_path read2_path read2_resolved_path fingerprint_method input_read1_name input_read2_name
    local -A seen_paths=() seen_samples=() seen_case=()
    while IFS=$'\t' read -r sample_id reads_layout reads_source read1_path read1_resolved_path read2_path read2_resolved_path _ _ _ _ fingerprint_method input_read1_name input_read2_name; do
        [[ -n "${sample_id}" ]] || { log_error "Empty sample_id in ${manifest_path}"; return 1; }
        if ! validate_sample_id "${sample_id}"; then
            log_error "Invalid sample_id '${sample_id}' in ${manifest_path}."
            return 1
        fi
        [[ "${reads_layout}" == "paired" ]] || { log_error "reads_layout must be paired in ${manifest_path}"; return 1; }
        [[ "${reads_source}" == "${expected_source}" ]] || { log_error "reads_source mismatch in ${manifest_path}: expected ${expected_source}, found ${reads_source}"; return 1; }
        [[ "${fingerprint_method}" == "size-mtime-ns" ]] || { log_error "Unsupported fingerprint_method in ${manifest_path}: ${fingerprint_method}"; return 1; }
        [[ -n "${input_read1_name}" && -n "${input_read2_name}" ]] || { log_error "input names cannot be empty in ${manifest_path}"; return 1; }
        [[ "${read1_path}" == /* && "${read2_path}" == /* ]] || { log_error "Read paths must be absolute in ${manifest_path}"; return 1; }
        [[ "${read1_path}" != "${read2_path}" ]] || { log_error "R1 and R2 paths must differ for sample ${sample_id}."; return 1; }
        [[ -f "${read1_resolved_path}" && -f "${read2_resolved_path}" ]] || { log_error "Resolved read paths must point to regular files for ${sample_id}."; return 1; }
        [[ -n "${seen_paths[${read1_resolved_path}]:-}" ]] && { log_error "Duplicate read1_resolved_path in ${manifest_path}: ${read1_resolved_path}"; return 1; }
        [[ -n "${seen_paths[${read2_resolved_path}]:-}" ]] && { log_error "Duplicate read2_resolved_path in ${manifest_path}: ${read2_resolved_path}"; return 1; }
        seen_paths["${read1_resolved_path}"]=1
        seen_paths["${read2_resolved_path}"]=1
        [[ -n "${seen_samples[${sample_id}]:-}" ]] && { log_error "Duplicate sample_id in ${manifest_path}: ${sample_id}"; return 1; }
        local sample_lower="${sample_id,,}"
        [[ -n "${seen_case[${sample_lower}]:-}" ]] && { log_error "Sample IDs must also be unique when compared case-insensitively: '${sample_id}' conflicts with '${seen_case[${sample_lower}]}'"; return 1; }
        seen_samples["${sample_id}"]=1
        seen_case["${sample_lower}"]="${sample_id}"
        (( expected_rows += 1 ))
    done < <(tail -n +2 "${manifest_path}")
    (( expected_rows > 0 )) || {
        log_error "Sample manifest contains no rows: ${manifest_path}"
        return 1
    }
    return 0
}

publish_sample_manifest() {
    local candidate_path="$1"
    local target_manifest="$2"
    local target_metadata="$3"
    local reads_source="$4"
    validate_sample_manifest "${candidate_path}" "${reads_source}" || return 1
    local dir
    dir="$(dirname "${target_manifest}")"
    mkdir -p "${dir}"
    local tmp_manifest
    tmp_manifest="$(mktemp "${dir}/.$(basename "${target_manifest}").XXXXXX")"
    cp "${candidate_path}" "${tmp_manifest}"
    mv -f "${tmp_manifest}" "${target_manifest}"
    write_sample_manifest_metadata "${target_manifest}" "${target_metadata}" "${reads_source}"
}

write_sample_manifest_metadata() {
    local manifest_path="$1"
    local metadata_path="$2"
    local reads_source="$3"
    local parent_manifest=""
    local parent_fingerprint=""
    local sample_manifest_fingerprint sample_manifest_location_fingerprint created_at
    created_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    if [[ "${reads_source}" == "fastp-clean" ]]; then
        parent_manifest="${RAW_SAMPLE_MANIFEST}"
        parent_fingerprint="$(python3 - "$parent_manifest" <<'PY'
import hashlib, pathlib, sys
p = pathlib.Path(sys.argv[1])
if not p.exists():
    print('missing')
    raise SystemExit(0)
text = p.read_text(encoding='utf-8')
# logical fingerprint based on manifest rows excluding paths and metadata fields
lines = [ln for ln in text.splitlines() if ln and not ln.startswith('sample_id\t')]
parts = []
for ln in lines:
    fields = ln.split('\t')
    if len(fields) < 14:
        continue
    sample_id, reads_layout, reads_source, _, _, _, _, _, _, read1_fp, read2_fp, fp_method, _, _ = fields[:14]
    parts.append("\t".join([sample_id, reads_layout, reads_source, read1_fp, read2_fp, fp_method]))
print(hashlib.sha256(('schema_version=1\n' + '\n'.join(sorted(parts)) + '\n').encode('utf-8')).hexdigest())
PY
)"
    else
        parent_manifest="null"
        parent_fingerprint="null"
    fi
    sample_manifest_fingerprint="$(python3 - "$manifest_path" <<'PY'
import hashlib, sys
from pathlib import Path
p = Path(sys.argv[1])
lines = p.read_text(encoding='utf-8').splitlines()
if not lines: raise SystemExit('missing')
header = lines[0]
rows = []
for line in lines[1:]:
    if not line.strip():
        continue
    fields = line.split('\t')
    if len(fields) < 14:
        continue
    sample_id, reads_layout, reads_source, _, _, _, _, _, _, read1_fp, read2_fp, fp_method, _, _ = fields[:14]
    rows.append((sample_id, reads_layout, reads_source, read1_fp, read2_fp, fp_method))
body = '\n'.join('|' .join((sample_id, reads_layout, reads_source, read1_fp, read2_fp, fp_method)) for sample_id, reads_layout, reads_source, read1_fp, read2_fp, fp_method in sorted(rows))
print(hashlib.sha256(("schema_version=1\n" + body + "\n").encode('utf-8')).hexdigest())
PY
)"
    sample_manifest_location_fingerprint="$(python3 - "$manifest_path" <<'PY'
import hashlib, sys
from pathlib import Path
p = Path(sys.argv[1])
lines = p.read_text(encoding='utf-8').splitlines()
if not lines:
    raise SystemExit('missing')
result = []
for line in lines[1:]:
    if not line.strip():
        continue
    fields = line.split('\t')
    if len(fields) < 14:
        continue
    sample_id, _, _, read1_path, read1_resolved_path, read2_path, read2_resolved_path, _, _, _, _, _, input_read1_name, input_read2_name = fields[:14]
    if len(fields) >= 14:
        input_read1_name = fields[12]
        input_read2_name = fields[13]
    result.append('|'.join([sample_id, read1_path, read1_resolved_path, read2_path, read2_resolved_path, input_read1_name, input_read2_name]))
print(hashlib.sha256(('location\n' + '\n'.join(sorted(result)) + '\n').encode('utf-8')).hexdigest())
PY
)"
    python3 - "$manifest_path" "$metadata_path" "$reads_source" "$created_at" "$sample_manifest_fingerprint" "$sample_manifest_location_fingerprint" "$parent_manifest" "$parent_fingerprint" <<'PY'
import json, os, sys
manifest_path, metadata_path, reads_source, created_at, sample_manifest_fingerprint, sample_manifest_location_fingerprint, parent_manifest, parent_fingerprint = sys.argv[1:9]
obj = {
    'schema_version': '1',
    'created_at': created_at,
    'reads_layout': 'paired',
    'reads_source': reads_source,
    'sample_count': sum(1 for line in open(manifest_path, encoding='utf-8') if line.strip() and not line.startswith('sample_id\t')),
    'fingerprint_method': 'size-mtime-ns',
    'sample_manifest_fingerprint': sample_manifest_fingerprint,
    'sample_manifest_location_fingerprint': sample_manifest_location_fingerprint,
    'parent_manifest': parent_manifest if parent_manifest != 'null' else None,
    'parent_manifest_fingerprint': parent_fingerprint if parent_fingerprint != 'null' else None,
    'generator': 'run_metagenomics_pipeline.sh',
    'generator_version': os.environ.get('PIPELINE_VERSION', 'unknown')
}
with open(metadata_path, 'w', encoding='utf-8') as fh:
    json.dump(obj, fh, indent=2, sort_keys=True)
    fh.write('\n')
PY
}

load_sample_manifest() {
    local manifest_path="$1"
    local expected_source="$2"
    local sample_name="$3"
    local r1_name="$4"
    local r2_name="$5"
    local -n sample_ref="${sample_name}"
    local -n r1_ref="${r1_name}"
    local -n r2_ref="${r2_name}"
    sample_ref=()
    r1_ref=()
    r2_ref=()
    local line sample_id reads_layout reads_source read1_path read2_path
    if [[ ! -s "${manifest_path}" ]]; then
        log_error "Sample manifest not found or empty: ${manifest_path}"
        return 1
    fi
    validate_sample_manifest "${manifest_path}" "${expected_source}" || return 1
    while IFS=$'\t' read -r sample_id reads_layout reads_source read1_path _ read2_path _ _ _ _ _ _ _ _; do
        [[ "${sample_id}" == "sample_id" ]] && continue
        [[ -n "${sample_id}" ]] || continue
        if ! validate_sample_id "${sample_id}"; then
            log_error "Manifest contains an invalid sample_id '${sample_id}' in ${manifest_path}."
            return 1
        fi
        sample_ref+=("${sample_id}")
        r1_ref+=("${read1_path}")
        r2_ref+=("${read2_path}")
    done < "${manifest_path}"
    (( ${#sample_ref[@]} > 0 )) || {
        log_error "Sample manifest loaded no samples: ${manifest_path}"
        return 1
    }
}

collect_files_by_pattern() {
    local directory="$1"
    local pattern="$2"
    local output_name="$3"
    # shellcheck disable=SC2178,SC2034
    local -n result_ref="${output_name}"
    result_ref=()
    [[ -d "${directory}" ]] || return 0
    # shellcheck disable=SC2034
    mapfile -d '' -t result_ref < <(
        find "${directory}" -maxdepth 1 -type f -name "${pattern}" -print0 | sort -z
    )
}

require_files_matching() {
    local directory="$1"
    local pattern="$2"
    local description="$3"
    local -a matches=()
    local file

    if [[ ! -d "${directory}" ]]; then
        log_error "[Validation] ${description} directory is missing: ${directory}"
        return 1
    fi
    collect_files_by_pattern "${directory}" "${pattern}" matches
    for file in "${matches[@]}"; do
        [[ -s "${file}" ]] && return 0
    done
    log_error "[Validation] No non-empty ${description} matching '${pattern}' in ${directory}"
    return 1
}

sample_names() {
    printf '%s\n' "${ACTIVE_SAMPLE_IDS[@]}"
}

expected_sample_names() {
    sample_names
}

stage_rank() {
    case "$1" in
        mapping) echo 1 ;;
        binning) echo 2 ;;
        refinement) echo 3 ;;
        reassembly) echo 4 ;;
        checkm2) echo 5 ;;
        mag-filtering) echo 6 ;;
        gtdbtk) echo 7 ;;
        coverm) echo 8 ;;
        taxonomy-summary) echo 9 ;;
        mag-summary) echo 10 ;;
        run-overview|summary) echo 11 ;;
        *) return 1 ;;
    esac
}

stage_in_scope() {
    local stage="$1"
    [[ -z "${STOP_AFTER}" ]] && return 0
    [[ "$(stage_rank "${stage}")" -le "$(stage_rank "${STOP_AFTER}")" ]]
}

stop_after_reached() {
    [[ -n "${STOP_AFTER}" && ( "$1" == "${STOP_AFTER}" || ( "$1" == run-overview && "${STOP_AFTER}" == summary ) ) ]]
}

validate_stage_dependencies() {
    local errors=0

    if [[ -n "${STOP_AFTER}" ]]; then
        if ! stage_rank "${STOP_AFTER}" >/dev/null; then
            log_error "Invalid --stop-after stage: ${STOP_AFTER}"
            (( errors += 1 ))
        elif ! stage_in_scope "${STOP_AFTER}"; then
            log_error "Invalid --stop-after stage: ${STOP_AFTER}"
            (( errors += 1 ))
        else
            case "${STOP_AFTER}" in
                mapping) [[ "${RUN_MAPPING}" == true ]] || { log_error "--stop-after mapping requires --mapping."; (( errors += 1 )); } ;;
                binning) [[ "${RUN_BINNING}" == true ]] || { log_error "--stop-after binning requires --binning."; (( errors += 1 )); } ;;
                refinement) [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] || { log_error "--stop-after refinement requires effective consensus refinement."; (( errors += 1 )); } ;;
                reassembly) [[ "${RUN_REASSEMBLY}" == true ]] || { log_error "--stop-after reassembly requires --reassembly."; (( errors += 1 )); } ;;
                checkm2) [[ "${RUN_CHECKM2}" == true ]] || { log_error "--stop-after checkm2 requires --checkm2."; (( errors += 1 )); } ;;
                mag-filtering) [[ "${RUN_MAG_FILTERING}" == true ]] || { log_error "--stop-after mag-filtering requires --mag-filtering."; (( errors += 1 )); } ;;
                gtdbtk) [[ "${RUN_GTDBTK}" == true ]] || { log_error "--stop-after gtdbtk requires --gtdbtk."; (( errors += 1 )); } ;;
                coverm) [[ "${RUN_COVERM}" == true ]] || { log_error "--stop-after coverm requires --coverm."; (( errors += 1 )); } ;;
                taxonomy-summary) [[ "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true ]] || { log_error "--stop-after taxonomy-summary requires --summary."; (( errors += 1 )); } ;;
                mag-summary) [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true ]] || { log_error "--stop-after mag-summary requires --summary."; (( errors += 1 )); } ;;
                run-overview|summary) [[ "${RUN_SUMMARY}" == true ]] || { log_error "--stop-after ${STOP_AFTER} requires --summary."; (( errors += 1 )); } ;;
            esac
        fi
    fi

    if [[ "${RUN_BINNING}" != true ]]; then
        if [[ "${RUN_CHECKM2}" == true || "${RUN_MAG_FILTERING}" == true || \
              "${RUN_GTDBTK}" == true || "${RUN_COVERM}" == true ]]; then
            log_error "CheckM2, MAG filtering, GTDB-Tk, and CoverM require RUN_BINNING=true; disable dependent stages explicitly."
            (( errors += 1 ))
        fi
        if [[ "${RUN_REASSEMBLY}" == true ]]; then
            log_error "Reassembly requires bins from an active binning/refinement path."
            (( errors += 1 ))
        fi
    fi
    if [[ "${RUN_REASSEMBLY}" == true && -z "${BINS_FOR_REASSEMBLY}" ]]; then
        log_error "Reassembly requires a resolved selected bin collection."
        (( errors += 1 ))
    fi
    if [[ "${RUN_CHECKM2}" == true && "${RUN_BINNING}" == true && -z "${BINS_FOR_QC}" ]]; then
        log_error "CheckM2 requires a resolved selected bin collection."
        (( errors += 1 ))
    fi
    if [[ "${RUN_MAG_FILTERING}" == true && "${RUN_CHECKM2}" != true ]]; then
        log_error "MAG filtering requires CheckM2; enable --checkm2 or disable --mag-filtering."
        (( errors += 1 ))
    fi
    if [[ "${RUN_GTDBTK}" == true && "${RUN_MAG_FILTERING}" != true ]]; then
        log_error "GTDB-Tk requires MAG filtering; enable --mag-filtering or disable --gtdbtk."
        (( errors += 1 ))
    fi
    if [[ "${RUN_COVERM}" == true && "${RUN_MAG_FILTERING}" != true ]]; then
        log_error "CoverM requires MAG filtering; enable --mag-filtering or disable --coverm."
        (( errors += 1 ))
    fi
    if [[ "${RUN_SUMMARY}" == true && \
          ( "${RUN_MAG_FILTERING}" != true || \
            ( "${RUN_CHECKM2}" != true && "${RUN_GTDBTK}" != true && "${RUN_COVERM}" != true ) ) ]]; then
        log_error "Summary requires selected MAGs and at least one active tabular source."
        (( errors += 1 ))
    fi

    if [[ "${errors}" -gt 0 ]]; then
        log_error "Stage configuration is incompatible (${errors} issue(s))."
        return 1
    fi
}

append_reads_manifest() {
    local step="$1" i file size mtime
    local -a r1_files=("${ACTIVE_R1_FILES[@]}")
    local -a r2_files=("${ACTIVE_R2_FILES[@]}")

    if [[ "${step}" == 01_fastp ]]; then
        r1_files=()
        r2_files=()
        local raw_r1
        collect_files_by_suffix "${RAW_READS_DIR}" "${READ_SUFFIX_R1}" r1_files
        for raw_r1 in "${r1_files[@]}"; do
            r2_files+=("${raw_r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
    fi
    printf 'active_reads_source=%s\nactive_reads_dir=%s\nr1_suffix=%s\nr2_suffix=%s\n' \
        "${ACTIVE_READS_SOURCE}" "$(realpath -m "${ACTIVE_READS_DIR}")" \
        "${READ_SUFFIX_R1}" "${READ_SUFFIX_R2}"
    printf 'sample_ids=%s\n' "$(IFS=,; echo "${ACTIVE_SAMPLE_IDS[*]}")"
    for i in "${!r1_files[@]}"; do
        for file in "${r1_files[$i]}" "${r2_files[$i]}"; do
            file="$(realpath -m "${file}")"
            size="$(stat -c '%s' "${file}" 2>/dev/null || printf missing)"
            mtime="$(stat -c '%y' "${file}" 2>/dev/null || printf missing)"
            printf 'read_file=%s|%s|%s\n' "${file}" "${size}" "${mtime}"
        done
    done
}

append_assembly_manifest() {
    local size mtime
    [[ "${RUN_MEGAHIT}" == false && "${NEEDS_ASSEMBLY}" == true && \
       -n "${ASSEMBLY_IN}" ]] || return 0
    size="$(stat -c '%s' "${ASSEMBLY_IN}")"
    mtime="$(stat -c '%y' "${ASSEMBLY_IN}")"
    printf 'external_assembly=%s|%s|%s\n' "${ASSEMBLY_IN}" "${size}" "${mtime}"
}

append_file_fingerprint() {
    local label="$1" path="$2" size mtime
    if [[ -f "${path}" ]]; then
        path="$(realpath "${path}")"
        size="$(stat -c '%s' "${path}")"
        mtime="$(stat -c '%y' "${path}")"
        printf '%s=%s|%s|%s\n' "${label}" "${path}" "${size}" "${mtime}"
    elif [[ -d "${path}" ]]; then
        path="$(realpath "${path}")"
        mtime="$(stat -c '%y' "${path}")"
        printf '%s=%s|directory|%s\n' "${label}" "${path}" "${mtime}"
    else
        printf '%s=missing\n' "${label}"
    fi
}

collect_binner_output_files() {
    local method="$1" output_name="$2" source_dir
    local -n output_files="${output_name}"
    output_files=()
    case "${method}" in
        vamb) source_dir="${VAMB_DIR}/vamb_out/bins" ;;
        metabinner) source_dir="${METABINNER_DIR}/metabinner_out/result_bins" ;;
        comebin) source_dir="${COMEBIN_DIR}/comebin_res/comebin_res_bins" ;;
        *) log_error "Unknown binner: ${method}"; return 1 ;;
    esac
    [[ -d "${source_dir}" ]] || return 0
    case "${method}" in
        vamb)
            # shellcheck disable=SC2034
            mapfile -d '' -t output_files < <(
                find "${source_dir}" -mindepth 1 -maxdepth 1 -type f \
                    \( -name '*.fa' -o -name '*.fna' \) -print0 | sort -z
            ) ;;
        metabinner|comebin)
            # shellcheck disable=SC2034
            mapfile -d '' -t output_files < <(
                find "${source_dir}" -mindepth 1 -maxdepth 1 -type f \
                    \( -name '*.fa' -o -name '*.fasta' \) -print0 | sort -z
            ) ;;
    esac
}

validate_bins_dir() {
    local bins_dir="$1" label="${2:-bins}" entry count=0
    local -A seen_names=()
    local -a entries=()
    if [[ ! -d "${bins_dir}" ]]; then
        log_error "[Validation] ${label} directory not found: ${bins_dir}"
        return 1
    fi
    mapfile -d '' -t entries < <(find "${bins_dir}" -mindepth 1 -maxdepth 1 -print0 | sort -z)
    for entry in "${entries[@]}"; do
        if [[ "$(basename "${entry}")" != *.fa ]]; then
            log_error "[Validation] Unexpected non-bin entry in ${label}: ${entry}"
            return 1
        fi
        if [[ ! -f "${entry}" || ! -s "${entry}" ]]; then
            log_error "[Validation] Missing, broken, or empty bin in ${label}: ${entry}"
            return 1
        fi
        if ! awk '/^>/ { header=1; next } header && /[^[:space:]]/ { sequence=1 } END { exit !(header && sequence) }' "${entry}"; then
            log_error "[Validation] Invalid FASTA bin in ${label}: ${entry}"
            return 1
        fi
        local bin_name
        bin_name="$(basename "${entry}")"
        if [[ -n "${seen_names[$bin_name]:-}" ]]; then
            log_error "[Validation] Duplicate bin name in ${label}: ${bin_name}"
            return 1
        fi
        seen_names["${bin_name}"]=1
        count=$((count + 1))
    done
    if (( count == 0 )); then
        log_error "[Validation] No FASTA bins found in ${label}: ${bins_dir}"
        return 1
    fi
}

normalize_binner_bins() {
    local method="$1" source_dir="$2" destination_dir original stem name relative
    local -a files=()
    local -A seen_names=()
    destination_dir="$(get_binner_bins_dir "${method}")" || return 1
    collect_binner_output_files "${method}" files
    if (( ${#files[@]} == 0 )); then
        log_error "[Binning: ${method}] No final bin FASTA files found in ${source_dir}."
        return 1
    fi
    mkdir -p "${destination_dir}"
    find "${destination_dir}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -delete
    for original in "${files[@]}"; do
        stem="$(basename "${original}")"
        stem="${stem%.*}"
        name="${stem}.fa"
        if [[ -n "${seen_names[$name]:-}" ]]; then
            log_error "[Binning: ${method}] Multiple source files normalize to ${name}."
            return 1
        fi
        seen_names["${name}"]=1
        relative="$(realpath --relative-to="${destination_dir}" "${original}")"
        ln -s -- "${relative}" "${destination_dir}/${name}"
    done
    validate_bins_dir "${destination_dir}" "${method} normalized bins"
}

append_bins_manifest() {
    local bins_dir="$1" bin origin size mtime
    local -a bins=()
    mapfile -d '' -t bins < <(find "${bins_dir}" -mindepth 1 -maxdepth 1 \
        \( -type f -o -type l \) -print0 2>/dev/null | sort -z)
    for bin in "${bins[@]}"; do
        origin="$(realpath "${bin}" 2>/dev/null || printf broken)"
        if [[ -f "${origin}" ]]; then
            size="$(stat -c '%s' "${origin}")"
            mtime="$(stat -c '%y' "${origin}")"
        else
            size=missing
            mtime=missing
        fi
        printf 'bin=%s|%s|%s|%s\n' "$(realpath -m "${bin}")" "${origin}" "${size}" "${mtime}"
    done
}

append_fasta_bins_manifest() {
    local label="$1" bins_dir="$2" bin target size mtime
    local -a bins=()
    mapfile -d '' -t bins < <(find "${bins_dir}" -mindepth 1 -maxdepth 1 \
        \( -type f -o -type l \) -name '*.fa' -print0 2>/dev/null | sort -z)
    for bin in "${bins[@]}"; do
        target="$(realpath "${bin}" 2>/dev/null || printf broken)"
        if [[ -f "${target}" ]]; then
            size="$(stat -c '%s' "${target}")"
            mtime="$(stat -c '%y' "${target}")"
        else
            size=missing
            mtime=missing
        fi
        printf '%s=%s|%s|%s|%s\n' "${label}" "$(realpath -m "${bin}")" \
            "${target}" "${size}" "${mtime}"
    done
}

append_vamb_abundance_manifest() {
    local header
    if [[ -s "${VAMB_ABUNDANCE_TSV}" ]]; then
        IFS= read -r header < "${VAMB_ABUNDANCE_TSV}"
        printf 'abundance_header=%q\n' "${header}"
        python3 - "${VAMB_ABUNDANCE_TSV}" <<'PY_VAMB_MANIFEST'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as stream:
    columns = stream.readline().rstrip("\r\n").split("\t")
print("abundance_samples_json=" + json.dumps(columns[1:], ensure_ascii=True, separators=(",", ":")))
PY_VAMB_MANIFEST
    else
        printf 'abundance_header=missing\n'
        printf 'abundance_samples_json=[]\n'
    fi
}

tsv_value() {
    local value="${1:-}"
    value="${value//$'\t'/_}"
    value="${value//$'\r'/_}"
    value="${value//$'\n'/_}"
    printf '%s' "${value}"
}

inventory_stage_tools() {
    case "$1" in
        01_fastp) printf '%s\n' fastp ;;
        02_megahit) printf '%s\n' megahit ;;
        03_kraken2) printf '%s\n' kraken2 bracken ;;
        04_mapping) printf '%s\n' bowtie2 samtools ;;
        05a_vamb) printf '%s\n' vamb ;;
        05b_metabinner) printf '%s\n' metabinner ;;
        05c_comebin) printf '%s\n' comebin ;;
        06_refinement|07_reassembly) printf '%s\n' metawrap ;;
        08_checkm2) printf '%s\n' checkm2 ;;
        10_gtdbtk) printf '%s\n' gtdbtk ;;
        11_coverm) printf '%s\n' coverm ;;
    esac
}

tool_runtime_spec() {
    local stage="$1" tool="$2"
    case "${stage}:${tool}" in
        01_fastp:fastp) printf '%s\t%s\t%s\t%s\n' "${SIF_FASTP}" "${SIF_FASTP_URI}" fastp "" ;;
        02_megahit:megahit) printf '%s\t%s\t%s\t%s\n' "${SIF_MEGAHIT}" "${SIF_MEGAHIT_URI}" megahit "" ;;
        03_kraken2:kraken2|03_kraken2:bracken) printf '%s\t%s\t%s\t%s\n' "${SIF_KRAKEN2}" "${SIF_KRAKEN2_URI}" "${tool}" "" ;;
        04_mapping:bowtie2|04_mapping:samtools) printf '%s\t%s\t%s\t%s\n' "${SIF_MAPPING}" "${SIF_MAPPING_URI}" "${tool}" "" ;;
        05a_vamb:vamb) printf '%s\t%s\t%s\t%s\n' "${SIF_VAMB}" "${SIF_VAMB_URI}" vamb "${ENV_VAMB}" ;;
        05b_metabinner:metabinner) printf '%s\t%s\t%s\t%s\n' "${SIF_METABINNER}" "${SIF_METABINNER_URI}" metabinner "${ENV_METABINNER}" ;;
        05c_comebin:comebin) printf '%s\t%s\t%s\t%s\n' "${SIF_COMEBIN}" "${SIF_COMEBIN_URI}" run_comebin.sh "${ENV_COMEBIN}" ;;
        06_refinement:metawrap|07_reassembly:metawrap) printf '%s\t%s\t%s\t%s\n' "${SIF_METAWRAP}" "${SIF_METAWRAP_URI}" metawrap "${ENV_METAWRAP}" ;;
        08_checkm2:checkm2) printf '%s\t%s\t%s\t%s\n' "${SIF_CHECKM2}" "${SIF_CHECKM2_URI}" checkm2 "${ENV_CHECKM2}" ;;
        10_gtdbtk:gtdbtk) printf '%s\t%s\t%s\t%s\n' "${SIF_GTDBTK}" "${SIF_GTDBTK_URI}" gtdbtk "${ENV_GTDBTK}" ;;
        11_coverm:coverm) printf '%s\t%s\t%s\t%s\n' "${SIF_COVERM}" "${SIF_COVERM_URI}" coverm "" ;;
    esac
}

tool_version_query() {
    case "$1" in
        fastp) printf '%s\n' '--version' ;;
        megahit|kraken2|bowtie2|samtools|vamb|checkm2|gtdbtk|coverm|Rscript) printf '%s\n' '--version' ;;
        bracken) printf '%s\n' '-v' ;;
        metabinner) printf '%s\n' '--version' ;;
        run_comebin.sh) printf '%s\n' '-V' ;;
        metawrap) printf '%s\n' '--version' ;;
        *) printf '%s\n' '--version' ;;
    esac
}

extract_reported_version() {
    python3 -c '
import re, sys
text = sys.stdin.read()
version = re.search(r"(?<![A-Za-z0-9])v?([0-9]+(?:\.[0-9A-Za-z][0-9A-Za-z+_-]*)+)", text)
print(version.group(1) if version else "unknown")
'
}

mamba_environment_snapshot() {
    python3 - "$1" <<'PY'
import os, sys
root = os.path.join(sys.argv[1], "conda-meta")
rows = []
try:
    for name in os.listdir(root):
        if not name.endswith(".json"):
            continue
        try:
            s = os.stat(os.path.join(root, name))
            rows.append((name.encode().hex(), s.st_size, s.st_mtime_ns, s.st_ctime_ns, s.st_dev, s.st_ino))
        except OSError:
            continue
except OSError:
    pass
for row in sorted(rows):
    print("\t".join(map(str, row)))
PY
}

runtime_identity_for_tool() {
    local stage="$1" tool="$2" spec sif uri executable env_name query command_used
    local mode="${EXEC_MODE}" runtime_name runtime_uri="not-applicable" runtime_path="unknown"
    local runtime_size="not-applicable" runtime_mtime="not-applicable"
    local method="unknown" fingerprint="unknown" version="unknown" version_source="unknown"
    local status=unknown output="" prefix="" key="" cache_value="" runtime_bin
    local command_presence="command -v \"\$1\" >/dev/null 2>&1"
    spec="$(tool_runtime_spec "${stage}" "${tool}")"
    IFS=$'\t' read -r sif uri executable env_name <<< "${spec}"
    query="$(tool_version_query "${executable}")"
    command_used="${executable} ${query}"
    runtime_name="${env_name:-${CONDA_DEFAULT_ENV:-host}}"

    if [[ "${mode}" == apptainer ]]; then
        runtime_name="${tool}"
        runtime_uri="${uri:-unknown}"
        if [[ ! -f "${sif}" || ! -r "${sif}" ]]; then
            status=unavailable
        elif ! runtime_bin="$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null)"; then
            status=unavailable
        else
            runtime_path="$(realpath "${sif}")"
            local sif_stat
            sif_stat="$(python3 -c 'import os,sys; s=os.stat(sys.argv[1]); print(f"{s.st_size}\t{s.st_mtime_ns}\t{s.st_ctime_ns}\t{s.st_dev}\t{s.st_ino}")' "${runtime_path}")"
            IFS=$'\t' read -r runtime_size runtime_mtime <<< "${sif_stat}"
            SIF_OBJECT_CACHE["${runtime_path}"]="$(cut -f3- <<< "${sif_stat}")"
            if [[ -z "${SIF_IDENTITY_CACHE[$runtime_path]:-}" ]]; then
                SIF_IDENTITY_CACHE["${runtime_path}"]="$(sha256sum "${runtime_path}" | awk '{print tolower($1)}')"
            fi
            method=sha256-file
            fingerprint="${SIF_IDENTITY_CACHE[$runtime_path]}"
            key="${mode}|${fingerprint}|${tool}|${command_used}"
            local is_present=false
            if [[ "${tool}" == metabinner || "${tool}" == comebin ]]; then
                if "${runtime_bin}" exec --cleanenv "${runtime_path}" bash -c \
                     "command -v \"${executable}\" >/dev/null 2>&1 || python3 -c \"import ${tool}\" >/dev/null 2>&1"; then
                    is_present=true
                fi
            elif "${runtime_bin}" exec --cleanenv "${runtime_path}" bash -c \
                     "${command_presence}" bash "${executable}"; then
                is_present=true
            fi

            if [[ "${is_present}" != true ]]; then
                status=unavailable
            elif [[ -n "${VERSION_CACHE[$key]:-}" ]]; then
                cache_value="${VERSION_CACHE[$key]}"
            else
                output="$("${runtime_bin}" exec --cleanenv "${runtime_path}" "${executable}" "${query}" 2>&1 || true)"
                version="$(printf '%s' "${output}" | extract_reported_version)"
                version_source="runtime-command"
                if [[ "${version}" == unknown && ( "${tool}" == metabinner || "${tool}" == comebin ) ]]; then
                    local package_names
                    package_names="$([[ "${tool}" == metabinner ]] && printf 'MetaBinner,metabinner' || printf 'COMEBin,comebin')"
                    local const_output
                    const_output="$("${runtime_bin}" exec --cleanenv "${runtime_path}" python3 -c '
import sys
for n in sys.argv[1].split(","):
    try:
        mod = __import__(n)
        v = getattr(mod, "__version__", None) or getattr(mod, "VERSION", None)
        if v:
            print(str(v))
            sys.exit(0)
    except Exception:
        pass
    try:
        import importlib.metadata as m
        v = m.version(n)
        if v:
            print(str(v))
            sys.exit(0)
    except Exception:
        pass
print("unknown")
' "${package_names}" 2>/dev/null || true)"
                    version="$(printf '%s' "${const_output}" | extract_reported_version)"
                    [[ "${version}" != unknown ]] && command_used="python3 -c import ${tool}"
                fi
                if [[ "${version}" == unknown ]]; then
                    version_source="unknown"
                    local label_json
                    label_json="$("${runtime_bin}" inspect --json "${runtime_path}" 2>/dev/null || true)"
                    local label_identity
                    label_identity="$(printf '%s' "${label_json}" | python3 -c '
import json,sys
try:
 d=json.load(sys.stdin)
 labels=d.get("data",{}).get("attributes",{}).get("labels",{})
 v=next((str(labels[k]) for k in ("org.opencontainers.image.version","VERSION","version") if labels.get(k)), "")
 c=next((str(labels[k]) for k in ("org.opencontainers.image.revision","org.opencontainers.image.vcs-ref","VCS_REF","COMMIT","GIT_COMMIT") if labels.get(k)), "")
 print(v + "\timage-label" if v else (c + "\tcommit" if c else "unknown\tunknown"))
except Exception: print("unknown\tunknown")
')"
                    IFS=$'\t' read -r version version_source <<< "${label_identity}"
                    if [[ "${version}" != unknown ]]; then
                        command_used="${runtime_bin} inspect --json"
                    elif [[ "${tool}" == metabinner || "${tool}" == comebin ]]; then
                        local commit_val
                        commit_val="$("${runtime_bin}" exec --cleanenv "${runtime_path}" sh -c 'cat /.commit /commit 2>/dev/null | head -n1' 2>/dev/null || true)"
                        commit_val="$(tr -d " \r\n\t" <<< "${commit_val}")"
                        if [[ -n "${commit_val}" && "${commit_val}" =~ ^[0-9a-fA-F]{7,40}$ ]]; then
                            version="${commit_val}"
                            version_source="commit"
                            command_used="cat /.commit"
                        fi
                    fi
                fi
                [[ "${version}" == unknown ]] && version_source="unknown"
                cache_value="${version}|${version_source}"
                VERSION_CACHE["$key"]="${cache_value}"
            fi
            if [[ "${status}" != unavailable ]]; then
                IFS='|' read -r version version_source <<< "${cache_value}"
                status=available
                [[ "${version}" == unknown ]] && status=unknown
            fi
        fi
    else
        if [[ -n "${env_name}" ]]; then
            runtime_name="${env_name}"
            if ! command -v mamba >/dev/null 2>&1; then
                status=unavailable
            else
                prefix="$(mamba run -n "${env_name}" python -c 'import sys; print(sys.prefix)' 2>/dev/null || true)"
                [[ -n "${prefix}" && -d "${prefix}" ]] || status=unavailable
            fi
        else
            prefix="${CONDA_PREFIX:-}"
            command -v "${executable}" >/dev/null 2>&1 || status=unavailable
        fi
        if [[ "${status}" != unavailable ]]; then
            [[ -z "${prefix}" ]] || runtime_path="$(realpath "${prefix}")"
            if [[ -n "${prefix}" && -d "${prefix}/conda-meta" ]]; then
                if [[ -z "${MAMBA_IDENTITY_CACHE[$runtime_path]:-}" ]]; then
                    MAMBA_IDENTITY_CACHE["${runtime_path}"]="$(python3 - "${runtime_path}" <<'PY'
import hashlib, json, os, sys
root = os.path.join(sys.argv[1], "conda-meta")
rows = []
invalid = False
for name in os.listdir(root):
    if not name.endswith(".json"):
        continue
    try:
        with open(os.path.join(root, name), encoding="utf-8") as stream:
            p = json.load(stream)
        channel = str(p.get("channel") or "")
        if channel and "://" in channel:
            from urllib.parse import urlsplit
            u = urlsplit(channel)
            channel = (u.hostname or "") + u.path.rstrip("/")
        if not channel and p.get("channel_url"):
            from urllib.parse import urlsplit
            u = urlsplit(str(p["channel_url"]))
            channel = (u.hostname or "") + (u.path.rstrip("/").rsplit("/", 1)[-1] if u.path else "")
        channel = channel.split("?", 1)[0].split("#", 1)[0]
        rows.append("\t".join(str(p.get(k) or "") for k in ("name", "version", "build")) + "\t" + channel)
    except (OSError, ValueError, TypeError):
        invalid = True
        continue
payload = ("\n".join(sorted(rows, key=lambda x: x.encode("utf-8"))) + "\n").encode()
print(hashlib.sha256(payload).hexdigest() if rows and not invalid else "unknown")
PY
)"
                fi
                fingerprint="${MAMBA_IDENTITY_CACHE[$runtime_path]}"
                method=conda-meta-v1
                [[ "${fingerprint}" != unknown ]] || method=unknown
                if [[ -z "${MAMBA_OBJECT_CACHE[$runtime_path]:-}" ]]; then
                    MAMBA_OBJECT_CACHE["${runtime_path}"]="$(mamba_environment_snapshot "${runtime_path}")"
                fi
            fi
            if [[ -n "${prefix}" || -z "${env_name}" ]]; then
                key="${mode}|${fingerprint}|${tool}|${command_used}"
                local is_present=false
                if [[ -n "${env_name}" ]]; then
                    if [[ "${tool}" == metabinner || "${tool}" == comebin ]]; then
                        if mamba run -n "${env_name}" bash -c \
                             "command -v \"${executable}\" >/dev/null 2>&1 || python -c \"import ${tool}\" >/dev/null 2>&1"; then
                            is_present=true
                        fi
                    elif mamba run -n "${env_name}" bash -c "${command_presence}" bash "${executable}"; then
                        is_present=true
                    fi
                else
                    if command -v "${executable}" >/dev/null 2>&1; then
                        is_present=true
                    fi
                fi

                if [[ "${is_present}" != true ]]; then
                    status=unavailable
                elif [[ -n "${VERSION_CACHE[$key]:-}" ]]; then
                    cache_value="${VERSION_CACHE[$key]}"
                else
                    if [[ -n "${env_name}" ]]; then
                        output="$(mamba run -n "${env_name}" "${executable}" "${query}" 2>&1 || true)"
                    else
                        output="$("${executable}" "${query}" 2>&1 || true)"
                    fi
                    version="$(printf '%s' "${output}" | extract_reported_version)"
                    version_source="runtime-command"
                    if [[ "${version}" == unknown && ( "${tool}" == metabinner || "${tool}" == comebin ) ]]; then
                        local package_names
                        package_names="$([[ "${tool}" == metabinner ]] && printf 'MetaBinner,metabinner' || printf 'COMEBin,comebin')"
                        local py_query='
import sys
for n in sys.argv[1].split(","):
    try:
        mod = __import__(n)
        v = getattr(mod, "__version__", None) or getattr(mod, "VERSION", None)
        if v:
            print(str(v))
            sys.exit(0)
    except Exception:
        pass
    try:
        import importlib.metadata as m
        v = m.version(n)
        if v:
            print(str(v))
            sys.exit(0)
    except Exception:
        pass
print("unknown")
'
                        if [[ -n "${env_name}" ]]; then
                            output="$(mamba run -n "${env_name}" python -c "${py_query}" "${package_names}" 2>/dev/null || true)"
                        else
                            output="$(python3 -c "${py_query}" "${package_names}" 2>/dev/null || true)"
                        fi
                        version="$(printf '%s' "${output}" | extract_reported_version)"
                        [[ "${version}" != unknown ]] && command_used="python -c import ${tool}"
                    fi
                    if [[ "${version}" == unknown && ( "${tool}" == metabinner || "${tool}" == comebin ) && -n "${prefix}" ]]; then
                        local commit_val
                        commit_val="$(python3 - "${prefix}" "${tool}" <<'PY_COMMIT'
import glob, json, os, sys
prefix, tool = sys.argv[1], sys.argv[2].lower()
found = ""
for f in glob.glob(os.path.join(prefix, "conda-meta", f"*{tool}*.json")):
    try:
        with open(f, encoding="utf-8") as stream:
            d = json.load(stream)
            c = d.get("commit") or d.get("git_url") or ""
            if c:
                found = str(c).split("#")[-1].split("/")[-1]
                break
    except Exception:
        pass
print(found if len(found) in (7, 40) else "")
PY_COMMIT
)"
                        if [[ -n "${commit_val}" ]]; then
                            version="${commit_val}"
                            version_source="commit"
                            command_used="conda-meta git revision"
                        fi
                    fi
                    [[ "${version}" == unknown ]] && version_source="unknown"
                    cache_value="${version}|${version_source}"
                    VERSION_CACHE["$key"]="${cache_value}"
                fi
                if [[ "${status}" != unavailable ]]; then
                    IFS='|' read -r version version_source <<< "${cache_value}"
                    status=available
                    [[ "${version}" == unknown ]] && status=unknown
                fi
            fi
            [[ "${runtime_path}" == unknown ]] && runtime_path=unknown
        fi
    fi

    if [[ "${status}" == unavailable ]]; then
        version=unavailable
        version_source="unknown"
    fi
    local identity_key="${stage}|${tool}"
    ID_VERSION["${identity_key}"]="$(tsv_value "${version}")"
    ID_VERSION_SOURCE["${identity_key}"]="${version_source}"
    ID_EXECUTION_MODE["${identity_key}"]="${mode}"
    ID_RUNTIME_NAME["${identity_key}"]="$(tsv_value "${runtime_name}")"
    ID_RUNTIME_URI["${identity_key}"]="$(tsv_value "${runtime_uri}")"
    ID_RUNTIME_PATH["${identity_key}"]="$(tsv_value "${runtime_path}")"
    ID_RUNTIME_SIZE["${identity_key}"]="${runtime_size}"
    ID_RUNTIME_MTIME["${identity_key}"]="${runtime_mtime}"
    ID_RUNTIME_METHOD["${identity_key}"]="${method}"
    ID_RUNTIME_FINGERPRINT["${identity_key}"]="${fingerprint}"
    ID_COMMAND["${identity_key}"]="$(tsv_value "${command_used}")"
    ID_STATUS["${identity_key}"]="${status}"
}

database_keyfile_rows() {
    local kind="$1" root="$2" output warning
    output="$(python3 - "${kind}" "${root}" <<'PY'
import os, sys
kind, root = sys.argv[1:]
candidates = {
 "kraken2-bracken": [
  "hash.k2d", "opts.k2d", "taxo.k2d", "seqid2taxid.map",
  "taxonomy/nodes.dmp", "taxonomy/names.dmp", "taxonomy/merged.dmp",
  "taxonomy/delnodes.dmp", "database.kraken",
 ]
}
paths = set(candidates.get(kind, []))
if kind == "kraken2-bracken":
    try:
        paths.update(n for n in os.listdir(root) if n.endswith(".kmer_distrib") or n.endswith(".kraken"))
    except OSError:
        pass
elif kind == "gtdbtk":
    paths.update([
      "VERSION", "release.txt", "gtdbtk.package_version",
      "metadata/metadata.tsv", "metadata/metadata.txt",
      "metadata/genome_metadata.tsv", "metadata/taxonomy.tsv",
      "taxonomy/gtdb_taxonomy.tsv",
    ])
rows = []
present = set()
for rel in paths:
    path = os.path.join(root, rel)
    try:
        st = os.stat(path)
        if os.path.isfile(path):
            rows.append((rel, st.st_size, st.st_mtime_ns))
            present.add(rel)
    except OSError:
        continue
for rel, size, mtime in sorted(rows, key=lambda x: x[0].encode("utf-8")):
    print(f"{rel}\t{size}\t{mtime}")
missing = sorted(set(candidates.get(kind, [])) - present)
if missing:
    print("#WARN\toptional keyfiles omitted: " + ", ".join(missing))
PY
)"
    warning="$(sed -n 's/^#WARN	//p' <<< "${output}")"
    [[ -z "${warning}" ]] || log_warn "[Inventory] ${kind} database: ${warning}."
    sed '/^#WARN	/d' <<< "${output}"
}

resolve_database_identity() {
    local stage="$1" kind="$2" path="$3"
    local key="${stage}|${kind}" resolved version=unknown version_source=unknown method=unknown fingerprint=unknown
    local size=0 mtime=0 status=unavailable rows="" release_file="" first_row
    [[ -n "${path}" && -e "${path}" ]] || {
        ID_DATABASE_PATH["$key"]="$(realpath -m "${path:-unknown}")"
        ID_DATABASE_RELEASE["$key"]=unknown
        ID_DATABASE_VERSION_SOURCE["$key"]=unknown
        ID_DATABASE_SIZE["$key"]=not-applicable
        ID_DATABASE_MTIME["$key"]=not-applicable
        ID_DATABASE_METHOD["$key"]=unknown
        ID_DATABASE_FINGERPRINT["$key"]=unknown
        ID_DATABASE_STATUS["$key"]=unavailable
        return 0
    }
    resolved="$(realpath "${path}")"
    if [[ "${kind}" == checkm2 ]]; then
        method=sha256-file
    else
        method=metadata-keyfiles-v1
    fi
    local cache_key="${kind}|${resolved}|${method}"
    if [[ -n "${DB_IDENTITY_CACHE[$cache_key]:-}" ]]; then
        IFS='|' read -r version version_source size mtime method fingerprint status <<< "${DB_IDENTITY_CACHE[$cache_key]}"
    elif [[ -f "${resolved}" && "${kind}" == checkm2 ]]; then
        local stat_fields
        stat_fields="$(python3 -c 'import os,sys; s=os.stat(sys.argv[1]); print(f"{s.st_size}\t{s.st_mtime_ns}")' "${resolved}")"
        DB_SNAPSHOT["$cache_key"]="${stat_fields}"
        IFS=$'\t' read -r size mtime <<< "${stat_fields}"
        fingerprint="$(sha256sum "${resolved}" | awk '{print tolower($1)}')"
        method=sha256-file
        status=available
        DB_IDENTITY_CACHE["$cache_key"]="${version}|${version_source}|${size}|${mtime}|${method}|${fingerprint}|${status}"
    elif [[ -d "${resolved}" && ( "${kind}" == kraken2-bracken || "${kind}" == gtdbtk ) ]]; then
        rows="$(database_keyfile_rows "${kind}" "${resolved}")"
        DB_SNAPSHOT["$cache_key"]="${rows}"
        if [[ -n "${rows}" ]]; then
            size=0
            mtime=0
            while IFS=$'\t' read -r _relative _size _mtime; do
                size=$((size + _size))
                (( _mtime > mtime )) && mtime="${_mtime}"
            done <<< "${rows}"
            fingerprint="$(printf '%s\n' "${rows}" | sha256sum | awk '{print tolower($1)}')"
            method=metadata-keyfiles-v1
            status=available
            local metadata_name
            for metadata_name in VERSION release.txt gtdbtk.package_version; do
                if [[ -f "${resolved}/${metadata_name}" ]]; then
                    release_file="${resolved}/${metadata_name}"
                    break
                fi
            done
            if [[ -n "${release_file}" ]]; then
                first_row="$(head -n 1 "${release_file}" | tr -d '\r\n')"
                if [[ -n "${first_row}" ]]; then
                    version="$(tsv_value "${first_row}")"
                    version_source=metadata
                fi
            fi
            DB_IDENTITY_CACHE["$cache_key"]="${version}|${version_source}|${size}|${mtime}|${method}|${fingerprint}|${status}"
        else
            status=unknown
            method=unknown
            fingerprint=unknown
            log_warn "[Inventory] No trusted keyfiles found for ${kind} database at ${resolved}."
            DB_IDENTITY_CACHE["$cache_key"]="${version}|${version_source}|${size}|${mtime}|${method}|${fingerprint}|${status}"
        fi
    else
        status=unknown
        method=unknown
        fingerprint=unknown
        DB_IDENTITY_CACHE["$cache_key"]="${version}|${version_source}|${size}|${mtime}|${method}|${fingerprint}|${status}"
    fi
    ID_DATABASE_PATH["$key"]="$(tsv_value "${resolved}")"
    ID_DATABASE_RELEASE["$key"]="$(tsv_value "${version}")"
    ID_DATABASE_VERSION_SOURCE["$key"]="${version_source}"
    ID_DATABASE_SIZE["$key"]="${size}"
    ID_DATABASE_MTIME["$key"]="${mtime}"
    ID_DATABASE_METHOD["$key"]="${method}"
    ID_DATABASE_FINGERPRINT["$key"]="${fingerprint}"
    ID_DATABASE_STATUS["$key"]="${status}"
}

resolve_all_current_identities() {
    local stage tool
    local -a tools=()
    for stage in "${STAGE_IDS[@]}"; do
        if ! stage_is_active "${stage}" || ! stage_in_scope "$(get_stage_public_name "${stage}")"; then
            continue
        fi
        mapfile -t tools < <(inventory_stage_tools "${stage}")
        for tool in "${tools[@]}"; do
            runtime_identity_for_tool "${stage}" "${tool}"
        done
        case "${stage}" in
            03_kraken2) resolve_database_identity "${stage}" kraken2-bracken "${KRAKEN2_DB}" ;;
            08_checkm2) resolve_database_identity "${stage}" checkm2 "${CHECKM2_DB}" ;;
            10_gtdbtk) resolve_database_identity "${stage}" gtdbtk "${GTDBTK_DATA_PATH}" ;;
        esac
    done
}

append_checkpoint_identity() {
    local stage="$1" tool key tools_joined="" source_joined="" mode="" name="" uri=""
    local path="" size="" mtime="" method="" fingerprint="" database_key=""
    local -a tools=()
    mapfile -t tools < <(inventory_stage_tools "${stage}")
    for tool in "${tools[@]}"; do
        key="${stage}|${tool}"
        [[ -n "${ID_STATUS[$key]:-}" ]] || continue
        [[ -n "${tools_joined}" ]] && tools_joined+=","
        [[ -n "${source_joined}" ]] && source_joined+=","
        tools_joined+="${tool}=${ID_VERSION[$key]}"
        source_joined+="${tool}=${ID_VERSION_SOURCE[$key]}"
        mode="${ID_EXECUTION_MODE[$key]}"
        name="${ID_RUNTIME_NAME[$key]}"
        uri="${ID_RUNTIME_URI[$key]}"
        path="${ID_RUNTIME_PATH[$key]}"
        size="${ID_RUNTIME_SIZE[$key]}"
        mtime="${ID_RUNTIME_MTIME[$key]}"
        method="${ID_RUNTIME_METHOD[$key]}"
        fingerprint="${ID_RUNTIME_FINGERPRINT[$key]}"
        printf 'resource_%s_version=%s\nresource_%s_version_source=%s\nresource_%s_runtime_fingerprint_method=%s\nresource_%s_runtime_fingerprint=%s\n' \
            "${tool}" "${ID_VERSION[$key]}" "${tool}" "${ID_VERSION_SOURCE[$key]}" \
            "${tool}" "${method}" "${tool}" "${fingerprint}"
        printf 'resource_%s_command_used=%s\n' "${tool}" "${ID_COMMAND[$key]}"
    done
    [[ -z "${tools_joined}" ]] || printf 'effective_tool_version=%s\ntool_version_source=%s\nexecution_mode=%s\nruntime_name=%s\nruntime_uri=%s\nruntime_path=%s\nruntime_size_bytes=%s\nruntime_mtime_ns=%s\nruntime_fingerprint_method=%s\nruntime_fingerprint=%s\n' \
        "${tools_joined}" "${source_joined}" "${mode}" "${name}" "${uri}" "${path}" "${size}" "${mtime}" "${method}" "${fingerprint}"
    case "${stage}" in
        03_kraken2) database_key="${stage}|kraken2-bracken" ;;
        08_checkm2) database_key="${stage}|checkm2" ;;
        10_gtdbtk) database_key="${stage}|gtdbtk" ;;
    esac
    if [[ -n "${database_key}" && -n "${ID_DATABASE_STATUS[$database_key]:-}" ]]; then
        printf 'database_release_or_version=%s\ndatabase_version_source=%s\ndatabase_path=%s\ndatabase_fingerprint_method=%s\ndatabase_fingerprint=%s\n' \
            "${ID_DATABASE_RELEASE[$database_key]}" "${ID_DATABASE_VERSION_SOURCE[$database_key]}" \
            "${ID_DATABASE_PATH[$database_key]}" "${ID_DATABASE_METHOD[$database_key]}" \
            "${ID_DATABASE_FINGERPRINT[$database_key]}"
    fi
}

build_inventory_rows() {
    local stage tool key role dbkey producer_file producer_version producer_method producer_fingerprint row
    local producer_name producer_uri producer_path producer_size producer_mtime producer_mode producer_command
    local producer_release producer_dbpath producer_dbmethod producer_dbfingerprint producer_dbsize producer_dbmtime
    SOFTWARE_INVENTORY_ROWS=()
    DATABASE_INVENTORY_ROWS=()
    INVENTORY_UNKNOWN_COUNT=0
    INVENTORY_UNAVAILABLE_COUNT=0
    for stage in "${STAGE_IDS[@]}"; do
        if ! stage_is_active "${stage}" || ! stage_in_scope "$(get_stage_public_name "${stage}")"; then
            continue
        fi
        case "${STAGE_ACTION[$stage]}" in
            run|force|restart) role="current-runtime" ;;
            reuse) role="current-runtime" ;;
            *) continue ;;
        esac
        while IFS= read -r tool; do
            [[ -n "${tool}" ]] || continue
            key="${stage}|${tool}"
            [[ -n "${ID_STATUS[$key]:-}" ]] || continue
            row="$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
                "${stage}" "${tool}" "${role}" "${ID_VERSION[$key]}" "${ID_VERSION_SOURCE[$key]}" \
                "${ID_EXECUTION_MODE[$key]}" "${ID_RUNTIME_NAME[$key]}" "${ID_RUNTIME_URI[$key]}" \
                "${ID_RUNTIME_PATH[$key]}" "${ID_RUNTIME_SIZE[$key]}" "${ID_RUNTIME_MTIME[$key]}" \
                "${ID_RUNTIME_METHOD[$key]}" "${ID_RUNTIME_FINGERPRINT[$key]}" "${ID_COMMAND[$key]}")"
            SOFTWARE_INVENTORY_ROWS+=("${row}"$'\t'"${ID_STATUS[$key]}")
            if [[ "${ID_STATUS[$key]}" == unknown ]]; then
                INVENTORY_UNKNOWN_COUNT=$((INVENTORY_UNKNOWN_COUNT + 1))
            elif [[ "${ID_STATUS[$key]}" == unavailable ]]; then
                INVENTORY_UNAVAILABLE_COUNT=$((INVENTORY_UNAVAILABLE_COUNT + 1))
            fi
        done < <(inventory_stage_tools "${stage}")
        if [[ "${STAGE_ACTION[$stage]}" == reuse ]]; then
            producer_file="${CHECKPOINT_DIR}/${stage}.config"
            producer_mode="$(awk -F= '$1=="engine" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            [[ "${producer_mode}" == apptainer || "${producer_mode}" == mamba ]] || producer_mode="${EXEC_MODE}"
            producer_name="$(awk -F= '$1=="runtime_name" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            producer_uri="$(awk -F= '$1=="runtime_uri" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            producer_path="$(awk -F= '$1=="runtime_path" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            producer_size="$(awk -F= '$1=="runtime_size_bytes" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            producer_mtime="$(awk -F= '$1=="runtime_mtime_ns" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
            [[ -n "${producer_name}" ]] || producer_name=unknown
            [[ -n "${producer_uri}" ]] || producer_uri=unknown
            [[ -n "${producer_path}" ]] || producer_path=unknown
            [[ -n "${producer_size}" ]] || producer_size=unknown
            [[ -n "${producer_mtime}" ]] || producer_mtime=unknown
            while IFS= read -r tool; do
                [[ -n "${tool}" ]] || continue
                producer_version="$(awk -F= -v key="resource_${tool}_version" '$1==key {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_method="$(awk -F= -v key="resource_${tool}_runtime_fingerprint_method" '$1==key {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_fingerprint="$(awk -F= -v key="resource_${tool}_runtime_fingerprint" '$1==key {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_command="$(awk -F= -v key="resource_${tool}_command_used" '$1==key {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                [[ -n "${producer_version}" ]] || producer_version=unknown
                [[ -n "${producer_method}" ]] || producer_method=unknown
                [[ -n "${producer_fingerprint}" ]] || producer_fingerprint=unknown
                [[ -n "${producer_command}" ]] || producer_command=unknown
                producer_version="$(tsv_value "${producer_version}")"
                producer_command="$(tsv_value "${producer_command}")"
                SOFTWARE_INVENTORY_ROWS+=("$(printf '%s\t%s\tproducer\t%s\tcheckpoint-manifest\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\treused' \
                    "${stage}" "${tool}" "${producer_version}" "${producer_mode}" \
                    "$(tsv_value "${producer_name}")" "$(tsv_value "${producer_uri}")" "$(tsv_value "${producer_path}")" "${producer_size}" \
                    "${producer_mtime}" "${producer_method}" "${producer_fingerprint}" "${producer_command}")")
                if [[ "${producer_version}" == unknown || "${producer_method}" == unknown || "${producer_fingerprint}" == unknown ]]; then
                    INVENTORY_UNKNOWN_COUNT=$((INVENTORY_UNKNOWN_COUNT + 1))
                fi
            done < <(inventory_stage_tools "${stage}")
        fi
        case "${stage}" in
            03_kraken2) dbkey="${stage}|kraken2-bracken" ;;
            08_checkm2) dbkey="${stage}|checkm2" ;;
            10_gtdbtk) dbkey="${stage}|gtdbtk" ;;
            *) dbkey="" ;;
        esac
        if [[ -n "${dbkey}" && -n "${ID_DATABASE_STATUS[$dbkey]:-}" ]]; then
            DATABASE_INVENTORY_ROWS+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
                "${stage}" "${dbkey#*|}" "${role}" "${ID_DATABASE_PATH[$dbkey]}" \
                "${ID_DATABASE_RELEASE[$dbkey]}" "${ID_DATABASE_VERSION_SOURCE[$dbkey]}" \
                "${ID_DATABASE_SIZE[$dbkey]}" "${ID_DATABASE_MTIME[$dbkey]}" \
                "${ID_DATABASE_METHOD[$dbkey]}" "${ID_DATABASE_FINGERPRINT[$dbkey]}" \
                "${ID_DATABASE_STATUS[$dbkey]}")")
            if [[ "${ID_DATABASE_STATUS[$dbkey]}" == unknown ]]; then
                INVENTORY_UNKNOWN_COUNT=$((INVENTORY_UNKNOWN_COUNT + 1))
            elif [[ "${ID_DATABASE_STATUS[$dbkey]}" == unavailable ]]; then
                INVENTORY_UNAVAILABLE_COUNT=$((INVENTORY_UNAVAILABLE_COUNT + 1))
            fi
            if [[ "${STAGE_ACTION[$stage]}" == reuse ]]; then
                producer_file="${CHECKPOINT_DIR}/${stage}.config"
                producer_release="$(awk -F= '$1=="database_release_or_version" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_dbpath="$(awk -F= '$1=="database_path" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_dbmethod="$(awk -F= '$1=="database_fingerprint_method" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                producer_dbfingerprint="$(awk -F= '$1=="database_fingerprint" {sub(/^[^=]*=/,""); print; exit}' "${producer_file}" 2>/dev/null)"
                [[ -n "${producer_release}" ]] || producer_release=unknown
                [[ -n "${producer_dbpath}" ]] || producer_dbpath=unknown
                [[ -n "${producer_dbmethod}" ]] || producer_dbmethod=unknown
                [[ -n "${producer_dbfingerprint}" ]] || producer_dbfingerprint=unknown
                producer_release="$(tsv_value "${producer_release}")"
                producer_dbpath="$(tsv_value "${producer_dbpath}")"
                producer_dbsize=unknown
                producer_dbmtime=unknown
                DATABASE_INVENTORY_ROWS+=("$(printf '%s\t%s\tproducer\t%s\t%s\tcheckpoint-manifest\t%s\t%s\t%s\t%s\treused' \
                    "${stage}" "${dbkey#*|}" "${producer_dbpath}" "${producer_release}" \
                    "${producer_dbsize}" "${producer_dbmtime}" "${producer_dbmethod}" "${producer_dbfingerprint}")")
                if [[ "${producer_release}" == unknown || "${producer_dbmethod}" == unknown || "${producer_dbfingerprint}" == unknown ]]; then
                    INVENTORY_UNKNOWN_COUNT=$((INVENTORY_UNKNOWN_COUNT + 1))
                fi
            fi
        fi
    done
}

validate_current_identity_plan() {
    local stage tool key dbkey failures=0
    for stage in "${STAGE_IDS[@]}"; do
        if ! stage_is_active "${stage}" || ! stage_in_scope "$(get_stage_public_name "${stage}")"; then
            continue
        fi
        stage_action_requires_execution "${stage}" || continue
        while IFS= read -r tool; do
            [[ -n "${tool}" ]] || continue
            key="${stage}|${tool}"
            if [[ "${ID_STATUS[$key]:-unavailable}" == unavailable ]]; then
                log_error "Required runtime for ${stage}/${tool} is unavailable."
                (( failures += 1 ))
            elif [[ "${ID_VERSION[$key]:-unknown}" == unknown && \
                    "${ID_RUNTIME_FINGERPRINT[$key]:-unknown}" == unknown ]]; then
                log_error "Runtime identity for ${stage}/${tool} has neither a known version nor a fingerprint."
                (( failures += 1 ))
            fi
        done < <(inventory_stage_tools "${stage}")
        case "${stage}" in
            03_kraken2) dbkey="${stage}|kraken2-bracken" ;;
            08_checkm2) dbkey="${stage}|checkm2" ;;
            10_gtdbtk) dbkey="${stage}|gtdbtk" ;;
            *) dbkey="" ;;
        esac
        if [[ -n "${dbkey}" && "${ID_DATABASE_STATUS[$dbkey]:-unavailable}" == unavailable ]]; then
            log_error "Required database for ${stage} is unavailable: ${ID_DATABASE_PATH[$dbkey]:-unknown}"
            (( failures += 1 ))
        fi
    done
    (( failures == 0 ))
}

software_inventory_header() {
    printf 'stage_id\ttool\trole\tversion\tversion_source\texecution_mode\truntime_name\truntime_uri\truntime_path\truntime_size_bytes\truntime_mtime_ns\truntime_fingerprint_method\truntime_fingerprint\tcommand_used\tstatus\n'
}

database_inventory_header() {
    printf 'stage_id\tdatabase\trole\tpath\trelease_or_version\tversion_source\tsize_bytes\tmtime_ns\tfingerprint_method\tfingerprint\tstatus\n'
}

validate_inventory_file() {
    local file="$1" kind="$2" header expected
    if [[ "${kind}" == software ]]; then
        expected="$(software_inventory_header)"
        awk -F'\t' 'NR>1 && NF!=15 {exit 1}' "${file}" || return 1
        awk -F'\t' 'NR>1 && ($3 !~ /^(producer|current-runtime)$/ || $5 !~ /^(runtime-command|checkpoint-manifest|image-label|commit|unknown)$/ || $6 !~ /^(apptainer|mamba)$/ || $12 !~ /^(sha256-file|conda-meta-v1|unknown)$/ || $15 !~ /^(available|unknown|unavailable|reused)$/) {exit 1}' "${file}" || return 1
        awk -F'\t' 'NR>1 {k=$1 FS $2 FS $3; if(seen[k]++) exit 1}' "${file}" || return 1
    else
        expected="$(database_inventory_header)"
        awk -F'\t' 'NR>1 && NF!=11 {exit 1}' "${file}" || return 1
        awk -F'\t' 'NR>1 && ($2 !~ /^(kraken2-bracken|checkm2|gtdbtk)$/ || $3 !~ /^(producer|current-runtime)$/ || $6 !~ /^(metadata|checkpoint-manifest|unknown)$/ || $9 !~ /^(sha256-file|metadata-keyfiles-v1|unknown)$/ || $11 !~ /^(available|unknown|unavailable|reused)$/) {exit 1}' "${file}" || return 1
        awk -F'\t' 'NR>1 {k=$1 FS $2 FS $3; if(seen[k]++) exit 1}' "${file}" || return 1
    fi
    IFS= read -r header < "${file}"
    [[ "${header}" == "${expected%$'\n'}" ]]
}

publish_inventory_generation() {
    local logs root build previous public target temp outdir_abs
    logs="$(realpath -m "${LOGS_DIR}")"
    outdir_abs="$(realpath -m "${OUTDIR}")"
    [[ "${logs}" == "${outdir_abs}/logs" || "${logs}" == "${outdir_abs}/logs/"* ]] || {
        log_error "Inventory logs path resolves outside OUTDIR: ${logs}"
        return 1
    }
    root="$(realpath -m "${logs}/inventory")"
    [[ "${root}" == "${logs}/inventory" || "${root}" == "${logs}/inventory/"* || \
       "${root}" == "${logs}/"* ]] || {
        log_error "Inventory directory resolves outside logs: ${root}"
        return 1
    }
    build="${root}/.build.$$"
    previous="${root}/.previous.$$"
    if [[ "${DRY_RUN}" == true ]]; then
        log_info "Software inventory      : ${logs}/software_versions.tsv (would write)"
        log_info "Database inventory      : ${logs}/database_versions.tsv (would write)"
        log_info "Software entries        : ${#SOFTWARE_INVENTORY_ROWS[@]}"
        log_info "Database entries        : ${#DATABASE_INVENTORY_ROWS[@]}"
        log_info "Unknown identities      : ${INVENTORY_UNKNOWN_COUNT}"
        log_info "Unavailable resources   : ${INVENTORY_UNAVAILABLE_COUNT}"
        local _inv_row
        for _inv_row in "${SOFTWARE_INVENTORY_ROWS[@]}"; do
            log_info "  [Software] ${_inv_row}"
        done
        for _inv_row in "${DATABASE_INVENTORY_ROWS[@]}"; do
            log_info "  [Database] ${_inv_row}"
        done
        return 0
    fi
    mkdir -p "${root}"
    [[ ! -e "${build}" && ! -e "${previous}" ]] || {
        log_error "Inventory generation path already exists for PID $$."
        return 1
    }
    mkdir "${build}"
    {
        software_inventory_header
        (( ${#SOFTWARE_INVENTORY_ROWS[@]} == 0 )) || printf '%s\n' "${SOFTWARE_INVENTORY_ROWS[@]}"
    } > "${build}/software_versions.tsv"
    {
        database_inventory_header
        (( ${#DATABASE_INVENTORY_ROWS[@]} == 0 )) || printf '%s\n' "${DATABASE_INVENTORY_ROWS[@]}"
    } > "${build}/database_versions.tsv"
    if ! validate_inventory_file "${build}/software_versions.tsv" software || \
       ! validate_inventory_file "${build}/database_versions.tsv" database; then
        rm -rf -- "${build}"
        log_error "Refusing to publish an invalid inventory generation."
        return 1
    fi
    if [[ -e "${root}/current" ]]; then
        mv -- "${root}/current" "${previous}" || return 1
    fi
    if ! mv -- "${build}" "${root}/current"; then
        [[ -e "${previous}" ]] && mv -- "${previous}" "${root}/current"
        return 1
    fi
    for public in software_versions.tsv database_versions.tsv; do
        target="inventory/current/${public}"
        temp="${logs}/.${public}.$$"
        rm -f -- "${temp}"
        ln -s "${target}" "${temp}" || {
            [[ -e "${previous}" ]] && { rm -rf -- "${root}/current"; mv -- "${previous}" "${root}/current"; }
            return 1
        }
        mv -Tf -- "${temp}" "${logs}/${public}" || {
            rm -f -- "${temp}"
            [[ -e "${previous}" ]] && { rm -rf -- "${root}/current"; mv -- "${previous}" "${root}/current"; }
            return 1
        }
        [[ "$(realpath -m "${logs}/${public}")" == "${root}/current/${public}" ]] || {
            log_error "Public inventory link escapes logs directory: ${logs}/${public}"
            return 1
        }
    done
    if ! validate_inventory_file "${root}/current/software_versions.tsv" software || \
       ! validate_inventory_file "${root}/current/database_versions.tsv" database; then
        log_error "Post-publication validation failed for inventory generation."
        if [[ -e "${previous}" ]]; then
            rm -rf -- "${root}/current"
            mv -- "${previous}" "${root}/current"
            for public in software_versions.tsv database_versions.tsv; do
                target="inventory/current/${public}"
                temp="${logs}/.${public}.$$"
                rm -f -- "${temp}"
                ln -s "${target}" "${temp}" && mv -Tf -- "${temp}" "${logs}/${public}" || true
            done
        fi
        return 1
    fi
    [[ -e "${previous}" ]] && rm -rf -- "${previous}"
    log_info "Software inventory      : ${logs}/software_versions.tsv"
    log_info "Database inventory      : ${logs}/database_versions.tsv"
    log_info "Software entries        : ${#SOFTWARE_INVENTORY_ROWS[@]}"
    log_info "Database entries        : ${#DATABASE_INVENTORY_ROWS[@]}"
    log_info "Unknown identities      : ${INVENTORY_UNKNOWN_COUNT}"
    log_info "Unavailable resources   : ${INVENTORY_UNAVAILABLE_COUNT}"
}

assert_runtime_identity_frozen() {
    local stage="$1" tool spec sif uri executable env_name path stat_fields size mtime object key
    local -a tools=()
    mapfile -t tools < <(inventory_stage_tools "${stage}")
    for tool in "${tools[@]}"; do
        key="${stage}|${tool}"
        [[ "${ID_STATUS[$key]:-unavailable}" == unavailable ]] && continue
        if [[ "${EXEC_MODE}" == mamba ]]; then
            path="${ID_RUNTIME_PATH[$key]:-unknown}"
            [[ "${path}" == unknown ]] && continue
            [[ -d "${path}" ]] || {
                log_error "Mamba runtime prefix disappeared after planning for ${stage}/${tool}."
                return 1
            }
            object="$(mamba_environment_snapshot "${path}")"
            if [[ "${object}" != "${MAMBA_OBJECT_CACHE[$path]:-}" ]]; then
                log_error "Mamba environment changed after the plan was frozen for ${stage}/${tool}."
                return 1
            fi
            continue
        fi
        spec="$(tool_runtime_spec "${stage}" "${tool}")"
        IFS=$'\t' read -r sif uri executable env_name <<< "${spec}"
        [[ -f "${sif}" && -r "${sif}" ]] || {
            log_error "Runtime image became unavailable after planning for ${stage}/${tool}."
            return 1
        }
        path="$(realpath "${sif}")"
        stat_fields="$(python3 -c 'import os,sys; s=os.stat(sys.argv[1]); print(f"{s.st_size}\t{s.st_mtime_ns}\t{s.st_ctime_ns}\t{s.st_dev}\t{s.st_ino}")' "${path}")"
        IFS=$'\t' read -r size mtime <<< "${stat_fields}"
        object="$(cut -f3- <<< "${stat_fields}")"
        if [[ "${path}" != "${ID_RUNTIME_PATH[$key]}" || \
              "${size}" != "${ID_RUNTIME_SIZE[$key]}" || \
              "${mtime}" != "${ID_RUNTIME_MTIME[$key]}" || \
              "${object}" != "${SIF_OBJECT_CACHE[$path]:-}" ]]; then
            log_error "Runtime image identity changed after the plan was frozen for ${stage}/${tool}; refusing to replan."
            return 1
        fi
    done
}

assert_database_identity_frozen() {
    local stage="$1" kind="" path="" key="" resolved="" cache_key="" snapshot="" current=""
    case "${stage}" in
        03_kraken2) kind=kraken2-bracken; path="${KRAKEN2_DB}" ;;
        08_checkm2) kind=checkm2; path="${CHECKM2_DB}" ;;
        10_gtdbtk) kind=gtdbtk; path="${GTDBTK_DATA_PATH}" ;;
        *) return 0 ;;
    esac
    key="${stage}|${kind}"
    [[ "${ID_DATABASE_STATUS[$key]:-unavailable}" == unavailable ]] && return 0
    [[ -e "${path}" ]] || {
        log_error "Database became unavailable after planning for ${stage}."
        return 1
    }
    resolved="$(realpath "${path}")"
    [[ "${resolved}" == "${ID_DATABASE_PATH[$key]}" ]] || {
        log_error "Database path changed after the plan was frozen for ${stage}."
        return 1
    }
    cache_key="${kind}|${resolved}"
    snapshot="${DB_SNAPSHOT[$cache_key]:-}"
    if [[ "${kind}" == checkm2 ]]; then
        current="$(python3 -c 'import os,sys; s=os.stat(sys.argv[1]); print(f"{s.st_size}\t{s.st_mtime_ns}")' "${resolved}")"
    elif [[ -d "${resolved}" ]]; then
        current="$(database_keyfile_rows "${kind}" "${resolved}")"
    fi
    if [[ "${current}" != "${snapshot}" ]]; then
        log_error "Database identity changed after the plan was frozen for ${stage}; refusing to replan."
        return 1
    fi
}

checkpoint_config_manifest() {
    local step="$1"
    {
        printf 'step=%s\nengine=%s\n' "${step}" "${EXEC_MODE}"
        append_checkpoint_identity "${step}"
        case "${step}" in
            01_fastp|02_megahit|03_kraken2|04_mapping|07_reassembly|11_coverm)
                append_reads_manifest "${step}" ;;
        esac
        case "${step}" in
            04_mapping|05a_vamb|05b_metabinner|05c_comebin)
                append_assembly_manifest ;;
        esac
        case "${step}" in
            01_fastp)
                printf 'enabled=%s\nthreads=%s\nparams=-5 20 -3 20 -r 20 -M 20 -e 20 -n 50\nimage=%s\n' \
                    "${RUN_FASTP}" "${THREADS}" "${SIF_FASTP_URI}" ;;
            02_megahit)
                printf 'threads=%s\nmemory=%s\nmin_contig_len=%s\nimage=%s\n' \
                    "${THREADS}" "${MEMORY}" "${MIN_CONTIG_LEN}" "${SIF_MEGAHIT_URI}" ;;
            03_kraken2)
                printf 'db=%s\nthreads=%s\nread_len=%s\ntax_level=%s\nthreshold=%s\nimage=%s\n' \
                    "$(realpath -m "${KRAKEN2_DB}")" "${THREADS}" "${BRACKEN_READ_LEN}" \
                    "${BRACKEN_TAX_LEVEL}" "${BRACKEN_THRESHOLD}" "${SIF_KRAKEN2_URI}" ;;
            04_mapping)
                printf 'enabled=%s\ncontigs=%s\nthreads=%s\nimage=%s\n' \
                    "${RUN_MAPPING}" "${CONTIGS}" "${THREADS}" "${SIF_MAPPING_URI}"
                append_file_fingerprint contigs "${CONTIGS}" ;;
            05a_vamb)
                printf 'coverage_source=%s\ncontigs=%s\ndepth=%s\nabundance_tsv=%s\noutput_dir=%s\nfinal_bins=%s\n' \
                    "${COVERAGE_SOURCE}" \
                    "${CONTIGS}" "${DEPTH_FILE}" "${VAMB_ABUNDANCE_TSV}" \
                    "${VAMB_DIR}/vamb_out/bins" "$(realpath -m "${VAMB_FINAL_BINS_DIR}")"
                printf 'vamb_version=%s\nvamb_runtime=%s\nimage_uri=%s\n' \
                    "${VAMB_VERSION}" "${EXEC_MODE}" "${SIF_VAMB_URI}"
                printf 'vamb_binsplit_separator=%q\nbinsplitting=%s\n' \
                    "${VAMB_BINSPLIT_SEPARATOR}" \
                    "$([[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]] && printf enabled || printf disabled)"
                printf 'vamb_minfasta=%s\nvamb_minfasta_meaning=minimum_bin_fasta_size_bp\n' \
                    "${VAMB_MINFASTA}"
                printf 'vamb_internal_min_contig_len=%s\nvamb_internal_min_contig_parameter=default_-m\n' \
                    "${VAMB_INTERNAL_MIN_CONTIG_LEN}"
                printf 'vamb_default_threads=%s\nvamb_cuda=false\nvamb_norefcheck=false\nvamb_seed=unset\n' \
                    "${VAMB_DEFAULT_THREADS}"
                printf 'vamb_parameters=bin_default,--minfasta=%s; no_explicit_-m_or_-p; no_cuda; no_seed\n' \
                    "${VAMB_MINFASTA}"
                append_file_fingerprint contigs "${CONTIGS}"
                append_file_fingerprint depth "${DEPTH_FILE}"
                append_file_fingerprint abundance_tsv "${VAMB_ABUNDANCE_TSV}"
                if [[ "${EXEC_MODE}" != apptainer ]]; then
                    printf 'mamba_environment=%s\n' "${ENV_VAMB}"
                fi
                append_vamb_abundance_manifest
                append_assembly_manifest
                append_bins_manifest "${VAMB_FINAL_BINS_DIR}" ;;
            05b_metabinner)
                printf 'coverage_source=%s\ncontigs=%s\ndepth=%s\noutput_dir=%s\nfinal_bins=%s\nthreads=%s\nmin_contig_len=%s\nscripts=%s\nrunner=%s\nimage=%s\n' \
                    "${COVERAGE_SOURCE}" "${CONTIGS}" "${DEPTH_FILE}" "${METABINNER_DIR}/metabinner_out" "${METABINNER_FINAL_BINS_DIR}" \
                    "${THREADS}" "${MIN_CONTIG_LEN}" "${METABINNER_SCRIPTS_DIR}" "${METABINNER_RUNNER}" "${SIF_METABINNER_URI}"
                append_file_fingerprint contigs "${CONTIGS}"
                append_file_fingerprint depth "${DEPTH_FILE}"
                append_assembly_manifest
                append_bins_manifest "${METABINNER_FINAL_BINS_DIR}" ;;
            05c_comebin)
                local comebin_bam_dir="${COMEBIN_BAM_DIR:-${BAM_DIR_ACTIVE}}"
                local comebin_output_dir
                comebin_output_dir="$(realpath -m "${COMEBIN_DIR}")"
                [[ "${COVERAGE_SOURCE}" == external && -z "${comebin_bam_dir}" ]] && \
                    comebin_bam_dir="${COMEBIN_DIR}/external_bams"
                [[ -n "${comebin_bam_dir}" ]] && \
                    comebin_bam_dir="$(realpath -m "${comebin_bam_dir}")"
                printf 'coverage_source=%s\ncontigs=%s\nbam_dir=%s\noutput_dir=%s\nfinal_bins=%s\n' \
                    "${COVERAGE_SOURCE}" "${CONTIGS}" "${comebin_bam_dir}" "${comebin_output_dir}/comebin_res/comebin_res_bins" \
                    "${COMEBIN_FINAL_BINS_DIR}"
                printf 'comebin_version=%s\ncomebin_device=%s\ncomebin_seed=%s\nthreads=%s\n' \
                    "${COMEBIN_VERSION}" "${COMEBIN_DEVICE}" "${COMEBIN_SEED}" "${THREADS}"
                printf 'pytorch_version=%s\npytorch_cuda_build=%s\npytorch_cuda_available=%s\ncuda_device_count=%s\n' \
                    "${COMEBIN_PYTORCH_VERSION}" "${COMEBIN_PYTORCH_CUDA_BUILD}" \
                    "${COMEBIN_PYTORCH_CUDA_AVAILABLE}" "${COMEBIN_CUDA_DEVICE_COUNT}"
                printf 'runtime=%s\nimage_uri=%s\ncomebin_parameters=%q\n' \
                    "${EXEC_MODE}" "${SIF_COMEBIN_URI}" \
                    "-a ${CONTIGS} -o ${comebin_output_dir} -p ${comebin_bam_dir} -s ${COMEBIN_SEED} -t ${THREADS} -d ${COMEBIN_DEVICE}; wrapper defaults: -n 6 -w same_as_-t -m 100 -E 1e-5 -l 0.07_if_N50_gt_10000_else_0.15 -e 2048 -c 2048 -b 1024"
                append_file_fingerprint contigs "${CONTIGS}"
                if [[ "${EXEC_MODE}" != apptainer ]]; then
                    printf 'mamba_environment=%s\n' "${ENV_COMEBIN}"
                fi
                append_assembly_manifest
                if [[ "${COVERAGE_SOURCE}" == external ]]; then
                    local coverage_index
                    for (( coverage_index=0; coverage_index<${#ACTIVE_BAMS[@]}; coverage_index++ )); do
                        printf 'external_bam_sample=%s\n' "${ACTIVE_BAM_SAMPLE_IDS[$coverage_index]}"
                        append_file_fingerprint external_bam "${ACTIVE_BAMS[$coverage_index]}"
                        append_file_fingerprint external_bam_index "${ACTIVE_BAM_INDEXES[$coverage_index]}"
                    done
                else
                    local bam
                    while IFS= read -r bam; do
                        append_file_fingerprint mapping "${bam}"
                        if [[ -f "${bam}.bai" ]]; then
                            append_file_fingerprint mapping_index "${bam}.bai"
                        elif [[ -f "${bam}.csi" ]]; then
                            append_file_fingerprint mapping_index "${bam}.csi"
                        fi
                    done \
                        < <(find "${MAPPING_DIR}" -maxdepth 1 -type f -name '*_sorted.bam' -print | sort)
                fi
                append_bins_manifest "${COMEBIN_FINAL_BINS_DIR}" ;;
            06_refinement)
                printf 'refinement_mode=%s\nrefinement_effective=%s\nselected_binners=%s\nthreads=%s\ncompleteness=%s\ncontamination=%s\nimage=%s\n' \
                    "${REFINEMENT_MODE}" "${RUN_REFINEMENT_EFFECTIVE}" "${SELECTED_BINNERS}" "${THREADS}" "${METAWRAP_COMPLETENESS}" \
                    "${METAWRAP_CONTAMINATION}" "${SIF_METAWRAP_URI}"
                local method bins_dir
                for method in "${SELECTED_BINNER_ARRAY[@]}"; do
                    bins_dir="$(get_binner_bins_dir "${method}")"
                    printf 'input_binner=%s\n' "${method}"
                    append_bins_manifest "${bins_dir}"
                done
                append_bins_manifest "${REFINED_BINS_DIR}" ;;
            07_reassembly)
                printf 'enabled=%s\nbins_for_reassembly=%s\nbins_for_qc=%s\nthreads=%s\nmemory=%s\nimage=%s\n' \
                    "${RUN_REASSEMBLY}" "${BINS_FOR_REASSEMBLY}" "${BINS_FOR_QC}" "${THREADS}" "${MEMORY}" "${SIF_METAWRAP_URI}"
                printf 'temporary_reads_strategy=%s\ntemporary_reads_strategy_version=%s\n' \
                    "${REASSEMBLY_TEMP_READS_STRATEGY}" "${REASSEMBLY_TEMP_READS_STRATEGY_VERSION}"
                printf 'temporary_reads_validation=%s\ntemporary_reads_validation_version=%s\n' \
                    "${REASSEMBLY_TEMP_READS_VALIDATION}" "${REASSEMBLY_TEMP_READS_VALIDATION_VERSION}"
                printf 'active_sample_manifest_fingerprint=%s\n' \
                    "$(sha256sum "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}" 2>/dev/null | awk '{print $1}')"
                printf 'active_sample_ids=%s\n' "${ACTIVE_SAMPLE_IDS[*]:-}"
                local _src_r1_count=0 _src_r1_bytes=0 _src_r2_count=0 _src_r2_bytes=0 _rsf
                for _rsf in "${ACTIVE_R1_FILES[@]}"; do
                    _src_r1_count=$(( _src_r1_count + 1 ))
                    _src_r1_bytes=$(( _src_r1_bytes + $(stat -c '%s' -- "${_rsf}" 2>/dev/null || echo 0) ))
                done
                for _rsf in "${ACTIVE_R2_FILES[@]}"; do
                    _src_r2_count=$(( _src_r2_count + 1 ))
                    _src_r2_bytes=$(( _src_r2_bytes + $(stat -c '%s' -- "${_rsf}" 2>/dev/null || echo 0) ))
                done
                printf 'source_read_pair_count=%s\nread1_total_size_bytes=%s\nread2_total_size_bytes=%s\n' \
                    "${_src_r1_count}" "${_src_r1_bytes}" "${_src_r2_bytes}"
                append_bins_manifest "${BINS_FOR_REASSEMBLY}" ;;
            08_checkm2)
                printf 'selected_binners=%s\nrefinement_mode=%s\nrefinement_effective=%s\nbins_source_requested=%s\nbins_source_resolved=%s\nbins=%s\ndb=%s\nthreads=%s\nimage=%s\n' \
                    "${SELECTED_BINNERS}" "${REFINEMENT_MODE}" "${RUN_REFINEMENT_EFFECTIVE}" \
                    "${BINS_SOURCE}" "${BINS_SOURCE_RESOLVED}" "$(realpath -m "${BINS_FOR_QC}")" \
                    "$(realpath -m "${CHECKM2_DB}")" "${THREADS}" "${SIF_CHECKM2_URI}"
                append_bins_manifest "${BINS_FOR_QC}" ;;
            09_filter_mags)
                printf 'enabled=%s\nmin_completeness=%s\nmax_contamination=%s\ncomparison_completeness=>=\ncomparison_contamination=<=\n' \
                    "${RUN_MAG_FILTERING}" "${MIN_COMPLETENESS}" "${MAX_CONTAMINATION}"
                printf 'checkm2_report_fingerprint=%s\ninput_bins_fingerprint=%s\nretained_mags_fingerprint=%s\nfilter_manifest_fingerprint=%s\n' \
                    "$(sha256sum "${CHECKM2_DIR}/quality_report.tsv" 2>/dev/null | awk '{print $1}')" \
                    "$(append_fasta_bins_manifest input_bin "${BINS_FOR_QC}" | sha256sum | awk '{print $1}')" \
                    "$(append_fasta_bins_manifest retained_mag "${MAGS_DIR}" | sha256sum | awk '{print $1}')" \
                    "$(sha256sum "${MAG_FILTER_MANIFEST}" 2>/dev/null | awk '{print $1}')" ;;
            10_gtdbtk)
                printf 'enabled=%s\nmags=%s\ndata=%s\nthreads=%s\nimage=%s\n' \
                    "${RUN_GTDBTK}" "${HQ_BINS_DIR}" "${GTDBTK_DATA_PATH}" "${THREADS}" "${SIF_GTDBTK_URI}" ;;
            11_coverm)
                printf 'enabled=%s\nmags=%s\nmethods=%s\nmin_covered_fraction=%s\nthreads=%s\nimage=%s\n' \
                    "${RUN_COVERM}" "${HQ_BINS_DIR}" "${COVERM_METHODS}" "${COVERM_MIN_COVERED_FRACTION}" "${THREADS}" "${SIF_COVERM_URI}" ;;
            12a_taxonomy_summary)
                printf 'schema_version=1\n'
                printf 'generator_version=%s\n' "$(sha256sum "${SUMMARY_GENERATOR}" 2>/dev/null | awk '{print $1}')"
                printf 'sample_manifest_fingerprint=%s\n' "$(sha256sum "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}" 2>/dev/null | awk '{print $1}')"
                printf 'kraken2_enabled=%s\n' "${RUN_KRAKEN2}"
                printf 'bracken_enabled=%s\n' "${RUN_BRACKEN:-true}"
                printf 'bracken_read_len=%s\n' "${BRACKEN_READ_LEN}"
                printf 'bracken_tax_level=%s\n' "${BRACKEN_TAX_LEVEL}"
                printf 'bracken_threshold=%s\n' "${BRACKEN_THRESHOLD}"
                if [[ "${RUN_KRAKEN2}" == true ]]; then
                    local sample
                    for sample in "${ACTIVE_SAMPLE_IDS[@]}"; do
                        append_file_fingerprint "kraken_report_${sample}" "${KRAKEN2_DIR}/kraken_report_${sample}.txt"
                        if [[ "${RUN_BRACKEN:-true}" == true ]]; then
                            append_file_fingerprint "bracken_${sample}" "${KRAKEN2_DIR}/bracken_${sample}.txt"
                        fi
                    done
                fi ;;
            12b_mag_summary)
                printf 'schema_version=1\n'
                printf 'generator_version=%s\n' "$(sha256sum "${SUMMARY_GENERATOR}" 2>/dev/null | awk '{print $1}')"
                printf 'checkm2_enabled=%s\n' "${RUN_CHECKM2}"
                printf 'filtering_enabled=%s\n' "${RUN_MAG_FILTERING}"
                printf 'gtdbtk_enabled=%s\n' "${RUN_GTDBTK}"
                printf 'coverm_enabled=%s\n' "${RUN_COVERM}"
                printf 'min_completeness=%s\nmax_contamination=%s\n' "${MIN_COMPLETENESS}" "${MAX_CONTAMINATION}"
                printf 'coverm_methods=%s\n' "${COVERM_METHODS}"
                if [[ "${RUN_CHECKM2}" == true ]]; then
                    append_file_fingerprint checkm2_report "${CHECKM2_DIR}/quality_report.tsv"
                fi
                if [[ "${RUN_MAG_FILTERING}" == true ]]; then
                    append_file_fingerprint mag_filter_manifest "${MAG_FILTER_MANIFEST}"
                    append_fasta_bins_manifest retained_mag "${MAGS_DIR}"
                fi
                if [[ "${RUN_GTDBTK}" == true ]]; then
                    append_file_fingerprint gtdbtk_bac "${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv"
                    append_file_fingerprint gtdbtk_ar "${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv"
                fi
                if [[ "${RUN_COVERM}" == true ]]; then
                    local sample
                    for sample in "${ACTIVE_SAMPLE_IDS[@]}"; do
                        append_file_fingerprint "coverm_abundance_${sample}" "${ABUNDANCE_DIR}/abundance_${sample}.tsv"
                    done
                    printf 'sample_manifest_fingerprint=%s\n' "$(sha256sum "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}" 2>/dev/null | awk '{print $1}')"
                fi ;;
            12c_run_overview)
                printf 'schema_version=1\n'
                printf 'generator_version=%s\n' "$(sha256sum "${SUMMARY_GENERATOR}" 2>/dev/null | awk '{print $1}')"
                append_file_fingerprint resume_plan "${RESUME_PLAN_FILE}"
                append_file_fingerprint run_status "${RUN_STATUS_FILE}"
                append_file_fingerprint software_inventory "${LOGS_DIR}/software_versions.tsv"
                append_file_fingerprint database_inventory "${LOGS_DIR}/database_versions.tsv"
                printf 'profile=%s\nthreads=%s\n' "${REQUESTED_PROFILE}" "${THREADS}"
                printf 'taxonomy_summary_status=%s\n' "$([[ -s "${SUMMARY_DIR}/taxonomy_summary.tsv" ]] && echo available || echo "${STAGE_FINAL_STATUS[12a_taxonomy_summary]:-not-applicable}")"
                printf 'mag_summary_status=%s\n' "$([[ -s "${SUMMARY_DIR}/mag_summary.tsv" ]] && echo available || echo "${STAGE_FINAL_STATUS[12b_mag_summary]:-not-applicable}")" ;;
        esac
    }
}

validate_step_products() {
    local step="$1"
    local -a samples=()
    local -a paths=()
    local sample path bam
    local found_data=0

    case "${step}" in
        01_fastp)
            for sample in "${ACTIVE_SAMPLE_IDS[@]}"; do
                require_nonempty_file "${QC_DIR}/trimmed/${sample}${READ_SUFFIX_R1}" "clean R1 read" || return 1
                require_nonempty_file "${QC_DIR}/trimmed/${sample}${READ_SUFFIX_R2}" "matching clean R2 read" || return 1
            done ;;
        02_megahit)
            require_nonempty_file "${ASSEMBLY_DIR}/final.contigs.fa" "MEGAHIT assembly" ;;
        03_kraken2)
            mapfile -t samples < <(expected_sample_names)
            [[ ${#samples[@]} -gt 0 ]] || {
                log_error "[Validation] No samples available for Kraken2."
                return 1
            }
            for sample in "${samples[@]}"; do
                require_nonempty_file "${KRAKEN2_DIR}/kraken_output_${sample}.txt" "Kraken2 output for ${sample}" || return 1
                require_nonempty_file "${KRAKEN2_DIR}/kraken_report_${sample}.txt" "Kraken2 report for ${sample}" || return 1
                require_nonempty_file "${KRAKEN2_DIR}/bracken_${sample}.txt" "Bracken output for ${sample}" || return 1
            done ;;
        04_mapping)
            collect_files_by_pattern "${MAPPING_DIR}" "*_sorted.bam" paths
            mapfile -t samples < <(expected_sample_names)
            [[ ${#samples[@]} -gt 0 && ${#paths[@]} -gt 0 ]] || {
                log_error "[Validation] No sorted BAM files found in ${MAPPING_DIR}."
                return 1
            }
            for sample in "${samples[@]}"; do
                bam="${MAPPING_DIR}/${sample}_sorted.bam"
                require_nonempty_file "${bam}" "sorted BAM" || return 1
                require_nonempty_file "${bam}.bai" "BAM index" || return 1
            done
            require_nonempty_file "${BINNING_DIR}/depth_matrix.txt" "JGI depth matrix" ;;
        05a_vamb)
            validate_vamb_fasta "${CONTIGS}" || return 1
            validate_vamb_abundance_tsv "${VAMB_ABUNDANCE_TSV}" "${CONTIGS}" || return 1
            validate_bins_dir "${VAMB_FINAL_BINS_DIR}" "Vamb normalized bins" ;;
        05b_metabinner)
            require_nonempty_file "${METABINNER_DIR}/assembly_${MIN_CONTIG_LEN}.fa" "MetaBinner filtered assembly" || return 1
            require_nonempty_file "${METABINNER_DIR}/assembly_${MIN_CONTIG_LEN}_kmer_4_f${MIN_CONTIG_LEN}.csv" "MetaBinner k-mer table" || return 1
            validate_bins_dir "${METABINNER_FINAL_BINS_DIR}" "MetaBinner normalized bins" ;;
        05c_comebin)
            validate_bins_dir "${COMEBIN_FINAL_BINS_DIR}" "COMEBin normalized bins" ;;
        06_refinement)
            validate_bins_dir "${REFINED_BINS_DIR}" "refined bins" ;;
        07_reassembly)
            validate_bins_dir "${REASSEMBLY_DIR}/reassembled_bins" "reassembled bins" ;;
        08_checkm2)
            require_nonempty_file "${CHECKM2_DIR}/quality_report.tsv" "CheckM2 quality report" || return 1
            if ! awk 'NR > 1 && NF { found=1 } END { exit !found }' "${CHECKM2_DIR}/quality_report.tsv"; then
                log_error "[Validation] CheckM2 quality_report.tsv has no data rows."
                return 1
            fi ;;
        09_filter_mags)
            require_files_matching "${MAGS_DIR}" "*.fa" "retained MAG FASTA" &&
                require_nonempty_file "${MAG_FILTER_MANIFEST}" "MAG filtering manifest" ;;
        10_gtdbtk)
            for path in "${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv" "${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv"; do
                if [[ -s "${path}" ]] && awk 'NR > 1 && NF { found=1 } END { exit !found }' "${path}"; then
                    found_data=1
                fi
            done
            if [[ "${found_data}" -ne 1 ]]; then
                log_error "[Validation] GTDB-Tk produced no non-empty summary with data rows."
                return 1
            fi ;;
        11_coverm)
            mapfile -t samples < <(expected_sample_names)
            [[ ${#samples[@]} -gt 0 ]] || {
                log_error "[Validation] No samples available for CoverM."
                return 1
            }
            for sample in "${samples[@]}"; do
                require_nonempty_file "${ABUNDANCE_DIR}/abundance_${sample}.tsv" "CoverM abundance for ${sample}" || return 1
            done ;;
        12a_taxonomy_summary)
            require_nonempty_file "${SUMMARY_DIR}/taxonomy_summary.tsv" "taxonomy summary" || return 1
            require_nonempty_file "${SUMMARY_DIR}/taxonomy_summary.meta.json" "taxonomy metadata" || return 1
            if ! awk 'NR > 1 && NF { found=1 } END { exit !found }' "${SUMMARY_DIR}/taxonomy_summary.tsv"; then
                log_error "[Validation] taxonomy_summary.tsv has no data rows."
                return 1
            fi ;;
        12b_mag_summary)
            require_nonempty_file "${SUMMARY_DIR}/mag_summary.tsv" "MAG summary" || return 1
            require_nonempty_file "${SUMMARY_DIR}/mag_summary.meta.json" "MAG metadata" || return 1
            if ! awk 'NR > 1 && NF { found=1 } END { exit !found }' "${SUMMARY_DIR}/mag_summary.tsv"; then
                log_error "[Validation] mag_summary.tsv has no data rows."
                return 1
            fi ;;
        12c_run_overview)
            require_nonempty_file "${SUMMARY_DIR}/run_overview.tsv" "run overview" || return 1
            require_nonempty_file "${SUMMARY_DIR}/run_overview.meta.json" "run overview metadata" || return 1
            if ! awk 'NR > 1 && NF { found=1 } END { exit !found }' "${SUMMARY_DIR}/run_overview.tsv"; then
                log_error "[Validation] run_overview.tsv has no data rows."
                return 1
            fi ;;
        *)
            log_error "[Validation] Unknown pipeline step: ${step}"
            return 1 ;;
    esac
}

describe_step_validations() {
    case "$1" in
        01_fastp) log_info "[Validation] Would require clean R1 files and matching non-empty clean R2 files." ;;
        02_megahit) log_info "[Validation] Would require a non-empty final.contigs.fa." ;;
        03_kraken2) log_info "[Validation] Would require non-empty Kraken2 output/report and Bracken output for every sample." ;;
        04_mapping) log_info "[Validation] Would require non-empty sorted BAMs, their indexes, and the JGI depth matrix." ;;
        05a_vamb) log_info "[Validation] Would require valid normalized Vamb bins in ${VAMB_FINAL_BINS_DIR}." ;;
        05b_metabinner) log_info "[Validation] Would require valid normalized MetaBinner bins in ${METABINNER_FINAL_BINS_DIR}." ;;
        05c_comebin) log_info "[Validation] Would require valid normalized COMEBin bins in ${COMEBIN_FINAL_BINS_DIR}." ;;
        06_refinement) log_info "[Validation] Would require valid refined bins in ${REFINED_BINS_DIR}." ;;
        07_reassembly) log_info "[Validation] Would require at least one non-empty reassembled *.fa bin." ;;
        08_checkm2) log_info "[Validation] Would require a non-empty quality_report.tsv with at least one data row." ;;
        09_filter_mags) log_info "[Validation] Would require retained MAG FASTAs and ${MAG_FILTER_MANIFEST}." ;;
        10_gtdbtk) log_info "[Validation] Would require a non-empty GTDB-Tk summary containing a data row." ;;
        11_coverm) log_info "[Validation] Would require a non-empty abundance TSV for every sample." ;;
        12a_taxonomy_summary) log_info "[Validation] Would require a non-empty taxonomy_summary.tsv and metadata." ;;
        12b_mag_summary) log_info "[Validation] Would require a non-empty mag_summary.tsv and metadata." ;;
        12c_run_overview) log_info "[Validation] Would require a non-empty run_overview.tsv and metadata." ;;
        *) log_error "[Validation] Unknown pipeline step: $1"; return 1 ;;
    esac
}

_checkpoint_dependents() {
    get_direct_dependents "$1"
}

checkpoint_stage_enabled() {
    case "$1" in
        01_fastp) [[ "${RUN_FASTP}" == true ]] ;;
        02_megahit) [[ "${RUN_MEGAHIT}" == true ]] ;;
        03_kraken2) [[ "${RUN_KRAKEN2}" == true ]] ;;
        04_mapping) [[ "${RUN_MAPPING}" == true ]] ;;
        05a_vamb) [[ "${RUN_BINNING}" == true ]] && binner_selected vamb ;;
        05b_metabinner) [[ "${RUN_BINNING}" == true ]] && binner_selected metabinner ;;
        05c_comebin) [[ "${RUN_BINNING}" == true ]] && binner_selected comebin ;;
        06_refinement) [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] ;;
        07_reassembly) [[ "${RUN_REASSEMBLY}" == true ]] ;;
        08_checkm2) [[ "${RUN_CHECKM2}" == true ]] ;;
        09_filter_mags) [[ "${RUN_MAG_FILTERING}" == true ]] ;;
        10_gtdbtk) [[ "${RUN_GTDBTK}" == true ]] ;;
        11_coverm) [[ "${RUN_COVERM}" == true ]] ;;
        12a_taxonomy_summary) [[ "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true && "${RUN_KRAKEN2}" == true ]] ;;
        12b_mag_summary) [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true && ( "${RUN_CHECKM2}" == true || "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || "${RUN_COVERM}" == true ) ]] ;;
        12c_run_overview) [[ "${RUN_SUMMARY}" == true && "${RUN_RUN_OVERVIEW}" == true ]] ;;
        *) return 1 ;;
    esac
}

stage_action_requires_execution() {
    local action="${STAGE_ACTION[$1]:-disabled}"
    [[ "${action}" == run || "${action}" == force || "${action}" == restart ]]
}

invalidate_checkpoint() {
    local step="$1"
    local reason="$2"
    local item
    local -a affected=("${step}")
    local -a _descendants=()
    mapfile -t _descendants < <(get_transitive_dependents "${step}")
    affected+=("${_descendants[@]}")
    local -A invalidated=()
    for item in "${affected[@]}"; do
        stage_is_active "${item}" || continue
        [[ -n "${invalidated[$item]:-}" ]] && continue
        invalidated["${item}"]=1
        if [[ "${DRY_RUN}" == true ]]; then
            log_warn "[Checkpoint] Dry-run: would invalidate ${item} because: ${reason}"
        else
            rm -f -- "${CHECKPOINT_DIR}/${item}.done" "${CHECKPOINT_DIR}/${item}.config"
            log_warn "[Checkpoint] Invalidated ${item}: ${reason}"
        fi
    done
}

is_done() {
    local step="$1"
    local current_config recorded_config comparable_current comparable_recorded
    local action="${STAGE_ACTION[$step]:-run}"
    if [[ "${action}" != reuse ]]; then
        if [[ -z "${PREPARED_STAGES[$step]:-}" ]]; then
            prepare_stage_execution "${step}"
            PREPARED_STAGES["${step}"]=1
        fi
        return 1
    fi
    assert_runtime_identity_frozen "${step}" || exit 1
    assert_database_identity_frozen "${step}" || exit 1
    if [[ -z "${STAGE_FROZEN_EXTERNAL[$step]:-}" || \
          "$(stage_external_fingerprint "${step}")" != "${STAGE_FROZEN_EXTERNAL[$step]}" ]]; then
        log_error "External inputs changed after the resume plan was frozen for ${step}; aborting without replanning."
        exit 1
    fi
    [[ -f "${CHECKPOINT_DIR}/${step}.done" && -f "${CHECKPOINT_DIR}/${step}.config" ]] || {
        log_error "Checkpoint for planned reuse stage ${step} disappeared after planning."
        exit 1
    }

    current_config=$(checkpoint_config_manifest "${step}")
    recorded_config="$(cat "${CHECKPOINT_DIR}/${step}.config")"
    comparable_current="$(comparable_checkpoint_manifest "${step}" "${current_config}")"
    comparable_recorded="$(comparable_checkpoint_manifest "${step}" "${recorded_config}")"
    if [[ "${comparable_recorded}" != "${comparable_current}" ]] || \
       ! checkpoint_identity_compatible "${step}" "${current_config}" "${recorded_config}"; then
        log_error "Checkpoint for ${step} no longer matches the frozen resume plan."
        exit 1
    fi
    if ! validate_step_products "${step}"; then
        log_error "Checkpoint products changed after the resume plan was frozen for ${step}."
        exit 1
    fi
    log_info "[Checkpoint] Step ${step} configuration and products validated."
    STAGE_FINAL_STATUS["${step}"]=reused
    update_run_status "${step}" reused
    return 0
}

comparable_checkpoint_manifest() {
    local step="$1" manifest="$2"
    printf '%s\n' "${manifest}" | sed \
        -e '/^effective_tool_version=/d' \
        -e '/^tool_version_source=/d' \
        -e '/^execution_mode=/d' \
        -e '/^runtime_name=/d' \
        -e '/^runtime_uri=/d' \
        -e '/^runtime_path=/d' \
        -e '/^runtime_size_bytes=/d' \
        -e '/^runtime_mtime_ns=/d' \
        -e '/^runtime_fingerprint_method=/d' \
        -e '/^runtime_fingerprint=/d' \
        -e '/^database_release_or_version=/d' \
        -e '/^database_version_source=/d' \
        -e '/^database_path=/d' \
        -e '/^database_fingerprint_method=/d' \
        -e '/^database_fingerprint=/d' \
        -e '/^resource_[^=]*=/d' \
        -e '/^image=/d' \
        -e '/^image_uri=/d' \
        -e '/^[[:alnum:]_]*_image=/d' \
        -e '/^vamb_version=/d' \
        -e '/^vamb_runtime=/d' \
        -e '/^comebin_version=/d' \
        -e '/^runtime=/d' \
        -e '/^mamba_environment=/d' \
        -e '/^db=/d' \
        -e '/^data=/d' \
        -e '/^kraken_db=/d' \
        -e '/^checkm2_db=/d' \
        -e '/^gtdbtk_data=/d' \
        -e '/^refinement_mode=/d' \
        -e '/^selected_binners=/d' \
        -e '/^bins_source_requested=/d'
}

checkpoint_identity_compatible() {
    local step="$1" current="$2" recorded="$3" tool key previous_version previous_method previous_fingerprint
    local current_release previous_release current_method current_fingerprint previous_mode
    local -a tools=()
    if ! grep -q '^resource_[^=]*_version=' <<< "${recorded}"; then
        log_warn "[Checkpoint] ${step} predates functional identity metadata; reusing under the non-destructive migration policy."
        return 0
    fi
    previous_mode="$(awk -F= '$1=="execution_mode" {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
    [[ -z "${previous_mode}" || "${previous_mode}" == "${EXEC_MODE}" ]] || return 1
    mapfile -t tools < <(inventory_stage_tools "${step}")
    for tool in "${tools[@]}"; do
        key="${step}|${tool}"
        previous_version="$(awk -F= -v key="resource_${tool}_version" '$1==key {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
        previous_method="$(awk -F= -v key="resource_${tool}_runtime_fingerprint_method" '$1==key {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
        previous_fingerprint="$(awk -F= -v key="resource_${tool}_runtime_fingerprint" '$1==key {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
        if [[ "${ID_STATUS[$key]:-unknown}" == unavailable ]]; then
            if [[ "${previous_fingerprint:-unknown}" == unknown && "${previous_version:-unknown}" == unknown ]]; then
                log_warn "[Checkpoint] ${step}/${tool} checkpoint has neither known version nor fingerprint."
                return 1
            fi
            continue
        fi
        current_method="${ID_RUNTIME_METHOD[$key]:-unknown}"
        current_fingerprint="${ID_RUNTIME_FINGERPRINT[$key]:-unknown}"
        if [[ -n "${previous_fingerprint}" && "${previous_fingerprint}" != unknown && \
              "${current_fingerprint}" != unknown ]]; then
            [[ "${previous_method}" == "${current_method}" && \
               "${previous_fingerprint}" == "${current_fingerprint}" ]] || return 1
        elif [[ -n "${previous_version}" && "${previous_version}" != unknown && \
                "${ID_VERSION[$key]:-unknown}" != unknown ]]; then
            [[ "${previous_version}" == "${ID_VERSION[$key]}" ]] || return 1
        else
            log_warn "[Checkpoint] ${step}/${tool} has no comparable known runtime identity."
            return 1
        fi
        if [[ -n "${previous_version}" && "${previous_version}" != unknown && \
              "${ID_VERSION[$key]:-unknown}" != unknown && \
              "${previous_version}" != "${ID_VERSION[$key]}" ]]; then
            return 1
        fi
    done
    case "${step}" in
        03_kraken2) key="${step}|kraken2-bracken" ;;
        08_checkm2) key="${step}|checkm2" ;;
        10_gtdbtk) key="${step}|gtdbtk" ;;
        *) return 0 ;;
    esac
    previous_fingerprint="$(awk -F= '$1=="database_fingerprint" {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
    previous_method="$(awk -F= '$1=="database_fingerprint_method" {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
    previous_release="$(awk -F= '$1=="database_release_or_version" {sub(/^[^=]*=/,""); print; exit}' <<< "${recorded}")"
    if [[ "${ID_DATABASE_STATUS[$key]:-unknown}" == unavailable ]]; then
        if [[ "${previous_fingerprint:-unknown}" == unknown && "${previous_release:-unknown}" == unknown ]]; then
            log_warn "[Checkpoint] ${step} database checkpoint has neither known release nor fingerprint."
            return 1
        fi
        return 0
    fi
    current_fingerprint="${ID_DATABASE_FINGERPRINT[$key]:-unknown}"
    current_method="${ID_DATABASE_METHOD[$key]:-unknown}"
    current_release="${ID_DATABASE_RELEASE[$key]:-unknown}"
    if [[ -n "${previous_fingerprint}" && "${previous_fingerprint}" != unknown && \
          "${current_fingerprint}" != unknown ]]; then
        [[ "${previous_method}" == "${current_method}" && \
           "${previous_fingerprint}" == "${current_fingerprint}" ]] || return 1
    elif [[ -n "${previous_release}" && "${previous_release}" != unknown && \
            "${current_release}" != unknown ]]; then
        [[ "${previous_release}" == "${current_release}" ]] || return 1
    else
        log_warn "[Checkpoint] ${step} has no comparable known database identity."
        return 1
    fi
    [[ "${previous_release}" == unknown || "${current_release}" == unknown || \
       "${previous_release}" == "${current_release}" ]]
}

stage_external_fingerprint() {
    local step="$1" path file sample
    local -a paths=()
    case "${step}" in
        01_fastp)
            if [[ -n "${RAW_READS_DIR}" && -d "${RAW_READS_DIR}" ]]; then
                paths+=("${RAW_READS_DIR}")
                while IFS= read -r -d '' path; do paths+=("${path}"); done \
                    < <(find "${RAW_READS_DIR}" -maxdepth 1 -type f \
                        \( -name "*${READ_SUFFIX_R1}" -o -name "*${READ_SUFFIX_R2}" \) -print0 | sort -z)
            fi ;;
        03_kraken2|02_megahit|04_mapping|07_reassembly|11_coverm)
            if [[ "${ACTIVE_READS_SOURCE}" == raw || "${ACTIVE_READS_SOURCE}" == external-clean ]]; then
                paths+=("${ACTIVE_READS_DIR}")
                paths+=("${ACTIVE_R1_FILES[@]}" "${ACTIVE_R2_FILES[@]}")
            fi
            [[ -n "${ASSEMBLY_IN}" && ( "${step}" == 04_mapping || \
                "${step}" == 05a_vamb || "${step}" == 05b_metabinner || "${step}" == 05c_comebin ) ]] && \
                paths+=("${ASSEMBLY_IN}")
            if [[ "${step}" == 05a_vamb || "${step}" == 05b_metabinner ]]; then
                [[ "${COVERAGE_SOURCE}" == external && -n "${JGI_DEPTH_FILE_ACTIVE}" ]] && \
                    paths+=("${JGI_DEPTH_FILE_ACTIVE}")
            fi
            if [[ "${step}" == 05c_comebin && "${COVERAGE_SOURCE}" == external ]]; then
                paths+=("${ACTIVE_BAMS[@]}" "${ACTIVE_BAM_INDEXES[@]}")
                [[ -n "${BAM_DIR_ACTIVE}" ]] && paths+=("${BAM_DIR_ACTIVE}")
                [[ -n "${BAM_MANIFEST_INPUT}" ]] && paths+=("${BAM_MANIFEST_ACTIVE}")
            fi ;;
        05a_vamb|05b_metabinner)
            [[ -n "${ASSEMBLY_IN}" ]] && paths+=("${ASSEMBLY_IN}")
            [[ "${COVERAGE_SOURCE}" == external && -n "${JGI_DEPTH_FILE_ACTIVE}" ]] && \
                paths+=("${JGI_DEPTH_FILE_ACTIVE}") ;;
        05c_comebin)
            [[ -n "${ASSEMBLY_IN}" ]] && paths+=("${ASSEMBLY_IN}")
            if [[ "${COVERAGE_SOURCE}" == external ]]; then
                paths+=("${ACTIVE_BAMS[@]}" "${ACTIVE_BAM_INDEXES[@]}")
                [[ -n "${BAM_DIR_ACTIVE}" ]] && paths+=("${BAM_DIR_ACTIVE}")
                [[ -n "${BAM_MANIFEST_INPUT}" ]] && paths+=("${BAM_MANIFEST_ACTIVE}")
            fi ;;
    esac
    if [[ "${step}" == 04_mapping && -n "${CONTIGS}" ]]; then
        paths+=("${CONTIGS}")
    fi
    for path in "${paths[@]}"; do
        if [[ -e "${path}" ]]; then
            file="$(realpath "${path}")"
            if [[ -f "${file}" ]]; then
                printf '%s|%s|%s\n' "${file}" "$(stat -c '%s' "${file}")" "$(stat -c '%y' "${file}")"
            else
                printf '%s|dir|%s\n' "${file}" "$(stat -c '%y' "${file}")"
            fi
        fi
    done | sort
}

prepare_stage_execution() {
    local step="$1" frozen="${STAGE_FROZEN_EXTERNAL[$1]:-}"
    local current
    current="$(stage_external_fingerprint "${step}")"
    if [[ "${current}" != "${frozen}" ]]; then
        log_error "External inputs changed after planning for ${step}; refusing to execute without replanning."
        exit 1
    fi
    [[ -n "${PREPARED_STAGES[$step]:-}" ]] && return 0
    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[Plan] ${STAGE_ACTION[$step]} ${step}: remove stale checkpoint and clean stage outputs immediately before execution."
        clean_stage_outputs "${step}"
        return 0
    fi
    rm -f -- "${CHECKPOINT_DIR}/${step}.done" "${CHECKPOINT_DIR}/${step}.config"
    clean_stage_outputs "${step}"
}

apply_restart_cleanup() {
    local id current
    for id in "${STAGE_IDS[@]}"; do
        [[ "${STAGE_ACTION[$id]}" == restart ]] || continue
        current="$(stage_external_fingerprint "${id}")"
        if [[ "${current}" != "${STAGE_FROZEN_EXTERNAL[$id]:-}" ]]; then
            log_error "External inputs changed after planning for ${id}; refusing restart cleanup."
            exit 1
        fi
    done
    for id in "${STAGE_IDS[@]}"; do
        [[ "${STAGE_ACTION[$id]}" == restart ]] || continue
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Plan] Restart would invalidate and clean ${id}."
            clean_stage_outputs "${id}"
            continue
        fi
        rm -f -- "${CHECKPOINT_DIR}/${id}.done" "${CHECKPOINT_DIR}/${id}.config"
        clean_stage_outputs "${id}"
        PREPARED_STAGES["${id}"]=1
        STAGE_FINAL_STATUS["${id}"]=not-reached
        update_run_status "${id}" not-reached
    done
}

stage_display_name() {
    local id="$1"
    case "${id}" in
        01_fastp) echo "Quality Control (fastp)" ;;
        02_megahit) echo "Assembly (MEGAHIT)" ;;
        03_kraken2) echo "Taxonomy (Kraken2/Bracken)" ;;
        04_mapping) echo "Read Mapping & Contig Depth" ;;
        05a_vamb) echo "Binning (Vamb)" ;;
        05b_metabinner) echo "Binning (MetaBinner)" ;;
        05c_comebin) echo "Binning (COMEBin)" ;;
        06_refinement) echo "Bin Refinement (MetaWRAP)" ;;
        07_reassembly) echo "Bin Reassembly" ;;
        08_checkm2) echo "Quality Assessment (CheckM2)" ;;
        09_filter_mags) echo "MAG Filtering" ;;
        10_gtdbtk) echo "Taxonomic Classification (GTDB-Tk)" ;;
        11_coverm) echo "MAG Coverage (CoverM)" ;;
        12a_taxonomy_summary) echo "Taxonomy Summary" ;;
        12b_mag_summary) echo "MAG Summary" ;;
        12c_run_overview) echo "Run Overview" ;;
        12_summary) echo "Consolidated Summary" ;;
        pipeline) echo "Preflight" ;;
        *) echo "${id}" ;;
    esac
}

stage_primary_log() {
    local id="$1"
    local log="unknown"
    if [[ -n "${LOGS_DIR:-}" ]]; then
        case "${id}" in
            01_fastp) log="${LOGS_DIR}/step_01_fastp.log" ;;
            02_megahit) log="${LOGS_DIR}/step_02_megahit.log" ;;
            03_kraken2) log="${LOGS_DIR}/step_03_kraken2.log" ;;
            04_mapping) log="${LOGS_DIR}/step_04_mapping.log" ;;
            05a_vamb) log="${LOGS_DIR}/step_05a_vamb.log" ;;
            05b_metabinner) log="${LOGS_DIR}/step_05b_metabinner_run.log" ;;
            05c_comebin) log="${LOGS_DIR}/step_05c_comebin.log" ;;
            06_refinement) log="${LOGS_DIR}/step_06_bin_refinement.log" ;;
            07_reassembly) log="${LOGS_DIR}/step_07_reassembly.log" ;;
            08_checkm2) log="${LOGS_DIR}/step_08_checkm2.log" ;;
            09_filter_mags) log="${LOGS_DIR}/step_09_filter_mags.log" ;;
            10_gtdbtk) log="${LOGS_DIR}/step_10_gtdbtk.log" ;;
            11_coverm) log="${LOGS_DIR}/step_11_coverm.log" ;;
            12a_taxonomy_summary) log="${LOGS_DIR}/step_12a_taxonomy_summary.log" ;;
            12b_mag_summary) log="${LOGS_DIR}/step_12b_mag_summary.log" ;;
            12c_run_overview) log="${LOGS_DIR}/step_12c_run_overview.log" ;;
            pipeline) log="${LOGS_DIR}/pipeline.log" ;;
            *) log="unknown" ;;
        esac
    fi
    echo "${log}"
}

stage_primary_output() {
    local id="$1"
    case "${id}" in
        01_fastp) echo "${FASTP_CLEAN_SAMPLE_MANIFEST:-${QC_DIR}/trimmed}" ;;
        02_megahit) echo "${CONTIGS:-${ASSEMBLY_DIR}/final.contigs.fa}" ;;
        03_kraken2) echo "${KRAKEN2_DIR}/reports" ;;
        04_mapping) echo "${DEPTH_FILE:-${BINNING_DIR}/depth_matrix.txt}" ;;
        05a_vamb) echo "${VAMB_FINAL_BINS_DIR}" ;;
        05b_metabinner) echo "${METABINNER_FINAL_BINS_DIR}" ;;
        05c_comebin) echo "${COMEBIN_FINAL_BINS_DIR}" ;;
        06_refinement) echo "${REFINED_BINS_DIR:-${REFINEMENT_DIR}}" ;;
        07_reassembly) echo "${REASSEMBLY_DIR}" ;;
        08_checkm2) echo "${CHECKM2_DIR}/quality_report.tsv" ;;
        09_filter_mags) echo "${MAGS_DIR}" ;;
        10_gtdbtk) echo "${GTDBTK_DIR}" ;;
        11_coverm) echo "${ABUNDANCE_DIR}" ;;
        12a_taxonomy_summary) echo "${SUMMARY_DIR}/taxonomy_summary.tsv" ;;
        12b_mag_summary) echo "${SUMMARY_DIR}/mag_summary.tsv" ;;
        12c_run_overview) echo "${SUMMARY_DIR}/run_overview.tsv" ;;
        *) echo "NA" ;;
    esac
}

stage_primary_output_fingerprint() {
    local out
    out="$(stage_primary_output "$1")"
    if [[ -f "${out}" ]]; then
        sha256sum "${out}" 2>/dev/null | awk '{print $1}' || echo "NA"
    elif [[ -d "${out}" ]]; then
        find "${out}" -maxdepth 2 -type f -exec stat -c '%n %s' {} + 2>/dev/null | sort | sha256sum 2>/dev/null | awk '{print $1}' || echo "NA"
    else
        echo "NA"
    fi
}

stage_primary_tool_version() {
    local stage="$1" tool
    local -a tools=()
    mapfile -t tools < <(inventory_stage_tools "${stage}" 2>/dev/null || true)
    for tool in "${tools[@]}"; do
        local v="${ID_VERSION[${stage}|${tool}]:-unknown}"
        [[ "${v}" != "unknown" && -n "${v}" ]] && { echo "${v}"; return 0; }
    done
    echo "NA"
}

stage_primary_runtime_fingerprint() {
    local stage="$1" tool
    local -a tools=()
    mapfile -t tools < <(inventory_stage_tools "${stage}" 2>/dev/null || true)
    for tool in "${tools[@]}"; do
        local fp="${ID_RUNTIME_FINGERPRINT[${stage}|${tool}]:-unknown}"
        [[ "${fp}" != "unknown" && -n "${fp}" ]] && { echo "${fp}"; return 0; }
    done
    echo "NA"
}

stage_primary_database_fingerprint() {
    local stage="$1" dbkey=""
    case "${stage}" in
        03_kraken2) dbkey="${stage}|kraken2-bracken" ;;
        08_checkm2) dbkey="${stage}|checkm2" ;;
        10_gtdbtk) dbkey="${stage}|gtdbtk" ;;
        *) echo "NA"; return 0 ;;
    esac
    local dfp="${DB_FINGERPRINT[${dbkey}]:-unknown}"
    if [[ "${dfp}" != "unknown" && -n "${dfp}" ]]; then
        echo "${dfp}"
    else
        echo "NA"
    fi
}

set_last_valid_checkpoint() {
    local chk_id="$1"
    local chk_path="${2:-none}"
    local prod="${3:-none}"
    local prod_fp="${4:-none}"
    local ts
    ts="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "unknown-time")"

    LAST_VALID_CHECKPOINT="${chk_id}"
    LAST_VALID_CHECKPOINT_PATH="${chk_path}"
    LAST_VALID_PRODUCT="${prod}"
    LAST_VALID_PRODUCT_FINGERPRINT="${prod_fp}"

    [[ "${DRY_RUN:-false}" == true ]] && return 0
    [[ -z "${OUTDIR:-}" || ! -d "${OUTDIR}" ]] && return 0

    local chk_state_dir="${LOGS_DIR}/checkpoint_state"
    mkdir -p "${chk_state_dir}"
    local run_chk_file="${chk_state_dir}/${RUN_ID}.tsv"
    local last_state_file="${LOGS_DIR}/last_valid_checkpoint.state"

    if [[ ! -s "${run_chk_file}" ]]; then
        printf 'run_id\tcheckpoint_id\tcheckpoint_path\tproduct\tproduct_fingerprint\ttimestamp_utc\n' > "${run_chk_file}"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${RUN_ID}" "${chk_id}" "${chk_path}" "${prod}" "${prod_fp}" "${ts}" >> "${run_chk_file}"

    local tmp_state
    tmp_state="$(mktemp "${LOGS_DIR}/.last_chk_state.XXXXXX")"
    printf 'run_id=%s\ncheckpoint_id=%s\ncheckpoint_path=%s\nproduct=%s\nproduct_fingerprint=%s\ntimestamp_utc=%s\n' \
        "${RUN_ID}" "${chk_id}" "${chk_path}" "${prod}" "${prod_fp}" "${ts}" > "${tmp_state}"
    mv -f "${tmp_state}" "${last_state_file}"
}

set_last_valid_product() {
    local prod="$1"
    local prod_fp="${2:-none}"
    LAST_VALID_PRODUCT="${prod}"
    LAST_VALID_PRODUCT_FINGERPRINT="${prod_fp}"
}

run_stage() {
    local id="$1" function_name="$2" action="${STAGE_ACTION[$1]:-disabled}"
    stage_is_active "${id}" || return 0
    if [[ "${PLAN_IS_FROZEN}" == true ]]; then
        verify_stage_preconditions "${id}" || return 1
    fi
    if stage_action_requires_execution "${id}"; then
        assert_runtime_identity_frozen "${id}" || return 1
        assert_database_identity_frozen "${id}" || return 1
    fi
    CURRENT_STAGE="${id}"
    CURRENT_STAGE_NAME="$(stage_display_name "${id}")"
    CURRENT_STAGE_ACTION="${action}"
    CURRENT_STAGE_LOG="$(stage_primary_log "${id}")"
    CURRENT_STAGE_RUNTIME="unknown"
    CURRENT_STAGE_COMMAND=""
    CURRENT_COMMAND_SOURCE_FILE=""
    CURRENT_COMMAND_SOURCE_LINE=""

    STAGE_START_TIME["${id}"]="$(date -Is)"
    STAGE_FINAL_STATUS["${id}"]=not-reached

    if [[ "${action}" == reuse ]]; then
        if ! is_done "${id}"; then
            log_error "[Step ${id}] Action is reuse but checkpoint is not valid or complete."
            return 1
        fi
        if ! validate_step_products "${id}"; then
            log_error "[Step ${id}] Action is reuse but products failed validation."
            return 1
        fi
        local chk_path="${CHECKPOINT_DIR}/${id}.done"
        set_last_valid_checkpoint "${id}" "${chk_path}" "reused" "none"
        STAGE_START_TIME["${id}"]="NA"
        STAGE_END_TIME["${id}"]="NA"
        STAGE_EXIT_CODE["${id}"]=0
        STAGE_FINAL_STATUS["${id}"]=reused
        if [[ "${DRY_RUN}" != true && "${RUN_STATUS_INITIALIZED:-false}" == true ]]; then
            python3 "${BENCHMARK_WRITER}" stage-finish \
                --stage-status-path "${STAGE_STATUS_FILE}" \
                --compat-run-status-path "${RUN_STATUS_FILE}" \
                --run-id "${RUN_ID}" \
                --stage-id "${id}" \
                --planned-action "reuse" \
                --checkpoint-path "${chk_path}" \
                --checkpoint-reason "${STAGE_REASON[$id]:-checkpoint_reused}" \
                --validation-status "passed" \
                --primary-log "${CURRENT_STAGE_LOG}" \
                --primary-output "$(stage_primary_output "${id}")" \
                --output-fingerprint "$(stage_primary_output_fingerprint "${id}")" \
                --runtime "${CURRENT_STAGE_RUNTIME:-apptainer}" \
                --message "reused from valid checkpoint" >/dev/null 2>&1 || true
        fi
    else
        local stage_start_utc stage_end_utc stage_start_mono stage_end_mono wall_seconds="NA"
        local stage_metrics_dir="${OUTDIR}/logs/.metrics_${id}_${RUN_ID}"
        CURRENT_STAGE_METRICS_DIR="${stage_metrics_dir}"

        stage_start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        stage_start_mono="$(python3 -c 'import time; print(f"{time.monotonic():.6f}")' 2>/dev/null || echo "0.0")"
        STAGE_START_TIME["${id}"]="${stage_start_utc}"

        if [[ "${DRY_RUN}" != true && "${RUN_STATUS_INITIALIZED:-false}" == true ]]; then
            mkdir -p "${stage_metrics_dir}"
            python3 "${BENCHMARK_WRITER}" stage-start \
                --stage-status-path "${STAGE_STATUS_FILE}" \
                --compat-run-status-path "${RUN_STATUS_FILE}" \
                --run-id "${RUN_ID}" \
                --stage-id "${id}" \
                --start-time-utc "${stage_start_utc}" >/dev/null 2>&1 || true
        fi

        # Immediate stage precheck for stages in DISK_PRECHECK_STAGES
        local is_precheck_stg=false stg_elem
        for stg_elem in "${DISK_PRECHECK_STAGES[@]}"; do
            if [[ "${stg_elem}" == "${id}" ]]; then is_precheck_stg=true; break; fi
        done
        if [[ "${is_precheck_stg}" == true ]]; then
            if ! run_stage_precheck "${id}"; then
                log_error "[Stage Precheck] Halting before execution of ${id} due to insufficient disk space."
                return 1
            fi
        fi

        "$function_name"
        local func_rc=$?

        # Record post-stage observation
        record_post_stage_observation "${id}"

        stage_end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        stage_end_mono="$(python3 -c 'import time; print(f"{time.monotonic():.6f}")' 2>/dev/null || echo "0.0")"
        wall_seconds="$(python3 -c "print(f'{(float(\"${stage_end_mono}\") - float(\"${stage_start_mono}\")):.3f}')" 2>/dev/null || echo "NA")"

        STAGE_END_TIME["${id}"]="${stage_end_utc}"
        STAGE_EXIT_CODE["${id}"]="${func_rc}"
        if [[ "${DRY_RUN}" == true ]]; then
            STAGE_FINAL_STATUS["${id}"]=dry-run-success
        else
            STAGE_FINAL_STATUS["${id}"]=success
        fi

        local chk_act="executed"
        [[ "${action}" == force ]] && chk_act="forced"
        [[ "${action}" == restart ]] && chk_act="restarted"

        if [[ "${DRY_RUN}" != true && "${RUN_STATUS_INITIALIZED:-false}" == true ]]; then
            local stg_finish_msg="stage completed successfully"
            if [[ -n "${CURRENT_STAGE_DISK_STATUS:-}" ]]; then
                stg_finish_msg="${stg_finish_msg};disk_check=${CURRENT_STAGE_DISK_STATUS};disk_required_bytes=${CURRENT_STAGE_DISK_REQ:-NA};disk_available_bytes=${CURRENT_STAGE_DISK_AVAIL:-NA}"
            fi
            python3 "${BENCHMARK_WRITER}" stage-finish \
                --stage-status-path "${STAGE_STATUS_FILE}" \
                --compat-run-status-path "${RUN_STATUS_FILE}" \
                --run-id "${RUN_ID}" \
                --stage-id "${id}" \
                --planned-action "${action}" \
                --checkpoint-action "${chk_act}" \
                --final-status "${STAGE_FINAL_STATUS[$id]}" \
                --start-time-utc "${stage_start_utc}" \
                --end-time-utc "${stage_end_utc}" \
                --wall-time-seconds "${wall_seconds}" \
                --metrics-dir "${stage_metrics_dir}" \
                --exit-code "${func_rc}" \
                --checkpoint-path "${CHECKPOINT_DIR}/${id}.done" \
                --checkpoint-status "$([[ "${STAGE_CHECKPOINT_WRITTEN[$id]:-false}" == true ]] && echo written || echo missing)" \
                --checkpoint-reason "${STAGE_REASON[$id]:-run_executed}" \
                --validation-status "passed" \
                --primary-log "${CURRENT_STAGE_LOG}" \
                --primary-output "$(stage_primary_output "${id}")" \
                --output-fingerprint "$(stage_primary_output_fingerprint "${id}")" \
                --runtime "${CURRENT_STAGE_RUNTIME:-apptainer}" \
                --tool-version "$(stage_primary_tool_version "${id}")" \
                --runtime-fingerprint "$(stage_primary_runtime_fingerprint "${id}")" \
                --database-fingerprint "$(stage_primary_database_fingerprint "${id}")" \
                --message "${stg_finish_msg}" >/dev/null 2>&1 || true
        fi
        rm -rf "${stage_metrics_dir}" 2>/dev/null || true
    fi

    CURRENT_STAGE=""
    CURRENT_STAGE_NAME=""
    CURRENT_STAGE_ACTION=""
    CURRENT_STAGE_LOG=""
    CURRENT_STAGE_RUNTIME="unknown"
    CURRENT_STAGE_COMMAND=""
    CURRENT_COMMAND_SOURCE_FILE=""
    CURRENT_COMMAND_SOURCE_LINE=""
}

record_pipeline_termination_state() {
    local status="$1"
    local exit_code="$2"
    local error_id="$3"

    local current_stg="${CURRENT_STAGE:-pipeline}"
    local end_ts
    end_ts="$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo "")"

    STAGE_FINAL_STATUS["${current_stg}"]="${status}"
    STAGE_END_TIME["${current_stg}"]="${end_ts}"
    STAGE_EXIT_CODE["${current_stg}"]="${exit_code}"

    local id
    for id in "${STAGE_IDS[@]+"${STAGE_IDS[@]}"}"; do
        if [[ "${STAGE_FINAL_STATUS[$id]:-}" == not-reached || -z "${STAGE_FINAL_STATUS[$id]:-}" ]]; then
            STAGE_FINAL_STATUS["${id}"]=not-reached
        fi
    done

    if [[ "${DRY_RUN:-false}" == true || "${RUN_STATUS_INITIALIZED:-false}" != true || ( ! -f "${RUN_STATUS_FILE:-}" && ! -f "${STAGE_STATUS_FILE:-}" ) ]]; then
        return 2
    fi

    if [[ -f "${RUN_STATUS_FILE:-}" ]]; then
        python3 - "${RUN_STATUS_FILE}" "${current_stg}" "${status}" "${end_ts}" "${exit_code}" "${error_id}" "${RUN_ID}" "${LAST_VALID_CHECKPOINT:-none}" "${OUTDIR:-.}/logs/last_error.json" "${CURRENT_STAGE_ACTION:-unknown}" <<'PY_STATUS_TERM'
import csv, sys
path, stage, status, end_time, exit_code, error_id, run_id, last_chk, err_rep, stg_action = sys.argv[1:11]
try:
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)
    if not rows:
        sys.exit(1)
    headers = rows[0]
    extra_cols = ["run_id", "error_id", "last_valid_checkpoint", "error_report", "failed_stage", "failed_stage_action"]
    for col in extra_cols:
        if col not in headers:
            headers.append(col)

    stg_idx = headers.index("stage_internal_id") if "stage_internal_id" in headers else 0
    status_idx = headers.index("final_status") if "final_status" in headers else 3
    end_idx = headers.index("end_time") if "end_time" in headers else 5
    code_idx = headers.index("exit_code") if "exit_code" in headers else 6

    out_rows = [headers]
    found = False
    for row in rows[1:]:
        if not row:
            continue
        while len(row) < len(headers):
            row.append("")
        stg_id = row[stg_idx]
        if stg_id == stage:
            found = True
            row[status_idx] = status
            row[end_idx] = end_time
            row[code_idx] = str(exit_code)
        elif row[status_idx] == "not-reached" or not row[status_idx]:
            row[status_idx] = "not-reached"

        if "run_id" in headers: row[headers.index("run_id")] = run_id
        if "error_id" in headers: row[headers.index("error_id")] = error_id
        if "last_valid_checkpoint" in headers: row[headers.index("last_valid_checkpoint")] = last_chk
        if "error_report" in headers: row[headers.index("error_report")] = err_rep
        if "failed_stage" in headers: row[headers.index("failed_stage")] = stage
        if "failed_stage_action" in headers: row[headers.index("failed_stage_action")] = stg_action
        out_rows.append(row)

    if not found and stage == "pipeline":
        new_row = ["" for _ in headers]
        new_row[stg_idx] = "pipeline"
        if "stage_public_name" in headers: new_row[headers.index("stage_public_name")] = "Preflight"
        if "planned_action" in headers: new_row[headers.index("planned_action")] = stg_action
        new_row[status_idx] = status
        new_row[end_idx] = end_time
        new_row[code_idx] = str(exit_code)
        if "run_id" in headers: new_row[headers.index("run_id")] = run_id
        if "error_id" in headers: new_row[headers.index("error_id")] = error_id
        if "last_valid_checkpoint" in headers: new_row[headers.index("last_valid_checkpoint")] = last_chk
        if "error_report" in headers: new_row[headers.index("error_report")] = err_rep
        if "failed_stage" in headers: new_row[headers.index("failed_stage")] = stage
        if "failed_stage_action" in headers: new_row[headers.index("failed_stage_action")] = stg_action
        out_rows.append(new_row)

    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(out_rows)
except Exception:
    sys.exit(1)
PY_STATUS_TERM
    fi

    local b_writer="${BENCHMARK_WRITER:-${PIPELINE_SCRIPTS_DIR:-./scripts}/write_benchmark_manifest.py}"
    if [[ -f "${STAGE_STATUS_FILE:-}" && -f "${b_writer}" ]]; then
        local end_mono
        end_mono="$(python3 -c 'import time; print(f"{time.monotonic():.6f}")' 2>/dev/null || echo "0.0")"
        local wall_sec="NA"
        if [[ -n "${PIPELINE_START_MONOTONIC:-}" ]]; then
            wall_sec="$(python3 -c "print(f'{(float(\"${end_mono}\") - float(\"${PIPELINE_START_MONOTONIC}\")):.3f}')" 2>/dev/null || echo "NA")"
        fi

        local manifest_status="failed"
        case "${status}" in
            interrupted) manifest_status="interrupted" ;;
            terminated) manifest_status="terminated" ;;
            *) manifest_status="failed" ;;
        esac

        if [[ "${status}" == "interrupted" || "${status}" == "terminated" ]]; then
            python3 "${b_writer}" stage-signal \
                --stage-status-path "${STAGE_STATUS_FILE}" \
                --run-id "${RUN_ID}" \
                --stage-id "${current_stg}" \
                --signal "$([[ "${status}" == interrupted ]] && echo INT || echo TERM)" \
                --exit-code "${exit_code}" \
                --end-time-utc "${end_ts}" \
                --error-id "${error_id}" \
                --message "Pipeline interrupted by signal: ${status}" >/dev/null 2>&1 || true
        else
            python3 "${b_writer}" stage-fail \
                --stage-status-path "${STAGE_STATUS_FILE}" \
                --run-id "${RUN_ID}" \
                --stage-id "${current_stg}" \
                --planned-action "${CURRENT_STAGE_ACTION:-run}" \
                --exit-code "${exit_code}" \
                --end-time-utc "${end_ts}" \
                --primary-log "${CURRENT_STAGE_LOG:-NA}" \
                --error-id "${error_id}" \
                --message "Pipeline stage failure in ${current_stg}" >/dev/null 2>&1 || true
        fi

        if [[ -f "${RUN_MANIFEST_FILE:-}" ]]; then
            local err_count=1
            if [[ -f "${LOGS_DIR:-}/errors.tsv" ]]; then
                err_count="$(awk 'NR>1 {count++} END {print count+0}' "${LOGS_DIR}/errors.tsv" 2>/dev/null || echo 1)"
            fi

            python3 "${b_writer}" finalize-run \
                --run-manifest-path "${RUN_MANIFEST_FILE}" \
                --run-id "${RUN_ID}" \
                --final-status "${manifest_status}" \
                --exit-code "${exit_code}" \
                --end-time-utc "${end_ts}" \
                --wall-time-seconds "${wall_sec}" \
                --warning-count "0" \
                --error-count "${err_count}" >/dev/null 2>&1 || true
        fi
    fi

    return 0
}

record_fallback_error() {
    local rc="$1"
    if (( rc == 0 )); then
        return 0
    fi
    if [[ "${DRY_RUN:-false}" == true || "${PLAN_MODE:-false}" == true ]]; then
        return 0
    fi

    local stg="${CURRENT_STAGE:-pipeline}"
    local stg_name="${CURRENT_STAGE_NAME:-$(stage_display_name "${stg}")}"
    local stg_action="${CURRENT_STAGE_ACTION:-unknown}"
    local stg_log="${CURRENT_STAGE_LOG:-$(stage_primary_log "${stg}")}"
    local runtime="${CURRENT_STAGE_RUNTIME:-unknown}"
    local final_cmd="${CURRENT_STAGE_COMMAND:-unknown}"

    local trap_src_file="${BASH_SOURCE[0]:-unknown}"
    local trap_src_line="${LINENO:-0}"
    local call_src_file="${CURRENT_COMMAND_SOURCE_FILE:-${trap_src_file}}"
    local call_src_line="${CURRENT_COMMAND_SOURCE_LINE:-${trap_src_line}}"

    local status="failed"
    local is_primary=false err_id="" err_key="" state=""

    if [[ ! -f "${WRITE_ERROR_RECORD}" ]]; then
        printf '[FATAL] %s [%s] Pipeline fallback failure: exit=%s line=%s\n' \
            "${RUN_ID}" "$$" "${rc}" "${trap_src_line}" >&2
        return 0
    fi

    while IFS='=' read -r k v; do
        case "${k}" in
            IS_PRIMARY) is_primary="${v}" ;;
            ERROR_ID) err_id="${v}" ;;
            ERROR_KEY) err_key="${v}" ;;
            STATE) state="${v}" ;;
        esac
    done < <(python3 "${WRITE_ERROR_RECORD}" reserve \
        --outdir "${OUTDIR}" \
        --run-id "${RUN_ID}" \
        --stage-id "${stg}" \
        --stage-name "${stg_name}" \
        --stage-action "${stg_action}" \
        --exit-code "${rc}" \
        --trap-source-file "${trap_src_file}" \
        --trap-source-line "${trap_src_line}" \
        --call-source-file "${call_src_file}" \
        --call-source-line "${call_src_line}" \
        --command "${final_cmd}" \
        --pid "$$")

    if [[ "${is_primary}" != true ]]; then
        return 0
    fi

    local run_status_update_status="not-available"
    record_pipeline_termination_state "${status}" "${rc}" "${err_id}"
    local rec_rc=$?
    if (( rec_rc == 0 )); then
        run_status_update_status="updated"
    elif (( rec_rc == 1 )); then
        run_status_update_status="failed"
    else
        run_status_update_status="not-available"
    fi

    cleanup_safe_temporaries "failure"
    python3 "${WRITE_ERROR_RECORD}" commit \
        --outdir "${OUTDIR}" \
        --run-id "${RUN_ID}" \
        --error-id "${err_id}" \
        --error-key "${err_key}" \
        --stage-id "${stg}" \
        --stage-name "${stg_name}" \
        --stage-action "${stg_action}" \
        --exit-code "${rc}" \
        --signal "none" \
        --trap-source-file "${trap_src_file}" \
        --trap-source-line "${trap_src_line}" \
        --call-source-file "${call_src_file}" \
        --call-source-line "${call_src_line}" \
        --function-name "record_fallback_error" \
        --command "${final_cmd}" \
        --log-path "${stg_log}" \
        --runtime "${runtime}" \
        --last-valid-checkpoint "${LAST_VALID_CHECKPOINT:-none}" \
        --last-valid-checkpoint-path "${LAST_VALID_CHECKPOINT_PATH:-none}" \
        --last-valid-product "${LAST_VALID_PRODUCT:-none}" \
        --last-valid-product-fingerprint "${LAST_VALID_PRODUCT_FINGERPRINT:-none}" \
        --message "Uncaught script termination fallback" \
        --status "${status}" \
        --run-status-update-status "${run_status_update_status}" \
        --pipeline-version "${PIPELINE_VERSION}" \
        --hostname "${PIPELINE_HOSTNAME}" \
        --pid "$$" \
        --temporary-directory "${PIPELINE_RUN_TMPDIR:-none}" \
        --temporary-files-preserved "${TEMPORARY_FILES_PRESERVED:-false}" \
        --temporary-cleanup-status "${TEMPORARY_CLEANUP_STATUS:-not-applicable}" >/dev/null 2>&1 || true

    printf '[FATAL] %s | stage=%s | action=%s | exit=%s | line=%s | log=%s | last_checkpoint=%s\n' \
        "${err_id}" "${stg}" "${stg_action}" "${rc}" "${trap_src_line}" "${stg_log}" "${LAST_VALID_CHECKPOINT:-none}" >&2
}

pipeline_error_handler() {
    local rc="$1"
    local line="$2"
    local trap_cmd="$3"

    trap - ERR
    set +e

    if (( ERROR_HANDLER_ACTIVE == 1 || ERROR_ALREADY_RECORDED == 1 )); then
        exit "${rc}"
    fi
    ERROR_HANDLER_ACTIVE=1

    local trap_src_file="${BASH_SOURCE[0]:-unknown}"
    local trap_src_line="${line}"
    local call_src_file="${CURRENT_COMMAND_SOURCE_FILE:-${trap_src_file}}"
    local call_src_line="${CURRENT_COMMAND_SOURCE_LINE:-${trap_src_line}}"
    local func_name="${FUNCNAME[1]:-${FUNCNAME[0]:-unknown}}"

    local stg="${CURRENT_STAGE:-pipeline}"
    local stg_name="${CURRENT_STAGE_NAME:-$(stage_display_name "${stg}")}"
    local stg_action="${CURRENT_STAGE_ACTION:-unknown}"
    local stg_log="${CURRENT_STAGE_LOG:-$(stage_primary_log "${stg}")}"
    local runtime="${CURRENT_STAGE_RUNTIME:-unknown}"

    local final_cmd="${CURRENT_STAGE_COMMAND:-}"
    if [[ -z "${final_cmd}" ]]; then
        final_cmd="$(sanitize_command "${trap_cmd}")"
    fi
    [[ -z "${final_cmd}" ]] && final_cmd="unknown"

    local status="failed"
    if [[ "${stg_action}" == "preflight" || "${stg}" == "pipeline" ]]; then
        status="preflight-failed"
    fi

    if [[ "${DRY_RUN:-false}" == true || "${IS_PLAN_MODE:-false}" == true ]]; then
        if [[ "${IS_PLAN_MODE:-false}" == true ]]; then
            printf '[FATAL] [PLAN] stage=%s | action=%s | exit=%s | line=%s | command=%s\n' \
                "${stg}" "${stg_action}" "${rc}" "${trap_src_line}" "${final_cmd}" >&2
        else
            printf '[FATAL] [DRY-RUN] stage=%s | action=%s | exit=%s | line=%s | command=%s | would write errors.tsv in %s/logs/\n' \
                "${stg}" "${stg_action}" "${rc}" "${trap_src_line}" "${final_cmd}" "${OUTDIR}" >&2
        fi
        PIPELINE_TERMINATION_RECORDED=true
        ERROR_ALREADY_RECORDED=1
        exit "${rc}"
    fi

    if [[ ! -f "${WRITE_ERROR_RECORD}" ]]; then
        printf '[FATAL] Pipeline error handler: exit=%s line=%s (error record manager unavailable)\n' \
            "${rc}" "${trap_src_line}" >&2
        PIPELINE_TERMINATION_RECORDED=true
        ERROR_ALREADY_RECORDED=1
        exit "${rc}"
    fi

    local is_primary=false err_id="" err_key="" state=""
    while IFS='=' read -r k v; do
        case "${k}" in
            IS_PRIMARY) is_primary="${v}" ;;
            ERROR_ID) err_id="${v}" ;;
            ERROR_KEY) err_key="${v}" ;;
            STATE) state="${v}" ;;
        esac
    done < <(python3 "${WRITE_ERROR_RECORD}" reserve \
        --outdir "${OUTDIR}" \
        --run-id "${RUN_ID}" \
        --stage-id "${stg}" \
        --stage-name "${stg_name}" \
        --stage-action "${stg_action}" \
        --exit-code "${rc}" \
        --trap-source-file "${trap_src_file}" \
        --trap-source-line "${trap_src_line}" \
        --call-source-file "${call_src_file}" \
        --call-source-line "${call_src_line}" \
        --command "${final_cmd}" \
        --pid "$$")

    if [[ "${is_primary}" == true ]]; then
        local run_status_update_status="not-available"
        record_pipeline_termination_state "${status}" "${rc}" "${err_id}"
        local rec_rc=$?
        if (( rec_rc == 0 )); then
            run_status_update_status="updated"
        elif (( rec_rc == 1 )); then
            run_status_update_status="failed"
        else
            run_status_update_status="not-available"
        fi

        cleanup_safe_temporaries "failure"
        python3 "${WRITE_ERROR_RECORD}" commit \
            --outdir "${OUTDIR}" \
            --run-id "${RUN_ID}" \
            --error-id "${err_id}" \
            --error-key "${err_key}" \
            --stage-id "${stg}" \
            --stage-name "${stg_name}" \
            --stage-action "${stg_action}" \
            --exit-code "${rc}" \
            --signal "none" \
            --trap-source-file "${trap_src_file}" \
            --trap-source-line "${trap_src_line}" \
            --call-source-file "${call_src_file}" \
            --call-source-line "${call_src_line}" \
            --function-name "${func_name}" \
            --command "${final_cmd}" \
            --log-path "${stg_log}" \
            --runtime "${runtime}" \
            --last-valid-checkpoint "${LAST_VALID_CHECKPOINT:-none}" \
            --last-valid-checkpoint-path "${LAST_VALID_CHECKPOINT_PATH:-none}" \
            --last-valid-product "${LAST_VALID_PRODUCT:-none}" \
            --last-valid-product-fingerprint "${LAST_VALID_PRODUCT_FINGERPRINT:-none}" \
            --message "Pipeline stage failure in ${stg}" \
            --status "${status}" \
            --run-status-update-status "${run_status_update_status}" \
            --pipeline-version "${PIPELINE_VERSION}" \
            --hostname "${PIPELINE_HOSTNAME}" \
            --pid "$$" \
            --temporary-directory "${PIPELINE_RUN_TMPDIR:-none}" \
            --temporary-files-preserved "${TEMPORARY_FILES_PRESERVED:-false}" \
            --temporary-cleanup-status "${TEMPORARY_CLEANUP_STATUS:-not-applicable}" >/dev/null 2>&1 || true

        printf '[FATAL] %s | stage=%s | action=%s | exit=%s | line=%s | log=%s | last_checkpoint=%s\n' \
            "${err_id}" "${stg}" "${stg_action}" "${rc}" "${trap_src_line}" "${stg_log}" "${LAST_VALID_CHECKPOINT:-none}" >&2
    fi

    PIPELINE_TERMINATION_RECORDED=true
    ERROR_ALREADY_RECORDED=1
    exit "${rc}"
}

pipeline_signal_handler() {
    local sig="$1"
    local rc="$2"

    trap - ERR INT TERM HUP
    set +e

    if (( ERROR_HANDLER_ACTIVE == 1 || ERROR_ALREADY_RECORDED == 1 )); then
        exit "${rc}"
    fi
    ERROR_HANDLER_ACTIVE=1

    local status="terminated"
    [[ "${sig}" == "INT" ]] && status="interrupted"

    local stg="${CURRENT_STAGE:-pipeline}"
    local stg_name="${CURRENT_STAGE_NAME:-$(stage_display_name "${stg}")}"
    local stg_action="${CURRENT_STAGE_ACTION:-unknown}"
    local stg_log="${CURRENT_STAGE_LOG:-$(stage_primary_log "${stg}")}"
    local runtime="${CURRENT_STAGE_RUNTIME:-unknown}"
    local final_cmd="${CURRENT_STAGE_COMMAND:-unknown}"

    local trap_src_file="${BASH_SOURCE[0]:-unknown}"
    local trap_src_line="${LINENO:-0}"
    local call_src_file="${CURRENT_COMMAND_SOURCE_FILE:-${trap_src_file}}"
    local call_src_line="${CURRENT_COMMAND_SOURCE_LINE:-${trap_src_line}}"

    if [[ "${DRY_RUN:-false}" == true ]]; then
        printf '[FATAL] [DRY-RUN] stage=%s | action=%s | exit=%s | signal=%s | line=%s\n' \
            "${stg}" "${stg_action}" "${rc}" "${sig}" "${trap_src_line}" >&2
        PIPELINE_TERMINATION_RECORDED=true
        ERROR_ALREADY_RECORDED=1
        exit "${rc}"
    fi

    if [[ ! -f "${WRITE_ERROR_RECORD}" ]]; then
        printf '[FATAL] Pipeline signal handler: %s (exit=%s)\n' "${sig}" "${rc}" >&2
        PIPELINE_TERMINATION_RECORDED=true
        ERROR_ALREADY_RECORDED=1
        exit "${rc}"
    fi

    local is_primary=false err_id="" err_key="" state=""
    while IFS='=' read -r k v; do
        case "${k}" in
            IS_PRIMARY) is_primary="${v}" ;;
            ERROR_ID) err_id="${v}" ;;
            ERROR_KEY) err_key="${v}" ;;
            STATE) state="${v}" ;;
        esac
    done < <(python3 "${WRITE_ERROR_RECORD}" reserve \
        --outdir "${OUTDIR}" \
        --run-id "${RUN_ID}" \
        --stage-id "${stg}" \
        --stage-name "${stg_name}" \
        --stage-action "${stg_action}" \
        --exit-code "${rc}" \
        --trap-source-file "${trap_src_file}" \
        --trap-source-line "${trap_src_line}" \
        --call-source-file "${call_src_file}" \
        --call-source-line "${call_src_line}" \
        --command "${final_cmd}" \
        --pid "$$")

    if [[ "${is_primary}" == true ]]; then
        local run_status_update_status="not-available"
        record_pipeline_termination_state "${status}" "${rc}" "${err_id}"
        local rec_rc=$?
        if (( rec_rc == 0 )); then
            run_status_update_status="updated"
        elif (( rec_rc == 1 )); then
            run_status_update_status="failed"
        else
            run_status_update_status="not-available"
        fi

        cleanup_safe_temporaries "signal"
        python3 "${WRITE_ERROR_RECORD}" commit \
            --outdir "${OUTDIR}" \
            --run-id "${RUN_ID}" \
            --error-id "${err_id}" \
            --error-key "${err_key}" \
            --stage-id "${stg}" \
            --stage-name "${stg_name}" \
            --stage-action "${stg_action}" \
            --exit-code "${rc}" \
            --signal "${sig}" \
            --trap-source-file "${trap_src_file}" \
            --trap-source-line "${trap_src_line}" \
            --call-source-file "${call_src_file}" \
            --call-source-line "${call_src_line}" \
            --function-name "pipeline_signal_handler" \
            --command "${final_cmd}" \
            --log-path "${stg_log}" \
            --runtime "${runtime}" \
            --last-valid-checkpoint "${LAST_VALID_CHECKPOINT:-none}" \
            --last-valid-checkpoint-path "${LAST_VALID_CHECKPOINT_PATH:-none}" \
            --last-valid-product "${LAST_VALID_PRODUCT:-none}" \
            --last-valid-product-fingerprint "${LAST_VALID_PRODUCT_FINGERPRINT:-none}" \
            --message "Pipeline interrupted by signal ${sig}" \
            --status "${status}" \
            --run-status-update-status "${run_status_update_status}" \
            --pipeline-version "${PIPELINE_VERSION}" \
            --hostname "${PIPELINE_HOSTNAME}" \
            --pid "$$" \
            --temporary-directory "${PIPELINE_RUN_TMPDIR:-none}" \
            --temporary-files-preserved "${TEMPORARY_FILES_PRESERVED:-false}" \
            --temporary-cleanup-status "${TEMPORARY_CLEANUP_STATUS:-not-applicable}" >/dev/null 2>&1 || true

        printf '[FATAL] %s | stage=%s | action=%s | exit=%s | signal=%s | line=%s | log=%s | last_checkpoint=%s\n' \
            "${err_id}" "${stg}" "${stg_action}" "${rc}" "${sig}" "${trap_src_line}" "${stg_log}" "${LAST_VALID_CHECKPOINT:-none}" >&2
    fi

    PIPELINE_TERMINATION_RECORDED=true
    ERROR_ALREADY_RECORDED=1
    exit "${rc}"
}

pipeline_exit_handler() {
    local rc="$1"
    if (( rc == 0 )); then
        return 0
    fi
    if [[ "${PLAN_MODE:-false}" == true ]]; then
        return 0
    fi
    if (( ERROR_HANDLER_ACTIVE == 1 || ERROR_ALREADY_RECORDED == 1 )); then
        return 0
    fi
    if [[ "${PIPELINE_TERMINATION_RECORDED}" == true ]]; then
        return 0
    fi

    record_fallback_error "${rc}"
}

install_pipeline_traps() {
    trap 'pipeline_error_handler "$?" "$LINENO" "$BASH_COMMAND"' ERR
    trap 'pipeline_signal_handler INT 130' INT
    trap 'pipeline_signal_handler TERM 143' TERM
    trap 'pipeline_signal_handler HUP 129' HUP
    trap 'pipeline_exit_handler "$?"' EXIT
}

mark_done() {
    local step="$1" tool key
    local -a stage_tools=()

    if [[ "${DRY_RUN}" == true ]]; then
        describe_step_validations "${step}"
        log_info "[Step ${step}] Dry-run: checkpoint not written."
        STAGE_FINAL_STATUS["${step}"]=dry-run
        return 0
    fi

    validate_step_products "${step}" || return 1

    mapfile -t stage_tools < <(inventory_stage_tools "${step}")
    for tool in "${stage_tools[@]}"; do
        key="${step}|${tool}"
        if [[ "${ID_VERSION[$key]:-unknown}" == unknown && \
              "${ID_RUNTIME_FINGERPRINT[$key]:-unknown}" == unknown ]]; then
            log_error "[Checkpoint] Refusing to write checkpoint for ${step}: ${tool} has neither a reliable version nor a functional fingerprint."
            return 1
        fi
    done

    checkpoint_config_manifest "${step}" > "${CHECKPOINT_DIR}/${step}.config"
    touch "${CHECKPOINT_DIR}/${step}.done"
    STAGE_CHECKPOINT_WRITTEN["${step}"]=true

    local main_prod="none" prod_fp="none"
    case "${step}" in
        01_fastp) main_prod="${FASTP_CLEAN_SAMPLE_MANIFEST:-none}" ;;
        02_megahit) main_prod="${CONTIGS:-none}" ;;
        03_kraken2) main_prod="${KRAKEN2_DIR}/reports" ;;
        04_mapping) main_prod="${DEPTH_FILE:-none}" ;;
        05a_vamb) main_prod="${VAMB_DIR}/final_bins" ;;
        05b_metabinner) main_prod="${METABINNER_DIR}/final_bins" ;;
        05c_comebin) main_prod="${COMEBIN_DIR}/final_bins" ;;
        06_refinement) main_prod="${REFINED_BINS_DIR:-none}" ;;
        07_reassembly) main_prod="${REASSEMBLY_DIR}" ;;
        08_checkm2) main_prod="${CHECKM2_DIR}/quality_report.tsv" ;;
        09_filter_mags) main_prod="${MAGS_DIR}" ;;
        10_gtdbtk) main_prod="${GTDBTK_DIR}" ;;
        11_coverm) main_prod="${ABUNDANCE_DIR}" ;;
        12a_taxonomy_summary) main_prod="${SUMMARY_DIR}/taxonomy_summary.tsv" ;;
        12b_mag_summary) main_prod="${SUMMARY_DIR}/mag_summary.tsv" ;;
        12c_run_overview) main_prod="${SUMMARY_DIR}/run_overview.tsv" ;;
        *) main_prod="none" ;;
    esac
    [[ -n "${main_prod}" && -f "${main_prod}" ]] && prod_fp="$(sha256sum "${main_prod}" 2>/dev/null | awk '{print $1}')"

    set_last_valid_checkpoint "${step}" "${CHECKPOINT_DIR}/${step}.done" "${main_prod}" "${prod_fp:-none}"
    set_last_valid_product "${main_prod}" "${prod_fp:-none}"
    log_info "[Step ${step}] Checkpoint written."
}

# ── Directory Layout ──────────────────────────────────────────────────────────
QC_DIR="${OUTDIR}/01_qc"
ASSEMBLY_DIR="${OUTDIR}/02_assembly"
KRAKEN2_DIR="${OUTDIR}/03_kraken2"
MAPPING_DIR="${OUTDIR}/04_mapping"
BINNING_DIR="${OUTDIR}/05_binning"
VAMB_DIR="${BINNING_DIR}/vamb"
VAMB_ABUNDANCE_TSV="${VAMB_DIR}/vamb_abundance.tsv"
METABINNER_DIR="${BINNING_DIR}/metabinner"
COMEBIN_DIR="${BINNING_DIR}/comebin"
REFINEMENT_DIR="${OUTDIR}/06_refinement"
VAMB_FINAL_BINS_DIR="${VAMB_DIR}/final_bins"
METABINNER_FINAL_BINS_DIR="${METABINNER_DIR}/final_bins"
COMEBIN_FINAL_BINS_DIR="${COMEBIN_DIR}/final_bins"
REFINED_BINS_DIR="${REFINEMENT_DIR}/metawrap_${METAWRAP_COMPLETENESS}_${METAWRAP_CONTAMINATION}_bins"
EXTERNAL_BINS_CACHE_ROOT="${OUTDIR}/00_inputs/external_bins"
EXTERNAL_BINS_CURRENT="${EXTERNAL_BINS_CACHE_ROOT}/current"
EXTERNAL_FINAL_BINS_DIR="${EXTERNAL_BINS_CURRENT}/final_bins"
EXTERNAL_BINS_MANIFEST="${EXTERNAL_BINS_CURRENT}/external_bins_manifest.tsv"
REASSEMBLY_DIR="${OUTDIR}/07_reassembly"
CHECKM2_DIR="${OUTDIR}/08_checkm2"
MAGS_DIR="${OUTDIR}/09_mags"
MAG_FILTER_MANIFEST="${MAGS_DIR}/mag_filter_manifest.tsv"
GTDBTK_DIR="${OUTDIR}/10_gtdbtk"
ABUNDANCE_DIR="${OUTDIR}/11_abundance"
SUMMARY_DIR="${OUTDIR}/12_summary"
LOGS_DIR="${OUTDIR}/logs"
CHECKPOINT_DIR="${OUTDIR}/.checkpoints"
TMPSCRIPTS_DIR="${OUTDIR}/logs/mamba_scripts"
DISK_SPACE_CURRENT_TSV="${LOGS_DIR}/disk_space.current.tsv"
DISK_SPACE_TSV="${LOGS_DIR}/disk_space.tsv"
DISK_PREFLIGHT_STATUS="pending"
DISK_WARNING_COUNT=0
OUTDIR_FS_ID=""
TMPDIR_FS_ID=""
SHARED_CAPACITY_FS="unknown"

# Install full traps now that directories and error system are initialized
install_pipeline_traps

# ── Global State Variables ────────────────────────────────────────────────────
CONTIGS=""
CLEAN_R1_FILES=()
CLEAN_R2_FILES=()
SORTED_BAMS=()
ACTIVE_BAMS=()
ACTIVE_BAM_INDEXES=()
ACTIVE_BAM_SAMPLE_IDS=()
INPUT_SAMPLE_IDS=()
INPUT_R1_FILES=()
INPUT_R2_FILES=()
DEPTH_FILE=""
BINS_FOR_QC=""
BINS_FOR_REASSEMBLY=""
SAMPLE_MANIFESTS_DIR="${OUTDIR}/00_inputs/manifests"
RAW_SAMPLE_MANIFEST="${SAMPLE_MANIFESTS_DIR}/raw_samples.tsv"
RAW_SAMPLE_METADATA="${SAMPLE_MANIFESTS_DIR}/raw_samples.meta.json"
EXTERNAL_CLEAN_SAMPLE_MANIFEST="${SAMPLE_MANIFESTS_DIR}/external_clean_samples.tsv"
EXTERNAL_CLEAN_SAMPLE_METADATA="${SAMPLE_MANIFESTS_DIR}/external_clean_samples.meta.json"
FASTP_CLEAN_SAMPLE_MANIFEST="${QC_DIR}/fastp_clean_samples.tsv"
FASTP_CLEAN_SAMPLE_METADATA="${QC_DIR}/fastp_clean_samples.meta.json"
HQ_BINS_DIR=""
BINNING_RESOLUTION_ERRORS=0
# Resource detection state (populated by detect_resources)
_PHYS_CORES=0
_PHYS_MEM_GB=0

STAGE_IDS=(01_fastp 03_kraken2 02_megahit 04_mapping 05a_vamb 05b_metabinner 05c_comebin
           06_refinement 07_reassembly 08_checkm2 09_filter_mags 10_gtdbtk 11_coverm
           12a_taxonomy_summary 12b_mag_summary 12c_run_overview)

resolve_stage_id() {
    case "$1" in
        fastp) printf '%s\n' 01_fastp ;;
        kraken2) printf '%s\n' 03_kraken2 ;;
        megahit) printf '%s\n' 02_megahit ;;
        mapping) printf '%s\n' 04_mapping ;;
        vamb) printf '%s\n' 05a_vamb ;;
        metabinner) printf '%s\n' 05b_metabinner ;;
        comebin) printf '%s\n' 05c_comebin ;;
        refinement) printf '%s\n' 06_refinement ;;
        reassembly) printf '%s\n' 07_reassembly ;;
        checkm2) printf '%s\n' 08_checkm2 ;;
        mag-filtering) printf '%s\n' 09_filter_mags ;;
        gtdbtk) printf '%s\n' 10_gtdbtk ;;
        coverm) printf '%s\n' 11_coverm ;;
        taxonomy-summary) printf '%s\n' 12a_taxonomy_summary ;;
        mag-summary) printf '%s\n' 12b_mag_summary ;;
        run-overview) printf '%s\n' 12c_run_overview ;;
        summary)
            log_error "The monolithic 'summary' stage has been replaced by 'taxonomy-summary', 'mag-summary', and 'run-overview'. Please specify one of these canonical stage names."
            return 1 ;;
        *)
            log_error "Unknown stage '$1'. Allowed stages: fastp, kraken2, megahit, mapping, vamb, metabinner, comebin, refinement, reassembly, checkm2, mag-filtering, gtdbtk, coverm, taxonomy-summary, mag-summary, run-overview."
            return 1 ;;
    esac
}

get_stage_public_name() {
    case "$1" in
        01_fastp) printf '%s\n' fastp ;;
        03_kraken2) printf '%s\n' kraken2 ;;
        02_megahit) printf '%s\n' megahit ;;
        04_mapping) printf '%s\n' mapping ;;
        05a_vamb) printf '%s\n' vamb ;;
        05b_metabinner) printf '%s\n' metabinner ;;
        05c_comebin) printf '%s\n' comebin ;;
        06_refinement) printf '%s\n' refinement ;;
        07_reassembly) printf '%s\n' reassembly ;;
        08_checkm2) printf '%s\n' checkm2 ;;
        09_filter_mags) printf '%s\n' mag-filtering ;;
        10_gtdbtk) printf '%s\n' gtdbtk ;;
        11_coverm) printf '%s\n' coverm ;;
        12a_taxonomy_summary) printf '%s\n' taxonomy-summary ;;
        12b_mag_summary) printf '%s\n' mag-summary ;;
        12c_run_overview) printf '%s\n' run-overview ;;
        12_summary) printf '%s\n' summary ;;
        *) return 1 ;;
    esac
}

get_stage_display_name() {
    case "$1" in
        01_fastp) printf '%s\n' "Read quality control with fastp" ;;
        03_kraken2) printf '%s\n' "Read taxonomic profiling with Kraken2/Bracken" ;;
        02_megahit) printf '%s\n' "Metagenome assembly with MEGAHIT" ;;
        04_mapping) printf '%s\n' "Read mapping and coverage estimation" ;;
        05a_vamb) printf '%s\n' "Genome binning with Vamb" ;;
        05b_metabinner) printf '%s\n' "Genome binning with MetaBinner" ;;
        05c_comebin) printf '%s\n' "Genome binning with COMEBin" ;;
        06_refinement) printf '%s\n' "Consensus bin refinement with MetaWRAP" ;;
        07_reassembly) printf '%s\n' "MAG reassembly with MetaWRAP" ;;
        08_checkm2) printf '%s\n' "MAG quality assessment with CheckM2" ;;
        09_filter_mags) printf '%s\n' "MAG quality filtering" ;;
        10_gtdbtk) printf '%s\n' "Taxonomic classification with GTDB-Tk" ;;
        11_coverm) printf '%s\n' "MAG abundance estimation with CoverM" ;;
        12a_taxonomy_summary) printf '%s\n' "Read-based taxonomic summary" ;;
        12b_mag_summary) printf '%s\n' "MAG quality, taxonomy, and abundance summary" ;;
        12c_run_overview) printf '%s\n' "Pipeline run overview and execution metrics" ;;
        12_summary) printf '%s\n' "Legacy summary" ;;
        *) return 1 ;;
    esac
}

stage_is_active() {
    [[ "${STAGE_STATE[$1]:-}" == active ]]
}

binner_requires_input() {
    local binner="$1" input="$2"
    case "${binner}:${input}" in
        vamb:assembly|vamb:jgi-depth) return 0 ;;
        metabinner:assembly|metabinner:jgi-depth) return 0 ;;
        comebin:assembly|comebin:bam|comebin:bam-index) return 0 ;;
        *) return 1 ;;
    esac
}

validate_binner_inputs() {
    local binner="$1" input
    for input in assembly bam bam-index jgi-depth vamb-abundance metabinner-abundance comebin-abundance; do
        binner_requires_input "${binner}" "${input}" || continue
        case "${input}" in
            assembly)
                [[ -n "${CONTIGS}" ]] || {
                    log_error "[Preflight] ${binner} requires an assembly; provide --assembly-in or enable MEGAHIT."
                    return 1
                } ;;
            bam|bam-index)
                if [[ "${COVERAGE_SOURCE}" == none || \
                      ( "${COVERAGE_SOURCE}" == external && ${#ACTIVE_BAMS[@]} -eq 0 ) ]]; then
                    log_error "[Preflight] ${binner} requires ${input} coverage; provide --bam-dir/--bam-manifest or enable internal mapping."
                    return 1
                fi ;;
            jgi-depth)
                if [[ "${COVERAGE_SOURCE}" == none || -z "${JGI_DEPTH_FILE_ACTIVE}" ]]; then
                    log_error "[Preflight] ${binner} requires JGI depth; provide --jgi-depth-file or enable internal mapping."
                    return 1
                fi ;;
        esac
    done
}

validate_vamb_fasta() {
    local fasta="$1"
    [[ -f "${fasta}" && -s "${fasta}" ]] || {
        log_error "[Vamb] Assembly FASTA does not exist or is empty: ${fasta}"
        return 1
    }
    python3 - "${fasta}" <<'PY_VAMB_FASTA'
import sys

path = sys.argv[1]
identifiers = set()
headers = 0
has_sequence = False
with open(path, encoding="utf-8") as fasta:
    for line_number, line in enumerate(fasta, 1):
        if line.startswith(">"):
            headers += 1
            fields = line[1:].split()
            if not fields:
                raise SystemExit(f"[Vamb] Empty FASTA identifier at line {line_number}: {path}")
            identifier = fields[0]
            if identifier in identifiers:
                raise SystemExit(f"[Vamb] Duplicate FASTA identifier '{identifier}': {path}")
            identifiers.add(identifier)
        elif headers > 0 and line.strip():
            has_sequence = True
if headers == 0:
    raise SystemExit(f"[Vamb] FASTA contains no headers: {path}")
if not has_sequence:
    raise SystemExit(f"[Vamb] FASTA contains no sequence: {path}")
print(f"[Vamb] FASTA validated: {headers} headers, at least one sequence.")
PY_VAMB_FASTA
}

validate_vamb_binsplit_headers() {
    local fasta="$1" separator="$2"
    python3 - "${fasta}" "${separator}" <<'PY_VAMB_BINSPLIT'
import sys

path, separator = sys.argv[1:]
if len(separator) != 1 or separator.isspace() or separator in "/>":
    raise SystemExit("[Vamb] Binsplit separator must be exactly one non-whitespace character other than '/' or '>'.")
checked = 0
first_invalid = None
with open(path, encoding="utf-8") as fasta:
    for line in fasta:
        if not line.startswith(">"):
            continue
        fields = line[1:].split()
        identifier = fields[0] if fields else ""
        checked += 1
        position = identifier.find(separator)
        if (position <= 0 or position == len(identifier) - 1) and first_invalid is None:
            first_invalid = identifier
if checked == 0:
    raise SystemExit(f"[Vamb] FASTA contains no headers; cannot validate binsplitting: {path}")
if first_invalid is not None:
    raise SystemExit(
        f"[Vamb] Invalid binsplit identifier '{first_invalid}': separator {separator!r} "
        f"must occur between non-empty sample and contig IDs; checked {checked} headers."
    )
print(f"[Vamb] Binsplit separator {separator!r} validated across {checked} FASTA headers.")
PY_VAMB_BINSPLIT
}

validate_vamb_abundance_tsv() {
    local tsv="$1" fasta="$2"
    [[ -f "${tsv}" && -s "${tsv}" ]] || {
        log_error "[Vamb] Abundance TSV does not exist or is empty: ${tsv}"
        return 1
    }
    python3 - "${tsv}" "${fasta}" <<'PY_VAMB_ABUNDANCE'
import math
import re
import sys

tsv_path, fasta_path = sys.argv[1:]
fasta_ids = []
seen_fasta = set()
with open(fasta_path, encoding="utf-8") as fasta:
    for line_number, line in enumerate(fasta, 1):
        if line.startswith(">"):
            fields = line[1:].split()
            if not fields:
                raise SystemExit(f"[Vamb] Empty FASTA identifier at line {line_number}: {fasta_path}")
            identifier = fields[0]
            if identifier in seen_fasta:
                raise SystemExit(f"[Vamb] Duplicate FASTA identifier '{identifier}': {fasta_path}")
            seen_fasta.add(identifier)
            fasta_ids.append(identifier)
if not fasta_ids:
    raise SystemExit(f"[Vamb] FASTA contains no headers: {fasta_path}")

with open(tsv_path, encoding="utf-8") as stream:
    header_line = stream.readline()
    if not header_line:
        raise SystemExit(f"[Vamb] Abundance TSV has no header: {tsv_path}")
    header = header_line.rstrip("\r\n").split("\t")
    if not header or header[0] != "contigname":
        raise SystemExit('[Vamb] Abundance TSV first column must be exactly "contigname".')
    if len(header) < 2:
        raise SystemExit("[Vamb] Abundance TSV must contain at least one sample column.")
    samples = header[1:]
    if any(sample == "" for sample in samples):
        raise SystemExit("[Vamb] Abundance TSV contains an empty sample name.")
    seen_samples = set()
    for sample in samples:
        if sample in seen_samples:
            raise SystemExit(f"[Vamb] Duplicate sample name in abundance TSV header: '{sample}'.")
        seen_samples.add(sample)

    abundance_ids = []
    seen_abundance = set()
    incompatible = []
    numeric_pattern = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")
    for line_number, line in enumerate(stream, 2):
        row = line.rstrip("\r\n").split("\t")
        if len(row) != len(header):
            raise SystemExit(
                f"[Vamb] Abundance TSV line {line_number} has {len(row)} columns; expected {len(header)}."
            )
        identifier = row[0]
        if not identifier:
            raise SystemExit(f"[Vamb] Empty contigname in abundance TSV line {line_number}.")
        if identifier in seen_abundance:
            raise SystemExit(f"[Vamb] Duplicate contig ID in abundance TSV: '{identifier}'.")
        seen_abundance.add(identifier)
        abundance_ids.append(identifier)
        if identifier not in seen_fasta:
            incompatible.append(identifier)
        for column, value in enumerate(row[1:], 2):
            if not numeric_pattern.fullmatch(value):
                raise SystemExit(
                    f"[Vamb] Non-numeric abundance at line {line_number}, column {column}: {value!r}."
                )
            number = float(value)
            if not math.isfinite(number):
                raise SystemExit(
                    f"[Vamb] Non-finite abundance at line {line_number}, column {column}: {value!r}."
                )
            if number < 0:
                raise SystemExit(
                    f"[Vamb] Negative abundance at line {line_number}, column {column}: {value!r}."
                )

if not abundance_ids:
    raise SystemExit(f"[Vamb] Abundance TSV has no data rows: {tsv_path}")
if incompatible:
    raise SystemExit(
        f"[Vamb] First abundance contig ID absent from FASTA: '{incompatible[0]}'; "
        f"total incompatible IDs: {len(incompatible)}."
    )

missing = [identifier for identifier in fasta_ids if identifier not in seen_abundance]
if missing:
    raise SystemExit(
        f"[Vamb] Vamb 5.0.4 requires abundance rows for the complete FASTA in the same order. "
        f"First FASTA ID absent from TSV: '{missing[0]}'; total missing IDs: {len(missing)}."
    )
if abundance_ids != fasta_ids:
    first = next(
        (index for index, (abundance_id, fasta_id) in enumerate(zip(abundance_ids, fasta_ids), 1)
         if abundance_id != fasta_id),
        None,
    )
    raise SystemExit(
        f"[Vamb] Abundance TSV contig order must match FASTA order; first mismatch at row {first}: "
        f"TSV='{abundance_ids[first - 1]}', FASTA='{fasta_ids[first - 1]}'."
    )
print(
    f"[Vamb] Abundance TSV validated: {len(abundance_ids)} contigs, "
    f"{len(samples)} samples; IDs match the complete FASTA in order."
)
PY_VAMB_ABUNDANCE
}

write_vamb_abundance_tsv() {
    local depth="$1" output="$2"
    awk 'BEGIN{OFS="\t"}
    NR==1{
        printf "contigname"
        for(i=4; i<=NF; i+=2) {
            n=$i; gsub(/.*\//,"",n); gsub(/_sorted\.bam$/,"",n)
            printf "\t%s", n
        }
        printf "\n"; next
    }
    {
        printf "%s", $1
        for(i=4; i<=NF; i+=2) printf "\t%s", $i
        printf "\n"
    }' "${depth}" > "${output}"
}

resolve_vamb_version() {
    local output help_output runtime="apptainer"
    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        if [[ ! -f "${SIF_VAMB}" ]]; then
            if [[ "${DRY_RUN}" == true ]]; then
                log_error "[Vamb] Dry-run cannot verify Vamb version because the image is unavailable: ${SIF_VAMB}"
                return 1
            fi
            ensure_sif_image "${SIF_VAMB}" "${SIF_VAMB_URI}" || return 1
        fi
        command -v apptainer &>/dev/null || runtime="singularity"
        output="$("${runtime}" exec --cleanenv "${SIF_VAMB}" vamb --version 2>&1)" || {
            log_error "[Vamb] Could not query Vamb version from ${SIF_VAMB}: ${output}"
            return 1
        }
        help_output="$("${runtime}" exec --cleanenv "${SIF_VAMB}" vamb bin default --help 2>&1)" || {
            log_error "[Vamb] Could not inspect 'vamb bin default --help' in ${SIF_VAMB}: ${help_output}"
            return 1
        }
    else
        output="$(mamba run -n "${ENV_VAMB}" vamb --version 2>&1)" || {
            log_error "[Vamb] Could not query Vamb version in Mamba environment ${ENV_VAMB}: ${output}"
            return 1
        }
        help_output="$(mamba run -n "${ENV_VAMB}" vamb bin default --help 2>&1)" || {
            log_error "[Vamb] Could not inspect 'vamb bin default --help' in ${ENV_VAMB}: ${help_output}"
            return 1
        }
    fi
    VAMB_VERSION="$(tr '\n' ' ' <<< "${output}" | sed 's/[[:space:]]*$//')"
    [[ -n "${VAMB_VERSION}" ]] || {
        log_error "[Vamb] Version query returned empty output."
        return 1
    }
    VAMB_INTERNAL_MIN_CONTIG_LEN="$(python3 -c '
import re, sys
for line in sys.stdin:
    match = re.match(r"^\s*-m\s+.*\[(\d+)\]", line)
    if match:
        print(match.group(1))
        break
' <<< "${help_output}")"
            VAMB_DEFAULT_THREADS="$(python3 -c '
        import re, sys
        for line in sys.stdin:
            match = re.match(r"^\s*-p\s+.*\[(\d+)\]", line)
            if match:
                print(match.group(1))
                break
        ' <<< "${help_output}")"
            if [[ -z "${VAMB_INTERNAL_MIN_CONTIG_LEN}" ]] || \
               [[ -z "${VAMB_DEFAULT_THREADS}" ]] || \
               ! grep -Eq -- '--minfasta[[:space:]]+.*Minimum bin size to output as fasta' <<< "${help_output}"; then
                log_error "[Vamb] Vamb help did not confirm the supported -m/-p defaults and --minfasta meaning."
                return 1
            fi
            log_info "[Vamb] Effective version: ${VAMB_VERSION}."
            log_info "[Vamb] Help confirms --minfasta is minimum bin FASTA size; -m ignores contigs shorter than ${VAMB_INTERNAL_MIN_CONTIG_LEN} bp; -p defaults to ${VAMB_DEFAULT_THREADS} threads."
}

build_vamb_args() {
    local output_name="$1"
    local -n vamb_args_ref="${output_name}"
    vamb_args_ref=(
        bin default
        --outdir "${VAMB_DIR}/vamb_out"
        --fasta "${CONTIGS}"
        --abundance_tsv "${VAMB_ABUNDANCE_TSV}"
        --minfasta "${VAMB_MINFASTA}"
    )
    if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
        vamb_args_ref+=(-o "${VAMB_BINSPLIT_SEPARATOR}")
    fi
}

write_vamb_script() {
    local script="$1" argument
    local -a vamb_args=()
    build_vamb_args vamb_args
    {
        printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'vamb_args=('
        for argument in "${vamb_args[@]}"; do
            printf '    %q\n' "${argument}"
        done
        printf '%s\n' ')' "vamb \"\${vamb_args[@]}\""
    } > "${script}"
    chmod +x "${script}"
}

validate_vamb_preflight() {
    local abundance_tmp="" argument
    local -a vamb_args=()
    [[ "${RUN_BINNING}" == true ]] && binner_selected vamb || return 0

    if [[ "${VAMB_BINSPLIT_SEPARATOR_WAS_SET}" == true && -z "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
        log_error "[Vamb] An explicitly supplied binsplit separator cannot be empty."
        return 1
    fi
    if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
        if ! python3 - "${VAMB_BINSPLIT_SEPARATOR}" <<'PY_VAMB_SEPARATOR'
import sys
separator = sys.argv[1]
if len(separator) != 1 or separator.isspace() or separator in "/>":
    raise SystemExit("[Vamb] Separator must be exactly one non-whitespace character other than '/' or '>'.")
PY_VAMB_SEPARATOR
        then
            return 1
        fi
    fi
    if [[ ! "${VAMB_MINFASTA}" =~ ^[0-9]+$ ]]; then
        log_error "[Vamb] VAMB_MINFASTA must be a non-negative integer (minimum output bin size in bp)."
        return 1
    fi

    if [[ -s "${CONTIGS}" ]]; then
        validate_vamb_fasta "${CONTIGS}" || return 1
        if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
            validate_vamb_binsplit_headers "${CONTIGS}" "${VAMB_BINSPLIT_SEPARATOR}" || return 1
        fi
        if [[ -n "${DEPTH_FILE}" && -s "${DEPTH_FILE}" ]]; then
            abundance_tmp="$(mktemp "${TMPDIR:-/tmp}/vamb-abundance.XXXXXX")" || return 1
            if ! write_vamb_abundance_tsv "${DEPTH_FILE}" "${abundance_tmp}" || \
               ! validate_vamb_abundance_tsv "${abundance_tmp}" "${CONTIGS}"; then
                rm -f -- "${abundance_tmp}"
                return 1
            fi
            rm -f -- "${abundance_tmp}"
        elif [[ -s "${VAMB_ABUNDANCE_TSV}" ]]; then
            validate_vamb_abundance_tsv "${VAMB_ABUNDANCE_TSV}" "${CONTIGS}" || return 1
        else
            log_info "[Vamb] Abundance validation is deferred until the active mapping stage creates its depth matrix."
        fi
    else
        if [[ "${RUN_MEGAHIT}" == true ]]; then
            log_info "[Vamb] FASTA, binsplit headers, and abundance validation are deferred until upstream assembly/mapping completes."
        else
            log_error "[Vamb] Active assembly is missing or empty: ${CONTIGS}"
            return 1
        fi
    fi

    resolve_vamb_version || return 1
    build_vamb_args vamb_args
    if [[ "${DRY_RUN}" == true ]]; then
        if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
            log_info "[Vamb] Binsplitting: enabled (separator=${VAMB_BINSPLIT_SEPARATOR})."
        else
            log_info "[Vamb] Binsplitting: disabled."
        fi
        printf '[DRY-RUN] Vamb command: vamb'
        for argument in "${vamb_args[@]}"; do printf ' %q' "${argument}"; done
        printf '\n'
    fi

    if [[ "${RUN_MEGAHIT}" == true ]]; then
        local delta=$(( VAMB_INTERNAL_MIN_CONTIG_LEN - MIN_CONTIG_LEN ))
        log_info "[Vamb] Pipeline assembly filter: contigs shorter than ${MIN_CONTIG_LEN} bp removed by MEGAHIT post-filter."
        log_info "[Vamb] Internal contig minimum: ${VAMB_INTERNAL_MIN_CONTIG_LEN} bp (default -m; separate from --minfasta=${VAMB_MINFASTA}, the minimum bin FASTA size)."
        log_info "[Vamb] Difference (Vamb internal minimum minus pipeline filter): ${delta} bp."
    else
        log_info "[Vamb] External assembly is used as supplied; no pipeline contig-length filter is recorded. Vamb internal minimum: ${VAMB_INTERNAL_MIN_CONTIG_LEN} bp."
    fi
}

build_comebin_args() {
    local output_name="$1"
    local output_dir
    # shellcheck disable=SC2178
    local -n output_args="${output_name}"
    local bam_dir="${COMEBIN_BAM_DIR:-${BAM_DIR_ACTIVE}}"
    output_dir="$(realpath -m "${COMEBIN_DIR}")"
    [[ "${COVERAGE_SOURCE}" == external && -z "${COMEBIN_BAM_DIR}" ]] && \
        bam_dir="${COMEBIN_DIR}/external_bams"
    [[ -n "${bam_dir}" ]] && bam_dir="$(realpath -m "${bam_dir}")"
    output_args=(
        -a "${CONTIGS}"
        -o "${output_dir}"
        -p "${bam_dir}"
        -s "${COMEBIN_SEED}"
        -t "${THREADS}"
        -d "${COMEBIN_DEVICE}"
    )
}

write_comebin_script() {
    local script="$1" argument quoted_wrapper
    local -a comebin_args=()
    build_comebin_args comebin_args
    printf -v quoted_wrapper '%q' "${COMEBIN_SCRIPT}"
    {
        printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'comebin_args=('
        for argument in "${comebin_args[@]}"; do
            printf '    %q\n' "${argument}"
        done
        printf '%s\n' ')' "bash ${quoted_wrapper} \"\${comebin_args[@]}\""
    } > "${script}"
    chmod +x "${script}"
}

show_comebin_command() {
    local argument
    local -a comebin_args=()
    build_comebin_args comebin_args
    printf '[DRY-RUN] COMEBin command: bash %q' "${COMEBIN_SCRIPT}"
    for argument in "${comebin_args[@]}"; do
        printf ' %q' "${argument}"
    done
    printf '\n'
}

resolve_coverage_inputs() {
    local errors=0
    BAM_DIR_ACTIVE=""
    BAM_MANIFEST_ACTIVE=""
    JGI_DEPTH_FILE_ACTIVE=""
    ACTIVE_BAMS=()
    ACTIVE_BAM_INDEXES=()
    ACTIVE_BAM_SAMPLE_IDS=()
    if [[ "${RUN_MAPPING}" == true ]]; then
        if [[ -n "${BAM_DIR_INPUT}" || -n "${BAM_MANIFEST_INPUT}" || -n "${JGI_DEPTH_FILE_INPUT}" ]]; then
            log_error "Internal mapping is enabled, but external coverage products were provided. Use --no-mapping to consume external inputs."
            return 1
        fi
        COVERAGE_SOURCE=internal
        BAM_DIR_ACTIVE="${MAPPING_DIR}"
        JGI_DEPTH_FILE_ACTIVE="${BINNING_DIR}/depth_matrix.txt"
        return 0
    fi
    if [[ -n "${BAM_DIR_INPUT}" || -n "${BAM_MANIFEST_INPUT}" || -n "${JGI_DEPTH_FILE_INPUT}" ]]; then
        COVERAGE_SOURCE=external
        if [[ -n "${BAM_DIR_INPUT}" && -n "${BAM_MANIFEST_INPUT}" ]]; then
            log_error "Use either --bam-dir or --bam-manifest, not both."
            errors=$((errors + 1))
        fi
        if [[ -n "${BAM_DIR_INPUT}" ]]; then
            if [[ ! -d "${BAM_DIR_INPUT}" ]]; then
                log_error "External BAM directory not found: ${BAM_DIR_INPUT}"
                errors=$((errors + 1))
            else
                BAM_DIR_ACTIVE="$(realpath "${BAM_DIR_INPUT}")"
            fi
        fi
        if [[ -n "${BAM_MANIFEST_INPUT}" ]]; then
            if [[ ! -f "${BAM_MANIFEST_INPUT}" || ! -s "${BAM_MANIFEST_INPUT}" ]]; then
                log_error "External BAM manifest not found or empty: ${BAM_MANIFEST_INPUT}"
                errors=$((errors + 1))
            else
                BAM_MANIFEST_ACTIVE="$(realpath "${BAM_MANIFEST_INPUT}")"
            fi
        fi
        if [[ -n "${JGI_DEPTH_FILE_INPUT}" ]]; then
            if [[ ! -f "${JGI_DEPTH_FILE_INPUT}" || ! -s "${JGI_DEPTH_FILE_INPUT}" ]]; then
                log_error "External JGI depth file not found or empty: ${JGI_DEPTH_FILE_INPUT}"
                errors=$((errors + 1))
            else
                JGI_DEPTH_FILE_ACTIVE="$(realpath "${JGI_DEPTH_FILE_INPUT}")"
            fi
        fi
    else
        COVERAGE_SOURCE=none
    fi
    (( errors == 0 ))
}

resolve_active_assembly() {
    if [[ "${RUN_MEGAHIT}" == true ]]; then
        CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
    elif [[ -n "${ASSEMBLY_IN}" ]]; then
        [[ -f "${ASSEMBLY_IN}" && -s "${ASSEMBLY_IN}" ]] || {
            log_error "External assembly file not found or empty: ${ASSEMBLY_IN}"
            return 1
        }
        ASSEMBLY_IN="$(realpath "${ASSEMBLY_IN}")"
        CONTIGS="${ASSEMBLY_IN}"
    else
        CONTIGS=""
    fi
}

_samtools_validation() {
    local mode="$1" bam="$2" output="" runtime="apptainer" parent
    if command -v samtools >/dev/null 2>&1; then
        if [[ "${mode}" == quickcheck ]]; then
            if output="$(samtools quickcheck -v "${bam}" 2>&1)"; then
                [[ -z "${output}" ]] || printf '%s\n' "${output}"
                return 0
            fi
        else
            if output="$(samtools view -H "${bam}" 2>&1)"; then
                printf '%s\n' "${output}"
                return 0
            fi
        fi
        log_warn "Host samtools could not validate ${bam}; trying the auxiliary mapping image."
    fi
    if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
        log_error "External BAM validation requires samtools on the host or the mapping container runtime."
        return 1
    fi
    command -v apptainer >/dev/null 2>&1 || runtime=singularity
    parent="$(dirname "$(realpath "${bam}")")"
    if [[ "${DRY_RUN}" == true && ! -f "${SIF_MAPPING}" ]]; then
        log_error "Cannot deeply validate external BAM during dry-run: host samtools failed and the mapping image is unavailable."
        return 1
    fi
    ensure_sif_image "${SIF_MAPPING}" "${SIF_MAPPING_URI}" || return 1
    log_info "[Validation] Using mapping image as auxiliary samtools validator; 04_mapping remains disabled."
    if [[ "${mode}" == quickcheck ]]; then
        output="$("${runtime}" exec --cleanenv -B "${parent}:${parent}" "${SIF_MAPPING}" samtools quickcheck -v "${bam}")"
    else
        output="$("${runtime}" exec --cleanenv -B "${parent}:${parent}" "${SIF_MAPPING}" samtools view -H "${bam}")"
    fi
    printf '%s\n' "${output}"
}

validate_external_bams() {
    local sample bam index line manifest_dir
    local -A seen_sample=() seen_path=()
    local -a required_samples=("${ACTIVE_SAMPLE_IDS[@]}")
    local -a rows=()
    [[ "${COVERAGE_SOURCE}" == external && ( -n "${BAM_DIR_ACTIVE}" || -n "${BAM_MANIFEST_ACTIVE}" ) ]] || return 0
    [[ "${RUN_BINNING}" == true ]] && binner_selected comebin || return 0
    if [[ -n "${BAM_MANIFEST_ACTIVE}" ]]; then
        manifest_dir="$(dirname "${BAM_MANIFEST_ACTIVE}")"
        mapfile -t rows < <(tail -n +2 "${BAM_MANIFEST_ACTIVE}")
        local header_sample header_bam header_index header_extra
        IFS=$'\t' read -r header_sample header_bam header_index header_extra < "${BAM_MANIFEST_ACTIVE}"
        if [[ "${header_sample}" != sample_id || "${header_bam}" != bam || \
              "${header_index}" != index || -n "${header_extra}" ]] || \
           ! awk -F '\t' 'NR == 1 { next } NF != 3 || $1 == "" || $2 == "" || $3 == "" { exit 1 }' \
                "${BAM_MANIFEST_ACTIVE}"; then
            log_error "BAM manifest header must begin with sample_id<TAB>bam<TAB>index."
            return 1
        fi
        for line in "${rows[@]}"; do
            IFS=$'\t' read -r sample bam index _ <<< "${line}"
            if [[ -z "${sample}" || -z "${bam}" || -z "${index}" ]]; then
                log_error "Malformed BAM manifest row: ${line}"
                return 1
            fi
            [[ "${bam}" = /* ]] || bam="${manifest_dir}/${bam}"
            [[ "${index}" = /* ]] || index="${manifest_dir}/${index}"
            bam="$(realpath "${bam}" 2>/dev/null)" || { log_error "BAM not found: ${bam}"; return 1; }
            index="$(realpath "${index}" 2>/dev/null)" || { log_error "BAM index not found: ${index}"; return 1; }
            [[ "${index}" == *.bai || "${index}" == *.csi ]] || {
                log_error "BAM index must use .bai or .csi: ${index}"
                return 1
            }
            if [[ -n "${seen_sample[$sample]:-}" || -n "${seen_path[$bam]:-}" || -n "${seen_path[$index]:-}" ]]; then
                log_error "Duplicate sample or path in BAM manifest: ${line}"
                return 1
            fi
            seen_sample["${sample}"]=1; seen_path["${bam}"]=1; seen_path["${index}"]=1
            ACTIVE_BAM_SAMPLE_IDS+=("${sample}")
            ACTIVE_BAMS+=("${bam}")
            ACTIVE_BAM_INDEXES+=("${index}")
        done
        for sample in "${required_samples[@]}"; do
            [[ -n "${seen_sample[$sample]:-}" ]] || {
                log_error "BAM manifest has no BAM for active sample ${sample}."
                return 1
            }
        done
        if (( ${#required_samples[@]} > 0 && ${#ACTIVE_BAMS[@]} != ${#required_samples[@]} )); then
            log_error "BAM manifest contains extra or missing sample BAMs."
            return 1
        fi
    else
        local -a files=()
        mapfile -d '' -t files < <(find "${BAM_DIR_ACTIVE}" -maxdepth 1 -type f -name '*.sorted.bam' -print0 | sort -z)
        if (( ${#files[@]} == 0 )); then
            log_error "No external <sample_id>.sorted.bam files found in ${BAM_DIR_ACTIVE}."
            return 1
        fi
        for bam in "${files[@]}"; do
            sample="$(basename "${bam}" .sorted.bam)"
            if [[ -n "${seen_sample[$sample]:-}" ]]; then
                log_error "Ambiguous external BAM association for sample ${sample}."
                return 1
            fi
            seen_sample["${sample}"]=1
            index=""
            if [[ -s "${bam}.bai" ]]; then index="${bam}.bai"
            elif [[ -s "${bam}.csi" ]]; then index="${bam}.csi"
            else
                log_error "External BAM has no BAI/CSI index: ${bam}"
                return 1
            fi
            ACTIVE_BAM_SAMPLE_IDS+=("${sample}")
            ACTIVE_BAMS+=("$(realpath "${bam}")")
            ACTIVE_BAM_INDEXES+=("$(realpath "${index}")")
        done
        for sample in "${required_samples[@]}"; do
            [[ -n "${seen_sample[$sample]:-}" ]] || {
                log_error "External BAM directory has no unambiguous BAM for active sample ${sample}."
                return 1
            }
        done
        if (( ${#required_samples[@]} > 0 && ${#ACTIVE_BAMS[@]} != ${#required_samples[@]} )); then
            log_error "External BAM directory contains extra or missing sample BAMs."
            return 1
        fi
    fi
    for (( line=0; line<${#ACTIVE_BAMS[@]}; line++ )); do
        bam="${ACTIVE_BAMS[$line]}"; index="${ACTIVE_BAM_INDEXES[$line]}"
        [[ -s "${bam}" && -s "${index}" ]] || { log_error "External BAM or index is empty: ${bam}"; return 1; }
        _samtools_validation quickcheck "${bam}" || { log_error "External BAM failed samtools quickcheck: ${bam}"; return 1; }
        local bam_header
        bam_header="$(_samtools_validation header "${bam}")" || return 1
        grep -Eq '^@HD.*SO:coordinate' <<< "${bam_header}" || {
            log_error "External BAM must be coordinate-sorted: ${bam}"
            return 1
        }
    done
    if [[ -z "${BAM_MANIFEST_INPUT}" ]]; then
        BAM_MANIFEST_ACTIVE="${OUTDIR}/logs/external_bam_manifest.tsv"
    fi
    log_warn "External BAM files cannot be proven automatically to match the active assembly. Compatibility is the user's responsibility."
}

persist_external_bam_manifest() {
    local index
    [[ -n "${ACTIVE_BAMS[*]:-}" && "${DRY_RUN}" != true ]] || return 0
    [[ -n "${BAM_MANIFEST_INPUT}" ]] && return 0
    mkdir -p "${LOGS_DIR}"
    {
        printf 'sample_id\tbam\tindex\n'
        for (( index=0; index<${#ACTIVE_BAMS[@]}; index++ )); do
            printf '%s\t%s\t%s\n' "${ACTIVE_BAM_SAMPLE_IDS[$index]}" \
                "${ACTIVE_BAMS[$index]}" "${ACTIVE_BAM_INDEXES[$index]}"
        done
    } > "${BAM_MANIFEST_ACTIVE}"
}

validate_external_jgi_depth() {
    local file="${JGI_DEPTH_FILE_ACTIVE}" assembly="${CONTIGS}" shared=0
    [[ "${COVERAGE_SOURCE}" == external && -n "${file}" ]] || return 0
    [[ -s "${file}" ]] || { log_error "External JGI depth file is empty: ${file}"; return 1; }
    awk -F '\t' '
        NR == 1 {
            if (NF < 4 || $1 != "contigName" || $2 != "contigLen" || $3 != "totalAvgDepth") exit 2
            columns=NF
            next
        }
        NF != columns || $1 == "" || seen[$1]++ { exit 3 }
        $2 !~ /^[0-9]+([.][0-9]+)?$/ || $3 !~ /^[0-9]+([.][0-9]+)?$/ { exit 4 }
        { for (i=4; i<=NF; i++) if ($i !~ /^[-+]?[0-9]*([.][0-9]+)?([eE][-+]?[0-9]+)?$/) exit 5 }
        END { if (NR < 2) exit 6 }
    ' "${file}" || { log_error "Malformed JGI depth table (header, duplicate IDs, columns, or numeric values): ${file}"; return 1; }
    if [[ ! -s "${assembly}" ]]; then
        log_info "[Preflight] JGI depth contig comparison deferred until the active assembly exists."
        return 0
    fi
    local missing_count shared
    missing_count="$(comm -23 <(cut -f1 "${file}" | tail -n +2 | sort -u) \
        <(awk '/^>/ { id=substr($0,2); sub(/[[:space:]].*/,"",id); print id }' "${assembly}" | sort -u) | wc -l)"
    if (( missing_count > 0 )); then
        log_error "JGI depth contains contig IDs absent from the active assembly."
        return 1
    fi
    local uncovered_count
    uncovered_count="$(comm -23 <(awk '/^>/ { id=substr($0,2); sub(/[[:space:]].*/,"",id); print id }' "${assembly}" | sort -u) \
        <(cut -f1 "${file}" | tail -n +2 | sort -u) | wc -l)"
    if (( uncovered_count > 0 )); then
        log_error "External JGI depth is missing active assembly contigs required by Vamb/MetaBinner."
        return 1
    fi
    shared="$(comm -12 <(cut -f1 "${file}" | tail -n +2 | sort -u) \
        <(awk '/^>/ { id=substr($0,2); sub(/[[:space:]].*/,"",id); print id }' "${assembly}" | sort -u) | wc -l)"
    if (( shared == 0 )); then
        log_error "External JGI depth and active assembly share zero contig IDs."
        return 1
    fi
}

classify_stages() {
    local id state
    for id in "${STAGE_IDS[@]}"; do
        state=disabled
        case "${id}" in
            01_fastp) [[ "${RUN_FASTP}" == true ]] && state=active ;;
            03_kraken2) [[ "${RUN_KRAKEN2}" == true ]] && state=active ;;
            02_megahit) [[ "${RUN_MEGAHIT}" == true ]] && state=active ;;
            04_mapping) [[ "${RUN_MAPPING}" == true ]] && state=active ;;
            05a_vamb) [[ "${RUN_BINNING}" == true ]] && binner_selected vamb && state=active ;;
            05b_metabinner) [[ "${RUN_BINNING}" == true ]] && binner_selected metabinner && state=active ;;
            05c_comebin) [[ "${RUN_BINNING}" == true ]] && binner_selected comebin && state=active ;;
            06_refinement)
                if [[ "${RUN_BINNING}" != true ]]; then state=disabled
                elif [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]]; then state=active
                else state=not-applicable; fi ;;
            07_reassembly) [[ "${RUN_REASSEMBLY}" == true ]] && state=active ;;
            08_checkm2) [[ "${RUN_CHECKM2}" == true ]] && state=active ;;
            09_filter_mags) [[ "${RUN_MAG_FILTERING}" == true ]] && state=active ;;
            10_gtdbtk) [[ "${RUN_GTDBTK}" == true ]] && state=active ;;
            11_coverm) [[ "${RUN_COVERM}" == true ]] && state=active ;;
            12a_taxonomy_summary)
                if [[ "${RUN_SUMMARY}" != true || "${RUN_TAXONOMY_SUMMARY}" != true ]]; then state=disabled
                elif [[ "${RUN_KRAKEN2}" == true ]]; then state=active
                else state=not-applicable; fi ;;
            12b_mag_summary)
                if [[ "${RUN_SUMMARY}" != true || "${RUN_MAG_SUMMARY}" != true ]]; then state=disabled
                elif [[ "${RUN_CHECKM2}" == true || "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || "${RUN_COVERM}" == true ]]; then state=active
                else state=not-applicable; fi ;;
            12c_run_overview)
                if [[ "${RUN_SUMMARY}" != true || "${RUN_RUN_OVERVIEW}" != true ]]; then state=disabled
                else state=active; fi ;;
        esac
        STAGE_STATE["${id}"]="${state}"
        STAGE_ACTION["${id}"]="${state}"
        STAGE_REASON["${id}"]="${state}"
        STAGE_FINAL_STATUS["${id}"]="${state}"
    done
}

build_effective_dependency_graph() {
    local id binner source_id
    GRAPH_CHILDREN=()
    GRAPH_PARENTS=()
    _graph_add_edge() {
        local upstream="$1" downstream="$2"
        stage_is_active "${upstream}" && stage_is_active "${downstream}" || return 0
        [[ " ${GRAPH_CHILDREN[$upstream]:-} " == *" ${downstream} "* ]] && return 0
        GRAPH_CHILDREN["${upstream}"]+="${downstream} "
        GRAPH_PARENTS["${downstream}"]+="${upstream} "
    }

    if [[ "${ACTIVE_READS_SOURCE}" == fastp-clean ]]; then
        [[ "${RUN_KRAKEN2}" == true ]] && _graph_add_edge 01_fastp 03_kraken2
        [[ "${RUN_MEGAHIT}" == true ]] && _graph_add_edge 01_fastp 02_megahit
        [[ "${RUN_MAPPING}" == true && "${COVERAGE_SOURCE}" == internal ]] && _graph_add_edge 01_fastp 04_mapping
        [[ "${RUN_REASSEMBLY}" == true ]] && _graph_add_edge 01_fastp 07_reassembly
        [[ "${RUN_COVERM}" == true ]] && _graph_add_edge 01_fastp 11_coverm
    fi
    [[ "${RUN_MEGAHIT}" == true && "${RUN_MAPPING}" == true ]] && _graph_add_edge 02_megahit 04_mapping
    if [[ "${RUN_BINNING}" == true ]]; then
        for binner in "${SELECTED_BINNER_ARRAY[@]}"; do
            id="$(resolve_stage_id "${binner}")"
            if [[ "${RUN_MAPPING}" == true && "${COVERAGE_SOURCE}" == internal ]] && \
                { binner_requires_input "${binner}" bam || binner_requires_input "${binner}" jgi-depth; }; then
                _graph_add_edge 04_mapping "${id}"
            elif [[ "${RUN_MEGAHIT}" == true ]]; then
                _graph_add_edge 02_megahit "${id}"
            fi
            [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] && _graph_add_edge "${id}" 06_refinement
        done
    fi
    if [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]]; then
        if [[ "${RUN_REASSEMBLY}" == true ]]; then
            _graph_add_edge 06_refinement 07_reassembly
        else
            _graph_add_edge 06_refinement 08_checkm2
        fi
    else
        case "${BINS_SOURCE_RESOLVED}" in
            vamb|metabinner|comebin)
                source_id="$(resolve_stage_id "${BINS_SOURCE_RESOLVED}")"
                if [[ "${RUN_REASSEMBLY}" == true ]]; then _graph_add_edge "${source_id}" 07_reassembly
                else _graph_add_edge "${source_id}" 08_checkm2; fi ;;
        esac
    fi
    if [[ "${RUN_REASSEMBLY}" == true ]]; then
        _graph_add_edge 07_reassembly 08_checkm2
    fi
    [[ "${RUN_CHECKM2}" == true ]] && [[ "${RUN_MAG_FILTERING}" == true ]] && _graph_add_edge 08_checkm2 09_filter_mags
    [[ "${RUN_MAG_FILTERING}" == true && "${RUN_GTDBTK}" == true ]] && _graph_add_edge 09_filter_mags 10_gtdbtk
    [[ "${RUN_MAG_FILTERING}" == true && "${RUN_COVERM}" == true ]] && _graph_add_edge 09_filter_mags 11_coverm
    if [[ "${RUN_KRAKEN2}" == true && "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true ]]; then
        _graph_add_edge 03_kraken2 12a_taxonomy_summary
    fi
    if [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true ]]; then
        [[ "${RUN_CHECKM2}" == true ]] && _graph_add_edge 08_checkm2 12b_mag_summary
        [[ "${RUN_MAG_FILTERING}" == true ]] && _graph_add_edge 09_filter_mags 12b_mag_summary
        [[ "${RUN_GTDBTK}" == true ]] && _graph_add_edge 10_gtdbtk 12b_mag_summary
        [[ "${RUN_COVERM}" == true ]] && _graph_add_edge 11_coverm 12b_mag_summary
    fi
    if [[ "${RUN_SUMMARY}" == true && "${RUN_RUN_OVERVIEW}" == true ]]; then
        stage_is_active 12a_taxonomy_summary && _graph_add_edge 12a_taxonomy_summary 12c_run_overview
        stage_is_active 12b_mag_summary && _graph_add_edge 12b_mag_summary 12c_run_overview
    fi
}

get_direct_dependents() {
    printf '%s\n' "${GRAPH_CHILDREN[$1]:-}" | xargs -r -n1
}

get_transitive_dependents() {
    local root="$1" node child
    local -a queue=() result=()
    local -A visited=()
    read -r -a queue <<< "${GRAPH_CHILDREN[$root]:-}"
    while (( ${#queue[@]} > 0 )); do
        node="${queue[0]}"; queue=("${queue[@]:1}")
        [[ -n "${visited[$node]:-}" ]] && continue
        visited["$node"]=1
        result+=("$node")
        for child in ${GRAPH_CHILDREN[$node]:-}; do queue+=("$child"); done
    done
    printf '%s\n' "${result[@]}"
}

stage_depends_on() {
    local downstream="$1" upstream="$2"
    [[ " ${GRAPH_PARENTS[$downstream]:-} " == *" ${upstream} "* ]]
}

assert_plan_mutable() {
    local helper_name="$1"
    if [[ "${COLLECTING_PLAN:-false}" != true || "${PLAN_IS_FROZEN:-false}" == true ]]; then
        echo "[ERROR] ${helper_name} may only be called inside collect_execution_plan before freeze." >&2
        exit 3
    fi
}

emit_plan_tsv() {
    local id
    for id in "${STAGE_IDS[@]}"; do
        local in_by diag deps down
        in_by="$(echo "${STAGE_INVALIDATED_BY[$id]:-}" | tr ' ' ',' | sed -e 's/^,//' -e 's/,$//' -e 's/,,*/,/g')"
        diag="$(echo "${STAGE_DIAGNOSTIC_CODES[$id]:-}" | tr ' ' ',' | sed -e 's/^,//' -e 's/,$//' -e 's/,,*/,/g')"
        deps="$(echo "${STAGE_DEPENDENCIES[$id]:-}" | tr ' ' ',' | sed -e 's/^,//' -e 's/,$//' -e 's/,,*/,/g')"
        down="$(echo "${STAGE_DOWNSTREAM[$id]:-}" | tr ' ' ',' | sed -e 's/^,//' -e 's/,$//' -e 's/,,*/,/g')"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${id}" \
            "${STAGE_STAGE_NAME[$id]:-$(get_stage_public_name "$id")}" \
            "${STAGE_CONFIGURED_ENABLED[$id]:-false}" \
            "${STAGE_REQUIRED[$id]:-false}" \
            "${STAGE_REQUESTED_ACTION[$id]:-normal}" \
            "${STAGE_APPLICABILITY[$id]:-not-applicable}" \
            "${STAGE_ACTION[$id]:-disabled}" \
            "${STAGE_CHECKPOINT_PATH[$id]:-${CHECKPOINT_DIR}/${id}.done}" \
            "${STAGE_CHECKPOINT_STATUS[$id]:-not-applicable}" \
            "${STAGE_CHECKPOINT_REASON[$id]:-}" \
            "${STAGE_INPUT_STATUS[$id]:-not-applicable}" \
            "${STAGE_RUNTIME_STATUS[$id]:-not-applicable}" \
            "${STAGE_DATABASE_STATUS[$id]:-not-applicable}" \
            "${STAGE_OUTPUT_STATUS[$id]:-not-applicable}" \
            "${STAGE_VALIDATION_STATUS[$id]:-not-applicable}" \
            "${in_by}" \
            "${STAGE_BLOCKING_STAGE[$id]:-}" \
            "${STAGE_REASON_CODE[$id]:-stage-disabled}" \
            "${STAGE_REASON_DETAIL[$id]:-}" \
            "${diag}" \
            "${deps}" \
            "${down}" \
            "${STAGE_PRIMARY_OUTPUT[$id]:-}" \
            "${STAGE_PRIMARY_LOG[$id]:-}"
    done
}

serialize_execution_plan_json() {
    emit_plan_tsv | python3 "${PIPELINE_SCRIPTS_DIR}/execution_plan.py" assemble --profile "${REQUESTED_PROFILE:-full}"
}

resolve_base_stage_actions() {
    assert_plan_mutable "resolve_base_stage_actions"
    check_legacy_summary_checkpoint
    classify_stages
    build_effective_dependency_graph

    local id pub req_act
    for id in "${STAGE_IDS[@]}"; do
        pub="$(get_stage_public_name "${id}")"
        STAGE_STAGE_NAME["${id}"]="${pub}"
        STAGE_CHECKPOINT_PATH["${id}"]="${CHECKPOINT_DIR}/${id}.done"
        STAGE_PRIMARY_OUTPUT["${id}"]="$(stage_primary_output "${id}")"
        STAGE_PRIMARY_LOG["${id}"]="$(stage_primary_log "${id}")"
        STAGE_INVALIDATED_BY["${id}"]=""
        STAGE_BLOCKING_STAGE["${id}"]=""
        STAGE_DIAGNOSTIC_CODES["${id}"]=""
        STAGE_REASON_DETAIL["${id}"]=""

        STAGE_DEPENDENCIES["${id}"]="${GRAPH_PARENTS[$id]:-}"
        STAGE_DOWNSTREAM["${id}"]="${GRAPH_CHILDREN[$id]:-}"

        case "${id}" in
            01_fastp) [[ "${RUN_FASTP}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            03_kraken2) [[ "${RUN_KRAKEN2}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            02_megahit) [[ "${RUN_MEGAHIT}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            04_mapping) [[ "${RUN_MAPPING}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            05a_vamb) [[ "${RUN_BINNING}" == true ]] && binner_selected vamb && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            05b_metabinner) [[ "${RUN_BINNING}" == true ]] && binner_selected metabinner && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            05c_comebin) [[ "${RUN_BINNING}" == true ]] && binner_selected comebin && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            06_refinement) [[ "${RUN_BINNING}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            07_reassembly) [[ "${RUN_REASSEMBLY}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            08_checkm2) [[ "${RUN_CHECKM2}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            09_filter_mags) [[ "${RUN_MAG_FILTERING}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            10_gtdbtk) [[ "${RUN_GTDBTK}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            11_coverm) [[ "${RUN_COVERM}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            12a_taxonomy_summary) [[ "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            12b_mag_summary) [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
            12c_run_overview) [[ "${RUN_SUMMARY}" == true && "${RUN_RUN_OVERVIEW}" == true ]] && STAGE_CONFIGURED_ENABLED["${id}"]=true || STAGE_CONFIGURED_ENABLED["${id}"]=false ;;
        esac

        req_act="normal"
        if [[ "${RESTART_REQUESTED}" == true && "${RESTART_FROM_NAME}" == "${pub}" ]]; then
            req_act="restart"
        elif (( ${#FORCE_STAGE_NAMES[@]} > 0 )); then
            local f_name
            for f_name in "${FORCE_STAGE_NAMES[@]}"; do
                if [[ "${f_name}" == "${pub}" ]]; then
                    req_act="force"
                    break
                fi
            done
        elif [[ "${RESUME_CONTROL}" == force ]]; then
            req_act="force"
        fi
        STAGE_REQUESTED_ACTION["${id}"]="${req_act}"

        if [[ "${STAGE_CONFIGURED_ENABLED[$id]}" == false ]]; then
            STAGE_APPLICABILITY["${id}"]="not-applicable"
            STAGE_REQUIRED["${id}"]=false
            STAGE_ACTION["${id}"]="disabled"
            STAGE_REASON_CODE["${id}"]="stage-disabled"
            STAGE_REASON_DETAIL["${id}"]="Stage is disabled in configuration"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-applicable"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_INPUT_STATUS["${id}"]="not-applicable"
            STAGE_RUNTIME_STATUS["${id}"]="not-applicable"
            STAGE_DATABASE_STATUS["${id}"]="not-applicable"
            STAGE_OUTPUT_STATUS["${id}"]="not-applicable"
            STAGE_VALIDATION_STATUS["${id}"]="not-applicable"
            STAGE_DEPENDENCIES["${id}"]=""
            STAGE_DOWNSTREAM["${id}"]=""
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        if ! stage_in_scope "${pub}"; then
            STAGE_APPLICABILITY["${id}"]="not-applicable"
            STAGE_REQUIRED["${id}"]=false
            STAGE_ACTION["${id}"]="not-applicable"
            STAGE_REASON_CODE["${id}"]="stage-not-applicable"
            STAGE_REASON_DETAIL["${id}"]="Stage is outside the scope of --stop-after ${STOP_AFTER}"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-applicable"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_INPUT_STATUS["${id}"]="not-applicable"
            STAGE_RUNTIME_STATUS["${id}"]="not-applicable"
            STAGE_DATABASE_STATUS["${id}"]="not-applicable"
            STAGE_OUTPUT_STATUS["${id}"]="not-applicable"
            STAGE_VALIDATION_STATUS["${id}"]="not-applicable"
            STAGE_DEPENDENCIES["${id}"]=""
            STAGE_DOWNSTREAM["${id}"]=""
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        if [[ "${id}" == 06_refinement && "${RUN_REFINEMENT_EFFECTIVE}" != true ]]; then
            STAGE_APPLICABILITY["${id}"]="not-applicable"
            STAGE_REQUIRED["${id}"]=false
            STAGE_ACTION["${id}"]="not-applicable"
            STAGE_REASON_CODE["${id}"]="stage-not-applicable"
            STAGE_REASON_DETAIL["${id}"]="Refinement is not effective with selected binners"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-applicable"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_INPUT_STATUS["${id}"]="not-applicable"
            STAGE_RUNTIME_STATUS["${id}"]="not-applicable"
            STAGE_DATABASE_STATUS["${id}"]="not-applicable"
            STAGE_OUTPUT_STATUS["${id}"]="not-applicable"
            STAGE_VALIDATION_STATUS["${id}"]="not-applicable"
            STAGE_DEPENDENCIES["${id}"]=""
            STAGE_DOWNSTREAM["${id}"]=""
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        if [[ "${id}" == 12a_taxonomy_summary && "${RUN_KRAKEN2}" != true ]]; then
            STAGE_APPLICABILITY["${id}"]="not-applicable"
            STAGE_REQUIRED["${id}"]=false
            STAGE_ACTION["${id}"]="not-applicable"
            STAGE_REASON_CODE["${id}"]="stage-not-applicable"
            STAGE_REASON_DETAIL["${id}"]="Kraken2 is disabled; taxonomy summary not applicable"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-applicable"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_INPUT_STATUS["${id}"]="not-applicable"
            STAGE_RUNTIME_STATUS["${id}"]="not-applicable"
            STAGE_DATABASE_STATUS["${id}"]="not-applicable"
            STAGE_OUTPUT_STATUS["${id}"]="not-applicable"
            STAGE_VALIDATION_STATUS["${id}"]="not-applicable"
            STAGE_DEPENDENCIES["${id}"]=""
            STAGE_DOWNSTREAM["${id}"]=""
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        if [[ "${id}" == 12b_mag_summary && "${RUN_CHECKM2}" != true && "${RUN_MAG_FILTERING}" != true && "${RUN_GTDBTK}" != true && "${RUN_COVERM}" != true ]]; then
            STAGE_APPLICABILITY["${id}"]="not-applicable"
            STAGE_REQUIRED["${id}"]=false
            STAGE_ACTION["${id}"]="not-applicable"
            STAGE_REASON_CODE["${id}"]="stage-not-applicable"
            STAGE_REASON_DETAIL["${id}"]="No MAG analytical stages are enabled; MAG summary not applicable"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-applicable"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_INPUT_STATUS["${id}"]="not-applicable"
            STAGE_RUNTIME_STATUS["${id}"]="not-applicable"
            STAGE_DATABASE_STATUS["${id}"]="not-applicable"
            STAGE_OUTPUT_STATUS["${id}"]="not-applicable"
            STAGE_VALIDATION_STATUS["${id}"]="not-applicable"
            STAGE_DEPENDENCIES["${id}"]=""
            STAGE_DOWNSTREAM["${id}"]=""
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        STAGE_APPLICABILITY["${id}"]="applicable"
        STAGE_REQUIRED["${id}"]=true

        local in_st="available" run_st="available" db_st="not-applicable" out_st="available"
        local blk_cause="" blk_detail=""

        case "${id}" in
            01_fastp)
                if [[ -z "${RAW_READS_DIR:-}" ]]; then
                    in_st="not-configured"; blk_cause="external-input-missing"; blk_detail="Raw reads directory (--reads-dir) not configured"
                elif [[ ! -d "${RAW_READS_DIR}" ]]; then
                    in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="Raw reads directory does not exist: ${RAW_READS_DIR}"
                elif [[ $(find "${RAW_READS_DIR}" -maxdepth 1 \( -name '*.fastq*' -o -name '*.fq*' \) 2>/dev/null | head -n 1 | wc -l) -eq 0 ]]; then
                    in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="No fastq/fq reads found in ${RAW_READS_DIR}"
                fi ;;
            02_megahit)
                if [[ "${RUN_FASTP}" != true ]]; then
                    if [[ -n "${CLEAN_READS_DIR:-}" ]]; then
                        if [[ ! -d "${CLEAN_READS_DIR}" ]]; then in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="Clean reads directory does not exist: ${CLEAN_READS_DIR}"; fi
                    elif [[ -n "${RAW_READS_DIR:-}" ]]; then
                        if [[ ! -d "${RAW_READS_DIR}" ]]; then in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="Raw reads directory does not exist: ${RAW_READS_DIR}"; fi
                    else
                        in_st="not-configured"; blk_cause="external-input-missing"; blk_detail="No reads directory provided"
                    fi
                fi ;;
            04_mapping)
                if [[ "${COVERAGE_SOURCE:-internal}" == external ]]; then
                    if [[ -n "${JGI_DEPTH_FILE:-}" && ! -s "${JGI_DEPTH_FILE}" ]]; then
                        in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="JGI depth file not found or empty: ${JGI_DEPTH_FILE}"
                    elif [[ -n "${BAM_DIR:-}" && ! -d "${BAM_DIR}" ]]; then
                        in_st="unavailable"; blk_cause="external-input-missing"; blk_detail="BAM directory not found: ${BAM_DIR}"
                    fi
                fi ;;
        esac

        if [[ "${EXEC_MODE}" == apptainer ]]; then
            if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
                run_st="unavailable"
                [[ -z "${blk_cause}" ]] && { blk_cause="runtime-unavailable"; blk_detail="Apptainer or Singularity executable not found on host"; }
            else
                local sif_img=""
                case "${id}" in
                    01_fastp) sif_img="fastp.sif" ;;
                    02_megahit) sif_img="megahit.sif" ;;
                    03_kraken2) sif_img="kraken2_bracken.sif" ;;
                    04_mapping) sif_img="mapping_depth.sif" ;;
                    05a_vamb) sif_img="vamb.sif" ;;
                    05b_metabinner) sif_img="metabinner.sif" ;;
                    05c_comebin) sif_img="comebin.sif" ;;
                    06_refinement|07_reassembly) sif_img="metawrap.sif" ;;
                    08_checkm2) sif_img="checkm2.sif" ;;
                    10_gtdbtk) sif_img="gtdbtk.sif" ;;
                    11_coverm) sif_img="coverm.sif" ;;
                esac
                if [[ -n "${sif_img}" ]]; then
                    local s_path="${SIF_DIR:-./containers}/${sif_img}"
                    if [[ ! -f "${s_path}" || ! -r "${s_path}" ]]; then
                        run_st="unavailable"
                        [[ -z "${blk_cause}" ]] && { blk_cause="runtime-unavailable"; blk_detail="Container image '${sif_img}' not found at ${s_path}"; }
                    fi
                fi
            fi
        fi

        case "${id}" in
            03_kraken2)
                if [[ -z "${KRAKEN2_DB:-}" ]]; then
                    db_st="not-configured"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="Kraken2 database path is not configured"; }
                elif [[ ! -d "${KRAKEN2_DB}" ]]; then
                    db_st="unavailable"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="Kraken2 database directory does not exist: ${KRAKEN2_DB}"; }
                elif [[ ! -f "${KRAKEN2_DB}/hash.k2d" || ! -f "${KRAKEN2_DB}/opts.k2d" || ! -f "${KRAKEN2_DB}/taxo.k2d" ]]; then
                    db_st="invalid"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="Kraken2 database missing required keyfiles"; }
                else
                    db_st="available"
                fi ;;
            08_checkm2)
                if [[ -z "${CHECKM2_DB:-}" ]]; then
                    db_st="not-configured"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="CheckM2 diamond database path is not configured"; }
                elif [[ ! -f "${CHECKM2_DB}" || ! -s "${CHECKM2_DB}" ]]; then
                    db_st="unavailable"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="CheckM2 diamond database file not found: ${CHECKM2_DB}"; }
                else
                    db_st="available"
                fi ;;
            10_gtdbtk)
                if [[ -z "${GTDBTK_DATA_PATH:-}" ]]; then
                    db_st="not-configured"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="GTDB-Tk reference data path is not configured"; }
                elif [[ ! -d "${GTDBTK_DATA_PATH}" ]]; then
                    db_st="unavailable"
                    [[ -z "${blk_cause}" ]] && { blk_cause="database-not-configured"; blk_detail="GTDB-Tk reference data directory does not exist: ${GTDBTK_DATA_PATH}"; }
                else
                    db_st="available"
                fi ;;
        esac

        local cur_p="${OUTDIR}"
        while [[ ! -d "${cur_p}" && "${cur_p}" != "/" && "${cur_p}" != "." ]]; do
            cur_p="$(dirname "${cur_p}")"
        done
        if [[ ! -w "${cur_p}" ]]; then
            out_st="invalid"
            [[ -z "${blk_cause}" ]] && { blk_cause="output-parent-not-writable"; blk_detail="Nearest parent directory '${cur_p}' is not writable"; }
        else
            if [[ -d "${OUTDIR}" && -e "${STAGE_PRIMARY_OUTPUT[$id]}" ]]; then
                out_st="available"
            else
                out_st="missing"
            fi
        fi

        STAGE_INPUT_STATUS["${id}"]="${in_st}"
        STAGE_RUNTIME_STATUS["${id}"]="${run_st}"
        STAGE_DATABASE_STATUS["${id}"]="${db_st}"
        STAGE_OUTPUT_STATUS["${id}"]="${out_st}"

        if [[ -n "${blk_cause}" ]]; then
            STAGE_ACTION["${id}"]="blocked"
            STAGE_REASON_CODE["${id}"]="${blk_cause}"
            STAGE_REASON_DETAIL["${id}"]="${blk_detail}"
            STAGE_CHECKPOINT_STATUS["${id}"]="not-inspected"
            STAGE_VALIDATION_STATUS["${id}"]="not-performed"
            if [[ "${req_act}" == force ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            elif [[ "${req_act}" == restart ]]; then
                STAGE_DIAGNOSTIC_CODES["${id}"]="restart-boundary"
            fi
            continue
        fi

        local chk="${CHECKPOINT_DIR}/${id}.done"
        local cfg="${CHECKPOINT_DIR}/${id}.config"
        if [[ ! -d "${OUTDIR}" || ! -f "${chk}" ]]; then
            STAGE_CHECKPOINT_STATUS["${id}"]="missing"
            STAGE_CHECKPOINT_REASON["${id}"]=""
            STAGE_VALIDATION_STATUS["${id}"]="not-performed"
            STAGE_ACTION["${id}"]="run"
            STAGE_REASON_CODE["${id}"]="checkpoint-missing"
            STAGE_REASON_DETAIL["${id}"]="Checkpoint does not exist"
        elif [[ ! -f "${cfg}" ]]; then
            STAGE_CHECKPOINT_STATUS["${id}"]="legacy"
            STAGE_CHECKPOINT_REASON["${id}"]="manifest missing"
            STAGE_VALIDATION_STATUS["${id}"]="not-performed"
            STAGE_ACTION["${id}"]="run"
            STAGE_REASON_CODE["${id}"]="checkpoint-legacy"
            STAGE_REASON_DETAIL["${id}"]="Checkpoint manifest missing"
        else
            local cur_manifest rec_manifest
            cur_manifest="$(checkpoint_config_manifest "${id}" 2>/dev/null || true)"
            rec_manifest="$(cat "${cfg}" 2>/dev/null || true)"
            if [[ "$(comparable_checkpoint_manifest "${id}" "${cur_manifest}")" != "$(comparable_checkpoint_manifest "${id}" "${rec_manifest}")" ]] || \
               ! checkpoint_identity_compatible "${id}" "${cur_manifest}" "${rec_manifest}" 2>/dev/null; then
                STAGE_CHECKPOINT_STATUS["${id}"]="stale"
                STAGE_CHECKPOINT_REASON["${id}"]="manifest changed"
                STAGE_VALIDATION_STATUS["${id}"]="not-performed"
                STAGE_ACTION["${id}"]="run"
                STAGE_REASON_CODE["${id}"]="configuration-changed"
                STAGE_REASON_DETAIL["${id}"]="Checkpoint configuration differs from current configuration"
            elif ! validate_step_products "${id}" >/dev/null 2>&1; then
                STAGE_CHECKPOINT_STATUS["${id}"]="invalid"
                STAGE_CHECKPOINT_REASON["${id}"]="products invalid or missing"
                STAGE_VALIDATION_STATUS["${id}"]="failed"
                STAGE_ACTION["${id}"]="run"
                STAGE_REASON_CODE["${id}"]="product-validation-failed"
                STAGE_REASON_DETAIL["${id}"]="Checkpoint products failed validation or are missing"
            else
                STAGE_CHECKPOINT_STATUS["${id}"]="valid"
                STAGE_CHECKPOINT_REASON["${id}"]="valid"
                STAGE_VALIDATION_STATUS["${id}"]="valid"
                STAGE_ACTION["${id}"]="reuse"
                STAGE_REASON_CODE["${id}"]="checkpoint-valid"
                STAGE_REASON_DETAIL["${id}"]="Checkpoint manifest and products are valid"
            fi
        fi
    done
}

apply_force_stage_to_plan() {
    assert_plan_mutable "apply_force_stage_to_plan"
    local id pub node child
    local -a queue=()
    local -A visited=()

    if [[ "${RESUME_CONTROL}" == force ]]; then
        for id in "${STAGE_IDS[@]}"; do
            [[ "${STAGE_APPLICABILITY[$id]}" == applicable ]] || continue
            if [[ "${STAGE_ACTION[$id]}" != blocked ]]; then
                STAGE_ACTION["${id}"]="force"
                STAGE_REASON_CODE["${id}"]="forced-by-user"
                STAGE_REASON_DETAIL["${id}"]="Global force rerun requested"
            else
                STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
            fi
        done
        return 0
    fi

    if (( ${#FORCE_STAGE_NAMES[@]} > 0 )); then
        for pub in "${FORCE_STAGE_NAMES[@]}"; do
            id="$(resolve_stage_id "${pub}")" || exit 2
            if [[ "${STAGE_APPLICABILITY[$id]}" == applicable ]]; then
                if [[ "${STAGE_ACTION[$id]}" != blocked ]]; then
                    STAGE_ACTION["${id}"]="force"
                    STAGE_REASON_CODE["${id}"]="forced-by-user"
                    STAGE_REASON_DETAIL["${id}"]="Stage forced by user"
                    queue+=("${id}")
                else
                    STAGE_DIAGNOSTIC_CODES["${id}"]="forced-by-user"
                fi
            fi
        done

        while (( ${#queue[@]} > 0 )); do
            node="${queue[0]}"
            queue=("${queue[@]:1}")
            [[ -n "${visited[$node]:-}" ]] && continue
            visited["${node}"]=1
            for child in ${GRAPH_CHILDREN[$node]:-}; do
                [[ "${STAGE_APPLICABILITY[$child]:-}" == applicable ]] || continue
                if [[ "${STAGE_ACTION[$child]}" != blocked && "${STAGE_ACTION[$child]}" != force ]]; then
                    STAGE_ACTION["${child}"]="run"
                    STAGE_REASON_CODE["${child}"]="upstream-forced"
                    STAGE_REASON_DETAIL["${child}"]="Upstream stage ${node} was forced by user"
                    if [[ -z "${STAGE_INVALIDATED_BY[$child]}" ]]; then
                        STAGE_INVALIDATED_BY["${child}"]="${node}"
                    elif [[ " ${STAGE_INVALIDATED_BY[$child]} " != *" ${node} "* ]]; then
                        STAGE_INVALIDATED_BY["${child}"]+=" ${node}"
                    fi
                    queue+=("${child}")
                fi
            done
        done
    fi
}

apply_restart_from_to_plan() {
    assert_plan_mutable "apply_restart_from_to_plan"
    [[ "${RESTART_REQUESTED}" == true ]] || return 0

    local root node child
    local -a queue=()
    local -A visited=()

    root="$(resolve_stage_id "${RESTART_FROM_NAME}")" || exit 2
    if [[ "${STAGE_APPLICABILITY[$root]}" == applicable ]]; then
        if [[ "${STAGE_ACTION[$root]}" != blocked ]]; then
            STAGE_ACTION["${root}"]="restart"
            STAGE_REASON_CODE["${root}"]="restart-boundary"
            STAGE_REASON_DETAIL["${root}"]="Restart boundary requested by user"
            queue+=("${root}")
        else
            STAGE_DIAGNOSTIC_CODES["${root}"]="restart-boundary"
        fi
    fi

    while (( ${#queue[@]} > 0 )); do
        node="${queue[0]}"
        queue=("${queue[@]:1}")
        [[ -n "${visited[$node]:-}" ]] && continue
        visited["${node}"]=1
        for child in ${GRAPH_CHILDREN[$node]:-}; do
            [[ "${STAGE_APPLICABILITY[$child]:-}" == applicable ]] || continue
            if [[ "${STAGE_ACTION[$child]}" != blocked ]]; then
                STAGE_ACTION["${child}"]="run"
                STAGE_REASON_CODE["${child}"]="downstream-of-restart"
                STAGE_REASON_DETAIL["${child}"]="Descendant of restart boundary ${RESTART_FROM_NAME}"
                if [[ -z "${STAGE_INVALIDATED_BY[$child]}" ]]; then
                    STAGE_INVALIDATED_BY["${child}"]="${root}"
                elif [[ " ${STAGE_INVALIDATED_BY[$child]} " != *" ${root} "* ]]; then
                    STAGE_INVALIDATED_BY["${child}"]+=" ${root}"
                fi
                queue+=("${child}")
            fi
        done
    done
}

propagate_plan_invalidation() {
    assert_plan_mutable "propagate_plan_invalidation"
    local node child
    local -a queue=()
    local -A visited=()

    # 1. Propagate blockages along effective graph
    for node in "${STAGE_IDS[@]}"; do
        if [[ "${STAGE_ACTION[$node]}" == blocked ]]; then
            queue+=("${node}")
        fi
    done

    while (( ${#queue[@]} > 0 )); do
        node="${queue[0]}"
        queue=("${queue[@]:1}")
        [[ -n "${visited[$node]:-}" ]] && continue
        visited["${node}"]=1
        for child in ${GRAPH_CHILDREN[$node]:-}; do
            [[ "${STAGE_APPLICABILITY[$child]:-}" == applicable ]] || continue
            if [[ "${STAGE_ACTION[$child]}" != blocked ]]; then
                STAGE_ACTION["${child}"]="blocked"
                STAGE_REASON_CODE["${child}"]="dependency-blocked"
                STAGE_BLOCKING_STAGE["${child}"]="${node}"
                STAGE_REASON_DETAIL["${child}"]="Upstream dependency ${node} is blocked"
                queue+=("${child}")
            fi
        done
    done

    # 2. Propagate rebuilds (run/force/restart invalidates downstream reuse)
    queue=()
    visited=()
    for node in "${STAGE_IDS[@]}"; do
        if [[ "${STAGE_ACTION[$node]}" == run || "${STAGE_ACTION[$node]}" == force || "${STAGE_ACTION[$node]}" == restart ]]; then
            queue+=("${node}")
        fi
    done

    while (( ${#queue[@]} > 0 )); do
        node="${queue[0]}"
        queue=("${queue[@]:1}")
        [[ -n "${visited[$node]:-}" ]] && continue
        visited["${node}"]=1
        for child in ${GRAPH_CHILDREN[$node]:-}; do
            [[ "${STAGE_APPLICABILITY[$child]:-}" == applicable ]] || continue
            if [[ "${STAGE_ACTION[$child]}" == reuse ]]; then
                STAGE_ACTION["${child}"]="run"
                STAGE_REASON_CODE["${child}"]="upstream-invalidated"
                STAGE_REASON_DETAIL["${child}"]="Upstream stage ${node} will be rebuilt"
                if [[ -z "${STAGE_INVALIDATED_BY[$child]}" ]]; then
                    STAGE_INVALIDATED_BY["${child}"]="${node}"
                elif [[ " ${STAGE_INVALIDATED_BY[$child]} " != *" ${node} "* ]]; then
                    STAGE_INVALIDATED_BY["${child}"]+=" ${node}"
                fi
                queue+=("${child}")
            fi
        done
    done
}

validate_execution_plan() {
    assert_plan_mutable "validate_execution_plan"
    local tmp_val
    tmp_val="$(mktemp "${TMPDIR:-/tmp}/metanexus_val.XXXXXX.json")"
    if ! serialize_execution_plan_json > "${tmp_val}"; then
        rm -f "${tmp_val}"
        echo "[ERROR] Failed to serialize plan for validation" >&2
        exit 3
    fi
    local val_rc=0
    python3 "${PIPELINE_SCRIPTS_DIR}/execution_plan.py" validate --input "${tmp_val}" || val_rc=$?
    rm -f "${tmp_val}"
    if (( val_rc != 0 )); then
        exit 3
    fi
    return 0
}

freeze_execution_plan() {
    assert_plan_mutable "freeze_execution_plan"
    PLAN_IS_FROZEN=true
    local id
    for id in "${STAGE_IDS[@]}"; do
        STAGE_REASON["${id}"]="${STAGE_REASON_CODE[$id]}"
        STAGE_PRECONDITIONS["${id}"]="$(printf '%s|%s|%s|%s|%s' \
            "${STAGE_REQUESTED_ACTION[$id]}" \
            "${STAGE_ACTION[$id]}" \
            "${STAGE_CHECKPOINT_STATUS[$id]}" \
            "${STAGE_INPUT_STATUS[$id]}" \
            "${STAGE_RUNTIME_STATUS[$id]}")"
    done
}

collect_execution_plan() {
    COLLECTING_PLAN=true
    PLAN_IS_FROZEN=false

    resolve_base_stage_actions || { COLLECTING_PLAN=false; return 1; }
    apply_force_stage_to_plan || { COLLECTING_PLAN=false; return 1; }
    apply_restart_from_to_plan || { COLLECTING_PLAN=false; return 1; }
    propagate_plan_invalidation || { COLLECTING_PLAN=false; return 1; }
    validate_execution_plan || { COLLECTING_PLAN=false; return 1; }
    freeze_execution_plan || { COLLECTING_PLAN=false; return 1; }

    COLLECTING_PLAN=false
    return 0
}

verify_stage_preconditions() {
    local id="$1"
    local recorded="${STAGE_PRECONDITIONS[$id]:-}"
    [[ -n "${recorded}" ]] || return 0
    local cur
    cur="$(printf '%s|%s|%s|%s|%s' \
        "${STAGE_REQUESTED_ACTION[$id]}" \
        "${STAGE_ACTION[$id]}" \
        "${STAGE_CHECKPOINT_STATUS[$id]}" \
        "${STAGE_INPUT_STATUS[$id]}" \
        "${STAGE_RUNTIME_STATUS[$id]}")"
    if [[ "${recorded}" != "${cur}" ]]; then
        log_error "Critical preconditions changed after planning for stage ${id} (recorded: ${recorded}, current: ${cur})."
        return 1
    fi
    return 0
}

write_resume_plan() {
    local id public display action active checkpoint manifest reason
    RESUME_PLAN_FILE="${OUTDIR}/resume_plan.tsv"
    printf 'stage_internal_id\tstage_public_name\tstage_display_name\tactive\tplanned_action\treason\tcheckpoint_path\tmanifest_path\n'
    for id in "${STAGE_IDS[@]}"; do
        public="$(get_stage_public_name "${id}")"
        display="$(get_stage_display_name "${id}")"
        active="${STAGE_STATE[$id]}"
        action="${STAGE_ACTION[$id]}"
        checkpoint="${CHECKPOINT_DIR}/${id}.done"
        manifest="${CHECKPOINT_DIR}/${id}.config"
        reason="${STAGE_REASON[$id]}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${id}" "${public}" "${display}" \
            "${active}" "${action}" "${reason}" "${checkpoint}" "${manifest}"
    done
}

persist_resume_plan() {
    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[Plan] Dry-run: resume_plan.tsv will not be written."
        return 0
    fi
    mkdir -p "${OUTDIR}"
    write_resume_plan > "${OUTDIR}/resume_plan.tsv"
    RESUME_PLAN_FILE="${OUTDIR}/resume_plan.tsv"
    log_info "Resume plan written to ${RESUME_PLAN_FILE}."
}

initialize_run_status() {
    RUN_MANIFEST_FILE="${LOGS_DIR}/run_manifest.tsv"
    STAGE_STATUS_FILE="${LOGS_DIR}/stage_status.tsv"
    RUN_STATUS_FILE="${OUTDIR}/run_status.tsv"

    [[ "${DRY_RUN}" == true ]] && return 0
    mkdir -p "${LOGS_DIR}" "${OUTDIR}"

    local cpu_req="NA"
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        cpu_req="${SLURM_CPUS_PER_TASK}"
    elif [[ -n "${PBS_NP:-}" ]]; then
        cpu_req="${PBS_NP}"
    fi

    local mem_lim="NA"
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        mem_lim="${SLURM_MEM_PER_NODE}"
    fi

    local -a keep_tmp_args=()
    [[ "${KEEP_TEMPORARY}" == true ]] && keep_tmp_args=(--keep-temporary)

    python3 "${BENCHMARK_WRITER}" init \
        --run-manifest-path "${RUN_MANIFEST_FILE}" \
        --stage-status-path "${STAGE_STATUS_FILE}" \
        --compat-run-status-path "${RUN_STATUS_FILE}" \
        --resume-plan-path "${RESUME_PLAN_FILE}" \
        --checkpoint-dir "${CHECKPOINT_DIR}" \
        --run-id "${RUN_ID}" \
        --pipeline-version "${PIPELINE_VERSION}" \
        --profile "${REQUESTED_PROFILE:-custom}" \
        --execution-mode "${EXEC_MODE}" \
        --start-time-utc "${PIPELINE_START_TIME_UTC}" \
        --outdir "${OUTDIR}" \
        --command-line "${PIPELINE_INVOCATION_COMMAND}" \
        --hostname "${PIPELINE_HOSTNAME:-unknown}" \
        --cpu-requested "${cpu_req}" \
        --threads-effective "${THREADS}" \
        --memory-limit-bytes "${mem_lim}" \
        --active-sample-count "${#ACTIVE_SAMPLE_IDS[@]}" \
        --active-reads-source "${ACTIVE_READS_SOURCE:-none}" \
        --sample-manifest-path "${SAMPLE_MANIFEST_ACTIVE:-NA}" \
        --sample-manifest-fingerprint "${ACTIVE_SAMPLE_MANIFEST_FINGERPRINT:-NA}" \
        --assembly-source "$([[ -n "${ASSEMBLY_IN}" ]] && echo external || echo megahit)" \
        --assembly-path "${CONTIGS:-NA}" \
        --assembly-fingerprint "${ASSEMBLY_FINGERPRINT:-NA}" \
        --bins-source-requested "${BINS_SOURCE:-NA}" \
        --bins-source-resolved "${BINS_SOURCE_RESOLVED:-NA}" \
        --bins-fingerprint "${EXTERNAL_BINS_FINGERPRINT:-NA}" \
        --min-completeness "${MIN_COMPLETENESS}" \
        --max-contamination "${MAX_CONTAMINATION}" \
        --tmpdir-source "${PIPELINE_TMPDIR_SOURCE:-fallback}" \
        --tmpdir-base "${PIPELINE_TMPDIR_BASE:-NA}" \
        "${keep_tmp_args[@]}" \
        --software-inventory-path "${LOGS_DIR}/software_versions.tsv" \
        --database-inventory-path "${LOGS_DIR}/database_versions.tsv" \
        --error-history-path "${LOGS_DIR}/errors.tsv" \
        --last-error-path "${LOGS_DIR}/last_error.json" \
        --taxonomy-summary-path "${SUMMARY_DIR}/taxonomy_summary.tsv" \
        --mag-summary-path "${SUMMARY_DIR}/mag_summary.tsv" \
        --run-overview-path "${SUMMARY_DIR}/run_overview.tsv" \
        --warning-count "0" \
        --error-count "0"

    RUN_STATUS_INITIALIZED=true
}

update_run_status() {
    local id="$1" status="$2"
    [[ "${DRY_RUN}" == true || "${RUN_STATUS_INITIALIZED}" != true ]] && return 0
    STAGE_FINAL_STATUS["${id}"]="${status}"
    : "${STAGE_START_TIME[$id]:-}" "${STAGE_END_TIME[$id]:-}" "${STAGE_EXIT_CODE[$id]:-}"
    python3 "${BENCHMARK_WRITER}" sync-compat \
        --stage-status-path "${STAGE_STATUS_FILE}" \
        --compat-run-status-path "${RUN_STATUS_FILE}" \
        --run-id "${RUN_ID}" >/dev/null 2>&1 || true
}

assert_path_within_outdir() {
    local path="$1" base resolved
    base="$(realpath -m "${OUTDIR}")"
    resolved="$(realpath -m "${path}")"
    [[ "${resolved}" == "${base}" || "${resolved}" == "${base}/"* ]] || {
        log_error "Refusing to modify path outside OUTDIR: ${path}"
        return 1
    }
}

clean_stage_outputs() {
    local id="$1" path
    local -a paths=()
    case "${id}" in
        01_fastp) paths=("${QC_DIR}") ;;
        03_kraken2) paths=("${KRAKEN2_DIR}") ;;
        02_megahit) paths=("${ASSEMBLY_DIR}") ;;
        04_mapping) paths=("${MAPPING_DIR}" "${BINNING_DIR}/depth_matrix.txt") ;;
        05a_vamb) paths=("${VAMB_DIR}/vamb_out" "${VAMB_FINAL_BINS_DIR}" "${VAMB_ABUNDANCE_TSV}") ;;
        05b_metabinner) paths=("${METABINNER_DIR}") ;;
        05c_comebin) paths=("${COMEBIN_DIR}") ;;
        06_refinement) paths=("${REFINEMENT_DIR}") ;;
        07_reassembly) paths=("${REASSEMBLY_DIR}") ;;
        08_checkm2) paths=("${CHECKM2_DIR}") ;;
        09_filter_mags) paths=("${MAGS_DIR}") ;;
        10_gtdbtk) paths=("${GTDBTK_DIR}") ;;
        11_coverm) paths=("${ABUNDANCE_DIR}") ;;
        12a_taxonomy_summary) paths=("${SUMMARY_DIR}/taxonomy_summary.tsv" "${SUMMARY_DIR}/taxonomy_summary.meta.json") ;;
        12b_mag_summary) paths=("${SUMMARY_DIR}/mag_summary.tsv" "${SUMMARY_DIR}/mag_summary.meta.json") ;;
        12c_run_overview) paths=("${SUMMARY_DIR}/run_overview.tsv" "${SUMMARY_DIR}/run_overview.meta.json") ;;
        12_summary) paths=("${SUMMARY_DIR}/final_summary.tsv") ;;
        *) log_error "Cannot clean unknown stage ID: ${id}"; return 1 ;;
    esac
    for path in "${paths[@]}"; do
        assert_path_within_outdir "${path}" || return 1
        local protected_path
        protected_path="$(realpath -m "${path}")"
        case "${protected_path}" in
            "$(realpath -m "${LOGS_DIR}/inventory")"|\
            "$(realpath -m "${LOGS_DIR}/inventory")"/*|\
            "$(realpath -m "${LOGS_DIR}/disk_space.current.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/disk_space.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/run_manifest.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/stage_status.tsv")"|\
            "$(realpath -m "${OUTDIR}/run_status.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/software_versions.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/database_versions.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/errors.tsv")"|\
            "$(realpath -m "${LOGS_DIR}/last_error.json")"|\
            "$(realpath -m "${LOGS_DIR}/checkpoint_state")"|\
            "$(realpath -m "${LOGS_DIR}/checkpoint_state")"/*|\
            "$(realpath -m "${LOGS_DIR}/last_valid_checkpoint.state")"|\
            "$(realpath -m "${LOGS_DIR}")"/.error_recorded.*)
                log_error "Refusing to clean protected persistent error or inventory path: ${path}"
                return 1 ;;
        esac
        if [[ "${id}" == 05a_vamb ]]; then
            local external_input external_input_abs output_abs
            for external_input in "${ASSEMBLY_IN}" "${JGI_DEPTH_FILE_ACTIVE}"; do
                [[ -n "${external_input}" && -e "${external_input}" ]] || continue
                external_input_abs="$(realpath -m "${external_input}")"
                output_abs="$(realpath -m "${path}")"
                if [[ "${external_input_abs}" == "${output_abs}" || \
                      "${external_input_abs}" == "${output_abs}/"* ]]; then
                    log_error "Refusing to clean Vamb outputs because an external input is inside a controlled output path: ${external_input}"
                    return 1
                fi
            done
        fi
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Plan] Would clean controlled output: ${path}"
        elif [[ -e "${path}" || -L "${path}" ]]; then
            rm -rf -- "${path}"
        fi
    done
}

calculate_stage_needs() {
    NEEDS_READS=false
    NEEDS_ASSEMBLY=false
    NEEDS_MAGS=false

    if [[ "${RUN_FASTP}" == true || "${RUN_KRAKEN2}" == true || \
          "${RUN_MEGAHIT}" == true || \
          ( "${RUN_MAPPING}" == true && "${COVERAGE_SOURCE:-internal}" == internal ) || \
          "${RUN_REASSEMBLY}" == true || "${RUN_COVERM}" == true ]]; then
        NEEDS_READS=true
    fi
    if [[ "${RUN_MEGAHIT}" == true || "${RUN_MAPPING}" == true || \
          "${RUN_BINNING}" == true || "${RUN_REFINEMENT_EFFECTIVE}" == true || \
          "${RUN_REASSEMBLY}" == true || "${RUN_CHECKM2}" == true || \
          "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || \
          "${RUN_COVERM}" == true || "${RUN_SUMMARY}" == true ]]; then
        NEEDS_ASSEMBLY=true
    fi
    if [[ "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || \
          "${RUN_COVERM}" == true || "${RUN_SUMMARY}" == true ]]; then
        NEEDS_MAGS=true
    fi
}

binner_selected() {
    local requested="$1" selected
    for selected in "${SELECTED_BINNER_ARRAY[@]}"; do
        [[ "${selected}" == "${requested}" ]] && return 0
    done
    return 1
}

get_binner_bins_dir() {
    case "$1" in
        vamb) printf '%s\n' "${VAMB_FINAL_BINS_DIR}" ;;
        metabinner) printf '%s\n' "${METABINNER_FINAL_BINS_DIR}" ;;
        comebin) printf '%s\n' "${COMEBIN_FINAL_BINS_DIR}" ;;
        *) log_error "Unknown binner: $1"; return 1 ;;
    esac
}

prepare_external_bins_cache() {
    local source_dir="$1" outdir_real="" external_real="" fingerprint=""
    local -a external_bins=()
    local entry="" name=""

    if [[ -z "${source_dir}" ]]; then
        log_error "--bins-dir requires a directory value."
        return 1
    fi
    if [[ ! -d "${source_dir}" ]]; then
        log_error "External bins directory not found: ${source_dir}"
        return 1
    fi

    outdir_real="$(realpath -m "${OUTDIR}")"
    external_real="$(realpath -- "${source_dir}")"

    if [[ "${external_real}" == "/" ]]; then
        log_error "External bins directory cannot be the filesystem root."
        return 1
    fi
    if [[ "${external_real}" == "${outdir_real}" || "${external_real}" == "${outdir_real}/"* || "${outdir_real}" == "${external_real}" || "${outdir_real}" == "${external_real}/"* ]]; then
        log_error "External bins directory cannot be OUTDIR or be inside/outside of each other in a way that cycles inputs and outputs: ${source_dir}"
        return 1
    fi

    mkdir -p "${EXTERNAL_BINS_CURRENT}"
    rm -rf -- "${EXTERNAL_FINAL_BINS_DIR}"
    mkdir -p "${EXTERNAL_FINAL_BINS_DIR}"

    mapfile -d '' -t external_bins < <(find "${external_real}" -mindepth 1 -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -print0 | sort -z)
    if (( ${#external_bins[@]} == 0 )); then
        log_error "External bins directory contains no final bin FASTA files: ${source_dir}"
        return 1
    fi

    for entry in "${external_bins[@]}"; do
        name="$(basename "${entry}")"
        name="${name%.*}.fa"
        cp -a -- "${entry}" "${EXTERNAL_FINAL_BINS_DIR}/${name}"
    done

    validate_bins_dir "${EXTERNAL_FINAL_BINS_DIR}" "external final bins" || return 1

    fingerprint="$(find "${EXTERNAL_FINAL_BINS_DIR}" -mindepth 1 -maxdepth 1 -type f -printf '%f\t%s\n' | sort | sha256sum | awk '{print $1}')"
    printf '%s\t%s\t%s\n' "${external_real}" "${fingerprint}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${EXTERNAL_BINS_MANIFEST}"
    BINS_FOR_QC="${EXTERNAL_FINAL_BINS_DIR}"
    BINS_FOR_REASSEMBLY="${EXTERNAL_FINAL_BINS_DIR}"
    BINS_SOURCE_RESOLVED="external"
    return 0
}

resolve_binning_configuration() {
    local errors=0 method token canonical=""
    local input="${SELECTED_BINNERS}"
    local -a parsed=()
    local -A selected=()

    case "${REFINEMENT_MODE}" in
        auto|always|never) ;;
        *) log_error "Invalid REFINEMENT_MODE '${REFINEMENT_MODE}' (expected auto, always, or never)."; errors=$((errors + 1)) ;;
    esac
    case "${BINS_SOURCE}" in
        auto|vamb|metabinner|comebin|refined|external) ;;
        *) log_error "Invalid BINS_SOURCE '${BINS_SOURCE}' (expected auto, vamb, metabinner, comebin, refined, or external)."; errors=$((errors + 1)) ;;
    esac

    if [[ -n "${BINS_DIR_INPUT}" ]]; then
        if [[ "${RUN_BINNING}" == true ]]; then
            log_error "External bins supplied with --bins-dir cannot be combined with internal binning in this version. Use --no-binning, or remove --bins-dir."
            errors=$((errors + 1))
        fi
        if [[ "${BINS_SOURCE}" != auto && "${BINS_SOURCE}" != external ]]; then
            log_error "External bins supplied with --bins-dir require --bins-source auto or external."
            errors=$((errors + 1))
        fi
        if (( errors == 0 )); then
            BINS_SOURCE="external"
            prepare_external_bins_cache "${BINS_DIR_INPUT}" || return 1
        fi
        BINNING_RESOLUTION_ERRORS="${errors}"
        return $(( errors == 0 ))
    fi

    SELECTED_BINNER_ARRAY=()
    if [[ -n "${input}" ]]; then
        if [[ "${input}" =~ ,[[:space:]]*$ || "${input}" =~ ^[[:space:]]*, ]]; then
            log_error "--binners contains an empty list item."
            errors=$((errors + 1))
        fi
        IFS=',' read -r -a parsed <<< "${input}"
        for token in "${parsed[@]}"; do
            token="${token#"${token%%[![:space:]]*}"}"
            token="${token%"${token##*[![:space:]]}"}"
            case "${token}" in
                vamb|metabinner|comebin) selected["${token}"]=1 ;;
                *) log_error "Unknown or empty binner in --binners: '${token}'."; errors=$((errors + 1)) ;;
            esac
        done
    elif [[ "${RUN_BINNING}" == true ]]; then
        log_error "--binners cannot be empty when RUN_BINNING=true."
        errors=$((errors + 1))
    fi
    for method in vamb metabinner comebin; do
        if [[ -n "${selected[$method]:-}" ]]; then
            SELECTED_BINNER_ARRAY+=("${method}")
            [[ -z "${canonical}" ]] || canonical+=","
            canonical+="${method}"
        fi
    done
    SELECTED_BINNERS="${canonical}"
    RUN_REFINEMENT_EFFECTIVE=false
    BINS_SOURCE_RESOLVED=""
    BINS_FOR_QC=""
    BINS_FOR_REASSEMBLY=""

    if [[ "${BINS_SOURCE}" == external && -z "${BINS_DIR_INPUT}" ]]; then
        log_error "--bins-source external requires --bins-dir DIR."
        errors=$((errors + 1))
    fi
    if [[ "${RUN_BINNING}" != true ]]; then
        if [[ "${BINS_SOURCE}" != auto ]]; then
            log_error "A non-auto --bins-source requires RUN_BINNING=true; external bins are not supported."
            errors=$((errors + 1))
        fi
    elif (( ${#SELECTED_BINNER_ARRAY[@]} == 1 )); then
        if [[ "${REFINEMENT_MODE}" == always ]]; then
            log_error "Consensus refinement needs at least two bin collections; use --refinement-mode auto or never."
            errors=$((errors + 1))
        fi
        method="${SELECTED_BINNER_ARRAY[0]}"
        case "${BINS_SOURCE}" in
            auto) BINS_SOURCE_RESOLVED="${method}" ;;
            refined)
                log_error "--bins-source refined requires at least two selected binners and active refinement."
                errors=$((errors + 1)) ;;
            "${method}") BINS_SOURCE_RESOLVED="${method}" ;;
            *)
                log_error "BINS_SOURCE '${BINS_SOURCE}' is not the selected binner '${method}'."
                errors=$((errors + 1)) ;;
        esac
    elif (( ${#SELECTED_BINNER_ARRAY[@]} >= 2 )); then
        case "${REFINEMENT_MODE}" in
            auto)
                case "${BINS_SOURCE}" in
                    auto|refined)
                        RUN_REFINEMENT_EFFECTIVE=true
                        BINS_SOURCE_RESOLVED=refined ;;
                    vamb|metabinner|comebin)
                        if binner_selected "${BINS_SOURCE}"; then
                            BINS_SOURCE_RESOLVED="${BINS_SOURCE}"
                        else
                            log_error "BINS_SOURCE '${BINS_SOURCE}' is not included in --binners (${SELECTED_BINNERS})."
                            errors=$((errors + 1))
                        fi ;;
                esac ;;
            always)
                case "${BINS_SOURCE}" in
                    auto|refined)
                        RUN_REFINEMENT_EFFECTIVE=true
                        BINS_SOURCE_RESOLVED=refined ;;
                    *)
                        if binner_selected "${BINS_SOURCE}"; then
                            log_error "Refinement was requested explicitly, but BINS_SOURCE selects direct bins; use --refinement-mode auto or never."
                        else
                            log_error "BINS_SOURCE '${BINS_SOURCE}' is not included in --binners (${SELECTED_BINNERS})."
                        fi
                        errors=$((errors + 1)) ;;
                esac ;;
            never)
                case "${BINS_SOURCE}" in
                    auto)
                        log_error "Multiple binners make BINS_SOURCE=auto ambiguous when refinement is disabled; select a binner or use --refinement-mode auto." 
                        errors=$((errors + 1)) ;;
                    refined)
                        log_error "BINS_SOURCE=refined requires refinement; use --refinement-mode auto or select a direct binner."
                        errors=$((errors + 1)) ;;
                    *)
                        if binner_selected "${BINS_SOURCE}"; then
                            BINS_SOURCE_RESOLVED="${BINS_SOURCE}"
                        else
                            log_error "BINS_SOURCE '${BINS_SOURCE}' is not included in --binners (${SELECTED_BINNERS})."
                            errors=$((errors + 1))
                        fi ;;
                esac ;;
        esac
    elif [[ "${RUN_BINNING}" == true ]]; then
        log_error "At least one binner must be selected when RUN_BINNING=true."
        errors=$((errors + 1))
    fi

    case "${BINS_SOURCE_RESOLVED}" in
        vamb) BINS_FOR_QC="${VAMB_FINAL_BINS_DIR}" ;;
        metabinner) BINS_FOR_QC="${METABINNER_FINAL_BINS_DIR}" ;;
        comebin) BINS_FOR_QC="${COMEBIN_FINAL_BINS_DIR}" ;;
        refined) BINS_FOR_QC="${REFINED_BINS_DIR}" ;;
    esac
    BINS_FOR_REASSEMBLY="${BINS_FOR_QC}"
    if [[ "${RUN_REASSEMBLY}" == true ]]; then
        BINS_FOR_QC="${REASSEMBLY_DIR}/reassembled_bins"
    fi
    BINNING_RESOLUTION_ERRORS="${errors}"
    (( errors == 0 ))
}

resolve_active_reads() {
    local source_dir="" source_kind="" current_manifest="" current_metadata="" r1 sample r2
    local -A seen_samples=()
    ACTIVE_R1_FILES=()
    ACTIVE_R2_FILES=()
    ACTIVE_SAMPLE_IDS=()
    INPUT_SAMPLE_IDS=()
    INPUT_R1_FILES=()
    INPUT_R2_FILES=()

    if [[ -n "${RAW_READS_DIR}" && -n "${CLEAN_READS_DIR}" ]]; then
        log_error "--reads-dir and --clean-reads-dir are mutually exclusive."
        return 1
    fi
    if [[ -n "${CLEAN_READS_DIR}" && "${RUN_FASTP}" == true ]]; then
        log_error "--clean-reads-dir cannot be combined with fastp; use --no-fastp."
        return 1
    fi

    if [[ -n "${CLEAN_READS_DIR}" ]]; then
        source_dir="$(realpath -m "${CLEAN_READS_DIR}")"
        source_kind="external-clean"
        ACTIVE_READS_DIR="${source_dir}"
        ACTIVE_READS_SOURCE="external-clean"
        current_manifest="${EXTERNAL_CLEAN_SAMPLE_MANIFEST}"
        current_metadata="${EXTERNAL_CLEAN_SAMPLE_METADATA}"
    elif [[ -n "${RAW_READS_DIR}" ]]; then
        source_dir="$(realpath -m "${RAW_READS_DIR}")"
        source_kind="raw"
        RAW_READS_DIR="${source_dir}"
        READS_DIR="${source_dir}"
        if [[ "${RUN_FASTP}" == true ]]; then
            ACTIVE_READS_SOURCE="fastp-clean"
            ACTIVE_READS_DIR="$(realpath -m "${QC_DIR}/trimmed")"
        else
            ACTIVE_READS_SOURCE="raw"
            ACTIVE_READS_DIR="${source_dir}"
            current_manifest="${RAW_SAMPLE_MANIFEST}"
            current_metadata="${RAW_SAMPLE_METADATA}"
        fi
    elif [[ "${NEEDS_READS}" == true ]]; then
        log_error "An active stage requires reads; provide --reads-dir or --clean-reads-dir."
        return 1
    else
        ACTIVE_READS_DIR=""
        ACTIVE_READS_SOURCE=""
        return 0
    fi

    if [[ ! -d "${source_dir}" ]]; then
        log_error "Read directory does not exist: ${source_dir}"
        return 1
    fi

    if [[ "${RUN_FASTP}" == true && "${ACTIVE_READS_SOURCE}" == fastp-clean ]]; then
        if [[ -s "${RAW_SAMPLE_MANIFEST}" && -f "${RAW_SAMPLE_MANIFEST}" ]]; then
            load_sample_manifest "${RAW_SAMPLE_MANIFEST}" "raw" INPUT_SAMPLE_IDS INPUT_R1_FILES INPUT_R2_FILES || return 1
        else
            local candidate_tmp
            candidate_tmp="$(mktemp "${TMPDIR:-/tmp}/.raw_sample_manifest.XXXXXX")"
            build_sample_manifest_candidate "${source_dir}" "raw" "${candidate_tmp}" || return 1
            if [[ "${DRY_RUN}" == true ]]; then
                validate_sample_manifest "${candidate_tmp}" "raw" || {
                    rm -f "${candidate_tmp}"
                    return 1
                }
                rm -f "${candidate_tmp}"
                log_info "[DRY-RUN] Candidate raw sample manifest validated for ${RAW_READS_DIR}."
                return 0
            fi
            publish_sample_manifest "${candidate_tmp}" "${RAW_SAMPLE_MANIFEST}" "${RAW_SAMPLE_METADATA}" "raw" || return 1
            rm -f "${candidate_tmp}"
            load_sample_manifest "${RAW_SAMPLE_MANIFEST}" "raw" INPUT_SAMPLE_IDS INPUT_R1_FILES INPUT_R2_FILES || return 1
        fi
        ACTIVE_READS_DIR="$(realpath -m "${QC_DIR}/trimmed")"
        return 0
    fi

    if [[ -s "${current_manifest:-}" ]]; then
        load_sample_manifest "${current_manifest}" "${source_kind}" ACTIVE_SAMPLE_IDS ACTIVE_R1_FILES ACTIVE_R2_FILES || return 1
    else
        local candidate_tmp
        candidate_tmp="$(mktemp "${TMPDIR:-/tmp}/.${source_kind}_sample_manifest.XXXXXX")"
        build_sample_manifest_candidate "${source_dir}" "${source_kind}" "${candidate_tmp}" || return 1
        if [[ "${DRY_RUN}" == true ]]; then
            validate_sample_manifest "${candidate_tmp}" "${source_kind}" || {
                rm -f "${candidate_tmp}"
                return 1
            }
            rm -f "${candidate_tmp}"
            log_info "[DRY-RUN] Candidate ${source_kind} sample manifest validated for ${source_dir}."
            return 0
        fi
        publish_sample_manifest "${candidate_tmp}" "${current_manifest}" "${current_metadata}" "${source_kind}" || return 1
        rm -f "${candidate_tmp}"
        load_sample_manifest "${current_manifest}" "${source_kind}" ACTIVE_SAMPLE_IDS ACTIVE_R1_FILES ACTIVE_R2_FILES || return 1
    fi
    if (( ${#ACTIVE_R1_FILES[@]} != ${#ACTIVE_SAMPLE_IDS[@]} || ${#ACTIVE_R2_FILES[@]} != ${#ACTIVE_SAMPLE_IDS[@]} )); then
        log_error "Manifest-loaded read arrays cannot be reconciled with sample IDs."
        return 1
    fi
    return 0
}

# ── Resource Detection ────────────────────────────────────────────────────────
detect_resources() {
    log_info "Detecting hardware resources..."

    # ── CPU detection ────────────────────────────────────────────────────────
    local total_cores
    if command -v nproc &>/dev/null; then
        total_cores=$(nproc)
    elif command -v getconf &>/dev/null; then
        total_cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
    elif [[ -r /proc/cpuinfo ]]; then
        total_cores=$(grep -c '^processor' /proc/cpuinfo)
    else
        total_cores=1
        log_warn "[Resources] Cannot detect CPU count; defaulting to 1 core."
    fi
    _PHYS_CORES="${total_cores}"

    if [[ "${THREADS}" == "auto" ]]; then
        if (( total_cores > 4 )); then
            THREADS=$(( total_cores - CPU_RESERVE ))
        else
            THREADS=$(( total_cores > 0 ? total_cores : 1 ))
        fi
        log_info "[Resources] CPU auto-detected: ${THREADS} of ${total_cores} logical cores assigned (reserve=${CPU_RESERVE})."
    else
        # Manual value — validate against hardware
        if (( THREADS > total_cores )); then
            log_warn "[Resources] THREADS=${THREADS} exceeds detected cores (${total_cores}). Consider lowering to avoid contention."
        fi
        log_info "[Resources] CPU (manual): ${THREADS} of ${total_cores} logical cores."
    fi

    # ── RAM detection ────────────────────────────────────────────────────────
    local mem_kb=0
    if [[ -r /proc/meminfo ]]; then
        # Prefer MemAvailable (free + reclaimable); fall back to MemTotal
        mem_kb=$(awk '/^MemAvailable:/{print $2; found=1} END{if(!found) exit 1}' \
                      /proc/meminfo 2>/dev/null) || \
        mem_kb=$(awk '/^MemTotal:/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
    else
        log_warn "[Resources] /proc/meminfo not readable; defaulting to 4 GB."
        mem_kb=4194304  # 4 GB in kB
    fi

    local mem_total_gb=$(( mem_kb / 1048576 ))   # kB → GB (integer)
    [[ "${mem_total_gb}" -lt 1 ]] && mem_total_gb=1
    _PHYS_MEM_GB="${mem_total_gb}"

    if [[ "${MEMORY}" == "auto" ]]; then
        MEMORY=$(( mem_total_gb * MEM_FRACTION / 100 ))
        [[ "${MEMORY}" -lt 4 ]] && MEMORY=4   # hard floor: 4 GB
        log_info "[Resources] RAM auto-detected: ${MEMORY} GB assigned from ${mem_total_gb} GB available (${MEM_FRACTION}% cap)."
    else
        # Manual value — validate
        if (( MEMORY > mem_total_gb )); then
            log_warn "[Resources] MEMORY=${MEMORY} GB exceeds detected RAM (${mem_total_gb} GB). OOM risk!"
        fi
        log_info "[Resources] RAM (manual): ${MEMORY} GB of ${mem_total_gb} GB available."
    fi
}

# ── Pre-Flight Validation ─────────────────────────────────────────────────────
preflight_check() {
    log_info "Running pre-flight validation..."
    local errors="${BINNING_RESOLUTION_ERRORS}"
    local need_tools=false
    if ! validate_percentage_threshold MIN_COMPLETENESS "${MIN_COMPLETENESS}"; then
        errors=$((errors + 1))
    fi
    if ! validate_percentage_threshold MAX_CONTAMINATION "${MAX_CONTAMINATION}"; then
        errors=$((errors + 1))
    fi
    if ! resolve_coverage_inputs; then errors=$((errors + 1)); fi
    calculate_stage_needs
    if ! validate_stage_dependencies; then
        errors=$((errors + 1))
    fi
    if ! resolve_active_assembly; then errors=$((errors + 1)); fi
    if ! resolve_active_reads; then
        (( errors += 1 ))
    fi
    if [[ "${COVERAGE_SOURCE}" == external ]]; then
        if [[ -n "${JGI_DEPTH_FILE_ACTIVE}" ]]; then
            DEPTH_FILE="${JGI_DEPTH_FILE_ACTIVE}"
            if [[ "${RUN_BINNING}" == true ]]; then
                local binner
                for binner in "${SELECTED_BINNER_ARRAY[@]}"; do
                    if binner_requires_input "${binner}" jgi-depth; then
                        validate_external_jgi_depth || errors=$((errors + 1))
                        break
                    fi
                done
            fi
        fi
        if ! validate_external_bams; then errors=$((errors + 1)); fi
    fi
    if [[ "${COVERAGE_SOURCE}" == internal || "${COVERAGE_SOURCE}" == external || \
          "${COVERAGE_SOURCE}" == none ]]; then
        DEPTH_FILE="${JGI_DEPTH_FILE_ACTIVE}"
        if [[ "${RUN_BINNING}" == true ]]; then
            local binner
            for binner in "${SELECTED_BINNER_ARRAY[@]}"; do
                if ! validate_binner_inputs "${binner}"; then errors=$((errors + 1)); fi
            done
        fi
    elif [[ "${RUN_BINNING}" == true ]]; then
        local binner
        for binner in "${SELECTED_BINNER_ARRAY[@]}"; do
            if ! validate_binner_inputs "${binner}"; then errors=$((errors + 1)); fi
        done
    fi
    if ! validate_vamb_preflight; then
        errors=$((errors + 1))
    fi
    if [[ "${RUN_BINNING}" == true ]] && binner_selected comebin; then
        if [[ ! "${COMEBIN_SEED}" =~ ^[0-9]+$ ]]; then
            log_error "COMEBIN_SEED must be a non-negative integer."
            errors=$((errors + 1))
        fi
        if ! validate_comebin_preflight; then
            errors=$((errors + 1))
        fi
    elif [[ "${COMEBIN_DEVICE_CLI_SET}" == true ]]; then
        if ! validate_comebin_device; then
            errors=$((errors + 1))
        fi
    fi
    [[ "${RUN_MAG_FILTERING}" == true ]] && HQ_BINS_DIR="${MAGS_DIR}"
    if ! collect_execution_plan; then
        errors=$((errors + 1))
    fi

    # ── Temporary workspace: resolve base and freeze per-run paths ────────────
    # Must run after collect_execution_plan so we know the reassembly action.
    if ! resolve_pipeline_tmpdir_base; then
        (( errors += 1 ))
    else
        validate_pipeline_tmpdir_base || (( errors += 1 )) || true
        # Freeze the per-run paths (after plan is frozen)
        PIPELINE_RUN_TMPDIR="${PIPELINE_TMP_PARENT}/${RUN_ID}"
        REASSEMBLY_TMPDIR="${PIPELINE_RUN_TMPDIR}/reassembly"
        PIPELINE_TMP_MARKER="${PIPELINE_RUN_TMPDIR}/.pipeline_tmp_owner"
        if [[ "${DRY_RUN}" == true ]]; then
            show_tmpdir_dry_run_info
        elif [[ "${RUN_REASSEMBLY}" == true ]] && stage_action_requires_execution 07_reassembly; then
            if ! check_tmpdir_space ACTIVE_R1_FILES ACTIVE_R2_FILES; then
                (( errors += 1 ))
            fi
        fi
    fi

    # External assemblies are opt-in and never disable MEGAHIT implicitly.
    if [[ "${RUN_MEGAHIT}" == true && -n "${ASSEMBLY_IN}" ]]; then
        log_error "--assembly-in cannot be combined with MEGAHIT; pass --no-megahit explicitly."
        (( errors += 1 ))
    fi
    if [[ "${RUN_MEGAHIT}" == false && "${NEEDS_ASSEMBLY}" == true && -z "${ASSEMBLY_IN}" ]]; then
        log_error "An active stage requires an assembly, but MEGAHIT is disabled and ASSEMBLY_IN is not set."
        (( errors += 1 ))
    fi
    if [[ "${RUN_MEGAHIT}" == false && -n "${ASSEMBLY_IN}" ]]; then
        if [[ ! -f "${ASSEMBLY_IN}" || ! -s "${ASSEMBLY_IN}" ]]; then
            log_error "External assembly file does not exist or is empty: ${ASSEMBLY_IN}"
            (( errors += 1 ))
        else
            ASSEMBLY_IN="$(realpath "${ASSEMBLY_IN}")"
            if ! awk '/^>/ { header=1; next } header && /[^[:space:]]/ { sequence=1 } END { exit !(header && sequence) }' \
                    "${ASSEMBLY_IN}"; then
                log_error "External assembly is not a valid FASTA with at least one header and sequence: ${ASSEMBLY_IN}"
                (( errors += 1 ))
            fi
        fi
    fi

    if [[ "${NEEDS_READS}" == true && -z "${ACTIVE_READS_DIR}" ]]; then
        log_error "An active stage requires a resolved source of reads."
        (( errors += 1 ))
    fi

    # Kraken2 DB
    if [[ "${RUN_KRAKEN2}" == true ]] && stage_action_requires_execution 03_kraken2; then
        if [[ ! -d "${KRAKEN2_DB}" ]]; then
            log_error "KRAKEN2_DB not set or not found: ${KRAKEN2_DB}"
            (( errors += 1 ))
        fi
    fi

    # CheckM2 DB
    if [[ "${RUN_CHECKM2}" == true ]] && stage_in_scope checkm2 && \
       stage_action_requires_execution 08_checkm2 && [[ ! -f "${CHECKM2_DB}" ]]; then
        log_error "CHECKM2_DB not set or not found: ${CHECKM2_DB}"
        (( errors += 1 ))
    fi

    # GTDB-Tk data
    if [[ "${RUN_GTDBTK}" == true ]] && stage_in_scope gtdbtk && \
       stage_action_requires_execution 10_gtdbtk; then
        if [[ ! -d "${GTDBTK_DATA_PATH}" ]]; then
            log_error "GTDBTK_DATA_PATH not set or not found: ${GTDBTK_DATA_PATH}"
            (( errors += 1 ))
        fi
    fi

    # MetaBinner scripts directory — validation differs by engine
    if [[ "${RUN_BINNING}" == true ]] && binner_selected metabinner && \
       stage_in_scope binning && stage_action_requires_execution 05b_metabinner; then
        if [[ "${EXEC_MODE}" == "mamba" ]]; then
            if [[ -z "${METABINNER_SCRIPTS_DIR}" || ! -d "${METABINNER_SCRIPTS_DIR}" ]]; then
                log_error "METABINNER_SCRIPTS_DIR not set or not found on host: ${METABINNER_SCRIPTS_DIR}"
                (( errors += 1 ))
            fi
        else
            if [[ -z "${METABINNER_SCRIPTS_DIR}" ]]; then
                log_error "METABINNER_SCRIPTS_DIR is empty."
                log_error "For EXEC_MODE=apptainer use an in-container path (default: /opt/metabinner/scripts)."
                (( errors += 1 ))
            elif [[ -f "${SIF_METABINNER}" ]] && [[ "${DRY_RUN}" == false ]]; then
                local _rt="apptainer"
                command -v apptainer &>/dev/null || _rt="singularity"
                log_info "[Preflight] Verifying METABINNER_SCRIPTS_DIR inside container image..."
                if ! "${_rt}" exec --cleanenv "${SIF_METABINNER}" \
                        test -d "${METABINNER_SCRIPTS_DIR}" 2>/dev/null; then
                    log_error "METABINNER_SCRIPTS_DIR '${METABINNER_SCRIPTS_DIR}' not found inside ${SIF_METABINNER}."
                    log_error "Check the container's internal layout or adjust METABINNER_SCRIPTS_DIR."
                    (( errors += 1 ))
                else
                    log_info "[Preflight] METABINNER_SCRIPTS_DIR OK inside container: ${METABINNER_SCRIPTS_DIR}"
                fi
            elif [[ ! -f "${SIF_METABINNER}" ]] && [[ "${AUTO_PULL_SIF}" == true ]]; then
                log_info "[Preflight] ${SIF_METABINNER} not yet downloaded; METABINNER_SCRIPTS_DIR will be validated after pull."
            fi
        fi
    fi

    # Required direct binaries are host dependencies only in mamba mode.
    if { [[ "${RUN_FASTP}" == true ]] && stage_action_requires_execution 01_fastp; } || \
       { [[ "${RUN_KRAKEN2}" == true ]] && stage_action_requires_execution 03_kraken2; } || \
       { [[ "${RUN_MEGAHIT}" == true ]] && stage_action_requires_execution 02_megahit; }; then
        need_tools=true
    fi
    local stage
    for stage in mapping binning refinement reassembly checkm2 gtdbtk coverm summary; do
        if stage_in_scope "${stage}"; then
            case "${stage}" in
                mapping) if [[ "${RUN_MAPPING}" == true ]] && stage_action_requires_execution 04_mapping; then need_tools=true; fi ;;
                binning)
                    if [[ "${RUN_BINNING}" == true ]] && \
                       { stage_action_requires_execution 05a_vamb || stage_action_requires_execution 05b_metabinner || \
                         stage_action_requires_execution 05c_comebin; }; then need_tools=true; fi ;;
                refinement) if [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] && stage_action_requires_execution 06_refinement; then need_tools=true; fi ;;
                reassembly) if [[ "${RUN_REASSEMBLY}" == true ]] && stage_action_requires_execution 07_reassembly; then need_tools=true; fi ;;
                checkm2) if [[ "${RUN_CHECKM2}" == true ]] && stage_action_requires_execution 08_checkm2; then need_tools=true; fi ;;
                gtdbtk) if [[ "${RUN_GTDBTK}" == true ]] && stage_action_requires_execution 10_gtdbtk; then need_tools=true; fi ;;
                coverm) if [[ "${RUN_COVERM}" == true ]] && stage_action_requires_execution 11_coverm; then need_tools=true; fi ;;
                summary|taxonomy-summary|mag-summary|run-overview) ;;
            esac
        fi
    done

    if [[ "${EXEC_MODE}" == "mamba" ]] && [[ "${need_tools}" == true ]]; then
        local required_bins=()
        if [[ "${RUN_MAPPING}" == true ]] && stage_in_scope mapping && stage_action_requires_execution 04_mapping; then
            required_bins+=(bowtie2 bowtie2-build samtools jgi_summarize_bam_contig_depths)
        fi
        if [[ "${RUN_COVERM}" == true ]] && stage_in_scope coverm && stage_action_requires_execution 11_coverm; then
            required_bins+=(coverm)
        fi
        if [[ "${RUN_FASTP}" == true ]] && stage_action_requires_execution 01_fastp; then required_bins+=(fastp); fi
        if [[ "${RUN_MEGAHIT}" == true ]] && stage_action_requires_execution 02_megahit; then required_bins+=(megahit); fi
        if [[ "${RUN_KRAKEN2}" == true ]] && stage_action_requires_execution 03_kraken2; then required_bins+=(kraken2 bracken); fi
        for tool in "${required_bins[@]}"; do
            if ! command -v "${tool}" &>/dev/null; then
                log_error "Required binary not found in PATH: ${tool}"
                (( errors += 1 ))
            fi
        done
    fi

    # Container / environment runtime
    if [[ "${EXEC_MODE}" == "apptainer" ]] && [[ "${need_tools}" == true ]]; then
        if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
            log_error "EXEC_MODE=apptainer but neither 'apptainer' nor 'singularity' found in PATH."
            (( errors += 1 ))
        else
            local _rt="apptainer"
            command -v apptainer &>/dev/null || _rt="singularity"
            log_info "Container runtime detected: ${_rt}"
        fi
        # Warn about SIF_DIR accessibility
        if [[ ! -d "${SIF_DIR}" ]]; then
            log_warn "SIF_DIR does not exist yet (${SIF_DIR}); it will be created on first pull."
        fi
    elif [[ "${EXEC_MODE}" == "mamba" ]] && [[ "${need_tools}" == true ]]; then
        if ! command -v mamba &>/dev/null; then
            log_error "'mamba' not found in PATH."
            (( errors += 1 ))
        else
            local -a required_envs=()
            local env_name
            if [[ "${RUN_BINNING}" == true ]] && stage_in_scope binning; then
                binner_selected vamb && stage_action_requires_execution 05a_vamb && required_envs+=("${ENV_VAMB}")
                binner_selected metabinner && stage_action_requires_execution 05b_metabinner && required_envs+=("${ENV_METABINNER}")
                binner_selected comebin && stage_action_requires_execution 05c_comebin && required_envs+=("${ENV_COMEBIN}")
            fi
            if [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] && stage_in_scope refinement; then
                stage_action_requires_execution 06_refinement && required_envs+=("${ENV_METAWRAP}")
            fi
            if [[ "${RUN_REASSEMBLY}" == true ]] && stage_in_scope reassembly; then
                stage_action_requires_execution 07_reassembly && required_envs+=("${ENV_METAWRAP}")
            fi
            if [[ "${RUN_CHECKM2}" == true ]] && stage_in_scope checkm2; then
                stage_action_requires_execution 08_checkm2 && required_envs+=("${ENV_CHECKM2}")
            fi
            if [[ "${RUN_GTDBTK}" == true ]] && stage_in_scope gtdbtk; then
                stage_action_requires_execution 10_gtdbtk && required_envs+=("${ENV_GTDBTK}")
            fi
            local env_list=""
            if (( ${#required_envs[@]} > 0 )); then
                if ! env_list="$(mamba env list 2>&1)"; then
                    log_error "Could not list Mamba environments: ${env_list}"
                    (( errors += 1 ))
                else
                    for env_name in "${required_envs[@]}"; do
                        if ! grep -Eq "^[[:space:]]*${env_name}[[:space:]]" <<< "${env_list}"; then
                            log_error "Required Mamba environment not found: ${env_name}"
                            (( errors += 1 ))
                        fi
                    done
                fi
            fi
        fi
    fi

    # Detect reads and verify each R1 has its matching R2.
    log_info "Requested profile: ${REQUESTED_PROFILE}; resulting profile: ${RESULTING_PROFILE}."
    log_info "Reads source: ${ACTIVE_READS_SOURCE:-<none>}; active reads directory: ${ACTIVE_READS_DIR:-<none>}."
    log_info "Coverage source: ${COVERAGE_SOURCE}; BAM directory: ${BAM_DIR_ACTIVE:-<none>}; JGI depth: ${JGI_DEPTH_FILE_ACTIVE:-<none>}."
    log_info "Selected binners: ${SELECTED_BINNERS:-<none>}."
    log_info "Refinement mode: ${REFINEMENT_MODE}; effective: ${RUN_REFINEMENT_EFFECTIVE}."
    log_info "Bins source requested: ${BINS_SOURCE}; resolved: ${BINS_SOURCE_RESOLVED:-<none>}."
    log_info "Bins directory for QC: ${BINS_FOR_QC:-<none>}."
    log_info "Minimum completeness  : ${MIN_COMPLETENESS}"
    log_info "Maximum contamination : ${MAX_CONTAMINATION}"
    log_info "Thresholds used       : ${MAG_FILTERING_APPLICABLE}"
    if [[ "${MAG_FILTERING_APPLICABLE}" == true ]]; then
        log_info "Retention rule        : completeness >= ${MIN_COMPLETENESS} AND contamination <= ${MAX_CONTAMINATION}"
    fi
    log_info "Derived requirements: NEEDS_READS=${NEEDS_READS}, NEEDS_ASSEMBLY=${NEEDS_ASSEMBLY}, NEEDS_MAGS=${NEEDS_MAGS}."
    if [[ "${RUN_KRAKEN2}" == true ]]; then
        if [[ "${NEEDS_ASSEMBLY}" == true ]]; then
            log_info "Active branches: taxonomy and assembly."
        else
            log_info "Active branches: taxonomy."
        fi
    elif [[ "${NEEDS_ASSEMBLY}" == true ]]; then
        log_info "Active branches: assembly."
    else
        log_info "Active branches: none."
    fi
    local -a stage_status=()
    local stage_name stage_flag
    for stage_name in fastp kraken2 megahit mapping binning refinement reassembly checkm2 mag-filtering gtdbtk coverm taxonomy-summary mag-summary run-overview; do
        case "${stage_name}" in
            fastp) stage_flag="${RUN_FASTP}" ;;
            kraken2) stage_flag="${RUN_KRAKEN2}" ;;
            megahit) stage_flag="${RUN_MEGAHIT}" ;;
            mapping) stage_flag="${RUN_MAPPING}" ;;
            binning) stage_flag="${RUN_BINNING}" ;;
            refinement) stage_flag="${RUN_REFINEMENT_EFFECTIVE}" ;;
            reassembly) stage_flag="${RUN_REASSEMBLY}" ;;
            checkm2) stage_flag="${RUN_CHECKM2}" ;;
            mag-filtering) stage_flag="${RUN_MAG_FILTERING}" ;;
            gtdbtk) stage_flag="${RUN_GTDBTK}" ;;
            coverm) stage_flag="${RUN_COVERM}" ;;
            taxonomy-summary) stage_flag=false; [[ "${RUN_SUMMARY}" == true && "${RUN_TAXONOMY_SUMMARY}" == true && "${RUN_KRAKEN2}" == true ]] && stage_flag=true ;;
            mag-summary) stage_flag=false; [[ "${RUN_SUMMARY}" == true && "${RUN_MAG_SUMMARY}" == true && ( "${RUN_CHECKM2}" == true || "${RUN_MAG_FILTERING}" == true || "${RUN_GTDBTK}" == true || "${RUN_COVERM}" == true ) ]] && stage_flag=true ;;
            run-overview) stage_flag=false; [[ "${RUN_SUMMARY}" == true && "${RUN_RUN_OVERVIEW}" == true ]] && stage_flag=true ;;
        esac
        if [[ "${stage_flag}" == true ]]; then
            stage_status+=("[RUN] ${stage_name}")
        else
            stage_status+=("[SKIP] ${stage_name}")
        fi
    done
    log_info "Final stage status: ${stage_status[*]}."
    if ! run_disk_preflight; then
        (( errors += 1 ))
    fi
    if [[ "${errors}" -gt 0 ]]; then
        log_error "Pre-flight failed with ${errors} error(s). Fix above and retry."
        exit 1
    fi
    log_info "Effective dependency graph:"
    local graph_parent graph_child graph_parent_name graph_child_name
    for graph_parent in "${STAGE_IDS[@]}"; do
        for graph_child in ${GRAPH_CHILDREN[$graph_parent]:-}; do
            graph_parent_name="$(get_stage_public_name "${graph_parent}")"
            graph_child_name="$(get_stage_public_name "${graph_child}")"
            log_info "  ${graph_parent_name} -> ${graph_child_name}"
        done
    done
    log_info "Planned stage actions:"
    for graph_parent in "${STAGE_IDS[@]}"; do
        log_info "  $(get_stage_public_name "${graph_parent}") [${STAGE_STATE[$graph_parent]}]: ${STAGE_ACTION[$graph_parent]} (${STAGE_REASON[$graph_parent]})"
    done
    # Detect GNU time for resource metrics
    if /usr/bin/time --version 2>&1 | grep -q "GNU Time"; then
        GNU_TIME_BIN="/usr/bin/time"
        GNU_TIME_AVAILABLE=true
        log_info "[Benchmark] GNU time detected: ${GNU_TIME_BIN}."
    else
        GNU_TIME_BIN=""
        GNU_TIME_AVAILABLE=false
        log_info "[Benchmark] GNU time not detected; CPU/RSS metrics will be recorded as NA."
    fi
    PIPELINE_START_MONOTONIC="$(python3 -c 'import time; print(f"{time.monotonic():.6f}")' 2>/dev/null || echo "0.0")"
    PIPELINE_INVOCATION_COMMAND="$(sanitize_command "${PIPELINE_INVOCATION_RAW[@]}")"

    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[DRY-RUN] Benchmark manifests planned:"
        log_info "  Run manifest   : would record ${LOGS_DIR}/run_manifest.tsv"
        log_info "  Stage status   : would record ${LOGS_DIR}/stage_status.tsv"
        log_info "  Compatibility  : would mirror to ${OUTDIR}/run_status.tsv"
    fi

    persist_resume_plan
    initialize_run_status
    log_info "Pre-flight OK."
}

###############################################################################
#  STEP FUNCTIONS                                                              #
###############################################################################

# ── Step 01: fastp QC ─────────────────────────────────────────────────────────
step_01_fastp() {
    if [[ "${RUN_FASTP}" == false ]]; then
        log_info "[SKIP] fastp stage disabled."
        return 0
    fi
    log_step "01" "Read quality control with fastp"

    if is_done "01_fastp"; then
        log_info "[Step 01] Checkpoint found — repopulating CLEAN arrays from existing trimmed files."
        if [[ -s "${FASTP_CLEAN_SAMPLE_MANIFEST}" ]]; then
            ACTIVE_READS_SOURCE="fastp-clean"
            load_sample_manifest "${FASTP_CLEAN_SAMPLE_MANIFEST}" "fastp-clean" ACTIVE_SAMPLE_IDS ACTIVE_R1_FILES ACTIVE_R2_FILES || return 1
        else
            collect_files_by_suffix "${QC_DIR}/trimmed" "${READ_SUFFIX_R1}" CLEAN_R1_FILES
            CLEAN_R2_FILES=()
            for r1 in "${CLEAN_R1_FILES[@]}"; do
                CLEAN_R2_FILES+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
            done
            ACTIVE_R1_FILES=("${CLEAN_R1_FILES[@]}")
            ACTIVE_R2_FILES=("${CLEAN_R2_FILES[@]}")
        fi
        log_info "[Step 01] Repopulated ${#ACTIVE_R1_FILES[@]} sample(s) from checkpoint."
        return 0
    fi

    mkdir -p "${QC_DIR}/trimmed" "${QC_DIR}/reports/html" "${QC_DIR}/reports/json"
    if [[ ${#INPUT_SAMPLE_IDS[@]} -eq 0 ]]; then
        load_sample_manifest "${RAW_SAMPLE_MANIFEST}" "raw" INPUT_SAMPLE_IDS INPUT_R1_FILES INPUT_R2_FILES || {
            discover_paired_reads "${RAW_READS_DIR}" INPUT_SAMPLE_IDS INPUT_R1_FILES INPUT_R2_FILES || exit 1
        }
    fi

    local fastp_threads=$(( THREADS > 16 ? 16 : THREADS ))
    CLEAN_R1_FILES=()
    CLEAN_R2_FILES=()
    ACTIVE_SAMPLE_IDS=()
    for i in "${!INPUT_SAMPLE_IDS[@]}"; do
        local sample="${INPUT_SAMPLE_IDS[$i]}"
        local r1="${INPUT_R1_FILES[$i]}"
        local r2="${INPUT_R2_FILES[$i]}"
        local out_r1="${QC_DIR}/trimmed/${sample}${READ_SUFFIX_R1}"
        local out_r2="${QC_DIR}/trimmed/${sample}${READ_SUFFIX_R2}"
        local html_report="${QC_DIR}/reports/html/${sample}.html"
        local json_report="${QC_DIR}/reports/json/${sample}.json"
        local step_log="${LOGS_DIR}/step_01_fastp_${sample}.log"

        log_info "[Step 01] Processing sample: ${sample}"
        run_cmd_with_tee "${step_log}" -- \
            run_tool "01" "${SIF_FASTP}" "${SIF_FASTP_URI}" fastp \
                -i "${r1}" \
                -I "${r2}" \
                -o "${out_r1}" \
                -O "${out_r2}" \
                -h "${html_report}" \
                -j "${json_report}" \
                -D \
                -5 20 \
                -3 20 \
                -r 20 \
                -M 20 \
                -e 20 \
                -n 50 \
                -w "${fastp_threads}"

        CLEAN_R1_FILES+=("${out_r1}")
        CLEAN_R2_FILES+=("${out_r2}")
        ACTIVE_SAMPLE_IDS+=("${sample}")
    done

    ACTIVE_R1_FILES=("${CLEAN_R1_FILES[@]}")
    ACTIVE_R2_FILES=("${CLEAN_R2_FILES[@]}")
    ACTIVE_READS_SOURCE="fastp-clean"
    ACTIVE_READS_DIR="$(realpath -m "${QC_DIR}/trimmed")"
    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[DRY-RUN] Would write fastp-clean sample manifest to ${FASTP_CLEAN_SAMPLE_MANIFEST}."
        return 0
    fi
    local candidate_manifest
    candidate_manifest="$(mktemp "${QC_DIR}/.fastp_clean_samples.XXXXXX")"
    build_sample_manifest_candidate "${QC_DIR}/trimmed" "fastp-clean" "${candidate_manifest}" || exit 1
    publish_sample_manifest "${candidate_manifest}" "${FASTP_CLEAN_SAMPLE_MANIFEST}" "${FASTP_CLEAN_SAMPLE_METADATA}" "fastp-clean" || exit 1
    rm -f "${candidate_manifest}"
    load_sample_manifest "${FASTP_CLEAN_SAMPLE_MANIFEST}" "fastp-clean" ACTIVE_SAMPLE_IDS ACTIVE_R1_FILES ACTIVE_R2_FILES || exit 1
    mark_done "01_fastp"
}

# ── Step 02: MEGAHIT Assembly ─────────────────────────────────────────────────
step_02_megahit() {
    if [[ "${RUN_MEGAHIT}" != true ]]; then
        log_info "[SKIP] MEGAHIT stage disabled."
        if [[ "${NEEDS_ASSEMBLY}" == true && -n "${ASSEMBLY_IN}" ]]; then
            CONTIGS="$(realpath -m "${ASSEMBLY_IN}")"
        fi
        return 0
    fi
    log_step "02" "Metagenome assembly with MEGAHIT"

    # Run MEGAHIT
    if is_done "02_megahit"; then
        log_info "[Step 02] Checkpoint found — skipping MEGAHIT."
        CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
        return 0
    fi

    local r1_list r2_list
    r1_list=$(IFS=','; echo "${ACTIVE_R1_FILES[*]}")
    r2_list=$(IFS=','; echo "${ACTIVE_R2_FILES[*]}")
    log_info "[Step 02] Using ${ACTIVE_READS_SOURCE} reads for assembly."

    local megahit_outdir="${ASSEMBLY_DIR}/megahit_raw"
    local mem_bytes=$(( MEMORY * 1073741824 ))
    local step_log="${LOGS_DIR}/step_02_megahit.log"

    # MEGAHIT will fail if output dir exists
    if [[ -d "${megahit_outdir}" ]]; then
        log_warn "[Step 02] Removing previous MEGAHIT output dir: ${megahit_outdir}"
        run_cmd rm -rf "${megahit_outdir}"
    fi

    run_cmd_with_tee "${step_log}" -- \
        run_tool "02" "${SIF_MEGAHIT}" "${SIF_MEGAHIT_URI}" megahit \
            -1 "${r1_list}" \
            -2 "${r2_list}" \
            -o "${megahit_outdir}" \
            -t "${THREADS}" \
            --min-contig-len "${MIN_CONTIG_LEN}" \
            -m "${mem_bytes}"

    local raw_contigs="${megahit_outdir}/final.contigs.fa"
    CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
    if [[ ! -s "${raw_contigs}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 02] Dry-run: FASTA filtering is deferred until MEGAHIT produces contigs."
            mark_done "02_megahit"
            return 0
        fi
        log_error "[Step 02] MEGAHIT contigs not found: ${raw_contigs}"
        exit 1
    fi

    log_info "[Step 02] Cleaning FASTA headers and filtering contigs < ${MIN_CONTIG_LEN} bp..."

    awk -v minlen="${MIN_CONTIG_LEN}" '
        BEGIN { header=""; seq="" }
        /^>/ {
            if (header != "" && length(seq) >= minlen) { print header; print seq }
            header = ">" $1; sub(/^>>/, ">", header); seq = ""; next
        }
        { seq = seq $0 }
        END { if (header != "" && length(seq) >= minlen) { print header; print seq } }
    ' "${raw_contigs}" > "${CONTIGS}"

    local contig_count
    contig_count=$(grep -c '^>' "${CONTIGS}" || true)
    log_info "[Step 02] Assembly complete. Contigs >= ${MIN_CONTIG_LEN} bp: ${contig_count}"

    mark_done "02_megahit"
}

# ── Step 03: Kraken2 + Bracken ────────────────────────────────────────────────
step_03_kraken2_bracken() {
    if [[ "${RUN_KRAKEN2}" == false ]]; then
        log_info "[SKIP] Kraken2/Bracken stage disabled."
        return 0
    fi
    log_step "03" "Read taxonomic profiling with Kraken2/Bracken"

    if is_done "03_kraken2"; then
        log_info "[Step 03] Checkpoint found — skipping Kraken2/Bracken."
        return 0
    fi

    local -a r1_sources=("${ACTIVE_R1_FILES[@]}")
    local -a r2_sources=("${ACTIVE_R2_FILES[@]}")
    log_info "[Step 03] Using ${ACTIVE_READS_SOURCE} reads for Kraken2."

    if [[ ${#r1_sources[@]} -eq 0 ]]; then
        log_error "[Step 03] No read files available for Kraken2."
        exit 1
    fi

    mkdir -p "${KRAKEN2_DIR}"

    for i in "${!r1_sources[@]}"; do
        local r1="${r1_sources[$i]}"
        local r2="${r2_sources[$i]}"
        local base
        base=$(basename "${r1}" "${READ_SUFFIX_R1}")
        local kraken_out="${KRAKEN2_DIR}/kraken_output_${base}.txt"
        local kraken_report="${KRAKEN2_DIR}/kraken_report_${base}.txt"
        local bracken_out="${KRAKEN2_DIR}/bracken_${base}.txt"
        local kraken_log="${LOGS_DIR}/step_03_kraken2_${base}.log"
        local bracken_log="${LOGS_DIR}/step_03_bracken_${base}.log"

        log_info "[Step 03] Running Kraken2 for sample: ${base}"
        run_cmd_with_tee "${kraken_log}" -- \
            run_tool "03" "${SIF_KRAKEN2}" "${SIF_KRAKEN2_URI}" kraken2 \
                --db "${KRAKEN2_DB}" \
                --threads "${THREADS}" \
                --paired "${r1}" "${r2}" \
                --output "${kraken_out}" \
                --report "${kraken_report}"

        log_info "[Step 03] Running Bracken for sample: ${base}"
        run_cmd_with_tee "${bracken_log}" -- \
            run_tool "03" "${SIF_KRAKEN2}" "${SIF_KRAKEN2_URI}" bracken \
                -d "${KRAKEN2_DB}" \
                -i "${kraken_report}" \
                -o "${bracken_out}" \
                -r "${BRACKEN_READ_LEN}" \
                -l "${BRACKEN_TAX_LEVEL}" \
                -t "${BRACKEN_THRESHOLD}"
    done

    mark_done "03_kraken2"
}

# ── Step 04: Bowtie2 Mapping + JGI Depth ─────────────────────────────────────
step_04_mapping_depth() {
    if [[ "${RUN_MAPPING}" != true ]]; then
        log_info "[SKIP] Mapping stage disabled."
        return 0
    fi
    log_step "04" "Mapping — Bowtie2 + SAMtools + JGI depth"

    if is_done "04_mapping"; then
        log_info "[Step 04] Checkpoint found — repopulating SORTED_BAMS and DEPTH_FILE."
        collect_files_by_pattern "${MAPPING_DIR}" "*_sorted.bam" SORTED_BAMS
        DEPTH_FILE="${BINNING_DIR}/depth_matrix.txt"
        BAM_DIR_ACTIVE="${MAPPING_DIR}"
        JGI_DEPTH_FILE_ACTIVE="${DEPTH_FILE}"
        log_info "[Step 04] Repopulated ${#SORTED_BAMS[@]} BAM(s) from checkpoint."
        return 0
    fi

    if [[ -z "${CONTIGS}" ]] || [[ ! -f "${CONTIGS}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 04] Dry-run: mapping is deferred until an assembly is available."
            return 0
        fi
        log_error "[Step 04] CONTIGS file not found: ${CONTIGS}"
        exit 1
    fi

    mkdir -p "${MAPPING_DIR}"

    local bt2_index="${MAPPING_DIR}/contigs_idx"
    local build_log="${LOGS_DIR}/step_04_bowtie2_build.log"

    log_info "[Step 04] Building Bowtie2 index..."
    run_cmd_with_tee "${build_log}" -- \
        run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" bowtie2-build \
            --threads "${THREADS}" \
            "${CONTIGS}" \
            "${bt2_index}"

    local -a r1_sources=("${ACTIVE_R1_FILES[@]}")
    local -a r2_sources=("${ACTIVE_R2_FILES[@]}")
    log_info "[Step 04] Using ${ACTIVE_READS_SOURCE} reads for mapping."

    if [[ ${#r1_sources[@]} -eq 0 ]]; then
        log_error "[Step 04] No read files available for mapping."
        exit 1
    fi

    SORTED_BAMS=()

    for i in "${!r1_sources[@]}"; do
        local r1="${r1_sources[$i]}"
        local r2="${r2_sources[$i]}"
        local base
        base=$(basename "${r1}" "${READ_SUFFIX_R1}")
        local bam_sorted="${MAPPING_DIR}/${base}_sorted.bam"
        local bt2_log="${LOGS_DIR}/step_04_bowtie2_${base}.log"
        local q_bt2_index q_r1 q_r2 q_bt2_log q_bam
        printf -v q_bt2_index '%q' "${bt2_index}"
        printf -v q_r1 '%q' "${r1}"
        printf -v q_r2 '%q' "${r2}"
        printf -v q_bt2_log '%q' "${bt2_log}"
        printf -v q_bam '%q' "${bam_sorted}"

        log_info "[Step 04] Mapping sample: ${base}"
        run_tool_shell "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" "bowtie2 --very-sensitive-local -p ${THREADS} --time --reorder \
            -x ${q_bt2_index} \
            -1 ${q_r1} \
            -2 ${q_r2} \
            2>>${q_bt2_log} \
            | samtools view -@ ${THREADS} -bS - \
            | samtools sort -@ ${THREADS} -o ${q_bam} -"

        log_info "[Step 04] Indexing BAM: ${bam_sorted}"
        run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" samtools index -@ "${THREADS}" "${bam_sorted}"

        SORTED_BAMS+=("${bam_sorted}")
    done

    DEPTH_FILE="${BINNING_DIR}/depth_matrix.txt"
    BAM_DIR_ACTIVE="${MAPPING_DIR}"
    JGI_DEPTH_FILE_ACTIVE="${DEPTH_FILE}"
    local jgi_log="${LOGS_DIR}/step_04_jgi_depth.log"

    log_info "[Step 04] Running jgi_summarize_bam_contig_depths..."
    run_cmd_with_tee "${jgi_log}" -- \
        run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" jgi_summarize_bam_contig_depths \
            --outputDepth "${DEPTH_FILE}" \
            "${SORTED_BAMS[@]}"

    mark_done "04_mapping"
}

# ── Step 05a: Vamb Binning ────────────────────────────────────────────────────
step_05a_vamb() {
    if [[ "${RUN_BINNING}" != true ]]; then
        log_info "[Binning: Vamb] Binning disabled; skipping."
        return 0
    fi
    if ! binner_selected vamb; then
        log_info "[Binning: Vamb] Method not selected; skipping."
        return 0
    fi
    log_step "05a" "Binning — Vamb"
    [[ "${COVERAGE_SOURCE}" != external ]] || validate_external_jgi_depth || exit 1

    if is_done "05a_vamb"; then
        log_info "[Step 05a] Checkpoint found — skipping Vamb."
        return 0
    fi

    if [[ -z "${DEPTH_FILE}" ]] || [[ ! -f "${DEPTH_FILE}" ]]; then
        log_error "[Step 05a] Depth file not found: ${DEPTH_FILE}"
        exit 1
    fi

    validate_vamb_fasta "${CONTIGS}" || exit 1
    if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
        validate_vamb_binsplit_headers "${CONTIGS}" "${VAMB_BINSPLIT_SEPARATOR}" || exit 1
    fi

    mkdir -p "${VAMB_DIR}"

    log_info "[Step 05a] Converting JGI depth matrix to Vamb abundance TSV..."
    write_vamb_abundance_tsv "${DEPTH_FILE}" "${VAMB_ABUNDANCE_TSV}"
    validate_vamb_abundance_tsv "${VAMB_ABUNDANCE_TSV}" "${CONTIGS}" || exit 1

    log_info "[Step 05a] Vamb abundance TSV created: ${VAMB_ABUNDANCE_TSV}"
    if [[ -n "${VAMB_BINSPLIT_SEPARATOR}" ]]; then
        log_info "[Step 05a] Binsplitting enabled with separator '${VAMB_BINSPLIT_SEPARATOR}'."
    else
        log_info "[Step 05a] Binsplitting disabled."
    fi

    local vamb_script="${TMPSCRIPTS_DIR}/step_05a_vamb.sh"
    write_vamb_script "${vamb_script}"

    run_cmd_with_tee "${LOGS_DIR}/step_05a_vamb.log" -- \
        run_mamba_file "05a" "${ENV_VAMB}" "${vamb_script}"

    local vamb_bins_dir="${VAMB_DIR}/vamb_out/bins"
    if [[ "${DRY_RUN}" == true ]]; then
        mark_done "05a_vamb"
        return 0
    fi
    normalize_binner_bins vamb "${vamb_bins_dir}"
    mark_done "05a_vamb"
}

# ── Step 05b: MetaBinner Binning ──────────────────────────────────────────────
step_05b_metabinner() {
    if [[ "${RUN_BINNING}" != true ]]; then
        log_info "[Binning: MetaBinner] Binning disabled; skipping."
        return 0
    fi
    if ! binner_selected metabinner; then
        log_info "[Binning: MetaBinner] Method not selected; skipping."
        return 0
    fi
    log_step "05b" "Binning — MetaBinner"
    [[ "${COVERAGE_SOURCE}" != external ]] || validate_external_jgi_depth || exit 1

    if is_done "05b_metabinner"; then
        log_info "[Step 05b] Checkpoint found — skipping MetaBinner."
        return 0
    fi

    if [[ -z "${DEPTH_FILE}" ]] || [[ ! -f "${DEPTH_FILE}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 05b] Dry-run: MetaBinner is deferred until mapping produces a depth file."
            return 0
        fi
        log_error "[Step 05b] Depth file not found: ${DEPTH_FILE}"
        exit 1
    fi

    mkdir -p "${METABINNER_DIR}"
    if [[ -d "${METABINNER_DIR}/metabinner_out" ]]; then
        run_cmd rm -rf -- "${METABINNER_DIR}/metabinner_out"
    fi

    # Step 1: Convert JGI depth matrix to MetaBinner coverage_profile.tsv
    log_info "[Step 05b] Converting JGI depth matrix to MetaBinner coverage_profile.tsv..."
    awk 'BEGIN{OFS="\t"}
    NR==1{
        printf "contig_name"
        for(i=4; i<=NF; i+=2) {
            n=$i; gsub(/.*\//,"",n); gsub(/_sorted\.bam$/,"",n)
            printf "\t%s", n
        }
        printf "\n"; next
    }
    {
        printf "%s", $1
        for(i=4; i<=NF; i+=2) printf "\t%s", $i
        printf "\n"
    }' "${DEPTH_FILE}" > "${METABINNER_DIR}/coverage_profile.tsv"

    # Step 2: Symlink assembly
    log_info "[Step 05b] Symlinking assembly to MetaBinner directory..."
    ln -sf "${CONTIGS}" "${METABINNER_DIR}/assembly.fa"

    # Step 3: Write and run prep script
    local mb_prep_script="${TMPSCRIPTS_DIR}/step_05b_metabinner_prep.sh"
    cat > "${mb_prep_script}" <<MB_PREP_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
cd "${METABINNER_DIR}"
echo "[MetaBinner] Filtering short contigs..."
python "${METABINNER_SCRIPTS_DIR}/Filter_tooshort.py" "${METABINNER_DIR}/assembly.fa" "${MIN_CONTIG_LEN}"
echo "[MetaBinner] Computing 4-mer frequencies..."
python "${METABINNER_SCRIPTS_DIR}/gen_kmer.py" "${METABINNER_DIR}/assembly.fa" "${MIN_CONTIG_LEN}" 4
MB_PREP_SCRIPT

    chmod +x "${mb_prep_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_05b_metabinner_prep.log" -- \
        run_mamba_file "05b" "${ENV_METABINNER}" "${mb_prep_script}"

    # Step 4: Write embedded Python sync script
    log_info "[Step 05b] Writing sync_metabinner.py..."
    cat > "${METABINNER_DIR}/sync_metabinner.py" <<'SYNC_PYTHON'
#!/usr/bin/env python3
"""
sync_metabinner.py — Synchronise contig IDs across MetaBinner input files.

Usage:
    python sync_metabinner.py <base_dir> <filtered_fasta_basename> \
                              <coverage_basename> <kmer_basename>
"""

import sys
import os

def main():
    if len(sys.argv) != 5:
        print("Usage: sync_metabinner.py <base_dir> <filtered_fasta_basename> "
              "<coverage_basename> <kmer_basename>", file=sys.stderr)
        sys.exit(1)

    base_dir             = sys.argv[1]
    filtered_fasta_name  = sys.argv[2]
    coverage_name        = sys.argv[3]
    kmer_name            = sys.argv[4]

    filtered_fasta_path  = os.path.join(base_dir, filtered_fasta_name)
    coverage_path        = os.path.join(base_dir, coverage_name)
    kmer_path            = os.path.join(base_dir, kmer_name)

    # 1. Read k-mer IDs (skip header, first column comma-split)
    kmer_ids = []
    kmer_data = {}
    with open(kmer_path, 'r') as fh:
        header = fh.readline()  # skip header
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split(',')
            seq_id = parts[0]
            kmer_ids.append(seq_id)
            kmer_data[seq_id] = line

    # 2. Read coverage IDs (skip header, first column tab-split)
    coverage_header = None
    coverage_ids = []
    coverage_data = {}
    with open(coverage_path, 'r') as fh:
        coverage_header = fh.readline().rstrip('\n')
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            seq_id = parts[0]
            coverage_ids.append(seq_id)
            coverage_data[seq_id] = line

    # 3. Compute intersection preserving order from kmer_ids
    coverage_set = set(coverage_ids)
    valid_ids = [sid for sid in kmer_ids if sid in coverage_set]
    valid_set = set(valid_ids)

    print(f"[sync_metabinner] k-mer IDs     : {len(kmer_ids)}")
    print(f"[sync_metabinner] Coverage IDs  : {len(coverage_ids)}")
    print(f"[sync_metabinner] Intersection  : {len(valid_ids)}")

    # 4. Write synced_kmer_no_header.csv (without header, only valid IDs)
    out_kmer = os.path.join(base_dir, 'synced_kmer_no_header.csv')
    with open(out_kmer, 'w') as fh:
        for sid in valid_ids:
            fh.write(kmer_data[sid] + '\n')
    print(f"[sync_metabinner] Written: {out_kmer}")

    # 5. Write synced_coverage.tsv (with header, only valid IDs)
    out_cov = os.path.join(base_dir, 'synced_coverage.tsv')
    with open(out_cov, 'w') as fh:
        fh.write(coverage_header + '\n')
        for sid in valid_ids:
            fh.write(coverage_data[sid] + '\n')
    print(f"[sync_metabinner] Written: {out_cov}")

    # 6. Write synced_contigs.fa (only valid IDs, header stripped to >seq_id)
    out_fa = os.path.join(base_dir, 'synced_contigs.fa')
    write_seq = False
    with open(filtered_fasta_path, 'r') as fh_in, open(out_fa, 'w') as fh_out:
        for line in fh_in:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                if seq_id in valid_set:
                    fh_out.write(f'>{seq_id}\n')
                    write_seq = True
                else:
                    write_seq = False
            else:
                if write_seq:
                    fh_out.write(line)
    print(f"[sync_metabinner] Written: {out_fa}")
    print("[sync_metabinner] Sync complete.")

if __name__ == '__main__':
    main()
SYNC_PYTHON

    # Step 5: Write and run MetaBinner run script
    local mb_run_script="${TMPSCRIPTS_DIR}/step_05b_metabinner_run.sh"
    cat > "${mb_run_script}" <<MB_RUN_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
cd "${METABINNER_DIR}"
echo "[MetaBinner] Syncing contig IDs..."
python "${METABINNER_DIR}/sync_metabinner.py" \\
    "${METABINNER_DIR}" \\
    "assembly_${MIN_CONTIG_LEN}.fa" \\
    "coverage_profile.tsv" \\
    "assembly_${MIN_CONTIG_LEN}_kmer_4_f${MIN_CONTIG_LEN}.csv"
echo "[MetaBinner] Re-attaching k-mer header..."
head -n 1 "assembly_${MIN_CONTIG_LEN}_kmer_4_f${MIN_CONTIG_LEN}.csv" > synced_kmer.csv
cat synced_kmer_no_header.csv >> synced_kmer.csv
rm synced_kmer_no_header.csv
echo "[MetaBinner] Running run_metabinner.sh..."
MBPATH=\$(dirname \$(which python))
bash "${METABINNER_RUNNER}" \\
    -a "${METABINNER_DIR}/synced_contigs.fa" \\
    -o "${METABINNER_DIR}/metabinner_out" \\
    -d "${METABINNER_DIR}/synced_coverage.tsv" \\
    -k "${METABINNER_DIR}/synced_kmer.csv" \\
    -p "\${MBPATH}" \\
    -t "${THREADS}"
MB_RUN_SCRIPT

    chmod +x "${mb_run_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_05b_metabinner_run.log" -- \
        run_mamba_file "05b" "${ENV_METABINNER}" "${mb_run_script}"

    local mb_bins_dir="${METABINNER_DIR}/metabinner_out/result_bins"
    if [[ "${DRY_RUN}" == true ]]; then
        mark_done "05b_metabinner"
        return 0
    fi
    normalize_binner_bins metabinner "${mb_bins_dir}"
    mark_done "05b_metabinner"
}

prepare_comebin_bam_dir() {
    local sample bam bam_index i destination="${COMEBIN_DIR}/external_bams"
    [[ "${COVERAGE_SOURCE}" == external ]] || { printf '%s\n' "${BAM_DIR_ACTIVE}"; return 0; }
    mkdir -p "${destination}"
    find "${destination}" -mindepth 1 -maxdepth 1 -type l -delete
    for (( i=0; i<${#ACTIVE_BAMS[@]}; i++ )); do
        sample="${ACTIVE_BAM_SAMPLE_IDS[$i]}"
        bam="${ACTIVE_BAMS[$i]}"
        bam_index="${ACTIVE_BAM_INDEXES[$i]}"
        ln -sfn -- "${bam}" "${destination}/${sample}_sorted.bam"
        ln -sfn -- "${bam_index}" "${destination}/${sample}_sorted.bam.$([[ "${bam_index}" == *.csi ]] && printf csi || printf bai)"
    done
    printf '%s\n' "${destination}"
}

# ── Step 05c: COMEBin Binning ─────────────────────────────────────────────────
step_05c_comebin() {
    if [[ "${RUN_BINNING}" != true ]]; then
        log_info "[Binning: COMEBin] Binning disabled; skipping."
        return 0
    fi
    if ! binner_selected comebin; then
        log_info "[Binning: COMEBin] Method not selected; skipping."
        return 0
    fi
    log_step "05c" "Binning — COMEBin"

    if is_done "05c_comebin"; then
        log_info "[Step 05c] Checkpoint found — skipping COMEBin."
        return 0
    fi

    if [[ ! -f "${CONTIGS}" || ! -f "${DEPTH_FILE}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 05c] Dry-run: COMEBin is deferred until assembly and mapping outputs are available."
            return 0
        fi
        log_error "[Step 05c] Assembly or mapping output is missing (CONTIGS=${CONTIGS}, DEPTH_FILE=${DEPTH_FILE})."
        exit 1
    fi

    mkdir -p "${COMEBIN_DIR}"
    if [[ -d "${COMEBIN_DIR}/comebin_res" ]]; then
        run_cmd rm -rf -- "${COMEBIN_DIR}/comebin_res"
    fi

    local comebin_script="${TMPSCRIPTS_DIR}/step_05c_comebin.sh"
    COMEBIN_BAM_DIR="$(prepare_comebin_bam_dir)"
    write_comebin_script "${comebin_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_05c_comebin.log" -- \
        run_mamba_file "05c" "${ENV_COMEBIN}" "${comebin_script}"

    local cb_bins_dir="${COMEBIN_DIR}/comebin_res/comebin_res_bins"
    if [[ "${DRY_RUN}" == true ]]; then
        mark_done "05c_comebin"
        return 0
    fi
    normalize_binner_bins comebin "${cb_bins_dir}"
    mark_done "05c_comebin"
}

# ── Step 06: MetaWrap Bin Refinement ──────────────────────────────────────────
step_06_bin_refinement() {
    if [[ "${RUN_REFINEMENT_EFFECTIVE}" != true ]]; then
        log_info "[SKIP] Refinement stage disabled."
        return 0
    fi
    log_step "06" "Refinement — MetaWrap bin_refinement"

    local method bins_dir
    local -a refinement_inputs=()
    for method in "${SELECTED_BINNER_ARRAY[@]}"; do
        bins_dir="$(get_binner_bins_dir "${method}")" || exit 1
        if [[ ! -d "${bins_dir}" ]]; then
            if [[ "${DRY_RUN}" == true ]]; then
                log_info "[Step 06] Dry-run: refinement is deferred until ${method} bins exist."
                return 0
            fi
            log_error "[Step 06] Selected ${method} bins are not available: ${bins_dir}"
            exit 1
        fi
        validate_bins_dir "${bins_dir}" "${method} normalized bins" || exit 1
        refinement_inputs+=("${bins_dir}")
    done

    if is_done "06_refinement"; then
        log_info "[Step 06] Checkpoint found — skipping bin refinement."
        return 0
    fi

    mkdir -p "${REFINEMENT_DIR}"
    if [[ -d "${REFINED_BINS_DIR}" ]]; then
        run_cmd rm -rf -- "${REFINED_BINS_DIR}"
    fi
    local refinement_script="${TMPSCRIPTS_DIR}/step_06_bin_refinement.sh"
    local -a refinement_command=(metawrap bin_refinement -o "${REFINEMENT_DIR}"
        -t "${THREADS}" -c "${METAWRAP_COMPLETENESS}" -x "${METAWRAP_CONTAMINATION}")
    local index=0 input_flag
    for bins_dir in "${refinement_inputs[@]}"; do
        case "${index}" in
            0) input_flag=-A ;;
            1) input_flag=-B ;;
            2) input_flag=-C ;;
            *) log_error "[Step 06] More than three input collections are unsupported."; exit 1 ;;
        esac
        refinement_command+=("${input_flag}" "${bins_dir}")
        index=$((index + 1))
    done
    refinement_command+=(--quick)
    {
        printf '#!/usr/bin/env bash\nset -euo pipefail\n'
        printf '%q ' "${refinement_command[@]}"
        printf '\n'
    } > "${refinement_script}"

    chmod +x "${refinement_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_06_bin_refinement.log" -- \
        run_mamba_file "06" "${ENV_METAWRAP}" "${refinement_script}"

    if [[ "${DRY_RUN}" == true ]]; then
        mark_done "06_refinement"
        return 0
    fi
    validate_bins_dir "${REFINED_BINS_DIR}" "refined bins" || exit 1
    if [[ "${BINS_SOURCE_RESOLVED}" != refined || "${BINS_FOR_QC}" != "${REFINED_BINS_DIR}" ]]; then
        log_error "[Step 06] Resolved bins source changed unexpectedly during refinement."
        exit 1
    fi
    mark_done "06_refinement"
}

# ── Step 07: MetaWrap Reassembly ──────────────────────────────────────────────
step_07_reassembly() {
    if [[ "${RUN_REASSEMBLY}" != true ]]; then
        log_info "[SKIP] Reassembly stage disabled."
        return 0
    fi
    log_step "07" "Reassembly — MetaWrap reassemble_bins"

    if is_done "07_reassembly"; then
        log_info "[Step 07] Checkpoint found — skipping reassembly."
        return 0
    fi

    if [[ -z "${BINS_FOR_REASSEMBLY}" || ! -d "${BINS_FOR_REASSEMBLY}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 07] Dry-run: reassembly is deferred until refined bins are available."
            return 0
        fi
        log_error "[Step 07] Selected bins directory not found: ${BINS_FOR_REASSEMBLY}"
        exit 1
    fi
    validate_bins_dir "${BINS_FOR_REASSEMBLY}" "MetaWRAP reassembly input bins" || exit 1

    mkdir -p "${REASSEMBLY_DIR}"

    # ── Create the per-run temporary workspace ────────────────────────────────
    if [[ -z "${REASSEMBLY_TMPDIR}" ]]; then
        log_error "[Step 07] REASSEMBLY_TMPDIR is not set. This is an internal error."
        exit 1
    fi
    create_pipeline_run_tmpdir || exit 1

    # ── Concatenate reads into REASSEMBLY_TMPDIR ──────────────────────────────
    local reassembly_r1="" reassembly_r2=""
    prepare_reassembly_tmpdir reassembly_r1 reassembly_r2 || exit 1

    log_info "[Step 07] Temporary reads:"
    log_info "          R1: ${reassembly_r1}"
    log_info "          R2: ${reassembly_r2}"
    log_info "[Step 07] Strategy      : ${REASSEMBLY_TEMP_READS_STRATEGY} v${REASSEMBLY_TEMP_READS_STRATEGY_VERSION}"
    log_info "[Step 07] Validation    : ${REASSEMBLY_TEMP_READS_VALIDATION} v${REASSEMBLY_TEMP_READS_VALIDATION_VERSION}"

    # ── Dispatch reassembly ───────────────────────────────────────────────────
    local reassembly_script="${TMPSCRIPTS_DIR}/step_07_reassembly.sh"

    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        # Build bind mounts; add the per-run tmpdir exclusively (read-write)
        local sif_path uri
        read -r sif_path uri < <(_env_to_sif "${ENV_METAWRAP}")
        local abs_sif
        abs_sif="$(realpath -m "${sif_path}" 2>/dev/null || echo "${sif_path}")"
        CURRENT_STAGE_RUNTIME="apptainer:${abs_sif}"

        ensure_sif_image "${sif_path}" "${uri}" || exit 1

        local -a bind_args=()
        build_bind_mounts bind_args
        # Add the isolated per-run tmpdir as an exclusive bind (same path inside)
        bind_args+=(-B "${PIPELINE_RUN_TMPDIR}:${PIPELINE_RUN_TMPDIR}:rw")

        local runtime="apptainer"
        command -v apptainer &>/dev/null || runtime="singularity"

        # Validate that the reads are visible and tmpdir is writable inside container
        log_info "[Step 07] Engine: Apptainer | Image: ${sif_path}"
        if [[ "${DRY_RUN}" != true ]]; then
            log_info "[Step 07] Validating container access to temporary workspace..."
            if ! "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" test -w "${REASSEMBLY_TMPDIR}"; then
                log_error "[Step 07] Container cannot write to temporary workspace: ${REASSEMBLY_TMPDIR}"
                exit 1
            fi
            if ! "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" test -r "${reassembly_r1}"; then
                log_error "[Step 07] Container cannot read temporary R1: ${reassembly_r1}"
                exit 1
            fi
            if ! "${runtime}" exec --cleanenv "${bind_args[@]}" "${sif_path}" test -r "${reassembly_r2}"; then
                log_error "[Step 07] Container cannot read temporary R2: ${reassembly_r2}"
                exit 1
            fi
        fi

        # Script references the in-container paths (identical due to same-path bind)
        cat > "${reassembly_script}" <<REASSEMBLY_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
metawrap reassemble_bins \\
    -o "${REASSEMBLY_DIR}" \\
    -b "${BINS_FOR_REASSEMBLY}" \\
    -1 "${reassembly_r1}" \\
    -2 "${reassembly_r2}" \\
    -t "${THREADS}" \\
    -m "${MEMORY}"
REASSEMBLY_SCRIPT
        chmod +x "${reassembly_script}"
        log_info "[Step 07] Script: ${reassembly_script}"

        CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[0]}"
        CURRENT_COMMAND_SOURCE_LINE="${LINENO}"
        CURRENT_STAGE_COMMAND="${runtime} exec --cleanenv ${bind_args[*]} ${sif_path} bash ${reassembly_script}"

        if [[ "${DRY_RUN}" == true ]]; then
            printf '[DRY-RUN] %q ' "${runtime}" exec --cleanenv "${bind_args[@]}" \
                --env "TMPDIR=${REASSEMBLY_TMPDIR}" \
                --env "TMP=${REASSEMBLY_TMPDIR}" \
                --env "TEMP=${REASSEMBLY_TMPDIR}" \
                "${sif_path}" bash "${reassembly_script}"
            printf '\n'
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${reassembly_script}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
            mark_done "07_reassembly"
            return 0
        fi

        local had_errexit=0 rc=0
        [[ $- == *e* ]] && had_errexit=1
        set +e
        run_cmd_with_tee "${LOGS_DIR}/step_07_reassembly.log" -- \
            env APPTAINERENV_TMPDIR="${REASSEMBLY_TMPDIR}" \
                APPTAINERENV_TMP="${REASSEMBLY_TMPDIR}" \
                APPTAINERENV_TEMP="${REASSEMBLY_TMPDIR}" \
                SINGULARITYENV_TMPDIR="${REASSEMBLY_TMPDIR}" \
                SINGULARITYENV_TMP="${REASSEMBLY_TMPDIR}" \
                SINGULARITYENV_TEMP="${REASSEMBLY_TMPDIR}" \
            "${runtime}" exec --cleanenv "${bind_args[@]}" \
                --env "TMPDIR=${REASSEMBLY_TMPDIR}" \
                --env "TMP=${REASSEMBLY_TMPDIR}" \
                --env "TEMP=${REASSEMBLY_TMPDIR}" \
                "${sif_path}" bash "${reassembly_script}"
        rc=$?
        (( had_errexit )) && set -e
        if (( rc != 0 )); then
            exit "${rc}"
        fi
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""

    else
        # Mamba mode: pass TMPDIR/TMP/TEMP only to the reassembly invocation
        CURRENT_STAGE_RUNTIME="mamba:${ENV_METAWRAP}"
        cat > "${reassembly_script}" <<REASSEMBLY_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
metawrap reassemble_bins \\
    -o "${REASSEMBLY_DIR}" \\
    -b "${BINS_FOR_REASSEMBLY}" \\
    -1 "${reassembly_r1}" \\
    -2 "${reassembly_r2}" \\
    -t "${THREADS}" \\
    -m "${MEMORY}"
REASSEMBLY_SCRIPT
        chmod +x "${reassembly_script}"
        log_info "[Step 07] Engine: mamba | Env: ${BOLD}${ENV_METAWRAP}${NC}"
        log_info "[Step 07] Script: ${reassembly_script}"

        CURRENT_COMMAND_SOURCE_FILE="${BASH_SOURCE[0]}"
        CURRENT_COMMAND_SOURCE_LINE="${LINENO}"
        CURRENT_STAGE_COMMAND="env TMPDIR=${REASSEMBLY_TMPDIR} TMP=${REASSEMBLY_TMPDIR} TEMP=${REASSEMBLY_TMPDIR} mamba run -n ${ENV_METAWRAP} bash ${reassembly_script}"

        if [[ "${DRY_RUN}" == true ]]; then
            echo -e "${YELLOW}[DRY-RUN]${NC} env TMPDIR=${REASSEMBLY_TMPDIR} TMP=${REASSEMBLY_TMPDIR} TEMP=${REASSEMBLY_TMPDIR} mamba run -n ${ENV_METAWRAP} bash ${reassembly_script}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${reassembly_script}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            CURRENT_STAGE_COMMAND=""
            CURRENT_COMMAND_SOURCE_FILE=""
            CURRENT_COMMAND_SOURCE_LINE=""
            mark_done "07_reassembly"
            return 0
        fi

        local had_errexit=0 rc=0
        [[ $- == *e* ]] && had_errexit=1
        set +e
        run_cmd_with_tee "${LOGS_DIR}/step_07_reassembly.log" -- \
            env TMPDIR="${REASSEMBLY_TMPDIR}" TMP="${REASSEMBLY_TMPDIR}" TEMP="${REASSEMBLY_TMPDIR}" \
                mamba run -n "${ENV_METAWRAP}" bash "${reassembly_script}"
        rc=$?
        (( had_errexit )) && set -e
        if (( rc != 0 )); then
            exit "${rc}"
        fi
        CURRENT_STAGE_COMMAND=""
        CURRENT_COMMAND_SOURCE_FILE=""
        CURRENT_COMMAND_SOURCE_LINE=""
    fi

    mark_done "07_reassembly"
    # Controlled success cleanup (after checkpoint written)
    cleanup_safe_temporaries "success"
}

# ── Step 08: CheckM2 Quality Assessment ───────────────────────────────────────
step_08_checkm2() {
    if [[ "${RUN_CHECKM2}" != true ]]; then
        log_info "[SKIP] CheckM2 stage disabled."
        return 0
    fi
    log_step "08" "Quality — CheckM2 predict"

    local expected_bins_dir
    case "${BINS_SOURCE_RESOLVED}" in
        vamb) expected_bins_dir="${VAMB_FINAL_BINS_DIR}" ;;
        metabinner) expected_bins_dir="${METABINNER_FINAL_BINS_DIR}" ;;
        comebin) expected_bins_dir="${COMEBIN_FINAL_BINS_DIR}" ;;
        refined) expected_bins_dir="${REFINED_BINS_DIR}" ;;
        *) log_error "[Step 08] Invalid resolved bins source: ${BINS_SOURCE_RESOLVED:-<empty>}"; exit 1 ;;
    esac
    [[ "${RUN_REASSEMBLY}" == true ]] && expected_bins_dir="${REASSEMBLY_DIR}/reassembled_bins"
    if [[ "${BINS_FOR_QC}" != "${expected_bins_dir}" ]]; then
        log_error "[Step 08] BINS_FOR_QC does not match the resolved bins source."
        exit 1
    fi
    if [[ ! -d "${BINS_FOR_QC}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 08] Dry-run: CheckM2 is deferred until selected bins are available."
            return 0
        fi
        log_error "[Step 08] Selected bins directory not found: ${BINS_FOR_QC}"
        exit 1
    fi
    validate_bins_dir "${BINS_FOR_QC}" "CheckM2 input bins" || exit 1

    if is_done "08_checkm2"; then
        log_info "[Step 08] Checkpoint found — skipping CheckM2."
        return 0
    fi

    if [[ ! -f "${CHECKM2_DB}" ]]; then
        log_error "[Step 08] CHECKM2_DB file not found: ${CHECKM2_DB}"
        exit 1
    fi

    mkdir -p "${CHECKM2_DIR}"

    local checkm2_script="${TMPSCRIPTS_DIR}/step_08_checkm2.sh"
    cat > "${checkm2_script}" <<CHECKM2_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
checkm2 predict \\
    --input "${BINS_FOR_QC}" \\
    -x fa \\
    --output-directory "${CHECKM2_DIR}" \\
    --threads "${THREADS}" \\
    --database_path "${CHECKM2_DB}"
CHECKM2_SCRIPT

    chmod +x "${checkm2_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_08_checkm2.log" -- \
        run_mamba_file "08" "${ENV_CHECKM2}" "${checkm2_script}"

    if [[ ! -f "${CHECKM2_DIR}/quality_report.tsv" ]]; then
        log_error "[Step 08] CheckM2 quality_report.tsv not found. CheckM2 may have failed."
        exit 1
    fi

    log_info "[Step 08] CheckM2 quality report: ${CHECKM2_DIR}/quality_report.tsv"
    mark_done "08_checkm2"
}

# ── Step 09: Filter MAGs ──────────────────────────────────────────────────────
step_09_filter_mags() {
    if [[ "${RUN_MAG_FILTERING}" != true ]]; then
        log_info "[SKIP] MAG filtering stage disabled."
        return 0
    fi
    log_step "09" "MAG filtering (completeness >= ${MIN_COMPLETENESS}%, contamination <= ${MAX_CONTAMINATION}%)"

    if is_done "09_filter_mags"; then
        log_info "[Step 09] Checkpoint found — skipping MAG filtering."
        HQ_BINS_DIR="${MAGS_DIR}"
        return 0
    fi

    if [[ ! -f "${CHECKM2_DIR}/quality_report.tsv" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 09] Dry-run: MAG filtering is deferred until CheckM2 produces a quality report."
            return 0
        fi
        log_error "[Step 09] CheckM2 quality_report.tsv not found: ${CHECKM2_DIR}/quality_report.tsv"
        exit 1
    fi

    mkdir -p "${MAGS_DIR}"

    local retained_count
    if ! retained_count="$(python3 - "${CHECKM2_DIR}/quality_report.tsv" \
            "${BINS_FOR_QC}" "${MAGS_DIR}" "${MAG_FILTER_MANIFEST}" \
            "${MIN_COMPLETENESS}" "${MAX_CONTAMINATION}" <<'PY_MAG_FILTER'
import csv
import logging
import os
import shutil
import sys
from decimal import Decimal, InvalidOperation

report, input_dir, output_dir, manifest_path, min_text, max_text = sys.argv[1:]
logging.basicConfig(level=logging.INFO, format="[MAG filtering] %(message)s")
min_completeness = Decimal(min_text)
max_contamination = Decimal(max_text)
required = ("Name", "Completeness", "Contamination")
records = []

try:
    with open(report, newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("CheckM2 quality report has a missing or duplicate header.")
        if any(name not in reader.fieldnames for name in required):
            raise ValueError("CheckM2 quality report must contain Name, Completeness, and Contamination columns.")
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"CheckM2 quality report has an inconsistent row at line {line_number}.")
            mag_id = (row.get("Name") or "").strip()
            comp_text = (row.get("Completeness") or "").strip()
            cont_text = (row.get("Contamination") or "").strip()
            completeness = None
            contamination = None
            reason = ""
            try:
                if not comp_text:
                    reason = "missing_completeness"
                elif not cont_text:
                    reason = "missing_contamination"
                else:
                    completeness = Decimal(comp_text)
                    contamination = Decimal(cont_text)
                    if (not completeness.is_finite() or not contamination.is_finite() or
                            completeness < 0 or completeness > 100 or
                            contamination < 0 or contamination > 100):
                        reason = "invalid_quality_values"
            except InvalidOperation:
                reason = "invalid_quality_values"

            if not reason:
                if completeness < min_completeness:
                    reason = "below_min_completeness"
                elif contamination > max_contamination:
                    reason = "above_max_contamination"
                else:
                    reason = "retained"

            source_path = os.path.join(input_dir, mag_id)
            if not source_path.endswith(".fa"):
                source_path += ".fa"
            retained = reason == "retained" and os.path.isfile(source_path)
            if reason == "retained" and not retained:
                reason = "missing_mag_file"
            output_path = os.path.join(output_dir, os.path.basename(source_path)) if retained else ""
            records.append({
                "mag_id": mag_id,
                "completeness": comp_text,
                "contamination": cont_text,
                "min_completeness": min_text,
                "max_contamination": max_text,
                "retained": "true" if retained else "false",
                "reason": reason,
                "source_path": source_path,
                "output_path": output_path,
            })
except (OSError, csv.Error, ValueError) as error:
    raise SystemExit(f"Invalid CheckM2 quality report: {error}")

if not records:
    raise SystemExit("Invalid CheckM2 quality report: no data rows.")

os.makedirs(output_dir, exist_ok=True)
for entry in os.scandir(output_dir):
    if entry.name.endswith(".fa") and (entry.is_file(follow_symlinks=False) or entry.is_symlink()):
        os.unlink(entry.path)
retained_count = 0
for record in records:
    if record["retained"] == "true":
        shutil.copy2(record["source_path"], record["output_path"])
        retained_count += 1
    elif record["reason"] in ("missing_completeness", "missing_contamination", "invalid_quality_values", "missing_mag_file"):
        logging.warning("%s excluded: %s", record["mag_id"], record["reason"])

with open(manifest_path, "w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=(
        "mag_id", "completeness", "contamination", "min_completeness",
        "max_contamination", "retained", "reason", "source_path", "output_path",
    ), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
print(retained_count)
PY_MAG_FILTER
    )"; then
        log_error "[Step 09] MAG filtering failed because the CheckM2 report structure or output operation is invalid."
        exit 1
    fi

    if [[ "${retained_count}" -eq 0 ]]; then
        log_error "[Step 09] Zero MAGs were retained (completeness >= ${MIN_COMPLETENESS}%, contamination <= ${MAX_CONTAMINATION}%). No successful filter checkpoint will be written."
        exit 1
    else
        log_info "[Step 09] ${retained_count} MAG(s) retained. Manifest: ${MAG_FILTER_MANIFEST}"
    fi

    HQ_BINS_DIR="${MAGS_DIR}"
    mark_done "09_filter_mags"
}

# ── Step 10: GTDB-Tk Taxonomy ─────────────────────────────────────────────────
step_10_gtdbtk() {
    if [[ "${RUN_GTDBTK}" != true ]]; then
        log_info "[SKIP] GTDB-Tk stage disabled."
        return 0
    fi
    log_step "10" "Taxonomy — GTDB-Tk classify_wf"

    if is_done "10_gtdbtk"; then
        log_info "[Step 10] Checkpoint found — skipping GTDB-Tk."
        return 0
    fi

    if [[ -z "${HQ_BINS_DIR}" || ! -d "${HQ_BINS_DIR}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 10] Dry-run: GTDB-Tk is deferred until selected MAGs are available."
            return 0
        fi
        log_error "[Step 10] Selected MAG directory not found: ${HQ_BINS_DIR}"
        exit 1
    fi

    # Count MAGs
    local mag_count
    mag_count=$(find "${HQ_BINS_DIR}" -name '*.fa' 2>/dev/null | wc -l)
    if [[ "${mag_count}" -eq 0 ]]; then
        log_warn "[Step 10] No MAGs found in ${HQ_BINS_DIR} — skipping GTDB-Tk."
        return 0
    fi
    log_info "[Step 10] Running GTDB-Tk on ${mag_count} MAG(s)."

    mkdir -p "${GTDBTK_DIR}/tmp"

    local gtdbtk_script="${TMPSCRIPTS_DIR}/step_10_gtdbtk.sh"
    cat > "${gtdbtk_script}" <<GTDBTK_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH}"
gtdbtk classify_wf \\
    --genome_dir "${HQ_BINS_DIR}" \\
    --out_dir    "${GTDBTK_DIR}" \\
    --extension  fa \\
    --tmpdir     "${GTDBTK_DIR}/tmp" \\
    --cpus       "${THREADS}" \\
    --pplacer_cpus 1
GTDBTK_SCRIPT

    chmod +x "${gtdbtk_script}"
    run_cmd_with_tee "${LOGS_DIR}/step_10_gtdbtk.log" -- \
        run_mamba_file "10" "${ENV_GTDBTK}" "${gtdbtk_script}"

    mark_done "10_gtdbtk"
}

# ── Step 11: CoverM Abundance ─────────────────────────────────────────────────
step_11_coverm() {
    if [[ "${RUN_COVERM}" != true ]]; then
        log_info "[SKIP] CoverM stage disabled."
        return 0
    fi
    log_step "11" "Abundance — CoverM genome"

    if is_done "11_coverm"; then
        log_info "[Step 11] Checkpoint found — skipping CoverM."
        return 0
    fi

    if [[ -z "${HQ_BINS_DIR}" || ! -d "${HQ_BINS_DIR}" ]]; then
        if [[ "${DRY_RUN}" == true ]]; then
            log_info "[Step 11] Dry-run: CoverM is deferred until selected MAGs are available."
            return 0
        fi
        log_error "[Step 11] Selected MAG directory not found: ${HQ_BINS_DIR}"
        exit 1
    fi

    # Count MAGs
    local mag_count
    mag_count=$(find "${HQ_BINS_DIR}" -name '*.fa' 2>/dev/null | wc -l)
    if [[ "${mag_count}" -eq 0 ]]; then
        log_warn "[Step 11] No MAGs found in ${HQ_BINS_DIR} — skipping CoverM."
        return 0
    fi
    log_info "[Step 11] Running CoverM on ${mag_count} MAG(s)."

    mkdir -p "${ABUNDANCE_DIR}"

    local -a r1_sources=("${ACTIVE_R1_FILES[@]}")
    local -a r2_sources=("${ACTIVE_R2_FILES[@]}")

    if [[ ${#r1_sources[@]} -eq 0 ]]; then
        log_error "[Step 11] No read files available for CoverM."
        exit 1
    fi

    for i in "${!r1_sources[@]}"; do
        local r1="${r1_sources[$i]}"
        local r2="${r2_sources[$i]}"
        local base
        base=$(basename "${r1}" "${READ_SUFFIX_R1}")
        local out_tsv="${ABUNDANCE_DIR}/abundance_${base}.tsv"
        local coverm_log="${LOGS_DIR}/step_11_coverm_${base}.log"

        # Build optional flags
        local -a methods_args=()
        local -a min_cov_args=()
        if [[ -n "${COVERM_METHODS}" ]]; then
            methods_args=(--methods "${COVERM_METHODS}")
        fi
        if [[ -n "${COVERM_MIN_COVERED_FRACTION}" ]]; then
            min_cov_args=(--min-covered-fraction "${COVERM_MIN_COVERED_FRACTION}")
        fi

        log_info "[Step 11] Running CoverM for sample: ${base}"
        run_cmd_with_tee "${coverm_log}" -- \
            run_tool "11" "${SIF_COVERM}" "${SIF_COVERM_URI}" coverm genome \
                --read1 "${r1}" \
                --read2 "${r2}" \
                --genome-fasta-directory "${HQ_BINS_DIR}" \
                --genome-fasta-extension fa \
                --threads "${THREADS}" \
                "${methods_args[@]}" \
                "${min_cov_args[@]}" \
                -o "${out_tsv}"
    done

    mark_done "11_coverm"
}

check_legacy_summary_checkpoint() {
    if [[ -f "${CHECKPOINT_DIR}/12_summary.done" || -f "${CHECKPOINT_DIR}/12_summary.config" ]]; then
        log_info "[Migration] Legacy monolithic 12_summary checkpoint detected. Marking as legacy; new summaries (12a_taxonomy_summary, 12b_mag_summary, 12c_run_overview) will be built independently without reusing the legacy checkpoint."
        touch "${CHECKPOINT_DIR}/12_summary.legacy" 2>/dev/null || true
    fi
}

# ── Step 12a: Taxonomy Summary ────────────────────────────────────────────────
step_12a_taxonomy_summary() {
    if [[ "${RUN_SUMMARY}" != true || "${RUN_TAXONOMY_SUMMARY}" != true ]]; then
        log_info "[SKIP] Taxonomy summary disabled."
        return 0
    fi
    if [[ "${RUN_KRAKEN2}" != true ]]; then
        log_info "[SKIP] Taxonomy summary not-applicable (Kraken2 disabled)."
        return 0
    fi
    log_step "12a" "Summary — Read-based taxonomy summary (Kraken2/Bracken)"

    if is_done "12a_taxonomy_summary"; then
        log_info "[Step 12a] Checkpoint found — skipping taxonomy summary."
        return 0
    fi

    mkdir -p "${SUMMARY_DIR}"
    local -a cmd=(
        python3 "${SUMMARY_GENERATOR}" taxonomy
        --outdir "${SUMMARY_DIR}"
        --pipeline-version "${PIPELINE_VERSION}"
        --sample-manifest "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}"
        --kraken2-dir "${KRAKEN2_DIR}"
        --kraken2-enabled
        --bracken-level "${BRACKEN_TAX_LEVEL:-S}"
    )
    if [[ "${RUN_BRACKEN:-true}" == true ]]; then
        cmd+=(--bracken-enabled)
    fi
    if [[ "${DRY_RUN}" == true ]]; then
        cmd+=(--dry-run)
        "${cmd[@]}" || true
        return 0
    fi

    local exit_code=0
    local had_errexit=0
    [[ $- == *e* ]] && had_errexit=1
    set +e
    "${cmd[@]}"
    exit_code=$?
    (( had_errexit )) && set -e

    if [[ "${exit_code}" -eq 0 ]]; then
        mark_done "12a_taxonomy_summary"
    elif [[ "${exit_code}" -eq 4 ]]; then
        log_info "[Step 12a] Taxonomy summary determined as not-applicable."
    else
        log_error "[Step 12a] Failed to generate taxonomy summary (exit code ${exit_code})."
        return "${exit_code}"
    fi
}

# ── Step 12b: MAG Summary ─────────────────────────────────────────────────────
step_12b_mag_summary() {
    if [[ "${RUN_SUMMARY}" != true || "${RUN_MAG_SUMMARY}" != true ]]; then
        log_info "[SKIP] MAG summary disabled."
        return 0
    fi
    if [[ "${RUN_CHECKM2}" != true && "${RUN_MAG_FILTERING}" != true && \
          "${RUN_GTDBTK}" != true && "${RUN_COVERM}" != true ]]; then
        log_info "[SKIP] MAG summary not-applicable (no MAG stages active)."
        return 0
    fi
    log_step "12b" "Summary — MAG quality, taxonomy, and abundance summary"

    if is_done "12b_mag_summary"; then
        log_info "[Step 12b] Checkpoint found — skipping MAG summary."
        return 0
    fi

    mkdir -p "${SUMMARY_DIR}"
    local -a cmd=(
        python3 "${SUMMARY_GENERATOR}" mag
        --outdir "${SUMMARY_DIR}"
        --pipeline-version "${PIPELINE_VERSION}"
        --bins-dir "${HQ_BINS_DIR:-${MAGS_DIR}}"
        --min-completeness "${MIN_COMPLETENESS}"
        --max-contamination "${MAX_CONTAMINATION}"
    )
    if [[ "${RUN_CHECKM2}" == true ]]; then
        cmd+=(--checkm2-enabled --checkm2-report "${CHECKM2_DIR}/quality_report.tsv")
    fi
    if [[ "${RUN_MAG_FILTERING}" == true ]]; then
        cmd+=(--filtering-enabled --filter-manifest "${MAG_FILTER_MANIFEST}")
    fi
    if [[ "${RUN_GTDBTK}" == true ]]; then
        cmd+=(--gtdbtk-enabled --gtdbtk-dir "${GTDBTK_DIR}")
    fi
    if [[ "${RUN_COVERM}" == true ]]; then
        cmd+=(
            --coverm-enabled
            --coverm-dir "${ABUNDANCE_DIR}"
            --sample-manifest "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}"
            --coverm-methods "${COVERM_METHODS:-relative_abundance}"
        )
    fi
    if [[ "${DRY_RUN}" == true ]]; then
        cmd+=(--dry-run)
        "${cmd[@]}" || true
        return 0
    fi

    local exit_code=0
    local had_errexit=0
    [[ $- == *e* ]] && had_errexit=1
    set +e
    "${cmd[@]}"
    exit_code=$?
    (( had_errexit )) && set -e

    if [[ "${exit_code}" -eq 0 ]]; then
        mark_done "12b_mag_summary"
    elif [[ "${exit_code}" -eq 4 ]]; then
        log_info "[Step 12b] MAG summary determined as not-applicable."
    else
        log_error "[Step 12b] Failed to generate MAG summary (exit code ${exit_code})."
        return "${exit_code}"
    fi
}

# ── Step 12c: Run Overview ────────────────────────────────────────────────────
step_12c_run_overview() {
    if [[ "${RUN_SUMMARY}" != true || "${RUN_RUN_OVERVIEW}" != true ]]; then
        log_info "[SKIP] Run overview disabled."
        return 0
    fi
    log_step "12c" "Summary — Pipeline run overview and execution metrics"

    if is_done "12c_run_overview"; then
        log_info "[Step 12c] Checkpoint found — skipping run overview."
        return 0
    fi

    mkdir -p "${SUMMARY_DIR}"
    local tax_stat="not-applicable"
    local mag_stat="not-applicable"
    [[ -s "${SUMMARY_DIR}/taxonomy_summary.tsv" ]] && tax_stat="available"
    [[ -s "${SUMMARY_DIR}/mag_summary.tsv" ]] && mag_stat="available"

    local -a cmd=(
        python3 "${SUMMARY_GENERATOR}" overview
        --outdir "${SUMMARY_DIR}"
        --pipeline-outdir "${OUTDIR}"
        --pipeline-version "${PIPELINE_VERSION}"
        --run-id "${RUN_ID}"
        --profile "${REQUESTED_PROFILE:-custom}"
        --execution-mode "${EXEC_MODE}"
        --resume-plan "${OUTDIR}/resume_plan.tsv"
        --run-status "${OUTDIR}/run_status.tsv"
        --run-manifest "${RUN_MANIFEST_FILE:-${LOGS_DIR}/run_manifest.tsv}"
        --stage-status "${STAGE_STATUS_FILE:-${LOGS_DIR}/stage_status.tsv}"
        --sample-manifest "${SAMPLE_MANIFEST_ACTIVE:-${FASTP_CLEAN_SAMPLE_MANIFEST:-${RAW_SAMPLE_MANIFEST}}}"
        --reads-source "${ACTIVE_READS_SOURCE:-raw}"
        --assembly-path "${CONTIGS:-}"
        --assembly-source "$([[ -n "${ASSEMBLY_IN}" ]] && echo external || echo megahit)"
        --bins-source-requested "${BINS_SOURCE:-}"
        --bins-source-resolved "${BINS_SOURCE_RESOLVED:-}"
        --external-bins-fingerprint "${EXTERNAL_BINS_FINGERPRINT:-}"
        --threads "${THREADS}"
        --min-completeness "${MIN_COMPLETENESS}"
        --max-contamination "${MAX_CONTAMINATION}"
        --comebin-device "${COMEBIN_DEVICE:-cpu}"
        --enabled-stages "${STAGE_IDS[*]}"
        --refinement-mode "${REFINEMENT_MODE:-auto}"
        --reassembly-enabled "${RUN_REASSEMBLY}"
        --checkpoint-dir "${CHECKPOINT_DIR}"
        --taxonomy-summary-status "${tax_stat}"
        --mag-summary-status "${mag_stat}"
        --software-inventory "${LOGS_DIR}/software_versions.tsv"
        --database-inventory "${LOGS_DIR}/database_versions.tsv"
        --log-reference "${LOGS_DIR}"
        --final-status "completed"
    )
    if [[ "${DRY_RUN}" == true ]]; then
        cmd+=(--dry-run)
        "${cmd[@]}" || true
        return 0
    fi

    local exit_code=0
    local had_errexit=0
    [[ $- == *e* ]] && had_errexit=1
    set +e
    "${cmd[@]}"
    exit_code=$?
    (( had_errexit )) && set -e

    if [[ "${exit_code}" -eq 0 ]]; then
        mark_done "12c_run_overview"
    else
        log_error "[Step 12c] Failed to generate run overview (exit code ${exit_code})."
        return "${exit_code}"
    fi
}

###############################################################################
#  MAIN                                                                        #
###############################################################################

main() {
    if [[ "${PLAN_MODE}" == true ]]; then
        resolve_binning_configuration || exit 2
        calculate_stage_needs
        validate_percentage_threshold MIN_COMPLETENESS "${MIN_COMPLETENESS}" || exit 2
        validate_percentage_threshold MAX_CONTAMINATION "${MAX_CONTAMINATION}" || exit 2
        if [[ -n "${STOP_AFTER}" ]] && ! stage_in_scope "${STOP_AFTER}"; then
            echo "[ERROR] Invalid --stop-after stage: ${STOP_AFTER}" >&2
            exit 2
        fi
        if [[ "${RUN_MEGAHIT}" == true && -n "${ASSEMBLY_IN}" ]]; then
            echo "[ERROR] --assembly-in cannot be combined with MEGAHIT; pass --no-megahit explicitly." >&2
            exit 2
        fi
        if [[ -n "${RAW_READS_DIR}" && -n "${CLEAN_READS_DIR}" ]]; then
            echo "[ERROR] --reads-dir and --clean-reads-dir are mutually exclusive." >&2
            exit 2
        fi
        if [[ -n "${CLEAN_READS_DIR}" && "${RUN_FASTP}" == true ]]; then
            echo "[ERROR] --clean-reads-dir cannot be combined with fastp; use --no-fastp." >&2
            exit 2
        fi

        collect_execution_plan || exit $?

        local tmp_plan
        tmp_plan="$(mktemp "${TMPDIR:-/tmp}/metanexus_plan.XXXXXX.json")"
        if ! serialize_execution_plan_json > "${tmp_plan}"; then
            rm -f "${tmp_plan}"
            echo "[ERROR] Failed to serialize execution plan to JSON." >&2
            exit 3
        fi

        local render_rc=0
        python3 "${PIPELINE_SCRIPTS_DIR}/execution_plan.py" render \
            --format "${PLAN_FORMAT}" \
            --details "${PLAN_DETAILS}" \
            --input "${tmp_plan}" || render_rc=$?

        if (( render_rc != 0 )); then
            rm -f "${tmp_plan}"
            exit "${render_rc}"
        fi

        local is_exec_rc=0
        python3 "${PIPELINE_SCRIPTS_DIR}/execution_plan.py" is-executable \
            --input "${tmp_plan}" || is_exec_rc=$?

        rm -f "${tmp_plan}"

        if [[ "${STRICT_PLAN}" == true ]]; then
            if (( is_exec_rc != 0 )); then
                exit 1
            fi
        fi
        exit 0
    fi

    resolve_binning_configuration || true
    MAG_FILTERING_APPLICABLE=false
    if [[ "${RUN_MAG_FILTERING}" == true ]] && stage_in_scope mag-filtering; then
        MAG_FILTERING_APPLICABLE=true
    fi
    # Auto-detect / validate compute resources before printing the banner
    detect_resources

    # Print run configuration banner
    echo -e "\n${CYAN}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}${BOLD}║              METANEXUS — PIPELINE CONFIGURATION              ║${NC}"
    echo -e "${CYAN}${BOLD}╠══════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${CYAN}${BOLD}║${NC}  Pipeline               : ${PIPELINE_NAME} v${PIPELINE_VERSION} (${PIPELINE_RELEASE_STATUS})"
    echo -e "${CYAN}${BOLD}║${NC}  PROFILE requested/result: ${REQUESTED_PROFILE} / ${RESULTING_PROFILE}"
    echo -e "${CYAN}${BOLD}║${NC}  READS_DIR (raw)        : ${READS_DIR:-<unset>}"
    echo -e "${CYAN}${BOLD}║${NC}  CLEAN_READS_DIR        : ${CLEAN_READS_DIR:-<unset>}"
    echo -e "${CYAN}${BOLD}║${NC}  ACTIVE_READS_SOURCE    : ${ACTIVE_READS_SOURCE:-<resolved in preflight>}"
    echo -e "${CYAN}${BOLD}║${NC}  OUTDIR                 : ${OUTDIR}"
    echo -e "${CYAN}${BOLD}║${NC}  EXTERNAL ASSEMBLY      : ${ASSEMBLY_IN:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  THREADS                : ${THREADS} (de ${_PHYS_CORES} disponibles)"
    echo -e "${CYAN}${BOLD}║${NC}  MEMORY (GB)            : ${MEMORY} GB (de ${_PHYS_MEM_GB} GB disponibles)"
    echo -e "${CYAN}${BOLD}║${NC}  CPU_RESERVE            : ${CPU_RESERVE} cores"
    echo -e "${CYAN}${BOLD}║${NC}  MEM_FRACTION           : ${MEM_FRACTION}%"
    echo -e "${CYAN}${BOLD}║${NC}  EXEC_MODE              : ${EXEC_MODE}"
    echo -e "${CYAN}${BOLD}║${NC}  SIF_DIR                : ${SIF_DIR}"
    echo -e "${CYAN}${BOLD}║${NC}  AUTO_PULL_SIF          : ${AUTO_PULL_SIF}"
    echo -e "${CYAN}${BOLD}║${NC}  MIN_CONTIG_LEN         : ${MIN_CONTIG_LEN} bp"
    echo -e "${CYAN}${BOLD}║${NC}  KRAKEN2_DB             : ${KRAKEN2_DB}"
    echo -e "${CYAN}${BOLD}║${NC}  CHECKM2_DB             : ${CHECKM2_DB}"
    echo -e "${CYAN}${BOLD}║${NC}  GTDBTK_DATA_PATH       : ${GTDBTK_DATA_PATH}"
    echo -e "${CYAN}${BOLD}║${NC}  Minimum completeness   : ${MIN_COMPLETENESS}"
    echo -e "${CYAN}${BOLD}║${NC}  Maximum contamination  : ${MAX_CONTAMINATION}"
    echo -e "${CYAN}${BOLD}║${NC}  Thresholds used        : ${MAG_FILTERING_APPLICABLE}"
    if [[ "${MAG_FILTERING_APPLICABLE}" == true ]]; then
        echo -e "${CYAN}${BOLD}║${NC}  Retention rule         : completeness >= ${MIN_COMPLETENESS} AND contamination <= ${MAX_CONTAMINATION}"
    fi
    echo -e "${CYAN}${BOLD}║${NC}  RUN_FASTP              : ${RUN_FASTP}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_MEGAHIT            : ${RUN_MEGAHIT}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_KRAKEN2            : ${RUN_KRAKEN2}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_REASSEMBLY         : ${RUN_REASSEMBLY}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_GTDBTK             : ${RUN_GTDBTK}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_MAPPING            : ${RUN_MAPPING}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_BINNING            : ${RUN_BINNING}"
    echo -e "${CYAN}${BOLD}║${NC}  SELECTED_BINNERS       : ${SELECTED_BINNERS:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  REFINEMENT_MODE        : ${REFINEMENT_MODE}"
    echo -e "${CYAN}${BOLD}║${NC}  REFINEMENT_EFFECTIVE   : ${RUN_REFINEMENT_EFFECTIVE}"
    echo -e "${CYAN}${BOLD}║${NC}  BINS_SOURCE requested  : ${BINS_SOURCE}"
    echo -e "${CYAN}${BOLD}║${NC}  BINS_SOURCE resolved   : ${BINS_SOURCE_RESOLVED:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  BINS_FOR_QC            : ${BINS_FOR_QC:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  BINS_FOR_REASSEMBLY    : ${BINS_FOR_REASSEMBLY:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_CHECKM2            : ${RUN_CHECKM2}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_MAG_FILTERING      : ${RUN_MAG_FILTERING}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_COVERM             : ${RUN_COVERM}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_SUMMARY            : ${RUN_SUMMARY}"
    echo -e "${CYAN}${BOLD}║${NC}  STOP_AFTER             : ${STOP_AFTER:-<none>}"
    echo -e "${CYAN}${BOLD}║${NC}  DRY_RUN                : ${DRY_RUN}"
    echo -e "${CYAN}${BOLD}║${NC}  FORCE_RERUN            : ${FORCE_RERUN}"
    echo -e "${CYAN}${BOLD}║${NC}  RESUME_CONTROL         : ${RESUME_CONTROL}"
    echo -e "${CYAN}${BOLD}║${NC}  TMPDIR_CLI             : ${PIPELINE_TMPDIR_CLI:-<unset>}"
    echo -e "${CYAN}${BOLD}║${NC}  KEEP_TEMPORARY         : ${KEEP_TEMPORARY}"
    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${CYAN}${BOLD}║${NC}  Error history          : would write on failure"
        echo -e "${CYAN}${BOLD}║${NC}  Last error report      : would write on failure"
        echo -e "${CYAN}${BOLD}║${NC}  Checkpoint state       : would track in normal execution"
        echo -e "${CYAN}${BOLD}║${NC}  Disk report            : would write"
    else
        echo -e "${CYAN}${BOLD}║${NC}  Error history          : ${LOGS_DIR}/errors.tsv"
        echo -e "${CYAN}${BOLD}║${NC}  Last error report      : ${LOGS_DIR}/last_error.json"
        echo -e "${CYAN}${BOLD}║${NC}  Checkpoint state       : ${LOGS_DIR}/checkpoint_state/${RUN_ID}.tsv"
        echo -e "${CYAN}${BOLD}║${NC}  Disk report            : pending"
        echo -e "${CYAN}${BOLD}║${NC}  Error log tail         : 20 lines"
        echo -e "${CYAN}${BOLD}║${NC}  Error schema           : 1"
    fi
    echo -e "${CYAN}${BOLD}╚══════════════════════════════════════════════════════════════╝${NC}\n"

    # Pre-flight check
    preflight_check
    build_inventory_rows
    if [[ "${DRY_RUN}" != true ]]; then
        mkdir -p "${OUTDIR}" "${LOGS_DIR}" "${TMPSCRIPTS_DIR}" "${CHECKPOINT_DIR}"
        [[ "${RUN_FASTP}" == true ]] && mkdir -p "${QC_DIR}"
        [[ "${RUN_MEGAHIT}" == true ]] && mkdir -p "${ASSEMBLY_DIR}"
        [[ "${RUN_KRAKEN2}" == true ]] && mkdir -p "${KRAKEN2_DIR}"
        [[ "${RUN_MAPPING}" == true ]] && mkdir -p "${MAPPING_DIR}" "${BINNING_DIR}"
        if [[ "${RUN_BINNING}" == true ]]; then
            mkdir -p "${BINNING_DIR}"
            binner_selected vamb && mkdir -p "${VAMB_DIR}"
            binner_selected metabinner && mkdir -p "${METABINNER_DIR}"
            binner_selected comebin && mkdir -p "${COMEBIN_DIR}"
        fi
        [[ "${RUN_REFINEMENT_EFFECTIVE}" == true ]] && mkdir -p "${REFINEMENT_DIR}"
        [[ "${RUN_REASSEMBLY}" == true ]] && mkdir -p "${REASSEMBLY_DIR}"
        [[ "${RUN_CHECKM2}" == true ]] && mkdir -p "${CHECKM2_DIR}"
        [[ "${RUN_MAG_FILTERING}" == true ]] && mkdir -p "${MAGS_DIR}"
        [[ "${RUN_GTDBTK}" == true ]] && mkdir -p "${GTDBTK_DIR}"
        [[ "${RUN_COVERM}" == true ]] && mkdir -p "${ABUNDANCE_DIR}"
        [[ "${RUN_SUMMARY}" == true ]] && mkdir -p "${SUMMARY_DIR}"
        persist_external_bam_manifest
    fi
    publish_inventory_generation
    apply_restart_cleanup

    if [[ "${DRY_RUN}" == true ]]; then
        local planned_id planned_action display internal_number
        for planned_id in "${STAGE_IDS[@]}"; do
            stage_is_active "${planned_id}" || continue
            planned_action="${STAGE_ACTION[$planned_id]}"
            display="$(get_stage_display_name "${planned_id}")"
            case "${planned_id}" in
                01_fastp) internal_number=01 ;;
                03_kraken2) internal_number=03 ;;
                02_megahit) internal_number=02 ;;
                04_mapping) internal_number=04 ;;
                05a_vamb) internal_number=05a ;;
                05b_metabinner) internal_number=05b ;;
                05c_comebin) internal_number=05c ;;
                06_refinement) internal_number=06 ;;
                07_reassembly) internal_number=07 ;;
                08_checkm2) internal_number=08 ;;
                09_filter_mags) internal_number=09 ;;
                10_gtdbtk) internal_number=10 ;;
                11_coverm) internal_number=11 ;;
                12a_taxonomy_summary) internal_number=12a ;;
                12b_mag_summary) internal_number=12b ;;
                12c_run_overview) internal_number=12c ;;
            esac
            log_step "${internal_number}" "[${planned_action}] ${display}"
            if [[ "${planned_action}" != reuse && "${planned_action}" != restart ]]; then
                log_info "[DRY-RUN] Would remove checkpoint and manifest for ${planned_id} immediately before execution."
                clean_stage_outputs "${planned_id}"
            fi
        done
        log_info "[DRY-RUN] No stage commands will execute and no result, checkpoint, plan, or status files will be changed."
        return 0
    fi

    # Execute all steps
    run_stage 01_fastp step_01_fastp
    run_stage 03_kraken2 step_03_kraken2_bracken
    run_stage 02_megahit step_02_megahit
    while true; do
        run_stage 04_mapping step_04_mapping_depth
        if stop_after_reached mapping; then log_info "Stopping successfully after mapping."; break; fi

        run_stage 05a_vamb step_05a_vamb
        run_stage 05b_metabinner step_05b_metabinner
        run_stage 05c_comebin step_05c_comebin
        if stop_after_reached binning; then log_info "Stopping successfully after binning."; break; fi

        run_stage 06_refinement step_06_bin_refinement
        if stop_after_reached refinement; then log_info "Stopping successfully after refinement."; break; fi

        run_stage 07_reassembly step_07_reassembly
        if stop_after_reached reassembly; then log_info "Stopping successfully after reassembly."; break; fi

        run_stage 08_checkm2 step_08_checkm2
        if stop_after_reached checkm2; then log_info "Stopping successfully after CheckM2."; break; fi

        run_stage 09_filter_mags step_09_filter_mags
        if stop_after_reached mag-filtering; then log_info "Stopping successfully after MAG filtering."; break; fi

        run_stage 10_gtdbtk step_10_gtdbtk
        if stop_after_reached gtdbtk; then log_info "Stopping successfully after GTDB-Tk."; break; fi

        run_stage 11_coverm step_11_coverm
        if stop_after_reached coverm; then log_info "Stopping successfully after CoverM."; break; fi

        run_stage 12a_taxonomy_summary step_12a_taxonomy_summary
        if stop_after_reached taxonomy-summary; then log_info "Stopping successfully after taxonomy summary."; break; fi

        run_stage 12b_mag_summary step_12b_mag_summary
        if stop_after_reached mag-summary; then log_info "Stopping successfully after MAG summary."; break; fi

        run_stage 12c_run_overview step_12c_run_overview
        if stop_after_reached run-overview || stop_after_reached summary; then log_info "Stopping successfully after run overview."; break; fi
        break
    done

    # Print final summary
    if [[ -n "${STOP_AFTER}" ]]; then
        echo -e "\n${GREEN}${BOLD}PIPELINE STOPPED SUCCESSFULLY AFTER ${STOP_AFTER} 🎉${NC}"
    else
        echo -e "\n${GREEN}${BOLD}PIPELINE COMPLETE 🎉${NC}"
    fi
    echo "Output root    : ${OUTDIR}"
    echo "Assembly       : ${CONTIGS}"
    echo "Mapped samples : ${#SORTED_BAMS[@]}"
    echo "Vamb bins      : $(find "${VAMB_FINAL_BINS_DIR}" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l)"
    echo "MetaBinner bins: $(find "${METABINNER_FINAL_BINS_DIR}" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l)"
    echo "COMEBin bins   : $(find "${COMEBIN_FINAL_BINS_DIR}" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l)"
    echo "Refined bins   : $(find "${REFINEMENT_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "Retained MAGs  : $(find "${MAGS_DIR}" -maxdepth 1 -name '*.fa' 2>/dev/null | wc -l)"
    
    local tax_status="not-applicable" mag_status="not-applicable" ov_status="pending"
    local tax_path="not-applicable" mag_path="not-applicable" ov_path="not-applicable"

    if [[ -s "${SUMMARY_DIR}/taxonomy_summary.tsv" ]]; then
        tax_status="available"
        tax_path="${SUMMARY_DIR}/taxonomy_summary.tsv"
    elif [[ "${STAGE_STATE[12a_taxonomy_summary]:-}" == "not-applicable" || "${RUN_TAXONOMY_SUMMARY}" != true || "${RUN_SUMMARY}" != true ]]; then
        tax_status="not-applicable"
        tax_path="not-applicable"
    elif [[ "${STAGE_FINAL_STATUS[12a_taxonomy_summary]:-}" == "failed" ]]; then
        tax_status="failed"
        tax_path="not-applicable"
    else
        tax_status="pending"
        tax_path="not-applicable"
    fi

    if [[ -s "${SUMMARY_DIR}/mag_summary.tsv" ]]; then
        mag_status="available"
        mag_path="${SUMMARY_DIR}/mag_summary.tsv"
    elif [[ "${STAGE_STATE[12b_mag_summary]:-}" == "not-applicable" || "${RUN_MAG_SUMMARY}" != true || "${RUN_SUMMARY}" != true ]]; then
        mag_status="not-applicable"
        mag_path="not-applicable"
    elif [[ "${STAGE_FINAL_STATUS[12b_mag_summary]:-}" == "failed" ]]; then
        mag_status="failed"
        mag_path="not-applicable"
    else
        mag_status="pending"
        mag_path="not-applicable"
    fi

    if [[ -s "${SUMMARY_DIR}/run_overview.tsv" ]]; then
        ov_status="available"
        ov_path="${SUMMARY_DIR}/run_overview.tsv"
    elif [[ "${STAGE_FINAL_STATUS[12c_run_overview]:-}" == "failed" ]]; then
        ov_status="failed"
        ov_path="${SUMMARY_DIR}/run_overview.tsv"
    else
        ov_status="pending"
        ov_path="${SUMMARY_DIR}/run_overview.tsv"
    fi

    echo "Taxonomy summary       : ${tax_status}"
    echo "MAG summary            : ${mag_status}"
    echo "Run overview           : ${ov_status}"
    echo "Taxonomy summary path  : ${tax_path}"
    echo "MAG summary path       : ${mag_path}"
    echo "Run overview path      : ${ov_path}"
    # Final disk space observations
    record_final_disk_observations

    local disk_sum_out
    disk_sum_out="$(get_latest_disk_status)"
    DISK_WARNING_COUNT="$(grep -E '^WARNING_COUNT=' <<< "${disk_sum_out}" | cut -d'=' -f2)"
    [[ -z "${DISK_WARNING_COUNT}" || ! "${DISK_WARNING_COUNT}" =~ ^[0-9]+$ ]] && DISK_WARNING_COUNT=0

    if [[ "${DRY_RUN}" != true && "${RUN_STATUS_INITIALIZED:-false}" == true ]]; then
        local total_end_mono
        total_end_mono="$(python3 -c 'import time; print(f"{time.monotonic():.6f}")' 2>/dev/null || echo "0.0")"
        local total_wall="NA"
        if [[ -n "${PIPELINE_START_MONOTONIC:-}" ]]; then
            total_wall="$(python3 -c "print(f'{(float(\"${total_end_mono}\") - float(\"${PIPELINE_START_MONOTONIC}\")):.3f}')" 2>/dev/null || echo "NA")"
        fi
        local err_count=0
        if [[ -f "${LOGS_DIR}/errors.tsv" ]]; then
            err_count="$(awk 'NR>1 {count++} END {print count+0}' "${LOGS_DIR}/errors.tsv" 2>/dev/null || echo 0)"
        fi
        python3 "${BENCHMARK_WRITER}" finalize-run \
            --run-manifest-path "${RUN_MANIFEST_FILE:-${LOGS_DIR}/run_manifest.tsv}" \
            --run-id "${RUN_ID}" \
            --final-status "success" \
            --exit-code 0 \
            --wall-time-seconds "${total_wall}" \
            --warning-count "${DISK_WARNING_COUNT}" \
            --error-count "${err_count}" \
            --taxonomy-summary-path "${SUMMARY_DIR}/taxonomy_summary.tsv" \
            --mag-summary-path "${SUMMARY_DIR}/mag_summary.tsv" \
            --run-overview-path "${SUMMARY_DIR}/run_overview.tsv" >/dev/null 2>&1 || true
    fi

    echo "Run manifest   : ${LOGS_DIR}/run_manifest.tsv"
    echo "Stage status   : ${LOGS_DIR}/stage_status.tsv"
    echo "Run status     : ${OUTDIR}/run_status.tsv"
    echo "Disk report    : ${DISK_SPACE_TSV}"
    echo "Logs           : ${LOGS_DIR}/"
    echo "Mamba scripts  : ${TMPSCRIPTS_DIR}/"
    echo "Error history  : ${LOGS_DIR}/errors.tsv"
    echo "Last error     : ${LOGS_DIR}/last_error.json"
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
