#!/usr/bin/env bash
###############################################################################
# run_metagenomics_pipeline.sh — Consolidated Shotgun Metagenomics Pipeline
#
# Workflow:
#   01  QC            fastp                          [optional: RUN_FASTP]
#   02  Assembly      MEGAHIT                        [optional: RUN_MEGAHIT]
#   03  Taxonomy      Kraken2 + Bracken (read-level) [optional: RUN_KRAKEN2]
#   04  Mapping       Bowtie2 + SAMtools + jgi_summarize_bam_contig_depths
#   05a Binning       Vamb (uses shared depth_matrix.txt)
#   05b Binning       MetaBinner (reuses Step 04 BAMs via JGI depth conversion)
#   05c Binning       COMEBin   (reuses Step 04 MAPPING_DIR)
#   06  Refinement    MetaWrap bin_refinement (-A Vamb -B COMEBin -C MetaBinner)
#   07  Reassembly    MetaWrap reassemble_bins       [optional: RUN_REASSEMBLY]
#   08  Quality       CheckM2 predict
#   09  Filter        HQ/MQ MAG selection
#   10  Taxonomy      GTDB-Tk classify_wf            [optional: RUN_GTDBTK]
#   11  Abundance     CoverM genome
#   12  Summary       R: left-join TSV (base R only, no tidyverse)
#
# Usage:
#   1. Edit all variables in the CONFIGURATION block below.
#   2. bash run_metagenomics_pipeline.sh [--dry-run] [--force] [-h|--help]
#
# Execution model:
#   Container steps: fastp, megahit, kraken2, bracken, bowtie2, samtools,
#                     jgi_summarize_bam_contig_depths, vamb, metabinner,
#                     comebin, metawrap, checkm2, gtdbtk, coverm, Rscript
#   mamba mode     : all tools retain their original host/mamba execution paths
###############################################################################

set -euo pipefail
trap 'echo -e "\n[FATAL] Pipeline aborted at line ${LINENO}. Check logs in ${OUTDIR}/logs/" >&2' ERR

###############################################################################
#  CONFIGURATION  —  Edit everything in this block before running             #
###############################################################################

# ── Input / Output ────────────────────────────────────────────────────────────
READS_DIR="${READS_DIR:-}"      # Dir containing paired-end reads
OUTDIR="${OUTDIR:-./metagenomics_out}"     # All pipeline output goes here
ASSEMBLY_IN=""                   # Pre-assembled FASTA. Empty = run MEGAHIT.

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

SIF_R="${SIF_DIR}/r_base.sif"
SIF_R_URI="docker://rocker/r-base:latest"

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

# ── Assembly Parameters ───────────────────────────────────────────────────────
MIN_CONTIG_LEN=1000    # bp — contigs shorter than this are discarded

# ── Kraken2 + Bracken Parameters ─────────────────────────────────────────────
BRACKEN_READ_LEN=150
BRACKEN_TAX_LEVEL="S"   # S=Species G=Genus F=Family O=Order C=Class P=Phylum
BRACKEN_THRESHOLD=10

# ── MAG Quality Filtering (CheckM2) ──────────────────────────────────────────
COMPLETENESS=50          # minimum % completeness
CONTAMINATION=10         # maximum % contamination

# ── MetaWrap Refinement Thresholds ────────────────────────────────────────────
METAWRAP_COMPLETENESS=50    # passed to bin_refinement -c
METAWRAP_CONTAMINATION=10   # passed to bin_refinement -x

# ── CoverM Parameters ─────────────────────────────────────────────────────────
COVERM_METHODS="relative_abundance"   # leave empty for CoverM built-in default
COVERM_MIN_COVERED_FRACTION="0.0"

# ── Pipeline Step Flags ───────────────────────────────────────────────────────
RUN_FASTP="${RUN_FASTP:-true}"                # Step 01
RUN_MEGAHIT="${RUN_MEGAHIT:-true}"            # Step 02 — false requires ASSEMBLY_IN to be set
RUN_KRAKEN2="${RUN_KRAKEN2:-true}"            # Step 03
RUN_REASSEMBLY="${RUN_REASSEMBLY:-false}"     # Step 07
RUN_GTDBTK="${RUN_GTDBTK:-true}"              # Step 10

# ── Execution Flags ───────────────────────────────────────────────────────────
DRY_RUN="${DRY_RUN:-false}"                   # Print commands without running them
FORCE_RERUN="${FORCE_RERUN:-false}"           # Ignore checkpoints and rerun every step

###############################################################################
#  END OF CONFIGURATION                                                       #
###############################################################################

# ── CLI Parsing ───────────────────────────────────────────────────────────────
_show_help() {
    cat <<'HELP'
run_metagenomics_pipeline.sh — Consolidated Shotgun Metagenomics Pipeline

Usage:
  bash run_metagenomics_pipeline.sh [OPTIONS]

Options:
  -r, --reads-dir <directory>  Directory containing paired-end reads
      --kraken2-db <directory> Kraken2 database directory
      --checkm2-db <file>      CheckM2 database file
      --gtdbtk-data-path <dir> GTDB-Tk reference data directory
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
  --megahit / --no-megahit     Enable/disable MEGAHIT (default: enabled)
  --kraken2 / --no-kraken2     Enable/disable Kraken2 + Bracken (default: enabled)
  --reassembly / --no-reassembly
                               Enable/disable MetaWrap reassembly (default: disabled)
  --gtdbtk / --no-gtdbtk       Enable/disable GTDB-Tk (default: enabled)
  --dry-run / --no-dry-run     Print commands without executing them (default: disabled)
  --force / --no-force         Ignore checkpoints (default: disabled)
  -h, --help                   Show this help message and exit

Example:
  bash run_metagenomics_pipeline.sh \\
    --reads-dir /data/reads \\
    --kraken2-db /databases/kraken2 \\
    --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \\
    --gtdbtk-data-path /databases/gtdbtk \\
    --outdir ./metagenomics_out

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

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--reads-dir)
            [[ -z "${2:-}" || "${2}" == -* ]] && {
                echo "[ERROR] $1 requires a directory value." >&2; exit 1; }
            READS_DIR="$2"; shift 2 ;;
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
        --reassembly)     RUN_REASSEMBLY=true;  shift ;;
        --no-reassembly)  RUN_REASSEMBLY=false; shift ;;
        --gtdbtk)     RUN_GTDBTK=true;    shift ;;
        --no-gtdbtk)  RUN_GTDBTK=false;   shift ;;
        --dry-run)     DRY_RUN=true;      shift ;;
        --no-dry-run)  DRY_RUN=false;     shift ;;
        --force)     FORCE_RERUN=true;   shift ;;
        --no-force)  FORCE_RERUN=false;  shift ;;
        -h|--help)  _show_help ;;
        *)  echo "[ERROR] Unknown option: $1. Use -h for help." >&2; exit 1 ;;
    esac
done

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

log_info()  { echo -e "${GREEN}[INFO]${NC}  $(date '+%Y-%m-%d %H:%M:%S')  $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $(date '+%Y-%m-%d %H:%M:%S')  $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S')  $*" >&2; }
log_step()  {
    local n="$1"; shift
    echo -e "\n${CYAN}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}${BOLD}║  [Step ${n}]  $*${NC}"
    echo -e "${CYAN}${BOLD}╚══════════════════════════════════════════════════════════════╝${NC}\n"
}

# ── Utility Functions ─────────────────────────────────────────────────────────
run_cmd() {
    if [[ "${DRY_RUN}" == true ]]; then
        printf '[DRY-RUN] '
        printf '%q ' "$@"
        printf '\n'
    else
        printf -v cmd_display '%q ' "$@"
        log_info "Running: ${cmd_display}"
        "$@"
    fi
}
# ── run_tool(step_label, sif_path, uri, command...) ───────────────────────────
# Runs a direct binary in its SIF image in Apptainer mode, preserving the
# original direct command path in mamba mode.
run_tool() {
    local step_label="$1"
    local sif_path="$2"
    local uri="$3"
    shift 3

    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        run_cmd "$@"
        return
    fi

    ensure_sif_image "${sif_path}" "${uri}" || return 1
    local runtime="apptainer"
    command -v apptainer &>/dev/null || runtime="singularity"
    local binds
    binds=$(build_bind_mounts)

    log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} ${runtime} exec --cleanenv ${binds} ${sif_path} $*"
        return 0
    fi

    # shellcheck disable=SC2086
    "${runtime}" exec --cleanenv ${binds} "${sif_path}" "$@"
}

# ── run_tool_shell(step_label, sif_path, uri, command) ─────────────────────────
# Runs compound commands (such as the mapping pipeline) inside one SIF image.
run_tool_shell() {
    local step_label="$1"
    local sif_path="$2"
    local uri="$3"
    local command="$4"

    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        run_cmd "${command}"
        return
    fi

    ensure_sif_image "${sif_path}" "${uri}" || return 1
    local runtime="apptainer"
    command -v apptainer &>/dev/null || runtime="singularity"
    local binds
    binds=$(build_bind_mounts)

    log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} ${runtime} exec --cleanenv ${binds} ${sif_path} bash -o pipefail -c ${command}"
        return 0
    fi

    # shellcheck disable=SC2086
    "${runtime}" exec --cleanenv ${binds} "${sif_path}" bash -o pipefail -c "${command}"
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

    # Determine available container runtime
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

# ── build_bind_mounts() ───────────────────────────────────────────────────────
# Returns a string of "-B /abs/path:/abs/path" arguments covering all
# directories that the containerised steps need to read or write.
# Uses the same absolute paths on host and inside the container so that
# any generated scripts that reference ${OUTDIR} etc. work without translation.
run_in_sif() {
    local sif_file="$1"
    shift
    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} [$(basename "${sif_file}")] $*"
        return 0
    fi
    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        local apptainer_bin
        apptainer_bin=$(command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null)
        local binds
        binds=$(build_bind_mounts)
        # Execute via bash -c so compound pipes (like bowtie2 | samtools) run entirely inside the container
        "${apptainer_bin}" exec --cleanenv ${binds} "${sif_file}" bash -c "$*"
    else
        eval "$@"
    fi
}


build_bind_mounts() {
    local -a bind_dirs=()

    _add_bind() {
        local raw_path="$1"
        [[ -z "${raw_path}" ]] && return
        # realpath -m resolves even if the path does not yet exist
        local abs
        abs=$(realpath -m "${raw_path}" 2>/dev/null || echo "${raw_path}")
        # Avoid duplicates
        local d
        for d in "${bind_dirs[@]+"${bind_dirs[@]}"}"; do
            [[ "${d}" == "${abs}" ]] && return
        done
        bind_dirs+=("${abs}")
    }

    _add_bind "${READS_DIR}"
    _add_bind "${OUTDIR}"
    _add_bind "$(dirname "${CHECKM2_DB}")"
    [[ "${RUN_KRAKEN2}" == true ]] && [[ -d "${KRAKEN2_DB}" ]] && _add_bind "${KRAKEN2_DB}"
    [[ "${RUN_GTDBTK}" == true ]] && [[ -d "${GTDBTK_DATA_PATH}" ]] && _add_bind "${GTDBTK_DATA_PATH}"
    [[ -d "${METABINNER_SCRIPTS_DIR}" ]] && _add_bind "${METABINNER_SCRIPTS_DIR}"
    [[ -n "${ASSEMBLY_IN}" ]] && _add_bind "$(dirname "${ASSEMBLY_IN}")"

    local bind_str=""
    local d
    for d in "${bind_dirs[@]+"${bind_dirs[@]}"}"; do
        bind_str+=" -B ${d}:${d}"
    done

    printf '%s' "${bind_str}"
}

# ── _env_to_sif(env_name) → sif_path / uri ───────────────────────────────────
# Maps a mamba environment name to its SIF_* variable pair.
_env_to_sif() {
    local env_name="$1"
    case "${env_name}" in
        "${ENV_VAMB}")       echo "${SIF_VAMB}"       "${SIF_VAMB_URI}" ;;
        "${ENV_METABINNER}") echo "${SIF_METABINNER}"  "${SIF_METABINNER_URI}" ;;
        "${ENV_COMEBIN}")    echo "${SIF_COMEBIN}"     "${SIF_COMEBIN_URI}" ;;
        "${ENV_METAWRAP}")   echo "${SIF_METAWRAP}"    "${SIF_METAWRAP_URI}" ;;
        "${ENV_CHECKM2}")    echo "${SIF_CHECKM2}"     "${SIF_CHECKM2_URI}" ;;
        "${ENV_GTDBTK}")     echo "${SIF_GTDBTK}"      "${SIF_GTDBTK_URI}" ;;
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

# ── run_mamba_file(step_label, env_name, script_file) ────────────────────────
# Hybrid dispatcher: routes execution through Apptainer or mamba depending on
# EXEC_MODE. Signature is unchanged so all step functions work without edits.
run_mamba_file() {
    local step_label="$1"
    local env_name="$2"
    local script_file="$3"

    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        # Resolve container image for this environment
        local sif_path uri
        read -r sif_path uri < <(_env_to_sif "${env_name}")

        log_info "[Step ${step_label}] Engine: Apptainer | Image: ${sif_path}"
        log_info "[Step ${step_label}] Script: ${script_file}"

        # Pull image if absent
        ensure_sif_image "${sif_path}" "${uri}" || return 1

        # Resolve bind-mount string
        local bind_args
        bind_args=$(build_bind_mounts)

        # Optional GPU flag (Vamb and COMEBin can use CUDA)
        local gpu_flag=""
        if { [[ "${env_name}" == "${ENV_VAMB}" ]] || \
             [[ "${env_name}" == "${ENV_COMEBIN}" ]]; } && _has_gpu; then
            gpu_flag="--nv"
            log_info "[Step ${step_label}] GPU detected — adding --nv flag."
        fi

        # Determine available runtime
        local runtime="apptainer"
        command -v apptainer &>/dev/null || runtime="singularity"

        if [[ "${DRY_RUN}" == true ]]; then
            echo -e "${YELLOW}[DRY-RUN]${NC} ${runtime} exec --cleanenv ${gpu_flag}${bind_args} ${sif_path} bash ${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            return 0
        fi

        # shellcheck disable=SC2086
        "${runtime}" exec --cleanenv ${gpu_flag} ${bind_args} \
            "${sif_path}" bash "${script_file}"

    else
        # ── mamba fallback ────────────────────────────────────────────────────
        log_info "[Step ${step_label}] Engine: mamba | Env: ${BOLD}${env_name}${NC}"
        log_info "[Step ${step_label}] Script: ${script_file}"
        if [[ "${DRY_RUN}" == true ]]; then
            echo -e "${YELLOW}[DRY-RUN]${NC} mamba run -n ${env_name} bash ${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
            cat "${script_file}"
            echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
            return 0
        fi
        mamba run -n "${env_name}" bash "${script_file}"
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

is_done() {
    local step="$1"
    [[ "${FORCE_RERUN}" == false ]] && [[ -f "${CHECKPOINT_DIR}/${step}.done" ]]
}

mark_done() {
    local step="$1"

    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[Step ${step}] Dry-run: checkpoint not written."
        return 0
    fi

    touch "${CHECKPOINT_DIR}/${step}.done"
    log_info "[Step ${step}] Checkpoint written."
}

# ── Directory Layout ──────────────────────────────────────────────────────────
QC_DIR="${OUTDIR}/01_qc"
ASSEMBLY_DIR="${OUTDIR}/02_assembly"
KRAKEN2_DIR="${OUTDIR}/03_kraken2"
MAPPING_DIR="${OUTDIR}/04_mapping"
BINNING_DIR="${OUTDIR}/05_binning"
VAMB_DIR="${BINNING_DIR}/vamb"
METABINNER_DIR="${BINNING_DIR}/metabinner"
COMEBIN_DIR="${BINNING_DIR}/comebin"
REFINEMENT_DIR="${OUTDIR}/06_refinement"
REASSEMBLY_DIR="${OUTDIR}/07_reassembly"
CHECKM2_DIR="${OUTDIR}/08_checkm2"
MAGS_DIR="${OUTDIR}/09_mags"
GTDBTK_DIR="${OUTDIR}/10_gtdbtk"
ABUNDANCE_DIR="${OUTDIR}/11_abundance"
SUMMARY_DIR="${OUTDIR}/12_summary"
LOGS_DIR="${OUTDIR}/logs"
CHECKPOINT_DIR="${OUTDIR}/.checkpoints"
TMPSCRIPTS_DIR="${OUTDIR}/logs/mamba_scripts"

# ── Global State Variables ────────────────────────────────────────────────────
CONTIGS=""
CLEAN_R1_FILES=()
CLEAN_R2_FILES=()
SORTED_BAMS=()
DEPTH_FILE=""
BINS_FOR_REFINEMENT_A=""
BINS_FOR_REFINEMENT_B=""
BINS_FOR_REFINEMENT_C=""
REFINED_BINS_DIR=""
BINS_FOR_QC=""
HQ_BINS_DIR=""
# Resource detection state (populated by detect_resources)
_PHYS_CORES=0
_PHYS_MEM_GB=0

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
    local errors=0

    # Reads directory must exist
    if [[ ! -d "${READS_DIR}" ]]; then
        log_error "READS_DIR does not exist: ${READS_DIR}"
        (( errors++ )) || true
    fi

    # Assembly logic
    if [[ "${RUN_MEGAHIT}" == false ]] && [[ -z "${ASSEMBLY_IN}" ]]; then
        log_error "RUN_MEGAHIT=false but ASSEMBLY_IN is not set."
        (( errors++ )) || true
    fi
    if [[ -n "${ASSEMBLY_IN}" ]] && [[ ! -f "${ASSEMBLY_IN}" ]]; then
        log_error "ASSEMBLY_IN file not found: ${ASSEMBLY_IN}"
        (( errors++ )) || true
    fi

    # Kraken2 DB
    if [[ "${RUN_KRAKEN2}" == true ]]; then
        if [[ ! -d "${KRAKEN2_DB}" ]]; then
            log_error "KRAKEN2_DB not set or not found: ${KRAKEN2_DB}"
            (( errors++ )) || true
        fi
    fi

    # CheckM2 DB
    if [[ ! -f "${CHECKM2_DB}" ]]; then
        log_error "CHECKM2_DB not set or not found: ${CHECKM2_DB}"
        (( errors++ )) || true
    fi

    # GTDB-Tk data
    if [[ "${RUN_GTDBTK}" == true ]]; then
        if [[ ! -d "${GTDBTK_DATA_PATH}" ]]; then
            log_error "GTDBTK_DATA_PATH not set or not found: ${GTDBTK_DATA_PATH}"
            (( errors++ )) || true
        fi
    fi

    # MetaBinner scripts directory — validation differs by engine
    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        # mamba: directory must exist on the host
        if [[ -z "${METABINNER_SCRIPTS_DIR}" || ! -d "${METABINNER_SCRIPTS_DIR}" ]]; then
            log_error "METABINNER_SCRIPTS_DIR not set or not found on host: ${METABINNER_SCRIPTS_DIR}"
            (( errors++ )) || true
        fi
    else
        # apptainer: directory lives inside the container
        if [[ -z "${METABINNER_SCRIPTS_DIR}" ]]; then
            log_error "METABINNER_SCRIPTS_DIR is empty."
            log_error "For EXEC_MODE=apptainer use an in-container path (default: /opt/metabinner/scripts)."
            (( errors++ )) || true
        elif [[ -f "${SIF_METABINNER}" ]] && [[ "${DRY_RUN}" == false ]]; then
            # .sif already cached → verify the path exists inside the container
            local _rt="apptainer"
            command -v apptainer &>/dev/null || _rt="singularity"
            log_info "[Preflight] Verifying METABINNER_SCRIPTS_DIR inside container image..."
            if ! "${_rt}" exec --cleanenv "${SIF_METABINNER}" \
                    test -d "${METABINNER_SCRIPTS_DIR}" 2>/dev/null; then
                log_error "METABINNER_SCRIPTS_DIR '${METABINNER_SCRIPTS_DIR}' not found inside ${SIF_METABINNER}."
                log_error "Check the container's internal layout or adjust METABINNER_SCRIPTS_DIR."
                (( errors++ )) || true
            else
                log_info "[Preflight] METABINNER_SCRIPTS_DIR OK inside container: ${METABINNER_SCRIPTS_DIR}"
            fi
        elif [[ ! -f "${SIF_METABINNER}" ]] && [[ "${AUTO_PULL_SIF}" == true ]]; then
            log_info "[Preflight] ${SIF_METABINNER} not yet downloaded; METABINNER_SCRIPTS_DIR will be validated after pull."
        fi
    fi

    # Required direct binaries are host dependencies only in mamba mode.
    if [[ "${EXEC_MODE}" == "mamba" ]]; then
        local required_bins=(bowtie2 bowtie2-build samtools jgi_summarize_bam_contig_depths coverm)
        [[ "${RUN_FASTP}" == true ]]   && required_bins+=(fastp)
        [[ "${RUN_MEGAHIT}" == true ]] && required_bins+=(megahit)
        [[ "${RUN_KRAKEN2}" == true ]] && required_bins+=(kraken2 bracken)
        for tool in "${required_bins[@]}"; do
            if ! command -v "${tool}" &>/dev/null; then
                log_error "Required binary not found in PATH: ${tool}"
                (( errors++ )) || true
            fi
        done
    fi

    # Container / environment runtime
    if [[ "${EXEC_MODE}" == "apptainer" ]]; then
        if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
            log_error "EXEC_MODE=apptainer but neither 'apptainer' nor 'singularity' found in PATH."
            (( errors++ )) || true
        else
            local _rt="apptainer"
            command -v apptainer &>/dev/null || _rt="singularity"
            log_info "Container runtime detected: ${_rt}"
        fi
        # Warn about SIF_DIR accessibility
        if [[ ! -d "${SIF_DIR}" ]]; then
            log_warn "SIF_DIR does not exist yet (${SIF_DIR}); it will be created on first pull."
        fi
    else
        # mamba mode
        if ! command -v mamba &>/dev/null; then
            log_error "'mamba' not found in PATH."
            (( errors++ )) || true
        fi
    fi

    # Detect reads and verify each R1 has its matching R2.
    if [[ -d "${READS_DIR}" ]]; then
        shopt -s nullglob
        local _probe=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        if [[ ${#_probe[@]} -eq 0 ]]; then
            log_error "No reads matching *${READ_SUFFIX_R1} found in ${READS_DIR}"
            (( errors++ )) || true
        else
            local _r1 _r2
            for _r1 in "${_probe[@]}"; do
                _r2="${_r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}"
                if [[ ! -f "${_r2}" ]]; then
                    log_error "R2 file not found for R1 read: ${_r1} (expected ${_r2})"
                    (( errors++ )) || true
                fi
            done
            log_info "Found ${#_probe[@]} R1 read file(s) in ${READS_DIR}"
        fi
    fi

    if [[ "${errors}" -gt 0 ]]; then
        log_error "Pre-flight failed with ${errors} error(s). Fix above and retry."
        exit 1
    fi
    log_info "Pre-flight OK."
}

###############################################################################
#  STEP FUNCTIONS                                                              #
###############################################################################

# ── Step 01: fastp QC ─────────────────────────────────────────────────────────
step_01_fastp() {
    log_step "01" "QC — fastp"

    if [[ "${RUN_FASTP}" == false ]]; then
        log_info "[Step 01] RUN_FASTP=false — skipping fastp."
        return 0
    fi

    if is_done "01_fastp"; then
        log_info "[Step 01] Checkpoint found — repopulating CLEAN arrays from existing trimmed files."
        shopt -s nullglob
        CLEAN_R1_FILES=("${QC_DIR}/trimmed"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        CLEAN_R2_FILES=()
        for r1 in "${CLEAN_R1_FILES[@]}"; do
            CLEAN_R2_FILES+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "[Step 01] Repopulated ${#CLEAN_R1_FILES[@]} sample(s) from checkpoint."
        return 0
    fi

    mkdir -p "${QC_DIR}/trimmed" "${QC_DIR}/reports/html" "${QC_DIR}/reports/json"

    local fastp_threads=$(( THREADS > 16 ? 16 : THREADS ))

    shopt -s nullglob
    local raw_r1_files=("${READS_DIR}"/*${READ_SUFFIX_R1})
    shopt -u nullglob

    if [[ ${#raw_r1_files[@]} -eq 0 ]]; then
        log_error "[Step 01] No R1 read files found in ${READS_DIR} with suffix ${READ_SUFFIX_R1}"
        exit 1
    fi

    CLEAN_R1_FILES=()
    CLEAN_R2_FILES=()

    for r1 in "${raw_r1_files[@]}"; do
        local base
        base=$(basename "${r1}" "${READ_SUFFIX_R1}")
        local r2="${READS_DIR}/${base}${READ_SUFFIX_R2}"

        if [[ ! -f "${r2}" ]]; then
            log_error "[Step 01] R2 file not found for sample ${base}: ${r2}"
            exit 1
        fi

        local out_r1="${QC_DIR}/trimmed/${base}${READ_SUFFIX_R1}"
        local out_r2="${QC_DIR}/trimmed/${base}${READ_SUFFIX_R2}"
        local html_report="${QC_DIR}/reports/html/${base}.html"
        local json_report="${QC_DIR}/reports/json/${base}.json"
        local step_log="${LOGS_DIR}/step_01_fastp_${base}.log"

        log_info "[Step 01] Processing sample: ${base}"
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
            -w "${fastp_threads}" \
            2>&1 | tee "${step_log}"

        CLEAN_R1_FILES+=("${out_r1}")
        CLEAN_R2_FILES+=("${out_r2}")
    done

    mark_done "01_fastp"
}

# ── Step 02: MEGAHIT Assembly ─────────────────────────────────────────────────
step_02_megahit() {
    log_step "02" "Assembly — MEGAHIT"

    # Case 1: User supplied a pre-assembled FASTA
    if [[ -n "${ASSEMBLY_IN}" ]]; then
        if is_done "02_megahit"; then
            log_info "[Step 02] Checkpoint found — using pre-filtered assembly."
            CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
            return 0
        fi
        log_info "[Step 02] ASSEMBLY_IN is set — filtering and cleaning headers from: ${ASSEMBLY_IN}"
        local raw_contigs="${ASSEMBLY_IN}"
        CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
        run_cmd awk \
            -v minlen="${MIN_CONTIG_LEN}" \
            "'BEGIN { header=\"\"; seq=\"\" }
            /^\>/ {
                if (header != \"\" && length(seq) >= minlen) { print header; print seq }
                header = \">\" \$1; sub(/^\>\>/, \">\", header); seq = \"\"; next
            }
            { seq = seq \$0 }
            END { if (header != \"\" && length(seq) >= minlen) { print header; print seq } }'" \
            "${raw_contigs}" ">" "${CONTIGS}"
        mark_done "02_megahit"
        return 0
    fi

    # Case 2: RUN_MEGAHIT=false (skip assembly, use existing)
    if [[ "${RUN_MEGAHIT}" == false ]]; then
        CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
        if [[ ! -f "${CONTIGS}" ]]; then
            log_error "[Step 02] RUN_MEGAHIT=false but assembly not found at: ${CONTIGS}"
            exit 1
        fi
        log_info "[Step 02] RUN_MEGAHIT=false — using existing assembly: ${CONTIGS}"
        return 0
    fi

    # Case 3: Run MEGAHIT
    if is_done "02_megahit"; then
        log_info "[Step 02] Checkpoint found — skipping MEGAHIT."
        CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"
        return 0
    fi

    # Determine input reads
    local r1_list r2_list
    if [[ ${#CLEAN_R1_FILES[@]} -gt 0 ]]; then
        r1_list=$(IFS=','; echo "${CLEAN_R1_FILES[*]}")
        r2_list=$(IFS=','; echo "${CLEAN_R2_FILES[*]}")
        log_info "[Step 02] Using fastp-trimmed reads for assembly."
    else
        shopt -s nullglob
        local raw_r1_arr=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        if [[ ${#raw_r1_arr[@]} -eq 0 ]]; then
            log_error "[Step 02] No raw R1 reads found for assembly in ${READS_DIR}"
            exit 1
        fi
        r1_list=$(IFS=','; echo "${raw_r1_arr[*]}")
        local raw_r2_arr=()
        for r1 in "${raw_r1_arr[@]}"; do
            raw_r2_arr+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        r2_list=$(IFS=','; echo "${raw_r2_arr[*]}")
        log_info "[Step 02] Using raw reads for assembly (fastp was skipped)."
    fi

    local megahit_outdir="${ASSEMBLY_DIR}/megahit_raw"
    local mem_bytes=$(( MEMORY * 1073741824 ))
    local step_log="${LOGS_DIR}/step_02_megahit.log"

    # MEGAHIT will fail if output dir exists
    if [[ -d "${megahit_outdir}" ]]; then
        log_warn "[Step 02] Removing previous MEGAHIT output dir: ${megahit_outdir}"
        run_cmd rm -rf "${megahit_outdir}"
    fi

    run_tool "02" "${SIF_MEGAHIT}" "${SIF_MEGAHIT_URI}" megahit \
        -1 "${r1_list}" \
        -2 "${r2_list}" \
        -o "${megahit_outdir}" \
        -t "${THREADS}" \
        --min-contig-len "${MIN_CONTIG_LEN}" \
        -m "${mem_bytes}" \
        2>&1 | tee "${step_log}"

    local raw_contigs="${megahit_outdir}/final.contigs.fa"
    CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"

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
    log_step "03" "Taxonomy — Kraken2 + Bracken"

    if [[ "${RUN_KRAKEN2}" == false ]]; then
        log_info "[Step 03] RUN_KRAKEN2=false — skipping."
        return 0
    fi

    if is_done "03_kraken2"; then
        log_info "[Step 03] Checkpoint found — skipping Kraken2/Bracken."
        return 0
    fi

    # Determine read sources
    local r1_sources=()
    local r2_sources=()
    if [[ ${#CLEAN_R1_FILES[@]} -gt 0 ]]; then
        r1_sources=("${CLEAN_R1_FILES[@]}")
        r2_sources=("${CLEAN_R2_FILES[@]}")
        log_info "[Step 03] Using fastp-trimmed reads for Kraken2."
    else
        shopt -s nullglob
        r1_sources=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        for r1 in "${r1_sources[@]}"; do
            r2_sources+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "[Step 03] Using raw reads for Kraken2 (fastp was skipped)."
    fi

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
        run_tool "03" "${SIF_KRAKEN2}" "${SIF_KRAKEN2_URI}" kraken2 \
            --db "${KRAKEN2_DB}" \
            --threads "${THREADS}" \
            --paired "${r1}" "${r2}" \
            --output "${kraken_out}" \
            --report "${kraken_report}" \
            2>&1 | tee "${kraken_log}"

        log_info "[Step 03] Running Bracken for sample: ${base}"
        run_tool "03" "${SIF_KRAKEN2}" "${SIF_KRAKEN2_URI}" bracken \
            -d "${KRAKEN2_DB}" \
            -i "${kraken_report}" \
            -o "${bracken_out}" \
            -r "${BRACKEN_READ_LEN}" \
            -l "${BRACKEN_TAX_LEVEL}" \
            -t "${BRACKEN_THRESHOLD}" \
            2>&1 | tee "${bracken_log}"
    done

    mark_done "03_kraken2"
}

# ── Step 04: Bowtie2 Mapping + JGI Depth ─────────────────────────────────────
step_04_mapping_depth() {
    log_step "04" "Mapping — Bowtie2 + SAMtools + JGI depth"

    if is_done "04_mapping"; then
        log_info "[Step 04] Checkpoint found — repopulating SORTED_BAMS and DEPTH_FILE."
        shopt -s nullglob
        SORTED_BAMS=("${MAPPING_DIR}"/*_sorted.bam)
        shopt -u nullglob
        DEPTH_FILE="${BINNING_DIR}/depth_matrix.txt"
        log_info "[Step 04] Repopulated ${#SORTED_BAMS[@]} BAM(s) from checkpoint."
        return 0
    fi

    if [[ -z "${CONTIGS}" ]] || [[ ! -f "${CONTIGS}" ]]; then
        log_error "[Step 04] CONTIGS file not found: ${CONTIGS}"
        exit 1
    fi

    mkdir -p "${MAPPING_DIR}"

    local bt2_index="${MAPPING_DIR}/contigs_idx"
    local build_log="${LOGS_DIR}/step_04_bowtie2_build.log"

    log_info "[Step 04] Building Bowtie2 index..."
    run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" bowtie2-build \
        --threads "${THREADS}" \
        "${CONTIGS}" \
        "${bt2_index}" \
        2>&1 | tee "${build_log}"

    # Determine read sources
    local r1_sources=()
    local r2_sources=()
    if [[ ${#CLEAN_R1_FILES[@]} -gt 0 ]]; then
        r1_sources=("${CLEAN_R1_FILES[@]}")
        r2_sources=("${CLEAN_R2_FILES[@]}")
        log_info "[Step 04] Using fastp-trimmed reads for mapping."
    else
        shopt -s nullglob
        r1_sources=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        for r1 in "${r1_sources[@]}"; do
            r2_sources+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "[Step 04] Using raw reads for mapping (fastp was skipped)."
    fi

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

        log_info "[Step 04] Mapping sample: ${base}"
        run_tool_shell "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" "bowtie2 --very-sensitive-local -p ${THREADS} --time --reorder \
            -x ${bt2_index} \
            -1 ${r1} \
            -2 ${r2} \
            2>>${bt2_log} \
            | samtools view -@ ${THREADS} -bS - \
            | samtools sort -@ ${THREADS} -o ${bam_sorted} -"

        log_info "[Step 04] Indexing BAM: ${bam_sorted}"
        run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" samtools index -@ "${THREADS}" "${bam_sorted}"

        SORTED_BAMS+=("${bam_sorted}")
    done

    DEPTH_FILE="${BINNING_DIR}/depth_matrix.txt"
    local jgi_log="${LOGS_DIR}/step_04_jgi_depth.log"

    log_info "[Step 04] Running jgi_summarize_bam_contig_depths..."
    run_tool "04" "${SIF_MAPPING}" "${SIF_MAPPING_URI}" jgi_summarize_bam_contig_depths \
        --outputDepth "${DEPTH_FILE}" \
        "${SORTED_BAMS[@]}" \
        2>&1 | tee "${jgi_log}"

    mark_done "04_mapping"
}

# ── Step 05a: Vamb Binning ────────────────────────────────────────────────────
step_05a_vamb() {
    log_step "05a" "Binning — Vamb"

    if is_done "05a_vamb"; then
        log_info "[Step 05a] Checkpoint found — skipping Vamb."
        BINS_FOR_REFINEMENT_A="${VAMB_DIR}/vamb_out/bins"
        return 0
    fi

    if [[ -z "${DEPTH_FILE}" ]] || [[ ! -f "${DEPTH_FILE}" ]]; then
        log_error "[Step 05a] Depth file not found: ${DEPTH_FILE}"
        exit 1
    fi

    mkdir -p "${VAMB_DIR}"

    log_info "[Step 05a] Converting JGI depth matrix to Vamb abundance TSV..."
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
    }' "${DEPTH_FILE}" > "${VAMB_DIR}/vamb_abundance.tsv"

    log_info "[Step 05a] Vamb abundance TSV created: ${VAMB_DIR}/vamb_abundance.tsv"

    local vamb_script="${TMPSCRIPTS_DIR}/step_05a_vamb.sh"
    cat > "${vamb_script}" <<VAMB_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
vamb bin default \\
    --outdir "${VAMB_DIR}/vamb_out" \\
    --fasta  "${CONTIGS}" \\
    --abundance_tsv "${VAMB_DIR}/vamb_abundance.tsv" \\
    -o C \\
    --minfasta ${MIN_CONTIG_LEN}
VAMB_SCRIPT

    chmod +x "${vamb_script}"

    run_mamba_file "05a" "${ENV_VAMB}" "${vamb_script}" 2>&1 | tee "${LOGS_DIR}/step_05a_vamb.log"

    # Rename .fna → .fa
    local vamb_bins_dir="${VAMB_DIR}/vamb_out/bins"
    if [[ -d "${vamb_bins_dir}" ]]; then
        local fna_count
        fna_count=$(find "${vamb_bins_dir}" -name '*.fna' | wc -l)
        if [[ "${fna_count}" -gt 0 ]]; then
            log_info "[Step 05a] Renaming ${fna_count} .fna file(s) to .fa in ${vamb_bins_dir}"
            find "${vamb_bins_dir}" -name '*.fna' | while read -r f; do
                mv "${f}" "${f%.fna}.fa"
            done
        fi
    else
        log_warn "[Step 05a] Vamb bins directory not found: ${vamb_bins_dir}"
    fi

    BINS_FOR_REFINEMENT_A="${VAMB_DIR}/vamb_out/bins"
    mark_done "05a_vamb"
}

# ── Step 05b: MetaBinner Binning ──────────────────────────────────────────────
step_05b_metabinner() {
    log_step "05b" "Binning — MetaBinner"

    if is_done "05b_metabinner"; then
        log_info "[Step 05b] Checkpoint found — skipping MetaBinner."
        BINS_FOR_REFINEMENT_C="${METABINNER_DIR}/metabinner_out/result_bins"
        return 0
    fi

    if [[ -z "${DEPTH_FILE}" ]] || [[ ! -f "${DEPTH_FILE}" ]]; then
        log_error "[Step 05b] Depth file not found: ${DEPTH_FILE}"
        exit 1
    fi

    mkdir -p "${METABINNER_DIR}"

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

    # Derive filenames
    local mb_filtered="${METABINNER_DIR}/assembly_${MIN_CONTIG_LEN}.fa"
    local mb_kmer="${METABINNER_DIR}/assembly_${MIN_CONTIG_LEN}_kmer_4_f${MIN_CONTIG_LEN}.csv"

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
    run_mamba_file "05b-prep" "${ENV_METABINNER}" "${mb_prep_script}" 2>&1 | tee "${LOGS_DIR}/step_05b_metabinner_prep.log"

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
    run_mamba_file "05b-run" "${ENV_METABINNER}" "${mb_run_script}" 2>&1 | tee "${LOGS_DIR}/step_05b_metabinner_run.log"

    # Rename .fasta → .fa
    local mb_bins_dir="${METABINNER_DIR}/metabinner_out/result_bins"
    if [[ -d "${mb_bins_dir}" ]]; then
        local fasta_count
        fasta_count=$(find "${mb_bins_dir}" -name '*.fasta' | wc -l)
        if [[ "${fasta_count}" -gt 0 ]]; then
            log_info "[Step 05b] Renaming ${fasta_count} .fasta file(s) to .fa in ${mb_bins_dir}"
            find "${mb_bins_dir}" -name '*.fasta' | while read -r f; do
                mv "${f}" "${f%.fasta}.fa"
            done
        fi
    else
        log_warn "[Step 05b] MetaBinner bins directory not found: ${mb_bins_dir}"
    fi

    BINS_FOR_REFINEMENT_C="${METABINNER_DIR}/metabinner_out/result_bins"
    mark_done "05b_metabinner"
}

# ── Step 05c: COMEBin Binning ─────────────────────────────────────────────────
step_05c_comebin() {
    log_step "05c" "Binning — COMEBin"

    if is_done "05c_comebin"; then
        log_info "[Step 05c] Checkpoint found — skipping COMEBin."
        BINS_FOR_REFINEMENT_B="${COMEBIN_DIR}/comebin_res/comebin_res_bins"
        return 0
    fi

    mkdir -p "${COMEBIN_DIR}"

    local comebin_script="${TMPSCRIPTS_DIR}/step_05c_comebin.sh"
    cat > "${comebin_script}" <<COMEBIN_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
bash "${COMEBIN_SCRIPT}" \\
    -a "${CONTIGS}" \\
    -o "${COMEBIN_DIR}" \\
    -p "${MAPPING_DIR}" \\
    -s 42 \\
    -t "${THREADS}"
COMEBIN_SCRIPT

    chmod +x "${comebin_script}"
    run_mamba_file "05c" "${ENV_COMEBIN}" "${comebin_script}" 2>&1 | tee "${LOGS_DIR}/step_05c_comebin.log"

    # Rename .fasta → .fa
    local cb_bins_dir="${COMEBIN_DIR}/comebin_res/comebin_res_bins"
    if [[ -d "${cb_bins_dir}" ]]; then
        local fasta_count
        fasta_count=$(find "${cb_bins_dir}" -name '*.fasta' | wc -l)
        if [[ "${fasta_count}" -gt 0 ]]; then
            log_info "[Step 05c] Renaming ${fasta_count} .fasta file(s) to .fa in ${cb_bins_dir}"
            find "${cb_bins_dir}" -name '*.fasta' | while read -r f; do
                mv "${f}" "${f%.fasta}.fa"
            done
        fi
    else
        log_warn "[Step 05c] COMEBin bins directory not found: ${cb_bins_dir}"
    fi

    BINS_FOR_REFINEMENT_B="${COMEBIN_DIR}/comebin_res/comebin_res_bins"
    mark_done "05c_comebin"
}

# ── Step 06: MetaWrap Bin Refinement ──────────────────────────────────────────
step_06_bin_refinement() {
    log_step "06" "Refinement — MetaWrap bin_refinement"

    if is_done "06_refinement"; then
        log_info "[Step 06] Checkpoint found — skipping bin refinement."
        REFINED_BINS_DIR="${REFINEMENT_DIR}/metawrap_${METAWRAP_COMPLETENESS}_${METAWRAP_CONTAMINATION}_bins"
        BINS_FOR_QC="${REFINED_BINS_DIR}"
        return 0
    fi

    if [[ -z "${BINS_FOR_REFINEMENT_A}" ]]; then
        log_error "[Step 06] BINS_FOR_REFINEMENT_A is not set (Vamb bins)."
        exit 1
    fi
    if [[ -z "${BINS_FOR_REFINEMENT_B}" ]]; then
        log_error "[Step 06] BINS_FOR_REFINEMENT_B is not set (COMEBin bins)."
        exit 1
    fi
    if [[ -z "${BINS_FOR_REFINEMENT_C}" ]]; then
        log_error "[Step 06] BINS_FOR_REFINEMENT_C is not set (MetaBinner bins)."
        exit 1
    fi

    mkdir -p "${REFINEMENT_DIR}"

    local refinement_script="${TMPSCRIPTS_DIR}/step_06_bin_refinement.sh"
    cat > "${refinement_script}" <<REFINEMENT_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
metawrap bin_refinement \\
    -o "${REFINEMENT_DIR}" \\
    -t "${THREADS}" \\
    -c "${METAWRAP_COMPLETENESS}" \\
    -x "${METAWRAP_CONTAMINATION}" \\
    -A "${BINS_FOR_REFINEMENT_A}" \\
    -B "${BINS_FOR_REFINEMENT_B}" \\
    -C "${BINS_FOR_REFINEMENT_C}" \\
    --quick
REFINEMENT_SCRIPT

    chmod +x "${refinement_script}"
    run_mamba_file "06" "${ENV_METAWRAP}" "${refinement_script}" 2>&1 | tee "${LOGS_DIR}/step_06_bin_refinement.log"

    REFINED_BINS_DIR="${REFINEMENT_DIR}/metawrap_${METAWRAP_COMPLETENESS}_${METAWRAP_CONTAMINATION}_bins"
    BINS_FOR_QC="${REFINED_BINS_DIR}"
    mark_done "06_refinement"
}

# ── Step 07: MetaWrap Reassembly ──────────────────────────────────────────────
step_07_reassembly() {
    log_step "07" "Reassembly — MetaWrap reassemble_bins"

    if [[ "${RUN_REASSEMBLY}" == false ]]; then
        log_info "[Step 07] RUN_REASSEMBLY=false — skipping reassembly."
        BINS_FOR_QC="${REFINED_BINS_DIR}"
        return 0
    fi

    if is_done "07_reassembly"; then
        log_info "[Step 07] Checkpoint found — skipping reassembly."
        BINS_FOR_QC="${REASSEMBLY_DIR}/reassembled_bins"
        return 0
    fi

    mkdir -p "${REASSEMBLY_DIR}"

    # Determine read sources
    local r1_sources=()
    local r2_sources=()
    if [[ ${#CLEAN_R1_FILES[@]} -gt 0 ]]; then
        r1_sources=("${CLEAN_R1_FILES[@]}")
        r2_sources=("${CLEAN_R2_FILES[@]}")
        log_info "[Step 07] Using fastp-trimmed reads for reassembly."
    else
        shopt -s nullglob
        r1_sources=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        for r1 in "${r1_sources[@]}"; do
            r2_sources+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "[Step 07] Using raw reads for reassembly (fastp was skipped)."
    fi

    local all_r1="${REASSEMBLY_DIR}/all_reads_R1.fastq.gz"
    local all_r2="${REASSEMBLY_DIR}/all_reads_R2.fastq.gz"

    log_info "[Step 07] Concatenating R1 reads into ${all_r1}..."
    run_cmd cat "${r1_sources[@]}" ">" "${all_r1}"

    log_info "[Step 07] Concatenating R2 reads into ${all_r2}..."
    run_cmd cat "${r2_sources[@]}" ">" "${all_r2}"

    local reassembly_script="${TMPSCRIPTS_DIR}/step_07_reassembly.sh"
    cat > "${reassembly_script}" <<REASSEMBLY_SCRIPT
#!/usr/bin/env bash
set -euo pipefail
metawrap reassemble_bins \\
    -o "${REASSEMBLY_DIR}" \\
    -b "${REFINED_BINS_DIR}" \\
    -1 "${REASSEMBLY_DIR}/all_reads_R1.fastq.gz" \\
    -2 "${REASSEMBLY_DIR}/all_reads_R2.fastq.gz" \\
    -t "${THREADS}" \\
    -m "${MEMORY}"
REASSEMBLY_SCRIPT

    chmod +x "${reassembly_script}"
    run_mamba_file "07" "${ENV_METAWRAP}" "${reassembly_script}" 2>&1 | tee "${LOGS_DIR}/step_07_reassembly.log"

    BINS_FOR_QC="${REASSEMBLY_DIR}/reassembled_bins"
    mark_done "07_reassembly"
}

# ── Step 08: CheckM2 Quality Assessment ───────────────────────────────────────
step_08_checkm2() {
    log_step "08" "Quality — CheckM2 predict"

    if is_done "08_checkm2"; then
        log_info "[Step 08] Checkpoint found — skipping CheckM2."
        return 0
    fi

    if [[ -z "${BINS_FOR_QC}" ]] || [[ ! -d "${BINS_FOR_QC}" ]]; then
        log_error "[Step 08] BINS_FOR_QC directory not found: ${BINS_FOR_QC}"
        exit 1
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
    run_mamba_file "08" "${ENV_CHECKM2}" "${checkm2_script}" 2>&1 | tee "${LOGS_DIR}/step_08_checkm2.log"

    if [[ ! -f "${CHECKM2_DIR}/quality_report.tsv" ]]; then
        log_error "[Step 08] CheckM2 quality_report.tsv not found. CheckM2 may have failed."
        exit 1
    fi

    log_info "[Step 08] CheckM2 quality report: ${CHECKM2_DIR}/quality_report.tsv"
    mark_done "08_checkm2"
}

# ── Step 09: Filter MAGs ──────────────────────────────────────────────────────
step_09_filter_mags() {
    log_step "09" "Filter — HQ/MQ MAG selection (Completeness >= ${COMPLETENESS}%, Contamination <= ${CONTAMINATION}%)"

    if is_done "09_filter_mags"; then
        log_info "[Step 09] Checkpoint found — skipping MAG filtering."
        HQ_BINS_DIR="${MAGS_DIR}"
        return 0
    fi

    if [[ ! -f "${CHECKM2_DIR}/quality_report.tsv" ]]; then
        log_error "[Step 09] CheckM2 quality_report.tsv not found: ${CHECKM2_DIR}/quality_report.tsv"
        exit 1
    fi

    mkdir -p "${MAGS_DIR}"

    local passing_bins=()

    # Read quality_report.tsv: skip header, apply float comparison
    while IFS=$'\t' read -r bin_name completeness contamination rest; do
        # Skip header line
        [[ "${bin_name}" == "Name" ]] && continue
        [[ -z "${bin_name}" ]] && continue

        # Use awk for float comparison
        local passes
        passes=$(awk -v comp="${completeness}" -v cont="${contamination}" \
                     -v min_comp="${COMPLETENESS}" -v max_cont="${CONTAMINATION}" \
                     'BEGIN { if (comp+0 >= min_comp+0 && cont+0 <= max_cont+0) print "yes"; else print "no" }')

        if [[ "${passes}" == "yes" ]]; then
            local bin_file="${BINS_FOR_QC}/${bin_name}.fa"
            if [[ -f "${bin_file}" ]]; then
                passing_bins+=("${bin_name}")
                cp "${bin_file}" "${MAGS_DIR}/${bin_name}.fa"
                log_info "[Step 09] Passing MAG: ${bin_name} (Completeness=${completeness}%, Contamination=${contamination}%)"
            else
                log_warn "[Step 09] Bin file not found, skipping: ${bin_file}"
            fi
        fi
    done < "${CHECKM2_DIR}/quality_report.tsv"

    if [[ ${#passing_bins[@]} -eq 0 ]]; then
        log_warn "[Step 09] WARNING: Zero bins passed the quality filters (Completeness >= ${COMPLETENESS}%, Contamination <= ${CONTAMINATION}%)."
        log_warn "[Step 09] Downstream steps (GTDB-Tk, CoverM, Summary) may produce empty results."
    else
        log_info "[Step 09] ${#passing_bins[@]} MAG(s) passed quality filters. Written to: ${MAGS_DIR}"
    fi

    HQ_BINS_DIR="${MAGS_DIR}"
    mark_done "09_filter_mags"
}

# ── Step 10: GTDB-Tk Taxonomy ─────────────────────────────────────────────────
step_10_gtdbtk() {
    log_step "10" "Taxonomy — GTDB-Tk classify_wf"

    if [[ "${RUN_GTDBTK}" == false ]]; then
        log_info "[Step 10] RUN_GTDBTK=false — skipping GTDB-Tk."
        return 0
    fi

    if is_done "10_gtdbtk"; then
        log_info "[Step 10] Checkpoint found — skipping GTDB-Tk."
        return 0
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
    run_mamba_file "10" "${ENV_GTDBTK}" "${gtdbtk_script}" 2>&1 | tee "${LOGS_DIR}/step_10_gtdbtk.log"

    mark_done "10_gtdbtk"
}

# ── Step 11: CoverM Abundance ─────────────────────────────────────────────────
step_11_coverm() {
    log_step "11" "Abundance — CoverM genome"

    if is_done "11_coverm"; then
        log_info "[Step 11] Checkpoint found — skipping CoverM."
        return 0
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

    # Determine read sources
    local r1_sources=()
    local r2_sources=()
    if [[ ${#CLEAN_R1_FILES[@]} -gt 0 ]]; then
        r1_sources=("${CLEAN_R1_FILES[@]}")
        r2_sources=("${CLEAN_R2_FILES[@]}")
        log_info "[Step 11] Using fastp-trimmed reads for CoverM."
    else
        shopt -s nullglob
        r1_sources=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        for r1 in "${r1_sources[@]}"; do
            r2_sources+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "[Step 11] Using raw reads for CoverM (fastp was skipped)."
    fi

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
        local methods_flag=""
        local min_cov_flag=""
        if [[ -n "${COVERM_METHODS}" ]]; then
            methods_flag="--methods ${COVERM_METHODS}"
        fi
        if [[ -n "${COVERM_MIN_COVERED_FRACTION}" ]]; then
            min_cov_flag="--min-covered-fraction ${COVERM_MIN_COVERED_FRACTION}"
        fi

        log_info "[Step 11] Running CoverM for sample: ${base}"
        run_tool "11" "${SIF_COVERM}" "${SIF_COVERM_URI}" coverm genome \
            --read1 "${r1}" \
            --read2 "${r2}" \
            --genome-fasta-directory "${HQ_BINS_DIR}" \
            --genome-fasta-extension fa \
            --threads "${THREADS}" \
            ${methods_flag} \
            ${min_cov_flag} \
            -o "${out_tsv}" \
            2>&1 | tee "${coverm_log}"
    done

    mark_done "11_coverm"
}

# ── Step 12: R Summary Table ───────────────────────────────────────────────────
step_12_summary_r() {
    log_step "12" "Summary — R left-join TSV (base R)"

    if is_done "12_summary"; then
        log_info "[Step 12] Checkpoint found — skipping summary."
        return 0
    fi

    mkdir -p "${SUMMARY_DIR}"

    local checkm2_tsv="${CHECKM2_DIR}/quality_report.tsv"
    local gtdbtk_bac="${GTDBTK_DIR}/gtdbtk.bac120.summary.tsv"
    local gtdbtk_ar="${GTDBTK_DIR}/gtdbtk.ar53.summary.tsv"
    local abundance_dir="${ABUNDANCE_DIR}"
    local out_tsv="${SUMMARY_DIR}/final_summary.tsv"

    if [[ ! -f "${checkm2_tsv}" ]]; then
        log_error "[Step 12] CheckM2 quality_report.tsv not found: ${checkm2_tsv}"
        exit 1
    fi

    log_info "[Step 12] Running R summary script..."

    run_tool "12" "${SIF_R}" "${SIF_R_URI}" Rscript - <<RSCRIPT
###############################################################################
# Final summary — base R only, no external packages
###############################################################################

checkm2_tsv  <- "${checkm2_tsv}"
gtdbtk_bac   <- "${gtdbtk_bac}"
gtdbtk_ar    <- "${gtdbtk_ar}"
abundance_dir <- "${abundance_dir}"
out_tsv      <- "${out_tsv}"

cat("[Step 12] Reading CheckM2 quality report...\n")
checkm2_df <- read.table(checkm2_tsv, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
# Rename first column to Bin_ID
colnames(checkm2_df)[1] <- "Bin_ID"
cat(sprintf("[Step 12] CheckM2: %d bins\n", nrow(checkm2_df)))

cat("[Step 12] Reading GTDB-Tk outputs...\n")
bac_df <- NULL
ar_df  <- NULL

if (file.exists(gtdbtk_bac)) {
    bac_df <- read.table(gtdbtk_bac, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
    colnames(bac_df)[colnames(bac_df) == "user_genome"] <- "Bin_ID"
    cat(sprintf("[Step 12] GTDB-Tk bac120: %d entries\n", nrow(bac_df)))
} else {
    cat("[Step 12] WARN: gtdbtk.bac120.summary.tsv not found.\n")
}

if (file.exists(gtdbtk_ar)) {
    ar_df <- read.table(gtdbtk_ar, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)
    colnames(ar_df)[colnames(ar_df) == "user_genome"] <- "Bin_ID"
    cat(sprintf("[Step 12] GTDB-Tk ar53: %d entries\n", nrow(ar_df)))
} else {
    cat("[Step 12] WARN: gtdbtk.ar53.summary.tsv not found.\n")
}

# Combine bac120 + ar53 by intersected columns
if (!is.null(bac_df) && !is.null(ar_df)) {
    shared_cols <- intersect(colnames(bac_df), colnames(ar_df))
    gtdbtk_df <- rbind(bac_df[, shared_cols, drop = FALSE],
                        ar_df[, shared_cols, drop = FALSE])
    cat(sprintf("[Step 12] GTDB-Tk combined: %d entries\n", nrow(gtdbtk_df)))
} else if (!is.null(bac_df)) {
    gtdbtk_df <- bac_df
} else if (!is.null(ar_df)) {
    gtdbtk_df <- ar_df
} else {
    gtdbtk_df <- NULL
    cat("[Step 12] WARN: No GTDB-Tk data available — adding NA column.\n")
}

cat("[Step 12] Reading CoverM abundance tables...\n")
abund_files <- list.files(abundance_dir, pattern = "^abundance_.*\\.tsv$",
                           full.names = TRUE)

if (length(abund_files) == 0) {
    cat("[Step 12] WARN: No CoverM abundance files found.\n")
    coverm_wide <- NULL
} else {
    abund_list <- lapply(abund_files, function(f) {
        sample_name <- sub("abundance_", "",
                           sub("\\.tsv$", "", basename(f)))
        df <- read.table(f, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
        colnames(df)[1] <- "Bin_ID"
        # Strip path and .fa from Bin_ID values
        df[["Bin_ID"]] <- sub("\\.fa$", "", basename(df[["Bin_ID"]]))
        # Rename remaining columns as <sample>__<metric>, sanitize spaces
        for (j in seq(2, ncol(df))) {
            orig <- colnames(df)[j]
            sanitized <- gsub(" ", "_", orig)
            colnames(df)[j] <- paste0(sample_name, "__", sanitized)
        }
        df
    })

    # Merge all CoverM tables wide by Bin_ID
    coverm_wide <- Reduce(function(a, b) {
        merge(a, b, by = "Bin_ID", all = TRUE)
    }, abund_list)
    cat(sprintf("[Step 12] CoverM wide table: %d bins x %d columns\n",
                nrow(coverm_wide), ncol(coverm_wide)))
}

cat("[Step 12] Merging all tables...\n")
# Start with CheckM2
summary_df <- checkm2_df

# Left-join GTDB-Tk
if (!is.null(gtdbtk_df)) {
    summary_df <- merge(summary_df, gtdbtk_df, by = "Bin_ID", all.x = TRUE)
} else {
    summary_df[["GTDB_classification"]] <- NA_character_
}

# Left-join CoverM wide
if (!is.null(coverm_wide)) {
    summary_df <- merge(summary_df, coverm_wide, by = "Bin_ID", all.x = TRUE)
}

cat(sprintf("[Step 12] Final summary: %d bins x %d columns\n",
            nrow(summary_df), ncol(summary_df)))

write.table(summary_df, file = out_tsv, sep = "\t",
            row.names = FALSE, quote = FALSE, na = "NA")

cat(sprintf("[Step 12] Summary written to: %s\n", out_tsv))
RSCRIPT

    if [[ ! -f "${out_tsv}" ]]; then
        log_error "[Step 12] Summary TSV not created: ${out_tsv}"
        exit 1
    fi

    mark_done "12_summary"
}

###############################################################################
#  MAIN                                                                        #
###############################################################################

main() {
    # Create all required directories
    mkdir -p "${OUTDIR}" "${LOGS_DIR}" "${CHECKPOINT_DIR}" "${TMPSCRIPTS_DIR}"
    for d in "${QC_DIR}" "${ASSEMBLY_DIR}" "${KRAKEN2_DIR}" "${MAPPING_DIR}" \
              "${BINNING_DIR}" "${VAMB_DIR}" "${METABINNER_DIR}" "${COMEBIN_DIR}" \
              "${REFINEMENT_DIR}" "${REASSEMBLY_DIR}" "${CHECKM2_DIR}" \
              "${MAGS_DIR}" "${GTDBTK_DIR}" "${ABUNDANCE_DIR}" "${SUMMARY_DIR}"; do
        mkdir -p "${d}"
    done

    # Auto-detect / validate compute resources before printing the banner
    detect_resources

    # Print run configuration banner
    echo -e "\n${CYAN}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}${BOLD}║         SHOTGUN METAGENOMICS PIPELINE — CONFIGURATION         ║${NC}"
    echo -e "${CYAN}${BOLD}╠══════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${CYAN}${BOLD}║${NC}  READS_DIR              : ${READS_DIR}"
    echo -e "${CYAN}${BOLD}║${NC}  OUTDIR                 : ${OUTDIR}"
    echo -e "${CYAN}${BOLD}║${NC}  ASSEMBLY_IN            : ${ASSEMBLY_IN:-<run MEGAHIT>}"
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
    echo -e "${CYAN}${BOLD}║${NC}  COMPLETENESS threshold : >= ${COMPLETENESS}%"
    echo -e "${CYAN}${BOLD}║${NC}  CONTAMINATION threshold: <= ${CONTAMINATION}%"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_FASTP              : ${RUN_FASTP}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_MEGAHIT            : ${RUN_MEGAHIT}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_KRAKEN2            : ${RUN_KRAKEN2}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_REASSEMBLY         : ${RUN_REASSEMBLY}"
    echo -e "${CYAN}${BOLD}║${NC}  RUN_GTDBTK             : ${RUN_GTDBTK}"
    echo -e "${CYAN}${BOLD}║${NC}  DRY_RUN                : ${DRY_RUN}"
    echo -e "${CYAN}${BOLD}║${NC}  FORCE_RERUN            : ${FORCE_RERUN}"
    echo -e "${CYAN}${BOLD}╚══════════════════════════════════════════════════════════════╝${NC}\n"

    # Pre-flight check
    preflight_check

    # Initialize read arrays:
    # If RUN_FASTP=false, pre-populate CLEAN_R1_FILES/CLEAN_R2_FILES with raw reads
    if [[ "${RUN_FASTP}" == false ]]; then
        shopt -s nullglob
        CLEAN_R1_FILES=("${READS_DIR}"/*${READ_SUFFIX_R1})
        shopt -u nullglob
        for r1 in "${CLEAN_R1_FILES[@]}"; do
            CLEAN_R2_FILES+=("${r1/${READ_SUFFIX_R1}/${READ_SUFFIX_R2}}")
        done
        log_info "RUN_FASTP=false — pre-populated read arrays with ${#CLEAN_R1_FILES[@]} raw sample(s)."
    fi

    # If ASSEMBLY_IN is set, pre-set CONTIGS path
    [[ -n "${ASSEMBLY_IN}" ]] && CONTIGS="${ASSEMBLY_DIR}/final.contigs.fa"

    # Execute all steps
    step_01_fastp
    step_02_megahit
    step_03_kraken2_bracken
    step_04_mapping_depth
    step_05a_vamb
    step_05b_metabinner
    step_05c_comebin
    step_06_bin_refinement
    step_07_reassembly
    # Ensure BINS_FOR_QC is set (fallback to refined bins)
    BINS_FOR_QC="${BINS_FOR_QC:-${REFINED_BINS_DIR}}"
    step_08_checkm2
    step_09_filter_mags
    step_10_gtdbtk
    step_11_coverm
    step_12_summary_r

    # Print final summary
    echo -e "\n${GREEN}${BOLD}PIPELINE COMPLETE 🎉${NC}"
    echo "Output root    : ${OUTDIR}"
    echo "Assembly       : ${CONTIGS}"
    echo "Mapped samples : ${#SORTED_BAMS[@]}"
    echo "Vamb bins      : $(find "${VAMB_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "MetaBinner bins: $(find "${METABINNER_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "COMEBin bins   : $(find "${COMEBIN_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "Refined bins   : $(find "${REFINEMENT_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "MAGs (HQ/MQ)   : $(find "${MAGS_DIR}" -name '*.fa' 2>/dev/null | wc -l)"
    echo "Final summary  : ${SUMMARY_DIR}/final_summary.tsv"
    echo "Logs           : ${LOGS_DIR}/"
    echo "Mamba scripts  : ${TMPSCRIPTS_DIR}/"
}

main "$@"
