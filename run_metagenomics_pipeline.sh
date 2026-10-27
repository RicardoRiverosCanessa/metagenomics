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
#   Direct binaries : fastp, megahit, kraken2, bracken, bowtie2, bowtie2-build,
#                     samtools, jgi_summarize_bam_contig_depths, coverm
#                     → called directly via $BIN_DIR prepended to $PATH
#   Isolated envs   : vamb, metabinner, comebin, metawrap, checkm2, gtdbtk
#                     → ALWAYS via  mamba run -n <ENV_NAME>  (NEVER conda activate)
###############################################################################

set -euo pipefail
trap 'echo -e "\n[FATAL] Pipeline aborted at line ${LINENO}. Check logs in ${OUTDIR}/logs/" >&2' ERR

###############################################################################
#  CONFIGURATION  —  Edit everything in this block before running             #
###############################################################################

# ── Input / Output ────────────────────────────────────────────────────────────
READS_DIR="/path/to/reads"      # Dir containing paired-end reads
OUTDIR="./metagenomics_out"     # All pipeline output goes here
ASSEMBLY_IN=""                   # Pre-assembled FASTA. Empty = run MEGAHIT.

# ── Read File Patterns ────────────────────────────────────────────────────────
READ_SUFFIX_R1="_R1.fastq.gz"
READ_SUFFIX_R2="_R2.fastq.gz"

# ── Compute Resources ─────────────────────────────────────────────────────────
THREADS=32
MEMORY=120      # GB — MEGAHIT memory cap (converted to bytes internally)

# ── Tool Binaries Directory ───────────────────────────────────────────────────
# Compiled binaries accessed directly. Leave empty ("") if already on $PATH.
BIN_DIR=""

# ── Mamba Environment Names ───────────────────────────────────────────────────
ENV_VAMB="vamb"
ENV_METABINNER="metabinner_env"
ENV_COMEBIN="comebin_env"
ENV_METAWRAP="metawrap_env"
ENV_CHECKM2="checkm2_env"
ENV_GTDBTK="gtdbtk_env"

# ── Database Paths (override by exporting env vars before running) ────────────
KRAKEN2_DB="${KRAKEN2_DB:-/path/to/kraken2_db}"
CHECKM2_DB="${CHECKM2_DB:-/path/to/uniref100.KO.1.dmnd}"
GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-/path/to/gtdbtk_data}"

# ── MetaBinner Script Paths ───────────────────────────────────────────────────
METABINNER_SCRIPTS_DIR="${METABINNER_SCRIPTS_DIR:-/path/to/metabinner/scripts}"
METABINNER_RUNNER="${METABINNER_RUNNER:-run_metabinner.sh}"

# ── COMEBin Script Path ───────────────────────────────────────────────────────
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
RUN_FASTP=true        # Step 01
RUN_MEGAHIT=true      # Step 02 — false requires ASSEMBLY_IN to be set
RUN_KRAKEN2=true      # Step 03
RUN_REASSEMBLY=false  # Step 07
RUN_GTDBTK=true       # Step 10

# ── Execution Flags ───────────────────────────────────────────────────────────
DRY_RUN=false         # Print commands without running them
FORCE_RERUN=false     # Ignore checkpoints and rerun every step

###############################################################################
#  END OF CONFIGURATION                                                       #
###############################################################################

# ── CLI Parsing ───────────────────────────────────────────────────────────────
_show_help() {
    grep '^#' "$0" | head -40 | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)  DRY_RUN=true;     shift ;;
        --force)    FORCE_RERUN=true;  shift ;;
        -h|--help)  _show_help ;;
        *)  echo "[ERROR] Unknown option: $1. Use -h for help." >&2; exit 1 ;;
    esac
done

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
        echo -e "${YELLOW}[DRY-RUN]${NC} $*"
    else
        log_info "Running: $*"
        eval "$@"
    fi
}

run_mamba_file() {
    local step_label="$1"
    local env_name="$2"
    local script_file="$3"
    log_info "[Step ${step_label}] Dispatching to Mamba env: ${BOLD}${env_name}${NC}"
    log_info "[Step ${step_label}] Script: ${script_file}"
    if [[ "${DRY_RUN}" == true ]]; then
        echo -e "${YELLOW}[DRY-RUN]${NC} mamba run -n ${env_name} bash ${script_file}"
        echo -e "${YELLOW}[DRY-RUN]${NC} --- script contents ---"
        cat "${script_file}"
        echo -e "${YELLOW}[DRY-RUN]${NC} --- end of script ---"
        return 0
    fi
    mamba run -n "${env_name}" bash "${script_file}"
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
        if [[ "${KRAKEN2_DB}" == *"/path/to/"* ]] || [[ ! -d "${KRAKEN2_DB}" ]]; then
            log_error "KRAKEN2_DB not set or not found: ${KRAKEN2_DB}"
            (( errors++ )) || true
        fi
    fi

    # CheckM2 DB
    if [[ "${CHECKM2_DB}" == *"/path/to/"* ]] || [[ ! -f "${CHECKM2_DB}" ]]; then
        log_error "CHECKM2_DB not set or not found: ${CHECKM2_DB}"
        (( errors++ )) || true
    fi

    # GTDB-Tk data
    if [[ "${RUN_GTDBTK}" == true ]]; then
        if [[ "${GTDBTK_DATA_PATH}" == *"/path/to/"* ]] || [[ ! -d "${GTDBTK_DATA_PATH}" ]]; then
            log_error "GTDBTK_DATA_PATH not set or not found: ${GTDBTK_DATA_PATH}"
            (( errors++ )) || true
        fi
    fi

    # MetaBinner scripts directory
    if [[ "${METABINNER_SCRIPTS_DIR}" == *"/path/to/"* ]] || [[ ! -d "${METABINNER_SCRIPTS_DIR}" ]]; then
        log_error "METABINNER_SCRIPTS_DIR not set or not found: ${METABINNER_SCRIPTS_DIR}"
        (( errors++ )) || true
    fi

    # Required direct binaries
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

    # mamba
    if ! command -v mamba &>/dev/null; then
        log_error "'mamba' not found in PATH."
        (( errors++ )) || true
    fi

    # Detect reads
    shopt -s nullglob
    local _probe=("${READS_DIR}"/*${READ_SUFFIX_R1})
    shopt -u nullglob
    if [[ ${#_probe[@]} -eq 0 ]]; then
        log_error "No reads matching *${READ_SUFFIX_R1} found in ${READS_DIR}"
        (( errors++ )) || true
    else
        log_info "Found ${#_probe[@]} read pair(s) in ${READS_DIR}"
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
        run_cmd fastp \
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

    run_cmd megahit \
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
        run_cmd kraken2 \
            --db "${KRAKEN2_DB}" \
            --threads "${THREADS}" \
            --paired "${r1}" "${r2}" \
            --output "${kraken_out}" \
            --report "${kraken_report}" \
            2>&1 | tee "${kraken_log}"

        log_info "[Step 03] Running Bracken for sample: ${base}"
        run_cmd bracken \
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
    run_cmd bowtie2-build \
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
        run_cmd "bowtie2 --very-sensitive-local -p ${THREADS} --time --reorder \
            -x ${bt2_index} \
            -1 ${r1} \
            -2 ${r2} \
            2>>${bt2_log} \
            | samtools view -@ ${THREADS} -bS - \
            | samtools sort -@ ${THREADS} -o ${bam_sorted} -"

        log_info "[Step 04] Indexing BAM: ${bam_sorted}"
        run_cmd samtools index -@ "${THREADS}" "${bam_sorted}"

        SORTED_BAMS+=("${bam_sorted}")
    done

    DEPTH_FILE="${BINNING_DIR}/depth_matrix.txt"
    local jgi_log="${LOGS_DIR}/step_04_jgi_depth.log"

    log_info "[Step 04] Running jgi_summarize_bam_contig_depths..."
    run_cmd jgi_summarize_bam_contig_depths \
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
metaWRAP reassemble_bins \\
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
        run_cmd coverm genome \
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

    run_cmd Rscript - <<RSCRIPT
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

    # Print run configuration banner
    echo -e "\n${CYAN}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}${BOLD}║         SHOTGUN METAGENOMICS PIPELINE — CONFIGURATION         ║${NC}"
    echo -e "${CYAN}${BOLD}╠══════════════════════════════════════════════════════════════╣${NC}"
    echo -e "${CYAN}${BOLD}║${NC}  READS_DIR              : ${READS_DIR}"
    echo -e "${CYAN}${BOLD}║${NC}  OUTDIR                 : ${OUTDIR}"
    echo -e "${CYAN}${BOLD}║${NC}  ASSEMBLY_IN            : ${ASSEMBLY_IN:-<run MEGAHIT>}"
    echo -e "${CYAN}${BOLD}║${NC}  THREADS                : ${THREADS}"
    echo -e "${CYAN}${BOLD}║${NC}  MEMORY (GB)            : ${MEMORY}"
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
