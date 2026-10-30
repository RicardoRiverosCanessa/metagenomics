# ==============================================================================
# Makefile - Autonomous Apptainer/Singularity (.sif) Container Management
# ==============================================================================

APPTAINER   ?= $(shell command -v apptainer 2>/dev/null || command -v singularity 2>/dev/null)
BUILD_FLAGS ?= --fakeroot
SIF_DIR     := ./containers
DEF_DIR     := ./defs
TMP_DIR     := $(CURDIR)/.apptainer_tmp
PIPELINE    := ./run_metagenomics_pipeline.sh
READS_DIR ?=
CLEAN_READS_DIR ?=
ASSEMBLY_IN ?=
RUN_MEGAHIT ?=
RUN_FASTP ?=
RUN_KRAKEN2 ?=
RUN_MAPPING ?= true
RUN_BINNING ?=
RUN_REASSEMBLY ?=
RUN_CHECKM2 ?=
RUN_MAG_FILTERING ?=
RUN_GTDBTK ?=
RUN_COVERM ?=
RUN_SUMMARY ?=
BINNERS ?= vamb,metabinner,comebin
COMEBIN_DEVICE ?= cpu
COMEBIN_SEED ?= 42
VAMB_BINSPLIT_SEPARATOR ?=
REFINEMENT_MODE ?= auto
BINS_SOURCE ?= auto
PROFILE ?= full
KRAKEN2_DB ?=
CHECKM2_DB ?=
GTDBTK_DATA_PATH ?=
OUTDIR ?= ./metagenomics_out_test
BAM_DIR ?=
BAM_MANIFEST ?=
JGI_DEPTH_FILE ?=
FORCE_STAGE ?=
RESTART_FROM ?=

READ_INPUT_ARGS = $(if $(strip $(READS_DIR)),--reads-dir "$(READS_DIR)") \
                  $(if $(strip $(CLEAN_READS_DIR)),--clean-reads-dir "$(CLEAN_READS_DIR)")
ASSEMBLY_ARGS = $(if $(strip $(ASSEMBLY_IN)),-a "$(ASSEMBLY_IN)")
MEGAHIT_ARGS = $(if $(filter no false 0,$(RUN_MEGAHIT)),--no-megahit) \
               $(if $(filter yes true 1,$(RUN_MEGAHIT)),--megahit)
FASTP_STAGE_ARGS = $(if $(filter no false 0,$(RUN_FASTP)),--no-fastp) \
                   $(if $(filter yes true 1,$(RUN_FASTP)),--fastp)
KRAKEN2_ARGS = $(if $(filter no false 0,$(RUN_KRAKEN2)),--no-kraken2) \
               $(if $(filter yes true 1,$(RUN_KRAKEN2)),--kraken2)
# Do not let the implicit true default override a profile such as taxonomy.
MAPPING_ARGS = $(if $(filter command% environment override,$(origin RUN_MAPPING)),\
                 $(if $(filter false,$(RUN_MAPPING)),--no-mapping,--mapping))
BINNING_ARGS = $(if $(filter no false 0,$(RUN_BINNING)),--no-binning) \
               $(if $(filter yes true 1,$(RUN_BINNING)),--binning)
REASSEMBLY_ARGS = $(if $(filter no false 0,$(RUN_REASSEMBLY)),--no-reassembly) \
                  $(if $(filter yes true 1,$(RUN_REASSEMBLY)),--reassembly)
CHECKM2_ARGS = $(if $(filter no false 0,$(RUN_CHECKM2)),--no-checkm2) \
               $(if $(filter yes true 1,$(RUN_CHECKM2)),--checkm2)
MAG_FILTERING_ARGS = $(if $(filter no false 0,$(RUN_MAG_FILTERING)),--no-mag-filtering) \
                     $(if $(filter yes true 1,$(RUN_MAG_FILTERING)),--mag-filtering)
GTDBTK_ARGS = $(if $(filter no false 0,$(RUN_GTDBTK)),--no-gtdbtk) \
              $(if $(filter yes true 1,$(RUN_GTDBTK)),--gtdbtk)
COVERM_ARGS = $(if $(filter no false 0,$(RUN_COVERM)),--no-coverm) \
              $(if $(filter yes true 1,$(RUN_COVERM)),--coverm)
SUMMARY_ARGS = $(if $(filter no false 0,$(RUN_SUMMARY)),--no-summary) \
               $(if $(filter yes true 1,$(RUN_SUMMARY)),--summary)
DATABASE_ARGS = $(if $(strip $(KRAKEN2_DB)),--kraken2-db "$(KRAKEN2_DB)") \
                $(if $(strip $(CHECKM2_DB)),--checkm2-db "$(CHECKM2_DB)") \
                $(if $(strip $(GTDBTK_DATA_PATH)),--gtdbtk-data-path "$(GTDBTK_DATA_PATH)")
FASTP_ARGS = $(if $(strip $(CLEAN_READS_DIR)),--no-fastp)
BAM_ARGS = $(if $(strip $(BAM_DIR)),--bam-dir "$(BAM_DIR)") \
           $(if $(strip $(BAM_MANIFEST)),--bam-manifest "$(BAM_MANIFEST)") \
           $(if $(strip $(JGI_DEPTH_FILE)),--jgi-depth-file "$(JGI_DEPTH_FILE)")
shell_quote = '$(subst ','"'"',$(1))'
VAMB_ARGS = $(if $(VAMB_BINSPLIT_SEPARATOR),\
             --vamb-binsplit-separator $(call shell_quote,$(VAMB_BINSPLIT_SEPARATOR)))
comma := ,
FORCE_STAGE_ARGS = $(foreach stage,$(subst $(comma), ,$(strip $(FORCE_STAGE))),--force-stage "$(stage)")
RESTART_FROM_ARGS = $(if $(strip $(RESTART_FROM)),--restart-from "$(RESTART_FROM)")
COMEBIN_ARGS = --comebin-device "$(COMEBIN_DEVICE)" \
               --comebin-seed "$(COMEBIN_SEED)"

# Redirect Apptainer/Singularity temp and cache dirs to local disk instead of /tmp RAM
export APPTAINER_TMPDIR     := $(TMP_DIR)/tmp
export APPTAINER_CACHEDIR   := $(TMP_DIR)/cache
export SINGULARITY_TMPDIR   := $(TMP_DIR)/tmp
export SINGULARITY_CACHEDIR := $(TMP_DIR)/cache

# Clean Python helper to query the latest active tag from Quay.io JSON API
GET_QUAY_TAG = python3 -c "import urllib.request, json, sys; u='https://quay.io/api/v1/repository/biocontainers/'+sys.argv[1]+'/tag/?limit=20&onlyActiveTags=true'; print(json.loads(urllib.request.urlopen(u).read())['tags'][0]['name'])"

LOCAL_SIFS := $(SIF_DIR)/kraken2_bracken.sif $(SIF_DIR)/mapping_depth.sif $(SIF_DIR)/metabinner.sif $(SIF_DIR)/comebin.sif $(SIF_DIR)/vamb.sif
REMOTE_SIFS := $(SIF_DIR)/fastp.sif $(SIF_DIR)/megahit.sif $(SIF_DIR)/metawrap.sif $(SIF_DIR)/checkm2.sif $(SIF_DIR)/gtdbtk.sif $(SIF_DIR)/coverm.sif $(SIF_DIR)/r_base.sif
ALL_SIFS := $(LOCAL_SIFS) $(REMOTE_SIFS)

.PHONY: all setup check-env pull-containers build-containers test-containers test-dry-run clean-tmp clean help

all: help

## setup: Check Apptainer installation, build/pull ALL .sif containers, and clean temp cache
setup: check-env $(ALL_SIFS) clean-tmp
	@echo ''
	@echo '[OK] All 12 (.sif) containers are ready in $(SIF_DIR)/'
	@echo 'Run make test-containers to verify binaries, or make test-dry-run to validate the pipeline.'

check-env:
	@if [ -z '$(APPTAINER)' ]; then echo '[ERROR] Neither apptainer nor singularity found in PATH.'; exit 1; fi
	@echo 'Container engine detected: $(APPTAINER)'
	@mkdir -p $(SIF_DIR) $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)

$(LOCAL_SIFS): $(SIF_DIR)/%.sif: $(DEF_DIR)/%.def | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@echo 'Building $@ from $<...'
	$(APPTAINER) build $(BUILD_FLAGS) $@ $<
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/fastp.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) fastp); echo "Pulling fastp:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/fastp:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/megahit.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) megahit); echo "Pulling megahit:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/megahit:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/metawrap.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) metawrap-mg); echo "Pulling metawrap-mg:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/metawrap-mg:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/checkm2.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) checkm2); echo "Pulling checkm2:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/checkm2:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/gtdbtk.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) gtdbtk); echo "Pulling gtdbtk:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/gtdbtk:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/coverm.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	@TAG=$$($(GET_QUAY_TAG) coverm); echo "Pulling coverm:$$TAG -> $@"; $(APPTAINER) pull $@ docker://quay.io/biocontainers/coverm:$$TAG
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

$(SIF_DIR)/r_base.sif: | check-env
	@mkdir -p $(APPTAINER_TMPDIR) $(APPTAINER_CACHEDIR)
	$(APPTAINER) pull $@ docker://rocker/r-base:4.3.2
	@rm -rf $(APPTAINER_TMPDIR)/* $(APPTAINER_CACHEDIR)/*

## build-containers: Build only the 5 local (.def) container images
build-containers: $(LOCAL_SIFS) clean-tmp

## pull-containers: Pull only the 7 remote BioContainers/DockerHub images
pull-containers: $(REMOTE_SIFS) clean-tmp

## test-containers: Verify that all 12 .sif images exist and their internal binaries respond
test-containers: check-env
	@echo '=========================================================='
	@echo 'Verifying all 12 Apptainer (.sif) containers...'
	@echo '=========================================================='
	@FAIL=0; \
	check_sif() { \
		SIF="$(SIF_DIR)/$$1"; CMD="$$2"; \
		if [ ! -f "$$SIF" ]; then \
			echo "[MISSING] $$SIF not found"; FAIL=1; \
		else \
			if $(APPTAINER) exec --cleanenv "$$SIF" bash -c "$$CMD" >/dev/null 2>&1; then \
				echo "[OK]      $$1 ($$CMD)"; \
			else \
				echo "[FAIL]    $$1 failed running: $$CMD"; FAIL=1; \
			fi; \
		fi; \
	}; \
	check_sif fastp.sif            'fastp --version'; \
	check_sif megahit.sif          'megahit --version'; \
	check_sif kraken2_bracken.sif  'kraken2 --version && which bracken'; \
	check_sif mapping_depth.sif    'bowtie2 --version && samtools --version && which jgi_summarize_bam_contig_depths'; \
	check_sif vamb.sif             'vamb --version'; \
	check_sif metabinner.sif       'which run_metabinner.sh && test -f /opt/metabinner/scripts/Filter_tooshort.py'; \
	check_sif comebin.sif          'python -c "import torch" && which run_comebin.sh'; \
	check_sif metawrap.sif         'metawrap --version || which metawrap'; \
	check_sif checkm2.sif          'checkm2 --version'; \
	check_sif gtdbtk.sif           'gtdbtk --version'; \
	check_sif coverm.sif           'coverm --version && samtools --version'; \
	check_sif r_base.sif           'Rscript --version'; \
	rm -rf $(TMP_DIR); \
	if [ "$$FAIL" -ne 0 ]; then \
		echo '----------------------------------------------------------'; \
		echo '[ERROR] One or more containers failed verification.'; \
		exit 1; \
	else \
		echo '----------------------------------------------------------'; \
		echo '[SUCCESS] All 12 containers are valid and ready!'; \
	fi

## clean-tmp: Remove temporary Apptainer build and cache directory (.apptainer_tmp)
clean-tmp:
	@echo 'Cleaning temporary build/cache directory ($(TMP_DIR))...'
	@rm -rf $(TMP_DIR)

## test-dry-run: Check Bash syntax and execute a pipeline dry-run (--dry-run)
test-dry-run:
	@case "$(RUN_MAPPING)" in true|false) ;; *) echo '[ERROR] RUN_MAPPING must be true or false.'; exit 1 ;; esac
	@case $(call shell_quote,$(COMEBIN_DEVICE)) in \
	  cpu|cuda) ;; \
	  cuda:*) device_index=$(call shell_quote,$(patsubst cuda:%,%,$(COMEBIN_DEVICE))); \
	    case "$$device_index" in ''|*[!0-9]*) echo '[ERROR] COMEBIN_DEVICE must be cpu, cuda, or cuda:N with N a non-negative integer.'; exit 1 ;; esac ;; \
	  *) echo '[ERROR] COMEBIN_DEVICE must be cpu, cuda, or cuda:N with N a non-negative integer.'; exit 1 ;; \
	esac
	@case $(call shell_quote,$(COMEBIN_SEED)) in ''|*[!0-9]*) echo '[ERROR] COMEBIN_SEED must be a non-negative integer.'; exit 1 ;; esac
	bash -n $(PIPELINE)
	bash $(PIPELINE) \
	  $(READ_INPUT_ARGS) \
	  --profile "$(PROFILE)" \
	  $(ASSEMBLY_ARGS) \
	  $(MEGAHIT_ARGS) \
	  $(FASTP_STAGE_ARGS) \
	  $(KRAKEN2_ARGS) \
	  $(MAPPING_ARGS) \
	  --binners "$(BINNERS)" \
	  $(COMEBIN_ARGS) \
	  $(VAMB_ARGS) \
	  --refinement-mode "$(REFINEMENT_MODE)" \
	  --bins-source "$(BINS_SOURCE)" \
	  $(BINNING_ARGS) \
	  $(REASSEMBLY_ARGS) \
	  $(CHECKM2_ARGS) \
	  $(MAG_FILTERING_ARGS) \
	  $(GTDBTK_ARGS) \
	  $(COVERM_ARGS) \
	  $(SUMMARY_ARGS) \
	  $(FASTP_ARGS) \
	  $(BAM_ARGS) \
	  $(FORCE_STAGE_ARGS) \
	  $(RESTART_FROM_ARGS) \
	  $(DATABASE_ARGS) \
	  --outdir "$(OUTDIR)" \
	  --dry-run

## clean: Remove all generated .sif containers and temp files
clean: clean-tmp
	rm -rf $(SIF_DIR)/*.sif

## help: Show available Makefile commands
help:
	@echo 'Usage: make <target>'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^## ' $(MAKEFILE_LIST) | sed 's/## /  /'
	@echo ''
	@echo 'Example:'
	@echo '  make test-dry-run \'
	@echo '    PROFILE=full READS_DIR=/data/raw_reads \'
	@echo '    KRAKEN2_DB=/databases/kraken2 \'
	@echo '    CHECKM2_DB=/databases/checkm2/uniref100.KO.1.dmnd \'
	@echo '    GTDBTK_DATA_PATH=/databases/gtdbtk'
	@echo ''
	@echo 'External assembly (explicitly disable MEGAHIT):'
	@echo '  make test-dry-run PROFILE=assembly RUN_MEGAHIT=no \'
	@echo '    ASSEMBLY_IN=/data/assembly/final.contigs.fa CLEAN_READS_DIR=/data/clean_reads'
	@echo ''
	@echo 'Binner selection dry-runs:'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb CLEAN_READS_DIR=/data/clean_reads  # binsplitting disabled by default'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb VAMB_BINSPLIT_SEPARATOR=C CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=metabinner CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=comebin CLEAN_READS_DIR=/data/clean_reads'
	@echo '  CPU:'
	@echo '  make test-dry-run \'
	@echo '    BINNERS=comebin \'
	@echo '    COMEBIN_DEVICE=cpu'
	@echo '  CUDA:'
	@echo '  make test-dry-run \'
	@echo '    BINNERS=comebin \'
	@echo '    COMEBIN_DEVICE=cuda:0'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb,metabinner REFINEMENT_MODE=auto BINS_SOURCE=auto CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb,metabinner,comebin REFINEMENT_MODE=auto BINS_SOURCE=auto CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly BINNERS=vamb,comebin REFINEMENT_MODE=never BINS_SOURCE=comebin CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=custom RUN_FASTP=no RUN_KRAKEN2=no RUN_MEGAHIT=no RUN_MAPPING=no RUN_BINNING=no RUN_CHECKM2=no RUN_MAG_FILTERING=no RUN_GTDBTK=no RUN_COVERM=no RUN_SUMMARY=no'
	@echo ''
	@echo 'External coverage and selective resume:'
	@echo '  make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads'
	@echo '  make test-dry-run PROFILE=assembly RUN_MAPPING=false BAM_DIR=/data/bams BINNERS=comebin'
	@echo '  make test-dry-run PROFILE=assembly RUN_MAPPING=false JGI_DEPTH_FILE=/data/depth.txt BINNERS=vamb'
	@echo '  make test-dry-run PROFILE=assembly RUN_MAPPING=false JGI_DEPTH_FILE=/data/depth.txt BINNERS=vamb FORCE_STAGE=vamb'
	@echo '  make test-dry-run PROFILE=assembly RUN_MAPPING=false BAM_MANIFEST=/data/bams.tsv BINNERS=comebin RESTART_FROM=comebin'
