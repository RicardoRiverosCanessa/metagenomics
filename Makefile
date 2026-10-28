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
KRAKEN2_DB ?=
CHECKM2_DB ?=
GTDBTK_DATA_PATH ?=
OUTDIR ?= ./metagenomics_out_test

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
	bash -n $(PIPELINE)
	@test -n "$(READS_DIR)" || { echo '[ERROR] Set READS_DIR.'; exit 1; }
	@test -n "$(KRAKEN2_DB)" || { echo '[ERROR] Set KRAKEN2_DB.'; exit 1; }
	@test -n "$(CHECKM2_DB)" || { echo '[ERROR] Set CHECKM2_DB.'; exit 1; }
	@test -n "$(GTDBTK_DATA_PATH)" || { echo '[ERROR] Set GTDBTK_DATA_PATH.'; exit 1; }
	bash $(PIPELINE) \
	  --reads-dir "$(READS_DIR)" \
	  --kraken2-db "$(KRAKEN2_DB)" \
	  --checkm2-db "$(CHECKM2_DB)" \
	  --gtdbtk-data-path "$(GTDBTK_DATA_PATH)" \
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
	@echo '    READS_DIR=/data/reads \'
	@echo '    KRAKEN2_DB=/databases/kraken2 \'
	@echo '    CHECKM2_DB=/databases/checkm2/uniref100.KO.1.dmnd \'
	@echo '    GTDBTK_DATA_PATH=/databases/gtdbtk'
