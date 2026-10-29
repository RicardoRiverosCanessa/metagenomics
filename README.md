# Shotgun Metagenomics Pipeline

A reproducible Bash pipeline for paired-end shotgun metagenomic reads, from quality control and assembly to metagenome-assembled genome (MAG) recovery, quality assessment, taxonomy, abundance estimation, and final summary generation.

The workflow supports two execution modes:

- **Apptainer/Singularity**: the default and recommended mode, using isolated `.sif` images.
- **Mamba**: an alternative mode using existing host binaries and named Conda/Mamba environments.

> **Status:** This pipeline is under active development. Validate it with representative test data before using it in production analyses.

## Workflow

| Step | Process | Main software | Optional |
|---:|---|---|:---:|
| 01 | Read quality control | fastp | Yes |
| 02 | Metagenome assembly | MEGAHIT | Yes |
| 03 | Read-level taxonomic profiling | Kraken2 + Bracken | Yes |
| 04 | Read mapping and contig depth | Bowtie2, SAMtools, JGI depth | Yes |
| 05a | Binning | Vamb | Yes |
| 05b | Binning | MetaBinner | Yes |
| 05c | Binning | COMEBin | Yes |
| 06 | Bin refinement | MetaWRAP | Yes |
| 07 | Bin reassembly | MetaWRAP | Yes |
| 08 | MAG quality assessment | CheckM2 | Yes |
| 09 | MAG quality filtering | Bash/AWK | Yes |
| 10 | MAG taxonomic classification | GTDB-Tk | Yes |
| 11 | MAG abundance estimation | CoverM | Yes |
| 12 | Integrated summary table | Base R | Yes |

All stages are enabled by default except MetaWRAP reassembly. Mapping, binning,
refinement, CheckM2, MAG filtering, CoverM, and summary can be disabled
independently when their dependencies are also disabled as described below.

## Repository structure

```text
.
├── run_metagenomics_pipeline.sh
├── Makefile
├── defs/
│   ├── kraken2_bracken.def
│   ├── mapping_depth.def
│   ├── metabinner.def
│   ├── comebin.def
│   └── vamb.def
└── containers/                 # Created automatically
```

The five local `.def` files are required by the Makefile targets that build the corresponding container images.

## Requirements

### Recommended: Apptainer mode

- Linux
- Bash
- GNU Make
- Python 3, used by the Makefile to query container tags
- Apptainer or Singularity
- Internet access during initial image build/download
- Sufficient disk space for container images, reference databases, intermediate BAM files, bins, and assemblies

Apptainer builds use `--fakeroot` by default. If `--fakeroot` is unavailable on your system, override the Make variable with an appropriate site-specific build option.

### Alternative: Mamba mode

Mamba mode requires:

- `mamba` available in `PATH`
- Host-accessible binaries for fastp, MEGAHIT, Kraken2, Bracken, Bowtie2, SAMtools, JGI depth, and CoverM as applicable
- Named environments configured in the script for Vamb, MetaBinner, COMEBin, MetaWRAP, CheckM2, and GTDB-Tk
- A host path to the MetaBinner scripts supplied with `--metabinner-scripts-dir`

## Input files

The reads directory must contain paired, gzip-compressed FASTQ files using these suffixes:

```text
sample_01_R1.fastq.gz
sample_01_R2.fastq.gz
sample_02_R1.fastq.gz
sample_02_R2.fastq.gz
```

For every `*_R1.fastq.gz` file, a corresponding `*_R2.fastq.gz` file must exist.

## Required databases

The pipeline accepts database paths as command-line arguments:

| Argument | Expected path | Required when |
|---|---|---|
| `--kraken2-db` | Kraken2/Bracken database directory | Kraken2 is enabled |
| `--checkm2-db` | CheckM2 DIAMOND database file, such as `uniref100.KO.1.dmnd` | CheckM2 is enabled and reached |
| `--gtdbtk-data-path` | GTDB-Tk reference-data directory | GTDB-Tk is enabled |

The Kraken2 database must contain the files needed by Kraken2 and the Bracken files matching the configured read length and taxonomic level.

## Installation

Clone the repository and enter its directory:

```bash
git clone https://github.com/RicardoRiverosCanessa/metagenomics
cd metagenomics
chmod +x run_metagenomics_pipeline.sh
```

Check that Apptainer or Singularity is available:

```bash
make check-env
```

Build or download all 12 container images:

```bash
make setup
```

Alternatively, manage the two image groups separately:

```bash
make build-containers   # Build the five images defined locally
make pull-containers    # Download the seven remote images
```

Verify the containers and their principal commands:

```bash
make test-containers
```

## Basic usage

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out
```

The script uses Apptainer mode by default. Missing remote images may be pulled automatically when `AUTO_PULL_SIF=true`.

## Command-line options

### Paths and resources

| Option | Description |
|---|---|
| `-r, --reads-dir DIR` | Directory containing paired-end reads |
| `--kraken2-db DIR` | Kraken2/Bracken database directory |
| `--checkm2-db FILE` | CheckM2 database file |
| `--gtdbtk-data-path DIR` | GTDB-Tk reference-data directory |
| `-o, --outdir DIR` | Output directory; default: `./metagenomics_out` |
| `--metabinner-scripts-dir DIR` | MetaBinner scripts path; mainly required in Mamba mode |
| `-t, --threads INT` | Number of CPU threads |
| `-m, --memory INT` | RAM limit in GB |
| `-e, --engine MODE` | Execution mode: `apptainer` or `mamba` |
| `-h, --help` | Display command-line help |

When threads or memory remain set to `auto`, the pipeline detects available resources. It reserves two logical CPU cores on hosts with more than four cores and uses up to 85% of available RAM.

### Enable or disable optional stages

| Enable | Disable | Default |
|---|---|---:|
| `--fastp` | `--no-fastp` | Enabled |
| `--megahit` | `--no-megahit` | Enabled |
| `--kraken2` | `--no-kraken2` | Enabled |
| `--reassembly` | `--no-reassembly` | Disabled |
| `--gtdbtk` | `--no-gtdbtk` | Enabled |
| `--mapping` | `--no-mapping` | Enabled |
| `--binning` | `--no-binning` | Enabled |
| `--refinement` | `--no-refinement` | Enabled |
| `--checkm2` | `--no-checkm2` | Enabled |
| `--mag-filtering` | `--no-mag-filtering` | Enabled |
| `--coverm` | `--no-coverm` | Enabled |
| `--summary` | `--no-summary` | Enabled |
| `--dry-run` | `--no-dry-run` | Disabled |
| `--force` | `--no-force` | Disabled |

Command-line arguments override environment variables and script defaults.

### Stage dependencies and stopping early

The pipeline checks the requested configuration during preflight and does not
silently enable or disable stages. The supported dependencies are:

- Binning requires mapping.
- Refinement requires binning.
- Reassembly requires refinement.
- CheckM2 requires refinement or reassembly.
- MAG filtering requires CheckM2.
- GTDB-Tk and CoverM each require MAG filtering.

Use `--stop-after <stage>` to finish successfully after a stage has completed
and its outputs have been validated. Accepted values are `mapping`, `binning`,
`refinement`, `reassembly`, `checkm2`, `mag-filtering`, `gtdbtk`, `coverm`, and
`summary`. The selected stage must be enabled and its dependencies must be
satisfied. For example:

```bash
bash run_metagenomics_pipeline.sh --reads-dir /data/reads \
  --kraken2-db /databases/kraken2 --stop-after binning
```

The final summary uses the available MAG FASTA list as its base and adds
CheckM2, GTDB-Tk, and CoverM columns when their result files are available.
Unavailable source fields remain `NA`; a summary is omitted with a warning
when no non-empty MAGs are available.

## Examples

### Validate a complete run without executing tools

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out_test \
  --dry-run
```

The dry run still performs pre-flight validation of reads, database paths, execution mode, and required configuration.

### Use the Makefile dry-run target

```bash
make test-dry-run \
  READS_DIR=/data/project/reads \
  KRAKEN2_DB=/databases/kraken2 \
  CHECKM2_DB=/databases/checkm2/uniref100.KO.1.dmnd \
  GTDBTK_DATA_PATH=/databases/gtdbtk \
  OUTDIR=./metagenomics_out_test
```

This target first runs `bash -n` and then launches the pipeline with `--dry-run`.

### Skip Kraken2 and GTDB-Tk

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --outdir ./metagenomics_out \
  --no-kraken2 \
  --no-gtdbtk
```

When these stages are disabled, their database paths are not required.

### Enable bin reassembly

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out \
  --reassembly
```

Reassembly concatenates the available paired reads before running MetaWRAP `reassemble_bins`, so ensure sufficient temporary and output storage.

### Use raw reads without fastp

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out \
  --no-fastp
```

### Force all stages to run again

```bash
bash run_metagenomics_pipeline.sh \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out \
  --force
```

### Run with Mamba

```bash
bash run_metagenomics_pipeline.sh \
  --engine mamba \
  --reads-dir /data/project/reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --metabinner-scripts-dir /opt/software/MetaBinner/scripts \
  --outdir ./metagenomics_out
```

## Output structure

```text
metagenomics_out/
├── 01_qc/
├── 02_assembly/
├── 03_kraken2/
├── 04_mapping/
├── 05_binning/
│   ├── vamb/
│   ├── metabinner/
│   └── comebin/
├── 06_refinement/
├── 07_reassembly/
├── 08_checkm2/
├── 09_mags/
├── 10_gtdbtk/
├── 11_abundance/
├── 12_summary/
│   └── final_summary.tsv
├── logs/
└── .checkpoints/
```

Important outputs include:

- `02_assembly/final.contigs.fa`: filtered assembly.
- `03_kraken2/`: per-sample Kraken2 and Bracken results.
- `04_mapping/`: sorted and indexed BAM files.
- `05_binning/`: bins from Vamb, MetaBinner, and COMEBin.
- `06_refinement/`: MetaWRAP-refined bins.
- `08_checkm2/quality_report.tsv`: MAG quality estimates.
- `09_mags/`: MAGs meeting the configured quality thresholds.
- `10_gtdbtk/`: GTDB-Tk classifications.
- `11_abundance/`: per-sample CoverM abundance tables.
- `12_summary/final_summary.tsv`: MAG-based summary with available CheckM2, GTDB-Tk, and CoverM data.
- `logs/`: step-specific execution logs and generated helper scripts.

## MAG filtering defaults

The pipeline retains bins meeting both thresholds:

```text
Completeness >= 50%
Contamination <= 10%
```

These values are currently configured inside the script as `COMPLETENESS` and `CONTAMINATION`.

## Checkpoints and reproducibility

After successful completion, each stage writes a marker under `.checkpoints/`. On later runs, completed stages are skipped and their expected outputs are reused.

Use `--force` to ignore all checkpoints and rerun the pipeline. Because force mode can overwrite or remove some previous stage outputs, use a new output directory when preserving an earlier analysis is important.

## Pre-flight validation

Before running the workflow, the script checks:

- existence of the reads directory;
- presence of at least one `*_R1.fastq.gz` file;
- availability of the matching R2 file for every R1 file;
- required database paths according to enabled stages;
- the selected Apptainer/Singularity or Mamba runtime;
- required host binaries in Mamba mode;
- MetaBinner script configuration;
- compatibility of assembly settings.

## Makefile targets

| Target | Purpose |
|---|---|
| `make help` | Show available targets |
| `make check-env` | Check Apptainer/Singularity and initialize directories |
| `make setup` | Build/download all containers and clean temporary cache |
| `make build-containers` | Build the five locally defined images |
| `make pull-containers` | Download the seven remote images |
| `make test-containers` | Verify all 12 images and their principal commands |
| `make test-dry-run ...` | Check Bash syntax and run pre-flight validation |
| `make clean-tmp` | Remove local container build/cache files |
| `make clean` | Remove generated `.sif` images and temporary files |

> `make clean` deletes all `.sif` files under `containers/`.

## Troubleshooting

### No reads detected

Confirm that the files are directly inside the supplied reads directory and follow the exact naming convention:

```text
*_R1.fastq.gz
*_R2.fastq.gz
```

### Database path rejected

- `--kraken2-db` expects a directory.
- `--checkm2-db` expects the CheckM2 `.dmnd` file itself when CheckM2 is enabled and reached.
- `--gtdbtk-data-path` expects the GTDB-Tk reference-data directory.

### MetaBinner scripts not found

In Apptainer mode, the default internal path is:

```text
/opt/metabinner/scripts
```

In Mamba mode, supply the host path explicitly:

```bash
--metabinner-scripts-dir /path/to/MetaBinner/scripts
```

### Assembly disabled

The current command-line interface can disable MEGAHIT with `--no-megahit`, but a valid alternative assembly must be available to the pipeline. Review the script configuration for `ASSEMBLY_IN` before using this mode.

### Container build consumes `/tmp`

The Makefile redirects Apptainer and Singularity temporary and cache files to:

```text
./.apptainer_tmp/
```

Run the following to remove those temporary files:

```bash
make clean-tmp
```

## Known implementation considerations

- The pipeline is designed for Linux and relies on Bash-specific features.
- The current read suffixes are fixed to `_R1.fastq.gz` and `_R2.fastq.gz`.
- Mapping, binning, refinement, CheckM2, MAG filtering, CoverM, and summary are enabled by default but can be disabled subject to their dependencies.
- The final summary uses a MAG list as its base and can include CheckM2, GTDB-Tk, and CoverM data independently.
- Container tags and tool compatibility should be validated before long production runs.

## Citation

If you use this workflow in a publication, cite the individual tools used in your selected configuration, including fastp, MEGAHIT, Kraken2, Bracken, Bowtie2, SAMtools, Vamb, MetaBinner, COMEBin, MetaWRAP, CheckM2, GTDB-Tk, CoverM, Apptainer, and R as applicable.


## License

No license is currently specified.

## Contributing

Bug reports and pull requests are welcome. When reporting an issue, include:

- the exact command used;
- the selected execution engine;
- the relevant log file from `OUTDIR/logs/`;
- the pipeline commit or release;
- software/container versions;
- a minimal description of the input naming structure, without sharing sensitive sequence data.
