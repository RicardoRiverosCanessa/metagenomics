# Shotgun Metagenomics Pipeline

A reproducible Bash pipeline for paired-end shotgun metagenomic reads, from quality control and assembly to metagenome-assembled genome (MAG) recovery, quality assessment, taxonomy, abundance estimation, and final summary generation.

The workflow supports two execution modes:

- **Apptainer/Singularity**: the default and recommended mode, using isolated `.sif` images.
- **Mamba**: an alternative mode using existing host binaries and named Conda/Mamba environments.

> **Status:** This pipeline is under active development. Validate it with representative test data before using it in production analyses.

## Workflow branches

```text
                         ┌─ Kraken2 ─ Bracken (taxonomy branch)
reads ── optional fastp ─┤
                         └─ MEGAHIT ─ mapping/depth ─ selected binners
                              ─ optional consensus refinement ─ optional reassembly
                              ─ CheckM2 ─ MAG filtering
                              ─ optional GTDB-Tk ─ CoverM ─ summary
```

The branches can run together, independently, or as a manually selected
combination. Reads are collected once and shared through the active-read arrays.
Visible progress numbering follows execution order: fastp, Kraken2/Bracken,
MEGAHIT, then the assembly stages.

Kraken2/Bracken intentionally runs before MEGAHIT so taxonomic results become
available early. Its internal checkpoint and output directory remain
`03_kraken2` for compatibility; MEGAHIT retains checkpoint `02_megahit` and
directory `02_assembly`.

### Execution profiles

| Profile | Enabled stages |
|---|---|
| `full` (default) | Both branches, all assembly stages, and fastp/reassembly/GTDB-Tk according to their configured values |
| `taxonomy` | fastp according to configuration, then Kraken2/Bracken only |
| `assembly` | Assembly branch only; fastp, reassembly, and GTDB-Tk follow their configured values |
| `custom` | Keeps stage defaults/configuration and applies individual stage options |

`--taxonomy-only` and `--assembly-only` are aliases for the corresponding
profiles. Individual `--stage`/`--no-stage` options override a selected profile,
regardless of argument order. The pipeline reports incompatible final
combinations during preflight rather than silently changing stage settings.

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
├── scripts/
│   ├── generate_pipeline_summaries.py
│   ├── validate_pipeline_outputs.py
│   ├── write_benchmark_manifest.py
│   └── write_error_record.py
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

## Read inputs

Both read-directory options use paired, gzip-compressed FASTQ files directly
inside the directory, with these exact suffixes:

```text
sample_01_R1.fastq.gz
sample_01_R2.fastq.gz
sample_02_R1.fastq.gz
sample_02_R2.fastq.gz
```

For every `*_R1.fastq.gz` file, a corresponding `*_R2.fastq.gz` file must exist.
`-r, --reads-dir DIR` supplies **raw reads**. With fastp enabled, fastp reads
from this directory and the remaining stages use `OUTDIR/01_qc/trimmed`; with
`--no-fastp`, they use the raw files directly. `--clean-reads-dir DIR` supplies
reads already cleaned outside the pipeline and requires `--no-fastp`. The
external files are consumed in place; the pipeline does not copy or rename
them. The two directory options are mutually exclusive.

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
  --profile full \
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
| `-r, --reads-dir DIR` | Directory containing raw paired-end reads |
| `--clean-reads-dir DIR` | Directory containing externally cleaned paired-end reads; use with `--no-fastp` |
| `--profile PROFILE` | `full`, `taxonomy`, `assembly`, or `custom`; default: `full` |
| `--taxonomy-only` | Alias for `--profile taxonomy` |
| `--assembly-only` | Alias for `--profile assembly` |
| `--kraken2-db DIR` | Kraken2/Bracken database directory |
| `--checkm2-db FILE` | CheckM2 database file |
| `--gtdbtk-data-path DIR` | GTDB-Tk reference-data directory |
| `-a, --assembly-in FILE` | External assembly FASTA; requires `--no-megahit` |
| `-o, --outdir DIR` | Output directory; default: `./metagenomics_out` |
| `--metabinner-scripts-dir DIR` | MetaBinner scripts path; mainly required in Mamba mode |
| `--binners LIST` | Comma-separated selection: `vamb`, `metabinner`, `comebin`; default: all three |
| `--comebin-device DEVICE` | COMEBin training device: `cpu`, `cuda`, or `cuda:N`; default: `cpu` |
| `--comebin-seed INTEGER` | Non-negative COMEBin random seed; default: `42` |
| `--min-completeness VALUE` | Minimum CheckM2 completeness to retain a MAG; inclusive, default: `70` |
| `--max-contamination VALUE` | Maximum CheckM2 contamination to retain a MAG; inclusive, default: `10` |
| `--vamb-binsplit-separator CHAR` | Advanced Vamb binsplitting separator; disabled unless explicitly supplied |
| `--refinement-mode MODE` | Consensus refinement policy: `auto`, `always`, or `never` |
| `--refinement` / `--no-refinement` | Aliases for `--refinement-mode always` / `never` |
| `--bins-dir DIR` | Directory containing one final external bin collection; requires `--no-binning` |
| `--bins-source SOURCE` | `auto`, `vamb`, `metabinner`, `comebin`, `refined`, or `external`; default: `auto` |
| `--bam-dir DIR` | External `<sample>.sorted.bam` files and `.bai`/`.csi` indexes |
| `--bam-manifest FILE` | External BAM manifest with `sample_id`, `bam`, and `index` columns |
| `--jgi-depth-file FILE` | External JGI `jgi_summarize_bam_contig_depths` table |
| `--tmpdir DIR` | Base directory for temporary pipeline workspaces (creates an isolated RUN_ID directory below) |
| `--keep-temporary` | Preserve current execution's temporary workspace after success or failure |
| `--no-keep-temporary` | Remove current execution's temporary workspace during controlled cleanup (default) |
| `--resume` | Reuse valid checkpoints (default) |
| `--force-stage STAGE` | Force one stage; repeatable, rebuilds its effective descendants |
| `--restart-from STAGE` | Destructively clears one stage and its effective descendants |
| `-t, --threads INT` | Number of CPU threads |
| `-m, --memory INT` | RAM limit in GB |
| `-e, --engine MODE` | Execution mode: `apptainer` or `mamba` |
| `-h, --help` | Display command-line help |

When threads or memory remain set to `auto`, the pipeline detects available resources. It reserves two logical CPU cores on hosts with more than four cores and uses up to 85% of available RAM.

### COMEBin training device

COMEBin defaults to CPU (`COMEBIN_DEVICE=cpu`, seed `COMEBIN_SEED=42`).
`--comebin-device DEVICE` accepts `cpu`, `cuda`, or `cuda:N`, where `N` is a
zero-based non-negative CUDA device index. `cpu` works with CPU-only or
CUDA-enabled PyTorch and does not require a GPU. `cuda` lets COMEBin choose the
default CUDA device; `cuda:N` requests a specific device and is checked against
the selected runtime's visible device count.

The pipeline passes the selection to COMEBin's short-option wrapper as
`run_comebin.sh -d DEVICE`. The wrapper parses options with `getopts`;
`--device` is not valid and must not be used. CUDA requires a CUDA-built
PyTorch, a visible GPU, and `torch.cuda.is_available()` to be true. In
Apptainer mode the pipeline adds `--nv` only to the COMEBin invocation when
`cuda` or `cuda:N` is requested; CPU execution never adds `--nv`. Mamba mode
checks the selected COMEBin environment directly.

The current `defs/comebin.def` builds PyTorch with `pytorch` and `cpuonly`.
That image supports CPU training, but it does not support CUDA requests. To use
CUDA, build and select a separate COMEBin image containing CUDA-enabled
PyTorch and compatible host GPU drivers.

Examples:

### CPU

```bash
bash run_metagenomics_pipeline.sh --binners comebin \
  --comebin-device cpu
```

### CUDA

```bash
bash run_metagenomics_pipeline.sh --binners comebin \
  --comebin-device cuda:0
```

The CUDA example requires a different COMEBin image built with CUDA-enabled
PyTorch. Both the device and seed are recorded in the COMEBin checkpoint; a
change invalidates COMEBin and its effective descendants, but not independent
stages such as mapping, Vamb, or MetaBinner.

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
| `--refinement` | `--no-refinement` | `auto` (consensus refinement is resolved from the selected binners and bins source) |
| `--checkm2` | `--no-checkm2` | Enabled |
| `--mag-filtering` | `--no-mag-filtering` | Enabled |
| `--coverm` | `--no-coverm` | Enabled |
| `--summary` | `--no-summary` | Enabled |
| `--taxonomy-summary` | `--no-taxonomy-summary` | Enabled |
| `--mag-summary` | `--no-mag-summary` | Enabled |
| `--run-overview` | `--no-run-overview` | Enabled |
| `--bracken` | `--no-bracken` | Enabled |
| `--keep-temporary` | `--no-keep-temporary` | Disabled (`--no-keep-temporary`) |
| `--dry-run` | `--no-dry-run` | Disabled |
| `--force` | `--no-force` | Disabled |

Priority is **script/environment defaults < profile < individual command-line
stage options**. The profile is resolved first, so for example
`--profile assembly --no-gtdbtk` and `--no-gtdbtk --profile assembly` have the
same result. If an individual option is repeated, its last occurrence wins.

### Stage dependencies and stopping early

The pipeline checks the requested configuration during preflight and does not
silently enable or disable stages. Binning, consensus refinement, CheckM2, and
MAG filtering are distinct operations: all selected binners run; refinement
combines two or more bin collections; CheckM2 evaluates the chosen bins; MAG
filtering selects bins using CheckM2 quality thresholds. The supported
dependencies are:

- Binning requires only the coverage products consumed by each selected binner;
  these can come from internal mapping or external inputs.
- Consensus refinement requires at least two selected binners.
- Reassembly requires effective consensus refinement.
- CheckM2 requires a resolved bin collection; it can consume direct bins or
  refined bins and does not itself perform refinement.
- MAG filtering requires CheckM2.
- GTDB-Tk and CoverM each require MAG filtering.
- Mapping requires reads and an assembly. When MEGAHIT is disabled, set the
  external FASTA with `-a, --assembly-in`; Kraken2 alone does not require one.
- Summary requires selected MAGs and at least one active tabular result source.

An external assembly is never selected implicitly: `--assembly-in FILE` must
be paired with `--no-megahit`, and `--megahit` with `--assembly-in` is rejected.
When no active stage needs an assembly, the file is optional. If provided, it is
validated as a non-empty FASTA. Its path is normalized, the parent directory is
bind-mounted in Apptainer, and the input file is used in place without copying
or modifying it.

### Mapping and coverage inputs

`RUN_MAPPING` controls only whether internal mapping runs. `COVERAGE_SOURCE` is
resolved independently as `internal`, `external`, or `none`:

| Source | Resolution | Mapping stage |
|---|---|---|
| `internal` | `RUN_MAPPING=true`; produces BAMs, indexes, and JGI depth | Active |
| `external` | `RUN_MAPPING=false` and external coverage options supplied | Disabled and absent from the graph |
| `none` | `RUN_MAPPING=false` and no external coverage options supplied | Disabled |

Do not pass external coverage options with internal mapping; use
`--no-mapping` to consume external inputs. External files and databases are
mounted read-only for consumers that need them and are never copied, renamed,
or deleted. Symlinks used by COMEBin are created inside the pipeline output
directory and point to those read-only inputs.

The requirements below are derived from the current commands and drive
preflight, bind mounts, manifests, graph edges, and invalidation:

| Binner | Assembly | BAM and index | JGI depth | Reads / abundance table |
|---|---:|---:|---:|---|
| Vamb | Yes | No | Yes | No |
| MetaBinner | Yes | No | Yes | No |
| COMEBin | Yes | Yes | No | No |

Vamb and MetaBinner can use internal mapping or `--jgi-depth-file`; COMEBin
uses internal mapping or BAM inputs. This implementation has no binner that
can run with `COVERAGE_SOURCE=none`, so disable binning for a no-coverage run.
No abundance-table options are exposed because none of the current binner
commands consumes an additional external abundance table.

Vamb's abundance TSV is generated by the pipeline from the selected JGI depth
table and validated against the assembly before Vamb runs.

`--bam-dir` uses `<sample_id>.sorted.bam` and its `<bam>.bai` or `<bam>.csi`
index. When these names associate unambiguously with active samples, the
pipeline generates an internal manifest. Otherwise, supply a tab-separated
`--bam-manifest` with exactly three columns:

```text
sample_id    bam    index
```

Manifest paths can be absolute or relative to the manifest. Each BAM must be
non-empty, pass `samtools quickcheck`, and be coordinate sorted. Validation
tries host `samtools` first and uses the mapping image only as an auxiliary
validator if needed; this does not enable mapping or create mapping outputs.
External BAM names cannot prove that they match the active assembly, so
compatibility remains the user's responsibility.

`--jgi-depth-file` accepts the JGI `jgi_summarize_bam_contig_depths` table with
`contigName`, `contigLen`, and `totalAvgDepth` leading columns, followed by
numeric coverage columns. IDs must be unique and rows well formed; depth IDs
are checked against the active assembly. The file is not modified and its
absolute path, size, and mtime are recorded.

With internal mapping, the graph adds `mapping -> binner` only when that binner
consumes a mapping product. With external coverage, mapping is not a graph
node; MEGAHIT connects directly to selected binners that consume the assembly.
External assemblies are fingerprinted by their consumers and do not create a
MEGAHIT edge.

### Selecting binners and bins

`--binners` accepts one or more exact names: `vamb`, `metabinner`, and
`comebin`. Whitespace around commas is ignored, duplicates are removed, and
the normalized list is always ordered Vamb, MetaBinner, COMEBin. Every selected
method runs; `--bins-source` only determines which collection proceeds to
CheckM2 and downstream MAG stages. The pipeline does not rank or compare
binners, choose a "best" method, or run per-binner benchmarks.

### Vamb binsplitting and inputs

Vamb binsplitting is **disabled by default**. In Vamb 5.0.4, `-o` means the
binsplit separator; it is not an alias for `--outdir`. The pipeline therefore
does not pass `-o` during a normal Vamb run. It never guesses a separator from
MEGAHIT IDs and never renames contigs.

The advanced option `--vamb-binsplit-separator CHAR` enables binsplitting only
when the FASTA first-token IDs encode the source sample in this form:

```text
<sample_id><separator><contig_id>
```

For example, `S1C19` and `S2C44` are compatible with separator `C`.
`k141_60693` is **not** compatible with separator `C`: `k141` is a MEGAHIT
contig prefix, not an inferred sample ID. A MEGAHIT coassembly should normally
be run without binsplitting. Use binsplitting only for a multi-split workflow
whose input headers were explicitly prepared to preserve sample origin.
Preflight checks every FASTA header and rejects empty, multi-character,
whitespace, `/`, and `>` separators, as well as IDs where the separator is at
either edge or absent.

Vamb 5.0.4 requires a tab-separated abundance table with a header of the form
`contigname<TAB><sample_name>...` and exactly one row per FASTA record. Vamb
verifies the ordered reference hash for contigs retained by its length filter;
the pipeline additionally enforces an exact full-FASTA ID order as a safeguard.
It fails rather than inventing or silently filling abundances. Preflight also
checks unique sample and contig IDs, row widths, finite non-negative numeric
values, and FASTA/table ID correspondence before running Vamb.

`--minfasta` has a separate meaning from contig filtering: it is Vamb's minimum
bin size to emit as a FASTA file. The pipeline setting `VAMB_MINFASTA` defaults
to 1000 bp. In Vamb 5.0.4, the contig minimum is a different option, `-m`, whose
internal default is 2000 bp; this pipeline does not override it. The pipeline's
`MIN_CONTIG_LEN` (1000 bp by default) is the post-MEGAHIT assembly filter, not
the meaning of `--minfasta`. The Vamb log reports both thresholds and their
difference when MEGAHIT created the assembly. For an external assembly, the
pipeline does not claim that it applied a contig-length filter.

Examples (binsplitting remains off in the first command):

```bash
make test-dry-run PROFILE=assembly BINNERS=vamb \
  CLEAN_READS_DIR=/data/clean_reads \
  RUN_CHECKM2=no RUN_MAG_FILTERING=no RUN_GTDBTK=no RUN_COVERM=no RUN_SUMMARY=no

make test-dry-run PROFILE=assembly BINNERS=vamb \
  VAMB_BINSPLIT_SEPARATOR=C CLEAN_READS_DIR=/data/clean_reads \
  RUN_CHECKM2=no RUN_MAG_FILTERING=no RUN_GTDBTK=no RUN_COVERM=no RUN_SUMMARY=no
```

The second command is valid only when **every** first-token assembly ID follows
the required sample/separator/contig format.

`--refinement-mode auto` refines by default when at least two binners are
selected and the resolved source is `auto` or `refined`. With one selected
binner it skips refinement. `always` requires at least two binners and a
refined source; `never` requires a direct source if multiple binners are
selected. A direct `--bins-source vamb|metabinner|comebin` skips refinement
under `auto`, but all selected binners still execute. Invalid combinations
fail preflight rather than changing the requested options.
Use `--bins-dir DIR` to point at one external final bin collection. The option
requires `--no-binning` and resolves the target to `external` automatically when
`--bins-source` is left at `auto`. A direct `--bins-source external` also works
when paired with `--bins-dir DIR`; internal and external collections are not
mixed in one run. With binning disabled, dependent MAG stages must be disabled
explicitly as needed.

Normalized bins are exposed at:

```text
05_binning/vamb/final_bins/
05_binning/metabinner/final_bins/
05_binning/comebin/final_bins/
```

These contain stable `*.fa` links to the binner's real final bins; source
outputs are not renamed or modified. `BINS_FOR_QC` resolves to one of these
directories or to the MetaWRAP refined-bin directory.

Examples:

```bash
# One method: no refinement
bash run_metagenomics_pipeline.sh --profile assembly --binners vamb \
  --refinement-mode auto --bins-source auto --clean-reads-dir /data/clean_reads

# Two methods and consensus refinement
bash run_metagenomics_pipeline.sh --profile assembly --binners vamb,metabinner \
  --refinement-mode auto --bins-source auto --clean-reads-dir /data/clean_reads

# Run both methods but use Vamb bins directly, without refinement
bash run_metagenomics_pipeline.sh --profile assembly --binners vamb,metabinner \
  --refinement-mode auto --bins-source vamb --clean-reads-dir /data/clean_reads

# Three methods, refinement forced
bash run_metagenomics_pipeline.sh --profile assembly \
  --binners vamb,metabinner,comebin --refinement-mode always \
  --bins-source refined --clean-reads-dir /data/clean_reads

# Several methods without refinement: choose the direct collection explicitly
bash run_metagenomics_pipeline.sh --profile assembly \
  --binners vamb,comebin --refinement-mode never --bins-source comebin \
  --clean-reads-dir /data/clean_reads

# No binning: turn off dependent MAG stages explicitly
bash run_metagenomics_pipeline.sh --profile custom --no-fastp --no-kraken2 \
  --no-megahit --no-mapping --no-binning --no-checkm2 --no-mag-filtering \
  --no-gtdbtk --no-coverm --no-summary
```

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

### Full pipeline with raw reads

```bash
bash run_metagenomics_pipeline.sh \
  --profile full \
  --reads-dir /data/raw_reads \
  --kraken2-db /databases/kraken2 \
  --checkm2-db /databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /databases/gtdbtk \
  --outdir ./metagenomics_out_test \
  --dry-run
```

The dry run still performs pre-flight validation of reads, database paths,
execution mode, and required configuration.

### Taxonomy-only with raw reads (fastp enabled)

```bash
bash run_metagenomics_pipeline.sh --profile taxonomy \
  --reads-dir /data/raw_reads --kraken2-db /db/kraken2
```

### Taxonomy-only with raw reads (fastp disabled)

```bash
bash run_metagenomics_pipeline.sh --profile taxonomy \
  --reads-dir /data/raw_reads --no-fastp --kraken2-db /db/kraken2
```

### Taxonomy-only with externally cleaned reads

```bash
bash run_metagenomics_pipeline.sh --profile taxonomy \
  --clean-reads-dir /data/clean_reads --no-fastp \
  --kraken2-db /db/kraken2
```

### Assembly-only with externally cleaned reads

```bash
bash run_metagenomics_pipeline.sh --profile assembly \
  --clean-reads-dir /data/clean_reads --no-fastp \
  --checkm2-db /db/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /db/gtdbtk
```

### Assembly-only with an external assembly

```bash
bash run_metagenomics_pipeline.sh --profile assembly \
  --no-megahit --assembly-in /data/assembly/final.contigs.fa \
  --clean-reads-dir /data/clean_reads --no-fastp \
  --checkm2-db /db/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-data-path /db/gtdbtk
```

### Custom stage selection

```bash
bash run_metagenomics_pipeline.sh --profile custom \
  --reads-dir /data/raw_reads --fastp \
  --no-kraken2 --megahit --mapping --binning --refinement \
  --checkm2 --mag-filtering --no-gtdbtk --coverm --summary \
  --checkm2-db /db/checkm2/uniref100.KO.1.dmnd
```

### Use the Makefile dry-run target

```bash
make test-dry-run \
  PROFILE=full \
  READS_DIR=/data/raw_reads \
  KRAKEN2_DB=/databases/kraken2 \
  CHECKM2_DB=/databases/checkm2/uniref100.KO.1.dmnd \
  GTDBTK_DATA_PATH=/databases/gtdbtk \
  OUTDIR=./metagenomics_out_test
```

This target first runs `bash -n` and then launches the pipeline with `--dry-run`.
Read directories and database arguments are passed only when their Make
variables are non-empty, so an inactive branch does not require its database.
When `CLEAN_READS_DIR` is set, the target also supplies `--no-fastp`.
`RUN_MAPPING` accepts only `true` or `false`. Its implicit Makefile default
does not override a selected profile; setting `RUN_MAPPING=true|false`
explicitly passes `--mapping` or `--no-mapping`. External coverage values are
passed only when their Make variables are non-empty. Comma-separated
`FORCE_STAGE` values become repeated `--force-stage` options.

The three profile-specific Make examples are:

```bash
make test-dry-run PROFILE=full READS_DIR=/data/raw_reads \
  KRAKEN2_DB=/db/kraken2 \
  CHECKM2_DB=/db/checkm2/uniref100.KO.1.dmnd \
  GTDBTK_DATA_PATH=/db/gtdbtk

make test-dry-run PROFILE=taxonomy CLEAN_READS_DIR=/data/clean_reads \
  KRAKEN2_DB=/db/kraken2

make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads \
  CHECKM2_DB=/db/checkm2/uniref100.KO.1.dmnd \
  GTDBTK_DATA_PATH=/db/gtdbtk

make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads \
  RUN_MAPPING=false BAM_DIR=/data/bams BINNERS=comebin

make test-dry-run \
  BINNERS=comebin \
  COMEBIN_DEVICE=cpu

make test-dry-run \
  BINNERS=comebin \
  COMEBIN_DEVICE=cuda:0

make test-dry-run \
  PROFILE=assembly \
  MIN_COMPLETENESS=80 \
  MAX_CONTAMINATION=5

make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads \
  RUN_MAPPING=false JGI_DEPTH_FILE=/data/coverage/depth.txt BINNERS=vamb

make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads \
  RUN_MAPPING=false JGI_DEPTH_FILE=/data/coverage/depth.txt \
  BINNERS=vamb FORCE_STAGE=vamb

make test-dry-run PROFILE=assembly CLEAN_READS_DIR=/data/clean_reads \
  RUN_MAPPING=false BAM_MANIFEST=/data/bams.tsv BINNERS=comebin \
  RESTART_FROM=comebin

make test-dry-run PROFILE=assembly RUN_MEGAHIT=no \
  ASSEMBLY_IN=/data/assembly/final.contigs.fa \
  CLEAN_READS_DIR=/data/clean_reads \
  CHECKM2_DB=/db/checkm2/uniref100.KO.1.dmnd \
  GTDBTK_DATA_PATH=/db/gtdbtk
```

For an external assembly, `RUN_MEGAHIT=no` explicitly passes `--no-megahit`;
setting `ASSEMBLY_IN` alone does not disable MEGAHIT.

Dry-run examples for binner selection:

```bash
make test-dry-run PROFILE=assembly BINNERS=vamb CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=assembly BINNERS=metabinner CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=assembly BINNERS=comebin CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=assembly BINNERS=vamb,metabinner \
  REFINEMENT_MODE=auto BINS_SOURCE=auto CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=assembly BINNERS=vamb,metabinner,comebin \
  REFINEMENT_MODE=auto BINS_SOURCE=auto CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=assembly BINNERS=vamb,comebin \
  REFINEMENT_MODE=never BINS_SOURCE=comebin CLEAN_READS_DIR=/data/clean_reads
make test-dry-run PROFILE=custom RUN_FASTP=no RUN_KRAKEN2=no RUN_MEGAHIT=no \
  RUN_MAPPING=no RUN_BINNING=no RUN_CHECKM2=no RUN_MAG_FILTERING=no \
  RUN_GTDBTK=no RUN_COVERM=no RUN_SUMMARY=no
```

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

Reassembly concatenates the available paired reads into files under its output
directory before running MetaWRAP `reassemble_bins`; it does not modify the
source FASTQ files. Ensure sufficient temporary and output storage.

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
│   ├── taxonomy_summary.tsv
│   ├── taxonomy_summary.meta.json
│   ├── mag_summary.tsv
│   ├── mag_summary.meta.json
│   ├── run_overview.tsv
│   └── run_overview.meta.json
├── logs/
└── .checkpoints/
```

Important outputs include:

- `02_assembly/final.contigs.fa`: filtered assembly.
- `03_kraken2/`: per-sample Kraken2 and Bracken results.
- `04_mapping/`: sorted and indexed BAM files.
- `05_binning/`: outputs from selected binners.
- `05_binning/<method>/final_bins/`: normalized bins from each selected method.
- `06_refinement/`: MetaWRAP-refined bins.
- `08_checkm2/quality_report.tsv`: MAG quality estimates.
- `09_mags/`: MAGs meeting the configured quality thresholds.
- `10_gtdbtk/`: GTDB-Tk classifications.
- `11_abundance/`: per-sample CoverM abundance tables.
- `12_summary/taxonomy_summary.tsv`: Read-based taxonomic classifications and abundances (Kraken2/Bracken).
- `12_summary/mag_summary.tsv`: MAG-based quality, taxonomy (GTDB-Tk), and sample coverage/abundance (CoverM).
- `12_summary/run_overview.tsv`: Technical execution overview, stage statuses, and environment inventories.
- `logs/`: step-specific execution logs and generated helper scripts.

## Configurable MAG filtering

MAG filtering defaults to `MIN_COMPLETENESS=70` and
`MAX_CONTAMINATION=10`. It retains a bin only when both inclusive comparisons
hold:

```text
completeness >= MIN_COMPLETENESS
contamination <= MAX_CONTAMINATION
```

Thresholds accept integers or decimal values between 0 and 100, such as `80`,
`70.0`, `82.5`, or `5.5`; values are not rounded. For example:

```bash
bash run_metagenomics_pipeline.sh --min-completeness 80 \
  --max-contamination 5
```

For the defaults, a bin with completeness `70` and contamination `10` is
retained; completeness `69.99` or contamination `10.01` is excluded.

Only MAG filtering and effective downstream consumers of the quality-filtered
MAG collection use these thresholds. CheckM2 evaluates the selected input bins
independently of the filtering thresholds. The thresholds have no effect on
taxonomy-only runs or upstream stages such as assembly, mapping, binning,
refinement, and reassembly. They are configurable pipeline selection criteria,
not a universal classification of MAG quality. Each filtering run writes
`09_mags/mag_filter_manifest.tsv`, including reported values, thresholds,
retention status, reason, and source/output paths.

Example manifest row:

```text
mag_id	completeness	contamination	min_completeness	max_contamination	retained	reason	source_path	output_path
bin_001	70	10	70	10	true	retained	/data/bins/bin_001.fa	/results/09_mags/bin_001.fa
```

## Summary Products: Taxonomy, MAG, and Run Overview

The pipeline decomposes final reporting into three conceptually and analytically independent products located in `OUTDIR/12_summary/`:

1. **`taxonomy_summary.tsv`**: Read-based taxonomic classifications and abundances from Kraken2 and Bracken.
2. **`mag_summary.tsv`**: MAG-based quality, GTDB-Tk taxonomy, and sample-level coverage/abundance from CoverM.
3. **`run_overview.tsv`**: Technical execution metadata, stage statuses, configuration, outputs, and runtime/database inventory links.

Each product includes an accompanying `.meta.json` document providing JSON-formatted provenance, timestamps, source file paths, SHA-256 fingerprints, and publication status.

### Fundamental Architectural Guarantees

To ensure scientific rigor and avoid cross-contamination of analytical concepts:
- **`taxonomy_summary.tsv` strictly does NOT contain MAGs**: It only reports reads classified directly against reference k-mer databases.
- **`mag_summary.tsv` strictly does NOT contain Kraken2 or Bracken read abundances**: It only evaluates assembled, binned, and filtered genomes.
- **`run_overview.tsv` is NOT a biological summary**: It strictly contains execution parameters, stage plans, runtimes, error counters, and technical inventories.
- **CoverM measures MAG coverage and relative depth**, not read-based taxonomic abundance.
- **Bracken estimates read redistribution**, but does NOT substitute GTDB-Tk genome-based taxonomic placement.
- **GTDB-Tk classifies assembled genomes**, but does NOT provide read-based community profiles.

---

### 1. `taxonomy_summary.tsv`

#### Exact TSV Header
```text
sample_id	source	taxonomic_level	taxon_name	taxonomy_id	assigned_reads	estimated_reads	fraction_total_reads	status
```

#### Granularity & Columns
- **Granularity**: Long format; one row per `(sample_id, source, taxonomic_level, taxon_name)`.
- **`sample_id`**: Canonical sample identifier from `sample_manifest.tsv`. Stray reports not in the active manifest are ignored.
- **`source`**: `kraken2` or `bracken`.
- **`taxonomic_level`**: Standard rank codes (e.g., `U`, `R`, `D`, `P`, `C`, `O`, `F`, `G`, `S`) for Kraken2; configured level (e.g., `S`) for Bracken. Intermediate clades are preserved without artificial flattening.
- **`taxon_name`**: Taxon name reported by the tool with leading hierarchical indentation stripped. Cell values are sanitised against embedded tabs and newlines.
- **`taxonomy_id`**: NCBI taxonomy ID (integer). For unclassified reads in Kraken2, `taxonomy_id=0`.
- **`assigned_reads`**: Direct read count assigned to the taxon. For Kraken2, reports `direct_reads` (clade reads for `unclassified` and `root`). For Bracken, reports `kraken_assigned_reads`.
- **`estimated_reads`**: Re-estimated read count from Bracken (`new_est_reads`). When source is `kraken2`, reports `NA`.
- **`fraction_total_reads`**: Fraction of total reads (0.0 to 1.0).
- **`status`**: `classified`, `unclassified`, `zero-count`, or `estimated`.

#### Applicability
- Applicable only when `RUN_KRAKEN2=true`.
- If Kraken2 is enabled and Bracken is enabled, includes both sources.
- If Kraken2 is enabled and Bracken is disabled, includes only Kraken2 rows (`estimated_reads=NA`).
- If Kraken2 is disabled, `taxonomy_summary` is marked `not-applicable` (exit code 4), produces no output TSV, and does not read residual files from prior runs.

---

### 2. `mag_summary.tsv`

#### Exact TSV Header
```text
mag_id	sample_id	completeness	contamination	retained	filter_reason	gtdb_classification	gtdb_domain	gtdb_class	gtdb_order	gtdb_family	gtdb_genus	gtdb_species	classification_method	closest_reference	ani	af	coverm_metric	coverm_value	status
```

#### Granularity & Columns
- **Granularity**: Long format by MAG, sample, and CoverM metric.
  - When CoverM is **active**: One row per `mag_id + sample_id + coverm_metric`.
  - When CoverM is **disabled**: One row per `mag_id`, with `sample_id=NA`, `coverm_metric=NA`, `coverm_value=NA`.
  - Filtered-out MAGs (`retained=false`): Excluded from downstream GTDB-Tk and CoverM evaluation. They appear as a single row with quality metrics, `sample_id=NA`, `coverm_metric=NA`, `coverm_value=NA`, and status `quality-reported;excluded`.
- **`mag_id`**: Canonical genome identifier without `.fa`/`.fasta` extensions.
- **`sample_id`**: Sample identifier from `sample_manifest.tsv` or `NA` when CoverM is not active.
- **`completeness` / `contamination`**: Estimated percentage (0..100) from CheckM2.
- **`retained`**: `true` or `false` when MAG filtering is active; `NA` when MAG filtering is disabled.
- **`filter_reason`**: Exact retention or exclusion rationale from `mag_filter_manifest.tsv` (e.g. `retained`, `completeness_too_low`, `contamination_too_high`). `NA` when filtering is disabled.
- **`gtdb_*`**: Lineage ranks parsed from GTDB-Tk (`domain`, `class`, `order`, `family`, `genus`, `species`) alongside full `gtdb_classification`. `NA` when GTDB-Tk is disabled.
- **`classification_method` / `closest_reference` / `ani` / `af`**: Placement method, reference accession, ANI (0..100), and alignment fraction (0..1) from GTDB-Tk. `NA` when GTDB-Tk is disabled.
- **`coverm_metric` / `coverm_value`**: Metric name (`relative_abundance`, `mean_coverage`, `covered_fraction`, etc.) and float value. If a sample was not mapped or is missing in CoverM tables, `coverm_value=NA`.
- **`status`**: Composite status flags separated by semicolons: `quality-reported`, `retained`, `excluded`, `checkm2-only`, `classified`, `unclassified`, `coverm-reported`, `coverm-zero`, `coverm-missing`.

#### Applicability
- Applicable when at least one MAG stage is active: CheckM2, MAG filtering, GTDB-Tk, or CoverM.
- If all MAG stages are disabled, `mag_summary` is marked `not-applicable` (exit code 4) and produces no TSV.
- Never reads residual outputs from disabled MAG stages.

---

### 3. `run_overview.tsv`

#### Exact TSV Header
```text
section	key	value	source	status
```

#### Granularity & Columns
- **Granularity**: Key-value table (one row per key).
- **`section`**: Permitted sections: `run`, `inputs`, `configuration`, `stages`, `outputs`, `software`, `databases`, `warnings`, `errors`.
- **`key`**: Stable unique identifier within each section (e.g. `pipeline_version`, `threads`, `01_fastp.action`, `01_fastp.status`, `01_fastp.duration_seconds`, `software_inventory.path`).
- **`value`**: Text serialized value, strictly sanitised against embedded tabs, carriage returns, and newlines.
- **`source`**: Provenance identifier: `configuration`, `resume_plan`, `run_status`, `sample_manifest`, `checkpoint`, `software_inventory`, `database_inventory`, or `summary_stage`.
- **`status`**: `available`, `unknown`, `unavailable`, `pending`, `not-applicable`, `warning`, or `failed`.

#### Applicability & Failure Handling
- Applicable in every execution that passes preflight checks.
- On successful runs, published with `final_status=completed`.
- On pipeline failure, `run_status.tsv` is updated directly and idempotently with `final_status=failed`. The fatal error handler intentionally does NOT call `generate_pipeline_summaries.py` or generate `run_overview.tsv` to avoid secondary failures, recursion, or partial overview corruption. Detailed diagnostics are recorded in `OUTDIR/logs/errors.tsv` and `OUTDIR/logs/last_error.json`.

---

### 4. Metadata JSON Files (`.meta.json`)

Each summary TSV is published atomically with a companion `.meta.json` file:
- `taxonomy_summary.meta.json`
- `mag_summary.meta.json`
- `run_overview.meta.json`

Common metadata fields include:
```json
{
  "schema_version": "1",
  "summary_type": "taxonomy|mag|run-overview",
  "created_at": "2026-10-31T15:30:00Z",
  "pipeline_version": "2.0.0",
  "row_count": 42,
  "source_files": ["..."],
  "source_fingerprints": {"file.tsv": "sha256..."},
  "summary_fingerprint": "sha256...",
  "status": "available|not-applicable|failed"
}
```

---

### 5. Independent Checkpoints and Selective Invalidation

The summary stage is divided into three distinct logical checkpoints under `OUTDIR/.checkpoints/`:
- `12a_taxonomy_summary.done`
- `12b_mag_summary.done`
- `12c_run_overview.done`

#### Dependency Graph and Invalidation Isolation
- **Kraken2 / Bracken changes**: Invalidate `12a_taxonomy_summary` and `12c_run_overview`. They **never** invalidate `12b_mag_summary`.
- **CheckM2 / Filter Thresholds / GTDB-Tk / CoverM changes**: Invalidate `12b_mag_summary` and `12c_run_overview`. They **never** invalidate `12a_taxonomy_summary`.
- **Inventory changes**: Only invalidate `12c_run_overview`.
- **Independent Execution**:
  - `--force-stage taxonomy-summary` reruns `12a_taxonomy_summary` and updates `12c_run_overview`, leaving `12b_mag_summary` intact.
  - `--force-stage mag-summary` reruns `12b_mag_summary` and updates `12c_run_overview`, leaving `12a_taxonomy_summary` intact.
  - `--force-stage run-overview` regenerates `12c_run_overview` without recalculating biological tables.
  - Requesting the legacy name `--force-stage summary` fails immediately with an informative error prompting the user to specify `taxonomy-summary`, `mag-summary`, or `run-overview`.

#### Legacy Migration Policy
When the pipeline detects an older monolithic checkpoint `12_summary.done` or `12_summary.config`:
- The legacy checkpoint is **not** deleted during preflight.
- It is renamed/marked as `12_summary.legacy`.
- It is **never** reused as compatible with `12a`, `12b`, or `12c`.
- The new summaries are constructed cleanly and independently.

---

### 6. Atomic Publishing and Security Contracts

- **No Partial Files**: Summary files are written to temporary staging folders (`.build.<summary>.<pid>/`), validated against expected column counts, and atomically renamed into `OUTDIR/12_summary/`.
- **Failure Safety**: If summary generation fails or aborts, previous valid versions of TSVs and metadata remain intact.
- **TSV Sanitization**: All tab, carriage return, and newline characters within metadata or taxonomic names are converted to single spaces.
- **Dry-Run Protection**: Invoking `--dry-run` performs validation and argument resolution without creating staging directories or writing files.

---

## Execution Tracking & Benchmark Manifests

To support reproducible benchmarking, granular performance tracking, and cross-run comparisons without interfering with analytical checkpoints, graph execution, or output validity, the pipeline creates two centralized, versioned manifest tables under `OUTDIR/logs/`:

1. `OUTDIR/logs/run_manifest.tsv` (1 row per `RUN_ID`)
2. `OUTDIR/logs/stage_status.tsv` (1 row per stage and `RUN_ID`)

### Security and Governance Contracts

- **Central Execution Tracking**: These manifests record the execution lifecycle, environment metadata, and stage metrics.
- **Not Stage Outputs**: They are not products of any analytical step.
- **Not Checkpoints**: They do not participate in checkpoint dependency graphs or compatibility decisions.
- **Immune to Invalidation**: They cannot be targets of `--force-stage` or `--restart-from`.
- **Protected Against Deletion**: `clean_stage_outputs` and selective clean routines are strictly barred from removing `run_manifest.tsv`, `stage_status.tsv`, or `run_status.tsv`.
- **Atomic Publishing**: Updates use atomic temporary writes (`.tmp_tsv_*`) and `os.replace` to prevent corrupted reads during concurrent inspections or abrupt interrupts.

---

### 1. `run_manifest.tsv`

Maintains execution-level identity, configurations, input fingerprints, allocated resources, final status, and references to all other manifests. Consecutive runs in the same `OUTDIR` append distinct rows indexed by their unique `RUN_ID`.

#### Exact TSV Header (43 Columns)
```text
schema_version	run_id	pipeline_version	profile	execution_mode	dry_run	start_time_utc	end_time_utc	wall_time_seconds	final_status	exit_code	outdir	command_line	hostname	cpu_requested	threads_effective	memory_limit_bytes	active_sample_count	active_reads_source	sample_manifest_path	sample_manifest_fingerprint	assembly_source	assembly_path	assembly_fingerprint	bins_source_requested	bins_source_resolved	bins_fingerprint	min_completeness	max_contamination	tmpdir_source	tmpdir_base	keep_temporary	resume_plan_path	stage_status_path	software_inventory_path	database_inventory_path	error_history_path	last_error_path	taxonomy_summary_path	mag_summary_path	run_overview_path	warning_count	error_count
```

#### Field Contracts & Permitted Values
- **`schema_version`**: Initial version `1`.
- **`run_id`**: Canonical `YYYYMMDDTHHMMSSZ_<PID>_<RAND>` execution identifier.
- **`execution_mode`**: `apptainer` or `mamba`.
- **`dry_run`**: `true` or `false`.
- **`final_status`**: `running`, `success`, `failed`, `interrupted`, `terminated`, `dry-run-success`, `dry-run-failed`.
- **`wall_time_seconds`**: Monotonic total execution duration formatted to 3 decimal places (or `NA` while running).
- **`command_line`**: Invocation command sanitized against tokens, passwords, secrets, credentials, and API keys.
- **`cpu_requested`**: Externally allocated CPU count from `SLURM_CPUS_PER_TASK` or `PBS_NP` (or `NA` if undetermined). Total physical host CPU is never substituted.
- **`memory_limit_bytes`**: Memory allocation limit from `SLURM_MEM_PER_NODE` or cgroup quotas (or `NA` if undetermined). Total host RAM is never substituted.
- **`warning_count` / `error_count`**: Derived strictly from structured records (`errors.tsv` row counts, validation checks), never by counting arbitrary regex matches.

---

### 2. `stage_status.tsv`

Canonical stage-level execution manifest. Records applicability, planned actions, observed checkpoint decisions, wall-clock duration (via monotonic clock), CPU and memory metrics (via GNU `/usr/bin/time -v`), exit codes, and output fingerprints.

#### Exact TSV Header (32 Columns)
```text
schema_version	run_id	stage_id	stage_name	applicability	planned_action	checkpoint_action	final_status	start_time_utc	end_time_utc	wall_time_seconds	user_cpu_seconds	system_cpu_seconds	cpu_utilization_percent	max_rss_kib	io_read_bytes	io_write_bytes	exit_code	checkpoint_path	checkpoint_status	checkpoint_reason	input_state	validation_status	primary_log	primary_output	output_fingerprint	runtime	tool_version	runtime_fingerprint	database_fingerprint	error_id	message
```

#### Stage Lifecycle & Metrics Separation Matrix

The pipeline strictly isolates planned decisions from observed runtime actions. Unexecuted stages never spoof zero-duration or zero-memory metrics:

| Field | Executed Stage (`run`, `force`, `restart`) | Reused Checkpoint (`reuse`) | Disabled Stage | Not-Applicable Stage | Unreached (After Upstream Failure) |
|---|---|---|---|---|---|
| **`applicability`** | `applicable` | `applicable` | `disabled` | `not-applicable` | `applicable` |
| **`planned_action`** | `run`, `force`, `restart` | `reuse` | `disabled` | `not-applicable` | `run`, `force`, `restart` |
| **`checkpoint_action`** | `executed`, `forced`, `restarted` | `reused` | `skipped-disabled` | `skipped-not-applicable` | `not-reached` |
| **`final_status`** | `success` / `failed` | `reused` | `disabled` | `not-applicable` | `not-reached` |
| **`start_time_utc`** | ISO-8601 UTC | `NA` | `NA` | `NA` | `NA` |
| **`end_time_utc`** | ISO-8601 UTC | `NA` | `NA` | `NA` | `NA` |
| **`wall_time_seconds`** | Measured monotonic seconds | `NA` | `NA` | `NA` | `NA` |
| **`user_cpu_seconds`** | Measured via GNU time | `NA` | `NA` | `NA` | `NA` |
| **`system_cpu_seconds`**| Measured via GNU time | `NA` | `NA` | `NA` | `NA` |
| **`max_rss_kib`** | Peak memory from GNU time | `NA` | `NA` | `NA` | `NA` |
| **`cpu_utilization_percent`** | Calculated % from GNU time | `NA` | `NA` | `NA` | `NA` |
| **`io_read_bytes`** | `NA` (block count not bytes) | `NA` | `NA` | `NA` | `NA` |
| **`io_write_bytes`** | `NA` (block count not bytes) | `NA` | `NA` | `NA` | `NA` |
| **`exit_code`** | `0` or command error code | `0` | `NA` | `NA` | `NA` |

---

### 3. Relationship with Other Pipeline Manifests

```
                      ┌────────────────────────────┐
                      │    resume_plan.tsv         │ (Pre-execution frozen plan)
                      └─────────────┬──────────────┘
                                    │ initializes
                                    ▼
                      ┌────────────────────────────┐
                      │    stage_status.tsv        │ (Canonical runtime stage state)
                      └─────────────┬──────────────┘
                                    │ generates
                                    ▼
                      ┌────────────────────────────┐
                      │    run_status.tsv          │ (Compatibility view)
                      └────────────────────────────┘
```

- **`resume_plan.tsv` (`OUTDIR/resume_plan.tsv`)**: Frozen pre-execution plan recording planned actions, graph parents, dependencies, and reasons. Remains immutable throughout execution.
- **`stage_status.tsv` (`OUTDIR/logs/stage_status.tsv`)**: The single canonical source of truth for observed stage state, durations, and resource metrics.
- **`run_status.tsv` (`OUTDIR/run_status.tsv`)**: Maintained as an automatic, idempotent compatibility projection generated directly from `stage_status.tsv`.
- **`run_overview.tsv` (`OUTDIR/12_summary/run_overview.tsv`)**: Summary table that references and consumes `run_manifest.tsv` and `stage_status.tsv` without independently recalculating execution parameters or stage metrics.

---

## Checkpoints and reproducibility

After successful completion, each stage writes a marker under `.checkpoints/`.
Read-consuming stages record the active read source and absolute directory,
suffixes, sample IDs, and each read's absolute path, size, and modification
time. This reproducible fingerprint avoids hashing large FASTQ files on every
run.

Binner checkpoints (`05a_vamb`, `05b_metabinner`, `05c_comebin`) are independent
and are consulted only for selected methods. Refinement (`06_refinement`)
depends on the selected input collections, while CheckM2 (`08_checkm2`) depends
on the effective bins collection. Changing the selected binners does not by
itself invalidate an unchanged binner's checkpoint, but changing refinement
inputs or the effective bin collection invalidates the appropriate downstream
MAG checkpoints. Unselected binner checkpoints and outputs are left untouched
and ignored. `--binners` configures the pipeline; it is not a comparison or
benchmark facility.

The MAG-filtering checkpoint records both thresholds, inclusive comparison
operators, the CheckM2 report fingerprint, input-bin fingerprint, retained-MAG
fingerprint, and filter-manifest fingerprint. Changing either threshold
repeats MAG filtering and its effective consumers (GTDB-Tk, CoverM, and/or
summary) while reusing a valid CheckM2 checkpoint and preserving upstream
stages. In taxonomy-only runs the threshold values are validated and reported
as unused; they do not activate MAG filtering or affect Kraken2/Bracken.

The Vamb checkpoint records its effective version and image URI/fingerprint,
assembly and generated abundance-table fingerprints, ordered abundance header
and sample names, binsplitting state/separator, `VAMB_MINFASTA` meaning and
value, other effective Vamb parameters, and the normalized final-bin path.
Changing any Vamb input or option rebuilds Vamb and only its effective
descendants (for example refinement and CheckM2 when they consume Vamb bins);
it does not invalidate Kraken2, MEGAHIT, mapping, or independent binners.
Failed or forced Vamb runs clean only Vamb-controlled outputs under `OUTDIR`.
Abundance inputs such as depth tables and assemblies are retained, partial
Vamb bins are not checkpointed, and downstream stages cannot consume them.

Kraken2 and MEGAHIT have independent checkpoints and dependencies: Kraken2
depends on the active reads, Kraken2/Bracken database, parameters, and image,
not on MEGAHIT or an assembly. MEGAHIT depends on the active reads and assembly
parameters, not on Kraken2. Invalidating Kraken2 does not invalidate the
assembly branch, and vice versa. A reads-source change invalidates every
enabled consumer of reads. Disabled stages do not inspect, write, or remove
their checkpoint.

When using `--assembly-in` with MEGAHIT disabled, assembly-dependent stage
manifests record the normalized FASTA path, file size, and modification time.
Changing that fingerprint invalidates the dependent assembly branch, never
Kraken2/Bracken. The disabled MEGAHIT stage does not inspect or write the
`02_megahit` checkpoint.

### Resume, force, and restart

The default is `--resume`: compatible checkpoints and valid products are reused;
missing or incompatible checkpoints, changed inputs, or invalid products cause
the stage to run again. Resume controls are mutually exclusive:

| Option | Effect |
|---|---|
| `--resume` | Explicitly select the default selective reuse behavior |
| `--force` | Force every active stage; each stage cleans its own outputs immediately before running |
| `--force-stage STAGE` | Repeatable; force each requested root and run only its effective descendants |
| `--restart-from STAGE` | Destructively clear the root and all effective descendants before running the chain |

Public stage names are `fastp`, `kraken2`, `megahit`, `mapping`, `vamb`,
`metabinner`, `comebin`, `refinement`, `reassembly`, `checkm2`,
`mag-filtering`, `gtdbtk`, `coverm`, and `summary`. A force-stage root is
cleaned just before execution; descendants clean themselves immediately before
their own execution. `--restart-from` is deliberately different: it
pre-invalidates checkpoints, manifests, and controlled outputs throughout the
descendant chain. Neither operation modifies independent branches.

Examples:

```bash
# Reuse everything still valid (the default)
bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads

# Recompute Kraken2 and the stages that consume its outputs; MEGAHIT is independent
bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads \
  --force-stage kraken2

# Force multiple roots in one graph traversal
bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads \
  --force-stage kraken2 --force-stage checkm2

# Restart mapping and its effective consumers only, using external coverage
bash run_metagenomics_pipeline.sh --profile assembly --no-mapping \
  --jgi-depth-file /data/coverage/depth.txt --binners vamb \
  --restart-from vamb
```

`resume_plan.tsv` records the frozen plan; `run_status.tsv` separately records
the real outcome (`completed`, `reused`, `failed`, `interrupted`, `dry-run`,
`not-reached`, `disabled`, or `not-applicable`). External reads, assemblies,
BAMs/indexes, depth tables, and databases use absolute-path/size/mtime
fingerprints rather than full FASTQ hashes. If an external input changes after
the plan is frozen, the run aborts instead of silently changing its plan.

Cleanup is restricted to stage-owned paths inside `OUTDIR`; a path resolving
outside it is rejected. External inputs, databases, containers, and `OUTDIR`
itself are never removed. Dry-run validates and prints the plan and prospective
cleanup paths without writing result, checkpoint, plan, or status files.

The controlled output roots are:

| Stage | Controlled outputs cleaned on rerun |
|---|---|
| fastp | `01_qc/` |
| kraken2 | `03_kraken2/` |
| megahit | `02_assembly/` |
| mapping | `04_mapping/`, `05_binning/depth_matrix.txt` |
| Vamb / MetaBinner / COMEBin | The selected method's `05_binning/<method>/` |
| refinement | `06_refinement/` |
| reassembly | `07_reassembly/` |
| CheckM2 | `08_checkm2/` |
| MAG filtering | `09_mags/` |
| GTDB-Tk | `10_gtdbtk/` |
| CoverM | `11_abundance/` |
| taxonomy-summary | `12_summary/taxonomy_summary.tsv`, `12_summary/taxonomy_summary.meta.json` |
| mag-summary | `12_summary/mag_summary.tsv`, `12_summary/mag_summary.meta.json` |
| run-overview | `12_summary/run_overview.tsv`, `12_summary/run_overview.meta.json` |

### Runtime and database identity inventories

Each planned run builds a coordinated inventory generation under
`logs/inventory/current/`. The public paths remain
`logs/software_versions.tsv` and `logs/database_versions.tsv`; they are
relative symbolic links to the current generation. A generation is validated
before publication, and the previous generation is retained until both new
TSVs and public links are checked. Inventory files are logs/traceability data,
not analytical outputs, checkpoints, or inputs to compatibility decisions.
The stage `.config` manifest remains the primary compatibility record.

The software inventory has this exact header and column order:

```text
stage_id	tool	role	version	version_source	execution_mode	runtime_name	runtime_uri	runtime_path	runtime_size_bytes	runtime_mtime_ns	runtime_fingerprint_method	runtime_fingerprint	command_used	status
```

The database inventory has this exact header and column order:

```text
stage_id	database	role	path	release_or_version	version_source	size_bytes	mtime_ns	fingerprint_method	fingerprint	status
```

Example rows (illustrative values):

```text
stage_id	tool	role	version	version_source	execution_mode	runtime_name	runtime_uri	runtime_path	runtime_size_bytes	runtime_mtime_ns	runtime_fingerprint_method	runtime_fingerprint	command_used	status
03_kraken2	kraken2	current-runtime	2.1.3	runtime-command	apptainer	kraken2	docker://staphb/kraken2:2.1.3	/opt/containers/kraken2.sif	123456789	1793415102948000000	sha256-file	0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef	kraken2 --version	available
```

```text
stage_id	database	role	path	release_or_version	version_source	size_bytes	mtime_ns	fingerprint_method	fingerprint	status
08_checkm2	checkm2	current-runtime	/databases/checkm2/uniref100.KO.1.dmnd	unknown	unknown	987654321	1793415102948000000	sha256-file	abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef	available
```

Software rows are emitted for selected active stages only. A `producer` row is
read from the historical checkpoint manifest; a `current-runtime` row records
the identity resolved for this invocation. A reused stage can therefore have
both rows. Database rows use the same role distinction. Valid role values are
`producer` and `current-runtime`; status values are `available`, `unknown`,
`unavailable`, and `reused`. `unknown` means the resource was present but some
identity could not be established; it is not the same as `unavailable`, and
`unknown` never compares equal as evidence of compatibility. A producer whose
old checkpoint has no identity metadata is recorded as `unknown`/`reused` and
is eligible for the non-destructive migration policy.

Version sources are `runtime-command`, `checkpoint-manifest`, `image-label`,
`commit`, or `unknown`. A current version query failure is not treated as an
absent tool. A tool/runtime that is missing for a stage that must execute is a
preflight error; reuse can proceed without the current runtime when the
checkpoint products remain valid. New checkpoints record the effective
runtime version, execution mode, runtime name/URI/path/size/mtime/fingerprint,
and database release/path/fingerprint fields where the stage consumes a
database. They do not record database fields in unrelated stages, nor
Apptainer-only values as if they applied to Mamba.

#### Current software queries

The executable is queried in the selected runtime. The command column contains
only the version query, not the analytical command.

| Tool | Query |
|---|---|
| fastp | `fastp --version` |
| MEGAHIT | `megahit --version` |
| Kraken2 | `kraken2 --version` |
| Bracken | `bracken -v` |
| Bowtie2 | `bowtie2 --version` |
| SAMtools | `samtools --version` |
| Vamb | `vamb --version` |
| MetaBinner | `metabinner --version`, then package metadata when available |
| COMEBin | `run_comebin.sh -V`, then package metadata when available |
| MetaWRAP | `metawrap --version` |
| CheckM2 | `checkm2 --version` |
| GTDB-Tk | `gtdbtk --version` |
| CoverM | `coverm --version` |
| Rscript summary runtime | `Rscript --version` |

When an Apptainer command does not report a version, an explicit image version
or revision label is used when present. Otherwise the version is `unknown`.
Mamba identity is calculated from the exact selected environment's
`conda-meta/*.json`: package name, version, build, and channel are sorted in
binary order and hashed as `conda-meta-v1`. Absolute paths, timestamps,
command history, credentials, and full environment exports are excluded.
`runtime_path` is the environment prefix when it can be resolved; size and
mtime are `not-applicable` for Mamba.

Apptainer images are resolved to absolute paths and identified with the
lowercase SHA-256 of the SIF (`sha256-file`). Size and nanosecond mtime are
also recorded for traceability. The SIF hash and version query are cached for
the invocation. Moving an unchanged SIF changes its recorded path but not its
functional identity. Once the plan is frozen, changes to its resolved path,
size, mtime, or filesystem object abort execution without replanning.

#### Database identities

`CHECKM2_DB` is a file and receives a full SHA-256 (`sha256-file`), size, and
nanosecond mtime. The release is taken only from trusted metadata when
available; it is otherwise `unknown`. A known file fingerprint is sufficient
to identify the database even when its release is unknown.

Kraken2/Bracken and GTDB-Tk databases are directories. They are never hashed
recursively. Their `metadata-keyfiles-v1` fingerprint is SHA-256 over sorted
`relative_path<TAB>size_bytes<TAB>mtime_ns` records with a final newline.
Only these known keyfiles are considered:

| Database | Metadata/keyfiles considered |
|---|---|
| Kraken2/Bracken | `hash.k2d`, `opts.k2d`, `taxo.k2d`, `seqid2taxid.map`, `taxonomy/nodes.dmp`, `taxonomy/names.dmp`, `taxonomy/merged.dmp`, `taxonomy/delnodes.dmp`, `database.kraken`, and top-level `*.kmer_distrib` / `*.kraken` files |
| GTDB-Tk | `VERSION`, `release.txt`, `gtdbtk.package_version`, `metadata/metadata.tsv`, `metadata/metadata.txt`, `metadata/genome_metadata.tsv`, `metadata/taxonomy.tsv`, and `taxonomy/gtdb_taxonomy.tsv` |

Missing optional keyfiles are omitted. If none of the known files is present,
the fingerprint and version remain `unknown`; the pipeline does not walk the
rest of a large database as a fallback and does not download metadata. For a
directory, `size_bytes` is the sum of only the included keyfiles and `mtime_ns`
is the maximum keyfile mtime; neither describes the total database size.
Database location is recorded for access and traceability but is not part of
functional identity.

#### Functional identity versus traceability

| Functional compatibility fields | Traceability/access fields |
|---|---|
| Runtime: execution mode, known runtime fingerprint method/hash, and known effective tool version | Runtime: logical name, configured URI, resolved path, size, mtime, and version-query command |
| Database: known fingerprint method/hash and known release/version | Database: absolute path, size, and mtime |

When a known fingerprint is available it is compared first; a known version or
release is also compared. Paths, URIs, size, and mtime do not invalidate
analytically identical resources. A missing historical identity follows the
legacy migration policy; a newly written checkpoint must have a known runtime
version or runtime fingerprint. An `unknown` value is never used as a matching
identity.

The order is: resolve active/applicable stages and current runtimes/databases;
calculate their identities; load producer identities from existing checkpoint
manifests; compare each stage's checkpoint; freeze the plan; build and publish
the inventory generation; then execute. Inventory TSVs are not read to make
reuse decisions. A changed functional identity rebuilds only its consuming
stage and effective descendants in the active dependency graph. For example,
Kraken2/Bracken changes do not rebuild the MAG branch, and CoverM changes do
not invalidate GTDB-Tk.

Dry-run performs only read-only resolution and queries, constructs inventory
rows in memory, and reports the paths and entry/unknown/unavailable counts it
would publish. It does not create the inventory directory, TSVs, checkpoints,
or analytical outputs. The inventories are protected from stage cleanup,
`--force-stage`, `--restart-from`, and selective invalidation.

Limitations: version options and labels vary between tool/image releases, so
some tools can legitimately be recorded as `unknown`; the keyfile sets are
deliberately conservative and do not prove byte-for-byte identity of every
large database file. Mamba host tools run outside a named environment use the
active Conda prefix when available; without a prefix, their version may be
known while their environment fingerprint is `unknown`.

## Pre-flight validation

Before running the workflow, the script checks only resources needed by active
stages, including:

- selected profile, resulting stage flags, active branch, and reads source;
- existence and pairing of reads when an active stage needs them;
- presence of at least one `*_R1.fastq.gz` file;
- availability of the matching R2 file for every R1 file;
- required database paths according to enabled stages;
- the selected Apptainer/Singularity or Mamba runtime;
- required host binaries in Mamba mode;
- MetaBinner script configuration;
- assembly availability only when an active assembly consumer needs one.

## Fatal Error Capture and Fault Recovery Architecture

The pipeline implements a structured, robust, and non-recursive fatal error capture system designed to ensure post-mortem reproducibility, process safety, and clear recovery paths for large-scale metagenomic workflows.

### 1. Dual-Phase Trap Architecture

The error handling system operates in two distinct phases:

```text
┌────────────────────────────────────────────────────────┐
│ Phase 1: Early Error Traps (Bootstrap & Preflight)     │
│ - Installed immediately after shebang (set -Eeuo pipefail)│
│ - Traps ERR, INT (130), TERM (143), HUP (129)          │
│ - Independent of OUTDIR, LOGS_DIR, Python, or disk writes│
│ - Single formatted error message to stderr             │
│ - Preserves original exit code without side effects    │
└───────────────────────────┬────────────────────────────┘
                            │ Argument parsing & directory layout complete
                            ▼
┌────────────────────────────────────────────────────────┐
│ Phase 2: Full Persistent Error System                  │
│ - Installed via install_pipeline_traps                 │
│ - Replaces early traps with canonical handlers:        │
│     * pipeline_error_handler                           │
│     * pipeline_signal_handler                          │
│     * pipeline_exit_handler                            │
│ - Process-safe deduplication via atomic file markers   │
│ - Thread/process-safe locked persistence in errors.tsv │
│ - Atomic materialized view in last_error.json          │
│ - Per-run checkpoint tracking in checkpoint_state/     │
│ - Idempotent run_status.tsv updates                    │
└────────────────────────────────────────────────────────┘
```

#### Early Traps vs Full System
| Property | Early Traps | Full Persistent Traps |
|---|---|---|
| **Activation Point** | Immediately after shebang (`set -Eeuo pipefail`) | After argument parsing and `OUTDIR/logs/` initialization |
| **Handlers** | `early_error_handler`, `early_signal_handler` | `pipeline_error_handler`, `pipeline_signal_handler`, `pipeline_exit_handler` |
| **Disk Writes** | None (read-only / stderr only) | `errors.tsv`, `last_error.json`, `run_status.tsv`, `checkpoint_state/` |
| **ERROR_ID** | Not assigned | Sequential `<RUN_ID>:E0001` assigned under `flock` |
| **Exit Code** | Strictly preserved | Strictly preserved |
| **Log Snippet** | None | Bounded log tail (max 50 lines, max 8 KiB) |
| **Scope** | Command-line errors, path resolution, early preflight | All pipeline stages, tool execution, and child workflows |

#### What Bash Traps Can and Cannot Capture
- **Captured**:
  - Unhandled non-zero exit codes from commands executed in Bash functions (`ERR`).
  - Process interruptions and terminations (`SIGINT` = 130, `SIGTERM` = 143, `SIGHUP` = 129).
  - Wrapper errors in command execution, container execution, and piped commands (`run_cmd_with_tee` preserves `PIPESTATUS[0]`).
  - Unexpected script exits with non-zero exit code (`EXIT` fallback).
- **Not Captured (Explicit Limitations)**:
  - Bash syntax errors occurring during initial script parsing (e.g. mismatched quotes or syntax errors before execution starts). These must be validated prior to execution using `bash -n run_metagenomics_pipeline.sh` or `shellcheck`.
  - Errors occurring before trap installation.
  - Certain fatal bash parameter expansions occurring during shell initialization.
  - Failures inside asynchronous background tasks spawned with `&` unless explicitly waited on.
  - Subprocesses spawned without trap propagation (child processes do not inherit bash error traps unless subshells execute within parent trap scopes).

### 2. Execution Identity and Context Variables

Every pipeline execution receives a unique, collision-resistant identifier:
```text
RUN_ID = YYYYMMDDTHHMMSSZ_<PID>_<RANDOM_SHORT>
Example: 20261031T212500Z_18432_a91f
```

Along with `RUN_ID`, execution provenance is tracked via:
- `PIPELINE_VERSION`: Semantic version string.
- `PIPELINE_START_TIME_UTC`: ISO 8601 UTC execution start timestamp.
- `PIPELINE_PID`: Process identifier of the root pipeline process.
- `PIPELINE_HOSTNAME`: Execution host system name.

During execution, wrappers maintain live stage context:
- `CURRENT_STAGE`: Canonical stage identifier (e.g. `02_megahit`, `05a_vamb`, `08_checkm2`).
- `CURRENT_STAGE_NAME`: Human-readable stage display name.
- `CURRENT_STAGE_ACTION`: Planned action (`run`, `reuse`, `force`).
- `CURRENT_STAGE_LOG`: Primary log file path associated with the stage.
- `CURRENT_STAGE_RUNTIME`: Active runtime environment (`apptainer`, `mamba`, or `host`).
- `CURRENT_STAGE_COMMAND`: Sanitized command string currently executing.
- `CURRENT_COMMAND_SOURCE_FILE`: File from which the execution wrapper was invoked.
- `CURRENT_COMMAND_SOURCE_LINE`: Line number in the calling script where execution was requested.
- `LAST_VALID_CHECKPOINT`: Name of the last successfully validated checkpoint in *this execution*.
- `LAST_VALID_CHECKPOINT_PATH`: Absolute path to the validated `.done` checkpoint marker.
- `LAST_VALID_PRODUCT`: Path to the primary output product validated in the last successful stage.
- `LAST_VALID_PRODUCT_FINGERPRINT`: SHA-256 / size / mtime fingerprint of the validated product.

### 3. Persistent Error Artifacts

All persistent diagnostic records are stored in `OUTDIR/logs/`:

| Artifact Path | Format | Role & Concurrency Semantics |
|---|---|---|
| `OUTDIR/logs/errors.tsv` | TSV (22 columns) | **Primary immutable history**. Append-only across pipeline runs. Writes are serialized using exclusive kernel advisory locks (`fcntl.flock`). |
| `OUTDIR/logs/last_error.json` | JSON | **Materialized view** of the most recent fatal termination. Written atomically using temporary file creation and POSIX rename (`os.replace`). |
| `OUTDIR/logs/checkpoint_state/<RUN_ID>.tsv` | TSV | **Execution-scoped checkpoint audit**. Contains only checkpoints validated during the specific `RUN_ID`, preventing state contamination across runs. |
| `OUTDIR/logs/last_valid_checkpoint.state` | Key-Value State | Atomic pointer recording the overall latest validated checkpoint and product fingerprint. |
| `OUTDIR/logs/.error_recorded.<RUN_ID>` | Marker File | **Deduplication guard**. Tracks `pending` and `committed` states to eliminate redundant writes across subshells, wrappers, and exit traps. |

#### Protection During Cleaning and Retries
None of the error logging artifacts (`errors.tsv`, `last_error.json`, `checkpoint_state/`, `last_valid_checkpoint.state`, `.error_recorded.*`) are analytical outputs or checkpoints. They are explicitly excluded from deletion by `clean_stage_outputs`, `--force-stage`, and `--restart-from`. Rerunning a failed pipeline retains the full historical log of previous errors.

### 4. Deduplication and Concurrency Safety

Metagenomic pipelines frequently run commands across subshells, background wrappers, and complex control structures. Without coordination, an error in a subshell can trigger the subshell's `ERR` trap, return non-zero to the parent wrapper, trigger the parent's `ERR` trap, and subsequently trigger the `EXIT` trap, resulting in four identical error entries and noisy logs.

To prevent duplication:
1. **In-Memory Guards**: The flags `ERROR_HANDLER_ACTIVE`, `ERROR_ALREADY_RECORDED`, and `PIPELINE_TERMINATION_RECORDED` prevent nested error handler recursion within the same shell process.
2. **Atomic File Marker State Machine**:
   - `scripts/write_error_record.py reserve` creates `OUTDIR/logs/.error_recorded.<RUN_ID>` with state `pending` and the active PID using atomic file creation (`O_CREAT | O_EXCL`).
   - If a marker already exists with `pending` or `committed`, subsequent handlers exit immediately with status 0 without writing duplicate records.
   - **Dead PID Detection & Recovery**: If an execution crashed or was killed (`SIGKILL`) while in `pending` state, subsequent invocations check whether the recorded PID is still alive using `os.kill(pid, 0)`. If the PID is dead, the stale pending marker is automatically cleaned and reclaimed.
3. **Single Persistent Record**: Only one record is committed per logical failure.
4. **Single Stderr Message**: Exactly one `[FATAL]` message block is printed to stderr.

### 5. Sanitization, Secret Redaction, and Log Truncation

Error records sanitize all sensitive and control data before writing:
- **Secret Redaction**: Passwords, API tokens, authorization headers, and credentials matching patterns (`password=...`, `token=...`, `api_key=...`, `http://user:pass@...`) are automatically redacted to `[REDACTED]`.
- **TSV Field Sanitization**: Tabs (`\t`), newlines (`\n`), and carriage returns (`\r`) in command strings or technical messages are replaced with single spaces to protect TSV table formatting.
- **Log Snippet Limits**: The fatal error recorder extracts the tail of the failed stage's primary log:
  - Maximum lines: 50 lines.
  - Maximum size: 8 KiB (8,192 bytes).
  - Binary file detection: If a log contains null bytes (`\x00`), log extraction is replaced with `[Binary log file suppressed]`.

### 6. Primary Log Mapping and Recovery Hints

Each stage is associated with a primary log file and a canonical recovery action:

| Stage ID | Stage Display Name | Primary Log Path | Canonical Recovery Action |
|---|---|---|---|
| `01_fastp` | QC / Adapter Trimming (fastp) | `OUTDIR/logs/01_fastp.log` | Check read files, permissions, and disk space; retry with `--force-stage fastp` |
| `02_megahit` | De Novo Assembly (MEGAHIT) | `OUTDIR/logs/02_megahit.log` | Check memory limits and CPU threads; rerun with `--force-stage megahit` |
| `03_kraken2` | Taxonomic Classification (Kraken2/Bracken) | `OUTDIR/logs/03_kraken2.log` | Verify Kraken2 database integrity; retry with `--force-stage kraken2` |
| `04_mapping` | Contig Read Mapping & Coverage (Bowtie2/SAMtools) | `OUTDIR/logs/04_mapping.log` | Check alignment memory and disk space; retry with `--force-stage mapping` |
| `05a_vamb` | Metagenomic Binning (Vamb) | `OUTDIR/logs/05a_vamb.log` | Verify contig depth matrix and Vamb runtime; retry with `--force-stage vamb` |
| `05b_metabinner` | Metagenomic Binning (MetaBinner) | `OUTDIR/logs/05b_metabinner.log` | Check MetaBinner scripts directory and dependencies; retry with `--force-stage metabinner` |
| `05c_comebin` | Metagenomic Binning (COMEBin) | `OUTDIR/logs/05c_comebin.log` | Check CUDA/CPU resources and COMEBin environment; retry with `--force-stage comebin` |
| `06_refinement` | Bin Consensus Refinement (MetaWRAP) | `OUTDIR/logs/06_refinement.log` | Inspect binner outputs and MetaWRAP environment; retry with `--force-stage refinement` |
| `07_reassembly` | MAG Reassembly (MetaWRAP) | `OUTDIR/logs/07_reassembly.log` | Check reassembly memory requirements; retry with `--force-stage reassembly` |
| `08_checkm2` | MAG Quality Assessment (CheckM2) | `OUTDIR/logs/08_checkm2.log` | Verify CheckM2 DIAMOND database path; retry with `--force-stage checkm2` |
| `09_mags` | MAG Filtering & Deduplication | `OUTDIR/logs/09_mags.log` | Verify completeness/contamination thresholds; retry with `--force-stage mag-filtering` |
| `10_gtdbtk` | Taxonomic Classification of MAGs (GTDB-Tk) | `OUTDIR/logs/10_gtdbtk.log` | Verify GTDB-Tk reference database path; retry with `--force-stage gtdbtk` |
| `11_abundance` | MAG Abundance Profiling (CoverM) | `OUTDIR/logs/11_abundance.log` | Check CoverM inputs and BAM alignments; retry with `--force-stage coverm` |
| `12a_taxonomy_summary` | Taxonomy Summary Report | `OUTDIR/logs/12_summary.log` | Check Bracken outputs; retry with `--force-stage taxonomy-summary` |
| `12b_mag_summary` | MAG Summary Report | `OUTDIR/logs/12_summary.log` | Check MAG filtering and CheckM2 outputs; retry with `--force-stage mag-summary` |
| `12c_run_overview` | Run Execution Overview | `OUTDIR/logs/12_summary.log` | Check run status logs; retry with `--force-stage run-overview` |

### 7. Signal Handling and Exit Semantics

- **Signals (`INT`, `TERM`, `HUP`)**:
  - `SIGINT` (Ctrl+C) terminates with exit code **130** and records `signal_name=INT`.
  - `SIGTERM` terminates with exit code **143** and records `signal_name=TERM`.
  - `SIGHUP` terminates with exit code **129** and records `signal_name=HUP`.
- **Exit Trap (`pipeline_exit_handler`)**:
  - Acts solely as a safety fallback. If the script exits with non-zero status and an error record was not yet committed by `ERR` or signal traps, the exit handler writes a fallback termination record.
  - If the script exits normally with status 0, the exit trap performs a clean no-op.
- **Strict Exit Code Preservation**:
  - The fatal error handler captures the exact exit code of the failing tool (e.g. 1, 2, 127, 137).
  - Wrapper utilities and helper scripts never substitute the exit code with handler exit codes.

### 8. Example Artifacts

#### Example `errors.tsv`
```text
run_id	error_id	timestamp_utc	stage_id	stage_display_name	planned_action	exit_code	signal_name	trap_source_file	trap_source_line	caller_source_file	caller_source_line	active_function	sanitized_command	primary_log	runtime_seconds	last_valid_checkpoint	last_valid_product	technical_message	recovery_hint	log_tail_snippet	status_updated
20261031T212500Z_18432_a91f	20261031T212500Z_18432_a91f:E0001	2026-10-31T21:28:45Z	02_megahit	De Novo Assembly (MEGAHIT)	run	137	NONE	run_metagenomics_pipeline.sh	3820	run_metagenomics_pipeline.sh	3815	run_megahit	megahit -1 ... -2 ... -o ...	/data/out/logs/02_megahit.log	225	01_fastp	/data/out/01_qc/trimmed/sample_01_R1.fastq.gz	Command exited with status 137	Check memory limits and CPU threads; rerun with --force-stage megahit	Error: out of memory allocating 34359738368 bytes	true
```

#### Example `last_error.json`
```json
{
  "schema_version": "1.0",
  "run_id": "20261031T212500Z_18432_a91f",
  "error_id": "20261031T212500Z_18432_a91f:E0001",
  "timestamp_utc": "2026-10-31T21:28:45Z",
  "stage_id": "02_megahit",
  "stage_display_name": "De Novo Assembly (MEGAHIT)",
  "planned_action": "run",
  "exit_code": 137,
  "signal_name": "NONE",
  "trap_source_file": "run_metagenomics_pipeline.sh",
  "trap_source_line": 3820,
  "caller_source_file": "run_metagenomics_pipeline.sh",
  "caller_source_line": 3815,
  "active_function": "run_megahit",
  "sanitized_command": "megahit -1 ... -2 ... -o ...",
  "primary_log": "/data/out/logs/02_megahit.log",
  "runtime_seconds": "225",
  "last_valid_checkpoint": "01_fastp",
  "last_valid_product": "/data/out/01_qc/trimmed/sample_01_R1.fastq.gz",
  "technical_message": "Command exited with status 137",
  "recovery_hint": "Check memory limits and CPU threads; rerun with --force-stage megahit",
  "log_tail_snippet": "Error: out of memory allocating 34359738368 bytes",
  "status_updated": true
}
```

#### Example `checkpoint_state/<RUN_ID>.tsv`
```text
run_id	stage_id	checkpoint_path	validated_product	product_fingerprint	timestamp_utc
20261031T212500Z_18432_a91f	01_fastp	/data/out/.checkpoints/01_fastp.done	/data/out/01_qc/trimmed/sample_01_R1.fastq.gz	size:1428570|mtime:1730410125	2026-10-31T21:26:15Z
```

### 9. Recovering and Resuming After a Failure

When a pipeline execution terminates with an error:
1. **Inspect the materialized error view**:
   ```bash
   cat OUTDIR/logs/last_error.json
   ```
2. **Review the technical message and recovery hint**:
   The `recovery_hint` field specifies whether the stage requires resource adjustments, database path corrections, or input verification.
3. **Resume the workflow**:
   - To resume reusing all previously validated checkpoints (the default):
     ```bash
     bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads --outdir OUTDIR
     ```
   - To rerun the failed stage specifically without invalidating unrelated branches:
     ```bash
     bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads --outdir OUTDIR --force-stage megahit
     ```
   - To clear the failed stage and all downstream consumers:
     ```bash
     bash run_metagenomics_pipeline.sh --profile full --reads-dir /data/reads --outdir OUTDIR --restart-from megahit
     ```

---

## Temporary Workspace Management

The pipeline provides secure, isolated, and reproducible management of temporary working directories for intermediate files—specifically the concatenated paired-end read files required by the MetaWRAP bin reassembly stage (`07_reassembly`).

### Command-Line Options

| Option | Argument | Description | Default |
|---|---|---|---|
| `--tmpdir` | `DIR` | Base directory for temporary pipeline workspaces. The pipeline creates an isolated `RUN_ID` directory beneath this path. | Resolved via precedence |
| `--keep-temporary` | None | Preserves the execution's temporary workspace after success or failure for manual inspection. | Disabled |
| `--no-keep-temporary` | None | Removes the execution's temporary workspace during controlled cleanup upon completion. | Enabled (default) |

### Base Resolution Precedence

The temporary base directory (`PIPELINE_TMPDIR_BASE`) is resolved during pre-flight using strict precedence:

| Priority | Tier | Source | Value Assigned | Location Rules |
|:---:|---|---|---|---|
| **1 (Highest)** | CLI | `--tmpdir DIR` | `PIPELINE_TMPDIR_SOURCE="cli"` | Canonical absolute path. **Must NOT** reside within `OUTDIR`. |
| **2** | Environment | `TMPDIR` | `PIPELINE_TMPDIR_SOURCE="environment"` | Captured from original `TMPDIR`. **Must NOT** reside within `OUTDIR`. |
| **3 (Lowest)** | Fallback | Default | `PIPELINE_TMPDIR_SOURCE="fallback"` | Exactly `${OUTDIR}/.tmp`. Permitted **only** when fallback is activated. |

> [!IMPORTANT]
> The path `${OUTDIR}/.tmp` is permitted **strictly and exclusively** when `PIPELINE_TMPDIR_SOURCE=fallback`. If a user manually supplies `--tmpdir ./metagenomics_out/.tmp` or sets `TMPDIR` to that path, the pipeline explicitly rejects it.

### Location Restrictions

The pipeline canonicalizes all candidate paths (resolving symlinks and eliminating `..`) and rejects:
- The filesystem root (`/`);
- The user's home directory (`$HOME`);
- The output directory (`OUTDIR`) itself;
- Any path situated inside `OUTDIR` (unless resolved via the fallback source);
- The pipeline repository / script root directory;
- Any path whose managed base is an existing symlink.

### Per-Run Isolated Workspace Layout

Temporary directories are never shared statically across executions. The pipeline creates an exclusive, per-execution directory tree using the current `RUN_ID`:

```text
PIPELINE_TMPDIR_BASE/
└── metagenomics_pipeline/
    └── <RUN_ID>/                                # Mode 0700, exclusive (no mkdir -p)
        ├── .pipeline_tmp_owner                  # Mode 0600, atomic ownership marker
        └── reassembly/                          # Mode 0700
            ├── concatenated_R1.fastq.gz.part    # Staging write
            ├── concatenated_R2.fastq.gz.part
            ├── concatenated_R1.fastq.gz         # Final validated concatenated reads
            └── concatenated_R2.fastq.gz
```

Workspaces are created with `umask 077` and `mkdir --mode=700`. If `<RUN_ID>` already exists, the pipeline aborts immediately rather than reusing or deleting existing contents, preventing race conditions between concurrent executions.

### Ownership Marker

Upon workspace initialization, an ownership marker (`.pipeline_tmp_owner`) is written atomically (via temporary file and `mv`) with mode `0600`:

```ini
schema_version=1
run_id=20261102T153000Z_12345_6789a
pipeline_pid=12345
created_at=2026-11-02T15:30:00Z
tmpdir_source=cli
tmpdir_base=/scratch/pipeline_tmp
run_tmpdir=/scratch/pipeline_tmp/metagenomics_pipeline/20261102T153000Z_12345_6789a
```

### Temporary Read Concatenation Strategy

MetaWRAP's `reassemble_bins` module requires all paired-end reads across all samples in the cohort to be provided as a single paired FASTQ set. The pipeline produces these files inside `REASSEMBLY_TMPDIR` using:
- **Strategy**: `gzip-member-concat` (Version `1`)
- **Validation**: `gzip-integrity-and-source-manifest` (Version `1`)

#### Multi-Member Gzip Concatenation Compatibility
The `gzip-member-concat` strategy leverages standard RFC 1952 multi-member gzip streams: concatenating individual `.fastq.gz` files via binary concatenation (`cat file1.fastq.gz file2.fastq.gz > concatenated.fastq.gz`). This stream is fully compliant with POSIX `gzip -dc`, `zcat`, `pigz`, Bowtie2, MEGAHIT, and the MetaWRAP runtime without requiring costly uncompress-recompress cycles.

#### Validation Scope
Concatenation is performed strictly in active manifest order to `.part` files, verified to be non-empty, validated with `gzip -t`, and renamed into place.
> [!NOTE]
> `gzip -t` confirms archive and decompression stream integrity. It does **not** validate mate-pairing synchronization between R1 and R2 reads. Pairing guarantees derive from the active sample manifest, which enforces identical sample order and paired integrity during pre-flight.

### Storage Space and Inode Verification

Prior to concatenation, the pipeline calculates total input size:
$$\text{concatenated\_input\_bytes} = \text{read1\_total\_size\_bytes} + \text{read2\_total\_size\_bytes}$$
$$\text{required\_concat\_bytes} = \lceil \text{concatenated\_input\_bytes} \times 1.10 \rceil$$

Using `df -Pk` and `df -Pi`, the pipeline verifies:
1. **Available Bytes**: At least $1.10\times$ the combined input read size is free on the target filesystem.
2. **Available Inodes**: At least 16 free inodes (`REQUIRED_FREE_INODES=16`) are available in the workspace.

> [!WARNING]
> This estimation covers pipeline-managed temporary read concatenation. It does **not** estimate or guarantee the internal scratch memory or disk requirements of MetaWRAP's internal SPAdes assembler.

### Apptainer and Mamba Runtime Integration

1. **Apptainer Container Bind**:
   - The per-run workspace `PIPELINE_RUN_TMPDIR` is bound into the container in read-write mode (`-B "${PIPELINE_RUN_TMPDIR}:${PIPELINE_RUN_TMPDIR}:rw"`).
   - The mount is strictly isolated to reassembly; the base directory is never bound pipeline-wide.
   - Container visibility and write permissions are tested before running MetaWRAP.
2. **Environment Variable Propagation**:
   - `TMPDIR`, `TMP`, and `TEMP` are passed exclusively to the reassembly process (`--env` and `APPTAINERENV_` for Apptainer; inline command prefix for Mamba).
   - `TMPDIR` is **never globally exported** across the pipeline host environment.

### Safe Ownership-Verified Cleanup

Removal of temporary workspaces is governed by `validate_run_tmpdir_for_cleanup`:
- Path is non-empty, absolute, and matches `${PIPELINE_TMPDIR_BASE}/metagenomics_pipeline/${RUN_ID}`;
- Not `/`, `$HOME`, `OUTDIR`, base, or parent directory;
- Canonical path matches resolved path;
- Neither workspace, reassembly subdirectory, nor concatenated files are symlinks;
- `.pipeline_tmp_owner` exists, is regular, and matches `RUN_ID`, `pipeline_pid`, `tmpdir_base`, and `run_tmpdir`.

#### Cleanup Contexts

| Trigger | Condition | Pipeline Action | Workspace Outcome |
|---|---|---|---|
| **Success** | `KEEP_TEMPORARY=false` | Controlled cleanup after checkpoint verified | Removed (`removed`) |
| **Success** | `KEEP_TEMPORARY=true` | Diagnostic preservation requested | Preserved (`preserved-by-flag`) |
| **Failure** | Child processes running | Safety guard: children may still access files | Preserved (`preserved-active-children`) |
| **Failure** | No child processes, valid owner | Controlled error cleanup | Removed (`removed`) |
| **Signal (INT/TERM/HUP)** | Any | Child processes cannot be guaranteed halted | Preserved (`preserved-after-signal`) |

Any diagnostic workspace preserved upon failure, signal, or via `--keep-temporary` can be safely inspected and manually deleted when no longer needed:
```bash
rm -rf /path/to/metagenomics_pipeline/<RUN_ID>
```

### Checkpoints and Stage Reusability

- The `07_reassembly.config` checkpoint records strategy name, strategy version, validation name, validation version, active manifest fingerprint, active sample IDs, source pair count, and total input read bytes.
- Temporary directory paths, file mtimes, and `--keep-temporary` flags are **not** recorded in checkpoints. Changing `--tmpdir` does not invalidate an existing valid reassembly checkpoint.
- If reassembly action is `reuse`, the stage completes immediately from existing products in `OUTDIR/07_reassembly/` without creating temporary directories or concatenating reads.
- Preserved temporary files are **never reused** as caches by subsequent runs. Every fresh run (`run`, `force`, `restart`) generates its own isolated workspace and re-concatenates active reads.

---

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
| `make test-error-record` | Run Python unit and integration test suite for fatal error capture |
| `make test` | Run all test suites (including error capture) |
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

MEGAHIT can be disabled with `--no-megahit`; provide a valid assembly with
`--assembly-in /path/to/assembly.fa` only if an active downstream stage needs
one. A taxonomy-only run does not require an assembly.

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
