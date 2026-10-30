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
| `--vamb-binsplit-separator CHAR` | Advanced Vamb binsplitting separator; disabled unless explicitly supplied |
| `--refinement-mode MODE` | Consensus refinement policy: `auto`, `always`, or `never` |
| `--refinement` / `--no-refinement` | Aliases for `--refinement-mode always` / `never` |
| `--bins-source SOURCE` | `auto`, `vamb`, `metabinner`, `comebin`, or `refined`; default: `auto` |
| `--bam-dir DIR` | External `<sample>.sorted.bam` files and `.bai`/`.csi` indexes |
| `--bam-manifest FILE` | External BAM manifest with `sample_id`, `bam`, and `index` columns |
| `--jgi-depth-file FILE` | External JGI `jgi_summarize_bam_contig_depths` table |
| `--resume` | Reuse valid checkpoints (default) |
| `--force-stage STAGE` | Force one stage; repeatable, rebuilds its effective descendants |
| `--restart-from STAGE` | Destructively clears one stage and its effective descendants |
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
| `--refinement` | `--no-refinement` | `auto` (consensus refinement is resolved from the selected binners and bins source) |
| `--checkm2` | `--no-checkm2` | Enabled |
| `--mag-filtering` | `--no-mag-filtering` | Enabled |
| `--coverm` | `--no-coverm` | Enabled |
| `--summary` | `--no-summary` | Enabled |
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
External bin directories are not supported as `--bins-source`; choose a
selected binner or `refined`. With binning disabled, dependent MAG stages must
also be disabled explicitly and `--bins-source` must remain `auto`.

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
│   └── final_summary.tsv
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
| summary | `12_summary/` |

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
