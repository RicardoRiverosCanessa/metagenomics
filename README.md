Shotgun Metagenomics Pipeline

A robust, modular, production-ready Bash pipeline for high-throughput shotgun metagenomic assembly, binning, refinement, taxonomic classification, and abundance profiling.The pipeline implements a hybrid execution architecture: compiled binaries (C/C++/Rust) are called directly via a custom $PATH, while complex Python-based environments and deep learning wrappers (Vamb, MetaBinner, COMEBin, metaWRAP, CheckM2, GTDB-Tk) are safely sandboxed and invoked via mamba run (completely avoiding dangerous conda activate calls inside subshells). It also features built-in checkpointing to allow seamless resumption after failures, a dry-run mode, and automated pre-flight validations.   

Workflow OverviewQC (fastp): 
1) Quality filtering and adapter trimming (optional).   
2) Assembly (MEGAHIT): De novo metagenomic assembly with contig header cleaning and minimum length filtering.
3) Taxonomy (Kraken2 + Bracken): Read-level taxonomic profiling and abundance re-estimation.
4) Mapping (Bowtie2 + SAMtools + jgi_summarize_bam_contig_depths): Read mapping, sorting, indexing, and generation of the shared depth matrix.
5) Binning (Multi-binner strategy):
  5a) Vamb (using the shared JGI depth matrix converted via awk).
  5b) MetaBinner (reusing Step 04 BAMs/depths with synchronized contig IDs).
  5c) COMEBin (reusing the Step 04 mapping directory).
6) Refinement (metaWRAP bin_refinement): Consensus binning combining Vamb (-A), COMEBin (-B), and MetaBinner (-C).
7) Reassembly (metaWRAP reassemble_bins): Optional bin reassembly using pooled clean reads.   
8) Quality Assessment (CheckM2): Completeness and contamination estimation for all bins.
9) Filtering: Selection of high-quality/medium-quality Metagenome-Assembled Genomes (MAGs) based on user-defined thresholds.
10) Taxonomic Classification (GTDB-Tk): Automated genome taxonomy assignment (optional).
11) Abundance Profiling (CoverM): Per-sample relative abundance estimation for filtered MAGs.
12) Summary (R): Base-R left-join merging CheckM2, GTDB-Tk, and CoverM outputs into a unified final .tsv matrix.   


#~~~~~~~~~~~~~~~~~~~~~~~~~

#Requirements

#~~~~~~~~~~~~~~~~~~~~~~~~~

System RequirementsOS: Linux (tested on Fedora 44).

Package Manager: Mamba / Conda installed.   

Core Tools / Binaries: fastp, megahit, kraken2, bracken, bowtie2, samtools, jgi_summarize_bam_contig_depths, coverm, Rscript.   

Mamba Environments: Separate isolated environments for vamb, metabinner_env, comebin_env, metawrap_env, checkm2_env, and gtdbtk_env.   


#~~~~~~~~~~~~~~~~~~~~~~~~~

#Installation

#~~~~~~~~~~~~~~~~~~~~~~~~~

Installation & ConfigurationClone or download the script to your working environment.   

Make the script executable:   Bashchmod +x run_metagenomics_pipeline.sh

Open run_metagenomics_pipeline.sh and edit the CONFIGURATION block at the top to configure your paths, threads, memory limits, and database locations. (Alternatively, database paths can be exported as environmental variables prior to execution, e.g., export KRAKEN2_DB="/path/to/db").   


#~~~~~~~~~~~~~~~~~~~~~~~~~

Usage

#~~~~~~~~~~~~~~~~~~~~~~~~~

Run the pipeline from your terminal:   Bash./run_metagenomics_pipeline.sh

Command-Line Flags--dry-run: Inspect all generated commands, sub-scripts, and logic blocks without executing them.   

--force: Ignore checkpoints and rerun all steps from scratch.   

-h, --help: Display inline documentation and usage guidelines.   
