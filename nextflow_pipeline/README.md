# Fiber-seq Nextflow Pipeline

A Nextflow pipeline for processing PacBio HiFi fiber-seq data through alignment, nucleosome calling, QC, pileup generation, bigWig conversion, heatmap visualization, and optional advanced analysis (k-mer phasing, FIRE peak calling).

## Requirements

- **Nextflow** >= 23.04.0
- **Singularity** or **Docker** (for containerized execution)
- **Reference genome** FASTA file
- **Chromosome sizes** file (for bigWig generation)

## Quick Start

```bash
nextflow run ./main.nf \
    -profile singularity \
    --input samples.tsv \
    --fasta /path/to/reference.fa \
    --chrom_sizes /path/to/chrom.sizes \
    --outdir results
```

## YAML Configuration

For complex analyses, use a YAML params file instead of command-line arguments:

```bash
# 1. Copy the template
cp params_template.yaml my_analysis.yaml

# 2. Edit parameters
vim my_analysis.yaml

# 3. Run with params file
nextflow run main.nf -profile singularity,minerva -params-file my_analysis.yaml

# 4. Override specific params if needed
nextflow run main.nf -profile singularity,minerva -params-file my_analysis.yaml --outdir different_dir
```

The `params_template.yaml` file contains all available parameters with documentation. Key sections:

- **Input/Output**: Required paths for samplesheet, output directory
- **Reference Genome**: FASTA, optional pre-built index, chromosome sizes
- **Step 1-6**: Part 1 processing options (alignment, nucleosomes, QC, pileup, bigWig, heatmaps)
- **Step 7-9**: Part 2 advanced analysis (chromosome filtering, k-mer phasing, FIRE)
- **Resource Limits**: Max CPUs, memory, and time

## Input Samplesheet

The pipeline requires a TSV file with 4 columns (no header):

```
sample_name    raw_bam_path    tissue    sequencer
```

Example (`samples.tsv`):
```
sample1    /path/to/sample1.hifi.bam    brain    Revio
sample2    /path/to/sample2.hifi.bam    liver    Sequel2
```

## Pipeline Steps

### Part 1: Core Processing

| Step | Process | Description |
|------|---------|-------------|
| 1 | PBMM2_INDEX | Create reference index (runs once) |
| 1 | PBMM2_ALIGN | Align reads with pbmm2 |
| 1b | SEQUENCING_QC | PacBio sequencing QC (optional) |
| 2 | FIBERTOOLS_ADD_NUCLEOSOMES | Call nucleosomes/MSPs with fibertools |
| 2b | SAMTOOLS_INDEX | Index ft.bam files |
| 3 | FIBERSEQ_QC | Fiber-seq QC metrics (optional) |
| 4 | FIBERTOOLS_PILEUP | Create pileup files |
| 5a | PILEUP_TO_BEDGRAPH | Convert pileups to bedGraph |
| 5b | BEDGRAPH_TO_BIGWIG | Convert bedGraph to bigWig |
| 6a | COMPUTE_MATRIX | deepTools matrix computation |
| 6b | PLOT_HEATMAP | Generate heatmaps |

### Part 2: Advanced Analysis (Optional)

| Step | Process | Description |
|------|---------|-------------|
| 7 | SAMTOOLS_FILTER_CHROMS | Filter BAM to specific chromosomes |
| 8 | KMER_PHASING | k-mer variant phasing (Snakemake wrapper) |
| 9 | FIRE | FIRE peak calling (Snakemake wrapper) |

## Execution DAG

Nextflow executes processes based on data dependencies, not source code order. Once Step 2b produces the indexed BAM, Part 1 and Part 2 branches run **in parallel**:

```
                          ch_samples (from samplesheet)
                                   │
                    ┌──────────────┴──────────────┐
                    ↓                             ↓
              PBMM2_ALIGN                   SEQUENCING_QC
                    │                        (optional)
                    ↓
      FIBERTOOLS_ADD_NUCLEOSOMES (Step 2)
                    │
                    ↓
             SAMTOOLS_INDEX (Step 2b)
                    │
                    ↓ indexed_bam
    ┌───────────────┼───────────────┐
    ↓               ↓               ↓
FIBERSEQ_QC   FIBERTOOLS_PILEUP   SAMTOOLS_FILTER_CHROMS
(optional)          │                     │
                    ↓                     ↓
           PILEUP_TO_BEDGRAPH      KMER_PHASING (optional)
                    │                     │
                    ↓                     ↓
           BEDGRAPH_TO_BIGWIG        FIRE (optional)
                    │
                    ↓
           COMPUTE_MATRIX (optional)
                    │
                    ↓
            PLOT_HEATMAP (optional)
```

**Key point:** Steps 3-6 (pileups, bigWigs, heatmaps) and Steps 7-9 (filtering, phasing, FIRE) execute concurrently once `indexed_bam` is available.

To visualize your specific run's DAG:
```bash
nextflow run main.nf [options] -preview
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet TSV |
| `--fasta` | Path to reference genome FASTA |
| `--outdir` | Output directory |

### Reference Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fasta_index` | null | Pre-built MMI index (skips indexing) |
| `--chrom_sizes` | null | Chromosome sizes file (required for bigWig) |

### Alignment Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_alignment` | false | Skip alignment step (BAMs in samplesheet are pre-aligned) |
| `--pbmm2_preset` | HIFI | Preset: SUBREAD, CCS, HIFI, ISOSEQ |

### Pileup Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pileup_region` | null | Optional region (chr:start-end) |
| `--ml_threshold` | null | ML threshold for ft pileup |
| `--marks` | m6a,nuc,cpg | Comma-separated marks for bigWig |

### Heatmap Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--regions_bed` | null | BED file for heatmaps |
| `--matrix_mode` | reference-point | Mode: reference-point or scale-regions |
| `--reference_point` | center | Reference point: center, TSS, TES |
| `--upstream` | 2000 | Upstream distance (bp) |
| `--downstream` | 2000 | Downstream distance (bp) |
| `--region_body_length` | 5000 | Body length for scale-regions mode |
| `--scale_upstream` | 3000 | Upstream for scale-regions mode |
| `--scale_downstream` | 3000 | Downstream for scale-regions mode |
| `--heatmap_height` | 15 | Heatmap height |
| `--heatmap_width` | 4 | Heatmap width |
| `--heatmap_colormap_m6a` | Greens | Colormap for m6A heatmaps |
| `--heatmap_colormap_nuc` | Blues | Colormap for nucleosome heatmaps |
| `--heatmap_colormap_cpg` | Reds | Colormap for CpG heatmaps |

### Part 2 Parameters

#### Chromosome Filtering (Step 7)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--filter_chromosomes` | null | Space-separated chromosome list |

#### K-mer Phasing (Step 8)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_kmer_phasing` | false | Enable k-mer phasing |
| `--kmer_use_filtered` | false | Use filtered BAM from step 7 |
| `--kmer_align` | false | Perform alignment in k-mer phasing |
| `--resume_kmer_phasing` | false | Resume incomplete Snakemake run |
| `--kmer_parental_mode` | false | Enable parental data for phasing |
| `--kmer_maternal_bam` | null | Path to maternal sequencing data |
| `--kmer_paternal_bam` | null | Path to paternal sequencing data |

#### FIRE Peak Calling (Step 9)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_fire` | false | Enable FIRE peak calling |
| `--fire_input_source` | nucleosomes | Input: "nucleosomes" or "kmer_phasing" |
| `--fire_ref_name` | hg38 | UCSC genome name for track hub |
| `--fire_ont` | false | ONT Fiber-seq mode (vs PacBio) |
| `--fire_max_threads` | null | Max threads for distributed steps |
| `--fire_keep_chromosomes` | null | Regex filter (e.g., "chr[0-9XY]+$") |
| `--fire_excludes` | null | BED files to exclude from null distribution |
| `--fire_min_contig_length` | 0 | Skip contigs smaller than this (bp) |

FIRE peak calling threshold parameters (null = use FIRE defaults):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fire_max_peak_fdr` | null | FDR for FIRE peaks (default: 0.05) |
| `--fire_min_fire_fdr` | null | FDR for individual elements (default: 0.10) |
| `--fire_min_coverage` | null | Min coverage for peaks (default: 4) |
| `--fire_coverage_within_n_sd` | null | Filter peaks beyond N SDs (default: 5) |
| `--fire_min_msp` | null | Min MSPs per read (default: 10) |
| `--fire_min_ave_msp_size` | null | Min average MSP size (default: 10) |

#### Snakemake Workflows

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kmer_pixi_manifest` | null | pixi.toml for k-mer phasing (k-mer-variant-phasing repo) |
| `--kmer_snakemake_profile` | null | Snakemake profile for k-mer phasing (e.g., profiles/lsf-executor) |
| `--fire_pixi_manifest` | null | pixi.toml for FIRE (FIRE repo) |
| `--fire_snakemake_profile` | null | Snakemake profile for FIRE (e.g., profiles/lsf-executor) |

### QC Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_sequencing_qc` | false | Run PacBio sequencing QC |
| `--run_fiberseq_qc` | false | Run fiberseq-qc suite |
| `--smrt_root` | null | SMRT Tools installation path |
| `--fiberseq_qc_root` | null | fiberseq-qc installation path |

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution |
| `singularity` | Enable Singularity containers |
| `docker` | Enable Docker containers |
| `lsf` | Generic LSF scheduler |
| `slurm` | Generic SLURM scheduler |
| `test` | Minimal resources for testing |

Combine profiles with comma separation:
```bash
-profile singularity,minerva
```

## Usage Examples

### Basic Run (Local)

```bash
nextflow run ./main.nf \
    -profile singularity \
    --input samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --outdir results
```

### Run on HPC with LSF

For HPC environments that require proxy settings for container downloads:

```bash
# 1. Edit run_fiberseq_pipeline.lsf to set your params file and project
# 2. Submit the job
cd nextflow_pipeline
bsub < run_fiberseq_pipeline.lsf

# 3. Monitor progress
bjobs
tail -f logs/fiberseq_pipeline_*.stdout
```

**Alternative: Interactive session**
```bash
# If your HPC requires proxy, set environment variables:
# export http_proxy=http://your-proxy.institution.edu:3128
# export https_proxy=http://your-proxy.institution.edu:3128

# Load modules and activate conda
module load singularity
source ~/miniforge3/etc/profile.d/conda.sh
conda activate nextflow

# Run pipeline
nextflow run ./main.nf \
    -profile singularity,lsf \
    -params-file my_analysis.yaml \
    -resume
```

### With Heatmaps

```bash
nextflow run ./main.nf \
    -profile singularity,minerva \
    --input samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --regions_bed promoters.bed \
    --matrix_mode reference-point \
    --upstream 3000 \
    --downstream 3000 \
    --outdir results
```

### With FIRE Peak Calling

```bash
nextflow run ./main.nf \
    -profile singularity,minerva \
    --input samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --filter_chromosomes "chr1 chr2 chr3 chr4 chr5" \
    --run_fire \
    --fire_ref_name hg38 \
    --fire_keep_chromosomes 'chr[0-9XY]+$' \
    --fire_pixi_manifest /path/to/FIRE/pixi.toml \
    --fire_snakemake_profile /path/to/FIRE/profiles/lsf-executor \
    --outdir results
```

### With K-mer Phasing and FIRE

```bash
# Run k-mer phasing, then use phased BAM for FIRE
nextflow run ./main.nf \
    -profile singularity,minerva \
    --input samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --run_kmer_phasing \
    --kmer_pixi_manifest /path/to/k-mer-variant-phasing/pixi.toml \
    --kmer_snakemake_profile /path/to/k-mer-variant-phasing/profiles/lsf-executor \
    --run_fire \
    --fire_input_source kmer_phasing \
    --fire_pixi_manifest /path/to/FIRE/pixi.toml \
    --fire_snakemake_profile /path/to/FIRE/profiles/lsf-executor \
    --outdir results
```

### Using YAML Params File

```bash
# Recommended for complex analyses
nextflow run ./main.nf \
    -profile singularity,minerva \
    -params-file my_analysis.yaml
```

### Starting from Pre-aligned BAMs

```bash
nextflow run ./main.nf \
    -profile singularity,minerva \
    --input aligned_samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --skip_alignment \
    --outdir results
```

**Note:** When using `--skip_alignment`, the `bam_path` column in your samplesheet should contain paths to pre-aligned BAM files (not raw HiFi BAMs). Sequencing QC is also skipped since it is designed for raw BAMs.

### Dry Run (Preview DAG)

```bash
nextflow run ./main.nf \
    -profile singularity \
    --input samples.tsv \
    --fasta reference.fa \
    --outdir test \
    -preview
```

### Resume Failed Run

```bash
nextflow run ./main.nf \
    -profile singularity,minerva \
    --input samples.tsv \
    --fasta reference.fa \
    --chrom_sizes chrom.sizes \
    --outdir results \
    -resume
```

## Output Directory Structure

```
results/
├── pipeline_info/
│   ├── execution_timeline_*.html
│   ├── execution_report_*.html
│   ├── execution_trace_*.txt
│   └── pipeline_dag_*.html
├── reference/
│   └── *.mmi                       # Reference index
└── {sample_name}/
    ├── 01_aligned/
    │   ├── {sample}.aligned.bam
    │   └── {sample}.aligned.bam.bai
    ├── 01_sequencing_qc/           # If --run_sequencing_qc
    │   └── {sample}.sequencing_qc_report.pdf
    ├── 02_nucleosomes/
    │   ├── {sample}.ft.bam
    │   └── {sample}.ft.bam.bai
    ├── 03_fiberseq_qc/             # If --run_fiberseq_qc
    │   └── qc_output/*
    ├── 04_pileups/
    │   ├── {sample}.genome_wide.pileup.bed
    │   ├── {sample}.*.m6a.bedgraph
    │   ├── {sample}.*.nuc.bedgraph
    │   └── {sample}.*.cpg.bedgraph
    ├── 05_bigwigs/
    │   ├── {sample}.*.m6a.bw
    │   ├── {sample}.*.nuc.bw
    │   └── {sample}.*.cpg.bw
    ├── 06_heatmaps/                # If --regions_bed provided
    │   ├── {sample}.*.matrix.gz
    │   ├── {sample}.*.heatmap.png
    │   └── {sample}.*.heatmap.pdf
    ├── 07_filtered/                # If --filter_chromosomes
    │   ├── {sample}.filtered.bam
    │   └── {sample}.filtered.bam.bai
    ├── 08_kmer_phasing/            # If --run_kmer_phasing
    │   └── output/**
    └── 09_fire/                    # If --run_fire
        └── output/**
```

## Resource Requirements

Default resource allocations (can be scaled on retry):

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| PBMM2_ALIGN | 12 | 64 GB | 8h |
| FIBERTOOLS_ADD_NUCLEOSOMES | 8 | 32 GB | 4h |
| FIBERTOOLS_PILEUP | 8 | 16 GB | 4h |
| COMPUTE_MATRIX | 8 | 16 GB | 2h |
| KMER_PHASING | 16 | 64 GB | 24h |
| FIRE | 16 | 64 GB | 24h |

Override max resources:
```bash
--max_cpus 24 --max_memory 128.GB --max_time 72.h
```

## Containers

The pipeline uses Biocontainers images:

| Process | Container |
|---------|-----------|
| pbmm2 | `biocontainers/pbmm2:1.14.99--h9ee0642_0` |
| fibertools | `biocontainers/fibertools-rs:0.7.1--h3b373d1_0` |
| samtools | `biocontainers/samtools:1.19.2--h50ea8bc_1` |
| bedGraphToBigWig | `biocontainers/ucsc-bedgraphtobigwig:447--h954228d_0` |
| deepTools | `biocontainers/deeptools:3.5.4--pyhdfd78af_1` |

## Troubleshooting

### Proxy connection refused on HPC

If you see errors like:
```
curl: (7) Failed to connect to proxy port: Connection refused
ERROR: Cannot download nextflow required file
```

You need to set proxy environment variables before running Nextflow. Either:
1. Use the LSF submission script (`run_fiberseq_pipeline.lsf`) which handles this automatically
2. Export proxy variables manually for your institution's proxy

### Pipeline fails to find container

Ensure Singularity is available and caching is configured:
```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/cache
```

### Out of memory errors

Increase max memory or let the pipeline retry with more resources:
```bash
--max_memory 256.GB
```

### Resume not working

Ensure you're running from the same directory with the same parameters. The work directory must be intact.

### LSF job submission issues

Check your project allocation is set correctly in your profile config or nextflow.config:
```groovy
clusterOptions = { "-P acc_YourProject" }
```

## Citation

If you use this pipeline, please cite:

- Nextflow: https://doi.org/10.1038/nbt.3820
- fibertools-rs: https://github.com/fiberseq/fibertools-rs
- FIRE: https://github.com/fiberseq/FIRE
