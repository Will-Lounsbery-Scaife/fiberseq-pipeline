# Fiber-seq Nextflow Pipeline

A Nextflow pipeline for processing PacBio HiFi fiber-seq data through alignment, nucleosome calling, pileup generation, bigWig conversion, heatmap visualization, and advanced analysis (k-mer phasing, FIRE peak calling).

## Requirements

- **Nextflow** >= 23.04.0
- **Singularity** or **Docker**
- **Reference genome** FASTA file
- **Chromosome sizes** file

## Input Samplesheet

TSV file with 4 columns (no header):

```
sample_name    raw_bam_path    tissue    sequencer
sample1        /path/to/sample1.hifi.bam    brain    Revio
```

## Execution DAG

Once Step 2b produces the indexed BAM, Part 1 and Part 2 branches run **in parallel**:

```
                          ch_samples (from samplesheet)
                                   │
                    ┌──────────────┴──────────────┐
                    ↓                             ↓
              PBMM2_ALIGN                   SEQUENCING_QC
                    │
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
                    │                     │
                    ↓                     ↓
           PILEUP_TO_BEDGRAPH       KMER_PHASING
                    │                     │
                    ↓                     ↓
           BEDGRAPH_TO_BIGWIG          FIRE
                    │
                    ↓
              COMPUTE_MATRIX
                    │
                    ↓
              PLOT_HEATMAP
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet TSV |
| `--fasta` | Path to reference genome FASTA |
| `--outdir` | Output directory |
| `--chrom_sizes` | Chromosome sizes file (required for bigWig) |

### Reference & Alignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fasta_index` | null | Pre-built MMI index (skips indexing) |
| `--skip_alignment` | false | Skip alignment (use pre-aligned BAMs) |
| `--pbmm2_preset` | HIFI | Preset: SUBREAD, CCS, HIFI, ISOSEQ |

### Pileup & BigWig

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pileup_region` | null | Region to pileup (chr:start-end) |
| `--marks` | m6a,nuc,cpg | Marks to extract for bigWig |

### Heatmaps

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--regions_bed` | null | BED file for heatmaps (enables heatmap generation) |
| `--matrix_mode` | reference-point | Mode: reference-point or scale-regions |
| `--reference_point` | center | Reference point: center, TSS, TES |
| `--upstream` | 2000 | Upstream distance (bp) |
| `--downstream` | 2000 | Downstream distance (bp) |

### Part 2: Chromosome Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--filter_chromosomes` | null | Space-separated chromosome list |

### Part 2: K-mer Phasing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_kmer_phasing` | false | Enable k-mer phasing |
| `--kmer_pixi_manifest` | null | pixi.toml from k-mer-variant-phasing repo |
| `--kmer_snakemake_profile` | null | Snakemake profile path |
| `--kmer_use_filtered` | false | Use filtered BAM from step 7 |

### Part 2: FIRE Peak Calling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_fire` | false | Enable FIRE peak calling |
| `--fire_pixi_manifest` | null | pixi.toml from FIRE repo |
| `--fire_snakemake_profile` | null | Snakemake profile path |
| `--fire_input_source` | nucleosomes | Input: "nucleosomes" or "kmer_phasing" |
| `--fire_ref_name` | hg38 | UCSC genome name for track hub |
| `--fire_keep_chromosomes` | null | Regex filter (e.g., "chr[0-9XY]+$") |

### QC

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_sequencing_qc` | false | Run PacBio sequencing QC |
| `--run_fiberseq_qc` | false | Run fiberseq-qc suite |
| `--smrt_root` | null | SMRT Tools installation path |
| `--fiberseq_qc_root` | null | fiberseq-qc installation path |

## Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Enable Singularity containers |
| `docker` | Enable Docker containers |
| `lsf` | Generic LSF scheduler |
| `slurm` | Generic SLURM scheduler |
| `test` | Minimal resources for testing |

Combine profiles: `-profile singularity,lsf`
