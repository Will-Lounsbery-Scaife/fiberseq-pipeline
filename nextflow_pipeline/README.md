# Fiber-seq Nextflow Pipeline

## Pipeline Steps

### Part 1: Core Processing

| Step | Process | Description |
|------|---------|-------------|
| 1 | PBMM2_INDEX | Create reference index (runs once) |
| 1 | PBMM2_ALIGN | Align reads with pbmm2 |
| 1b | SEQUENCING_QC | PacBio sequencing QC |
| 2 | FIBERTOOLS_ADD_NUCLEOSOMES | Call nucleosomes/MSPs with fibertools |
| 2b | SAMTOOLS_INDEX | Index ft.bam files |
| 3 | FIBERSEQ_QC | Fiber-seq QC metrics |
| 4 | FIBERTOOLS_PILEUP | Create pileup files |
| 5a | PILEUP_TO_BEDGRAPH | Convert pileups to bedGraph |
| 5b | BEDGRAPH_TO_BIGWIG | Convert bedGraph to bigWig |
| 6a | COMPUTE_MATRIX | deepTools matrix computation |
| 6b | PLOT_HEATMAP | Generate heatmaps |

### Part 2: Advanced Analysis

| Step | Process | Description |
|------|---------|-------------|
| 7 | SAMTOOLS_FILTER_CHROMS | Filter BAM to specific chromosomes |
| 8 | KMER_PHASING | k-mer variant phasing (Snakemake wrapper) |
| 9 | FIRE | FIRE peak calling (Snakemake wrapper) |

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
