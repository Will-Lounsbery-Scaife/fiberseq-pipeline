# Fiberseq Pipeline

A complete shell-based pipeline for processing and analyzing PacBio HiFi fiber-seq data, from raw BAM files to publication-quality heatmaps.

---

## Biological Background

### What is Fiber-seq?

Fiber-seq is a single-molecule long-read sequencing method that enables mapping chromatin structure/accessibility at the resolution of individual DNA fibers. It works by treating nuclei with a N6-adenine methyltransferase (such as Hia5). The MTase acts like a stencil, methylating adenine bases in accessible regions while leaving nucleosome-bound or TF-bound DNA unmodified. The pattern of methylated bases along individual DNA strands shows where nucleosomes were positioned and where the DNA was open and active.

Original paper: https://www.science.org/doi/abs/10.1126/science.aaz1646

Fibertools: https://fiberseq.github.io/index.html

### What is FIRE?

FIRE ("Fiber-seq Inferred Regulatory Elements") is a peak-calling and annotation method that takes Fiber-seq single-molecule data and classifies which methyltransferase-sensitive patches (MSPs; stretches of m6A-labeled, exposed DNA on a single read) correspond to regulatory element "actuation" versus ordinary internucleosomal linker-like accessibility. 

https://www.biorxiv.org/content/10.1101/2024.06.14.599122v2
https://fiberseq.github.io/fire/fire.html


---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Pipeline Structure](#pipeline-structure)
4. [Configuration](#configuration)
5. [Running the Pipeline](#running-the-pipeline)
   - [Part 1: Core Processing](#part-1-core-processing)
     - [Step 1: Alignment and Sequencing QC](#step-1-alignment-and-sequencing-qc)
     - [Step 2: Add Nucleosomes](#step-2-add-nucleosomes)
     - [Step 3: Fiber-seq QC](#step-3-fiber-seq-qc)
     - [Step 4: Create Pileups](#step-4-create-pileups)
     - [Step 5: Pileup to bigWig](#step-5-pileup-to-bigwig)
     - [Step 6: Create Heatmaps](#step-6-create-heatmaps)
   - [Part 2: Phasing and peak calling](#part-2-phasing-and-peak-calling)
     - [Step 7: Filter Chromosomes](#step-7-filter-chromosomes)
     - [Step 8: K-mer Variant Phasing](#step-8-k-mer-variant-phasing)
     - [Step 9: FIRE Peak Calling](#step-9-fire-peak-calling)
6. [Output Structure](#output-structure)

---

## Overview

This pipeline processes fiber-seq data through the following steps:

**Part 1: Core Processing**
1. **Alignment and QC**: Align reads to reference genome and generate QC reports
2. **Add Nucleosomes**: Call methylation-sensitive patches (MSPs) and nucleosomes
3. **Fiber-seq QC**: Generate fiber-seq-specific quality metrics
4. **Create Pileups**: Generate coverage pileups for genomic regions
5. **Pileup to bigWig**: Convert pileups to bigWig format for visualization
6. **Create Heatmaps**: Generate heatmaps centered on genomic features

**Part 2: Haplotype Phasing and FIRE**
7. **Filter Chromosomes**: Filter BAMs to specific chromosomes (e.g., canonical human chromosomes)
8. **K-mer Variant Phasing**: Phase variants (Snakemake)
9. **FIRE Peak Calling**: Identify Fiber-seq Inferred Regulatory Elements using FDR-based peak calling (Snakemake)

---

## Prerequisites

### Required Files

- Raw unaligned fiber-seq BAMs files (PacBio HiFi reads)
- Reference genome FASTA
- Chromosome sizes
- BED file with regions of interest (for heatmaps). I recommend using promoters.

### Required Software

- **PacBio SMRT Tools** (pbmm2, pbindex, dataset, runqc-reports)
- **samtools**
- **tcsh** (for fiberseq-qc)
- **fiberseq-qc**
- **fibertools** (ft)
- **bedGraphToBigWig** (UCSC utils)
- **deepTools** (computeMatrix, plotHeatmap)

---

## Pipeline Structure

```
pipeline/
├── part1_scripts/               # Steps 1-6: Preprocessing, QC, initial heatmaps
│   ├── 01_align_and_qc_wrapper.sh
│   ├── 01_align_and_qc_main.sh
│   ├── 02_add_nucleosomes_wrapper.sh
│   ├── 02_add_nucleosomes_main.sh
│   ├── 03_fiberseq_qc_wrapper.sh
│   ├── 03_fiberseq_qc_main.sh
│   ├── 04_create_pileups_wrapper.sh
│   ├── 04_create_pileups_main.sh
│   ├── 05_pileup_to_bigwig_wrapper.sh
│   ├── 05_pileup_to_bigwig_main.sh
│   ├── 06_create_heatmaps_wrapper.sh
│   └── 06_create_heatmaps_main.sh
├── part2_scripts/               # Steps 7-9: phasing, peak calling
│   ├── 07_filter_chroms_wrapper.sh
│   ├── 07_filter_chroms_main.sh
│   ├── 08_kmer_phasing_wrapper.sh
│   ├── 08_kmer_phasing_main.sh
│   ├── 09_fire_wrapper.sh
│   ├── 09_fire_main.sh
│   └── config/                  # Generated YAML configs
└── README.md
```

---

## Configuration

### 1. Create Sample Sheet

Create a tab-delimited file at `config/samples.tsv`:

```tsv
sample_name	raw_bam_path	tissue	sequencer
sample1	/path/to/sample1.bam	brain	revio
sample2	/path/to/sample2.bam	brain	revio
sample3	/path/to/sample3.bam	brain	revio
```

**Columns:**
- `sample_name`: Unique identifier for each sample
- `raw_bam_path`: Full path to raw unaligned BAM file
- `tissue`: Tissue type (optional metadata)
- `sequencer`: Sequencer used (optional metadata)

### 2. Prepare Reference Files

**Reference genome:**
```bash
# You should already have your reference FASTA
REFERENCE_FASTA="/path/to/reference.fa"
```

**Chromosome sizes file:**
```bash
# Create from reference FASTA
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > reference.chrom.sizes
```

**Regions BED file (for heatmaps):**
```bash
# Example: Gene TSS regions
# Format: chr  start  end  name  score  strand
chr1    1000000    1000001    gene1    .    +
chr1    2000000    2000001    gene2    .    -
```

### 3. Make Scripts Executable

```bash
chmod +x part1_scripts/*.sh part2_scripts/*.sh
```

---

## Running the Pipeline

### Part 1: Core Processing

#### Step 1: Alignment and Sequencing QC

**Purpose:** Align raw BAM files to reference genome and generate sequencing QC reports.

**Configuration:**

Edit `part1_scripts/01_align_and_qc_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
REFERENCE_FASTA="/path/to/reference.fa"
WORKDIR="/path/to/output/directory"
SMRT_ROOT="/path/to/smrttools"
THREADS=48
PRESET="HIFI"
```

**Run:**

```bash
bash part1_scripts/01_align_and_qc_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
├── 01_aligned/
│   ├── {sample_name}.aligned.bam
│   └── {sample_name}.aligned.bam.pbi
└── 02_sequencing_qc/
    ├── {sample_name}.alignmentset.xml
    └── {sample_name}.pbqc.report.pdf
```

---

#### Step 2: Add Nucleosomes

**Purpose:** Call methylation-sensitive patches (MSPs) and add nucleosome annotations using fibertools.

**Configuration:**

Edit `part1_scripts/02_add_nucleosomes_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"  # Same as Step 1
THREADS=48
FT_OPTIONS=""  # Optional fibertools parameters
```

**Run:**

```bash
# Activate fibertools environment
conda activate fibertools_env

# Run the script
bash part1_scripts/02_add_nucleosomes_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 03_nucleosomes/
    ├── {sample_name}.ft.bam
    ├── {sample_name}.ft.bam.bai
    ├── {sample_name}.ft.bam.pbi
    └── {sample_name}.ft.stats.txt
```

---

#### Step 3: Fiber-seq QC

**Purpose:** Generate fiber-seq-specific quality control metrics and visualizations.

**Configuration:**

Edit `part1_scripts/03_fiberseq_qc_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"  # Same as previous steps
FIBERSEQ_QC_ROOT="/path/to/fiberseq-qc"
THREADS=48
```

**Run:**

```bash
# Ensure tcsh is installed
which tcsh

# Run the script
bash part1_scripts/03_fiberseq_qc_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 04_fiberseq_qc/
    ├── qc_summary.txt
    ├── plots/
    ├── {sample_name}.qc_run_summary.txt
    └── [other QC outputs]
```

---

#### Step 4: Create Pileups

**Purpose:** Generate coverage pileups aggregating fiber-seq data across genomic regions.

**Configuration:**

Edit `part1_scripts/04_create_pileups_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"  # Same as previous steps
REGIONS_BED=""  # Optional: leave empty for genome-wide, or specify BED file
THREADS=48
FT_PILEUP_OPTIONS=""  # Optional fibertools pileup parameters
```

**Run:**

```bash
# Activate fibertools environment
conda activate fibertools_env

# Run the script
bash part1_scripts/04_create_pileups_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 05_pileups/
    ├── {sample_name}.genome_wide.pileup.bed
    └── {sample_name}.pileup_stats.txt
```

---

#### Step 5: Pileup to bigWig

**Purpose:** Split pileups by mark type (coverage, 5mC, nucleosomes) and convert to bigWig format.

**Configuration:**

Edit `part1_scripts/05_pileup_to_bigwig_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"  # Same as previous steps
GENOME_CHROMSIZES="/path/to/reference.chrom.sizes"
BEDGRAPHTOBIGWIG="/path/to/bedGraphToBigWig"

# Configure which columns to extract
declare -A PILEUP_COLUMNS=(
    ["coverage"]=4      # Column 4: coverage
    ["perc5mc"]=5       # Column 5: percent 5mC
    ["nucleosome"]=6    # Column 6: nucleosome occupancy
)
```

**Run:**

```bash
bash part1_scripts/05_pileup_to_bigwig_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 06_bigwigs/
    ├── {sample_name}.genome_wide.coverage.bw
    ├── {sample_name}.genome_wide.perc5mc.bw
    ├── {sample_name}.genome_wide.nucleosome.bw
    └── {sample_name}.bigwig_conversion_summary.txt
```

---

#### Step 6: Create Heatmaps

**Purpose:** Generate publication-quality heatmaps showing fiber-seq signals around genomic features.

**Configuration:**

Edit `part1_scripts/06_create_heatmaps_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"  # Same as previous steps
REGIONS_BED="/path/to/regions_of_interest.bed"
THREADS=48

# Matrix parameters
MATRIX_MODE="reference-point"  # or "scale-regions"
REFERENCE_POINT="center"       # "TSS", "TES", or "center"
UPSTREAM=3000
DOWNSTREAM=3000

# Heatmap parameters
COLORMAP="RdYlBu_r"
HEATMAP_HEIGHT=20
HEATMAP_WIDTH=10

# Which marks to plot
MARKS_TO_PLOT=("coverage" "perc5mc" "nucleosome")
```

**Run:**

```bash
# Activate deepTools environment
conda activate deeptools_env

# Run the script
bash part1_scripts/06_create_heatmaps_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 07_heatmaps/
    ├── {sample}.{regions}.coverage.matrix.gz
    ├── {sample}.{regions}.coverage.heatmap.png
    ├── {sample}.{regions}.coverage.heatmap.pdf
    ├── {sample}.{regions}.perc5mc.matrix.gz
    ├── {sample}.{regions}.perc5mc.heatmap.png
    ├── {sample}.{regions}.nucleosome.matrix.gz
    ├── {sample}.{regions}.nucleosome.heatmap.png
    ├── {sample}.{regions}.all_marks.matrix.gz
    ├── {sample}.{regions}.all_marks.heatmap.png
    ├── {sample}.{regions}.all_marks.heatmap.pdf
    └── {sample}.heatmap_summary.txt
```

---

### Part 2: Phasing and peak calling

Part 2 steps are optional and require additional setup. They use Snakemake workflows managed by pixi.

**Prerequisites for Part 2:**
- Singularity/Apptainer: `ml singularity/3.6.4`
- k-mer-variant-phasing repository (for Step 8)
- FIRE repository (for Step 9)

---

#### Step 7: Filter Chromosomes

**Purpose:** Filter nucleosome BAMs to specific chromosomes (e.g., remove contigs, keep only canonical chromosomes).

**Configuration:**

Edit `part2_scripts/07_filter_chroms_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
CHROMOSOMES="chr1 chr2 chr3 ... chrX chrY"  # Space-separated list
FILTER_THREADS=8
```

**Run:**

```bash
module load samtools/1.21
bash part2_scripts/07_filter_chroms_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 07_filter_chroms/
    ├── {sample_name}.filtered.bam
    └── {sample_name}.filtered.bam.bai
```

---

#### Step 8: K-mer Variant Phasing

**Purpose:** Phase variants using k-mer-based approach with the k-mer-variant-phasing Snakemake workflow.

**Configuration:**

Edit `part2_scripts/08_kmer_phasing_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
REFERENCE_FASTA="/path/to/reference.fa"
KMER_PHASING_ROOT="/path/to/k-mer-variant-phasing"

ALIGN="false"        # Set to "true" if alignment not done in Step 1
FILTERED_CHR="false" # Set to "true" to use filtered BAMs from Step 7
RESUME="false"       # Set to "true" to resume incomplete run

# Optional: Parental data for trio phasing
PARENTAL_MODE="false"
MATERNAL_BAM=""
PATERNAL_BAM=""
```

**Run:**

```bash
ml singularity/3.6.4
bash part2_scripts/08_kmer_phasing_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 08_kmer_phasing/
    └── results/{sample_name}/
        └── hiphase/
            ├── {sample_name}.bam
            └── {sample_name}.bam.bai
```

---

#### Step 9: FIRE Peak Calling

**Purpose:** Identify Fiber-seq Inferred Regulatory Elements (FIRE) using FDR-based peak calling.

**Configuration:**

Edit `part2_scripts/09_fire_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
REFERENCE_FASTA="/path/to/reference.fa"
FIRE_ROOT="/path/to/FIRE"
REF_NAME="hg38"                # UCSC genome name for track hub
INPUT_SOURCE="kmer_phasing"    # "kmer_phasing" or "nucleosomes"
KEEP_CHROMOSOMES=""            # Regex filter (e.g., "chr[0-9XY]+$")

# Peak calling thresholds (optional, uses FIRE defaults if empty)
MAX_PEAK_FDR=""         # Default: 0.05
MIN_FIRE_FDR=""         # Default: 0.10
MIN_COVERAGE=""         # Default: 4
```

See `config/FIRE_config_readme.md` for detailed configuration options.

**Run:**

```bash
ml singularity/3.6.4
bash part2_scripts/09_fire_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 09_FIRE/
    ├── {sample_name}.fire.bed
    ├── {sample_name}.fire.bw
    ├── {sample_name}.fire.bb
    └── [additional FIRE outputs]
```

---

## Output Structure

Complete output directory structure for each sample:

```
WORKDIR/
└── {sample_name}/
    ├── 01_aligned/
    │   ├── {sample_name}.aligned.bam
    │   └── {sample_name}.aligned.bam.pbi
    ├── 02_sequencing_qc/
    │   ├── {sample_name}.alignmentset.xml
    │   └── {sample_name}.pbqc.report.pdf
    ├── 03_nucleosomes/
    │   ├── {sample_name}.ft.bam
    │   ├── {sample_name}.ft.bam.bai
    │   ├── {sample_name}.ft.bam.pbi
    │   └── {sample_name}.ft.stats.txt
    ├── 04_fiberseq_qc/
    │   ├── qc_summary.txt
    │   ├── plots/
    │   └── {sample_name}.qc_run_summary.txt
    ├── 05_pileups/
    │   ├── {sample_name}.genome_wide.pileup.bed
    │   └── {sample_name}.pileup_stats.txt
    ├── 06_bigwigs/
    │   ├── {sample_name}.genome_wide.coverage.bw
    │   ├── {sample_name}.genome_wide.perc5mc.bw
    │   ├── {sample_name}.genome_wide.nucleosome.bw
    │   └── {sample_name}.bigwig_conversion_summary.txt
    ├── 07_heatmaps/
    │   ├── {sample}.{regions}.coverage.heatmap.png
    │   ├── {sample}.{regions}.perc5mc.heatmap.png
    │   ├── {sample}.{regions}.nucleosome.heatmap.png
    │   ├── {sample}.{regions}.all_marks.heatmap.png
    │   └── {sample}.heatmap_summary.txt
    ├── 07_filter_chroms/           # Part 2 (optional)
    │   ├── {sample_name}.filtered.bam
    │   └── {sample_name}.filtered.bam.bai
    ├── 08_kmer_phasing/            # Part 2 (optional)
    │   └── results/{sample_name}/hiphase/
    │       ├── {sample_name}.bam
    │       └── {sample_name}.bam.bai
    ├── 09_FIRE/                    # Part 2 (optional)
        ├── {sample_name}.fire.bed
        ├── {sample_name}.fire.bw
        └── {sample_name}.fire.bb
```

---
