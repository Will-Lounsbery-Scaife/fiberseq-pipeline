# Fiberseq Pipeline

A complete shell-based pipeline for processing and analyzing PacBio HiFi fiber-seq data, from raw BAM files to publication-quality heatmaps.

---

## Biological Background

### What is Fiber-seq?

Fiber-seq is a single-molecule sequencing method that allows researchers to map chromatin structure and accessibility with high resolution on individual DNA fibers. The technique works by treating cell nuclei with a non-specific N6-adenine methyltransferase (such as Hia5). This enzyme acts like a "stencil," methylating adenine bases within open, accessible regions of the chromatin while leaving DNA wrapped around nucleosomes or bound by transcription factors unmodified. The genomic DNA is then sequenced using PacBio long-read sequencing, which can directly detect these N6-methyladenine (m6A) marks. By reading the pattern of methylated bases along long strands of DNA, the method creates a digital footprint of the chromatin architecture, showing exactly where nucleosomes were positioned and where the DNA was open and active on that specific molecule.

### Why Use Fiber-seq?

Fiber-seq is useful due to its ability to capture the exact regulatory state of single DNA molecules, rather than an average signal from millions of cells (like ATAC-seq). This enables observing how multiple regulatory elements, such as promoters and enhancers, coordinate their activity on the same strand of DNA (co-actuation). It's particularly useful for resolving chromatin structure in complex, repetitive regions of the genome. It can also be useful for phasing regulatory information to specific maternal or paternal haplotypes.

### What is FIRE?

FIRE ("Fiber-seq Inferred Regulatory Elements") is a computational peak-calling and annotation framework that takes Fiber-seq single-molecule data and classifies which methyltransferase-sensitive patches (MSPs; stretches of m6A-labeled, exposed DNA on a single read) correspond to bona fide regulatory element "actuation" versus ordinary internucleosomal linker-like accessibility. Concretely, FIRE computes a set of per-MSP features along each long molecule and applies a semi-supervised machine-learning classifier (implemented with Mokapot + XGBoost) to assign each candidate MSP an estimated precision/probability of being a true regulatory element; high-precision calls are labeled FIREs, while lower-precision MSPs are treated as background accessibility typical of linker DNA. FIRE then aggregates single-molecule calls across reads to generate genome-wide FIRE "peaks" with an empirical FDR strategy, and it naturally yields a quantitative accessibility metric often described as percent chromatin actuation (the fraction of molecules showing a FIRE at a locus), enabling allele/haplotype-resolved and heterogeneity-aware regulatory element maps from Fiber-seq alone.

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Pipeline Structure](#pipeline-structure)
5. [Configuration](#configuration)
6. [Running the Pipeline](#running-the-pipeline)
   - [Part 1: Core Processing](#part-1-core-processing)
     - [Step 1: Alignment and Sequencing QC](#step-1-alignment-and-sequencing-qc)
     - [Step 2: Add Nucleosomes](#step-2-add-nucleosomes)
     - [Step 3: Fiber-seq QC](#step-3-fiber-seq-qc)
     - [Step 4: Create Pileups](#step-4-create-pileups)
     - [Step 5: Pileup to bigWig](#step-5-pileup-to-bigwig)
     - [Step 6: Create Heatmaps](#step-6-create-heatmaps)
   - [Part 2: Advanced Analysis](#part-2-advanced-analysis)
     - [Step 7: Filter Chromosomes](#step-7-filter-chromosomes)
     - [Step 8: K-mer Variant Phasing](#step-8-k-mer-variant-phasing)
     - [Step 9: FIRE Peak Calling](#step-9-fire-peak-calling)
   - [Part 3: TF Footprinting Analysis](#part-3-tf-footprinting-analysis)
     - [Step 10: Footprint Detection](#step-10-footprint-detection)
     - [Step 11: Fiber Classification](#step-11-fiber-classification)
     - [Step 12: Footprint Visualization](#step-12-footprint-visualization)
7. [Output Structure](#output-structure)
8. [Troubleshooting](#troubleshooting)

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

**Part 2: Advanced Analysis (Optional)**
7. **Filter Chromosomes**: Filter BAMs to specific chromosomes (e.g., canonical human chromosomes)
8. **K-mer Variant Phasing**: Phase variants using k-mer-based approach (Snakemake workflow)
9. **FIRE Peak Calling**: Identify Fiber-seq Inferred Regulatory Elements using FDR-based peak calling (Snakemake workflow)

---

## Prerequisites

### Required Software

- **PacBio SMRT Tools** (pbmm2, pbindex, dataset, runqc-reports)
- **fibertools** (ft)
- **fiberseq-qc**
- **bedGraphToBigWig** (UCSC tools)
- **deepTools** (computeMatrix, plotHeatmap)
- **samtools** (optional, for indexing)
- **tcsh** (required for fiberseq-qc)

### Required Files

- Raw unaligned BAM files (PacBio HiFi reads)
- Reference genome FASTA file
- Genome chromosome sizes file
- BED file with regions of interest (for heatmaps)

---

## Installation

### 1. Install SMRT Tools

```bash
# Download and install SMRT Tools
# See: https://www.pacb.com/support/software-downloads/

# Example installation path
SMRT_ROOT="/path/to/smrttools"
```

### 2. Install fibertools (conda environment)

```bash
conda create -n fibertools_env -c bioconda fibertools
conda activate fibertools_env
```

### 3. Install fiberseq-qc

```bash
git clone https://github.com/fiberseq/fiberseq-qc.git
# Follow installation instructions in the repository
```

### 4. Install bedGraphToBigWig

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig
# Move to a directory in your PATH or note the full path
```

### 5. Install deepTools

```bash
conda create -n deeptools_env python=3.9
conda activate deeptools_env
conda install -c bioconda deeptools
```

---

## Pipeline Structure

```
pipeline/
├── config/
│   ├── samples.tsv              # Sample metadata (default)
│   ├── samples_example.tsv      # Template sample manifest
│   ├── samples_HIV.tsv          # HIV genome samples
│   ├── samples_chr21.tsv        # Chromosome 21 only (testing)
│   └── FIRE_config_readme.md    # FIRE configuration documentation
├── part1_scripts/               # Steps 1-6: Core fiber-seq processing
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
├── part2_scripts/               # Steps 7-9: Advanced analysis
│   ├── 07_filter_chroms_wrapper.sh
│   ├── 07_filter_chroms_main.sh
│   ├── 08_kmer_phasing_wrapper.sh
│   ├── 08_kmer_phasing_main.sh
│   ├── 09_fire_wrapper.sh
│   ├── 09_fire_main.sh
│   └── config/                  # Generated YAML configs
├── part3_scripts/               # Steps 10-12: TF Footprinting
│   ├── 10_footprint_detection_wrapper.sh
│   ├── 10_footprint_detection_main.sh
│   ├── 11_fiber_classification_wrapper.sh
│   ├── 11_fiber_classification_main.sh
│   ├── 12_footprint_viz_wrapper.sh
│   ├── 12_footprint_viz_main.sh
│   └── parse_footprints.py      # Python script for fiber classification
└── README.md
```

---

## Configuration

### 1. Create Sample Sheet

Create a tab-delimited file at `config/samples.tsv`:

```tsv
sample_name	raw_bam_path	tissue	sequencer
SA271_deep_PBG4150	/path/to/SA271_deep_PBG4150.bam	lymphocyte	revio
SA273_deep_1_PBG4134	/path/to/SA273_deep_1_PBG4134.bam	lymphocyte	revio
SA273_deep_2_PBG4134	/path/to/SA273_deep_2_PBG4134.bam	lymphocyte	revio
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

### Part 2: Advanced Analysis

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

### Part 3: TF Footprinting Analysis

Part 3 provides transcription factor (TF) footprinting analysis using fiber-seq data. It detects where TFs are actively bound to DNA at the single-molecule level by identifying "footprints" (protected regions within accessible chromatin).

**Prerequisites for Part 3:**
- Completed Step 2 (nucleosome calling with fibertools)
- TF motifs BED file (user-provided, e.g., from FIMO motif scanning)
- Python 3 (standard library only)
- deepTools (via Singularity container)

---

#### Step 10: Footprint Detection

**Purpose:** Detect TF footprints at motif sites using `ft footprint`. Optionally scan for motifs in ChIP peaks using FIMO.

**Configuration:**

Edit `part3_scripts/10_footprint_detection_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
REFERENCE_FASTA="/path/to/reference.fa"

# Motif source - provide pre-computed motifs BED file
MOTIFS_BED="/path/to/tf_motifs.bed"

# Optional: Scan motifs from ChIP peaks using FIMO
SCAN_MOTIFS="false"    # Set to "true" to enable FIMO
PEAKS_BED=""           # ChIP-seq peaks BED (if SCAN_MOTIFS=true)
MOTIF_FILE=""          # MEME format motif file (if SCAN_MOTIFS=true)
FIMO_THRESH="1e-4"     # FIMO p-value threshold

THREADS=48
```

**Run:**

```bash
conda activate fibertools_env
module load samtools/1.21
# If using FIMO: module load meme/5.5.1 bedtools
bash part3_scripts/10_footprint_detection_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 10_footprinting/
    ├── {sample_name}.footprints.tsv
    ├── {sample_name}.footprint.yaml
    ├── {sample_name}.footprint_summary.txt
    └── motif_scan/                    # Only if SCAN_MOTIFS=true
        ├── peaks.fa
        ├── fimo_out/
        └── scanned_motifs.bed
```

---

#### Step 11: Fiber Classification

**Purpose:** Parse footprinting results to classify fibers as TF-bound (footprinted) or accessible (unbound), then create separate BAM files for each category.

**Configuration:**

Edit `part3_scripts/11_fiber_classification_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
THREADS=48
```

**Run:**

```bash
module load samtools/1.21
bash part3_scripts/11_fiber_classification_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 11_fiber_classification/
    ├── fiber_lists/
    │   ├── footprinted_fibers.txt
    │   └── accessible_fibers.txt
    ├── filtered_bams/
    │   ├── footprinted.bam
    │   ├── footprinted.bam.bai
    │   ├── accessible.bam
    │   └── accessible.bam.bai
    └── {sample_name}.classification_summary.txt
```

---

#### Step 12: Footprint Visualization

**Purpose:** Generate CPM-normalized m6A bigWig files and deepTools profile plots comparing TF-bound vs accessible fibers.

**Configuration:**

Edit `part3_scripts/12_footprint_viz_wrapper.sh`:

```bash
SAMPLESHEET="/path/to/config/samples.tsv"
WORKDIR="/path/to/output/directory"
GENOME_CHROMSIZES="/path/to/reference.chrom.sizes"
MOTIFS_BED="/path/to/tf_motifs.bed"      # Same as Step 10
SINGULARITY_IMAGE="/path/to/deeptools_container.sif"

# BigWig parameters
ML_THRESHOLD=220       # m6A ML score threshold
THREADS=48

# Plot parameters
UPSTREAM=1000          # bp upstream of motif center
DOWNSTREAM=1000        # bp downstream of motif center
```

**Run:**

```bash
conda activate fibertools_env
module load ucsc-utils/2023-10-17
module load singularity/3.6.4
bash part3_scripts/12_footprint_viz_wrapper.sh
```

**Output:**
```
WORKDIR/{sample_name}/
└── 12_footprint_viz/
    ├── pileups/
    │   ├── footprinted.pileup.bed
    │   ├── accessible.pileup.bed
    │   ├── footprinted.m6a.bedgraph
    │   └── accessible.m6a.bedgraph
    ├── bigwigs/
    │   ├── footprinted.m6a_cpm.bw
    │   └── accessible.m6a_cpm.bw
    ├── plots/
    │   ├── {sample_name}.footprint_matrix.gz
    │   ├── {sample_name}.footprint_profile.svg
    │   ├── {sample_name}.footprint_profile.png
    │   ├── {sample_name}.footprint_profile.tab
    │   └── {sample_name}.footprint_heatmap.svg
    └── {sample_name}.footprint_viz_summary.txt
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
    │   ├── {sample_name}.fire.bed
    │   ├── {sample_name}.fire.bw
    │   └── {sample_name}.fire.bb
    ├── 10_footprinting/            # Part 3 (optional)
    │   ├── {sample_name}.footprints.tsv
    │   └── {sample_name}.footprint_summary.txt
    ├── 11_fiber_classification/    # Part 3 (optional)
    │   ├── fiber_lists/
    │   └── filtered_bams/
    └── 12_footprint_viz/           # Part 3 (optional)
        ├── bigwigs/
        └── plots/
```

---

## Troubleshooting

### General Issues

**Scripts not executable:**
```bash
chmod +x part1_scripts/*.sh part2_scripts/*.sh
```

**Sample directory not found:**
- Ensure previous steps completed successfully
- Check that WORKDIR is the same across all steps
- Verify sample names in samples.tsv match output directories

### Step 1: Alignment and QC

**pbmm2 not found:**
```bash
# Verify SMRT_ROOT is set correctly
ls $SMRT_ROOT/smrtcmds/bin/pbmm2

# Update SMRT_ROOT in wrapper script
```

**Reference index creation fails:**
- Check reference FASTA file exists and is readable
- Ensure sufficient disk space
- Verify file permissions

### Step 2: Add Nucleosomes

**fibertools not found:**
```bash
# Activate conda environment
conda activate fibertools_env

# Verify installation
ft --version
```

**Input BAM not found:**
- Ensure Step 1 completed successfully
- Check that aligned BAM exists in `01_aligned/`

### Step 3: Fiber-seq QC

**tcsh not found:**
```bash
# Install tcsh
# Ubuntu/Debian:
sudo apt-get install tcsh

# CentOS/RHEL:
sudo yum install tcsh
```

**fiberseq-qc script not found:**
- Verify FIBERSEQ_QC_ROOT points to correct directory
- Check that `src/runall-qc.tcsh` exists in that directory

### Step 4: Create Pileups

**Empty pileup file:**
- Verify input BAM has fibertools annotations from Step 2
- Check that regions in BED file match chromosome names in BAM
- Ensure there's data in the specified regions

### Step 5: Pileup to bigWig

**bedGraphToBigWig fails:**
```bash
# Common issue: chromosome name mismatch
# Check pileup chromosomes
cut -f1 pileup.bed | sort -u

# Check chrom.sizes chromosomes
cut -f1 reference.chrom.sizes | sort -u

# They must match (both use "chr1" or both use "1")
```

**Empty bigWig files:**
- Verify PILEUP_COLUMNS are correct for your pileup format
- Check pileup file has data: `head pileup.bed`

### Step 6: Create Heatmaps

**computeMatrix fails:**
```bash
# Check chromosome naming consistency
bigWigInfo file.bw | grep chrom
cut -f1 regions.bed | sort -u

# Fix if needed (add/remove 'chr' prefix)
sed 's/^/chr/' regions.bed > regions_with_chr.bed
```

**deepTools not found:**
```bash
# Activate correct environment
conda activate deeptools_env

# Verify installation
computeMatrix --version
plotHeatmap --version
```

### Memory Issues

**Large files causing memory errors:**
- Increase available memory/RAM
- Process samples one at a time
- Reduce number of regions
- Decrease upstream/downstream distances

### Disk Space Issues

**Running out of disk space:**
```bash
# Check disk usage
df -h $WORKDIR
df -h $TMPDIR

# Clean up intermediate files if needed
# (Only after confirming outputs are correct)
```

---

## Support and References

### Documentation
- **fibertools**: https://fiberseq.github.io/fibertools/
- **PacBio SMRT Tools**: https://www.pacb.com/support/software-downloads/
- **deepTools**: https://deeptools.readthedocs.io/
- **fiberseq-qc**: https://github.com/fiberseq/fiberseq-qc

### Citation
If you use this pipeline, please cite the relevant tools:
- PacBio SMRT Tools
- fibertools
- deepTools
- UCSC Genome Browser tools

---

## Quick Start Example

```bash
# 1. Clone or copy the pipeline
cd /path/to/your/project
# Copy pipeline scripts to your project directory

# 2. Create sample sheet
cat > config/samples.tsv << EOF
sample_name	raw_bam_path	tissue	sequencer
sample1	/path/to/sample1.bam	tissue1	revio
EOF

# 3. Configure wrapper scripts
# Edit the CUSTOMIZE section in each wrapper script with your paths

# 4. Run Part 1 pipeline steps sequentially

# Step 1: Alignment
bash part1_scripts/01_align_and_qc_wrapper.sh

# Step 2-4: Nucleosomes, QC, Pileups
conda activate fibertools_env
bash part1_scripts/02_add_nucleosomes_wrapper.sh
bash part1_scripts/03_fiberseq_qc_wrapper.sh
bash part1_scripts/04_create_pileups_wrapper.sh
conda deactivate

# Step 5: bigWig conversion
bash part1_scripts/05_pileup_to_bigwig_wrapper.sh

# Step 6: Heatmaps
conda activate deeptools_env
bash part1_scripts/06_create_heatmaps_wrapper.sh
conda deactivate

# 5. (Optional) Run Part 2 advanced analysis
ml singularity/3.6.4
bash part2_scripts/07_filter_chroms_wrapper.sh
bash part2_scripts/08_kmer_phasing_wrapper.sh
bash part2_scripts/09_fire_wrapper.sh
```

---

**Pipeline Version:** 1.0
**Last Updated:** 2026