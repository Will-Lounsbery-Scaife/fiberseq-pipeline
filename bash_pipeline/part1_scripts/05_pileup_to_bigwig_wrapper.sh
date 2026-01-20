#!/bin/bash
# Step 5 Wrapper: Pileup to bigWig
# Prerequisites: UCSC tools (bedGraphToBigWig) and GNU Parallel
# On HPC, load modules: ml ucsc-utils gnuParallel

# Load required modules (adjust for your HPC environment)
ml ucsc-utils 2>/dev/null || true
ml gnuParallel 2>/dev/null || true

BEDGRAPHTOBIGWIG="$(command -v bedGraphToBigWig 2>/dev/null)"

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
GENOME_CHROMSIZES=""  # CONFIGURE: Path to chromosome sizes file
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback
THREADS=48            # Number of parallel jobs for bedGraph/bigWig conversion

# Pileup columns: chrom, start, end, coverage, fire_coverage, score, nuc_coverage, msp_coverage, m6a_coverage, cpg_coverage
# Output will be percentages: mark_coverage / total_coverage (column 4)
declare -A PILEUP_COLUMNS=(["perc_nuc"]=7 ["perc_m6a"]=9 ["perc_cpg"]=10)

### VALIDATION ###
if [[ -z "$SAMPLESHEET" || -z "$WORKDIR" || -z "$GENOME_CHROMSIZES" ]]; then
    echo "ERROR: Required variables not configured. Edit this script and set:"
    echo "  - SAMPLESHEET: Path to your samplesheet TSV"
    echo "  - WORKDIR: Base output directory"
    echo "  - GENOME_CHROMSIZES: Path to chromosome sizes file"
    exit 1
fi

if [[ -z "$BEDGRAPHTOBIGWIG" ]]; then
    echo "ERROR: bedGraphToBigWig not found in PATH"
    echo "On HPC, try: ml ucsc-utils"
    exit 1
fi

### AUTO ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/05_pileup_to_bigwig_main.sh"

mkdir -p "$TMPDIR"
export TMPDIR

# Build columns string
COLUMNS_STR=""
for mark in "${!PILEUP_COLUMNS[@]}"; do
    COLUMNS_STR="${COLUMNS_STR}${mark}:${PILEUP_COLUMNS[$mark]},"
done
COLUMNS_STR="${COLUMNS_STR%,}"

echo "=== Pileup to bigWig Pipeline ==="
echo "Samplesheet: $SAMPLESHEET"
echo "Workdir: $WORKDIR"
echo "Marks: ${!PILEUP_COLUMNS[*]}"
echo "Threads: $THREADS"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR" "$GENOME_CHROMSIZES" "$BEDGRAPHTOBIGWIG" "$COLUMNS_STR" "$THREADS"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

ml purge 2>/dev/null || true
echo "=== Pipeline complete ==="
exit 0
