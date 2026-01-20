#!/bin/bash
# Step 5 Wrapper: Pileup to bigWig
# Prerequisites: UCSC tools (bedGraphToBigWig) and GNU Parallel
# On HPC, load modules: ml ucsc-utils gnuParallel

# Load required modules (adjust for your HPC environment)
ml ucsc-utils 2>/dev/null || true
ml gnuParallel 2>/dev/null || true

BEDGRAPHTOBIGWIG="$(command -v bedGraphToBigWig 2>/dev/null)"

### CONFIGURE ###
SAMPLESHEET=""
WORKDIR=""
GENOME_CHROMSIZES=""
TMPDIR="${TMPDIR:-/tmp}"
THREADS=48

# Pileup columns: chrom, start, end, coverage, fire_coverage, score, nuc_coverage, msp_coverage, m6a_coverage, cpg_coverage
declare -A PILEUP_COLUMNS=(["perc_nuc"]=7 ["perc_m6a"]=9 ["perc_cpg"]=10)

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

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR" "$GENOME_CHROMSIZES" "$BEDGRAPHTOBIGWIG" "$COLUMNS_STR" "$THREADS"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

ml purge 2>/dev/null || true
echo "Pipeline complete!"
exit 0
