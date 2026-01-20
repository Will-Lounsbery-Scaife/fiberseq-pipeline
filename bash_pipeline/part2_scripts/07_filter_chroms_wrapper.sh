#!/bin/bash
# Step 7 Wrapper: Filter chromosomes from nucleosome BAMs
# Prerequisites: samtools (load via module or conda)

### CONFIGURE ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback

# Space-separated list of chromosomes to keep (empty string keeps all)
# Default: canonical human chromosomes (chr1-chr22, chrX, chrY)
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
FILTER_THREADS=8


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/07_filter_chroms_main.sh"


mkdir -p "$TMPDIR"
export TMPDIR


tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "Processing: $sample_name"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR" "$CHROMOSOMES" "$FILTER_THREADS"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "Filter Chromosomes complete"
exit 0
