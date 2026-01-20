#!/bin/bash
# Step 2 Wrapper: Add Nucleosomes
# Prerequisites: conda environment with fibertools-rs and samtools
# Usage: bash 02_add_nucleosomes_wrapper.sh

### CONFIGURE ###
SAMPLESHEET=""
WORKDIR=""
THREADS=48
TMPDIR="${TMPDIR:-/tmp}"
FT_OPTIONS=""


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/02_add_nucleosomes_main.sh"

mkdir -p "$TMPDIR"
export TMPDIR

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR" "$THREADS" "$FT_OPTIONS"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
