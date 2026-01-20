#!/bin/bash
# Step 4 Wrapper: Create Pileups
# Prerequisites: conda environment with fibertools-rs

### CONFIGURE ###
SAMPLESHEET=""
WORKDIR=""
REGION=""
ML_THRESHOLD=220
THREADS=48
TMPDIR="${TMPDIR:-/tmp}"
FT_PILEUP_OPTIONS="-v"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/04_create_pileups_main.sh"

command -v ft &>/dev/null || { echo "ERROR: fibertools (ft) not in PATH"; exit 1; }

mkdir -p "$TMPDIR"
export TMPDIR


tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "${WORKDIR}/${sample_name}" "$REGION" "$THREADS" "$FT_PILEUP_OPTIONS" "$ML_THRESHOLD"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
