#!/bin/bash
# Step 4 Wrapper: Create Pileups
# Prerequisites: conda environment with fibertools-rs

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
REGION=""             # Optional: "chr1:1-100000" for subset, empty for genome-wide
ML_THRESHOLD=220      # ML score threshold for m6A/5mC calls
THREADS=48
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback
FT_PILEUP_OPTIONS="-v"    # Additional ft pileup options

### VALIDATION ###
if [[ -z "$SAMPLESHEET" || -z "$WORKDIR" ]]; then
    echo "ERROR: Required variables not configured. Edit this script and set:"
    echo "  - SAMPLESHEET: Path to your samplesheet TSV"
    echo "  - WORKDIR: Base output directory"
    exit 1
fi

### AUTO ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/04_create_pileups_main.sh"

command -v ft &>/dev/null || { echo "ERROR: fibertools (ft) not in PATH"; exit 1; }

mkdir -p "$TMPDIR"
export TMPDIR

echo "=== Create Pileups Pipeline ==="
echo "Samplesheet: $SAMPLESHEET"
echo "Workdir: $WORKDIR"
echo "Region: ${REGION:-genome-wide}"
echo "ML threshold: $ML_THRESHOLD"
echo "ft version: $(ft --version 2>&1)"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "${WORKDIR}/${sample_name}" "$REGION" "$THREADS" "$FT_PILEUP_OPTIONS" "$ML_THRESHOLD"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
