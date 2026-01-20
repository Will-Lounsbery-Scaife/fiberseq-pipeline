#!/bin/bash
# Step 2 Wrapper: Add Nucleosomes
# Prerequisites: conda environment with fibertools-rs and samtools
# Usage: bash 02_add_nucleosomes_wrapper.sh

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
THREADS=48
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback
FT_OPTIONS=""         # Optional: e.g., "--min-ml-score 150"

### VALIDATION ###
if [[ -z "$SAMPLESHEET" || -z "$WORKDIR" ]]; then
    echo "ERROR: Required variables not configured. Edit this script and set:"
    echo "  - SAMPLESHEET: Path to your samplesheet TSV"
    echo "  - WORKDIR: Base output directory"
    exit 1
fi

### AUTO ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/02_add_nucleosomes_main.sh"

mkdir -p "$TMPDIR"
export TMPDIR

echo "=== Add Nucleosomes Pipeline ==="
echo "Samplesheet: $SAMPLESHEET"
echo "Workdir: $WORKDIR"
echo "Threads: $THREADS"
echo "ft version: $(ft --version 2>&1)"

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
