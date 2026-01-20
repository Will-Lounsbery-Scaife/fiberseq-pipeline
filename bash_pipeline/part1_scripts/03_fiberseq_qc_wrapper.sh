#!/bin/bash
# Step 3 Wrapper: Fiber-seq QC
# Prerequisites: conda environment with fiberseq-qc tools

### CONFIGURE ###
SAMPLESHEET=""
WORKDIR=""
FIBERSEQ_QC_ROOT=""
TMPDIR="${TMPDIR:-/tmp}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/03_fiberseq_qc_main.sh"

mkdir -p "$TMPDIR"
export TMPDIR FIBERSEQ_QC_ROOT

export CONDA_DEFAULT_ENV="fiberseq-qc"


tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "Processing $sample_name"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
