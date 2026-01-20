#!/bin/bash
# Step 3 Wrapper: Fiber-seq QC
# Prerequisites: conda environment with fiberseq-qc tools

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
FIBERSEQ_QC_ROOT=""   # CONFIGURE: Path to fiberseq-qc installation
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback

### VALIDATION ###
if [[ -z "$SAMPLESHEET" || -z "$WORKDIR" || -z "$FIBERSEQ_QC_ROOT" ]]; then
    echo "ERROR: Required variables not configured. Edit this script and set:"
    echo "  - SAMPLESHEET: Path to your samplesheet TSV"
    echo "  - WORKDIR: Base output directory"
    echo "  - FIBERSEQ_QC_ROOT: Path to fiberseq-qc installation"
    exit 1
fi

### AUTO ###
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/03_fiberseq_qc_main.sh"

mkdir -p "$TMPDIR"
export TMPDIR FIBERSEQ_QC_ROOT

# Override CONDA_DEFAULT_ENV to satisfy runall-qc.tcsh's environment check
export CONDA_DEFAULT_ENV="fiberseq-qc"

echo "=== Fiber-seq QC Pipeline ==="
echo "Samplesheet: $SAMPLESHEET"
echo "Workdir: $WORKDIR"
echo "QC root: $FIBERSEQ_QC_ROOT"

# Note: fiberseq-qc does not support multithreading
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
