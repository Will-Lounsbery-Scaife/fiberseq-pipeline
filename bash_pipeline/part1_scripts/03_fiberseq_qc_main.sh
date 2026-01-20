#!/bin/bash
# Step 3: Fiber-seq QC

SAMPLE_NAME=$1
SAMPLE_DIR=$2

echo "[Step 3] Fiber-seq QC - ${SAMPLE_NAME} ($(date))"


RUNALL_QC="${FIBERSEQ_QC_ROOT}/src/runall-qc.tcsh"

# Locate input files
INPUT_BAM="$SAMPLE_DIR/02_nucleosomes/${SAMPLE_NAME}.ft.bam"

# Setup output
QC_DIR="$SAMPLE_DIR/03_fiberseq_qc"
    mkdir -p "$QC_DIR"

# Run fiberseq-qc
echo "Input BAM: $INPUT_BAM"
echo "Output dir: $QC_DIR"

cd "$QC_DIR"
tcsh "$RUNALL_QC" . "$SAMPLE_NAME" "$INPUT_BAM"
QC_EXIT_CODE=$?
cd - >/dev/null

[[ $QC_EXIT_CODE -ne 0 ]] && { echo "ERROR: fiberseq-qc failed (exit code: $QC_EXIT_CODE)"; exit 1; }

# Generate summary
SUMMARY_FILE="$QC_DIR/${SAMPLE_NAME}.qc_run_summary.txt"
{
    echo "Fiber-seq QC Summary"
    echo "Sample: ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Input BAM: $INPUT_BAM"
    echo "Output directory: $QC_DIR"
    echo "Exit code: $QC_EXIT_CODE"
    echo ""
    echo "Files generated: $(find "$QC_DIR" -type f | wc -l)"
    echo ""
    ls -lh "$QC_DIR" | tail -n +2
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  QC output: $QC_DIR"
echo "End time: $(date)"

exit 0
