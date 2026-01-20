#!/bin/bash
# Step 9: FIRE Main Script
# Runs snakemake FIRE workflow for a single sample

SAMPLE_NAME=$1
OUTPUT_DIR=$2
CONFIG_FILE=$3
PIXI_MANIFEST=$4
PROFILE_PATH=$5

echo "[FIRE] ${SAMPLE_NAME} ($(date))"

command -v pixi &>/dev/null || { echo "ERROR: pixi not found in PATH"; exit 1; }
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "ERROR: Cannot cd to $OUTPUT_DIR"; exit 1; }



# Run pixi fire task with --manifest-path to use the FIRE environment
pixi run --manifest-path "$PIXI_MANIFEST" fire \
    --profile "$PROFILE_PATH" \
    --configfile "$CONFIG_FILE"

SNAKEMAKE_EXIT=$?

# Generate summary
SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.fire_summary.txt"
{
    echo "FIRE Summary"
    echo "Sample: ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Config: $CONFIG_FILE"
    echo "Exit code: $SNAKEMAKE_EXIT"
    echo ""
    echo "Output files:"
    find "$OUTPUT_DIR" -type f \( -name "*.bed" -o -name "*.bw" -o -name "*.bb" -o -name "*.cram" \) 2>/dev/null | head -30
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Output: $OUTPUT_DIR"
echo "End time: $(date)"

exit 0
