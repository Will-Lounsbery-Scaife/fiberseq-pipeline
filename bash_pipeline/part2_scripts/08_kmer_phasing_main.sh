#!/bin/bash
# Step 8: k-mer Variant Phasing Main Script
# Runs snakemake k-mer-variant-phasing workflow for a single sample

SAMPLE_NAME=$1
OUTPUT_DIR=$2
CONFIG_FILE=$3
PIXI_MANIFEST=$4
PROFILE_PATH=$5
RESUME=${6:-"false"}

echo "[k-mer Phasing] ${SAMPLE_NAME} ($(date))"
echo "  Resume mode: $RESUME"

# Check pixi is available
command -v pixi &>/dev/null || { echo "ERROR: pixi not found in PATH"; exit 1; }

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "ERROR: Cannot cd to $OUTPUT_DIR"; exit 1; }

# Build snakemake arguments
SNAKEMAKE_ARGS=(
    --profile "$PROFILE_PATH"
    --configfile "$CONFIG_FILE"
    --verbose
)
if [[ "$RESUME" == "true" ]]; then
    SNAKEMAKE_ARGS+=(--rerun-incomplete)
fi

# Run pixi snakemake with --manifest-path to use the k-mer-variant-phasing environment
pixi run --manifest-path "$PIXI_MANIFEST" snakemake "${SNAKEMAKE_ARGS[@]}"

SNAKEMAKE_EXIT=$?

# Generate summary
SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.kmer_phasing_summary.txt"
{
    echo "k-mer Variant Phasing Summary"
    echo "Sample: ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Config: $CONFIG_FILE"
    echo "Exit code: $SNAKEMAKE_EXIT"
    echo ""
    echo "Output files:"
    find "$OUTPUT_DIR" -type f -name "*.vcf*" -o -name "*.bam" -o -name "*.tsv" 2>/dev/null | head -20
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Output: $OUTPUT_DIR"
echo "End time: $(date)"

exit 0
