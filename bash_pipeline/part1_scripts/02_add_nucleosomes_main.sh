#!/bin/bash
# Step 2: Add Nucleosomes (MSPs) using fibertools

SAMPLE_NAME=$1
SAMPLE_DIR=$2
THREADS=$3
FT_OPTIONS=$4

echo "[Step 2] Add Nucleosomes - ${SAMPLE_NAME} ($(date))"

# Locate input files
INPUT_BAM="$SAMPLE_DIR/01_aligned/${SAMPLE_NAME}.aligned.bam"

# Setup output
NUCLEOSOME_DIR="$SAMPLE_DIR/02_nucleosomes"
    mkdir -p "$NUCLEOSOME_DIR"
OUTPUT_BAM="$NUCLEOSOME_DIR/${SAMPLE_NAME}.ft.bam"

# Run fibertools add-nucleosomes
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_BAM"

FT_CMD="ft add-nucleosomes -t $THREADS"
[[ -n "$FT_OPTIONS" ]] && FT_CMD="$FT_CMD $FT_OPTIONS"
FT_CMD="$FT_CMD $INPUT_BAM $OUTPUT_BAM"

echo "Running: $FT_CMD"
eval $FT_CMD || { echo "ERROR: ft add-nucleosomes failed"; exit 1; }

# Index output BAM
echo "Indexing output BAM..."
if command -v samtools &>/dev/null; then
    samtools index -@ $THREADS "$OUTPUT_BAM" && echo "Created .bai index"
    fi

if command -v pbindex &>/dev/null; then
    pbindex "$OUTPUT_BAM" && echo "Created .pbi index"
elif [[ -n "$SMRT_ROOT" && -f "$SMRT_ROOT/smrtcmds/bin/pbindex" ]]; then
    "$SMRT_ROOT/smrtcmds/bin/pbindex" "$OUTPUT_BAM" && echo "Created .pbi index"
        fi

# Generate summary stats
STATS_FILE="$NUCLEOSOME_DIR/${SAMPLE_NAME}.ft.stats.txt"
{
    echo "Fibertools Add-Nucleosomes Summary"
    echo "Sample: ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Input BAM: $INPUT_BAM"
    echo "Output BAM: $OUTPUT_BAM"
    echo ""
    if command -v samtools &>/dev/null; then
        echo "Total reads: $(samtools view -c "$OUTPUT_BAM" 2>/dev/null || echo 'N/A')"
    fi
    echo "Input size: $(du -h "$INPUT_BAM" | cut -f1)"
    echo "Output size: $(du -h "$OUTPUT_BAM" | cut -f1)"
} > "$STATS_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Output BAM: $OUTPUT_BAM"
echo "End time: $(date)"

exit 0
