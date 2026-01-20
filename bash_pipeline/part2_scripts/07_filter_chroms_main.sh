#!/bin/bash
# Step 7: Filter chromosomes from nucleosome BAMs

SAMPLE_NAME=$1
SAMPLE_DIR=$2
CHROMOSOMES=$3
FILTER_THREADS=${4:-8}

SOURCE_BAM="${SAMPLE_DIR}/02_nucleosomes/${SAMPLE_NAME}.ft.bam"

OUTPUT_DIR="${SAMPLE_DIR}/07_filter_chroms"
mkdir -p "$OUTPUT_DIR"
FILTERED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.filtered.bam"


# Filter BAM to specified chromosomes (empty CHROMOSOMES keeps all)
if [[ -n "$CHROMOSOMES" ]]; then
    TEMP_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}.tmp.bam"
    samtools view -@ "$FILTER_THREADS" -b -h "$SOURCE_BAM" $CHROMOSOMES > "$TEMP_BAM"
    
    FILTER_EXIT=$?
    if [[ $FILTER_EXIT -ne 0 ]]; then
        echo "ERROR: samtools view filtering failed (exit: $FILTER_EXIT)"
        rm -f "$TEMP_BAM"
        exit 1
    fi
    
    FULL_HEADER="${OUTPUT_DIR}/${SAMPLE_NAME}.full_header.sam"
    CLEAN_HEADER="${OUTPUT_DIR}/${SAMPLE_NAME}.clean_header.sam"
    samtools view -H "$TEMP_BAM" > "$FULL_HEADER"

    grep -v "^@SQ" "$FULL_HEADER" > "$CLEAN_HEADER"

    for chrom in $CHROMOSOMES; do
        grep -P "\tSN:${chrom}\t" "$FULL_HEADER" >> "$CLEAN_HEADER"
    done

    samtools reheader "$CLEAN_HEADER" "$TEMP_BAM" > "$FILTERED_BAM"
    
    rm -f "$TEMP_BAM" "$FULL_HEADER" "$CLEAN_HEADER"

else
    # If no chromosomes provided, just copy with samtools view
    samtools view -@ "$FILTER_THREADS" -b -h "$SOURCE_BAM" > "$FILTERED_BAM"
fi

if [[ ! -f "$FILTERED_BAM" ]]; then
    echo "ERROR: Output BAM was not created successfully."
    exit 1
fi

echo "Indexing filtered BAM..."
samtools index -@ "$FILTER_THREADS" "$FILTERED_BAM"
INDEX_EXIT=$?

if [[ $INDEX_EXIT -ne 0 ]]; then
    echo "ERROR: samtools index failed (exit: $INDEX_EXIT)"
    rm -f "$FILTERED_BAM" "${FILTERED_BAM}.bai"
    exit 1
fi

# Summary
SUMMARY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.filtered_summary.txt"
{
    echo "Filtered Chromosomes Summary"
    echo "Sample: ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Input BAM: $SOURCE_BAM"
    echo "Output BAM: $FILTERED_BAM"
    echo "Chromosomes: ${CHROMOSOMES:-all}"
    echo "Threads: $FILTER_THREADS"
    echo "Input size: $(du -h "$SOURCE_BAM" | cut -f1)"
    echo "Output size: $(du -h "$FILTERED_BAM" | cut -f1)"
    echo "Header Check (SQ count): $(samtools view -H "$FILTERED_BAM" | grep -c "^@SQ")"
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Filtered BAM: $FILTERED_BAM"
echo "End time: $(date)"

exit 0