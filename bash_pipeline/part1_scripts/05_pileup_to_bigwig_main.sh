#!/bin/bash
# Step 5: Pileup to bedGraph to bigWig

SAMPLE_NAME=$1
SAMPLE_DIR=$2
GENOME_CHROMSIZES=$3
BEDGRAPHTOBIGWIG=$4
COLUMNS_STR=$5
THREADS=${6:-4}

echo "[Step 5] Pileup to bigWig - ${SAMPLE_NAME} ($(date))"

# Parse columns string
declare -A PILEUP_COLUMNS
IFS=',' read -ra PAIRS <<< "$COLUMNS_STR"
for pair in "${PAIRS[@]}"; do
    IFS=':' read -r mark col <<< "$pair"
    PILEUP_COLUMNS[$mark]=$col
done

PILEUP_DIR="$SAMPLE_DIR/04_pileups"
[[ ! -d "$PILEUP_DIR" ]] && { echo "ERROR: Pileups directory not found: $PILEUP_DIR"; exit 1; }

mapfile -t PILEUP_FILES < <(find "$PILEUP_DIR" -name "*.pileup.bed" -type f)
[[ ${#PILEUP_FILES[@]} -eq 0 ]] && { echo "ERROR: No pileup files found"; exit 1; }

echo "Found ${#PILEUP_FILES[@]} pileup file(s)"

BIGWIG_DIR="$SAMPLE_DIR/05_bigwigs"
mkdir -p "$BIGWIG_DIR"

RESULTS_FILE=$(mktemp)
trap "rm -f $RESULTS_FILE" EXIT

# Function to process a single mark from a pileup file
process_mark() {
    local PILEUP_FILE="$1"
    local mark="$2"
    local COL_NUM="$3"
    local PILEUP_DIR="$4"
    local BIGWIG_DIR="$5"
    local GENOME_CHROMSIZES="$6"
    local BEDGRAPHTOBIGWIG="$7"
    local RESULTS_FILE="$8"
    
    local PILEUP_BASENAME=$(basename "$PILEUP_FILE" .pileup.bed)
    local BEDGRAPH_FILE="$PILEUP_DIR/${PILEUP_BASENAME}.${mark}.bedgraph"
    local BIGWIG_FILE="$BIGWIG_DIR/${PILEUP_BASENAME}.${mark}.bw"
    
    # Calculate percentage (mark_coverage / total_coverage) and create bedGraph
    awk -F'\t' -v col="$COL_NUM" \
        '/^#/ {next} 
         $4 > 0 {printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, $col/$4}' \
        "$PILEUP_FILE" > "$BEDGRAPH_FILE"
    
    if [[ ! -s "$BEDGRAPH_FILE" ]]; then
        echo "  Skipping empty bedGraph for $mark" >&2
        return 1
    fi
    
    "$BEDGRAPHTOBIGWIG" "$BEDGRAPH_FILE" "$GENOME_CHROMSIZES" "$BIGWIG_FILE" 2>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "  ERROR: bedGraphToBigWig failed for $mark" >&2
        return 1
    fi
    
    echo "  Created: $(basename $BEDGRAPH_FILE) ($(du -h "$BEDGRAPH_FILE" | cut -f1))"
    echo "  Created: $(basename $BIGWIG_FILE) ($(du -h "$BIGWIG_FILE" | cut -f1))"
    
    echo "$mark" >> "$RESULTS_FILE"
    return 0
}

export -f process_mark

# Build list of jobs
JOBS=()
for PILEUP_FILE in "${PILEUP_FILES[@]}"; do
    [[ ! -s "$PILEUP_FILE" ]] && continue
    for mark in "${!PILEUP_COLUMNS[@]}"; do
        COL_NUM=${PILEUP_COLUMNS[$mark]}
        JOBS+=("$PILEUP_FILE|$mark|$COL_NUM")
    done
done

echo "Running ${#JOBS[@]} conversion jobs with up to $THREADS parallel workers..."

# Run jobs in parallel with GNU parallel
printf '%s\n' "${JOBS[@]}" | parallel -j "$THREADS" --colsep '\|' \
    process_mark {1} {2} {3} "$PILEUP_DIR" "$BIGWIG_DIR" "$GENOME_CHROMSIZES" "$BEDGRAPHTOBIGWIG" "$RESULTS_FILE"

TOTAL_BIGWIGS=$(wc -l < "$RESULTS_FILE" 2>/dev/null || echo 0)
MARKS_SEEN=$(sort -u "$RESULTS_FILE" 2>/dev/null | tr '\n' ' ')

# Generate summary
SUMMARY_FILE="$BIGWIG_DIR/${SAMPLE_NAME}.bigwig_summary.txt"
{
    echo "bigWig Conversion Summary - ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Pileup files: ${#PILEUP_FILES[@]}"
    echo "Parallel threads used: $THREADS"
    echo "Marks processed: $MARKS_SEEN"
    echo "bigWig files created: $TOTAL_BIGWIGS"
    echo ""
    echo "Note: Values are percentages (mark_coverage / total_coverage)"
    echo ""
    echo "bedGraph files (saved in 04_pileups):"
    for bg in $(find "$PILEUP_DIR" -name "*.bedgraph" -type f 2>/dev/null); do
        echo "  $(basename $bg) ($(du -h "$bg" | cut -f1))"
    done
    echo ""
    echo "bigWig files:"
    for bw in $(find "$BIGWIG_DIR" -name "*.bw" -type f 2>/dev/null); do
        echo "  $(basename $bw) ($(du -h "$bw" | cut -f1))"
    done
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  bedGraph files: $TOTAL_BIGWIGS (saved in 04_pileups)"
echo "  bigWig files: $TOTAL_BIGWIGS"
echo "End time: $(date)"

exit 0
