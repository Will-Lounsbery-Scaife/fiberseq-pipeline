#!/bin/bash
# Step 6: Create Heatmaps using deepTools

SAMPLE_NAME=$1
SAMPLE_DIR=$2
REGIONS_BED=$3
MATRIX_MODE=$4
REFERENCE_POINT=$5
UPSTREAM=$6
DOWNSTREAM=$7
REGION_BODY_LENGTH=$8
SCALE_UPSTREAM=$9
SCALE_DOWNSTREAM=${10}
COMBINED_COLORMAP=${11}
HEATMAP_HEIGHT=${12}
HEATMAP_WIDTH=${13}
THREADS=${14}
MARKS_STR=${15}

echo "[Step 6] Create Heatmaps - ${SAMPLE_NAME} ($(date))"

# Helper functions
get_label() {
    case "$1" in
        m6a|perc_m6a) echo "% 6mA";;
        cpg|perc_cpg|perc5mc) echo "% 5mC";;
        nuc|perc_nuc|nucleosome) echo "% Nucleosome";;
        coverage) echo "Coverage";;
        *) echo "$1";;
    esac
}

get_colormap() {
    case "$1" in
        m6a|perc_m6a) echo "${M6A_COLORMAP:-Greens}";;
        nuc|perc_nuc|nucleosome) echo "${NUCLEOSOME_COLORMAP:-Blues}";;
        cpg|perc_cpg|perc5mc) echo "${CPG_COLORMAP:-Reds}";;
        *) echo "${COMBINED_COLORMAP:-RdYlBu_r}";;
    esac
}

# Parse marks
IFS=',' read -ra MARKS <<< "$MARKS_STR"
echo "Marks: ${MARKS[*]}"

# Locate bigWig files
BIGWIG_DIR="$SAMPLE_DIR/05_bigwigs"
[[ ! -d "$BIGWIG_DIR" ]] && { echo "ERROR: bigWig directory not found: $BIGWIG_DIR"; exit 1; }

declare -A BIGWIG_FILES
for mark in "${MARKS[@]}"; do
    files=($(ls "${BIGWIG_DIR}"/*.${mark}.bw 2>/dev/null))
    if [[ ${#files[@]} -gt 0 ]]; then
        BIGWIG_FILES[$mark]="${files[0]}"
        echo "Found: $(basename ${files[0]})"
    else
        echo "WARNING: No bigWig for $mark"
    fi
done

[[ ${#BIGWIG_FILES[@]} -eq 0 ]] && { echo "ERROR: No bigWig files found"; exit 1; }

# Setup output
HEATMAP_DIR="$SAMPLE_DIR/06_heatmaps"
mkdir -p "$HEATMAP_DIR"
REGIONS_BASENAME=$(basename "$REGIONS_BED" .bed)

# Track sorted regions from first mark (m6A) to maintain consistent row order across all heatmaps
MASTER_SORTED_BED=""

# Build computeMatrix args based on mode
if [[ "$MATRIX_MODE" == "reference-point" ]]; then
    MATRIX_ARGS="reference-point --referencePoint $REFERENCE_POINT -a $DOWNSTREAM -b $UPSTREAM"
else
    MATRIX_ARGS="scale-regions --regionBodyLength $REGION_BODY_LENGTH -a $SCALE_DOWNSTREAM -b $SCALE_UPSTREAM"
fi

# Process each mark
for mark in "${MARKS[@]}"; do
    BIGWIG_FILE="${BIGWIG_FILES[$mark]}"
    [[ -z "$BIGWIG_FILE" ]] && continue
    
    echo "Processing mark: $mark"
    
    # Determine which regions file to use and whether to filter zeros:
    # - First mark (m6A): use original REGIONS_BED, sort by signal, skip all-zero rows
    # - Subsequent marks: use sorted/filtered regions from m6A to maintain same row order
    if [[ -n "$MASTER_SORTED_BED" && -f "$MASTER_SORTED_BED" ]]; then
        CURRENT_REGIONS="$MASTER_SORTED_BED"
        SORT_FLAG="--sortRegions no"
        SKIP_ZEROS_FLAG=""
        echo "  Using m6A-sorted regions for consistent row order"
    else
        CURRENT_REGIONS="$REGIONS_BED"
        SORT_FLAG=""
        SKIP_ZEROS_FLAG="--skipZeros"
        echo "  Using original regions with --skipZeros (will set sort order for other marks)"
    fi
    
    MATRIX_FILE="$HEATMAP_DIR/${SAMPLE_NAME}.${REGIONS_BASENAME}.${mark}.matrix.gz"
    HEATMAP_PNG="$HEATMAP_DIR/${SAMPLE_NAME}.${REGIONS_BASENAME}.${mark}.heatmap.png"
    HEATMAP_PDF="$HEATMAP_DIR/${SAMPLE_NAME}.${REGIONS_BASENAME}.${mark}.heatmap.pdf"
    SORTED_REGIONS_FILE="${HEATMAP_PNG%.png}.sorted_regions.bed"
    
    COLORMAP=$(get_colormap "$mark")
    LABEL=$(get_label "$mark")
    
    # Compute matrix using appropriate regions file
    $DEEPTOOLS_PREFIX computeMatrix $MATRIX_ARGS \
        -S "$BIGWIG_FILE" -R "$CURRENT_REGIONS" -p $THREADS \
        --missingDataAsZero $SKIP_ZEROS_FLAG -o "$MATRIX_FILE" || { echo "ERROR: computeMatrix failed for $mark"; continue; }
    
    # Plot heatmaps (with optional sort flag for non-m6A marks)
    $DEEPTOOLS_PREFIX plotHeatmap -m "$MATRIX_FILE" -o "$HEATMAP_PNG" \
        --colorMap "$COLORMAP" --heatmapHeight $HEATMAP_HEIGHT --heatmapWidth $HEATMAP_WIDTH \
        --plotTitle "${SAMPLE_NAME} - ${LABEL}" --legendLocation best \
        --outFileSortedRegions "$SORTED_REGIONS_FILE" $SORT_FLAG || echo "WARNING: plotHeatmap PNG failed"
    
    $DEEPTOOLS_PREFIX plotHeatmap -m "$MATRIX_FILE" -o "$HEATMAP_PDF" \
        --colorMap "$COLORMAP" --heatmapHeight $HEATMAP_HEIGHT --heatmapWidth $HEATMAP_WIDTH \
        --plotTitle "${SAMPLE_NAME} - ${LABEL}" --legendLocation best $SORT_FLAG || echo "WARNING: plotHeatmap PDF failed"
    
    # Capture sorted regions from first mark (m6A) to use for subsequent marks
    if [[ -z "$MASTER_SORTED_BED" && -f "$SORTED_REGIONS_FILE" ]]; then
        MASTER_SORTED_BED="$SORTED_REGIONS_FILE"
        echo "  Captured m6A sort order: $(basename $MASTER_SORTED_BED)"
    fi
    
    echo "  Created: $(basename $HEATMAP_PNG)"
done

# Generate summary
SUMMARY_FILE="$HEATMAP_DIR/${SAMPLE_NAME}.heatmap_summary.txt"
{
    echo "Heatmap Summary - ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Regions: $REGIONS_BED ($(wc -l < "$REGIONS_BED") regions)"
    echo "Mode: $MATRIX_MODE"
    echo "Marks: ${!BIGWIG_FILES[*]}"
    echo ""
    echo "Output files:"
    for f in $(find "$HEATMAP_DIR" -name "*.png" -o -name "*.pdf" | sort); do
        echo "  $(basename $f)"
    done
} > "$SUMMARY_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Heatmaps: $HEATMAP_DIR"
echo "End time: $(date)"

exit 0
