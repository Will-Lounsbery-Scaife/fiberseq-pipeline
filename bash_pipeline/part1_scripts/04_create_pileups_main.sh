#!/bin/bash
# Step 4: Create Pileups using fibertools

SAMPLE_NAME=$1
SAMPLE_DIR=$2
REGION=$3
THREADS=$4
FT_PILEUP_OPTIONS=$5
ML_THRESHOLD=$6

echo "[Step 4] Create Pileups - ${SAMPLE_NAME} ($(date))"

# Locate input files
INPUT_BAM="$SAMPLE_DIR/02_nucleosomes/${SAMPLE_NAME}.ft.bam"
[[ ! -f "$INPUT_BAM" ]] && { echo "ERROR: Nucleosome BAM not found: $INPUT_BAM"; exit 1; }

if [[ ! -f "$INPUT_BAM.bai" && ! -f "$INPUT_BAM.pbi" ]]; then
    echo "ERROR: No BAM index found (.bai or .pbi)"
    exit 1
fi

# Setup output
PILEUP_DIR="$SAMPLE_DIR/04_pileups"
    mkdir -p "$PILEUP_DIR"

if [[ -n "$REGION" ]]; then
    region_safe=$(echo "$REGION" | tr ':,-' '_')
    OUTPUT_PILEUP="$PILEUP_DIR/${SAMPLE_NAME}.${region_safe}.pileup.bed"
else
    OUTPUT_PILEUP="$PILEUP_DIR/${SAMPLE_NAME}.genome_wide.pileup.bed"
fi

# Build ft pileup command using array for safety
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_PILEUP"
[[ -n "$REGION" ]] && echo "Region: $REGION" || echo "Region: genome-wide"

ft_args=("pileup" "-t" "$THREADS" "-m" "-c")
[[ -n "$FT_PILEUP_OPTIONS" ]] && read -ra opts <<< "$FT_PILEUP_OPTIONS" && ft_args+=("${opts[@]}")
[[ -n "$ML_THRESHOLD" ]] && ft_args+=("--ml" "$ML_THRESHOLD")
ft_args+=("$INPUT_BAM")
[[ -n "$REGION" ]] && ft_args+=("$REGION")
ft_args+=("-o" "$OUTPUT_PILEUP")

echo "Running: ft ${ft_args[*]}"
ft "${ft_args[@]}" || { echo "ERROR: ft pileup failed"; exit 1; }

# Verify output
[[ ! -f "$OUTPUT_PILEUP" ]] && { echo "ERROR: Pileup file not created"; exit 1; }

num_lines=$(wc -l < "$OUTPUT_PILEUP")
echo "Pileup: $(du -h "$OUTPUT_PILEUP" | cut -f1), $num_lines lines"

[[ "$num_lines" -eq 0 ]] && echo "WARNING: Pileup file is empty"

# Generate summary statistics
STATS_FILE="$PILEUP_DIR/${SAMPLE_NAME}.pileup_stats.txt"
{
    echo "Pileup Summary - ${SAMPLE_NAME}"
    echo "Date: $(date)"
    echo "Input: $INPUT_BAM"
    echo "Output: $OUTPUT_PILEUP"
    echo "Region: ${REGION:-genome-wide}"
    echo "File size: $(du -h "$OUTPUT_PILEUP" | cut -f1)"
    echo "Total lines: $num_lines"
echo ""
    echo "Positions per chromosome:"
    cut -f1 "$OUTPUT_PILEUP" | sort | uniq -c | sort -k2 -V | awk '{printf "  %-20s %s\n", $2, $1}'
echo ""
    echo "Column stats (first 10000 lines):"
    head -n 10000 "$OUTPUT_PILEUP" | awk -F'\t' '
    NR==1 {next}
    {
        for (i=4; i<=NF && i<=10; i++) {
            if ($i ~ /^[0-9.]+$/) {
                sum[i]+=$i; count[i]++
                if (count[i]==1 || $i<min[i]) min[i]=$i
                if (count[i]==1 || $i>max[i]) max[i]=$i
            }
        }
    }
    END {
        n[4]="coverage"; n[5]="fire"; n[6]="score"; n[7]="nuc"; n[8]="msp"; n[9]="m6a"; n[10]="cpg"
        for (i=4; i<=10; i++) if (count[i]>0) printf "  %s: mean=%.2f min=%.2f max=%.2f\n", n[i], sum[i]/count[i], min[i], max[i]
    }'
} > "$STATS_FILE"

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Pileup: $OUTPUT_PILEUP"
echo "End time: $(date)"

exit 0
