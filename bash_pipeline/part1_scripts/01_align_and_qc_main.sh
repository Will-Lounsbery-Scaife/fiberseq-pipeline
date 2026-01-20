#!/bin/bash
# Step 1: Alignment and Sequencing QC

SAMPLE_NAME=$1
SAMPLE_DIR=$2
RAW_BAM=$3
REFERENCE_FASTA=$4
THREADS=$5
PRESET=$6

echo "[Step 1] Alignment and QC - ${SAMPLE_NAME} ($(date))"

PBMM2="${SMRT_ROOT}/smrtcmds/bin/pbmm2"
PBINDEX="${SMRT_ROOT}/smrtcmds/bin/pbindex"
DATASET="${SMRT_ROOT}/smrtcmds/bin/dataset"
RUNQC="${SMRT_ROOT}/smrtcmds/bin/runqc-reports"

for tool in "$PBMM2" "$PBINDEX" "$DATASET" "$RUNQC"; do
    [[ ! -f "$tool" ]] && { echo "ERROR: Tool not found: $tool"; exit 1; }
done

# Create reference index if it doesn't exist
REFERENCE_INDEX="${REFERENCE_FASTA%.fasta}.mmi"
[[ "$REFERENCE_INDEX" == "$REFERENCE_FASTA" ]] && REFERENCE_INDEX="${REFERENCE_FASTA%.fa}.mmi"

if [[ ! -f "$REFERENCE_INDEX" ]]; then
    echo "Creating reference index: $REFERENCE_INDEX"
    $PBMM2 index --preset $PRESET -j $THREADS "$REFERENCE_FASTA" "$REFERENCE_INDEX" || { echo "ERROR: pbmm2 index failed"; exit 1; }
fi
    
mkdir -p "$SAMPLE_DIR"
ALIGN_DIR="$SAMPLE_DIR/01_aligned"
QC_DIR="$SAMPLE_DIR/01_sequencing_qc"
mkdir -p "$ALIGN_DIR" "$QC_DIR"

# Sequencing QC

# Index raw bam
if [[ ! -f "${RAW_BAM}.pbi" ]]; then
    $PBINDEX -j $THREADS "$RAW_BAM" || { echo "ERROR: pbindex on raw bam failed"; exit 1; }
fi

DATASET_XML="$QC_DIR/${SAMPLE_NAME}.consensusreadset.xml"
QC_REPORT_PDF="$QC_DIR/${SAMPLE_NAME}.sequencing_qc_report.pdf"

echo "Creating ConsensusReadSet..."
$DATASET create --type ConsensusReadSet --force --name "${SAMPLE_NAME}" "$DATASET_XML" "$RAW_BAM" || { echo "ERROR: dataset create failed"; exit 1; }

echo "Generating Sequencing QC Report..."
CURRENT_DIR=$(pwd)
cd "$QC_DIR" || { echo "ERROR: Could not cd to $QC_DIR"; exit 1; }

$RUNQC -b --pdf-report "$QC_REPORT_PDF" "$DATASET_XML" || echo "WARNING: Sequencing QC PDF generation failed (likely missing 'reportlab' module). Proceeding to alignment..."

cd "$CURRENT_DIR"
echo "QC Step Finished"

# Alignment
OUTPUT_BAM="$ALIGN_DIR/${SAMPLE_NAME}.aligned.bam"
echo "Aligning: $RAW_BAM -> $OUTPUT_BAM"

$PBMM2 align --preset $PRESET --sort -j $THREADS --log-level INFO \
    "$REFERENCE_INDEX" "$RAW_BAM" "$OUTPUT_BAM" || { echo "ERROR: pbmm2 align failed"; exit 1; }

echo "Indexing aligned BAM..."
$PBINDEX -j $THREADS "$OUTPUT_BAM" || { echo "ERROR: pbindex on aligned bam failed"; exit 1; }

echo "COMPLETED: ${SAMPLE_NAME}"
echo "  Aligned BAM: $OUTPUT_BAM"
echo "End time: $(date)"

exit 0