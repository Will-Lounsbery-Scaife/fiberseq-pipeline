#!/bin/bash
# Step 6 Wrapper: Create Heatmaps
# Prerequisites: Singularity/Apptainer with deepTools container

# Load singularity module (adjust for your HPC)
module load singularity 2>/dev/null || ml singularity 2>/dev/null || true

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
WORKDIR=""            # CONFIGURE: Base output directory (same as step 1)
REGIONS_BED=""        # CONFIGURE: BED file with regions for heatmaps (e.g., promoters)
SINGULARITY_IMAGE=""  # CONFIGURE: Path to deepTools Singularity container (.sif)
THREADS=24
TMPDIR="${TMPDIR:-/tmp}"

# Matrix parameters
MATRIX_MODE="reference-point"  # or "scale-regions"
REFERENCE_POINT="center"       # TSS, TES, or center
UPSTREAM=2000
DOWNSTREAM=2000
REGION_BODY_LENGTH=5000        # for scale-regions mode
SCALE_UPSTREAM=3000
SCALE_DOWNSTREAM=3000

# Heatmap parameters
HEATMAP_HEIGHT=20
HEATMAP_WIDTH=10
M6A_COLORMAP="Greens"
NUCLEOSOME_COLORMAP="Blues"
CPG_COLORMAP="Reds"
COMBINED_COLORMAP="RdYlBu_r"
MARKS_TO_PLOT=("perc_m6a" "perc_nuc" "perc_cpg")


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/06_create_heatmaps_main.sh"

mkdir -p "$TMPDIR"
DEEPTOOLS_PREFIX="singularity exec $SINGULARITY_IMAGE"

export TMPDIR DEEPTOOLS_PREFIX SINGULARITY_IMAGE
export M6A_COLORMAP NUCLEOSOME_COLORMAP CPG_COLORMAP COMBINED_COLORMAP

MARKS_STR=$(IFS=,; echo "${MARKS_TO_PLOT[*]}")

echo "=== Create Heatmaps Pipeline ==="
echo "Samplesheet: $SAMPLESHEET"
echo "Workdir: $WORKDIR"
echo "Regions: $REGIONS_BED"
echo "Mode: $MATRIX_MODE"
echo "Marks: ${MARKS_TO_PLOT[*]}"

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    echo "--- Processing: $sample_name ---"
    bash "$MAIN" "$sample_name" "$SAMPLE_DIR" "$REGIONS_BED" "$MATRIX_MODE" \
        "$REFERENCE_POINT" "$UPSTREAM" "$DOWNSTREAM" "$REGION_BODY_LENGTH" \
        "$SCALE_UPSTREAM" "$SCALE_DOWNSTREAM" "$COMBINED_COLORMAP" \
        "$HEATMAP_HEIGHT" "$HEATMAP_WIDTH" "$THREADS" "$MARKS_STR"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
