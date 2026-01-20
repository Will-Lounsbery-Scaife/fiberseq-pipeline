#!/bin/bash
# Step 8: k-mer Variant Phasing Wrapper
# Prerequisites: Singularity/Apptainer

# Load singularity module (adjust for your HPC)
ml singularity 2>/dev/null || module load singularity 2>/dev/null || true

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""
WORKDIR=""
REFERENCE_FASTA=""
KMER_PHASING_ROOT=""
TMPDIR="${TMPDIR:-/tmp}"
SNAKEMAKE_CONDA_PREFIX=""
APPTAINER_CACHEDIR=""

ALIGN="false"
FILTERED_CHR="false"
RESUME="false"

# Parental data options (set to empty string "" if no parental data available)
PARENTAL_MODE="false"
MATERNAL_BAM=""
PATERNAL_BAM=""


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/08_kmer_phasing_main.sh"
CONFIG_DIR="${SCRIPT_DIR}/config"
PROFILE_PATH="${KMER_PHASING_ROOT}/profiles/lsf-executor"
PIXI_MANIFEST="${KMER_PHASING_ROOT}/pixi.toml"

# Set default cache directories if not specified
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-${KMER_PHASING_ROOT}/.snakemake_conda_envs}"
APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${KMER_PHASING_ROOT}/.apptainer_cache}"

mkdir -p "$TMPDIR" "$CONFIG_DIR" "$SNAKEMAKE_CONDA_PREFIX" "$APPTAINER_CACHEDIR"
export TMPDIR SNAKEMAKE_CONDA_PREFIX APPTAINER_CACHEDIR

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    # Locate input BAM (from nucleosomes step)
    SOURCE_BAM="${SAMPLE_DIR}/02_nucleosomes/${sample_name}.ft.bam"
    [[ ! -f "$SOURCE_BAM" ]] && { echo "ERROR: Nucleosome BAM not found: $SOURCE_BAM"; continue; }
    [[ ! -f "${SOURCE_BAM}.bai" ]] && { echo "ERROR: BAM index not found: ${SOURCE_BAM}.bai"; continue; }

    # Setup output directory
    OUTPUT_DIR="${SAMPLE_DIR}/08_kmer_phasing"
    mkdir -p "$OUTPUT_DIR"
    FILTERED_BAM="${SAMPLE_DIR}/07_filter_chroms/${sample_name}.filtered.bam"

    if [[ "$FILTERED_CHR" == "true" ]]; then
        INPUT_BAM="$FILTERED_BAM"
        [[ ! -f "$INPUT_BAM" ]] && { echo "ERROR: Filtered BAM not found: $INPUT_BAM (run 07_filter_chroms_wrapper.sh or set FILTERED_CHR=false)"; continue; }
        [[ ! -f "${INPUT_BAM}.bai" ]] && { echo "ERROR: Filtered BAM index not found: ${INPUT_BAM}.bai"; continue; }
    else
        INPUT_BAM="$SOURCE_BAM"
    fi

    # Generate per-sample config YAML
    CONFIG_FILE="${CONFIG_DIR}/${sample_name}_kmer_phasing.yaml"
    echo "sample: \"${sample_name}\"" > "$CONFIG_FILE"
    echo "reference: ${REFERENCE_FASTA}" >> "$CONFIG_FILE"
    echo "hifi_bam: ${INPUT_BAM}" >> "$CONFIG_FILE"
    echo "align: ${ALIGN}" >> "$CONFIG_FILE"
    echo "variant_callers: false" >> "$CONFIG_FILE"  # Disable SV callers (pbsv/sniffles) to avoid contig mismatch issues

    # Add parental data if available
    if [[ "$PARENTAL_MODE" == "true" ]]; then
        echo "maternal: ${MATERNAL_BAM}" >> "$CONFIG_FILE"
        echo "paternal: ${PATERNAL_BAM}" >> "$CONFIG_FILE"
    fi

    echo "--- Processing: $sample_name ---"
    echo "Input BAM: $INPUT_BAM"
    echo "Output dir: $OUTPUT_DIR"
    echo "Config: $CONFIG_FILE"

    bash "$MAIN" "$sample_name" "$OUTPUT_DIR" "$CONFIG_FILE" "$PIXI_MANIFEST" "$PROFILE_PATH" "$RESUME"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
