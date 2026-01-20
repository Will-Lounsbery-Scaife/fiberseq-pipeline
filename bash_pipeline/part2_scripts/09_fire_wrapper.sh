#!/bin/bash
# Step 9: FIRE (Fiber-seq Inferred Regulatory Elements) Wrapper
# Prerequisites: Singularity/Apptainer
# Usage: Run from an interactive session or submit via bsub
# Note: FIRE requires whole genome data for FDR calibration - use KEEP_CHROMOSOMES to limit analysis scope

# Load singularity module (adjust for your HPC)
ml singularity 2>/dev/null || module load singularity 2>/dev/null || true

### CONFIGURE ###
SAMPLESHEET=""
WORKDIR=""
REFERENCE_FASTA=""
FIRE_ROOT=""
TMPDIR="${TMPDIR:-/tmp}"
SNAKEMAKE_CONDA_PREFIX=""
KEEP_CHROMOSOMES=""
REF_NAME="hg38"
INPUT_SOURCE="kmer_phasing"

# Peak calling thresholds (leave blank for FIRE defaults)
MAX_PEAK_FDR=""             # FDR threshold for FIRE peaks (default: 0.05); increase for more peaks
MIN_FIRE_FDR=""             # FDR for individual FIRE elements (default: 0.10); increase for more elements
MIN_COVERAGE=""             # Minimum coverage for peak calling (default: 4); decrease for low-coverage regions
COVERAGE_WITHIN_N_SD=""     # Filter peaks beyond N SDs from mean coverage (default: 5)
MIN_MSP=""                  # Minimum MSPs per Fiber-seq read (default: 10)
MIN_AVE_MSP_SIZE=""         # Minimum average MSP size (default: 10)


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/09_fire_main.sh"
CONFIG_DIR="${SCRIPT_DIR}/config"
MANIFEST_DIR="${SCRIPT_DIR}/config/manifests"
PROFILE_PATH="${FIRE_ROOT}/profiles/lsf-executor"
PIXI_MANIFEST="${FIRE_ROOT}/pixi.toml"

# Set default cache directory if not specified
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-${FIRE_ROOT}/.snakemake_conda_envs}"

mkdir -p "$TMPDIR" "$CONFIG_DIR" "$MANIFEST_DIR" "$SNAKEMAKE_CONDA_PREFIX"
export TMPDIR SNAKEMAKE_CONDA_PREFIX


tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    SAMPLE_DIR="${WORKDIR}/${sample_name}"
    [[ ! -d "$SAMPLE_DIR" ]] && { echo "ERROR: Sample dir not found: $SAMPLE_DIR"; continue; }

    # Locate input BAM based on input source
    if [[ "$INPUT_SOURCE" == "kmer_phasing" ]]; then
        # k-mer phasing outputs phased BAM in hiphase subdirectory
        INPUT_BAM="${SAMPLE_DIR}/08_kmer_phasing/results/${sample_name}/hiphase/${sample_name}.bam"
    else
        INPUT_BAM="${SAMPLE_DIR}/02_nucleosomes/${sample_name}.ft.bam"
    fi

    [[ ! -f "$INPUT_BAM" ]] && { echo "ERROR: Input BAM not found: $INPUT_BAM"; continue; }
    [[ ! -f "${INPUT_BAM}.bai" ]] && { echo "ERROR: BAM index not found: ${INPUT_BAM}.bai"; continue; }

    # Setup output directory
    OUTPUT_DIR="${SAMPLE_DIR}/09_FIRE"
    mkdir -p "$OUTPUT_DIR"

    # Generate per-sample manifest file (FIRE requires sample\tbam format)
    MANIFEST_FILE="${MANIFEST_DIR}/${sample_name}_fire_manifest.tbl"
    {
        echo -e "sample\tbam"
        echo -e "${sample_name}\t${INPUT_BAM}"
    } > "$MANIFEST_FILE"

    # Generate per-sample config YAML
    CONFIG_FILE="${CONFIG_DIR}/${sample_name}_fire_config.yaml"
    {
        echo "ref: ${REFERENCE_FASTA}"
        echo "ref_name: ${REF_NAME}"
        echo "manifest: ${MANIFEST_FILE}"
        # Optional parameters (only added if non-empty)
        [[ -n "$KEEP_CHROMOSOMES" ]] && echo "keep_chromosomes: \"${KEEP_CHROMOSOMES}\""
        [[ -n "$MAX_PEAK_FDR" ]] && echo "max_peak_fdr: ${MAX_PEAK_FDR}"
        [[ -n "$MIN_FIRE_FDR" ]] && echo "min_fire_fdr: ${MIN_FIRE_FDR}"
        [[ -n "$MIN_COVERAGE" ]] && echo "min_coverage: ${MIN_COVERAGE}"
        [[ -n "$COVERAGE_WITHIN_N_SD" ]] && echo "coverage_within_n_sd: ${COVERAGE_WITHIN_N_SD}"
        [[ -n "$MIN_MSP" ]] && echo "min_msp: ${MIN_MSP}"
        [[ -n "$MIN_AVE_MSP_SIZE" ]] && echo "min_ave_msp_size: ${MIN_AVE_MSP_SIZE}"
    } > "$CONFIG_FILE"


    bash "$MAIN" "$sample_name" "$OUTPUT_DIR" "$CONFIG_FILE" "$PIXI_MANIFEST" "$PROFILE_PATH"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
