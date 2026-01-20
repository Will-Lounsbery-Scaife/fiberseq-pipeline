#!/bin/bash
# =============================================================================
# PRE-PULL SINGULARITY CONTAINERS
# =============================================================================
# Run this script from a LOGIN NODE (not compute node) before running the pipeline.
# This downloads all required containers to a shared cache directory.
#
# Usage:
#   bash pull_containers.sh
# =============================================================================

set -e

# CONFIGURE: If your HPC requires a proxy for internet access, set it here
# Example for institutions with proxy servers:
# export http_proxy=http://your-proxy.institution.edu:3128
# export https_proxy=http://your-proxy.institution.edu:3128

# CONFIGURE: Set your Singularity cache directory
# Recommended: shared project storage accessible from compute nodes
CACHE_DIR="${NXF_SINGULARITY_CACHEDIR:-${HOME}/.singularity/cache}"
mkdir -p "${CACHE_DIR}"

# Clear any conflicting singularity env vars and set fresh
unset SINGULARITY_PULLFOLDER
unset SINGULARITY_LOCALCACHEDIR
export SINGULARITY_CACHEDIR="${CACHE_DIR}"
export SINGULARITY_TMPDIR="${CACHE_DIR}"

# Load singularity module (adjust for your HPC)
ml singularity 2>/dev/null || module load singularity 2>/dev/null || true

echo "=============================================="
echo "Pulling Singularity containers to: ${CACHE_DIR}"
echo "=============================================="

# List of containers used by the pipeline (Galaxy depot URLs - no Docker required)
CONTAINERS=(
    "https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0"
    "https://depot.galaxyproject.org/singularity/fibertools-rs:0.7.1--h3b373d1_0"
    "https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1"
    "https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0"
)

for container in "${CONTAINERS[@]}"; do
    echo ""
    echo "Pulling: ${container}"
    echo "----------------------------------------------"

    # Extract image name for the .img file (Nextflow naming convention)
    # Galaxy depot URL: https://depot.galaxyproject.org/singularity/tool:version
    img_name=$(echo "$container" | sed 's|https://||' | tr '/:' '-')
    img_file="${CACHE_DIR}/${img_name}.img"

    if [[ -f "${img_file}" ]]; then
        echo "Already cached: ${img_file}"
    else
        singularity pull "${img_file}" "${container}"
        echo "Downloaded: ${img_file}"
    fi
done

echo ""
echo "=============================================="
echo "All containers pulled successfully!"
echo "Cache directory: ${CACHE_DIR}"
echo ""
echo "Contents:"
ls -lh "${CACHE_DIR}"/*.img 2>/dev/null || echo "No .img files found"
ls -lh "${CACHE_DIR}"/*.sif 2>/dev/null || echo "No .sif files found"
echo "=============================================="
