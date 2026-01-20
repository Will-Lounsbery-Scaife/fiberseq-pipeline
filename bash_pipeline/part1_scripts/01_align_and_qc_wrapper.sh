#!/bin/bash
# Step 1 Wrapper: Alignment and Sequencing QC
# Usage: bash 01_align_and_qc_wrapper.sh
#
# Example LSF submission (adjust resources as needed):
# bsub -P ${LSF_PROJECT} -q ${LSF_QUEUE} -n 12 -W 36:00 -R rusage[mem=64000] -Is /bin/bash

### CONFIGURE: Set these paths for your environment ###
SAMPLESHEET=""        # CONFIGURE: Path to your samplesheet TSV file
REFERENCE_FASTA=""    # CONFIGURE: Path to your reference genome FASTA
WORKDIR=""            # CONFIGURE: Base output directory for results
SMRT_ROOT=""          # CONFIGURE: Path to SMRT Tools installation (optional, for QC)
THREADS=12
PRESET="HIFI"         # Options: SUBREAD, CCS, HIFI, ISOSEQ, UNROLLED
TMPDIR="${TMPDIR:-/tmp}"  # Uses system TMPDIR or /tmp as fallback

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN="${SCRIPT_DIR}/01_align_and_qc_main.sh"

mkdir -p "$WORKDIR" "$TMPDIR"
export TMPDIR SMRT_ROOT

tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r sample_name raw_bam_path tissue sequencer; do
    [[ -z "$sample_name" ]] && continue

    echo "Processing: $sample_name"
    [[ ! -f "$raw_bam_path" ]] && { echo "ERROR: BAM not found: $raw_bam_path"; continue; }

    bash "$MAIN" "$sample_name" "${WORKDIR}/${sample_name}" "$raw_bam_path" "$REFERENCE_FASTA" "$THREADS" "$PRESET"
    [[ $? -ne 0 ]] && echo "ERROR: Failed for $sample_name" || echo "Done: $sample_name"
done

echo "=== Pipeline complete ==="
exit 0
