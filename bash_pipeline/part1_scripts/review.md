# Pipeline Code Review

**Reviewed by:** Automated Code Review  
**Date:** November 25, 2025  
**Scope:** `pipeline/scripts/` (excluding `scripts_v1`)

---

## Summary

| Severity | Count |
|----------|-------|
| 🔴 Critical | 0 (2 fixed) |
| 🟠 Medium | 4 |
| 🟡 Low | 4 |

---

## 🔴 Critical Issues (FIXED)

### 1. ✅ FIXED: Incomplete "Generating Summary Statistics" Section in `04_create_pileups_main.sh`

**Location:** Lines 183-199

**Description:** The "Generating summary statistics" section was incomplete.

**Fix Applied:** Added comprehensive statistics generation including:
- File size and total line count
- Chromosome coverage breakdown (positions per chromosome)
- Column statistics (mean, min, max) for coverage and mark columns
- Statistics are written to `${SAMPLE_NAME}.pileup_stats.txt`

---

### 2. ✅ FIXED: `THREADS` Variable Not Used in `03_fiberseq_qc_main.sh`

**Location:** Line 6 (parameter), entire script

**Description:** The `THREADS` parameter was passed but never used.

**Fix Applied:** Removed the unused `THREADS` parameter from both:
- `03_fiberseq_qc_main.sh` (removed parameter parsing)
- `03_fiberseq_qc_wrapper.sh` (removed THREADS variable and argument passing)

**Note:** fiberseq-qc does not support multithreading, so this parameter was unnecessary.

---

## 🟠 Medium Issues

### 3. Missing Validation for `BEDGRAPHTOBIGWIG` in `05_pileup_to_bigwig_wrapper.sh`

**Location:** Lines 7-8

**Description:** The `BEDGRAPHTOBIGWIG` variable is assigned from `command -v` with `|| true`, which means it will be empty if the tool isn't found. There's no validation before passing it to the main script.

```bash
ml ucsc-utils/2023-10-17 # provides bedGraphToBigWig
BEDGRAPHTOBIGWIG="$(command -v bedGraphToBigWig 2>/dev/null || true)"
# No check if BEDGRAPHTOBIGWIG is empty!
```

**Impact:** The main script will fail with a confusing error message if `bedGraphToBigWig` isn't in PATH.

**Recommendation:** Add validation:
```bash
if [ -z "$BEDGRAPHTOBIGWIG" ]; then
    echo "ERROR: bedGraphToBigWig not found in PATH"
    echo "Please load the ucsc-utils module"
    exit 1
fi
```

---

### 4. Reference Index Extension Handling May Fail for Non-Standard Extensions

**Location:** `01_align_and_qc_main.sh`, Lines 78-82

**Description:** The logic to determine the reference index filename only handles `.fasta` and `.fa` extensions. References with other common extensions (`.fna`, `.fas`, `.fa.gz`, etc.) will result in incorrect index paths.

```bash
REFERENCE_INDEX="${REFERENCE_FASTA%.fasta}.mmi"
if [[ "$REFERENCE_INDEX" == "$REFERENCE_FASTA" ]]; then
    # .fasta wasn't found, try .fa
    REFERENCE_INDEX="${REFERENCE_FASTA%.fa}.mmi"
fi
# No handling for .fna, .fas, or compressed files
```

**Impact:** Pipeline will fail to find or create index for references with non-standard extensions.

**Recommendation:** Use a more robust approach:
```bash
# Remove any .gz first, then handle common extensions
BASE="${REFERENCE_FASTA%.gz}"
for ext in .fasta .fa .fna .fas; do
    if [[ "$BASE" == *"$ext" ]]; then
        REFERENCE_INDEX="${BASE%$ext}.mmi"
        break
    fi
done
```

---

### 5. Potential Race Condition with Shared Reference Index

**Location:** `01_align_and_qc_main.sh`, Lines 84-105

**Description:** When running multiple samples in parallel, multiple jobs could attempt to create the reference index simultaneously if it doesn't exist, potentially causing corruption or race conditions.

```bash
if [ ! -f "$REFERENCE_INDEX" ]; then
    echo "Creating pbmm2 index..."
    $PBMM2 index ...
fi
```

**Impact:** Data corruption or failed indexing if samples are processed in parallel.

**Recommendation:** Use file locking or create the index in a separate preparation step before running samples:
```bash
LOCK_FILE="${REFERENCE_INDEX}.lock"
if [ ! -f "$REFERENCE_INDEX" ]; then
    (
        flock -x 200
        # Double-check after acquiring lock
        if [ ! -f "$REFERENCE_INDEX" ]; then
            $PBMM2 index ...
        fi
    ) 200>"$LOCK_FILE"
fi
```

---

### 6. Module Load Commands May Fail Silently

**Location:** 
- `05_pileup_to_bigwig_wrapper.sh`, Line 7
- `06_create_heatmaps_wrapper.sh`, Line 8

**Description:** The `module load` / `ml` commands don't check for success. If the module doesn't exist or can't be loaded, the script continues and will fail later with confusing errors.

```bash
ml ucsc-utils/2023-10-17 # provides bedGraphToBigWig
# or
module load macs/2.1.0 singularity/3.6.4
```

**Impact:** Confusing downstream errors when required tools aren't available.

**Recommendation:** Check module load success:
```bash
if ! ml ucsc-utils/2023-10-17 2>/dev/null; then
    echo "ERROR: Failed to load ucsc-utils module"
    exit 1
fi
```

---

## 🟡 Low Issues

### 7. Hardcoded Expected Output Files in `03_fiberseq_qc_main.sh`

**Location:** Lines 146-149

**Description:** The script checks for `qc_summary.txt` and `plots` as expected outputs, but these may not match actual fiberseq-qc outputs.

```bash
EXPECTED_FILES=(
    "qc_summary.txt"
    "plots"
)
```

**Impact:** False warnings about missing files; doesn't verify actual expected outputs.

**Recommendation:** Verify actual fiberseq-qc outputs and update the expected file list accordingly.

---

### 8. While Loop in Subshell (All Wrapper Scripts)

**Location:** All `*_wrapper.sh` files

**Description:** Using `tail | while read` creates a subshell, meaning any variables set inside the loop don't persist outside.

```bash
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r ...; do
    # Variables set here won't persist
done
```

**Impact:** Can't track overall success/failure counts or aggregate data across samples.

**Recommendation:** For this use case it's acceptable, but if you need to track data:
```bash
while IFS=$'\t' read -r ...; do
    ...
done < <(tail -n +2 "$SAMPLESHEET")
```

---

### 9. Inconsistent Quoting in File Paths

**Location:** Multiple scripts

**Description:** Some file paths are not properly quoted, which could cause issues with paths containing spaces or special characters.

Examples:
- `05_pileup_to_bigwig_main.sh`, Line 50: `PILEUP_FILES=($(find "$PILEUP_DIR" ...))`
- `05_pileup_to_bigwig_main.sh`, Line 147: `$BEDGRAPHTOBIGWIG "$BEDGRAPH_FILE" ...`

**Impact:** Scripts will fail if any paths contain spaces.

**Recommendation:** Use consistent quoting and consider using `mapfile`/`readarray` for array population from `find`.

---

### 10. Directory Numbering Confusion

**Location:** All scripts

**Description:** Output directory numbers don't align with step numbers, which could cause confusion:

| Step | Script | Output Directory |
|------|--------|------------------|
| 1 | 01_align_and_qc | `01_aligned/`, `02_sequencing_qc/` |
| 2 | 02_add_nucleosomes | `03_nucleosomes/` |
| 3 | 03_fiberseq_qc | `04_fiberseq_qc/` |
| 4 | 04_create_pileups | `05_pileups/` |
| 5 | 05_pileup_to_bigwig | `06_bigwigs/` |
| 6 | 06_create_heatmaps | `07_heatmaps/` |

**Impact:** Cognitive overhead when navigating outputs; potential for mistakes when manually inspecting results.

**Recommendation:** Either:
- Match directory numbers to step numbers (e.g., step 2 → `02_nucleosomes`)
- Or document this numbering scheme clearly in a README

---

## Additional Notes

### Positive Observations

1. **Good error handling**: Most scripts include comprehensive validation checks at the start
2. **Helpful logging**: Timestamps and progress messages aid debugging
3. **Modular design**: Separation of wrapper and main scripts allows for flexibility
4. **Clean exit codes**: Proper use of `exit 0` and `exit 1`

### Suggested Improvements

1. **Add a pipeline README**: Document the overall workflow, dependencies, and expected inputs/outputs
2. **Create a master wrapper**: A single script to run all steps sequentially with proper error handling
3. **Add `set -euo pipefail`**: At the top of scripts to catch errors early (with careful handling of expected failures)
4. **Version tracking**: Log tool versions used for reproducibility

---

*End of Review*

