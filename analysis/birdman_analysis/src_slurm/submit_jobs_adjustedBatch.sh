#!/bin/bash

set -euo pipefail

if [ -z "${1:-}" ]; then
  echo "Usage: $0 <variable_name>"
  exit 1
fi

VARIABLE_NAME="$1"

# Covariates to adjust for (no spaces, Patsy RHS style)
BASE_COVARS="batch"

# Build formula RHS with no spaces
FORMULA="${VARIABLE_NAME}+${BASE_COVARS}"

OUT_DIR="./outputs_adjustedBatch/$VARIABLE_NAME"

echo "Running BIRDMAn for variable: $VARIABLE_NAME"
echo "Formula RHS: $FORMULA"
echo "Output dir: $OUT_DIR"

job_id1=$(
  birdman-cli run \
    -i ./data/table-filt.biom \
    -o "$OUT_DIR" \
    -m ./data/metadata-filt.tsv \
    -f "$FORMULA" \
    -e YOUR_EMAIL | grep -oP 'Submitted batch job \K\d+'
)

if [ -z "$job_id1" ]; then
  echo "Failed to capture job ID for the first job."
  exit 1
fi

echo "First job submitted with ID: $job_id1"

second_job=$(
  sbatch --dependency=afterok:$job_id1 \
    --wrap="birdman-cli summarize -i '$OUT_DIR' && { rm -r '$OUT_DIR/inferences/' '$OUT_DIR/tmp/' '$OUT_DIR/logs/' '$OUT_DIR/slurm_out/' 2>/dev/null || true; }"
)

if [ $? -ne 0 ]; then
  echo "Failed to submit the second job."
  exit 1
fi

echo "Second job submitted: $second_job"

