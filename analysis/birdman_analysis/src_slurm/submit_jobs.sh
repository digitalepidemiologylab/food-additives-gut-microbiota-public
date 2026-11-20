#!/bin/bash

# Check if a variable is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <variable_name>"
  exit 1
fi

# Assign the first argument to a variable
VARIABLE_NAME=$1

# Submit the first job and capture the job ID
job_id1=$(birdman-cli run -i ./data/table-filt.biom -o ./outputs/$VARIABLE_NAME -m ./data/metadata-filt.tsv -f "$VARIABLE_NAME" -e YOUR_EMAIL | grep -oP 'Submitted batch job \K\d+')

if [ -z "$job_id1" ]; then
  echo "Failed to capture job ID for the first job."
  exit 1
fi

echo "First job submitted with ID: $job_id1"

# Submit the second job with dependency on the first job
second_job=$(sbatch --dependency=afterok:$job_id1 --wrap="birdman-cli summarize -i outputs/$VARIABLE_NAME && { rm -r outputs/$VARIABLE_NAME/inferences/; rm -r outputs/$VARIABLE_NAME/tmp/; rm -r outputs/$VARIABLE_NAME/logs/; rm -r outputs/$VARIABLE_NAME/slurm_out/; }")

if [ $? -ne 0 ]; then
  echo "Failed to submit the second job."
  exit 1
fi

echo "Second job submitted: $second_job"
