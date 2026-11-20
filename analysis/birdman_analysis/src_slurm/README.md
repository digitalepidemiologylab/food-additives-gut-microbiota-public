# BIRDMAn differential abundance analysis (backup scripts)

This directory contains a backup of scripts used to run BIRDMAn for identifying differentially abundant taxa under different dietary and additive exposure models.

Two model types are supported:

- **Unadjusted models**: models without covariate adjustment  
- **Adjusted models**: models including covariates such as age, BMI, diet quality, and other metadata variables

### Scripts overview

- `filter_biom_and_metadata.py`  
  Filters the metadata to match samples present in the BIOM feature table for **unadjusted** models.

- `filter_biom_and_metadata_with_covariate_adjustment.py`  
  Filters the metadata to match samples present in the BIOM feature table and prepares covariates for **adjusted** models.

- `submit_jobs.sh`  
  Submits BIRDMAn jobs for unadjusted models.

- `submit_jobs_adjustedModel.sh`  
  Submits BIRDMAn jobs for adjusted models.

### Execution order

Regardless of the model type, the workflow must be run in the following order:

1. Run the appropriate metadata filtering script:
   - Use `filter_biom_and_metadata.py` for unadjusted models
   - Use `filter_biom_and_metadata_with_covariate_adjustment.py` for adjusted models

2. Submit the corresponding BIRDMAn jobs:
   - Use `submit_jobs.sh` for unadjusted models
   - Use `submit_jobs_adjustedModel.sh` for adjusted models

### Required input files

Before running any script, the following files must be present in the `data/` directory:

- `feature-table.biom`  
  BIOM-format feature table containing microbial counts

- `meta_zscores.csv`  
  Metadata file containing exposure variables and covariates

### Output structure

- Results from **unadjusted models** are written to:
  - `outputs/`


- Results from **adjusted models** are written to:
  - `outputs_adjusted/`


### Execution safeguard

If results for a given variable of interest already exist in the corresponding output directory, execution of the BIRDMAn job is automatically terminated. This prevents accidental overwriting of previously generated results.

