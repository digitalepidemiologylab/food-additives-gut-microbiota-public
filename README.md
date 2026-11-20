# Food additive exposure and gut microbiota diversity

This repository contains the analysis code for the manuscript:

**Food additive exposure associated with reduction in gut diversity**
Rohan Singh, Daniel McDonald, Rob Knight, Marcel Salathé

---

## Analysis pipeline

Run the notebooks in the following order to generate additive exposure variables, mixtures, and composite metrics used across all downstream analyses:

1. `processing/additive_mapping_and_preprocessing.ipynb`
   Generates participant-level food additive exposure variables from MyFoodRepo barcode logs using OpenFoodFacts ingredient data.

2. `analysis/mixtures_analysis/mixtures_generation.ipynb`
   Identifies additive co-consumption patterns using non-negative matrix factorization and outputs additive mixture scores per participant.

3. `processing/zscore_and_load_difference_prep.ipynb`
   Z-standardizes dietary and additive variables and computes composite contrast metrics (e.g., vegetables minus additives, HEI minus preservatives) used in downstream analyses.

All subsequent analyses and plotting takes outputs from these steps.

---

## Setup

This repository contains shared Python utilities in `src/additives_microbiota/` that are imported by notebooks and scripts across `analysis/` and `processing/`.
To make these imports work reliably, install the repository in editable mode.

### 1) Create environment

Use your preferred environment manager (conda, venv, uv), then install:

```bash
pip install -e .
```

---

## Running analyses

* Analysis notebooks are located in `analysis/` and data processing notebooks in `processing/`.
* Figures and tables are written to `analysis/figures/` and `analysis/outputs/` (some notebooks may write to subdirectories).
* Large intermediate outputs (e.g., BIRDMAn HPC results) are not tracked in git.

Detailed descriptions of individual analyses and figure-level mappings are provided in `analysis/README.md`.

---

## License

The code in this repository is released under the MIT License (see `LICENSE`).
