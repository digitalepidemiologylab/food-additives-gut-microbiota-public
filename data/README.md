# Data directory

This directory contains a subset of data files required to reproduce specific analyses and visualizations in this repository. Only non-sensitive and derived files are included here.

## Included files

### Additive annotation files
- additive_name_mapping_v3.csv  
  Mapping between E codes, additive names, and additive categories.  
  Used for network visualization and food additive annotation.

- additive_name_mapping_v3_merged.csv  
  Additive annotation table with combined additive classes for downstream analyses.

### Additive co occurrence network data
- additive_analysis_files/json/  
  JSON files containing additive co occurrence counts within specific food groups.  
  These files are used by the interactive additive network visualization located in  
  analysis/food_group_associations/network/.

### Aggregated additive summaries
- additive_analysis_files/additive_prevalence_across_users.csv  
  Prevalence of each food additive across study participants (Supplementary Data 1).

- additive_analysis_files/summed_additives_across_food_categories.csv  
  Number of unique barcoded products containing each additive (E code), per food group.

- additive_analysis_files/summed_additives_types_across_food_categories.csv  
  The same counts summed by additive category, per food group.

## Not included in this repository

The following data are intentionally excluded:

### Food Intake data : 
- *.pkl files  
  Large intermediate objects generated during data processing and analysis.

### Metadata
- metadata_files/  
  Metadata containing demographic and nutritional variables.

These files are excluded to reduce repository size and because they contain sensitive participant information.

## Access to raw and sensitive data

Sequencing data from this study have been deposited in:
- European Nucleotide Archive (ENA) under project accession PRJEB85942
- Qiita under study ID 15880

Metadata containing clinical, demographic, and nutritional variables cannot be deposited publicly due to participant privacy and ethical restrictions. Access to this metadata can be requested by contacting the corresponding author, subject to institutional ethical compliance.

## Notes

Users attempting to reproduce the full analysis pipeline will need to obtain the restricted metadata separately from the corresponding author (Marcel Salathé).

## Supplementary materials

Supplementary tables and files associated with this study are available in the analysis folder: [`analysis/outputs`](../analysis/outputs)

