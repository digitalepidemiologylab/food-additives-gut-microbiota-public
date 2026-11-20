# MaAsLin2 setup for PICRUSt associations

# install.packages("Maaslin2") using version:1.14.1
library(Maaslin2)

# Paths
meta_path <- "../../data/metadata_files/meta_zscores.csv"

picrust_data_type <- "pathway"  # pathway | ko | ec

picrust_pathway_path <- "~/WORK/qiime2_2024/qiime_processed_files/exported_rep_seqs/pathways_out/path_abun_unstrat_descrip.tsv"
picrust_ko_path      <- "~/WORK/qiime2_2024/qiime_processed_files/exported_rep_seqs/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv"
picrust_ec_path      <- "~/WORK/qiime2_2024/qiime_processed_files/exported_rep_seqs/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv"

out_dir_base <- file.path("results", paste0("maaslin2_", picrust_data_type))

# Read metadata
meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)

# MaAsLin2 expects sample IDs as rownames
stopifnot("sample-id" %in% colnames(meta))
rownames(meta) <- meta[["sample-id"]]

# Ensure categorical covariates are factors
meta$gender <- as.factor(meta$gender)
meta$antibiotics_treatment <- as.factor(meta$antibiotics_treatment)
meta$smoking <- factor(meta$smoking, levels = c("non-smoker", setdiff(unique(meta$smoking), "non-smoker")))
meta$batch <- as.factor(meta$batch)
meta$language <- as.factor(meta$language)
meta$urbanity <- as.factor(meta$urbanity)

read_picrust <- function(path, type) {
  df <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (type == "pathway") {
    desc_map <- setNames(df$description, df$pathway)
    rownames(df) <- df$pathway
    df <- df[, setdiff(colnames(df), c("pathway", "description")), drop = FALSE]
  } else {
    desc_map <- setNames(df$description, df[["function"]])
    rownames(df) <- df[["function"]]
    df <- df[, setdiff(colnames(df), c("function", "description")), drop = FALSE]
  }
  
  sample_by_feat <- as.data.frame(t(df), check.names = FALSE)
  list(data = sample_by_feat, desc = desc_map)
}


pic <- switch(
  picrust_data_type,
  "pathway" = read_picrust(picrust_pathway_path, "pathway"),
  "ko"      = read_picrust(picrust_ko_path, "ko"),
  "ec"      = read_picrust(picrust_ec_path, "ec"),
  stop("Invalid picrust_data_type")
)

features <- pic$data
desc_map <- pic$desc

desc_map_maaslin <- desc_map
names(desc_map_maaslin) <- make.names(names(desc_map_maaslin))

# Align samples between metadata and features
common_ids <- intersect(rownames(meta), rownames(features))
meta <- meta[common_ids, , drop = FALSE]
features <- features[common_ids, , drop = FALSE]

# Optional: prevalence + variance filtering similar to your Python
min_prev <- 0.05
min_var  <- 1e-8

# prevalence on nonzero
keep_prev <- colMeans(features != 0, na.rm = TRUE) >= min_prev
features_f <- features[, keep_prev, drop = FALSE]

# variance filter
keep_var <- apply(features_f, 2, var, na.rm = TRUE) > min_var
features_f <- features_f[, keep_var, drop = FALSE]

# Define common covariates 
base_covariates <- c(
  "urbanity", "language",
  "gender", "age", "bmi",
  "antibiotics_treatment",
  "smoking",
  "defecate_quantity_per_day",
  "batch"
  # "vegetables_fruits"
)

additive_vars <- c( "thickeners_stabilizers_emulsifiers_gelling_agents",
                    "sweeteners_glazing_agents",
                    "preservatives_antioxidants",
                    "ph_regulators_anti_caking_leavening_agents",
                    "flavor_enhancers",
                    "colorants")
mixtures <- c("Mixture1","Mixture2","Mixture3")
exposures <- c(additive_vars, mixtures)

# Run MaAsLin2 per exposure, so each exposure is the "main effect"
# and all other covariates are adjusted for.
run_one_exposure <- function(exposure, features_mat, meta_df, out_dir) {
  # Ensure output folder unique per exposure
  out_path <- file.path(out_dir, exposure)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  fixed_effects <- c(base_covariates, exposure)
  fixed_effects <- fixed_effects[fixed_effects %in% colnames(meta_df)]

  # Drop rows with missing in any model variables
  keep_rows <- complete.cases(meta_df[, fixed_effects, drop = FALSE])
  meta_use <- meta_df[keep_rows, , drop = FALSE]
  feat_use <- features_mat[rownames(meta_use), , drop = FALSE]

  fit <- Maaslin2(
    input_data = feat_use,
    input_metadata = meta_use,
    output = out_path,
    fixed_effects = fixed_effects,
    random_effects = c(),
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = 1.0,
    standardize = FALSE
  )

  invisible(fit)
}

for (exposure in exposures) {
  if (!(exposure %in% colnames(meta))) next
  run_one_exposure(exposure, features_f, meta, out_dir_base)
}

# Combine all per-exposure results and do a global BH across all feature-exposure pairs 
collect_results <- function(out_dir) {
  exp_dirs <- list.dirs(out_dir, full.names = TRUE, recursive = FALSE)
  all <- list()

  for (d in exp_dirs) {
    f <- file.path(d, "all_results.tsv")
    if (!file.exists(f)) next
    res <- read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    exposure <- basename(d)
    res$food <- exposure

    all[[length(all) + 1]] <- res
  }

  if (length(all) == 0) return(NULL)
  out <- do.call(rbind, all)
  out
}

res_all <- collect_results(out_dir_base)

if (!is.null(res_all)) {
  # Keep only the exposure term for each run
  res_all_main <- res_all[res_all$metadata == res_all$food, , drop = FALSE]

  # Global BH across all main terms
  res_all_main$qval_bh_global <- p.adjust(res_all_main$pval, method = "BH")
  res_all_main$signif_fdr05_global <- res_all_main$qval_bh_global < 0.05

  if (!is.null(desc_map)) {
    res_all_main$description <- unname(desc_map[res_all_main$feature])
    na_idx <- is.na(res_all_main$description)
    res_all_main$description[na_idx] <-
      unname(desc_map_maaslin[make.names(res_all_main$feature[na_idx])])
  }
  

  write.csv(res_all_main,
            file = file.path(out_dir_base, paste0("maaslin2_associations_", picrust_data_type, "_globalBH.csv")),
            row.names = FALSE)
  
  write.csv(hits,
            file = file.path(out_dir_base, paste0("maaslin2_hits_", picrust_data_type, "_globalBH.csv")),
            row.names = FALSE)
  
}
