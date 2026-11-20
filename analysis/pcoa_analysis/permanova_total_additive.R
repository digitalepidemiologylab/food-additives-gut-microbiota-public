library(dplyr)
library(vegan)

# =========================
# PERMANOVA script (no PCoA)
# additive_intake quartiles
# - unadjusted adonis2
# - betadisper + permutest
# - adjusted adonis2 with covariates
# =========================

# -------------------------
# Paths
# -------------------------
meta_path <- "../../data/metadata_files/meta_zscores.csv"
distance_matrix_path <- "~/WORK/dietary-microbiota-paper/qiime/unweighted_unifrac_distance_matrix/distance-matrix.tsv"

add_col <- "ph_regulators_anti_caking_leavening_agents"

out_dir_res <- "./results"
if (!dir.exists(out_dir_res)) dir.create(out_dir_res, recursive = TRUE)

out_txt <- file.path(out_dir_res, paste0("permanova_", add_col, "_unweighted_unifrac.txt"))
out_csv_unadj <- file.path(out_dir_res, paste0("permanova_unadjusted_", add_col, ".csv"))
out_csv_adj <- file.path(out_dir_res, paste0("permanova_adjusted_", add_col, ".csv"))
out_csv_betadisper_overall <- file.path(out_dir_res, paste0("betadisper_overall_", add_col, ".csv"))
out_csv_betadisper_pairwise <- file.path(out_dir_res, paste0("betadisper_pairwise_", add_col, ".csv"))

# -------------------------
# Load metadata
# -------------------------
meta <- read.csv(meta_path, check.names = FALSE)
meta <- meta[, names(meta) != "", drop = FALSE]
meta$sample.id <- as.character(meta[["sample-id"]])

# -------------------------
# Load distance matrix
# -------------------------
beta_dist_tbl <- read.table(distance_matrix_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
beta_dist <- as.dist(beta_dist_tbl)
dist_ids <- rownames(as.matrix(beta_dist))

# -------------------------
# Align samples (strict)
# -------------------------
shared_ids <- intersect(meta$sample.id, dist_ids)
meta_sub <- meta[meta$sample.id %in% shared_ids, , drop = FALSE]
meta_sub <- meta_sub[match(shared_ids, meta_sub$sample.id), , drop = FALSE]

dist_mat <- as.matrix(beta_dist)
dist_sub <- as.dist(dist_mat[shared_ids, shared_ids])

# -------------------------
# Create quartiles for additive feature
# -------------------------
if (!(add_col %in% colnames(meta_sub))) stop("Missing additive feature in metadata.")

qs <- quantile(meta_sub[[add_col]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# If ties cause non-unique breaks, jitter slightly
if (length(unique(qs)) < length(qs)) {
  vals2 <- meta_sub[[add_col]] + rnorm(nrow(meta_sub), mean = 0, sd = 1e-10)
  meta_sub[[add_col]] <- vals2
  qs <- quantile(meta_sub[[add_col]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
}

meta_sub$additive_quartile <- cut(
  meta_sub[[add_col]],
  breaks = qs,
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)
meta_sub$additive_quartile <- factor(meta_sub$additive_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))

# -------------------------
# Helper: coerce factors
# -------------------------
to_factor_if_present <- function(df, col) {
  if (col %in% colnames(df)) df[[col]] <- factor(df[[col]])
  df
}

meta_sub <- to_factor_if_present(meta_sub, "urbanity")
meta_sub <- to_factor_if_present(meta_sub, "language")
meta_sub <- to_factor_if_present(meta_sub, "gender")
meta_sub <- to_factor_if_present(meta_sub, "antibiotics_treatment")
meta_sub <- to_factor_if_present(meta_sub, "smoking")
meta_sub <- to_factor_if_present(meta_sub, "batch")

if ("smoking" %in% colnames(meta_sub)) {
  meta_sub$smoking <- relevel(factor(meta_sub$smoking), ref = "non-smoker")
}

# -------------------------
# Run analyses and capture output
# -------------------------
sink(out_txt)
cat("PERMANOVA (unweighted UniFrac) for additive quartiles\n")
cat("Samples in distance matrix:", length(dist_ids), "\n")
cat("Samples in metadata:", nrow(meta), "\n")
cat("Samples used (intersection):", length(shared_ids), "\n\n")

cat("Unadjusted PERMANOVA: dist ~ additive_quartile\n")
perm_unadj <- adonis2(
  dist_sub ~ additive_quartile,
  data = meta_sub,
  permutations = 999
)
print(perm_unadj)
cat("\n")

cat("Betadisper overall + permutest\n")
bd <- betadisper(dist_sub, meta_sub$additive_quartile)
bd_overall <- permutest(bd, permutations = 999)
print(bd_overall)
cat("\n")

cat("Betadisper pairwise + permutest\n")
bd_pairwise <- permutest(bd, pairwise = TRUE, permutations = 999)
print(bd_pairwise)
cat("\n")

needed_cols <- c(
  #"additive_quartile",
  add_col,
  "urbanity", "language",
  "gender", "age", "bmi",
  "antibiotics_treatment",
  "smoking",
  "HEI",
  "alcohol_eaten",
  "energy_kcal_eaten",
  "defecate_quantity_per_day",
  "batch"
)
missing_cols <- setdiff(needed_cols, colnames(meta_sub))

if (length(missing_cols) > 0) {
  cat("WARNING: Missing columns for adjusted model:\n")
  print(missing_cols)
  perm_adj <- NULL
} else {
  cat("Adjusted PERMANOVA (marginal effects): dist ~ additive intake + covariates\n")
  form <- as.formula(
    paste(
      "dist_sub ~",
      paste(
        c(
          add_col,
          "urbanity", "language",
          "gender", "age", "bmi",
          "antibiotics_treatment",
          "smoking",
          "HEI",
          "alcohol_eaten",
          "energy_kcal_eaten",
          "defecate_quantity_per_day",
          "batch"
        ),
        collapse = " + "
      )
    )
  )
  
  perm_adj <- adonis2(
    form,
    data = meta_sub,
    permutations = 999,
    by = "margin"
  )
  print(perm_adj)
}
sink()

# -------------------------
# Save tables (CSV)
# -------------------------
unadj_df <- as.data.frame(perm_unadj)
unadj_df$term <- rownames(unadj_df)
rownames(unadj_df) <- NULL
write.csv(unadj_df, out_csv_unadj, row.names = FALSE)

bd_overall_df <- as.data.frame(bd_overall$tab)
bd_overall_df$term <- rownames(bd_overall_df)
rownames(bd_overall_df) <- NULL
write.csv(bd_overall_df, out_csv_betadisper_overall, row.names = FALSE)

# pair_mat <- bd_pairwise$pairwise
# if (!is.null(pair_mat)) {
#   pw <- as.data.frame(as.table(pair_mat))
#   colnames(pw) <- c("group1", "group2", "p_value")
#   pw <- pw[pw$group1 != pw$group2, , drop = FALSE]
#   write.csv(pw, out_csv_betadisper_pairwise, row.names = FALSE)
# }
# 
# if (!is.null(perm_adj)) {
#   adj_df <- as.data.frame(perm_adj)
#   adj_df$term <- rownames(adj_df)
#   rownames(adj_df) <- NULL
#   write.csv(adj_df, out_csv_adj, row.names = FALSE)
# }
