library(ggplot2)
library(dplyr)
library(ape)
library(vegan)
library(ggrepel)

# =========================
# PCoA + envfit origin arrows (Option B)
# - Unconstrained PCoA from distance matrix (eg unweighted UniFrac)
# - envfit fits continuous variables as vectors in the ordination
# - vectors are drawn from the origin (0,0) for directionality
# =========================

# -------------------------
# Paths
# -------------------------
meta_path <- "../../data/metadata_files/meta_zscores.csv"
distance_matrix_path <- "~/WORK/dietary-microbiota-paper/qiime/unweighted_unifrac_distance_matrix/distance-matrix.tsv"

out_dir_fig <- "../figures/Supplementary"
out_dir_res <- "./results"

if (!dir.exists(out_dir_fig)) dir.create(out_dir_fig, recursive = TRUE)
if (!dir.exists(out_dir_res)) dir.create(out_dir_res, recursive = TRUE)

out_fig <- file.path(out_dir_fig, "pcoa_envfit_origin_arrows.png")
out_csv_vectors <- file.path(out_dir_res, "pcoa_envfit_vectors.csv")
out_txt_envfit <- file.path(out_dir_res, "pcoa_envfit_results.txt")

# -------------------------
# Choose variables to fit as arrows (must exist in meta, preferably numeric)
# Add or remove here
# -------------------------
vars_for_vectors <- c(
  "thickeners_stabilizers_emulsifiers_gelling_agents",
  "sweeteners_glazing_agents",
  "preservatives_antioxidants",
  "ph_regulators_anti_caking_leavening_agents",
  # "other_additives",
  # "flavor_enhancers",
  # "colorants",
  # "total_additive_intake",
  "Mixture1",
  "Mixture2",
  "Mixture3",
  "HEI",
  "fruits_fg_eaten", 
  "vegetables_fg_eaten",
  "meat_fg_eaten",
  "fast_food_fg_eaten",
  "alcohol_fg_eaten"
)

# -------------------------
# Load metadata
# -------------------------
meta <- read.csv(meta_path, check.names = FALSE)
meta <- meta[, names(meta) != "", drop = FALSE]
meta$sample.id <- as.character(meta[["sample-id"]])

# -------------------------
# Load distance matrix + PCoA
# -------------------------
beta_dist_tbl <- read.table(distance_matrix_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
beta_dist <- as.dist(beta_dist_tbl)

pcoa_results <- pcoa(beta_dist)
pcoa_df <- as.data.frame(pcoa_results$vectors)
pcoa_df$sample.id <- rownames(pcoa_df)

# -------------------------
# Merge PCoA with metadata and keep only valid rows for envfit
# -------------------------
plot_df <- pcoa_df %>%
  inner_join(meta, by = "sample.id")

if (nrow(plot_df) == 0) stop("No overlapping samples between PCoA and metadata. Check sample IDs.")

# Keep only variables present in metadata
vars_present <- vars_for_vectors[vars_for_vectors %in% colnames(plot_df)]
vars_missing <- setdiff(vars_for_vectors, vars_present)

if (length(vars_missing) > 0) {
  cat("WARNING: Missing variables (will be skipped):\n")
  print(vars_missing)
}
if (length(vars_present) == 0) stop("None of vars_for_vectors are present in metadata after merge.")

# Coerce selected vars to numeric and drop rows with missing values
for (v in vars_present) plot_df[[v]] <- suppressWarnings(as.numeric(plot_df[[v]]))

keep <- is.finite(plot_df$Axis.1) & is.finite(plot_df$Axis.2)
for (v in vars_present) keep <- keep & is.finite(plot_df[[v]])
plot_fit <- plot_df[keep, , drop = FALSE]

if (nrow(plot_fit) < 10) stop("Too few complete cases after filtering for envfit. Check missingness in vars.")

# -------------------------
# envfit: fit variables as vectors in PCoA space
# -------------------------
ord_mat <- as.matrix(plot_fit[, c("Axis.1", "Axis.2")])
fit_mat <- as.data.frame(plot_fit[, vars_present, drop = FALSE])

ef <- envfit(ord_mat, fit_mat, permutations = 999)

# Save envfit output
sink(out_txt_envfit)
print(ef)
sink()

# Extract vectors (these are already coordinates in ordination space)
vec_mat <- scores(ef, display = "vectors")
vec_df <- as.data.frame(vec_mat)
colnames(vec_df)[1:2] <- c("PCoA1", "PCoA2")
vec_df$var <- rownames(vec_df)
vec_df$p <- ef$vectors$pvals
vec_df$r2 <- ef$vectors$r

# Rename variables for nicer labels
var_labels <- c(
  "thickeners_stabilizers_emulsifiers_gelling_agents" = "TSEGs",
  "sweeteners_glazing_agents" = "SG",
  "preservatives_antioxidants" = "PA",
  "ph_regulators_anti_caking_leavening_agents" = "pH/AC/L",
  "other_additives" = "Other additives",
  "flavor_enhancers" = "Flavor enhancers",
  "colorants" = "Colorants",
  "total_additive_intake" = "Total additives",
  "fruits_fg_eaten" = "Fruits",
  "vegetables_fg_eaten" = "Vegetables",
  "meat_fg_eaten" = "Meats",
  "fast_food_fg_eaten" = "Fast foods",
  "alcohol_fg_eaten" = "Alcohol"
)

vec_df$var <- ifelse(vec_df$var %in% names(var_labels), 
                     var_labels[vec_df$var], 
                     vec_df$var)

# Save vectors
write.csv(vec_df, out_csv_vectors, row.names = FALSE)

# -------------------------
# Plot: points + origin arrows
# -------------------------
p <- ggplot(plot_fit, aes(x = Axis.1, y = Axis.2)) +
  geom_point(size = 1.2, alpha = 0.3, color = "gray60") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.title.x = element_text(face = "bold", colour = "black", size = 12),
    axis.title.y = element_text(face = "bold", colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank()
  ) +
  labs(
    x = paste0("PCoA 1 (", 100 * round(pcoa_results$values$Relative_eig[1], 4), "%)"),
    y = paste0("PCoA 2 (", 100 * round(pcoa_results$values$Relative_eig[2], 4), "%)")
  )

# Visual scaling for arrows. Increase mult if arrows too short
mult <- 1

p2 <- p +
  geom_segment(
    data = vec_df,
    aes(x = 0, y = 0, xend = PCoA1 * mult, yend = PCoA2 * mult),
    inherit.aes = FALSE,
    linewidth = 0.4,
    arrow = arrow(length = unit(0.1, "inches"))
  ) +
  geom_text_repel(
    data = vec_df,
    aes(x = PCoA1 * mult, y = PCoA2 * mult, label = var),
    inherit.aes = FALSE,
    size = 5,
    point.padding = 0.3
  )

ggsave(
  filename = out_fig,
  plot = p2,
  width = 5.1,
  height = 5,
  dpi = 300,
  units = "in"
)

print(p2)
