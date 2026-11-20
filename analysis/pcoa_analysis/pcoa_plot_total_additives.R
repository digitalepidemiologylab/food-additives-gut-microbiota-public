library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -------------------------
# Paths
# -------------------------
meta_path <- "../../data/metadata_files/meta_zscores.csv"
distance_matrix_path <- "~/WORK/dietary-microbiota-paper/qiime/unweighted_unifrac_distance_matrix/distance-matrix.tsv"

out_dir_fig <- "../figures/Supplementary"
if (!dir.exists(out_dir_fig)) dir.create(out_dir_fig, recursive = TRUE)

out_fig_pcoa <- file.path(out_dir_fig, "pcoa_total_additive_quartiles.png")

# -------------------------
# Load metadata
# -------------------------
meta <- read.csv(meta_path, check.names = FALSE)
meta <- meta[, names(meta) != "", drop = FALSE]
meta$sample.id <- as.character(meta[["sample-id"]])

# -------------------------
# Load distance matrix + compute PCoA
# -------------------------
beta_dist_tbl <- read.table(distance_matrix_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
beta_dist <- as.dist(beta_dist_tbl)

pcoa_results <- ape::pcoa(beta_dist)
pcoa_df <- as.data.frame(pcoa_results$vectors)
pcoa_df$sample.id <- rownames(pcoa_df)

# -------------------------
# Merge and create quartiles
# -------------------------
plot_df <- pcoa_df %>%
  inner_join(meta, by = "sample.id")

if (nrow(plot_df) == 0) {
  stop("No overlapping samples between PCoA and metadata. Check sample IDs.")
}

add_col <- "total_additive_intake"
if (!(add_col %in% colnames(plot_df))) stop("Missing total_additive_intake in metadata.")

qs <- quantile(plot_df[[add_col]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# If ties cause non-unique breaks, jitter slightly
if (length(unique(qs)) < length(qs)) {
  vals2 <- plot_df[[add_col]] + rnorm(nrow(plot_df), mean = 0, sd = 1e-10)
  plot_df[[add_col]] <- vals2
  qs <- quantile(plot_df[[add_col]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
}

plot_df$additive_quartile <- cut(
  plot_df[[add_col]],
  breaks = qs,
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)
plot_df$additive_quartile <- factor(plot_df$additive_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))

######################
library(vegan)
ids <- intersect(rownames(beta_dist_tbl), plot_df$sample.id)
dist_sub <- as.dist(beta_dist_tbl[ids, ids])

meta_sub <- plot_df[match(ids, plot_df$sample.id), , drop = FALSE]

set.seed(1)
perm_res <- adonis2(dist_sub ~ additive_quartile, data = meta_sub, permutations = 999)
print(perm_res)

# perm_res_adj <- adonis2(dist_sub ~ additive_quartile + age + BMI + HEI, data = meta_sub, permutations = 999)
# print(perm_res_adj)

######################
bd <- betadisper(dist_sub, meta_sub$additive_quartile)
anova(bd)                 # overall test
permutest(bd, permutations = 999)  # permutation test
######################


medians <- plot_df %>%
  group_by(additive_quartile) %>%
  summarise(
    Axis.1 = median(Axis.1, na.rm = TRUE),
    Axis.2 = median(Axis.2, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# Plot
# -------------------------
p_a <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, color = additive_quartile)) +
  geom_point() +
  stat_ellipse(aes(fill = additive_quartile), geom = "polygon", alpha = 0.15, linetype = 3) +
  scale_color_brewer(type = "seq", palette = "OrRd", direction = 1) +
  scale_fill_brewer(type = "seq", palette = "OrRd", direction = 1) +
  geom_point(
    data = medians,
    aes(x = Axis.1, y = Axis.2, fill = additive_quartile),
    size = 3, shape = 21, color = "black", stroke = 0.5, inherit.aes = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.title.x = element_text(face = "bold", colour = "black", size = 12),
    axis.title.y = element_text(face = "bold", colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 10),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank()
  ) +
  labs(
    x = paste0("PCoA 1 (", 100 * round(pcoa_results$values$Relative_eig[1], 4), "%)"),
    y = paste0("PCoA 2 (", 100 * round(pcoa_results$values$Relative_eig[2], 4), "%)"),
    color = "Additive\nQuartile",
    fill = "Additive\nQuartile"
  )

# -------------------------
# Save
# -------------------------
ggsave(
  filename = out_fig_pcoa,
  plot = p_a,
  width = 6.1,
  height = 5,
  dpi = 300,
  units = "in"
)
