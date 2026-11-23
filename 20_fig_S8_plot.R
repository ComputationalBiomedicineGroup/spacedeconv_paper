library(spacedeconv)
library(ComplexHeatmap)
library(grid)

# Run from the repository root after 19_fig_S8_S9_S10_prepare.R.
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

prepared <- readRDS(file.path(objects_dir, "fig_S8_S9_S10_prepared.rds"))
correlations <- prepared$correlations
columns <- colnames(correlations)
metadata <- prepared$metadata[match(columns, prepared$metadata$column), , drop = FALSE]
fractions <- prepared$fractions[, columns, drop = FALSE]
stopifnot(identical(metadata$column, columns), identical(rownames(correlations), columns))

spe <- SpatialExperiment::SpatialExperiment(
  colData = S4Vectors::DataFrame(fractions, check.names = FALSE)
)
category_colors <- setNames(
  grDevices::palette.colors(10, "Tableau 10")[c(1, 8, 2, 3, 6, 7, 4, 5, 9, 10)],
  c("B cells", "Plasmablasts", "CAFs", "Cancer Epithelial", "Normal Epithelial",
    "Other", "Endothelial", "Myeloid", "PVL", "T/NK cells")
)
categories <- setNames(metadata$category, columns)
labels <- setNames(paste(metadata$method_label, metadata$category, sep = " | "), columns)

plot_heatmap <- function() {
  on.exit(invisible(dev.off()))
  old_options <- ht_opt()[c("HEATMAP_LEGEND_PADDING", "ANNOTATION_LEGEND_PADDING", "legend_gap")]
  on.exit(ht_opt(old_options), add = TRUE)
  ht_opt(HEATMAP_LEGEND_PADDING = unit(8, "mm"),
         ANNOTATION_LEGEND_PADDING = unit(8, "mm"),
         legend_gap = unit(c(8, 18), "mm"))

  result <- spacedeconv::spatialcorr(
    spe, column = columns, categories = categories,
    category_colors = category_colors[unique(categories)],
    labels = labels, label_title = "Method | Cell type",
    label_title_gp = gpar(fontsize = 20, fontface = "bold"),
    row_names_gp = gpar(fontsize = 17),
    width = unit(11.2, "inches"), height = unit(11.2 * (1 + sqrt(5)) / 2, "inches"),
    title = NULL,
    draw_param = list(heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                      padding = unit(c(5, 2, 2, 10), "mm")),
    annotation_legend_param = list(
      nrow = 2, by_row = TRUE, title_position = "topleft",
      column_gap = unit(5, "mm"), row_gap = unit(2, "mm"),
      grid_width = unit(5, "mm"), grid_height = unit(5, "mm"),
      title_gp = gpar(fontsize = 26, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    ),
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1), direction = "horizontal",
      legend_width = unit(6, "cm"), title_position = "topcenter",
      title_gp = gpar(fontsize = 26, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    )
  )
  stopifnot(isTRUE(all.equal(result$corr, correlations, tolerance = 1e-12)),
            identical(result$used_spots, rownames(fractions)))
}

pdf(file.path(figures_dir, "fig_S8.pdf"), width = 15.65, height = 21.25, useDingbats = FALSE)
plot_heatmap()
png(file.path(figures_dir, "fig_S8.png"), width = 4695, height = 6375, res = 300)
plot_heatmap()
