library(spacedeconv)
library(ggplot2)
library(patchwork)

# Run from the repository root after 19_fig_S8_S9_S10_prepare.R.
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

prepared <- readRDS(file.path(objects_dir, "fig_S8_S9_S10_prepared.rds"))
labels <- c(rectangle = "Rectangle", cell2location = "cell2location",
            spatialdwls = "spatialDWLS", rctd = "RCTD")
caf_fractions <- prepared$caf_fractions[, paste0(names(labels), "_CAFs"), drop = FALSE]
spe <- SpatialExperiment::SpatialExperiment(
  colData = S4Vectors::DataFrame(caf_fractions, check.names = FALSE)
)

pairs <- combn(names(labels), 2, simplify = FALSE)
plots <- lapply(pairs, function(pair) {
  spacedeconv::plot_scatter(
    spe, value1 = paste0(pair[1], "_CAFs"), value2 = paste0(pair[2], "_CAFs"),
    xlab = labels[[pair[1]]], ylab = labels[[pair[2]]], title = NULL,
    point_alpha = 0.2, point_size = 2, stats_size = 18,
    fix_coords = FALSE, coord_range = c(0, 1)
  ) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.2)) +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 17))
})

combined <- wrap_plots(plots, ncol = 2)
ggsave(file.path(figures_dir, "fig_S9.png"), combined,
       width = 14, height = 19, dpi = 600, bg = "white")
