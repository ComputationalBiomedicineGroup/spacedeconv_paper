library(spacedeconv)
library(ggplot2)

# Run from the repository root after 19_fig_S8_S9_S10_prepare.R.
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

prepared <- readRDS(file.path(objects_dir, "fig_S8_S9_S10_prepared.rds"))
metadata <- prepared$metadata
columns <- metadata$column
spe <- SpatialExperiment::SpatialExperiment(
  colData = S4Vectors::DataFrame(prepared$fractions, check.names = FALSE)
)
groups <- setNames(metadata$category, columns)
sources <- setNames(metadata$method_label, columns)
colors <- setNames(grDevices::palette.colors(10, "Tableau 10")[seq_len(9)],
  c("spatialDWLS", "cell2location", "SPOTlight", "RCTD", "CARD", "DOT",
    "Rectangle", "EPIC", "quanTIseq"))

plot <- spacedeconv::plot_boxplot(
  spe, column = columns, groups = groups, sources = sources,
  colors = colors, point_size = 5, ylab = "Mean estimated relative abundance"
) +
  scale_x_discrete(labels = c("Cancer Epithelial" = "Cancer\nEpithelial",
                              "Endothelial" = "Endo-\nthelial",
                              "Plasmablasts" = "Plasma-\nblasts")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.1)) +
  labs(color = "Method") +
  theme(
    text = element_text(family = "sans"),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 28),
    legend.text = element_text(size = 26),
    legend.key.height = grid::unit(0.9, "cm")
  )

expected_means <- setNames(prepared$boxplot_means$fraction, prepared$boxplot_means$column)
stopifnot(isTRUE(all.equal(plot$data$mean, unname(expected_means[plot$data$column]),
                          tolerance = 1e-12)))

ggsave(file.path(figures_dir, "fig_S10.png"), plot,
       width = 22, height = 10, dpi = 600, bg = "white")
