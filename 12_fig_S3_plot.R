library(ggplot2)
library(gridExtra)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


plots_3a <- readRDS("./export/objects/fig_S3A_panels.rds")
plots_3b <- readRDS("./export/objects/fig_S3B_panels.rds")

plots <- c(plots_3a, plots_3b)
plots <- plots[-c(5, 6, 7, 11)]
plots <- plots[-10] # Remove the COL18A1-ITGB1 panel (formerly J).

add_enum <- function(p, label) {
  p +
    labs(tag = paste0(label)) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1),
      plot.title.position = "plot",
      plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 24),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 6, 6, 6)
    )
}

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)

layout_matrix <- rbind(
  rep(1:4, each = 3),
  rep(5:7, each = 4),
  rep(8:10, each = 4),
  rep(11:13, each = 4)
)

final_plot <- grid.arrange(
  grobs = plots,
  layout_matrix = layout_matrix
)

ggplot2::ggsave(
  filename = "./export/figures/fig_S3.png",
  plot = final_plot,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600,
  bg = "white"
)
