library(spacedeconv)
library(SpatialExperiment)
library(gridExtra)
library(ggplot2)
library(grid)
library(ggpubr)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


allresults_minor_1142243F <- readRDS("./data/allresults_minor_1142243F.rds")
allresults_minor_1160920F <- readRDS("./data/allresults_minor_1160920F.rds")
allresults_minor_4465 <- readRDS("./data/allresults_minor_4465.rds")
allresults_minor_44971 <- readRDS("./data/allresults_minor_44971.rds")
allresults_minor_4290 <- readRDS("./data/allresults_minor_4290.rds") # lumA
allresults_minor_4535 <- readRDS("./data/allresults_minor_4535.rds") # lumB

title_size <- 0 #not needed
font_size <- 18
legend_size <- 20

# most abundant cancer cell types, with density
x1 <- plot_spatial(
  allresults_minor_1142243F,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x2 <- plot_spatial(
  allresults_minor_1142243F,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x3 <- plot_spatial(
  allresults_minor_1142243F,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x4 <- plot_spatial(
  allresults_minor_1142243F,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x5 <- plot_spatial(
  allresults_minor_1160920F,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x6 <- plot_spatial(
  allresults_minor_1160920F,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x7 <- plot_spatial(
  allresults_minor_1160920F,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x8 <- plot_spatial(
  allresults_minor_1160920F,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x9 <- plot_spatial(
  allresults_minor_4465,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x10 <- plot_spatial(
  allresults_minor_4465,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x11 <- plot_spatial(
  allresults_minor_4465,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x12 <- plot_spatial(
  allresults_minor_4465,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x13 <- plot_spatial(
  allresults_minor_44971,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x14 <- plot_spatial(
  allresults_minor_44971,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x15 <- plot_spatial(
  allresults_minor_44971,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x16 <- plot_spatial(
  allresults_minor_44971,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x17 <- plot_spatial(
  allresults_minor_4290,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x18 <- plot_spatial(
  allresults_minor_4290,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x19 <- plot_spatial(
  allresults_minor_4290,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x20 <- plot_spatial(
  allresults_minor_4290,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x21 <- plot_spatial(
  allresults_minor_4535,
  result = "rctd_Cancer.Basal.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)
x22 <- plot_spatial(
  allresults_minor_4535,
  result = "rctd_Cancer.LumA.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x23 <- plot_spatial(
  allresults_minor_4535,
  result = "rctd_Cancer.LumB.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

x24 <- plot_spatial(
  allresults_minor_4535,
  result = "rctd_Cancer.Her2.SC",
  density = FALSE,
  title_size = title_size,
  font_size = font_size,
  legend_size = legend_size,
  smooth = TRUE,
  limits = c(0, 1),
  title = ""
)

plots <- list(
  x1, x2, x3, x4,
  x5, x6, x7, x8,
  x9, x10, x11, x12,
  x13, x14, x15, x16,
  x17, x18, x19, x20,
  x21, x22, x23, x24
)

column_labels <- c("Cancer Basal SC", "Cancer LumA SC", "Cancer LumB SC", "Cancer Her2 SC")
sample_ids <- c("1142243F", "1160920F", "CID4465", "CID44971", "CID4290", "CID4535")
sample_subtypes <- c("Basal", "Basal", "Basal", "Basal", "LumA", "LumB")

column_headers <- lapply(column_labels, function(label) {
  textGrob(label, x = 0.45, gp = gpar(fontsize = 32))
})
row_headers <- Map(function(id, subtype) {
  grobTree(
    textGrob(id, x = 0.37, rot = 90, gp = gpar(fontsize = 32, fontface = "bold")),
    textGrob(paste0("(", subtype, ")"), x = 0.67, rot = 90, gp = gpar(fontsize = 32))
  )
}, sample_ids, sample_subtypes)

grid_layout <- rbind(
  c(NA, 25:28),
  cbind(29:34, matrix(seq_along(plots), nrow = 6, byrow = TRUE))
)

# Build the layout without opening a file-backed graphics device.
pdf(NULL)
grid_1 <- tryCatch(
  arrangeGrob(
    grobs = c(plots, column_headers, row_headers),
    layout_matrix = grid_layout,
    widths = unit.c(unit(1.8, "in"), unit(rep(1, 4), "null")),
    heights = unit.c(unit(0.6, "in"), unit(rep(1, 6), "null"))
  ),
  finally = dev.off()
)

ggsave(
  "./export/figures/fig_S3.png",
  plot = grid_1,
  dpi = 600,
  width = 21.8,
  height = 30.6,
  units = "in",
  bg = "white"
)
