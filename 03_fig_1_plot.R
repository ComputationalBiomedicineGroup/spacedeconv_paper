library(ggplot2)
library(gridExtra)
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)


plots <- readRDS("./export/objects/fig_1_brain_panels.rds")
plots <- plots[-5]

add_enum <- function(p, label) {
  p +
    labs(tag = paste0(label)) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 1),
      plot.title.position = "plot",
      plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 24),
      plot.tag.position = c(0, 1),
      plot.margin = margin(0, 0, 0, 0)
    )
}

add_enum_list <- function(plots) {
  Map(function(p, i) add_enum(p, LETTERS[i+1]), plots, seq_along(plots))
}

plots <- add_enum_list(plots)
plots[[1]] <- plots[[1]] + theme(plot.margin = margin(0, 250, 0, 0))

grid_layout <- rbind(
    c(1, 1, 1, 1, 2, 2, 3, 3),
    c(4, 4, 5, 5, 6, 6, 7, 7)
)

# Build the layout without creating an Rplots.pdf.
pdf(NULL)
final <- tryCatch(
    arrangeGrob(grobs = plots, layout_matrix = grid_layout),
    finally = dev.off()
)

# Assemble the overview (A) and brain panels (B–H).
convert_image <- function(args) {
    if (system2("convert", shQuote(args)) != 0L) {
        stop("ImageMagick could not assemble Figure 1.", call. = FALSE)
    }
}

brain_file <- tempfile("fig_1_brain_", fileext = ".png")
overview_file <- tempfile("fig_1_overview_resized_", fileext = ".png")
gap_file <- tempfile("fig_1_gap_", fileext = ".png")
tryCatch({
    ggsave(
        filename = brain_file,
        plot = final,
        dpi = 600,
        width = 20,
        height = 10,
        units = "in",
        bg = "white"
    )

    figure_width <- system2("identify", shQuote(c("-format", "%w", brain_file)), stdout = TRUE)
    if (length(figure_width) != 1L || !grepl("^[0-9]+$", figure_width)) {
        stop("Could not read the brain panel image width.", call. = FALSE)
    }

    convert_image(c("data/spacedeconv_overview.png", "-resize", paste0(figure_width, "x"),
                    "-background", "white", "-alpha", "remove", "-alpha", "off", overview_file))
    convert_image(c(overview_file, "-font", "DejaVu-Sans-Bold", "-pointsize", "200",
                    "-fill", "black", "-gravity", "northwest", "-annotate", "+0+0", "A", overview_file))
    convert_image(c("-size", paste0(figure_width, "x550"), "xc:white", gap_file))
    convert_image(c(overview_file, gap_file, brain_file, "-background", "white",
                    "-gravity", "center", "-append", file.path(figures_dir, "fig_1.png")))
}, finally = unlink(c(brain_file, overview_file, gap_file)))
