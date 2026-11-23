library(SpatialExperiment)
library(spacedeconv)

# Run from the repository root.
objects_dir <- "./export/objects"
second_gen <- c("spatialdwls", "cell2location", "spotlight", "rctd", "card",
                "dot", "rectangle")
methods <- c(second_gen, "epic", "quantiseq")
method_labels <- setNames(c("spatialDWLS", "cell2location", "SPOTlight", "RCTD",
                            "CARD", "DOT", "Rectangle", "EPIC", "quanTIseq"), methods)
paths <- setNames(file.path(objects_dir, paste0("fig_S8_S9_S10_", methods, ".rds")), methods)
if (any(!file.exists(paths))) {
  stop("Missing deconvolution results for: ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "), call. = FALSE)
}

# Define the expected spots and the Wu annotations
# to map each minor type to one major type.
spots <- colnames(readRDS("data/allresults_minor_andersson.rds"))
reference <- readRDS("data/wu_reference_filtered.rds")
types <- unique(data.frame(minor = as.character(reference$celltype_minor),
                           major = as.character(reference$celltype_major)))
rm(reference)
stopifnot(length(spots) == 3798L, !anyNA(spots), !anyDuplicated(spots),
          !anyNA(types), !anyDuplicated(make.names(types$minor)))
types$minor <- make.names(types$minor)
types$major <- make.names(types$major)

categories <- c(B.cells = "B cells", CAFs = "CAFs",
                Cancer.Epithelial = "Cancer Epithelial", Endothelial = "Endothelial",
                Myeloid = "Myeloid", Plasmablasts = "Plasmablasts",
                PVL = "PVL", T.cells = "T/NK cells")
stopifnot(setequal(types$major, names(categories)))

read_result <- function(method) {
  result <- readRDS(paths[[method]])
  if (!is.matrix(result) || !is.numeric(result) || any(!is.finite(result))) {
    stop(method, ": expected a finite numeric result matrix.", call. = FALSE)
  }
  spot_ids <- rownames(result)
  if (anyDuplicated(spot_ids) || !setequal(spot_ids, spots)) {
    stop(method, ": missing, unexpected or duplicate spot IDs.", call. = FALSE)
  }
  if (is.null(colnames(result)) || anyNA(colnames(result)) ||
      anyDuplicated(colnames(result))) {
    stop(method, ": missing or duplicate cell-type column names.", call. = FALSE)
  }
  if (any(result < -1e-8)) {
    stop(method, ": negative deconvolution estimates.", call. = FALSE)
  }
  result[result < 0] <- 0
  result[spots, , drop = FALSE]
}

aggregate_groups <- function(result, groups, method) {
  missing <- setdiff(unlist(groups, use.names = FALSE), colnames(result))
  if (length(missing)) {
    stop(method, ": missing result columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  spe <- SpatialExperiment(colData = S4Vectors::DataFrame(result, check.names = FALSE))
  values <- lapply(seq_along(groups), function(i) {
    columns <- groups[[i]]
    if (length(columns) == 1L) {
      return(result[, columns])
    }
    group <- spacedeconv::aggregate_results(
      spe, cell_types = columns, name = ".aggregated"
    )
    colData(group)[[".aggregated"]]
  })
  out <- do.call(cbind, values)
  dimnames(out) <- list(rownames(result), names(groups))
  out
}

second_fractions <- lapply(second_gen, function(method) {
  result <- read_result(method)
  groups <- split(paste0(method, "_", types$minor),
                  factor(types$major, levels = names(categories)))
  expected <- c(unlist(groups, use.names = FALSE),
                if (method == "rectangle") "rectangle_Unknown")
  extra <- setdiff(colnames(result), expected)
  if (length(extra)) {
    stop(method, ": unexpected result columns: ", paste(extra, collapse = ", "),
         call. = FALSE)
  }
  missing <- setdiff(expected, colnames(result))
  if (length(missing)) {
    stop(method, ": missing result columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  totals <- rowSums(result) # Includes Rectangle Unknown.
  if (any(!is.finite(totals) | totals <= 0)) {
    stop(method, ": invalid total abundance for at least one spot.", call. = FALSE)
  }
  names(groups) <- paste0(method, "_", names(groups))
  aggregate_groups(result, groups, method) / totals
})

# Match only the predefined signature populations used in the figures.
epic_groups <- list(
  epic_B.cells = "epic_B cell",
  epic_CAFs = "epic_Cancer associated fibroblast",
  epic_Endothelial = "epic_Endothelial cell",
  epic_T.cells = c("epic_T cell CD4+", "epic_T cell CD8+", "epic_NK cell")
)
quantiseq_groups <- list(
  quantiseq_B.cells = "quantiseq_B cell",
  quantiseq_Myeloid = c("quantiseq_Macrophage M1", "quantiseq_Macrophage M2",
    "quantiseq_Monocyte", "quantiseq_Neutrophil", "quantiseq_Myeloid dendritic cell"),
  quantiseq_T.cells = c("quantiseq_T cell CD4+ (non-regulatory)", "quantiseq_T cell CD8+",
    "quantiseq_T cell regulatory (Tregs)", "quantiseq_NK cell")
)
first_groups <- list(epic = epic_groups, quantiseq = quantiseq_groups)
first_fractions <- lapply(names(first_groups), function(method) {
  result <- read_result(method)
  if (any(abs(rowSums(result) - 1) > 1e-6)) {
    stop(method, ": native fractions do not sum to one.", call. = FALSE)
  }
  aggregate_groups(result, first_groups[[method]], method)
})

fractions <- do.call(cbind, c(second_fractions, first_fractions))
columns <- colnames(fractions)
method <- sub("_.*$", "", columns)
category <- unname(categories[sub("^[^_]+_", "", columns)])
stopifnot(!anyNA(category))
metadata <- data.frame(column = columns, method = method,
                       method_label = unname(method_labels[method]),
                       category = category)
metadata <- metadata[order(match(metadata$category, categories),
                           match(metadata$method, methods)), , drop = FALSE]
fractions <- fractions[, metadata$column, drop = FALSE]

# S8: Pearson correlations (undefined for constant columns).
variable <- apply(fractions, 2, function(x) any(x != x[1]))
if (any(!variable)) {
  warning("Excluded constant columns from the heatmap: ",
          paste(colnames(fractions)[!variable], collapse = ", "), call. = FALSE)
}
heatmap_columns <- colnames(fractions)[variable]
correlations <- cor(fractions[, heatmap_columns, drop = FALSE], method = "pearson")

# S9: CAF fractions for the scatterplot.
caf_columns <- paste0(c("rectangle", "cell2location", "spatialdwls", "rctd"), "_CAFs")
caf_fractions <- fractions[, caf_columns, drop = FALSE]

# S10: mean spot fractions per method and category.
boxplot_means <- transform(metadata, fraction = unname(colMeans(fractions)))

saveRDS(list(fractions = fractions, metadata = metadata, correlations = correlations,
             caf_fractions = caf_fractions, boxplot_means = boxplot_means),
        file.path(objects_dir, "fig_S8_S9_S10_prepared.rds"))
