library(SpatialExperiment)

set.seed(42)

# Run from the repository root.
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)

check_deconvolution_results <- function(result, spatial_obj, method) {
  if (!identical(rownames(result), colnames(spatial_obj))) {
    warning(method, ": unexpected spot IDs or order.",
            call. = FALSE, immediate. = TRUE)
  }
  if (any(!is.finite(result))) {
    warning(method, ": non-finite estimates.",
            call. = FALSE, immediate. = TRUE)
  }
}

# Load the same Visium section used for Fig. S8.
visium <- readRDS("data/allresults_minor_andersson.rds")
spe <- SpatialExperiment(
  assays = list(counts = assay(visium, "counts")),
  spatialCoords = spatialCoords(visium),
  imgData = imgData(visium)
)
rm(visium)

# Filtered Wu et al. single-cell reference.
sce <- readRDS("data/wu_reference_filtered.rds")

# Match gene symbols and retain genes detected in both datasets.
rownames(spe) <- toupper(rownames(spe))
stopifnot(
  !anyDuplicated(rownames(spe)),
  !anyDuplicated(rownames(sce)),
  !anyNA(sce$orig.ident),
  !anyNA(sce$celltype_minor),
  !anyNA(sce$celltype_major)
)

genes <- intersect(rownames(spe), rownames(sce))
detected_in_both <- Matrix::rowSums(assay(spe[genes, ], "counts")) > 0 &
  Matrix::rowSums(assay(sce[genes, ], "counts")) > 0
genes <- genes[detected_in_both]

# Keep the full spatial counts in spe for EPIC and quanTIseq.
spe_shared <- spe[genes, ]
sce <- sce[genes, ]
stopifnot(
  identical(rownames(spe_shared), rownames(sce)),
  all(Matrix::colSums(assay(spe, "counts")) > 0),
  all(Matrix::colSums(assay(sce, "counts")) > 0)
)

sce$celltype_minor <- make.names(as.character(sce$celltype_minor))

# Following the SPOTlight tutorial, select features on the full reference and then subsample cells.
# Do this before loading spacedeconv to avoid a scran method-dispatch conflict
# in the R environment.
spotlight_features <- local({
  sce_log <- scuttle::logNormCounts(sce)
  eligible <- !grepl("^(RPL|RPS|MT-)", rownames(sce_log))
  hvg <- scran::getTopHVGs(
    scran::modelGeneVar(sce_log, subset.row = eligible), n = 3000
  )
  scores <- scran::scoreMarkers(
    assay(sce_log, "logcounts"), groups = sce_log$celltype_minor,
    subset.row = eligible, BPPARAM = BiocParallel::SerialParam()
  )
  markers <- do.call(rbind, lapply(names(scores), function(cell_type) {
    ranked <- scores[[cell_type]]
    ranked <- ranked[is.finite(ranked$mean.AUC) & ranked$mean.AUC > 0.8, ]
    ranked <- ranked[order(ranked$mean.AUC, decreasing = TRUE), ]
    data.frame(gene = rownames(ranked), cluster = rep(cell_type, nrow(ranked)),
               mean.AUC = ranked$mean.AUC)
  }))
  list(markers = markers, hvg = hvg)
})
stopifnot(setequal(spotlight_features$markers$cluster, sce$celltype_minor))

library(spacedeconv)

# Compute unlogged CPM using all Visium genes for both signature-based methods.
spe_first <- spacedeconv::normalize(spe, method = "cpm")

# EPIC
epic <- as.matrix(spacedeconv::deconvolute(
  spe_first, method = "epic", assay_sp = "cpm",
  tumor = TRUE, scale_mrna = TRUE, return_object = FALSE
))
check_deconvolution_results(epic, spe_first, "epic")
saveRDS(epic, file.path(objects_dir, "fig_S8_S9_S10_epic.rds"))

# quanTIseq
quantiseq <- as.matrix(spacedeconv::deconvolute(
  spe_first, method = "quantiseq", assay_sp = "cpm",
  tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, return_object = FALSE
))
check_deconvolution_results(quantiseq, spe_first, "quantiseq")
saveRDS(quantiseq, file.path(objects_dir, "fig_S8_S9_S10_quantiseq.rds"))

# RCTD
rctd <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "rctd",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", n_cores = 16, return_object = FALSE
))
check_deconvolution_results(rctd, spe_shared, "rctd")
saveRDS(rctd, file.path(objects_dir, "fig_S8_S9_S10_rctd.rds"))

# CARD
card <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "card",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", return_object = FALSE
))
check_deconvolution_results(card, spe_shared, "card")
saveRDS(card, file.path(objects_dir, "fig_S8_S9_S10_card.rds"))

# DOT
dot <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "dot",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", max_spot_size = 20, ratios_weight = 0,
  return_object = FALSE
))
check_deconvolution_results(dot, spe_shared, "dot")
saveRDS(dot, file.path(objects_dir, "fig_S8_S9_S10_dot.rds"))

# Rectangle
rectangle_model <- spacedeconv::build_model(
  sce, spatial_obj = spe_shared, method = "rectangle",
  cell_type_col = "celltype_minor", n_cpus = 16
)
rectangle <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = rectangle_model, method = "rectangle",
  assay_sp = "counts", n_cpus = 16, return_object = FALSE
))
check_deconvolution_results(rectangle, spe_shared, "rectangle")
saveRDS(rectangle, file.path(objects_dir, "fig_S8_S9_S10_rectangle.rds"))

# spatialDWLS (build_model uses scran for marker selection).
spatialdwls_model <- spacedeconv::build_model(
  sce, method = "spatialdwls", cell_type_col = "celltype_minor",
  topNgenes = 100L
)
spatialdwls <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = spatialdwls_model, method = "spatialdwls",
  assay_sp = "counts", return_object = FALSE
))
check_deconvolution_results(spatialdwls, spe_shared, "spatialdwls")
saveRDS(spatialdwls, file.path(objects_dir, "fig_S8_S9_S10_spatialdwls.rds"))

# cell2location
cell2location_model <- spacedeconv::build_model(
  sce, method = "cell2location", cell_type_col = "celltype_minor",
  batch_id_col = "orig.ident", epochs = 250L, posterior_samples = 1000L,
  gpu = TRUE
)
cell2location <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = cell2location_model, method = "cell2location",
  epochs = 30000L, n_cell = 10, alpha = 20, posterior_samples = 1000L,
  values = "absolute", gpu = TRUE, return_object = FALSE
))
check_deconvolution_results(cell2location, spe_shared, "cell2location")
saveRDS(cell2location, file.path(objects_dir, "fig_S8_S9_S10_cell2location.rds"))

# SPOTlight
spotlight_cells <- unlist(lapply(split(seq_len(ncol(sce)), sce$celltype_minor),
  function(cells) sample(cells, min(100L, length(cells)))), use.names = FALSE)
spotlight_model <- spacedeconv::build_model(
  sce[, spotlight_cells], spatial_obj = spe_shared, method = "spotlight",
  cell_type_col = "celltype_minor", markers = spotlight_features$markers,
  hvg = spotlight_features$hvg, nrun = 1L, .options = list(parallel = FALSE)
)
spotlight <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = spotlight_model, method = "spotlight",
  assay_sp = "counts", return_object = FALSE
))
check_deconvolution_results(spotlight, spe_shared, "spotlight")
saveRDS(spotlight, file.path(objects_dir, "fig_S8_S9_S10_spotlight.rds"))
