library(SpatialExperiment)

# Run from the repository root.
figures_dir <- "./export/figures"
objects_dir <- "./export/objects"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)

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

# Compute unlogged CPM using all Visium genes for both signature-based methods.
spe_first <- spacedeconv::normalize(spe, method = "cpm")

# EPIC
epic <- as.matrix(spacedeconv::deconvolute(
  spe_first, method = "epic", assay_sp = "cpm",
  tumor = TRUE, scale_mrna = TRUE, return_object = FALSE
))
stopifnot(identical(rownames(epic), colnames(spe_first)), all(is.finite(epic)))
saveRDS(epic, file.path(objects_dir, "fig_S8_epic.rds"))

# quanTIseq
quantiseq <- as.matrix(spacedeconv::deconvolute(
  spe_first, method = "quantiseq", assay_sp = "cpm",
  tumor = TRUE, arrays = FALSE, scale_mrna = TRUE, return_object = FALSE
))
stopifnot(identical(rownames(quantiseq), colnames(spe_first)),
          all(is.finite(quantiseq)))
saveRDS(quantiseq, file.path(objects_dir, "fig_S8_quantiseq.rds"))

sce$celltype_minor <- make.names(as.character(sce$celltype_minor))

# RCTD
set.seed(42)
rctd <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "rctd",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", n_cores = 16, return_object = FALSE
))
saveRDS(rctd, file.path(objects_dir, "fig_S8_rctd.rds"))

# CARD
set.seed(42)
card <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "card",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", return_object = FALSE
))
saveRDS(card, file.path(objects_dir, "fig_S8_card.rds"))

# DOT
set.seed(42)
dot <- as.matrix(spacedeconv::deconvolute(
  spe_shared, single_cell_obj = sce, method = "dot",
  cell_type_col = "celltype_minor", batch_id_col = "orig.ident",
  assay_sp = "counts", max_spot_size = 20, ratios_weight = 0,
  return_object = FALSE
))
saveRDS(dot, file.path(objects_dir, "fig_S8_dot.rds"))

# Rectangle
set.seed(42)
rectangle_model <- spacedeconv::build_model(
  sce, spatial_obj = spe_shared, method = "rectangle",
  cell_type_col = "celltype_minor", n_cpus = 32
)
rectangle <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = rectangle_model, method = "rectangle",
  assay_sp = "counts", n_cpus = 32, return_object = FALSE
))
saveRDS(rectangle, file.path(objects_dir, "fig_S8_rectangle.rds"))

# spatialDWLS (build_model uses scran for marker selection).
set.seed(42)
spatialdwls_model <- spacedeconv::build_model(
  sce, method = "spatialdwls", cell_type_col = "celltype_minor",
  topNgenes = 100L
)
spatialdwls <- as.matrix(spacedeconv::deconvolute(
  spe_shared, signature = spatialdwls_model, method = "spatialdwls",
  assay_sp = "counts", return_object = FALSE
))
saveRDS(spatialdwls, file.path(objects_dir, "fig_S8_spatialdwls.rds"))
