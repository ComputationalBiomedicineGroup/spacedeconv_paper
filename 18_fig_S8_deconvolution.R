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
