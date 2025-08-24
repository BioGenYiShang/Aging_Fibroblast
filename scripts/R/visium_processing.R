#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(SummarizedExperiment)
  library(scater); library(scran); library(ggplot2); library(readr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript visium_processing.R <spaceranger_out> <sample_name> <outdir>")
in_dir <- args[1]; sample <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seu <- Load10X_Spatial(in_dir)
sce <- as.SingleCellExperiment(seu)
qc <- perCellQCMetrics(sce)
keep <- qc$detected > 200 & sce$sum > 500 & qc$subsets_Mito_percent < 15
sce <- sce[, keep]
sce <- logNormCounts(sce)
saveRDS(sce, file.path(outdir, paste0(sample, "_visium_sce.rds")))
pdf(file.path(outdir, paste0(sample, "_visium_qc.pdf")), width=6, height=5)
print(plotColData(sce, "sum") + ggtitle("UMI per spot"))
dev.off()
