#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(CellChat); library(patchwork); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript cellchat_analysis.R <seurat_rds> <group_by_age_csv> <outdir>")
seu <- readRDS(args[1]); group_csv <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
meta_age <- read.csv(group_csv)
seu$AgeGroup <- meta_age$age_group[match(seu$orig.ident, meta_age$sample)]
mat <- GetAssayData(seu, assay="RNA", slot="data")
meta <- data.frame(group=seu$celltype, row.names=colnames(seu))
cellchat <- createCellChat(object = mat, meta = meta, group.by = "group")
CellChatDB <- CellChatDB.human; cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat); cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat); cellchat <- filterCommunication(cellchat, min.cells=10)
cellchat <- computeCommunProbPathway(cellchat); cellchat <- aggregateNet(cellchat)
pdf(file.path(outdir, "net_circular.pdf"), width=7, height=7)
netVisual_circle(cellchat@net$count, weight.edge = TRUE, label.edge= FALSE)
dev.off()
saveRDS(cellchat, file.path(outdir, "cellchat.rds"))
