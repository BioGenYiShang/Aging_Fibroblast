#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(stringr); library(ggplot2); library(patchwork); library(plyr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript seurat_annotation.R <input_mtx_or_h5> <sample_name> <outdir>")
input <- args[1]; sample <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

obj <- if (grepl("\.h5$", input)) Read10X_h5(input) else Read10X(input)
seu <- CreateSeuratObject(counts = obj, project = sample, min.cells = 3, min.features = 200)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 10)
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
seu <- FindNeighbors(seu, dims = 1:30) %>% FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:30)

markers <- list(
  Keratinocyte = c("KRT14","KRT5"),
  Fibroblast   = c("COL1A1","COL1A2","DCN"),
  Endothelial  = c("PECAM1","VWF"),
  Tcell        = c("CD3D","TRAC"),
  Bcell        = c("MS4A1","CD79A"),
  Myeloid      = c("LYZ","CSF1R")
)
Idents(seu) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.25)
write_csv(cluster_markers, file.path(outdir, paste0(sample, "_cluster_markers.csv")))

cluster_to_label <- function(markers, cm) {
  labs <- rep("Unknown", length(unique(cm$cluster)))
  for (cid in unique(cm$cluster)) {
    top <- cm %>% filter(cluster == cid) %>% arrange(desc(avg_log2FC)) %>% head(200) %>% pull(gene)
    for (lbl in names(markers)) {
      if (length(intersect(top, markers[[lbl]]))>0) { labs[as.integer(cid)+1] <- lbl; break }
    }
  }
  setNames(labs, unique(cm$cluster))
}
labs <- cluster_to_label(markers, cluster_markers)
seu$celltype <- plyr::mapvalues(as.character(seu$seurat_clusters), from=names(labs), to=as.vector(labs))
saveRDS(seu, file.path(outdir, paste0(sample, "_seurat.rds")))
pdf(file.path(outdir, paste0(sample, "_umap_celltype.pdf")), width=6, height=5)
print(DimPlot(seu, group.by="celltype", label=TRUE) + ggtitle(sample))
dev.off()
