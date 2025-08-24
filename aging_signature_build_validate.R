#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript aging_signature_build_validate.R <merged_seurat_rds> <metadata_csv_with_age> <outdir> <supp_table_out>")
seu_path <- args[1]; meta_csv <- args[2]; outdir <- args[3]; supp_out <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seu <- readRDS(seu_path)
meta <- readr::read_csv(meta_csv)
seu$donor_age <- meta$age[match(seu$orig.ident, meta$sample)]
Idents(seu) <- "celltype"

sig_list <- list()
for (ct in unique(seu$celltype)) {
  sub <- subset(seu, idents=ct)
  young <- sub$donor_age <= median(sub$donor_age, na.rm=TRUE)
  sub$age_bin <- ifelse(young, "Young", "Old")
  Idents(sub) <- "age_bin"
  de <- FindMarkers(sub, ident.1="Old", ident.2="Young", test.use="wilcox", logfc.threshold=0.1)
  de$gene <- rownames(de); de$celltype <- ct
  sig <- dplyr::filter(de, p_val_adj < 0.05)
  sig$direction <- ifelse(sig$avg_log2FC>0, "up_with_age", "down_with_age")
  sig_list[[ct]] <- sig
}
sig_df <- dplyr::bind_rows(sig_list)
readr::write_csv(sig_df, supp_out)

seu <- AddModuleScore(seu, features = list(unique(sig_df$gene)), name = "AgingSignature")
p <- ggplot(seu@meta.data, aes(x=donor_age, y=AgingSignature1)) + geom_point(alpha=.3) + geom_smooth(method="lm") +
  labs(x="Chronological age", y="Aging signature score")
ggsave(file.path(outdir, "aging_signature_vs_age.pdf"), p, width=6, height=5)
