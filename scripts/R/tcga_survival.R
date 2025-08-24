#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(survival); library(survminer) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript tcga_survival.R <expr_csv> <clinical_csv> <outdir> <cancer_label>")
expr <- read_csv(args[1])  # patient, score
clin <- read_csv(args[2])  # patient, time, status(1=death)
outdir <- args[3]; cancer <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
merged <- inner_join(expr, clin, by="patient")
merged$group <- ifelse(merged$score >= median(merged$score, na.rm=TRUE), "High", "Low")
fit <- survfit(Surv(time, status) ~ group, data=merged)
p <- ggsurvplot(fit, data=merged, risk.table=TRUE, pval=TRUE, ggtheme=theme_minimal(), title=paste0("TCGA ", cancer, ": Aging signature"))
ggsave(file.path(outdir, paste0("KM_", cancer, ".pdf")), p$plot, width=6, height=5)
