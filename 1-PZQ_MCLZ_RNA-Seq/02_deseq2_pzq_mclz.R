# 02_deseq2_pzq_mclz.R - DESeq2, PZQ vs DMSO and MCLZ vs DMSO (one model per cohort)
# in : output/pzq_counts_raw.rds, mclz_counts_raw.rds, ../data_files/sample_metadata.csv
# out: output/pzq_vivo_adult.PZQ_DMSO.DEG.csv (Supp File 2, sheet 1)
#      output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv (Supp File 2, sheet 2)

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "1-PZQ_MCLZ_RNA-Seq"))
data_dir <- file.path(code_dir, "data_files")

library(tidyverse)
library(DESeq2)

run_deseq <- function(stub, cohort, drug) {
  samples <- read.csv(file.path(data_dir, "sample_metadata.csv")) %>%
    filter(project_batch == cohort) %>%
    select(sample, condition, project_batch, stage) %>%
    column_to_rownames("sample")
  samples$condition <- relevel(factor(samples$condition), ref = "DMSO")

  gene_count <- readRDS(paste0("output/", stub, "_counts_raw.rds")) %>%
    filter(sample %in% rownames(samples)) %>%
    pivot_wider(names_from = sample, values_from = expression) %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  gene_count <- gene_count[, rownames(samples)]
  stopifnot(identical(colnames(gene_count), rownames(samples)))

  dds  <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ condition)
  keep <- rowSums(counts(dds) >= 10) >= 5          # >= 10 reads in >= 5 samples
  dds  <- dds[keep, ]
  dds  <- DESeq(dds)

  res <- results(dds, contrast = c("condition", drug, "DMSO")) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)
  write.csv(res, paste0("output/", cohort, ".", drug, "_DMSO.DEG.csv"), quote = FALSE, row.names = FALSE)

  deg <- filter(res, !is.na(padj), padj < 0.01, abs(log2FoldChange) > 1)
  message(drug, ": genes tested ", nrow(dds), "; DEGs (padj<0.01, |log2FC|>1) ", nrow(deg),
          ", up ", sum(deg$log2FoldChange > 0), ", down ", sum(deg$log2FoldChange < 0))
}

run_deseq("pzq",  "pzq_vivo_adult",  "PZQ")
run_deseq("mclz", "mclz_vivo_adult", "MCLZ")
