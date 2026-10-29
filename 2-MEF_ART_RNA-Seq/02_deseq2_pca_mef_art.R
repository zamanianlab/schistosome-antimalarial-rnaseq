# 02_deseq2_pca_mef_art.R - DESeq2, MEF and ART vs DMSO; Supp Fig 1 PCA; Fig 3A Venn diagram
# in : output/mef_art_counts_raw.rds, ../data_files/sample_metadata.csv
# out: output/malarial_vivo_adult.MEF_DMSO.DEG.csv, .ART_DMSO.DEG.csv (Supp File 2, sheets 3-4),
#      SuppFig1_PCA.pdf, Fig3A_venn_MEF_ART.pdf

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "2-MEF_ART_RNA-Seq"))
data_dir <- file.path(code_dir, "data_files")
cohort   <- "malarial_vivo_adult"

library(tidyverse)
library(DESeq2)
library(ggVennDiagram)
library(sf)

samples <- read.csv(file.path(data_dir, "sample_metadata.csv")) %>%
  filter(project_batch == cohort) %>%
  select(sample, condition, project_batch) %>%
  column_to_rownames("sample")
samples$condition <- relevel(factor(samples$condition), ref = "DMSO")

gene_count <- readRDS("output/mef_art_counts_raw.rds") %>%
  filter(sample %in% rownames(samples)) %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
gene_count <- gene_count[, rownames(samples)]
stopifnot(identical(colnames(gene_count), rownames(samples)))

dds  <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 5          # >= 10 reads in >= 5 samples
dds  <- dds[keep, ]
message("genes tested: ", nrow(dds))
dds  <- DESeq(dds)

## DEG tables
get_res <- function(drug) {
  results(dds, contrast = c("condition", drug, "DMSO")) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)
}
res <- list(MEF = get_res("MEF"), ART = get_res("ART"))
for (d in names(res)) {
  write.csv(res[[d]], paste0("output/malarial_vivo_adult.", d, "_DMSO.DEG.csv"),
            quote = FALSE, row.names = FALSE)
  deg <- filter(res[[d]], !is.na(padj), padj < 0.01, abs(log2FoldChange) > 1)
  message(d, " DEGs (padj<0.01, |log2FC|>1): ", nrow(deg),
          "  up ", sum(deg$log2FoldChange > 0), "  down ", sum(deg$log2FoldChange < 0))
}

## Supp Fig 1: PCA
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(size = 2, alpha = 0.75, aes(shape = condition)) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 8), axis.text = element_text(size = 7.5),
        legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
        legend.margin = margin(-0.3, -0.5, 0, -0.5, "cm"))
ggsave("output/SuppFig1_PCA.pdf",  p_pca, width = 3.5, height = 3.7, units = "in")

## Fig 3A: Venn diagram
up   <- function(d) d %>% filter(!is.na(padj), padj < 0.01, log2FoldChange >=  1) %>% pull(gene_id)
down <- function(d) d %>% filter(!is.na(padj), padj < 0.01, log2FoldChange <= -1) %>% pull(gene_id)
venn_list <- list("MEF up" = up(res$MEF), "ART up" = up(res$ART),
                  "MEF down" = down(res$MEF), "ART down" = down(res$ART))

venn.plot <- ggVennDiagram(venn_list, label = "both", label_alpha = 0,
                           set_size = 4, edge_size = 0.6) +
  scale_fill_gradient(low = "white", high = "grey") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +   # room for set labels
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "pt"))
ggsave("output/Fig3A_venn_MEF_ART.pdf",  venn.plot, width = 5, height = 5, units = "in")
