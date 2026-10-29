# fig1b_deg_heatmap.R - Fig 1B heat map of DEGs from all four treatments
# in : DESeq2 tables from folders 1-2, ../data_files/WBP19_gene_annotation.csv
# out: output/DEG_log2FC_all_drugs.csv, Fig1B_DEG_heatmap.pdf, DEG_gene_clusters_h5_annot.csv, DEG_cluster_<k>_genes.txt

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "3-Hierarchical_clustering"))
data_dir <- file.path(code_dir, "data_files")
deg_files <- c(
  ART  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.ART_DMSO.DEG.csv",
  MEF  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.MEF_DMSO.DEG.csv",
  PZQ  = "../1-PZQ_MCLZ_RNA-Seq/output/pzq_vivo_adult.PZQ_DMSO.DEG.csv",
  MCLZ = "../1-PZQ_MCLZ_RNA-Seq/output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv")
h_cut <- 5          # height at which the gene dendrogram is cut into clusters
min_drugs_tested <- 2
dir.create("output", showWarnings = FALSE)

library(tidyverse)
library(pheatmap)
library(scico)

## 1. log2FC matrix, union of DEGs
res <- map(deg_files, read.csv)
degs <- map(res, ~ filter(.x, !is.na(padj), padj < 0.01, abs(log2FoldChange) > 1)$gene_id)
deg_union <- sort(unique(unlist(degs)))
message("DEGs per drug: ", paste(names(degs), lengths(degs), collapse = ", "),
        "; union: ", length(deg_union))

deg_all <- imap(res, ~ .x %>%
                  filter(gene_id %in% deg_union) %>%
                  transmute(gene_id, !!paste0(.y, "_log2FC") := log2FoldChange)) %>%
  reduce(full_join, by = "gene_id") %>%
  arrange(gene_id)
write.csv(deg_all, "output/DEG_log2FC_all_drugs.csv", row.names = FALSE, quote = FALSE)

## 2. z-scale per drug, cluster
mat <- t(as.matrix(column_to_rownames(deg_all, "gene_id")))   # rows = drugs
mat_scaled <- t(scale(t(mat)))                                 # NA ignored by scale()
mat_scaled[!is.finite(mat_scaled)] <- NA
mat_scaled <- mat_scaled[rowSums(!is.na(mat_scaled)) > 0, , drop = FALSE]
mat_scaled <- mat_scaled[, colSums(!is.na(mat_scaled)) > 0, drop = FALSE]

## order genes by variance; genes with a log2FC for only one drug have no variance and are dropped
n_tested <- colSums(!is.na(mat_scaled))
message("DEGs tested for < ", min_drugs_tested, " drugs (not shown): ", sum(n_tested < min_drugs_tested))
mat_scaled <- mat_scaled[, n_tested >= min_drugs_tested, drop = FALSE]
gene_var   <- apply(mat_scaled, 2, var, na.rm = TRUE)
mat_scaled <- mat_scaled[, names(sort(gene_var, decreasing = TRUE)), drop = FALSE]
message("genes in heat map: ", ncol(mat_scaled))

mat_clust <- mat_scaled
mat_clust[is.na(mat_clust)] <- 0
rowhc <- hclust(dist(mat_clust))        # drugs
colhc <- hclust(dist(t(mat_clust)))     # genes

gene_clusters <- cutree(colhc, h = h_cut)
gene_cluster_df <- data.frame(gene_id = names(gene_clusters), cluster = gene_clusters)

annot <- read.csv(file.path(data_dir, "WBP19_gene_annotation.csv"),
                  fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
write.csv(left_join(gene_cluster_df, annot, by = "gene_id"),
          paste0("output/DEG_gene_clusters_h", h_cut, "_annot.csv"), row.names = FALSE)
for (cl in sort(unique(gene_clusters)))
  write.table(gene_cluster_df$gene_id[gene_cluster_df$cluster == cl],
              paste0("output/DEG_cluster_", cl, "_genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

## 3. heat map
annotation_col <- data.frame(Cluster = factor(gene_clusters), row.names = names(gene_clusters))
gaps_col <- which(diff(gene_clusters[colhc$order]) != 0)
pal <- scico(200, palette = "vik")

draw <- function(file, ...) pheatmap(
  mat_scaled, cluster_rows = rowhc, cluster_cols = colhc,
  annotation_col = annotation_col, gaps_col = gaps_col,
  show_rownames = TRUE, show_colnames = FALSE, fontsize_row = 12,
  color = pal, na_col = "grey70", filename = file, ...)
draw("output/Fig1B_DEG_heatmap.pdf",  width = 12, height = 4)
