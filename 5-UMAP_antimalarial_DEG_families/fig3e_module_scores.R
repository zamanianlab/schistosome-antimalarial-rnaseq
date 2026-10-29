# fig3e_module_scores.R - Fig 3E gene-set module scores on the adult single-cell atlas UMAP
# in : GSE146736_adult_scseq_seurat.rds (GEO GSE146736; 3.9 GB, not in repo; set seurat_rds below)
# out: output/<set>/<set>_module_score_UMAP.pdf, _zscore_mean_by_cluster.csv, _module_score_by_cluster.csv

## settings
code_dir   <- "path/to/code"
seurat_rds <- "path/to/GSE146736_adult_scseq_seurat.rds"  # <- edit
setwd(file.path(code_dir, "5-UMAP_antimalarial_DEG_families"))
cluster_col <- "integrated_snn_res.3"

gene_sets <- list(
  ART_saposins = c("Smp_327810", "Smp_028870", "Smp_181070", "Smp_014570", "Smp_181060",
                   "Smp_105420", "Smp_105450", "Smp_137700", "Smp_194910", "Smp_137720",
                   "Smp_130100", "Smp_016490", "Smp_335370", "Smp_202610"),
  ART_cathepsins = c("Smp_210500", "Smp_141610", "Smp_019030", "Smp_105370", "Smp_085180",
                     "Smp_344940", "Smp_085010", "Smp_187140", "Smp_139240", "Smp_067060",
                     "Smp_103610", "Smp_343260"),
  ART_bax_inhibitor = c("Smp_345180", "Smp_331230", "Smp_131490", "Smp_315230", "Smp_331240",
                        "Smp_331250", "Smp_131470", "Smp_121930", "Smp_324250", "Smp_004420",
                        "Smp_317430", "Smp_336290"),
  MEF_eggshell = c("Smp_112450", "Smp_000390", "Smp_077890", "Smp_000270", "Smp_332440",
                   "Smp_000280", "Smp_000420", "Smp_326610", "Smp_326540", "Smp_340010",
                   "Smp_307900"))

library(Seurat)
library(tidyverse)

## load atlas
schistodata <- readRDS(seurat_rds)
umap    <- schistodata@reductions$umap@cell.embeddings
rna_mat <- schistodata@assays$RNA@data                  # log-normalized
common_cells <- intersect(colnames(rna_mat), rownames(umap))
umap     <- umap[common_cells, , drop = FALSE]
clusters <- schistodata@meta.data[common_cells, cluster_col]
rm(schistodata); invisible(gc())

## one gene set
module_score_set <- function(set_name, genes) {
  out <- file.path("output", set_name); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  genes_hyphen <- gsub("_", "-", genes)                  # atlas uses Smp-XXXXXX
  genes_found  <- intersect(genes_hyphen, rownames(rna_mat))
  if (length(setdiff(genes_hyphen, genes_found)))
    message(set_name, ": not in atlas: ", paste(setdiff(genes_hyphen, genes_found), collapse = ", "))

  expr_scaled <- t(scale(t(as.matrix(rna_mat[genes_found, common_cells, drop = FALSE]))))
  expr_scaled[is.na(expr_scaled)] <- 0                   # constant genes
  module_score <- colMeans(expr_scaled)

  ## UMAP, high-score cells drawn last
  mod_df <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], score = module_score)
  mod_df <- mod_df[order(mod_df$score), ]
  p <- ggplot(mod_df, aes(UMAP_1, UMAP_2, color = score)) +
    geom_point(size = 0.6, alpha = 0.9) +
    scale_color_gradientn(colours = c("lavender", "#600000")) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          legend.title = element_text(size = 12), legend.text = element_text(size = 10),
          legend.key.width = unit(0.4, "cm"), legend.key.height = unit(0.4, "cm")) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Module\nscore")
  ggsave(file.path(out, paste0(set_name, "_module_score_UMAP.pdf")), p, width = 6, height = 5)

  ## mean z-score per gene per atlas cluster
  cl_levels <- sort(unique(clusters))
  z_by_cl <- sapply(cl_levels, function(cl) rowMeans(expr_scaled[, clusters == cl, drop = FALSE]))
  colnames(z_by_cl) <- cl_levels
  write.csv(z_by_cl, file.path(out, paste0(set_name, "_zscore_mean_by_cluster.csv")), quote = FALSE)

  ## module-score summary per cluster
  data.frame(cluster = clusters, score = module_score) %>%
    group_by(cluster) %>%
    summarise(n_cells = n(), mean_score = mean(score), median_score = median(score),
              sd_score = sd(score), .groups = "drop") %>%
    arrange(desc(mean_score)) %>%
    write.csv(file.path(out, paste0(set_name, "_module_score_by_cluster.csv")), row.names = FALSE)
  message(set_name, ": ", length(genes_found), " genes scored")
}

iwalk(gene_sets, ~ module_score_set(.y, .x))
