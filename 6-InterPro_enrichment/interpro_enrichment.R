# interpro_enrichment.R - InterPro domain enrichment in the DEGs shared by both drug classes
# in : DESeq2 tables from folders 1-2, ../data_files/WBP19_interpro_domains.csv
# out: output/InterPro_enrichment_overlap_both_classes.csv

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "6-InterPro_enrichment"))

padj_cut   <- 0.01  # DESeq2 adjusted P for calling a DEG
lfc_cut    <- 1     # |log2 fold change| for calling a DEG
min_family <- 5     # smallest InterPro family that is tested
min_hits   <- 2     # smallest number of family members in the gene set

library(tidyverse)

## InterPro annotation, one row per gene-domain pair
anno <- read.csv("../data_files/WBP19_interpro_domains.csv", check.names = FALSE) %>%
  select(gene = `Gene stable ID`, term = `InterPro ID`,
         description = `InterPro short description`) %>%
  filter(!is.na(term), term != "") %>%
  distinct(gene, term, .keep_all = TRUE)

## DESeq2 tables, one per treatment
deg_files <- c(
  PZQ  = "../1-PZQ_MCLZ_RNA-Seq/output/pzq_vivo_adult.PZQ_DMSO.DEG.csv",
  MCLZ = "../1-PZQ_MCLZ_RNA-Seq/output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv",
  MEF  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.MEF_DMSO.DEG.csv",
  ART  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.ART_DMSO.DEG.csv")

deg <- imap_dfr(deg_files, ~ read.csv(.x) %>%
                  transmute(gene = gene_id, drug = .y,
                            lfc = log2FoldChange, padj = padj))

## background: genes with a usable result in all four treatments
background <- deg %>% filter(!is.na(padj)) %>% count(gene) %>%
  filter(n == length(deg_files)) %>% pull(gene)

## the test universe is the annotated part of that background
universe <- intersect(background, unique(anno$gene))
message(sprintf("background %d genes, annotated universe %d genes",
                length(background), length(universe)))

## gene sets: changed in the same direction by at least one drug of each class
deg_in <- function(drugs, direction) {
  deg %>% filter(drug %in% drugs, gene %in% background, !is.na(padj), padj < padj_cut,
                 if (direction == "up") lfc >= lfc_cut else lfc <= -lfc_cut) %>%
    pull(gene) %>% unique()
}
overlap_set <- function(direction)
  intersect(deg_in(c("PZQ", "MCLZ"), direction), deg_in(c("MEF", "ART"), direction))

set_up   <- overlap_set("up")
set_down <- overlap_set("down")
message(sprintf("shared by both classes: %d up, %d down", length(set_up), length(set_down)))

## one-sided Fisher exact test per InterPro term, Benjamini-Hochberg adjusted
interpro_enrichment <- function(genes, universe, annotation,
                                min_family = 5, min_hits = 2) {
  universe <- unique(universe)
  genes    <- intersect(unique(genes), universe)
  n_all    <- length(universe)
  n_genes  <- length(genes)

  annotation %>%
    filter(gene %in% universe) %>%
    group_by(term, description) %>%
    summarise(n_universe = n_distinct(gene),
              n_set      = n_distinct(gene[gene %in% genes]), .groups = "drop") %>%
    filter(n_universe >= min_family, n_set >= min_hits) %>%
    mutate(fold = (n_set / n_genes) / (n_universe / n_all),
           p = map2_dbl(n_set, n_universe, ~ fisher.test(
             matrix(c(.x, .y - .x, n_genes - .x, n_all - .y - n_genes + .x), nrow = 2),
             alternative = "greater")$p.value),
           padj = p.adjust(p, method = "BH")) %>%
    arrange(padj, p) %>%
    select(term, description, n_universe, n_set, fold, p, padj)
}

## run both directions and write the table
enrichment <- bind_rows(
  interpro_enrichment(set_up,   universe, anno, min_family, min_hits) %>%
    mutate(direction = "up", .before = 1),
  interpro_enrichment(set_down, universe, anno, min_family, min_hits) %>%
    mutate(direction = "down", .before = 1))

write.csv(enrichment, "output/InterPro_enrichment_overlap_both_classes.csv", row.names = FALSE)
message(sprintf("%d terms tested, %d with padj < 0.05",
                nrow(enrichment), sum(enrichment$padj < 0.05)))
