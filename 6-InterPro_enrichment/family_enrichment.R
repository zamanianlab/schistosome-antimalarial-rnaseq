# family_enrichment.R - enrichment of three gene families that carry no InterPro domain
# in : DESeq2 tables from folders 1-2, ../data_files/gene_families_MEG_PSI.csv,
#      ../data_files/WBP19_interpro_domains.csv
# out: output/Family_enrichment_MEG_PSI_CD59.csv
#
# The micro-exon (MEG) and plexin-containing (PSI) families are invisible to
# interpro_enrichment.R because no member carries an InterPro domain. They are tested
# here as named gene sets against the full background of genes tested in all four
# comparisons, rather than against the annotated universe used for the domain analysis.

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "6-InterPro_enrichment"))

padj_cut <- 0.01  # DESeq2 adjusted P for calling a DEG
lfc_cut  <- 1     # |log2 fold change| for calling a DEG

library(tidyverse)

## DESeq2 tables, one per treatment
deg_files <- c(
  PZQ  = "../1-PZQ_MCLZ_RNA-Seq/output/pzq_vivo_adult.PZQ_DMSO.DEG.csv",
  MCLZ = "../1-PZQ_MCLZ_RNA-Seq/output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv",
  MEF  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.MEF_DMSO.DEG.csv",
  ART  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.ART_DMSO.DEG.csv")

deg <- imap_dfr(deg_files, ~ read.csv(.x) %>%
                  transmute(gene = gene_id, drug = .y,
                            lfc = log2FoldChange, padj = padj))

background <- deg %>% filter(!is.na(padj)) %>% count(gene) %>%
  filter(n == length(deg_files)) %>% pull(gene)

## families. MEG and PSI are the curated lists used for the Figure 3 volcano plots, which
## ship as an explicit file because neither family is recoverable from gene descriptions
## or domain annotation. CD59-like is taken from the InterPro table and is included as the
## comparator that the domain analysis can also see.
fam <- read.csv("../data_files/gene_families_MEG_PSI.csv")
families <- list(
  "Micro-exon genes (MEG)"       = fam$gene_id[fam$family == "MEG"],
  "Plexin-containing (PSI)"      = fam$gene_id[fam$family == "PSI"],
  "CD59-like (snake toxin-like)" = read.csv("../data_files/WBP19_interpro_domains.csv",
      check.names = FALSE) %>% filter(`InterPro ID` == "IPR045860") %>%
      pull(`Gene stable ID`) %>% unique())

## gene sets: changed by at least one drug of each class, by direction and pooled
deg_in <- function(drugs, mode) {
  deg %>% filter(drug %in% drugs, gene %in% background, !is.na(padj), padj < padj_cut,
                 switch(mode, up = lfc >= lfc_cut, down = lfc <= -lfc_cut,
                        any = abs(lfc) >= lfc_cut)) %>%
    pull(gene) %>% unique()
}
gene_sets <- set_names(map(c("up", "down", "any"),
  ~ intersect(deg_in(c("PZQ", "MCLZ"), .x), deg_in(c("MEF", "ART"), .x))),
  c("up", "down", "any"))

## one-sided Fisher exact test per family
family_enrichment <- function(genes, universe, family) {
  genes  <- intersect(unique(genes), universe)
  family <- intersect(unique(family), universe)
  k <- length(intersect(family, genes))
  n_all <- length(universe); n_genes <- length(genes); n_fam <- length(family)
  tibble(n_set = k, n_family = n_fam,
         fold = if (k > 0) (k / n_genes) / (n_fam / n_all) else 0,
         p = fisher.test(matrix(c(k, n_fam - k, n_genes - k,
                                  n_all - n_fam - n_genes + k), nrow = 2),
                         alternative = "greater")$p.value,
         genes = paste(sort(intersect(family, genes)), collapse = " "))
}

enrichment <- imap_dfr(gene_sets, function(gs, set_name)
  imap_dfr(families, ~ family_enrichment(gs, background, .x) %>%
             mutate(gene_set = set_name, family = .y, .before = 1))) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  relocate(genes, .after = padj)

write.csv(enrichment, "output/Family_enrichment_MEG_PSI_CD59.csv", row.names = FALSE)
print(as.data.frame(enrichment %>% select(-genes)))
