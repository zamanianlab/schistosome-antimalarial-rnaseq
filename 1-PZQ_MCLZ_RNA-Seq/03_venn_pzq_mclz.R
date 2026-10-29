# 03_venn_pzq_mclz.R - Fig 2A Venn diagram, PZQ and MCLZ DEGs
# in : output/pzq_vivo_adult.PZQ_DMSO.DEG.csv, mclz_vivo_adult.MCLZ_DMSO.DEG.csv
# out: output/Fig2A_venn_PZQ_MCLZ.pdf

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "1-PZQ_MCLZ_RNA-Seq"))

library(tidyverse)
library(ggVennDiagram)
library(sf)

pzq  <- read.csv("output/pzq_vivo_adult.PZQ_DMSO.DEG.csv")
mclz <- read.csv("output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv")

up   <- function(d) d %>% filter(!is.na(padj), padj < 0.01, log2FoldChange >=  1) %>% pull(gene_id)
down <- function(d) d %>% filter(!is.na(padj), padj < 0.01, log2FoldChange <= -1) %>% pull(gene_id)

venn_list <- list(PZQ_up = up(pzq), MCLZ_up = up(mclz), PZQ_down = down(pzq), MCLZ_down = down(mclz))

venn.plot <- ggVennDiagram(venn_list, label = "both", label_alpha = 0,
                           set_size = 5, edge_size = 0.6) +
  scale_fill_gradient(low = "white", high = "grey60") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +   # room for set labels
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0, "pt"))

ggsave("output/Fig2A_venn_PZQ_MCLZ.pdf",  venn.plot, width = 5, height = 5, units = "in")
