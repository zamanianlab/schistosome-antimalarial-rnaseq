# supp_table1_cell_markers.R - Supp Table 1, % of single-cell cluster markers (Wendt et al. 2020) that are DEGs
# in : ../data_files/abb7709_tables1.xlsx (Wendt et al. 2020 Table S1), DESeq2 tables from folders 1-2
# out: output/Supp_Table1_markers_DEG.csv, Supp_Table1_marker_gene_status.csv,
#      Supplemental_Table_1_heatmap_data.csv, Supplemental_Table_1_heatmap_magma.pdf

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "4-Cell_markers_DEGs"))
data_dir <- file.path(code_dir, "data_files")
marker_file <- "abb7709_tables1.xlsx"
if (!file.exists(file.path(data_dir, marker_file)))
  stop("Download Table S1 from doi 10.1126/science.abb7709 to data_files/", marker_file)
deg_files <- c(
  PZQ  = "../1-PZQ_MCLZ_RNA-Seq/output/pzq_vivo_adult.PZQ_DMSO.DEG.csv",
  MCLZ = "../1-PZQ_MCLZ_RNA-Seq/output/mclz_vivo_adult.MCLZ_DMSO.DEG.csv",
  MEF  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.MEF_DMSO.DEG.csv",
  ART  = "../2-MEF_ART_RNA-Seq/output/malarial_vivo_adult.ART_DMSO.DEG.csv")
drugs <- names(deg_files)
dir.create("output", showWarnings = FALSE)

library(tidyverse)
library(viridis)
library(readxl)

## 1. markers and DEG status
markers <- read_excel(file.path(data_dir, marker_file), skip = 1) %>%    # row 1 is the table title
  rename(gene_id = gene) %>%
  filter(p_val_adj < 0.01) %>%
  distinct(gene_id, cluster, .keep_all = TRUE)

status <- imap_dfr(deg_files, ~ read.csv(.x) %>%
  transmute(gene_id, drug = .y,
            status = case_when(!is.na(padj) & padj < 0.01 & log2FoldChange >  1 ~ "up",
                               !is.na(padj) & padj < 0.01 & log2FoldChange < -1 ~ "down",
                               TRUE ~ "not DEG")))

gene_status <- markers %>%
  select(gene_id, cluster, avg_logFC, pct.1, pct.2, p_val_adj) %>%
  crossing(drug = drugs) %>%
  left_join(status, by = c("gene_id", "drug")) %>%
  mutate(status = replace_na(status, "not tested"))
write.csv(pivot_wider(gene_status, names_from = drug, values_from = status, names_glue = "{drug}_status"),
          "output/Supp_Table1_marker_gene_status.csv", row.names = FALSE)

tab <- gene_status %>%
  group_by(cluster, drug) %>%
  summarise(markers = n(), tested = sum(status != "not tested"),
            up = sum(status == "up"), down = sum(status == "down"), .groups = "drop") %>%
  mutate(DEG_n = up + down, DEG_pct = round(100 * DEG_n / markers, 3))
write.csv(pivot_wider(tab, names_from = drug, values_from = c(tested, up, down, DEG_n, DEG_pct),
                      names_glue = "{drug}_{.value}"),
          "output/Supp_Table1_markers_DEG.csv", row.names = FALSE)

hm <- tab %>%
  filter(cluster != "unknown") %>%
  select(cluster, markers, drug, DEG_pct) %>%
  pivot_wider(names_from = drug, values_from = DEG_pct) %>%
  select(cluster, markers, all_of(drugs)) %>%
  arrange(tolower(cluster))         # case-insensitive alphabetical, as published
write.csv(hm, "output/Supplemental_Table_1_heatmap_data.csv", row.names = FALSE, quote = FALSE)

## 2. heat map
if (!exists("theme_zlab_white")) theme_zlab_white <- function(base_size = 12, base_family = "Helvetica") { }
font <- "Helvetica"
magma_direction <- -1
out_stub <- "output/Supplemental_Table_1_heatmap_magma"
n_rows <- nrow(hm)
cluster_levels <- rev(hm$cluster)
pal <- scales::gradient_n_pal(viridis::viridis_pal(option = "magma", direction = magma_direction)(6))

srgb_lin <- function(u) ifelse(u <= 0.03928, u / 12.92, ((u + 0.055) / 1.055)^2.4)
rel_lum  <- function(hex) { rgb <- col2rgb(hex) / 255
  as.numeric(0.2126 * srgb_lin(rgb[1, ]) + 0.7152 * srgb_lin(rgb[2, ]) + 0.0722 * srgb_lin(rgb[3, ])) }
text_for <- function(pct) { L <- rel_lum(pal(pct / 100))
  if_else((1.05 / (L + 0.05)) >= ((L + 0.05) / 0.05), "white", "black") }

long <- hm %>%
  pivot_longer(all_of(drugs), names_to = "drug", values_to = "pct") %>%
  mutate(drug = factor(drug, levels = drugs),
         cluster = factor(cluster, levels = cluster_levels),
         label = floor(pct + 0.5),                      # round half up
         txt_col = text_for(pct))
n_col <- hm %>% mutate(cluster = factor(cluster, levels = cluster_levels))

p <- ggplot() +
  geom_tile(data = long, aes(drug, cluster, fill = pct), colour = "grey55", linewidth = 0.15) +
  geom_text(data = long, aes(drug, cluster, label = label, colour = txt_col), size = 2.5, family = font) +
  geom_text(data = n_col, aes(x = 0.42, y = cluster, label = markers),
            hjust = 1, size = 2.5, colour = "grey30", family = font) +
  annotate("text", x = 0.42, y = n_rows + 1.1, label = "Markers (n)",
           hjust = 1, size = 2.8, fontface = "bold", colour = "grey30", family = font) +
  annotate("text", x = seq_along(drugs), y = n_rows + 1.1, label = drugs,
           size = 2.8, fontface = "bold", family = font) +
  scale_colour_identity() +
  scale_fill_viridis(option = "magma", direction = magma_direction, limits = c(0, 100),
                     breaks = seq(0, 100, 25), name = "% of cluster markers that are DEGs") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, length(drugs) + 0.5), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_family = font) + theme_zlab_white() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 0, size = 7, colour = "black"),
        axis.ticks = element_blank(), plot.margin = margin(6, 6, 6, 6),
        legend.position = "bottom", legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = unit(1.8, "in"), barheight = unit(0.08, "in"),
                               ticks.colour = "white", frame.colour = "grey55", frame.linewidth = 0.2))

ggsave(paste0(out_stub, ".pdf"),  p, width = 5.4, height = 10.4, units = "in", bg = "white")
