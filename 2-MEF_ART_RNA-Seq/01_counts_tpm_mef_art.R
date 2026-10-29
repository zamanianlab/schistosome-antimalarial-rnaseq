# 01_counts_tpm_mef_art.R - MEF/ART cohort count and TPM tables from STAR ReadsPerGene files
# in : ../data_files/star_counts/, Sma.geneset.csv, sample_metadata.csv
# out: output/mef_art_counts_raw.rds, mef_art_counts_tpm.rds (Supp File 1, sheets 5-6)

## settings
code_dir <- "path/to/code"
setwd(file.path(code_dir, "2-MEF_ART_RNA-Seq"))
data_dir <- file.path(code_dir, "data_files")
cohort   <- "malarial_vivo_adult"
dir.create("output", showWarnings = FALSE)

library(tidyverse)

## Sma.geneset.csv was made from the WBPS19 canonical GTF with:
## awk '$3 ~ /^gene/' schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf |
##   awk '{print $1 "," $4 "," $5 "," $7 "," $10}' | sed 's/\"//g; s/;//g' > Sma.geneset.csv
gtf <- read.csv(file.path(data_dir, "Sma.geneset.csv"), header = FALSE,
                col.names = c("chr", "start", "end", "strand", "gene_id")) %>%
  mutate(gene_len = abs(end - start)) %>%
  select(gene_id, gene_len)

samples <- read.csv(file.path(data_dir, "sample_metadata.csv")) %>%
  filter(project_batch == cohort)

## one STAR ReadsPerGene file: drop the 4 N_* summary lines, keep the unstranded column
read_star <- function(s) {
  read.delim(file.path(data_dir, "star_counts", paste0(s, ".ReadsPerGene.tab")),
             header = FALSE, col.names = c("gene_id", "count")) %>%
    slice(-(1:4)) %>%
    left_join(gtf, by = "gene_id") %>%
    mutate(rpk = count / (1000 * gene_len),          # constant factor cancels in TPM
           tpm = 1e6 * rpk / sum(rpk, na.rm = TRUE),
           sample = s) %>%
    select(gene_id, sample, count, tpm)
}
long <- map_dfr(samples$sample, read_star)

counts.raw <- long %>% transmute(gene_id, sample = factor(sample), expression = count)
counts.tpm <- long %>% transmute(gene_id, sample = factor(sample), expression = tpm)
saveRDS(counts.raw, "output/mef_art_counts_raw.rds")
saveRDS(counts.tpm, "output/mef_art_counts_tpm.rds")
