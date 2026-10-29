# schistosome-antimalarial-rnaseq
R scripts for comparative RNA-seq of Schistosoma mansoni treated with artemether and mefloquine

# A. Directories

Run folders 1 and 2 before folders 3, 4 and 6. Read mapping was done on a computing cluster, not in R. 1-PZQ_MCLZ_RNA-Seq/00_map_reads_Core_RNAseq-nf_CHTC.sh is included as a record, and its output is data_files/star_counts/.
	
1-PZQ_MCLZ_RNA-Seq

Input	data_files/star_counts/	
Output	PZQ and MCLZ counts and TPM, DESeq2 tables

2-MEF_ART_RNA-Seq

Input	data_files/star_counts/	
Output	MEF and ART counts and TPM, DESeq2 tables, PCA plot

3-Hierarchical_clustering	

Input	DESeq2 tables from folders 1–2	
Output 	Clustering shown in Fig. 1B

4-Cell_markers_DEGs	

Input	DESeq2 tables from folders 1–2, Wendt et al 2020 supplemental table 1 (abb7709_tables1)	
Output	Supp. Table 1

5-UMAP_antimalarial_DEG_families	

Input	GSE146736 Seurat object (set seurat_rds)
Output	Fig. 3E

6-InterPro_enrichment	

Input	DESeq2 tables from folders 1–2, data_files/WBP19_interpro_domains.csv,
	data_files/gene_families_MEG_PSI.csv
Output	interpro_enrichment.R - InterPro domain enrichment among the DEGs changed by both drug classes
	family_enrichment.R - the same test for the micro-exon, plexin-containing and CD59-like
	families, which carry no InterPro domain

# B. Download instructions for the publicly available datasets

ART, MEF and vehicle RNA-seq (this study)

Source	NCBI BioProject PRJNA1544251

PZQ RNA-seq

Source	NCBI BioProject PRJNA597909

MCLZ RNA-seq

Source	NCBI BioProject PRJNA1150591

Sample sheet

Files	data_files/sample_metadata.csv

Protein domain annotations

Source	WormBase ParaSite release 19 BioMart, S. mansoni PRJEA36577
Files	data_files/WBP19_interpro_domains.csv

Gene family lists

Source	Curated micro-exon and plexin-domain gene lists used for the Fig. 3 volcano plots
Files	data_files/gene_families_MEG_PSI.csv

S. mansoni gene set and gene descriptions

Source	WormBase ParaSite release 19, S. mansoni PRJEA36577
Files	data_files/Sma.geneset.csv, data_files/WBP19_gene_annotation.csv

Single-cell cluster markers

Source	Wendt et al. 2020, Science 369:1644, Table S1 (doi 10.1126/science.abb7709)

Single-cell atlas

Source	GEO GSE146736
