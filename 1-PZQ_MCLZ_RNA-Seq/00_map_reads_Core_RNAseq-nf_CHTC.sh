#!/bin/bash
# =============================================================================
# Read trimming (fastp) and alignment (STAR) with the Zamanian lab Nextflow
# pipeline Core_RNAseq-nf  (https://github.com/zamanianlab/Core_RNAseq-nf)
#
# Run on the UW-Madison CHTC HTCondor pool. This is the job wrapper used for
# the paired-end S. mansoni libraries in this study (PZQ adults re-mapped from
# PRJNA597909, SAMN13691171-80; ART / MEF / control adults, PRJNA1544251).
# The recorded run for the ART/MEF cohort (Nextflow 20.07.1, 26 Aug 2025) was:
#
#   nextflow run Core_RNAseq-nf/WB-pe.nf -w work -c Core_RNAseq-nf/chtc.config \
#       --dir Schisto_S3 --star --release WBPS19 \
#       --species schistosoma_mansoni --prjn PRJEA36577 --rlen 50
#
# WB-pe.nf downloads the S. mansoni v10 genome and canonical gene set from
# WormBase ParaSite release 19, builds the STAR index, trims reads with fastp,
# aligns with STAR (--quantMode GeneCounts) and writes <sample>.ReadsPerGene.tab
# (gene_id <tab> unstranded count). Those files are in ../data_files/star_counts/
# and are the input to 01_counts_tpm_pzq_mclz.R and ../2-MEF_ART_RNA-Seq/01_counts_tpm_mef_art.R.
#
# FASTQ download: PZQ     - NCBI BioProject PRJNA597909 (BioSamples SAMN13691171-80)
#                 ART/MEF - NCBI BioProject PRJNA1544251 (SRR41195616-SRR41195630)
# =============================================================================

export HOME=$PWD
mkdir -p input work output

# copy and extract the input tar file(s) staged for this job
for sample in "$@"; do
    echo "Processing sample: $sample"
    cp -r /staging/groups/zamanian_group/input/${sample}.tar input/
    (cd input && tar -xvf ${sample}.tar && rm ${sample}.tar)
done

git clone https://github.com/zamanianlab/Core_RNAseq-nf.git

export NXF_OPTS='-Xms1g -Xmx8g'
nextflow run Core_RNAseq-nf/WB-pe.nf -w work -c Core_RNAseq-nf/chtc.config --dir input \
  --star --release "WBPS19" --species "schistosoma_mansoni" --prjn "PRJEA36577" --rlen "50"

# return outputs (counts/, star/, trim_stats/) as one tarball
rm -r work input
output_name="combined_inputs_$(date +%Y%m%d_%H%M%S)"
(cd output && tar -cvf ${output_name}.tar ./*)
STAGING_OUT="/staging/groups/zamanian_group/output"
rm -f ${STAGING_OUT}/${output_name}.tar
mv output/${output_name}.tar ${STAGING_OUT}/
