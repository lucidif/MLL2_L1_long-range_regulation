#!/usr/bin/env bash

REF_NAME="chr6_Cxcl12"

# TSV="sheets/VirtualC_ref_coordinates.tsv"
# NF_INPUT_ALL="sheets/NFSS_Lara_all_microC.csv"

# MM10_FA="/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa"
# MM10_FAI="/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa.fai"
# MM10_FAI_SRC="/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fai"

PIPELINE_PATH="git/nf-core-microc"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"

export NXF_VER=23.10.0

# BIN=5000
# THREADS=12
# TOPN1=10
# TOPN2=20

# Main folders (keep everything under outs/virtual4C for the test)
OUTDIR_MAIN="outs/virtual4C"
CUSTOM_REF_ROOT="in/customReferences"

# ----------------------------
# Derived paths
# ----------------------------
 mainprj_folder="$(pwd)"

 mkdir -p "${OUTDIR_MAIN}" "${CUSTOM_REF_ROOT}"

# OUTBED="${OUTDIR_MAIN}/customRef/${REF_NAME}.bed"
# REFDIR="${CUSTOM_REF_ROOT}/${REF_NAME}"
# FASTA_PATH="${REFDIR}/genome.fa"

source "${DOCKER_ENV}"

bash git/nf-core-microc/bin/virtual4C_make_reference.sh \
  --ref-name chr6_Cxcl12 \
  --tsv sheets/VirtualC_ref_coordinates.tsv \
  --mm10-fa /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
  --mm10-fai-src /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fai \
  --custom-ref-root in/customReferences \
  --outdir-main outs/virtual4C \
  --docker-env git/Lara_MLL2/bin/docker_env.sh

#=========================================
#       D4 wt
#=========================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_WT_day4"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_WT_day4"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_WT_day4_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false
#  -resume

git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

#mkdir -p outs/virtual4C/bw_tracks_chr6_Cxcl12_WT_day4
#mkdir -p outs/virtual4C/bw_tracks_chr6_Cxcl12_WT_day4/1500
#mkdir -p outs/virtual4C/bw_tracks_chr6_Cxcl12_WT_day4/1500/WT_day4/

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr6_Cxcl12_WT_day4 \
  --viewpoint 117168656 \
  --window 1500 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix WT_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0


#change viewpoint

sudo git/nf-core-microc/bin/v4c_recalc_window_from_pairs_and_bam.sh \
  --ref-name chr6_Cxcl12_WT_day4 \
  --sam-prefix WT_day4 \
  --viewpoint 117168574 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --bam-mm10 outs/virtual4C/bw_tracks_chr6_Cxcl12_WT_day4/1500/WT_day4/WT_day4.pooled.viewpoint.filtered.lifted.mm10.bam \
  --outroot outs/virtual4C/bw_new_viewp \
  --smooth 0

#=========================================
#       D4 ko
#=========================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_KO_day4"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_KO_day4"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_KO_day4_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

sudo -v
git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr6_Cxcl12_KO_day4 \
  --viewpoint 117168574 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/bw_new_viewp \
  --sam-prefix KO_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

#=============================================
#     D0 WT
#=============================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_WT_day0"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_WT_day0"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_WT_day0_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

sudo -v
git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr6_Cxcl12_WT_day0 \
  --viewpoint 117168574 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/bw_new_viewp \
  --sam-prefix WT_day0 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0


#=============================================
#     D0 KO
#=============================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_KO_day0"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_KO_day0"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_KO_day0_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr6_Cxcl12" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

sudo -v
git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr6_Cxcl12_KO_day0 \
  --viewpoint 117168574 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/bw_new_viewp \
  --sam-prefix KO_day0 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0
