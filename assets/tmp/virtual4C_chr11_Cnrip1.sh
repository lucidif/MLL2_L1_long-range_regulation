#!/usr/bin/env bash

REF_NAME="chr11_Cnrip1"

PIPELINE_PATH="git/nf-core-microc"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"

export NXF_VER=23.10.0

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


# ----------------------------
# make references
# ----------------------------

bash git/nf-core-microc/bin/virtual4C_make_reference.sh \
  --ref-name chr11_Cnrip1 \
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
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr11_Cnrip1/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr11_Cnrip1" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

echo "git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}"
git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr11_Cnrip1_WT_day4 \
  --viewpoint 17052170 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix WT_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

# git/nf-core-microc/bin/v4c_recalc_window_from_pairs_and_bam.sh \
#   --ref-name chr11_Cnrip1_WT_day4 \
#   --sam-prefix WT_day4 \
#   --viewpoint 17052170 \
#   --window 1000 \
#   --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
#   --nf-out "${NF_OUT}" \
#   --bam-mm10 outs/virtual4C/bw_tracks_chr11_Cnrip1_WT_day4/1500/WT_day4/WT_day4.pooled.viewpoint.filtered.lifted.mm10.bam \
#   --outroot outs/virtual4C \
#   --smooth 0

git/nf-core-microc/bin/v4c_recalc_window_from_pairs_and_bam.sh \
  --ref-name chr11_Cnrip1_WT_day4 \
  --sam-prefix WT_day4 \
  --viewpoint 17051818 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --bam-mm10 outs/virtual4C/bw_tracks_chr11_Cnrip1_WT_day4/1500/WT_day4/WT_day4.pooled.viewpoint.filtered.lifted.mm10.bam \
  --outroot outs/virtual4C/bw_new_viewp \
  --smooth 0

#===========================================
#   D4 ko
#===========================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_KO_day4"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_KO_day4"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_KO_day4_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr11_Cnrip1/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr11_Cnrip1" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

echo "git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}"
git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

sudo -v
git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr11_Cnrip1_KO_day4 \
  --viewpoint 17051818 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/bw_new_viewp \
  --sam-prefix KO_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0
