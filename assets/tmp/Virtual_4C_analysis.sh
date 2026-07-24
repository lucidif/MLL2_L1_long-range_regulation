#!/usr/bin/env bash
#set -euo pipefail

#==========================test with generalized scripts
REF_NAME="chr10_Rfx6"

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

#NF_INPUT="sheets/NFSS_Lara_WTd0_microC.csv"

source "${DOCKER_ENV}"

# # Where we keep test cool/hic/bedpe artifacts
# TEST_OUT="${OUTDIR_MAIN}/test_virtual4C_${REF_NAME}_WT_day4"
# COOL_DIR="${TEST_OUT}/custom_cool"
# ONEFILE_DIR="${TEST_OUT}/_onefile"
# LIFTED_DIR="${TEST_OUT}/_lifted_mm10"
# REBUILT_DIR="${TEST_OUT}/_rebuilt_mm10"
# HIC_DIR="${TEST_OUT}/_hic_mm10"
# TMP_DIR="${TEST_OUT}/_tmp_rebuild"
# BEDPE_FROM_COOL_DIR="${TEST_OUT}/bedpe_from_cool"



bash git/nf-core-microc/bin/virtual4C_make_reference.sh \
  --ref-name chr10_Rfx6 \
  --tsv sheets/VirtualC_ref_coordinates.tsv \
  --mm10-fa /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
  --mm10-fai-src /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fai \
  --custom-ref-root in/customReferences \
  --outdir-main outs/virtual4C \
  --docker-env git/Lara_MLL2/bin/docker_env.sh

#==============================
# take only specific samples
#==============================

# bash git/Lara_MLL2/bin/virtual4C_filter_samplesheet.sh \
#   --in  sheets/NFSS_Lara_all_microC.csv \
#   --out sheets/NFSS_Lara_WT_day0_microC.csv \
#   --match 'WT' \
#   --match 'day0|d0' \
#   --dry-run


# ----------------------------
# 5) Run nf-core/microc (single reference, WT day4 only)
# ----------------------------

# echo "[INFO] Running nf-core/microc for ${REF_NAME} (WT day0 only)"
# mkdir -p "${NF_OUT}" "${NF_WORK}"

#=============== 

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_WT_day4"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_WT_day4"


sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_WT_day0_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false \
  --removeUnmapped true \
  -resume

#==============================
# take only WT day4
#==============================

NF_OUT="outs/virtual4C/nfout/chr10_Rfx6_WT_day4"

bash git/Lara_MLL2/bin/virtual4C_filter_samplesheet.sh \
  --in  sheets/NFSS_Lara_all_microC.csv \
  --out sheets/NFSS_Lara_WT_day4_microC.csv \
  --match 'WT' \
  --match 'day4|d4' \
  --dry-run

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_WT_day4_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false \
  -resume

#chmod +x git/nf-core-microc/bin/virtual4C_rebuild_stdref.sh

# git/nf-core-microc/bin/virtual4C_rebuild_stdref.sh \
#     --nf-out "${NF_OUT}" \
#     --bin 5000 \
#     --outdir outs/virtual4C/96kb_virtual4C_chr10_Rfx6_WT_day4 \
#     --docker-env git/Lara_MLL2/bin/docker_env.sh \
#     --lift-script git/nf-core-microc/bin/lift_and_export_hic.sh \
#     --chrom chr10 --offset 1234567 \
#     --assembly mm10 \
#     --threads 12

chmod +x git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
    --ref-name chr10_Rfx6 \
    --viewpoint 51677810 \
    --window 1500 \
    --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
    --nf-out "${NF_OUT}" \
    --outroot outs/virtual4C/96kb \
    --sam-prefix WT_day4 \
    --docker-env git/Lara_MLL2/bin/docker_env.sh \
    --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh









# # ----------------------------
# # 5) Run nf-core/microc (single reference, WT day4 only)
# # ----------------------------
# echo "[INFO] Running nf-core/microc for ${REF_NAME} (WT day0 only)"
# mkdir -p "${NF_OUT}" "${NF_WORK}"

# sudo nextflow run "${PIPELINE_PATH}" -profile docker \
#   -work-dir "${NF_WORK}" \
#   --input "sheets/NFSS_Lara_WT_day0_microC.csv" \
#   --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6/genome.fa" \
#   --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Rfx6" \
#   --outdir "${NF_OUT}" \
#   --exclude_lcextrap false \
#   --removeUnmapped true \
#   -resume













