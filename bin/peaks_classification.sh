#!/usr/bin/env bash
set -euo pipefail

# ============================
# Configuration (edit as needed)
# ============================

# Root folder for outputs
ROOT="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream/otherouts"

# Location of source BEDs to copy in (proximal/distal and the DOWN set). Adjust if different.
SRC_HEATMAPS="${ROOT}/../deeptools_heatmaps"

# Working directory where all steps run and outputs are written
WORK="${ROOT}/bedtools_window"

# NMI BED file (must be present in WORK or available to copy from SRC_HEATMAPS)
NMI_BED="NMIs_mm10_mESC_afterLiftover_from_mm9_to_mm10.bed"

# K4me3 DOWN BED filename (expected to exist under SRC_HEATMAPS; will be copied into WORK)
DOWN_BED="fold2_th0_05_K4me3_Double_KO_vs_F_F_DBA_DESEQ2_DOWN.bed"

# Docker image for bedtools
BEDTOOLS_IMG="quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"

# ============================
# Helper to run bedtools in Docker
# ============================
docker_bedtools() {
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v "${WORK}:/work" \
    -w /work \
    "${BEDTOOLS_IMG}" \
    bedtools "$@"
}

# ============================
# Prepare folders and inputs
# ============================
mkdir -p "${WORK}"

# NOTE: Original filenames were spelled "preaks". This script keeps that spelling
# to match your files. If your files are actually named "*peaks.bed", update the lines below.
cp -f "${SRC_HEATMAPS}/proximal_preaks.bed" "${WORK}/"
cp -f "${SRC_HEATMAPS}/distal_preaks.bed"   "${WORK}/"

# Bring the NMI file and the K4me3 DOWN file into WORK if available at the source
for src in "${NMI_BED}" "${DOWN_BED}"; do
  if [[ -f "${SRC_HEATMAPS}/${src}" && ! -f "${WORK}/${src}" ]]; then
    cp -f "${SRC_HEATMAPS}/${src}" "${WORK}/"
  fi
done

cd "${WORK}"

# Sanity checks
for f in proximal_preaks.bed distal_preaks.bed "${NMI_BED}" "${DOWN_BED}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing required file: $f"; exit 1; }
done

echo "[INFO] Counting lines in proximal_preaks.bed"
wc -l proximal_preaks.bed

# ============================
# Proximal CpG +
# ============================
docker_bedtools window -w 1000 \
  -a proximal_preaks.bed -b "${NMI_BED}" \
  > proximal_CpG_plus.bed

# Keep unique intervals by first 3 columns
Rscript -e 'x<-read.table("proximal_CpG_plus.bed")[,1:3]; x<-x[!duplicated(x),]; write.table(x,"proximal_CpG_plus_unique.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")'

echo "[INFO] Counting lines in proximal_CpG_plus_unique.bed"
wc -l proximal_CpG_plus_unique.bed

# ============================
# Proximal CpG -
# ============================
docker_bedtools subtract -A \
  -a proximal_preaks.bed -b proximal_CpG_plus_unique.bed \
  > proximal_CpG_minus.bed

echo "[INFO] Counting lines in proximal_CpG_minus.bed"
wc -l proximal_CpG_minus.bed

# Diagnostic intersection (should be empty if plus/minus are disjoint)
docker_bedtools intersect \
  -a proximal_CpG_plus_unique.bed -b proximal_CpG_minus.bed \
  > proximal_CpG_intersect_plusminus.bed

echo "[INFO] Counting lines in proximal_CpG_intersect_plusminus.bed"
wc -l proximal_CpG_intersect_plusminus.bed

# ============================
# Distal set
# ============================
echo "[INFO] Counting lines in distal_preaks.bed"
wc -l distal_preaks.bed

# Distal CpG +
docker_bedtools window -w 1000 \
  -a distal_preaks.bed -b "${NMI_BED}" \
  > distal_CpG_plus.bed

Rscript -e 'x<-read.table("distal_CpG_plus.bed")[,1:3]; x<-x[!duplicated(x),]; write.table(x,"distal_CpG_plus_unique.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")'

echo "[INFO] Counting lines in distal_CpG_plus_unique.bed"
wc -l distal_CpG_plus_unique.bed

# Distal CpG -
docker_bedtools subtract -A \
  -a distal_preaks.bed -b distal_CpG_plus_unique.bed \
  > distal_CpG_minus.bed

echo "[INFO] Counting lines in distal_CpG_minus.bed"
wc -l distal_CpG_minus.bed

# ============================================
# K4me3 DOWN subsets (integrated properly)
# ============================================

echo "[INFO] Creating K4me3 DOWN subsets using ${DOWN_BED}"

# 1) Proximal/distal peaks intersecting the DOWN set
docker_bedtools intersect -wa \
  -a proximal_preaks.bed \
  -b "${DOWN_BED}" \
  > K4me3_proximal_Double_KO_vs_F_F.bed

docker_bedtools intersect -wa \
  -a distal_preaks.bed \
  -b "${DOWN_BED}" \
  > K4me3_distal_Double_KO_vs_F_F.bed

echo "[INFO] Lines in K4me3_proximal_Double_KO_vs_F_F.bed"
wc -l K4me3_proximal_Double_KO_vs_F_F.bed
echo "[INFO] Lines in K4me3_distal_Double_KO_vs_F_F.bed"
wc -l K4me3_distal_Double_KO_vs_F_F.bed

# 2) Proximal CpG +/- intersecting the DOWN set
docker_bedtools intersect -wa \
  -a proximal_CpG_plus_unique.bed \
  -b "${DOWN_BED}" \
  > K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed

docker_bedtools intersect -wa \
  -a proximal_CpG_minus.bed \
  -b "${DOWN_BED}" \
  > K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed

echo "[INFO] Lines in K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed"
wc -l K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed
echo "[INFO] Lines in K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed"
wc -l K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed

# 3) Distal CpG +/- intersecting the DOWN set
docker_bedtools intersect -wa \
  -a distal_CpG_plus_unique.bed \
  -b "${DOWN_BED}" \
  > K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed

docker_bedtools intersect -wa \
  -a distal_CpG_minus.bed \
  -b "${DOWN_BED}" \
  > K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed

echo "[INFO] Lines in K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed"
wc -l K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed
echo "[INFO] Lines in K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
wc -l K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed

echo "[DONE] All outputs written to: ${WORK}"
