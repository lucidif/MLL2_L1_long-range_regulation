# ============================================================
# virtual4C_chr10_Smpdl3a.sh
#
# Master script for virtual 4C analysis at the Smpdl3a locus
# (chr10, viewpoint mm10:57,794,374).
#
# The script runs the full virtual 4C pipeline for all four
# conditions: WT day 0, WT day 4, double-KO day 0, double-KO day 4.
#
# For each condition, the following steps are executed:
#   1. Build a locus-specific reference genome for chr10_Smpdl3a
#      (extracted from mm10, indexed with bwa and samtools).
#   2. Align Micro-C reads to the custom reference using the
#      nf-core-microc Nextflow pipeline (Dovetail protocol,
#      Nextflow v23.10.0).
#   3. Generate per-sample virtual 4C profiles and a pooled
#      library-size-normalized BigWig track (--pairs-counts-tsv
#      total_nodups.tsv; scaleFactor = 1e6 / sum(total_nodups)).
#   4. Count pairs in the viewpoint window (±1,000 bp around
#      mm10:57,794,374) for viewpoint-local normalization.
#   5. Generate a second pooled BigWig track with viewpoint-local
#      normalization (--pairs-counts-tsv total_nodups.window_
#      vp57794374_w1000.tsv; scaleFactor = 1e6 / sum(pairs_in_
#      window)); output written to outs/virtual4C/vpnorm/.
#
# The viewpoint-locally normalized tracks (step 5) are used for
# WT vs double-KO comparisons, as they correct for potential
# differences in the fraction of reads originating from the
# viewpoint locus between conditions.
# ============================================================

REF_NAME="chr10_Smpdl3a"

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
  --ref-name chr10_Smpdl3a \
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
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

echo "git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}"
git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix WT_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

bash git/nf-core-microc/bin/virtual4c_reads_counting.sh \
  ${NF_OUT}/pairix \
  git/Lara_MLL2/bin/docker_env.sh \
  57794374 \
  1000

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/vpnorm \
  --sam-prefix WT_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.window_vp57794374_w1000.tsv \
  --smoothLength 0


#============================================
#   D4 ko
#============================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_KO_day4"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_KO_day4"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_KO_day4_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix KO_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

#recalculate counts on window only

bash git/nf-core-microc/bin/virtual4c_reads_counting.sh \
  ${NF_OUT}/pairix \
  git/Lara_MLL2/bin/docker_env.sh \
  57794374 \
  1000


git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C/vpnorm \
  --sam-prefix KO_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.window_vp57794374_w1000.tsv \
  --smoothLength 0

#==========================================
# D0 WT
#==========================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_WT_day0"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_WT_day0"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_WT_day0_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

echo "git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}"
git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix WT_day0 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

#==========================================
# D0 KO
#==========================================

NF_OUT="${OUTDIR_MAIN}/nfout/${REF_NAME}_KO_day0"
NF_WORK="${OUTDIR_MAIN}/work_${REF_NAME}_KO_day0"

sudo nextflow run "${PIPELINE_PATH}" -profile docker \
  -work-dir "${NF_WORK}" \
  --input "sheets/NFSS_Lara_KO_day0_microC.csv" \
  --fasta "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a/genome.fa" \
  --index "/home/lucio/wkdir/projects/MLL2_L1_regulation/in/customReferences/chr10_Smpdl3a" \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false

echo "git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}"
git/nf-core-microc/bin/counts_pairs.sh ${NF_OUT}/pairix ${DOCKER_ENV}

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix KO_day0 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv ${NF_OUT}/pairix/total_nodups.tsv \
  --smoothLength 0

