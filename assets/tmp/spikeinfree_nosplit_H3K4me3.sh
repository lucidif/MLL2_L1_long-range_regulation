#!/usr/bin/env bash

MAIN="${PWD}"
SUBANDIR="${MAIN}/outs/quantile_normalization_analysis/nosplit_H3K4me3"
INBAM="/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp/inbam"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
GENOME_SIZES="./outs/Lara_multiomic_analysis/in/mm10.sizes"
source "${DOCKER_ENV}"



# METAFILES=(
#     "${SUBANDIR}/D0_H3K27ac.tsv"
#     "${SUBANDIR}/D0_H3K27me3.tsv"
#     "${SUBANDIR}/D0_H3K4me1.tsv"
#     "${SUBANDIR}/D0_H3K4me2.tsv"
#     "${SUBANDIR}/D0_H3K4me3.tsv"
#     "${SUBANDIR}/D0_H3K9me3.tsv"
#     "${SUBANDIR}/D0_H4K16ac.tsv"
#     "${SUBANDIR}/D4_H3K27ac.tsv"
#     "${SUBANDIR}/D4_H3K27me3.tsv"
#     "${SUBANDIR}/D4_H3K4me2.tsv"
#     "${SUBANDIR}/D4_H3K4me3.tsv"
#     "${SUBANDIR}/D4_H3K9me3.tsv"
#     "${SUBANDIR}/D4_H4K16ac.tsv"
# )


#"${SUBANDIR}/NOsplit_H3K4me3.tsv"
METAFILE="${MAIN}/sheets/spikeinfree_NOSPLIT_all_samples.tsv"
BASENAME=$(basename "$METAFILE" .tsv)   # es. D0_H3K4me3
OUTDIR="${SUBANDIR}"
PREFIX="${BASENAME}_spikeinFree"

R_SCRIPT=$(cat <<EOF
library("ChIPseqSpikeInFree")
setwd("${OUTDIR}")
meta.info <- read.table("${METAFILE}", sep="\t", header=TRUE)
bams <- file.path("${INBAM}", meta.info\$ID)
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = "${METAFILE}", prefix = "${PREFIX}")
quit()
EOF
)

    echo ">>> Lanciando Docker per ChIPseqSpikeInFree..."
    echo "$R_SCRIPT" | sudo docker run "${DOCKER_ARGS[@]}" --rm -i \
        -w "${OUTDIR}" \
        -v "${SUBANDIR}:${SUBANDIR}" \
        --user "$(id -u)":"$(id -g)" \
        lucidif/chipseqspikeinfree:1.2.4 R --no-save

    SF_FILE="${OUTDIR}/${PREFIX}_SF.txt"
    if [[ ! -f "$SF_FILE" ]]; then
        echo "ERRORE: SF file non trovato: $SF_FILE — skip bigwig step"
        continue
    fi

    # --- Step 3: genera bigwig ---
    BW_OUTDIR="${OUTDIR}/spikeinfree_bw"
    echo ">>> Generando bigwigs..."
    git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
        -s "$SF_FILE" \
        -b "${INBAM}/" \
        -g "${GENOME_SIZES}" \
        -o "${BW_OUTDIR}" \
        -t 100000000

    # --- Step 4: copia tutti i bw in tmp/ per averaging ---
    mkdir -p "${BW_OUTDIR}/tmp/"
    cp "${BW_OUTDIR}"/*/*.bw "${BW_OUTDIR}/tmp/" 2>/dev/null

    echo ">>> Done: ${BASENAME}"
    echo ""


#======================================
#
#======================================

#======================================
# Step 5: media dei bigwig per pattern (replicati)
#======================================

AVG_OUTDIR="${BW_OUTDIR}/average_bw"
mkdir -p "${AVG_OUTDIR}"

echo "============================================"
echo "Averaging: ${BASENAME}"
echo "============================================"

if [[ ! -d "${BW_OUTDIR}/tmp" ]]; then
    echo "ERRORE: ${BW_OUTDIR}/tmp non trovata — skip"
    exit 1
fi

PATTERNS=$(
    ls "${BW_OUTDIR}/tmp/"*.bw 2>/dev/null \
    | xargs -n1 basename \
    | sed -E 's/\.spikeinfree\.bw$//' \
    | sed -E 's/_rep[AB]$//' \
    | sort -u
)

if [[ -z "$PATTERNS" ]]; then
    echo "WARN: nessun .bw in ${BW_OUTDIR}/tmp/ — skip"
    exit 1
fi

echo ">>> Patterns trovati:"
echo "$PATTERNS"
echo ""

while IFS= read -r PATTERN; do
    MATCHING=( "${BW_OUTDIR}/tmp/${PATTERN}"*.bw )

    if [[ ! -e "${MATCHING[0]}" ]]; then
        echo "WARN: nessun file per pattern ${PATTERN} — skip"
        continue
    fi

    OUT_BW="${AVG_OUTDIR}/${PATTERN}_average.bw"
    echo ">>> Averaging: ${PATTERN} (${#MATCHING[@]} file)"

    sudo docker run "${DOCKER_ARGS[@]}" --rm \
        -v "${SUBANDIR}:${SUBANDIR}" \
        --user "$(id -u)":"$(id -g)" \
        quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
        bigwigAverage \
            -b "${MATCHING[@]}" \
            -o "${OUT_BW}" \
            -p 8 </dev/null

    echo ">>> Salvato: ${OUT_BW}"
done <<< "$PATTERNS"

echo ">>> Done: ${BASENAME}"
echo "=== Averaging completato ==="


#=======================================
#  make heatmaps
#=======================================

main_folder="$PWD"
an_folder="${main_folder}/outs/Lara_multiomic_analysis/"
wk_folder="${main_folder}/outs/quantile_normalization_analysis/nosplit_H3K4me3"

avr_path="${wk_folder}/spikeinfree_bw/average_bw/"
newrun_path="${avr_path}"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

mkdir -p "${wk_folder}"

outname="nosplit_l1_d4mll2"
macs_peaks="${an_folder}/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"

cut -f 1-$(($(head -1 $macs_peaks | awk '{print NF}') - 1)) $macs_peaks > "${wk_folder}/unstranded.bed"

#====================================================================================
# D4/D0 heatmap - nosplit spikein-normalized tracks
#====================================================================================

# RNA-seq tracks, presi tal quali dal path originale (eccezione, non spostati)
doav_dko="${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP1.markdup.sorted.bigWig \
    ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP2.markdup.sorted.bigWig \
        ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP3.markdup.sorted.bigWig"
doav_wt="${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP1.markdup.sorted.bigWig \
    ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP2.markdup.sorted.bigWig \
        ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP1.markdup.sorted.bigWig"

samples="${an_folder}/in/ucsc/gc5Base.bw \
${main_folder}/outs/test_chipseq_dowstream/otherouts/parallel_averageBigwig_D4/D4_WT_MLL2__average.bw \
${main_folder}/outs/test_chipseq_dowstream/otherouts/parallel_averageBigwig_D4/D4_WT_RbBP5__average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw \
${avr_path}/D4_WT_H3K27ac_average.bw \
${avr_path}/D4_double-KO_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${avr_path}/D4_WT_H3K27me3_average.bw \
${avr_path}/D4_double-KO_H3K27me3_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
${avr_path}/D4_WT_H3K4me3_average.bw \
${avr_path}/D4_double-KO_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw \
${avr_path}/D4_WT_H3K4me2_average.bw \
${avr_path}/D4_double-KO_H3K4me2_average.bw
"

# 19 tracce totali: gc5Base, MLL2, RbBP5, poi 4x(H3K27ac, H3K27me3, H3K4me3, H3K4me2)
thrsholds="60 5 10 40 40 40 40 80 80 80 80 80 80 80 80 70 70 70 70"
mins="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"

sudo docker run "${DOCKER_ARGS[@]}" -w "${wk_folder}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${wk_folder}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName "${wk_folder}/${outname}_deeptools_matrix.gzip" --referencePoint "center" \
--beforeRegionStartLength 1000 --afterRegionStartLength 1500

sudo docker run "${DOCKER_ARGS[@]}" -w "${wk_folder}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${wk_folder}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile "${wk_folder}/${outname}_deeptools_matrix.gzip" \
--outFileName "${wk_folder}/${outname}_plotHeatmap.pdf" \
--outFileNameMatrix "${wk_folder}/${outname}_plotHeatmap.mat.tab" \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 12


#===============================================
#
#===============================================


#inducedFC3
mainprj="${PWD}"
prj_folder="${mainprj}/outs/Lara_multiomic_analysis"

wk_folder="${mainprj}/outs/quantile_normalization_analysis/nosplit_H3K4me3"
avr_path="${wk_folder}/spikeinfree_bw/average_bw/"

outpath="${wk_folder}/dcp_lossk4me3_near_k27ac_heatmap/"
mkdir -p "${outpath}"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

#================================
# macs_peaks resta invariato: eccezione, non è output della normalizzazione spikein
macs_peaks="${prj_folder}/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed"
wc -l $macs_peaks
#396 /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/Lara_multiomic_analysis/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed

samples="${prj_folder}/in/test_chipseq_downstream/parallel_averageBigwig_D4/D4_WT_MLL2__average.bw \
${prj_folder}/in/test_chipseq_downstream/parallel_averageBigwig_D4/D4_WT_RbBP5__average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw \
${avr_path}/D4_WT_H3K27ac_average.bw \
${avr_path}/D4_double-KO_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${avr_path}/D4_WT_H3K27me3_average.bw \
${avr_path}/D4_double-KO_H3K27me3_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
${avr_path}/D4_WT_H3K4me3_average.bw \
${avr_path}/D4_double-KO_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw \
${avr_path}/D4_WT_H3K4me2_average.bw \
${avr_path}/D4_double-KO_H3K4me2_average.bw  \
"

# 18 tracce totali: MLL2, RbBP5, poi 4x(H3K27ac, H3K27me3, H3K4me3, H3K4me2)
thrsholds="4 4 80 80 80 80 250 250 250 250 150 150 150 150 150 150 150 150"
mins="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"

outname="nosplit_dcp_gainh3k27ac_in_d4"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName "${outpath}/${outname}_deeptools_matrix.gzip" --referencePoint "center" \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile "${outpath}/${outname}_deeptools_matrix.gzip" \
--outFileName "${outpath}/${outname}_plotHeatmap.pdf" \
--outFileNameMatrix "${outpath}/${outname}_plotHeatmap.mat.tab" \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 5 --sortRegions "descend"


#===============================================================
#
#===============================================================

#inducedFC3
mainprj="${PWD}"
prj_folder="${mainprj}/outs/Lara_multiomic_analysis"

wk_folder="${mainprj}/outs/quantile_normalization_analysis/nosplit_H3K4me3"
avr_path="${wk_folder}/spikeinfree_bw/average_bw/"

outpath="${wk_folder}/D4_H3K27me3_broad_narrow/"
mkdir -p "${outpath}"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

#================================
# macs_peaks resta invariato: eccezione, non è output della normalizzazione spikein
macs_peaks="${prj_folder}/outs/D4_H3K27me3_broad_narrow/induced_broad.bed"
wc -l $macs_peaks
#627 outs/D4_H3K27me3_broad_narrow/induced_broad.bed
# macs_peaks_FC2="${prj_folder}/outs/D4_H3K27me3_broad_narrow/inducedFC2_broad.bed"
# wc -l $macs_peaks_FC2
# #449 outs/D4_H3K27me3_broad_narrow/inducedFC2_broad.bed
# macs_peaks_FC3="${prj_folder}/outs/D4_H3K27me3_broad_narrow/inducedFC3_broad.bed"
# wc -l $macs_peaks_FC3
# #303 outs/D4_H3K27me3_broad_narrow/inducedFC3_broad.bed

samples="${prj_folder}/in/test_chipseq_downstream/parallel_averageBigwig_D4/D4_WT_MLL2__average.bw \
${prj_folder}/in/test_chipseq_downstream/parallel_averageBigwig_D4/D4_WT_RbBP5__average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw \
${avr_path}/D4_WT_H3K27ac_average.bw \
${avr_path}/D4_double-KO_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${avr_path}/D4_WT_H3K27me3_average.bw \
${avr_path}/D4_double-KO_H3K27me3_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
${avr_path}/D4_WT_H3K4me3_average.bw \
${avr_path}/D4_double-KO_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw \
${avr_path}/D4_WT_H3K4me2_average.bw \
${avr_path}/D4_double-KO_H3K4me2_average.bw  \
"

# 18 tracce totali: MLL2, RbBP5, poi 4x(H3K27ac, H3K27me3, H3K4me3, H3K4me2)
thrsholds="5 5 70 70 70 70 80 80 80 80 160 160 160 160 150 150 150 150"
mins="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"

outname="nosplit_induced_genes_broad_sorted_h3k4me3"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName "${outpath}/${outname}_deeptools_matrix.gzip" --referencePoint "center" \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile "${outpath}/${outname}_deeptools_matrix.gzip" \
--outFileName "${outpath}/${outname}_plotHeatmap.pdf" \
--outFileNameMatrix "${outpath}/${outname}_plotHeatmap.mat.tab" \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 13 --sortRegions "descend"

wc -l ${prj_folder}/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed
#396 /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/Lara_multiomic_analysis/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed