#!/usr/bin/env bash

main_folder="$PWD"
an_folder="${main_folder}/outs/Lara_multiomic_analysis/outs/"
wk_folder=${an_folder}

#avr_path="${main_folder}/outs/quantile_normalization_analysis/average_bw"
#old1_avr_path="${main_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"

avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"

newrun_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

mkdir -p "${wk_folder}"

#sudo docker run -v `pwd`:`pwd` -w `pwd` abralab/wigtobigwig:v4 wigToBigWig in/ucsc/UCSC_Main_on_Mouse_gc5BaseBw.wig in/ucsc/chrNameLength.txt in/ucsc/UCSC_Main_on_Mouse_gc5BaseBw.bw

outname="TE_heatmap"
macs_peaks="${wk_folder}/LINE.bed ${wk_folder}/LTR.bed ${wk_folder}/SINE.bed"

#cut -f 1-$(($(head -1 $macs_peaks | awk '{print NF}') - 1)) $macs_peaks > unstranded.bed

samples="/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/Lara_multiomic_analysis/in/ucsc/gc5Base.bw \
    ${old1_avr_path}/Anti-GFP_average.bw \
    ${avr_path}/D0_WT_H3K4me3_average.bw \
    ${avr_path}/D0_double-KO_H3K4me3_average.bw \
    ${avr_path}/D0_WT_H3K27ac_average.bw \
    ${avr_path}/D0_double-KO_H3K27ac_average.bw \
    ${avr_path}/D0_WT_H3K27me3_average.bw \
    ${avr_path}/D0_double-KO_H3K27me3_average.bw \
    ${avr_path}/D0_WT_H3K9me3_average.bw \
    ${avr_path}/D0_double-KO_H3K9me3_average.bw \
    "

#k3k9 as d0 cpg 

thrsholds="60 10 80 80 50 50 80 80 30 30"

#thrsholds="60 10 40 40 10 10 10 10 6 6"
#thrsholds="60 10 300 300 40 40 30 30 30 30"

mins="40 0 0 0 0 0 0 0 0 0"

sudo docker run "${DOCKER_ARGS[@]}" -w "${wk_folder}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${wk_folder}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName ${outname}_deeptools_matrix.gzip --referencePoint "center" \
--beforeRegionStartLength 1000 --afterRegionStartLength 1500

sudo docker run "${DOCKER_ARGS[@]}" -w "${wk_folder}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${wk_folder}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile ${outname}_deeptools_matrix.gzip \
--outFileName ${outname}_plotHeatmap.pdf \
--outFileNameMatrix ${outname}_plotHeatmap.mat.tab \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 3

