%%script bash

#inducedFC3

mainprj="${PWD}"
prj_folder="${mainprj}/outs/Lara_multiomic_analysis"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"

outpath="${prj_folder}/outs/D4_H3K27me3_broad_narrow/"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"


#================================

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

thrsholds="10 15 70 70 70 70 80 80 80 80 160 160 160 160 150 150 150 150"

#mins="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
mins="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"

outname="induced_genes_broad_sorted_h3k4me3"
sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName ${outname}_deeptools_matrix.gzip --referencePoint "center" \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --missingDataAsZero

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile ${outname}_deeptools_matrix.gzip \
--outFileName ${outname}_plotHeatmap.pdf \
--outFileNameMatrix ${outname}_plotHeatmap.mat.tab \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 13 --sortRegions "descend"



wc -l ${prj_folder}/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.be
d
#396 /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/Lara_multiomic_analysis/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed

