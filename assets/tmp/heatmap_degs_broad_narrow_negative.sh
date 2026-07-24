

#=====================================
#same scale
#=====================================

%%script bash

#all track

#ref_regions="outs/downGenesHeatmap/D0dkoVSwt.downgenes.bed"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"


DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

ref_regions="${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed \
${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed \
${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed"

echo "broad.downgenes n:"
wc -l ${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed

echo "narrow.downgenes n:"
wc -l ${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed

echo "negative.downgenes n:"
wc -l ${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed


outname="fig_1I"

outpath="outs/Lara_multiomic_analysis/"
plabels="broad narrow negative"


#thrsholds="10 150 150 20 20 200 200"
thrsholds="10 150 150 20 20 160 160"

samples="${old1_avr_path}/Anti-GFP_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/D0_WT_RING1B__average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/D0_KO_RING1B__average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
"

slabels="Anti-GFP_av \
WT_K27me3_av \
KO_K27me3_av \
WT_RING1B_av \
KO_RING1B_av \
WT_K4me3_av \
KO_K4me3_av \
"

sudo docker run "${DOCKER_ARGS[@]}" -w "${prj_folder}/${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -R $ref_regions -S $samples --sortUsingSamples 2 --numberOfProcessors 6 --scale 1 --binSize 10 --averageTypeBins "median" --outFileName ${prj_folder}/outs/Lara_multiomic_analysis/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -w "${prj_folder}/${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
plotHeatmap --colorMap 'seismic' --missingDataColor 0.6 --matrixFile outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip --zMax ${thrsholds} --outFileName outs/downGenesHeatmap/${outname}_sorted_by_k4me3_plotHeatmap.pdf --sortUsingSamples 6 --regionsLabel ${plabels} --samplesLabel ${slabels} --outFileNameMatrix outs/downGenesHeatmap/${outname}_plotHeatmap.mat.tab


