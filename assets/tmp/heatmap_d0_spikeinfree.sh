#!/usr/bin/env bash

prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
outpath="${PWD}/outs/quantile_normalization_analysis/d0_heatmaps"

${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/


mkdir -p $outpath

#inpath="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream/otherouts/deeptools_heatmaps/"
#wkdir="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream"

plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"
macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/proximal_CpG_plus_unique.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/proximal_CpG_minus.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/distal_CpG_plus_unique.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/distal_CpG_minus.bed"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

#panel A

outname="fig1b"

samples="${old1_avr_path}/Anti-GFP_average.bw \
${old1_avr_path}/Anti-Mll1_average.bw \
${old1_avr_path}/Mll2_KO_Mll1_average.bw \
${old1_avr_path}/Anti-RbBp5_average.bw \
${old1_avr_path}/Double_KO_RbBP5_average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_Mll1-KO_H3K27ac_average.bw \
${avr_path}/D0_Mll2-KO_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K27me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me2_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw"

slabels="D0_WT_MLL2 \
D0_WT_MLL1 \
D0_Mll2-KO_MLL1 \
D0_WT_RBBP5 \
D0_double-KO_RBBP5 \
D0_WT_H3K27ac \
D0_Mll1-KO_H3K27ac \
D0_Mll2-KO_H3K27ac \
D0_double-KO_H3K27ac \
D0_WT_H3K27me3 \
D0_Mll1-KO_H3K27me3 \
D0_Mll2-KO_H3K27me3 \
D0_double-KO_H3K27me3 \
D0_WT_H3K4me3 \
D0_Mll1-KO_H3K4me3 \
D0_Mll2-KO_H3K4me3 \
D0_double-KO_H3K4me3 \
D0_WT_H3K4me2 \
D0_Mll1-KO_H3K4me2 \
D0_Mll2-KO_H3K4me2 \
D0_double-KO_H3K4me2"


#thrsholds="10 10 10 10 15 15 10 10 10 10 10 10 10 10 40 40 40 40 30 30 30 30 10 10 10 10"
#thrsholds="10 10 10 15 15 10 10 10 10 10 10 10 10 40 40 40 40 30 30 30 30"
#thrsholds="10 10 10 15 15 70 70 70 70 35 35 35 35 400 400 400 400 150 150 150 150"

# fig1b e fig1c
thrsholds="10 10 10 15 15 70 70 70 70 80 80 80 80 160 160 160 160 150 150 150 150"

# fig2a
fig2a_thrsholds="10 80 80 30 30 160 160"

# figS2a — H3K27ac: 70 (già corretto), H3K4me1 invariato
s2_thrholds="10 70 70 70 70 70"

# figS2b — H3K4me3: 40 → 200, H3K4me1 invariato
s2b_thrholds="10 160 70 70 70 70"


sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile $outpath/${outname}_deeptools_matrix.gzip \
--zMax ${thrsholds} --outFileName $outpath/${outname}_plotHeatmap.pdf --sortUsingSamples 6 \
--samplesLabel ${slabels} --regionsLabel ${plabels} \
--outFileNameMatrix $outpath/${outname}_plotHeatmap.mat.tab

#====================================
#   fig S2a
#====================================

sup_outname="figS2a"

s2_samples="${old1_avr_path}/Anti-GFP_average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K4me1_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me1_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me1_average.bw \
${avr_path}/D0_double-KO_H3K4me1_average.bw"

s2_slabels="D0_WT_MLL2 \
D0_WT_H3K27ac \
D0_WT_H3K4me1 \
D0_Mll1-KO_H3K4me1 \
D0_Mll2-KO_H3K4me1 \
D0_double-KO_H3K4me1"


sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${s2_samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${sup_outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile $outpath/${sup_outname}_deeptools_matrix.gzip \
--zMax ${s2_thrholds} --outFileName $outpath/${sup_outname}_plotHeatmap.pdf --sortUsingSamples 2 \
--samplesLabel ${s2_slabels} --regionsLabel ${plabels} \
--outFileNameMatrix $outpath/${sup_outname}_plotHeatmap.mat.tab


#====================================
#       fig1c
#====================================

fig1c_outname="fig1c"

fig1c_macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"

plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

#thrsholds="10 10 10 15 15 70 70 70 70 35 35 35 35 400 400 400 400 150 150 150 150"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $fig1c_macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${fig1c_outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile $outpath/${fig1c_outname}_deeptools_matrix.gzip \
--zMax ${thrsholds} --outFileName $outpath/${fig1c_outname}_plotHeatmap.pdf --sortUsingSamples 14 \
--samplesLabel ${slabels} --regionsLabel ${plabels} \
--outFileNameMatrix $outpath/${fig1c_outname}_plotHeatmap.mat.tab

#====================================
#       fig2A
#====================================

fig2a_outname="fig2a"
plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

fig2a_macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"

fig2a_samples="${old1_avr_path}/Anti-GFP_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${prj_folder}/outs/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps/D0_WT_RING1B__average.bw \
${prj_folder}/outs/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps/D0_KO_RING1B__average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
"

#fig2a_thrsholds="10 35 35 20 20 150 150"
#fig2a_thrsholds="10 35 35 20 20 400 400"

fig2a_slabels="D0_WT_MLL2 \
D0_WT_H3K27me3 \
D0_double-KO_H3K27me3 \
D0_WT_RING1B_av \
D0_KO_RING1B_av \
D0_WT_H3K4me3 \
D0_double-KO_H3K4me3 \
"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${fig2a_samples} -R ${fig2a_macs_peaks} \
-b 1000 --sortUsingSamples 6 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${fig2a_outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile $outpath/${fig2a_outname}_deeptools_matrix.gzip \
--zMax ${fig2a_thrsholds} --outFileName $outpath/${fig2a_outname}_plotHeatmap.pdf --sortUsingSamples 6 \
--samplesLabel ${fig2a_slabels} --regionsLabel ${plabels} \
--outFileNameMatrix $outpath/${fig2a_outname}_plotHeatmap.mat.tab




#====================================
#   fig S2b
#====================================

sup_outname="figS2b"

s2b_samples="${old1_avr_path}/Anti-GFP_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me1_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me1_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me1_average.bw \
${avr_path}/D0_double-KO_H3K4me1_average.bw"

s2b_slabels="D0_WT_MLL2 \
D0_WT_H3K4me3 \
D0_WT_H3K4me1 \
D0_Mll1-KO_H3K4me1 \
D0_Mll2-KO_H3K4me1 \
D0_double-KO_H3K4me1"

#s2b_thrholds="10 400 100 100 100 100"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${s2b_samples} -R $fig1c_macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${sup_outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile $outpath/${sup_outname}_deeptools_matrix.gzip \
--zMax ${s2b_thrholds} --outFileName $outpath/${sup_outname}_plotHeatmap.pdf --sortUsingSamples 2 \
--samplesLabel ${s2b_slabels} --regionsLabel ${plabels} \
--outFileNameMatrix $outpath/${sup_outname}_plotHeatmap.mat.tab




