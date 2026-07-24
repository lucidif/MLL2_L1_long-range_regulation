#!/usr/bin/env bash

inpath="${PWD}/outs/Lara_multiomic_analysis"

prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
outpath="${PWD}/outs/Lara_multiomic_analysis/outs/gc_content_heatmap/"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

# macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed \
# ${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed \
# ${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed \
# ${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"

#plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"



#=============K4me3================= two groups

title="K4me3"
outname="2grps_K4me3_DOWN_Double_KO_vs_F_F_${title}_only"
macs_peaks="${PWD}/outs/Lara_multiomic_analysis/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"
#plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

samples="${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw"

slabels="WT DKO"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#ffb00e" --plotTitle "${title} profile l1" --samplesLabel ${slabels} 




#=============K9 =================

title="H3K9me3"
outname="K4me3_DOWN_Double_KO_vs_F_F_${title}_only"
#macs_peaks="outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"
#plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

samples="${avr_path}/D0_WT_H3K9me3_average.bw \
${avr_path}/D0_double-KO_H3K9me3_average.bw"


slabels="WT DKO"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#ffb00e" --plotTitle "${title} profile l1" --samplesLabel ${slabels} 


#=========== h3k27ac

# title="H3K27ac"
# outname="K4me3_DOWN_Double_KO_vs_F_F_${title}_only"

# #macs_peaks="outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"
# #plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

# #samples="./in/test_chipseq_downstream/parallel_averageBigwig/D0_WT_H3K27ac__average.bw ./in/test_chipseq_downstream/parallel_averageBigwig/D0_DoubleKO_H3K27ac__average.bw"

# samples="${avr_path}/D0_WT_H3K27ac_average.bw \
# ${avr_path}/D0_double-KO_H3K27ac_average.bw"


# slabels="WT DKO"

# sudo docker run \
# "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $macs_peaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run \
# "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#ffb00e" --plotTitle "${title} profile l1" --samplesLabel ${slabels} 

