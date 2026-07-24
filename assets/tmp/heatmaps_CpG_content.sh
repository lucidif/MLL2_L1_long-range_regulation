#!/usr/bin/env bash


#===============================================================================
#main heatmap
#===============================================================================

main_folder="$PWD"
an_folder="${main_folder}/outs/Lara_multiomic_analysis/"
wk_folder="${an_folder}/outs/gc_content_heatmap/"

#avr_path="${main_folder}/outs/quantile_normalization_analysis/average_bw"
#old1_avr_path="${main_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"

avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"

newrun_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

mkdir -p "${wk_folder}"

#sudo docker run -v `pwd`:`pwd` -w `pwd` abralab/wigtobigwig:v4 wigToBigWig in/ucsc/UCSC_Main_on_Mouse_gc5BaseBw.wig in/ucsc/chrNameLength.txt in/ucsc/UCSC_Main_on_Mouse_gc5BaseBw.bw

macs_peaks="${wk_folder}/recentered_distfil500_DKO_K4me3_dcm.l1.bed"


#===================================================================

outname="fig3E_CpG_content"


cut -f 1-$(($(head -1 $macs_peaks | awk '{print NF}') - 1)) $macs_peaks > unstranded.bed

ls ${an_folder}/in/test_chipseq_downstream/deeptools_heatmap_tmp/
#ls ${an_folder}/in/test_chipseq_downstream/parallel_averageBigwig

samples="${an_folder}/in/ucsc/gc5Base.bw \
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

# fig3E_CpG_content
thrsholds="60 10 80 80 40 40 80 80 30 30"

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

#--sortUsingSamples 7
#--samplesLabel ${slabels} --regionsLabel ${plabels} \
#--zMax ${thrsholds}

#========================================================
#supplementary all tracks heatmap
#========================================================

# outname="CpG_content_all_tracks"
# macs_peaks="outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"

# samples="in/ucsc/gc5Base.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Anti-GFP_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Anti-Menin_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Anti-Mll1_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Mll2_KO_Mll1_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Anti-RbBp5_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_RbBP5_average.bw \
# ./in/test_chipseq_downstream/parallel_averageBigwig/D0_WT_H3K27ac__average.bw \
# ./in/test_chipseq_downstream/parallel_averageBigwig/D0_Mll1KO_H3K27ac__average.bw \
# ./in/test_chipseq_downstream/parallel_averageBigwig/D0_Mll2KO_H3K27ac__average.bw \
# ./in/test_chipseq_downstream/parallel_averageBigwig/D0_DoubleKO_H3K27ac__average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K27me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Mll1-KO_K27me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/FC_FC_K27me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K27me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K4me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Mll1-KO_K4me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/FC_FC_K4me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K4me3_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K4me2_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Mll1-KO_K4me2_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/FC_FC_K4me2_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K4me2_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K4me1_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Mll1-KO_K4me1_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/FC_FC_K4me1_average.bw \
# ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K4me1_average.bw \
# "

# slabels="Anti-GFP_av Anti-Menin_av Anti-Mll1_av Mll2_KO_Mll1_av Anti-RbBp5_av Double_KO_RbBP5_av \
# WT_H3K27ac_av Mll2KO_H3K27ac_av Mll1KO_H3K27ac_av DoubleKO_H3K27ac_av \
# F_F_K27me3_av Mll1-KO_K27me3_av FC_FC_K27me3_av Double_KO_K27me3_av \
# F_F_K4me3_av Mll1-KO_K4me3_av FC_FC_K4me3_av Double_KO_K4me3_av \
# F_F_K4me2_av Mll1-KO_K4me2_av FC_FC_K4me2_av Double_KO_K4me2_av \
# F_F_K4me1_av Mll1-KO_K4me1_av FC_FC_K4me1_av Double_KO_K4me1_av \
# "

# thrsholds="60 10 10 10 10 15 15 10 10 10 10 10 10 10 10 40 40 40 40 30 30 30 30 10 10 10 10"
# mins="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"


# sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
# -b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
# --outFileName outs/gc_content_heatmap/${outname}_deeptools_matrix.gzip --referencePoint "center" \
# --beforeRegionStartLength 1000 --afterRegionStartLength 1500

# sudo docker run -v $outpath:$outpath -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
# --missingDataColor 0.6 --matrixFile outs/gc_content_heatmap/${outname}_deeptools_matrix.gzip \
# --outFileName outs/gc_content_heatmap/${outname}_plotHeatmap.pdf \
# --outFileNameMatrix outs/gc_content_heatmap/${outname}_plotHeatmap.mat.tab \
# --zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 16

#====================================================================================
# D4/D0 heatmap
#====================================================================================

cd ./outs/Lara_multiomic_analysis

cp /mnt/datawk1/analysis/Lara/test_chipseq_dowstream/nfout/D4_MLL2_broad_0_5/star/mergedLibrary/macs2/broadPeak/D4_MLL2B_broad_0_5_thr3.bed ./in/test_chipseq_downstream/macs_broadpeaks/
 
cp /mnt/datawk1/analysis/Lara/test_chipseq_dowstream/nfout/D4_MLL2_broad_std/star/mergedLibrary/macs2/broadPeak/D4_MLL2B_broad_std_thr3.bed ./in/test_chipseq_downstream/macs_broadpeaks/


#outname="l1_d4mll2"
outname="Dsplit_l1_d4mll2"
macs_peaks="${an_folder}/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"

doav_dko="${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP1.markdup.sorted.bigWig \
    ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP2.markdup.sorted.bigWig \
        ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4DoubleKO_REP3.markdup.sorted.bigWig"
doav_wt="${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP1.markdup.sorted.bigWig \
    ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP2.markdup.sorted.bigWig \
        ${an_folder}/in/build38_DEseq2_RNAseq/bw/D4WT_REP1.markdup.sorted.bigWig"



#rna seq removed by the showed tracks

mkdir -p outs/l1_d4_heatmap_v2


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
${avr_path}/D4_double-KO_H3K4me2_average.bw \
${avr_path}/D4_WT_H4K16ac_average.bw \
${avr_path}/D4_double-KO_H4K16ac_average.bw
"


#thrsholds="60 5 15 50 50 50 50 16 16 16 16 60 60 60 60 150 150 150 150 40 40 40 40"
#thrsholds="60 5 10 25 25 25 25 16 16 16 16 60 60 60 60 50 50 50 50 40 40"
thrsholds="60 5 10 40 40 40 40 80 80 80 80 80 80 80 80 70 70 70 70 40 40"

# thrsholds="4 4 13 13 13 13 7 7 7 7 16 16 16 16 30 30 30 30 8 8"

# mins="0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"



# slabels="gc5Base.bw \
# D4_WT_MLL2 \
# "

#thrsholds="60 5 10 8 8 8 8 10 10 10 10 20 20 20 20 10 10 5 5 5 5 0.4 0.4 0.4 0.4"
#mins="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
mins="40 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"

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
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 12


#======================================================================================
#  all other tracks two plor old s6b and s6c
#======================================================================================

# samples="${old1_avr_path}/Anti-GFP_average.bw \
# ${old1_avr_path}/Anti-Mll1_average.bw \
# ${old1_avr_path}/Mll2_KO_Mll1_average.bw \
# ${old1_avr_path}/Anti-RbBp5_average.bw \
# ${old1_avr_path}/Double_KO_RbBP5_average.bw \
# ${avr_path}/D0_WT_H3K27ac_average.bw \
# ${avr_path}/D0_Mll1-KO_H3K27ac_average.bw \
# ${avr_path}/D0_Mll2-KO_H3K27ac_average.bw \
# ${avr_path}/D0_double-KO_H3K27ac_average.bw \
# ${avr_path}/D0_WT_H3K27me3_average.bw \
# ${avr_path}/D0_Mll1-KO_H3K27me3_average.bw \
# ${avr_path}/D0_Mll2-KO_H3K27me3_average.bw \
# ${avr_path}/D0_double-KO_H3K27me3_average.bw \
# ${avr_path}/D0_WT_H3K4me3_average.bw \
# ${avr_path}/D0_Mll1-KO_H3K4me3_average.bw \
# ${avr_path}/D0_Mll2-KO_H3K4me3_average.bw \
# ${avr_path}/D0_double-KO_H3K4me3_average.bw \
# ${avr_path}/D0_WT_H3K4me2_average.bw \
# ${avr_path}/D0_Mll1-KO_H3K4me2_average.bw \
# ${avr_path}/D0_Mll2-KO_H3K4me2_average.bw \
# ${avr_path}/D0_double-KO_H3K4me2_average.bw"

    # ${avr_path}/D0_WT_H3K27ac_average.bw \
    # ${avr_path}/D0_double-KO_H3K27ac_average.bw \
    # ${avr_path}/D0_WT_H3K27me3_average.bw \
    # ${avr_path}/D0_double-KO_H3K27me3_average.bw \
    # ${avr_path}/D0_WT_H3K9me3_average.bw \
    # ${avr_path}/D0_double-KO_H3K9me3_average.bw \

outname="figS6B_CpG_content"

macs_peaks="${an_folder}/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"

samples="${an_folder}/in/ucsc/gc5Base.bw \
    ${old1_avr_path}/Anti-GFP_average.bw \
    ${old1_avr_path}/Anti-Mll1_average.bw \
    ${old1_avr_path}/Mll2_KO_Mll1_average.bw \
    ${old1_avr_path}/Anti-RbBp5_average.bw \
    ${old1_avr_path}/Double_KO_RbBP5_average.bw \
    ${avr_path}/D0_WT_H3K4me3_average.bw \
    ${avr_path}/D0_Mll1-KO_H3K4me3_average.bw \
    ${avr_path}/D0_Mll2-KO_H3K4me3_average.bw \
    ${avr_path}/D0_double-KO_H3K4me3_average.bw \
    "

# fig3E_CpG_content
thrsholds="60 10 10 10 15 15 80 80 80 80"

#thrsholds="60 10 40 40 10 10 10 10 6 6"
#thrsholds="60 10 300 300 40 40 30 30 30 30"

mins="40 0 0 0 0 0 0 0 0 0 0"

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
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 7

outname="figS6C_CpG_content"

samples="${an_folder}/in/ucsc/gc5Base.bw \
    ${old1_avr_path}/Anti-GFP_average.bw \
    ${avr_path}/D0_WT_H3K27ac_average.bw \
    ${avr_path}/D0_Mll1-KO_H3K27ac_average.bw \
    ${avr_path}/D0_Mll2-KO_H3K27ac_average.bw \
    ${avr_path}/D0_double-KO_H3K27ac_average.bw \
    ${avr_path}/D0_WT_H3K27me3_average.bw \
    ${avr_path}/D0_Mll1-KO_H3K27me3_average.bw \
    ${avr_path}/D0_Mll2-KO_H3K27me3_average.bw \
    ${avr_path}/D0_double-KO_H3K27me3_average.bw \
    ${avr_path}/D0_WT_H3K4me3_average.bw \
    "

# fig3E_CpG_content
thrsholds="60 10 40 40 40 40 80 80 80 80 80"

#thrsholds="60 10 40 40 10 10 10 10 6 6"
#thrsholds="60 10 300 300 40 40 30 30 30 30"

mins="40 0 0 0 0 0 0 0 0 0 0"

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
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 11







