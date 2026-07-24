## plot profile separately
## --perGroup make one image per BED file instead of per bigWig file
## you must split the in matrix by modification

#all preaks file

prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
outpath="${PWD}/outs/quantile_normalization_analysis/d0_heatmaps"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"


macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"

plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

#allpeaks="input/fold2_th0_05_K4me2_Double_KO_vs_F_F_DBA_DESEQ2_DOWN.bed"


#====================K3k27ac

samples="${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_Mll1-KO_H3K27ac_average.bw \
${avr_path}/D0_Mll2-KO_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw"

slabels="WT Mll1KO Mll2KO DKO"


outname="K4me3_DOWN_Double_KO_vs_F_F_H3K27ac_only"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run \
# -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $allpeaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/allpeaks_${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip \
-out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" \
--plotTitle "H3K27ac loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/allpeaks_${outname}_deeptools_matrix.gzip -out ${outpath}/allpeaks_${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "H3K27ac loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel "all preaks"


# data_d0[grep("D0DoubleKO_",data_d0[,1]),2]<-"#ffb00e" #"DKO"
# data_d0[grep("D0Mll1KO_",data_d0[,1]),2]<-"#de217d" #"Mll1KO"
# data_d0[grep("D0Mll2KO_",data_d0[,1]),2]<-"#ff5f01" #"Mll2KO"
# data_d0[grep("D0WTA_",data_d0[,1]),2]<-"#715eee" #"WTA"

#=================

#=============H3K27me3=================

title="H3K27me3"
outname="K4me3_DOWN_Double_KO_vs_F_F_${title}_only"
#macs_peaks="K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

# samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K27me3_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K27me3_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K27me3_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K27me3_average.bw"

samples="${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K27me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw"

slabels="WT Mll1KO Mll2KO DKO"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $allpeaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/allpeaks_${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}

# sudo docker "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \ 
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/allpeaks_${outname}_deeptools_matrix.gzip -out ${outpath}/allpeaks_${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel "all preaks"



#=============K4me3=================

title="H3K4me3"
outname="K4me3_DOWN_Double_KO_vs_F_F_${title}_only"
#macs_peaks="K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

samples="${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw"


slabels="WT Mll1KO Mll2KO DKO"


sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run \
# "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $allpeaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/allpeaks_${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/allpeaks_${outname}_deeptools_matrix.gzip -out ${outpath}/allpeaks_${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel "all preaks"

#=============K4me2=================

title="H3K4me2"
outname="H3K4me2_DOWN_Double_KO_vs_F_F_${title}_only"
#macs_peaks="H3K4me2_proximal_CpG_plus_Double_KO_vs_F_F.bed H3K4me2_proximal_CpG_minus_Double_KO_vs_F_F.bed H3K4me2_distal_CpG_plus_Double_KO_vs_F_F.bed H3K4me2_distal_CpG_minus_Double_KO_vs_F_F.bed"
plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

samples="${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me2_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw"

#samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K4me2_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K4me2_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K4me2_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K4me2_average.bw"

slabels="WT Mll1KO Mll2KO DKO"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $allpeaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/allpeaks_${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run \
"${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) -e MPLCONFIGDIR=${prj_folder}/${outpath}/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/allpeaks_${outname}_deeptools_matrix.gzip -out ${outpath}/allpeaks_${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel "all preaks"

#=============K4me1=================

# title="K4me1"
# outname="K4me3_DOWN_Double_KO_vs_F_F_${title}_only"
# macs_peaks="K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
# plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

# samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K4me1_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K4me1_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K4me1_average.bw ./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K4me1_average.bw"

# slabels="WT Mll1KO Mll2KO DKO"



# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $macs_peaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
# quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
# computeMatrix reference-point -S ${samples} -R $allpeaks \
# -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
# --averageTypeBins "median" --outFileName $outpath/allpeaks_${outname}_deeptools_matrix.gzip \
# --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/${outname}_deeptools_matrix.gzip -out ${outpath}/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}

# sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m $outpath/allpeaks_${outname}_deeptools_matrix.gzip -out ${outpath}/allpeaks_${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} lose K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel "all preaks"
