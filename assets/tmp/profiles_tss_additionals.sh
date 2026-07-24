#!/usr/bin/env bash

inpath="${PWD}/outs/Lara_multiomic_analysis"

prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
#old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
#outpath="${PWD}/outs/Lara_multiomic_analysis/outs/gc_content_heatmap/"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"



outpath="/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis"


#=============RIG1B=================



samples="./in/test_chipseq_downstream/deeptools_heatmaps/D0_WT_RING1B__average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/D0_KO_RING1B__average.bw \
"

title="RIG1B"
outname="DEGs_${title}_only"

macs_peaks="outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed"
plabels="broad narrow negative"
#allpeaks="input/fold2_th0_05_K4me2_Double_KO_vs_F_F_DBA_DESEQ2_DOWN.bed"


slabels="WT DKO"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
-e MPLCONFIGDIR=$outpath/.matplotlib \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) -e MPLCONFIGDIR=$outpath/.matplotlib quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip -out ${outpath}/outs/downGenesHeatmap/${outname}_Profile.png --perGroup --colors "#715eee""#ffb00e" --plotTitle "${title} loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}



# data_d0[grep("D0DoubleKO_",data_d0[,1]),2]<-"#ffb00e" #"DKO"
# data_d0[grep("D0Mll1KO_",data_d0[,1]),2]<-"#de217d" #"Mll1KO"
# data_d0[grep("D0Mll2KO_",data_d0[,1]),2]<-"#ff5f01" #"Mll2KO"
# data_d0[grep("D0WTA_",data_d0[,1]),2]<-"#715eee" #"WTA"

#=================

#=============H3K27me3=================


title="H3K27me3"
outname="DEGs_${title}_only"
macs_peaks="outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed"
plabels="broad narrow negative"

samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K27me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K27me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K27me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K27me3_average.bw \
"

slabels="WT Mll1KO Mll2KO DKO"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip -out ${outpath}/outs/downGenesHeatmap/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}




#=============K4me3=================

title="K4me3"
outname="DEGs_${title}_only"
macs_peaks="outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed"
plabels="broad narrow negative"

samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K4me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K4me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K4me3_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K4me3_average.bw \
"

slabels="WT Mll1KO Mll2KO DKO"


sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip -out ${outpath}/outs/downGenesHeatmap/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}


#=============K4me2=================

title="K4me2"
outname="DEGs_${title}_only"
macs_peaks="outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed \
outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed"
plabels="broad narrow negative"

samples="./in/test_chipseq_downstream/deeptools_heatmaps/F_F_K4me2_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Mll1-KO_K4me2_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/FC_FC_K4me2_average.bw \
./in/test_chipseq_downstream/deeptools_heatmaps/Double_KO_K4me2_average.bw \
"

slabels="WT Mll1KO Mll2KO DKO"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R $macs_peaks \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotProfile -m ${outpath}/outs/downGenesHeatmap/${outname}_deeptools_matrix.gzip -out ${outpath}/outs/downGenesHeatmap/${outname}_Profile.png --perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" --plotTitle "${title} loss K4me3 in dko vs wt" --samplesLabel ${slabels} --regionsLabel ${plabels}



#===========================RIG1B
