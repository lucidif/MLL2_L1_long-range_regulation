#!/usr/bin/env bash

prj_folder=${PWD}

mkdir -p outs/gc_content_heatmap


outname="CpG_content"
macs_peaks="outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"

cut -f 1-$(($(head -1 $macs_peaks | awk '{print NF}') - 1)) $macs_peaks > unstranded.bed

doav_dko="in/build38_DEseq2_RNAseq/bw/D0DoubleKO_REP1.markdup.sorted.bigWig \
    in/build38_DEseq2_RNAseq/bw/D0DoubleKO_REP2.markdup.sorted.bigWig \
        in/build38_DEseq2_RNAseq/bw/D0DoubleKO_REP3.markdup.sorted.bigWig"
doav_wt="in/build38_DEseq2_RNAseq/bw/D0WTA_REP1.markdup.sorted.bigWig \
    in/build38_DEseq2_RNAseq/bw/D0WTA_REP2.markdup.sorted.bigWig \
        in/build38_DEseq2_RNAseq/bw/D0WTA_REP1.markdup.sorted.bigWig"


#rna seq removed by the showed tracks

sudo docker run -v /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis -w `pwd` quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b $doav_dko -o outs/gc_content_heatmap/dko_rnaseq.bw

sudo docker run -v /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis -w `pwd` quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b $doav_wt -o outs/gc_content_heatmap/wt_rnaseq.bw

#plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"

#slabels="Anti-GFP_av Anti-Menin_av Anti-Mll1_av Mll2_KO_Mll1_av Anti-RbBp5_av Double_KO_RbBP5_av \
#WT_H3K27ac_C Mll2KO_H3K27ac_C Mll1KO_H3K27ac_C DoubleKO_H3K27ac_C \
#F_F_K27me3_av FC_FC_K27me3_av Mll1-KO_K27me3_av Double_KO_K27me3_av \
#F_F_K4me3_av FC_FC_K4me3_av Mll1-KO_K4me3_av Double_KO_K4me3_av \
#F_F_K4me2_av FC_FC_K4me2_av Mll1-KO_K4me2_av Double_KO_K4me2_av \
#F_F_K4me1_av FC_FC_K4me1_av Mll1-KO_K4me1_av Double_KO_K4me1_av \
#"

#to add
#H3K27ac #H3K27me3 #H3K9me3 to 0-6 

ls ./in/test_chipseq_downstream/deeptools_heatmap_tmp/
ls ./in/test_chipseq_downstream/parallel_averageBigwig

samples="in/ucsc/gc5Base.bw \
    ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Anti-GFP_average.bw \
    ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K4me3_average.bw \
    ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K4me3_average.bw \
    ./in/test_chipseq_downstream/parallel_averageBigwig/D0_WT_H3K27ac__average.bw \
    ./in/test_chipseq_downstream/parallel_averageBigwig/D0_DoubleKO_H3K27ac__average.bw \
    ./in/test_chipseq_downstream/deeptools_heatmap_tmp/F_F_K27me3_average.bw \
    ./in/test_chipseq_downstream/deeptools_heatmap_tmp/Double_KO_K27me3_average.bw \
    ./in/test_chipseq_downstream/parallel_averageBigwig/D0_WT_H3K9me3__average.bw \
    ./in/test_chipseq_downstream/parallel_averageBigwig/D0_KO_H3K9me3__average.bw \
    "


thrsholds="60 10 40 40 10 10 10 10 6 6"
mins="40 0 0 0 0 0 0 0 0 0"

sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S ${samples} -R ${macs_peaks} \
-b 1000 --averageTypeBins "median" --numberOfProcessors 8 --scale 1 --binSize 10 \
--outFileName outs/gc_content_heatmap/${outname}_deeptools_matrix.gzip --referencePoint "center" \
--beforeRegionStartLength 1000 --afterRegionStartLength 1500

sudo docker run -v $outpath:$outpath -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'seismic' \
--missingDataColor 0.6 --matrixFile outs/gc_content_heatmap/${outname}_deeptools_matrix.gzip \
--outFileName outs/gc_content_heatmap/${outname}_plotHeatmap.pdf \
--outFileNameMatrix outs/gc_content_heatmap/${outname}_plotHeatmap.mat.tab \
--zMax ${thrsholds} --zMin ${mins} --sortUsingSamples 3

#--sortUsingSamples 7
#--samplesLabel ${slabels} --regionsLabel ${plabels} \
#--zMax ${thrsholds}