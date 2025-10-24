#!/bin/bash

bigwig=$1
macs_peaks=$2
outpath=$3
inpath=$4
outname=$5

bedir=$(dirname "$bedfile")

sed 's/ /\t/g' ${macs_peaks} | cut -f 1-6 > ${outpath}/coordinate.bed

echo "bigwig= $bigwig"
echo "macs_peaks= $macs_peaks"
echo "outpath= $outpath"
echo "inpath= $inpath"
echo "sampleLabel= $sampleLabel"


#echo "sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 computeMatrix reference-point -S $inpath/$bigwig -R ${outpath}/coordinate.bed -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 --averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip --beforeRegionStartLength 5000 --afterRegionStartLength 5000" --referencePoint "center"

#sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 computeMatrix reference-point -S $inpath/$bigwig -R ${outpath}/coordinate.bed -b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 --averageTypeBins "median" --outFileName $outpath/${outname}_deeptools_matrix.gzip --beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 plotHeatmap --colorMap 'viridis' --missingDataColor 0.6 --matrixFile $outpath/${outname}_deeptools_matrix.gzip --outFileName $outpath/plotHeatmap.pdf --outFileNameMatrix $outpath/plotHeatmap.mat.tab


#--sortUsingSamples
#List of sample numbers (order as in matrix), that are used for sorting by –sortUsing, no value uses all samples, example: –sortUsingSamples 1 3

#--samplesLabel
#Labels for the samples. This will then be passed to plotHeatmap and plotProfile. The default is to use the file name of the sample. The sample labels should be separated by spaces and quoted if a label itselfcontains a space E.g. –samplesLabel label-1 “label 2”

#--scale
#If set, all values are multiplied by this number. (Default: 1)

#--numberOfProcessors, -p
#Number of processors to use. Type “max/2” to use half the maximum number of processors or “max” to use all available processors. (Default: 1)

#--binSize, -bs
#Length, in bases, of the non-overlapping bins for averaging the score over the regions length. (Default: 10)

#--averageTypeBins
#Possible choices: mean, median, min, max, std, sum

#--minThreshold
#Numeric value. Any region containing a value that is less than or equal to this will be skipped. This is useful to skip, for example, genes where the read count is zero for any of the bins. This could be the result of unmappable areas and can bias the overall results. (Default: None)

#--maxThreshold
#Numeric value. Any region containing a value greater than or equal to this will be skipped. The maxThreshold is useful to skip those few regions with very high read counts (e.g. micro satellites) that may bias the average values. (Default: None)

