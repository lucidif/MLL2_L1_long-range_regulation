###########################################################################################################################
## Script for performing Quantile Normalization of read counts, after RPKM calculation, for bigwig generation in MOUSE
## Bins of 100 bps considered for qnorm and bigwig generation
## In addition, only information from Chr 1:19 + X and Y is taken into account
############################################################################################################################
#!/bin/bash
  ##Tracking time spent
  start=`date +%s`

echo "Way of running script:"
echo "1st - Run from folder where you want the bigwigs to be generated!!"
echo "2nd - Running like: bash ~/pathToQuantileNormalizedScript PathTo/bam1.bam PathTo/bamN.bam..."

#path list
outpath="outs/quantile_normalization_analysis/tmp/"

#support functions
source "git/base/bin/dockerhub_helpers.sh"

#libraries and images
IMAGE_featuresCounts="quay.io/biocontainers/subread:2.1.1--h577a1d6_0"
IMAGE_quantNorm="quantilenormalization:edger4.8.1"
IMAGE_bedgraphtobigwig="quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0"

ensure_image "$IMAGE_featuresCounts"
ensure_image "$IMAGE_quantNorm"
ensure_image "$IMAGE_bedgraphtobigwig"


# Gather all BAM files passed as arguments
#bam_files="$@"
bam_files=("$@")
echo "Bam files considered for Quantile Normalization: " $bam_files

##Loading environment for feature counts execution
##This required to initialize conda in Stardust
#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate deeptools2024

#featureCounts -a reference_peaks.saf -F SAF -o reference_peaks.countmatrix ../bams/*.bam -p -O
#Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...

# The one working, commented for testing
#featureCounts -p -O -F SAF -a ~/Documents/genomesData/safFiles/mm10_100bpWindows_onlyTargetChr.saf -o bins.countmatrix $bam_files
#featureCounts -p -O -F SAF -a in/mm10_100bpWindows_onlyTargetChr.saf -o bins.countmatrix $bam_files
sudo docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -v  /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  "$IMAGE_featuresCounts" \
  featureCounts -p -O -F SAF \
    -a in/mm10_100bpWindows_onlyTargetChr.saf \
    -o ${outpath}bins.countmatrix \
    "${bam_files[@]}"

#conda deactivate

#############################################################################
#Running R quantile normalization and generating normalized bedgraphs
#outputDir=$(pwd)
#Rscript git/Lara_MLL2/bin/quantileNorm_counts_pipelineScript.R bins.countmatrix ${outpath}

sudo docker run --rm \
  -u "$(id -u)":"$(id -g)" \
  -v /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  "$IMAGE_quantNorm" \
  Rscript git/Lara_MLL2/bin/quantileNorm_counts_pipelineScript.R bins.countmatrix "${outpath}"


###################################
##Converting bedgraphs to bigwigs

# Loop through all .bdg files in current folder
for bdg_file in ${outpath}/*.bdg; do
  # Extract the basename (filename without extension)
  base_name="${bdg_file%.bdg}"  # Correctly removes the .bdg extension

  # Generate the output bigWig filename (same base name with .bw extension)
  bw_file="${base_name}_QNORM-RPKM.bw"

  # Execute the conversion command with adapted filenames
  sudo docker run --rm \
      -u "$(id -u)":"$(id -g)" \
      -v "$PWD":"$PWD" \
      -v  /media/lucio/easystore:/media/lucio/easystore \
      -w "$PWD" \
      "$IMAGE_bedgraphtobigwig" \
      bedGraphToBigWig "$bdg_file" in/mm10.onlyTargetChrom.sizes "$bw_file"

    echo "Converted $bdg_file to $bw_file"
done

##Remove bedgraphs and count matrix to avoid high storage accumulation
#rm ./bins.countmatrix
#rm ./bins.countmatrix.summary
#rm ./*.bdg

#########################
##Tracking time spent
end=`date +%s`

runtime=$((end-start))

echo ''
echo 'Total time spent: '
echo $runtime