sudo apt install python3-deeptools

chrom_sizes="outs/Lara_multiomic_analysis/in/mm10.sizes"
infolder="outs/geo_sub/out/MLL2_chipSEQ_processed/"
outfolder="outs/quantile_normalization_analysis/tmp/"

#binize referece
bedtools makewindows -g ${chrom_sizes} -w 50000 -s 50000 > ${outfolder}mm10.50kb.bed

multiBamSummary BED-file –BED ${outfolder}mm10.50kb.bed -b bam_list -out results.npz –outRawCounts output_rawCount.txt

mkdir ${outfolder}

multiBigwigSummary bins \
  --binSize 50 \
  --bwfiles \
    ${infolder}*.bigWig \
  --outFileName ${outfolder}allMarks_allSamples_bins50.npz \
  --outRawCounts ${outfolder}allMarks_allSamples_bins50.tsv



