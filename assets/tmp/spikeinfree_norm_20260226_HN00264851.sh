#!/usr/bin/env bash

main_prj_folder=${PWD}

#dos2unix sheets/GEO_rename_chipSEQ_bamFiles.tsv

BAMPATH="${main_prj_folder}/outs/20260226_HN00264851_alln/nfout/random/star/mergedLibrary/"

ls ${BAMPATH}

SAMPLESHEET="sheets/GEO_rename_chipSEQ_bamFiles.tsv"
OUTDIR="${main_prj_folder}/outs/quantile_normalization_analysis/tmp/"
METAFILE="${OUTDIR}/meta_spikeinfree.tsv"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"

source "${DOCKER_ENV}"

add_run_samples="D4_double-KO_H3K9me3_repA.mLb.mkD.sorted.bam \
D4_double-KO_H3K9me3_repB.mLb.mkD.sorted.bam \
D4_double-KO_H4K16ac_repA.mLb.mkD.sorted.bam \
D4_double-KO_H4K16ac_repB.mLb.mkD.sorted.bam \
D4_double-KO_input.mLb.mkD.sorted.bam \
D0_cdMll2_input.mLb.mkD.sorted.bam \
D0_cdMll2_H3K4me3_repA.mLb.mkD.sorted.bam \
D4_WT_H3K9me3.mLb.mkD.sorted.bam \
D4_WT_H3K9me3_2.mLb.mkD.sorted.bam \
WT_H4K16ac_repA.mLb.mkD.sorted.bam \
WT_H4K16ac_repB.mLb.mkD.sorted.bam \
D4_WT_input.mLb.mkD.sorted.bam"

for i in $add_run_samples ; do echo "bam copy: ${i}" ; cp ${BAMPATH}/$i ${OUTDIR}/inbam/ ; done

for i in $add_run_samples ; do echo "bai copy: ${i}.bai" ; cp ${BAMPATH}/${i}.bai ${OUTDIR}/inbam/ ; done

#import other bam files in inbam folder

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
  --user "$(id -u)":"$(id -g)" \
  lucidif/chipseqspikeinfree:1.2.4 R

library("ChIPseqSpikeInFree")
main_folder=getwd()

setwd("./outs/quantile_normalization_analysis/tmp/20260226_HN00264851_spikeinFree")

metaFile <- paste0(main_folder,"/sheets/spikeinfree_20260226_HN00264851.tsv")

meta.info<-read.table(metaFile, sep="\t", header=TRUE)

#bams <-meta.info$ID
bams <- file.path("..", "inbam", meta.info$ID) 

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "test")

#TODO this for see if something change

setwd(paste0(main_folder,"/outs/quantile_normalization_analysis/tmp/Dsplit_20260226_HN00264851_spikeinFree"))
metaFile <- paste0(main_folder,"/sheets/Dsplit_spikeinfree_20260226_HN00264851.tsv")

meta.info<-read.table(metaFile, sep="\t", header=TRUE)

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "Dsplit")

quit()

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s ./outs/quantile_normalization_analysis/tmp/20260226_HN00264851_spikeinFree/test_SF.txt \
  -b ./outs/quantile_normalization_analysis/tmp/inbam \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o ./outs/quantile_normalization_analysis/tmp/20260226_HN00264851_spikeinFree \
  -t 100000000

#dsplit

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s ./outs/quantile_normalization_analysis/tmp/Dsplit_20260226_HN00264851_spikeinFree/Dsplit_SF.txt \
  -b ./outs/quantile_normalization_analysis/tmp/inbam \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o ./outs/quantile_normalization_analysis/tmp/Dsplit_20260226_HN00264851_spikeinFree \
  -t 100000000

#====================

INBW="${PWD}/outs/quantile_normalization_analysis/tmp/20260226_HN00264851_spikeinFree/"
#META="${PWD}/outs/quantile_normalization_analysis/tmp/meta_spikeinfree.tsv" #METAFILE
OUT="${PWD}/outs/quantile_normalization_analysis/tmp/20260226_HN00264851_spikeinFree/average_bw"

mkdir -p $INBW/tmp/
cp $INBW/*/*.bw $INBW/tmp/


PATTERNS="$(
  awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++) if($i=="ID") c=i; next}
    {print $c}
  ' "$METAFILE" \
  | sed -E 's/\.bam$//' \
  | sed -E 's/_rep[AB]$//' \
  | sort -u \
  | tr '\n' ' ' \
  | sed 's/[[:space:]]*$//'
)"

echo "PATTERNS: $PATTERNS"

mkdir -p "$OUT"

bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "$INBW/tmp" \
  "$PATTERNS" \
  "$OUT"

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b ${INBW}/H4K16ac/D4_double-KO_H4K16ac_rep*.bw -o "${INBW}/average_bw/D4_double-KO_H4K16ac_average.bw" -p 8

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b ${INBW}/H4K16ac/D4_double-KO_H4K16ac_rep*.bw -o "${INBW}/average_bw/D4_double-KO_H4K16ac_average.bw" -p 8

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b ${INBW}/H4K16ac/D4_WT_H4K16ac_rep*.bw -o "${INBW}/average_bw/D4_WT_H4K16ac_average.bw" -p 8

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b ${INBW}/tmp/D4_WT_H3K9me3*.mLb.mkD.sorted.spikeinfree.bw -o "${INBW}/average_bw/D4_WT_H3K9me3_average.bw" -p 8

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b ${INBW}/tmp/D4_double-KO_H3K9me3_rep*.mLb.mkD.sorted.spikeinfree.bw -o "${INBW}/average_bw/D4_double-KO_H3K9me3_average.bw" -p 8

