#!/usr/bin/env bash
#find files in folder and save in temp folder to analyze


#generate spikein free metadata 

#SAMPLESHEET="sheets/GEO_rename_chipSEQ_bamFiles.tsv"
MAIN="${PWD}"
SUBANDIR="${MAIN}/outs/quantile_normalization_analysis/tmp/"
OUTDIR="${SUBANDIR}/Dsplit_first_batch"
METAFILE="${SUBANDIR}/meta_spikeinfree.tsv"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"

source "${DOCKER_ENV}"


# printf "ID\tANTIBODY\tGROUP\n" > "$METAFILE"
# tail -n +2 "$SAMPLESHEET" \
# | awk -F'\t' -v OFS='\t' '{print $1".bam", $4, $5}' >> "$METAFILE"

NEWMETA="${OUTDIR}/meta_spikeinfree_withday.tsv"

awk 'BEGIN{OFS="\t"} NR==1{print; next} {
    day = $1; sub(/_.*/, "", day)
    $2 = day "_" $2
    print
}' "$METAFILE" > "$NEWMETA"

echo "Salvato in: $NEWMETA"

cat $NEWMETA

#===========================================================

sudo docker run "${DOCKER_ARGS[@]}" --rm -it \
  -w ${OUTDIR} \
  --user "$(id -u)":"$(id -g)" \
  lucidif/chipseqspikeinfree:1.2.4 R

library("ChIPseqSpikeInFree")
main_folder=getwd()

#setwd("./outs/quantile_normalization_analysis/tmp/spikeinFree_wk")

metaFile <- "meta_spikeinfree_withday.tsv"

meta.info<-read.table(metaFile, sep="\t", header=TRUE)

#bams <-meta.info$ID
bams <- file.path("..", "inbam", meta.info$ID) 

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "Dsplit_first_batch_test")

quit()

#/media/lucio/easystore/Lucio/Analysis/Lara/quantile_normalization_analysis/tmp/spikeinFree_wk

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp//Dsplit_first_batch/Dsplit_first_batch_test_SF.txt \
  -b /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp//inbam/ \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp//Dsplit_first_batch/spikeinfree_bw \
  -t 100000000

INBW="${OUTDIR}/spikeinfree_bw"
#META="${PWD}/outs/quantile_normalization_analysis/tmp/meta_spikeinfree.tsv" #METAFILE
OUT="${OUTDIR}/spikeinfree_bw/average_bw"

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

