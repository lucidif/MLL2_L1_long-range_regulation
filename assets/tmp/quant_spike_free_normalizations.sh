#import bam files

#H3K4me3

#D0

#D0_WT_H3K4me3_repA
#D0_WT_H3K4me3_repB
#D0_Mll1-KO_H3K4me3_repA
#D0_Mll1-KO_H3K4me3_repB
#D0_Mll2-KO_H3K4me3_repA
#D0_Mll2-KO_H3K4me3_repB
#D0_double-KO_H3K4me3_repA
#D0_double-KO_H3K4me3_repB

rep1names="D0_WT_H3K4me3_repA D0_Mll1-KO_H3K4me3_repA D0_Mll2-KO_H3K4me3_repA D0_double-KO_H3K4me3_repA"

# rep1samples="outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K4me3_A.mLb.mkD.sorted.bam \
# outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K4me3.mLb.mkD.sorted.bam \
# outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K4me3_A.mLb.mkD.sorted.bam \
# outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K4me3_A.mLb.mkD.sorted.bam
# "

rep1samples=(
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K4me3_A.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K4me3.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K4me3_A.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K4me3_A.mLb.mkD.sorted.bam
)

bash git/Lara_MLL2/bin/quantileNormalizingBigwigs_Mouse.sh "${rep1samples[@]}"
mv ./outs/quantile_normalization_analysis/tmp/*.bw ./outs/quantile_normalization_analysis/bw/

rep2names="D0_WT_H3K4me3_repB D0_Mll1-KO_H3K4me3_repB D0_Mll2-KO_H3K4me3_repB D0_double-KO_H3K4me3_repB"

rep2samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K4me3_B.mLb.mkD.sorted.bam
)

bash git/Lara_MLL2/bin/quantileNormalizingBigwigs_Mouse.sh "${rep2samples[@]}"
mv ./outs/quantile_normalization_analysis/tmp/*.bw ./outs/quantile_normalization_analysis/bw/

#renames files
bwdir="./outs/quantile_normalization_analysis/bw"

# nomi nuovi (array)
rep1names=(D0_WT_H3K4me3_repA D0_Mll1-KO_H3K4me3_repA D0_Mll2-KO_H3K4me3_repA D0_double-KO_H3K4me3_repA)
rep2names=(D0_WT_H3K4me3_repB D0_Mll1-KO_H3K4me3_repB D0_Mll2-KO_H3K4me3_repB D0_double-KO_H3K4me3_repB)

# file vecchi (array) - in ordine coerente con i nomi sopra
repA_files=(
  "$bwdir/F_F_K4me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Mll1-KO_K4me3.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/FC_FC_K4me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Double_KO_K4me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
)

repB_files=(
  "$bwdir/F_F_K4me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Mll1-KO_K4me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/FC_FC_K4me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Double_KO_K4me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
)

# rinomina (mv -n = non sovrascrive se esiste già)
for i in "${!repA_files[@]}"; do
  src="${repA_files[$i]}"
  dst="$bwdir/${rep1names[$i]}.bw"
  [[ -e "$src" ]] || { echo "Manca: $src"; exit 1; }
  mv -n -- "$src" "$dst"
  echo "OK: $(basename "$src") -> $(basename "$dst")"
done

for i in "${!repB_files[@]}"; do
  src="${repB_files[$i]}"
  dst="$bwdir/${rep2names[$i]}.bw"
  [[ -e "$src" ]] || { echo "Manca: $src"; exit 1; }
  mv -n -- "$src" "$dst"
  echo "OK: $(basename "$src") -> $(basename "$dst")"
done

#H3K27me3

#D0

#"outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/"

#TODO :

rep1names="D0_WT_H3K27me3_repA D0_Mll1-KO_H3K27me3_repA D0_Mll2-KO_H3K27me3_repA D0_double-KO_H3K27me3_repA" 


rep1samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K27me3_A.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K27me3.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K27me3_A.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K27me3_A.mLb.mkD.sorted.bam
)

bash git/Lara_MLL2/bin/quantileNormalizingBigwigs_Mouse.sh "${rep1samples[@]}"
mv ./outs/quantile_normalization_analysis/tmp/*.bw ./outs/quantile_normalization_analysis/bw/


rep2names="D0_WT_H3K27me3_repB D0_Mll1-KO_H3K27me3_repB D0_Mll2-KO_H3K27me3_repB D0_double-KO_H3K27me3_repB"

rep2samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K27me3_B.mLb.mkD.sorted.bam
)

bash git/Lara_MLL2/bin/quantileNormalizingBigwigs_Mouse.sh "${rep2samples[@]}"
mv ./outs/quantile_normalization_analysis/tmp/*.bw ./outs/quantile_normalization_analysis/bw/

# renames files (H3K27me3)
bwdir="./outs/quantile_normalization_analysis/bw"


rep1names=(D0_WT_H3K27me3_repA D0_Mll1-KO_H3K27me3_repA D0_Mll2-KO_H3K27me3_repA D0_double-KO_H3K27me3_repA)
rep2names=(D0_WT_H3K27me3_repB D0_Mll1-KO_H3K27me3_repB D0_Mll2-KO_H3K27me3_repB D0_double-KO_H3K27me3_repB)

# file vecchi (array) - in ordine coerente con i nomi sopra
repA_files=(
  "$bwdir/F_F_K27me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Mll1-KO_K27me3.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/FC_FC_K27me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Double_KO_K27me3_A.mLb.mkD.sorted_QNORM-RPKM.bw"
)

repB_files=(
  "$bwdir/F_F_K27me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Mll1-KO_K27me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/FC_FC_K27me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
  "$bwdir/Double_KO_K27me3_B.mLb.mkD.sorted_QNORM-RPKM.bw"
)

# rinomina (mv -n = non sovrascrive se esiste già)
for i in "${!repA_files[@]}"; do
  src="${repA_files[$i]}"
  dst="$bwdir/${rep1names[$i]}.bw"
  [[ -e "$src" ]] || { echo "Manca: $src"; exit 1; }
  mv -n -- "$src" "$dst"
  echo "OK: $(basename "$src") -> $(basename "$dst")"
done

for i in "${!repB_files[@]}"; do
  src="${repB_files[$i]}"
  dst="$bwdir/${rep2names[$i]}.bw"
  [[ -e "$src" ]] || { echo "Manca: $src"; exit 1; }
  mv -n -- "$src" "$dst"
  echo "OK: $(basename "$src") -> $(basename "$dst")"
done

#make average bigwig

# bash ./git/chipseq_downstream_macs/bin/parallel_averageBigwig.sh \
# "./outs/quantile_normalization_analysis/bw/" \
# "D0_WT_H3K4me3 D0_Mll1-KO_H3K4me3 D0_Mll2-KO_H3K4me3 D0_double-KO_H3K4me3 D0_WT_H3K27me3 D0_Mll1-KO_H3K27me3 D0_Mll2-KO_H3K27me3 D0_double-KO_H3K27me3" \
# ./outs/quantile_normalization_analysis/average_bw

bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "${PWD}/outs/quantile_normalization_analysis/bw/" \
  "D0_WT_H3K4me3 D0_Mll1-KO_H3K4me3 D0_Mll2-KO_H3K4me3 D0_double-KO_H3K4me3 D0_WT_H3K27me3 D0_Mll1-KO_H3K27me3 D0_Mll2-KO_H3K27me3 D0_double-KO_H3K27me3" \
  ${PWD}/outs/quantile_normalization_analysis/average_bw


  #=====================
  #   spike in free normalization
  #=====================

copy_and_rename() {
  local names_str="$1"; shift
  local -a samples=("$@")

  # split della stringa nomi in array
  read -r -a names <<< "$names_str"

  if [[ ${#names[@]} -ne ${#samples[@]} ]]; then
    echo "[ERROR] mismatch: ${#names[@]} names vs ${#samples[@]} samples" >&2
    return 1
  fi

  for i in "${!samples[@]}"; do
    src="${samples[$i]}"
    dst="${outdir}/${names[$i]}.bam"

    if [[ ! -f "$src" ]]; then
      echo "[ERROR] file non trovato: $src" >&2
      return 1
    fi

    echo "[INFO] cp $src -> $dst"
    cp -a "$src" "$dst"
  done
}


outdir="./outs/quantile_normalization_analysis/tmp"

H3K4me3_rep1names="D0_WT_H3K4me3_repA D0_Mll1-KO_H3K4me3_repA D0_Mll2-KO_H3K4me3_repA D0_double-KO_H3K4me3_repA"

H3K4me3_rep1samples=(
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K4me3_A.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K4me3.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K4me3_A.mLb.mkD.sorted.bam
  outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K4me3_A.mLb.mkD.sorted.bam
)

H3K4me3_rep2names="D0_WT_H3K4me3_repB D0_Mll1-KO_H3K4me3_repB D0_Mll2-KO_H3K4me3_repB D0_double-KO_H3K4me3_repB"

H3K4me3_rep2samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K4me3_B.mLb.mkD.sorted.bam 
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K4me3_B.mLb.mkD.sorted.bam
)

H3K27me3_rep1names="D0_WT_H3K27me3_repA D0_Mll1-KO_H3K27me3_repA D0_Mll2-KO_H3K27me3_repA D0_double-KO_H3K27me3_repA" 

H3K27me3_rep1samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K27me3_A.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K27me3.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K27me3_A.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K27me3_A.mLb.mkD.sorted.bam
)

H3K27me3_rep2names="D0_WT_H3K27me3_repB D0_Mll1-KO_H3K27me3_repB D0_Mll2-KO_H3K27me3_repB D0_double-KO_H3K27me3_repB"

H3K27me3_rep2samples=(
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/F_F_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Mll1-KO_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/FC_FC_K27me3_B.mLb.mkD.sorted.bam
    outs/Lara_ChIP_with_spikein_A1/star/mergedLibrary/Double_KO_K27me3_B.mLb.mkD.sorted.bam
)

# --------- ESECUZIONE ---------

copy_and_rename "$H3K4me3_rep1names" "${H3K4me3_rep1samples[@]}"
copy_and_rename "$H3K4me3_rep2names" "${H3K4me3_rep2samples[@]}"
copy_and_rename "$H3K27me3_rep1names" "${H3K27me3_rep1samples[@]}"
copy_and_rename "$H3K27me3_rep2names" "${H3K27me3_rep2samples[@]}"

cp ./sheets/SpikeinFree_wkss.txt ./outs/quantile_normalization_analysis/tmp/

sudo docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -v  /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  lucidif/chipseqspikeinfree:1.2.4 R

library("ChIPseqSpikeInFree")

setwd("./outs/quantile_normalization_analysis/tmp")

metaFile <- "./SpikeinFree_wkss.txt"
bams <- c("D0_WT_H3K4me3_repA.bam",
          "D0_Mll1-KO_H3K4me3_repA.bam",
          "D0_Mll2-KO_H3K4me3_repA.bam",
          "D0_double-KO_H3K4me3_repA.bam",
          "D0_WT_H3K4me3_repB.bam",
          "D0_Mll1-KO_H3K4me3_repB.bam",
          "D0_Mll2-KO_H3K4me3_repB.bam",
          "D0_double-KO_H3K4me3_repB.bam",
          "D0_WT_H3K27me3_repA.bam",
          "D0_Mll1-KO_H3K27me3_repA.bam",
          "D0_Mll2-KO_H3K27me3_repA.bam",
          "D0_double-KO_H3K27me3_repA.bam",
          "D0_WT_H3K27me3_repB.bam",
          "D0_Mll1-KO_H3K27me3_repB.bam",
          "D0_Mll2-KO_H3K27me3_repB.bam",
          "D0_double-KO_H3K27me3_repB.bam"  
          )

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "test")

quit()


sudo docker run --rm -it \
  -u "$(id -u)":"$(id -g)" \
  -v /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  quantilenormalization:edger4.8.1 R

library("edgeR")

# SF + metadata
dat <- read.table("./outs/quantile_normalization_analysis/tmp/test_SF.txt",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE,
                  quote="", check.names=FALSE, fill=TRUE)

# rawCounts 
rc <- read.table("./outs/quantile_normalization_analysis/tmp/test_rawCounts.txt",
                 sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

coord_cols <- 1:3
counts_all <- as.matrix(rc[, -coord_cols])
rownames(counts_all) <- paste(rc[[1]], rc[[2]], rc[[3]], sep=":")

# match names (perche?)
dat$key <- sub("\\.bam$", "", dat$ID)
cnt_key <- sub("\\.bam$", "", colnames(counts_all))

dge_list <- list()

for (ab in unique(dat$ANTIBODY)) {

  d <- subset(dat, ANTIBODY == ab & QC == "pass" & !is.na(SF))

  # quali colonne della counts appartengono a questo anticorpo?
  keep <- cnt_key %in% d$key
  counts <- counts_all[, keep, drop=FALSE]

  # riallinea d all’ordine delle colonne
  d2 <- d[match(sub("\\.bam$", "", colnames(counts)), d$key), ]
  stopifnot(!any(is.na(d2$key)))

  SF    <- as.numeric(d2$SF)
  GROUP <- factor(d2$GROUP)

  dge <- DGEList(counts=counts, group=GROUP, norm.factors=SF)

  dge_list[[ab]] <- list(dge=dge, meta=d2)
}

quit()

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s ./outs/quantile_normalization_analysis/tmp/test_SF.txt \
  -b ./outs/quantile_normalization_analysis/tmp \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o ./outs/quantile_normalization_analysis/tmp \
  -t 100000000

bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "${PWD}/outs/quantile_normalization_analysis/tmp/H3K4me3/" \
  "D0_WT_H3K4me3 D0_Mll1-KO_H3K4me3 D0_Mll2-KO_H3K4me3 D0_double-KO_H3K4me3" \
  "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinfree_tracks_average_bw"


bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "${PWD}/outs/quantile_normalization_analysis/tmp/H3K27me3/" \
  "D0_WT_H3K27me3 D0_Mll1-KO_H3K27me3 D0_Mll2-KO_H3K27me3 D0_double-KO_H3K27me3" \
  "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinfree_tracks_average_bw"

inpath="${PWD}/outs/quantile_normalization_analysis/tmp/H3K27me3"
i="D0_Mll1-KO_H3K27me3"

infiles=( "$inpath"/"$i"*.bw )
printf "N=%s\n" "${#infiles[@]}"
printf "%s\n" "${infiles[@]}"

#test average
sudo docker run --rm \
  -e MPLCONFIGDIR="$mpldir" \
  -e TMPDIR="$tmpdir" -e TMP="$tmpdir" -e TEMP="$tmpdir" \
  -v /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  -u "$(id -u)":"$(id -g)" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bigwigAverage -b ${PWD}/outs/quantile_normalization_analysis/tmp/H3K4me3/D0_Mll1-KO_H3K4me3_*.bw -o "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinfree_tracks_average_bw/D0_Mll1-KO_H3K4me3_average.bw" -p 8

in1="${PWD}/outs/quantile_normalization_analysis/tmp/H3K4me3/D0_Mll1-KO_H3K4me3_repA.spikeinfree.bw"
in2="${PWD}/outs/quantile_normalization_analysis/tmp/H3K4me3/D0_Mll1-KO_H3K4me3_repB.spikeinfree.bw"
out="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinfree_tracks_average_bw/D0_Mll1-KO_H3K4me3_average_TEST.bw"

rm -f "$out"

sudo docker run --rm \
  -e MPLCONFIGDIR="$mpldir" \
  -e TMPDIR="$tmpdir" -e TMP="$tmpdir" -e TEMP="$tmpdir" \
  -v "$PWD":"$PWD" \
  -v /media/lucio/easystore:/media/lucio/easystore \
  -u "$(id -u)":"$(id -g)" \
  -w "$PWD" \
  -u "$(id -u)":"$(id -g)" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bigwigAverage -b "$in1" "$in2" -o "${out}" -p 8

echo "exit=$?"
ls -lh "$out"


bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "${PWD}/outs/quantile_normalization_analysis/tmp/H3K4me3/" \
  "D0_Mll1-KO_H3K4me3" \
  "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinfree_tracks_average_bw"



D0_Mll1-KO_H3K4me3_repA.bam
