  #=====================
  #   spike in free normalization
  #=====================


#=======load functions
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
#======================================

#find files in folder and save in temp folder to analyze

chmod +x git/Lara_MLL2/bin/spikein_free_search_samples.sh

dos2unix sheets/GEO_rename_chipSEQ_bamFiles.tsv

./git/Lara_MLL2/bin/spikein_free_search_samples.sh sheets/GEO_rename_chipSEQ_bamFiles.tsv > bam_existence_report.tsv

awk -F'\t' '$1=="MISSING"{print $2 "\t" $5}' bam_existence_report.tsv

OUTDIR="./outs/quantile_normalization_analysis/tmp/inbam"
REPORT="bam_existence_report.tsv"

mkdir -p "$OUTDIR"

awk -F'\t' 'NR>1 && $1=="OK"{print $2 "\t" $5}' "$REPORT" \
| while IFS=$'\t' read -r geoname srcbam; do
    base="${geoname%.bigWig}"
    dstbam="${OUTDIR}/${base}.bam"

    echo "[CP] $srcbam -> $dstbam"
    cp -a "$srcbam" "$dstbam"

    # copia indice (se esiste)
    if [[ -f "${srcbam}.bai" ]]; then
      echo "[CP] ${srcbam}.bai -> ${dstbam}.bai"
      cp -a "${srcbam}.bai" "${dstbam}.bai"
    elif [[ -f "${srcbam%.bam}.bai" ]]; then
      echo "[CP] ${srcbam%.bam}.bai -> ${dstbam}.bai"
      cp -a "${srcbam%.bam}.bai" "${dstbam}.bai"
    fi
  done



#==========================================================

#generate spikein free metadata 

SAMPLESHEET="sheets/GEO_rename_chipSEQ_bamFiles.tsv"
OUTDIR="./outs/quantile_normalization_analysis/tmp/"
METAFILE="${OUTDIR}/meta_spikeinfree.tsv"

printf "ID\tANTIBODY\tGROUP\n" > "$METAFILE"
tail -n +2 "$SAMPLESHEET" \
| awk -F'\t' -v OFS='\t' '{print $1".mLb.mkD.sorted.bam", $4, $5}' >> "$METAFILE"


#===========================================================


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
