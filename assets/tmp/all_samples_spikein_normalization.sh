#!/usr/bin/env bash
#find files in folder and save in temp folder to analyze

main_prj_folder=${PWD}

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

cp /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/HN00226631_HN00220667_chip/nfout/random/star/mergedLibrary/D4_double-KO_H3K27me3_repB.mLb.mkD.sorted.bam ./outs/quantile_normalization_analysis/tmp/inbam/D4_double-KO_H3K27me3_repB.bam

cp /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/HN00226631_HN00220667_chip/nfout/random/star/mergedLibrary/D4_double-KO_H3K27me3_repB.mLb.mkD.sorted.bam.bai ./outs/quantile_normalization_analysis/tmp/inbam/D4_double-KO_H3K27me3_repB.bam.bai

#==========================================================

#generate spikein free metadata 

SAMPLESHEET="sheets/GEO_rename_chipSEQ_bamFiles.tsv"
OUTDIR="./outs/quantile_normalization_analysis/tmp/"
METAFILE="${OUTDIR}/meta_spikeinfree.tsv"

printf "ID\tANTIBODY\tGROUP\n" > "$METAFILE"
tail -n +2 "$SAMPLESHEET" \
| awk -F'\t' -v OFS='\t' '{print $1".bam", $4, $5}' >> "$METAFILE"


#===========================================================

sudo docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -v  /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  lucidif/chipseqspikeinfree:1.2.4 R

library("ChIPseqSpikeInFree")

setwd("./outs/quantile_normalization_analysis/tmp/spikeinFree_wk")

metaFile <- "../meta_spikeinfree.tsv"

meta.info<-read.table(metaFile, sep="\t", header=TRUE)

#bams <-meta.info$ID
bams <- file.path("..", "inbam", meta.info$ID) 

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "test")

quit()

/media/lucio/easystore/Lucio/Analysis/Lara/quantile_normalization_analysis/tmp/spikeinFree_wk

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/test_SF.txt \
  -b ./outs/quantile_normalization_analysis/tmp/inbam/ \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/ \
  -t 100000000

#=========================================make average

# bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
#   "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/" \
#   "D0_WT_H3K4me3 D0_Mll1-KO_H3K4me3 D0_Mll2-KO_H3K4me3 D0_double-KO_H3K4me3" \
#   "${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw"


#WK="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk"
INBAM="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/"
META="${PWD}/outs/quantile_normalization_analysis/tmp/meta_spikeinfree.tsv"
OUT="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw"

PATTERNS="$(
  awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++) if($i=="ID") c=i; next}
    {print $c}
  ' "$META" \
  | sed -E 's/\.bam$//' \
  | sed -E 's/_rep[AB]$//' \
  | sort -u \
  | tr '\n' ' ' \
  | sed 's/[[:space:]]*$//'
)"

echo "PATTERNS: $PATTERNS"

mkdir -p "$OUT"

bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "$INBAM/" \
  "$PATTERNS" \
  "$OUT"

#compare with previous average

INPATH="${PWD}/outs/geo_sub/out/MLL2_chipSEQ_processed"
OUTPATH="${PWD}/outs/geo_sub/out/MLL2_chipSEQ_processed/average_bw"

PATTERNS=$(find "$INPATH" -maxdepth 1 -type f -name "*.bigWig" -printf "%f\n" \
  | sed -E 's/_rep[AB]\.bigWig$//; s/\.bigWig$//' \
  | sort -u \
  | tr '\n' ' ')

echo "PATTERNS: $PATTERNS"

bash git/chipseq_downstream_macs/bin/averageBigwig.sh \
  "$INPATH" \
  "$PATTERNS" \
  "$OUTPATH"

ls ${PWD}/outs/geo_sub/out/MLL2_chipSEQ_processed/average_bw/*.bw
ls ${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw/*.bw

BW_FREE="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw"
BW_CLASSIC="${PWD}/outs/geo_sub/out/MLL2_chipSEQ_processed/average_bw"
OUTDIR="${PWD}/quantile_normalization_analysis/tmp/compare_norm_methods"

mkdir -p "$OUTDIR"

#2) Costruisci una lista “pairwise” solo dei file presenti in entrambe
comm -12 \
  <(ls -1 "$BW_FREE"/*.bw    | xargs -n1 basename | sort) \
  <(ls -1 "$BW_CLASSIC"/*.bw | xargs -n1 basename | sort) \
  > "$OUTDIR/common_bw_names.txt"


# 3) Crea la lista completa (free + classic) per correlation matrix

awk -v A="$BW_FREE" -v B="$BW_CLASSIC" '{
  print A"/"$0"\n"B"/"$0
}' "$OUTDIR/common_bw_names.txt" > "$OUTDIR/bw_list.txt"


bash git/Lara_MLL2/bin/compare_normalization.sh

#==============file missed

awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || ($2=="MLL2" && $5=="pass" && $7!="NA" && $7!="")' \
  ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/test_SF.txt \
  > ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/test_SF_MLL2_only.txt

#rereun mll2 samples

SAMPLESHEET="sheets/GEO_rename_mll2_chipSEQ_bam.tsv"
OUTDIR="./outs/quantile_normalization_analysis/tmp/"
METAFILE="${OUTDIR}/mll2_meta_spikeinfree.tsv"

printf "ID\tANTIBODY\tGROUP\n" > "$METAFILE"
tail -n +2 "$SAMPLESHEET" \
| awk -F'\t' -v OFS='\t' '{print $1".bam", $4, $5}' >> "$METAFILE"

sudo docker run --rm -it \
  --user "$(id -u)":"$(id -g)" \
  -v  /media/lucio/easystore:/media/lucio/easystore \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  lucidif/chipseqspikeinfree:1.2.4 R

library("ChIPseqSpikeInFree")

setwd("./outs/quantile_normalization_analysis/tmp/spikeinFree_wk")

metaFile <- "../mll2_meta_spikeinfree.tsv"

meta.info<-read.table(metaFile, sep="\t", header=TRUE)

#bams <-meta.info$ID
bams <- file.path("..", "inbam", meta.info$ID) 

ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "mll2",cutoff_QC = 0.70)

quit()

/media/lucio/easystore/Lucio/Analysis/Lara/quantile_normalization_analysis/tmp/spikeinFree_wk

git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
  -s ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2_SF.txt \
  -b ./outs/quantile_normalization_analysis/tmp/inbam/ \
  -g ./outs/Lara_multiomic_analysis/in/mm10.sizes \
  -o ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2only/ \
  -t 1000000000

DIR="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2only/MLL2"
cd "$DIR"

#sudo chown -R "$(id -u)":"$(id -g)" "$DIR"
#chmod -R u+rwX "$DIR"

# D0_WT_MLL2 average
sudo docker run --rm \
  -v "$DIR":"$DIR" \
  -v /media/lucio/easystore:/media/lucio/easystore \
  -w "$DIR" \
  -u "$(id -u)":"$(id -g)" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bash -lc '
    set -euo pipefail
    out="D0_WT_MLL2_average.bw.tmp"
    bigwigAverage -p 8 -b D0_WT_MLL2_repA.spikeinfree.bw D0_WT_MLL2_repB.spikeinfree.bw -o "$out"
    test -s "$out"
    mv "$out" D0_WT_MLL2_average.bw
  '

# D4_WT_MLL2 average
sudo docker run --rm \
  -v "$DIR":"$DIR" \
  -w "$DIR" \
  -u "$(id -u)":"$(id -g)" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bash -lc '
    set -euo pipefail
    out="D4_WT_MLL2_average.bw.tmp"
    bigwigAverage -p 8 -b D4_WT_MLL2_repA.spikeinfree.bw D4_WT_MLL2_repB.spikeinfree.bw -o "$out"
    test -s "$out"
    mv "$out" D4_WT_MLL2_average.bw
  '

ls -lh *_average.bw

DIR_SPIKE="$PWD"  # qui sei in .../mll2only/MLL2 con gli *_average.bw spikeinfree
DIR_CLASSIC="${PWD}/../../../../../geo_sub/out/MLL2_chipSEQ_processed/average_bw"
BED="${PWD}/../../../../../test_chipseq_dowstream/otherouts/deeptools_heatmaps/coordinate.bed"

OUTDIR="$PWD/compare_spike_vs_classic_BED"
mkdir -p "$OUTDIR"

SAMPLES="D0_WT_MLL2 D4_WT_MLL2"

sudo docker run --rm \
  -v "$DIR_SPIKE":"$DIR_SPIKE" \
  -v "$DIR_CLASSIC":"$DIR_CLASSIC" \
  -v "$BED":"$BED" \
  -v "$OUTDIR":"$OUTDIR" \
  -w "$OUTDIR" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bash -lc '
    set -euo pipefail

    echo -e "sample\tbins_used\tpearson\tspearman" > pairwise_spike_vs_classic_BED.tsv
    : > failed_samples.tsv

    DIR_SPIKE="'"$DIR_SPIKE"'"
    DIR_CLASSIC="'"$DIR_CLASSIC"'"
    BED="'"$BED"'"
    SAMPLES="'"$SAMPLES"'"

    rank_avg_ties() {
      # ranks 1..n con average per ties; input: lista di float su stdin, output: lista di rank su stdout
      python3 - "$@" << "PY"
import sys
vals=[float(x) for x in sys.stdin.read().strip().split()]
n=len(vals)
idx=sorted(range(n), key=lambda i: vals[i])
r=[0.0]*n
i=0
rank=1
while i<n:
    j=i
    while j+1<n and vals[idx[j+1]]==vals[idx[i]]:
        j+=1
    # average rank for ties: ranks rank..rank+(j-i)
    avg=(rank + (rank + (j-i)))/2.0
    for k in range(i, j+1):
        r[idx[k]]=avg
    rank += (j-i+1)
    i=j+1
print(" ".join(str(x) for x in r))
PY
    }

    for sample in $SAMPLES; do
      bwA="$DIR_SPIKE/${sample}_average.bw"
      bwB="$DIR_CLASSIC/${sample}_average.bw"

      echo "[INFO] $sample" >&2

      if [[ ! -s "$bwA" ]]; then
        echo -e "${sample}\tmissing_spike\t$bwA" >> failed_samples.tsv
        continue
      fi
      if [[ ! -s "$bwB" ]]; then
        echo -e "${sample}\tmissing_classic\t$bwB" >> failed_samples.tsv
        continue
      fi

      multiBigwigSummary BED-file \
        -b "$bwA" "$bwB" \
        --BED "$BED" \
        -o "${sample}.npz" \
        --outRawCounts "${sample}.counts.tsv" \
        2> "${sample}.multibigwigsummary.log" || {
          echo -e "${sample}\tmultiBigwigSummary_failed\tsee_log" >> failed_samples.tsv
          continue
        }

      # calcolo correlazioni dal counts.tsv senza pandas
      python3 - "$sample" >> pairwise_spike_vs_classic_BED.tsv << "PY"
import sys, math, csv
sample=sys.argv[1]
fn=f"{sample}.counts.tsv"

a=[]; b=[]
with open(fn, newline="") as f:
    r=csv.reader(f, delimiter="\t")
    header=next(r, None)
    for row in r:
        if not row: 
            continue
        # multiBigwigSummary: chr start end <bw1> <bw2>
        try:
            x=float(row[3]); y=float(row[4])
        except Exception:
            continue
        if math.isnan(x) or math.isnan(y):
            continue
        if x==0.0 and y==0.0:
            continue
        a.append(x); b.append(y)

n=len(a)
if n < 10:
    print(f"{sample}\t{n}\tNA\tNA")
    sys.exit(0)

# Pearson
mx=sum(a)/n; my=sum(b)/n
vx=sum((x-mx)**2 for x in a)
vy=sum((y-my)**2 for y in b)
if vx==0.0 or vy==0.0:
    pear="NA"
else:
    cov=sum((a[i]-mx)*(b[i]-my) for i in range(n))
    pear=cov/math.sqrt(vx*vy)

# Spearman = Pearson sui rank (average ties)
# rankdata con ties (metodo semplice)
def rank_avg(vals):
    idx=sorted(range(len(vals)), key=lambda i: vals[i])
    r=[0.0]*len(vals)
    i=0; rank=1
    while i<len(vals):
        j=i
        while j+1<len(vals) and vals[idx[j+1]]==vals[idx[i]]:
            j+=1
        avg=(rank + (rank + (j-i)))/2.0
        for k in range(i, j+1):
            r[idx[k]]=avg
        rank += (j-i+1)
        i=j+1
    return r

ra=rank_avg(a); rb=rank_avg(b)
mrx=sum(ra)/n; mry=sum(rb)/n
vx=sum((x-mrx)**2 for x in ra)
vy=sum((y-mry)**2 for y in rb)
if vx==0.0 or vy==0.0:
    spear="NA"
else:
    cov=sum((ra[i]-mrx)*(rb[i]-mry) for i in range(n))
    spear=cov/math.sqrt(vx*vy)

def fmt(x):
    if isinstance(x,str): return x
    return f"{x:.6f}"

print(f"{sample}\t{n}\t{fmt(pear)}\t{fmt(spear)}")
PY

    done

    echo "[DONE] wrote: pairwise_spike_vs_classic_BED.tsv and failed_samples.tsv" >&2
  '

#DIR_SPIKE="$PWD"  # cartella dove sei ora (mll2only/MLL2) con *_average.bw spikeinfree
DIR_CLASSIC="${PWD}/../../../../../geo_sub/out/MLL2_chipSEQ_processed/average_bw"
BED="${PWD}/../../../../../test_chipseq_dowstream/otherouts/deeptools_heatmaps/coordinate.bed"

OUTDIR="$PWD/compare_spike_vs_classic_BED"
mkdir -p "$OUTDIR"

SAMPLES="D0_WT_MLL2 D4_WT_MLL2"

sudo docker run --rm \
  -v "$DIR_SPIKE":"$DIR_SPIKE" \
  -v "$DIR_CLASSIC":"$DIR_CLASSIC" \
  -v "$BED":"$BED" \
  -v "$OUTDIR":"$OUTDIR" \
  -w "$OUTDIR" \
  quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
  bash -lc '
    set -euo pipefail

    echo -e "sample\tbins_used\tpearson\tspearman" > pairwise_spike_vs_classic_BED.tsv
    : > failed_samples.tsv

    DIR_SPIKE="'"$DIR_SPIKE"'"
    DIR_CLASSIC="'"$DIR_CLASSIC"'"
    BED="'"$BED"'"
    SAMPLES="'"$SAMPLES"'"

    for sample in $SAMPLES; do
      bwA="$DIR_SPIKE/${sample}_average.bw"
      bwB="$DIR_CLASSIC/${sample}_average.bw"

      echo "[INFO] $sample" >&2

      if [[ ! -s "$bwA" ]]; then echo -e "${sample}\tmissing_spike\t$bwA" >> failed_samples.tsv; continue; fi
      if [[ ! -s "$bwB" ]]; then echo -e "${sample}\tmissing_classic\t$bwB" >> failed_samples.tsv; continue; fi

      multiBigwigSummary BED-file \
        -b "$bwA" "$bwB" \
        --BED "$BED" \
        -o "${sample}.npz" \
        --outRawCounts "${sample}.counts.tsv" \
        2> "${sample}.multibigwigsummary.log" || {
          echo -e "${sample}\tmultiBigwigSummary_failed\tsee_log" >> failed_samples.tsv
          continue
        }

      # conta punti utili (no NaN, non entrambi zero)
      n=$(python3 - "$sample" << "PY"
import sys, math, csv
s=sys.argv[1]
fn=f"{s}.counts.tsv"
cnt=0
with open(fn, newline="") as f:
    r=csv.reader(f, delimiter="\t")
    next(r, None)
    for row in r:
        if not row: 
            continue
        try:
            x=float(row[3]); y=float(row[4])
        except:
            continue
        if math.isnan(x) or math.isnan(y): 
            continue
        if x==0.0 and y==0.0:
            continue
        cnt += 1
print(cnt)
PY
)

      if [[ "$n" -lt 10 ]]; then
        echo -e "${sample}\t${n}\tNA\tNA" >> pairwise_spike_vs_classic_BED.tsv
        echo -e "${sample}\ttoo_few_points\t${n}" >> failed_samples.tsv
        continue
      fi

      # scatterplot come compare_normalization.sh
      plotCorrelation \
        --corData "${sample}.npz" \
        --corMethod pearson \
        --whatToPlot scatterplot \
        --labels spikeinfree classic \
        --skipZeros \
        --removeOutliers \
        --plotFile "${sample}_scatter_pearson.pdf" \
        2> "${sample}.plotCorrelation.stderr.log" || {
          echo -e "${sample}\tplotCorrelation_failed\tsee_log" >> failed_samples.tsv
          continue
        }

      # (opzionale) aggiungi anche Spearman scatter
      plotCorrelation \
        --corData "${sample}.npz" \
        --corMethod spearman \
        --whatToPlot scatterplot \
        --labels spikeinfree classic \
        --skipZeros \
        --removeOutliers \
        --plotFile "${sample}_scatter_spearman.pdf" \
        2> "${sample}.plotCorrelation_spearman.stderr.log" || true

      # se vuoi anche heatmap 2x2 (non clustering pesante)
      plotCorrelation \
        --corData "${sample}.npz" \
        --corMethod pearson \
        --whatToPlot heatmap \
        --labels spikeinfree classic \
        --skipZeros \
        --removeOutliers \
        --plotNumbers \
        --plotFile "${sample}_heatmap_pearson.pdf" \
        2> "${sample}.plotCorrelation_heatmap.stderr.log" || true

      # qui puoi anche appendere le correlazioni numeriche (se vuoi)
      # per ora lasciamo solo i PDF + eventuali log
      echo -e "${sample}\t${n}\tNA\tNA" >> pairwise_spike_vs_classic_BED.tsv
    done

    echo "[DONE] output in: '"$OUTDIR"'" >&2
  '

#=============save out files

cd $main_prj_folder

mkdir ./outs/quantile_normalization_analysis/average_bw/noinput
mkdir ./outs/quantile_normalization_analysis/average_bw/scatter
mkdir ./outs/quantile_normalization_analysis/spikeinfree_corrected_bw
mv ./outs/quantile_normalization_analysis/average_bw/*.bw ./outs/quantile_normalization_analysis/average_bw/noinput

cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/*.bw ./outs/quantile_normalization_analysis/spikeinfree_corrected_bw/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/test_* ./outs/quantile_normalization_analysis/spikeinfree_corrected_bw/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2only/MLL2/*.spikeinfree.bw ./outs/quantile_normalization_analysis/spikeinfree_corrected_bw/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2_*.txt ./outs/quantile_normalization_analysis/spikeinfree_corrected_bw/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw/*average.bw ./outs/quantile_normalization_analysis/average_bw/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2only/MLL2/*average.bw ./outs/quantile_normalization_analysis/average_bw/

cp ./outs/quantile_normalization_analysis/tmp/compare_norm_methods/pairwise/* ./outs/quantile_normalization_analysis/average_bw/scatter/
cp ./outs/quantile_normalization_analysis/tmp/spikeinFree_wk/mll2only/MLL2/compare_spike_vs_classic_BED/* ./outs/quantile_normalization_analysis/average_bw/scatter/