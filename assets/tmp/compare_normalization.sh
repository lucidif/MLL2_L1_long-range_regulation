#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DIR_A="${PWD}/outs/quantile_normalization_analysis/tmp/spikeinFree_wk/average_bw"   # spikeinfree
DIR_B="${PWD}/outs/geo_sub/out/MLL2_chipSEQ_processed/average_bw"                  # classico

OUTDIR="${PWD}/outs/quantile_normalization_analysis/tmp/compare_norm_methods/pairwise"
mkdir -p "$OUTDIR"

# lista dei sample comuni (basename senza path)
comm -12 \
  <(ls -1 "$DIR_A"/*_average.bw 2>/dev/null | xargs -n1 basename | sort) \
  <(ls -1 "$DIR_B"/*_average.bw 2>/dev/null | xargs -n1 basename | sort) \
  > "$OUTDIR/common_samples.txt"

echo -e "sample\tbins_used\tpearson\tspearman" > "$OUTDIR/pairwise_correlations.tsv"

# qui scegli binSize: più grande => meno NaN/zeri; 50kb spesso è più stabile
BINSIZE=50000

FAIL="$OUTDIR/failed_samples.tsv"
: > "$FAIL"

while read -r bw; do
  sample="${bw%_average.bw}"
  bwA="$DIR_A/${sample}_average.bw"
  bwB="$DIR_B/${sample}_average.bw"

  echo "[INFO] $sample"

  # sanity
  if [[ ! -s "$bwA" || ! -s "$bwB" ]]; then
    echo -e "${sample}\tmissing_input_bw" >> "$FAIL"
    continue
  fi

  # 1) summary (se fallisce, continua)
  if ! multiBigwigSummary BED-file \
      -b "$bwA" "$bwB" \
      --BED ./outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/coordinate.bed \
      -o "$OUTDIR/${sample}.npz" \
      --outRawCounts "$OUTDIR/${sample}.counts.tsv" ; then
    echo -e "${sample}\tmultiBigwigSummary_failed" >> "$FAIL"
    continue
  fi

  # 2) check: quanti punti “validi” abbiamo? (no NaN e non entrambi 0)
  #    uso le ultime 2 colonne (robusto anche se BED-file aggiunge colonne davanti)
  n_ok=$(awk 'NR>1{
      a=$(NF-1); b=$NF;
      if(a!="nan" && b!="nan" && !(a==0 && b==0)) c++
    } END{print c+0}' "$OUTDIR/${sample}.counts.tsv")

  if (( n_ok < 10 )); then
    echo -e "${sample}\ttoo_few_points\t${n_ok}" >> "$FAIL"
    continue
  fi

  # 3) scatter (se fallisce, continua)
  if ! plotCorrelation \
      --corData "$OUTDIR/${sample}.npz" \
      --corMethod pearson \
      --whatToPlot scatterplot \
      --skipZeros \
      --removeOutliers \
      --labels spikeinfree classic \
      --plotFile "$OUTDIR/${sample}_scatter_pearson.pdf" ; then
    echo -e "${sample}\tplotCorrelation_failed" >> "$FAIL"
    continue
  fi

done < "$OUTDIR/common_samples.txt"


echo "[DONE] risultati: $OUTDIR/pairwise_correlations.tsv"
