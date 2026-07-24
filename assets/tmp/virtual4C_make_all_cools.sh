#!/usr/bin/env bash
set -euo pipefail

BASE="outs/Lara_multiomic_analysis/nfout/virtual4C_chr10_Rfx6"
PAIRIX_DIR="$BASE/pairix"
OUTDIR="outs/Lara_multiomic_analysis/outs/test_virtual4C/custom_cool"          # output finale dei cool
CS="$BASE/chromsizes/genome.fa.sizes"

#BASE="outs/Lara_multiomic_analysis/nfout/virtual4C_chr10_Rfx6"
#PAIRIX_DIR="$BASE/pairix"
MERGEDIR="outs/Lara_multiomic_analysis/outs/test_virtual4C/merge_lp_for_cool"

#OUTDIR="$BASE/cool_mergedLP"

TMPROOT="$MERGEDIR/tmp"
CS="$BASE/chromsizes/genome.fa.sizes"

# risoluzioni
BINSIZES=(5000 15000 25000)

mkdir -p "$MERGEDIR" "$OUTDIR" "$TMPROOT"

DOCKER_ARGS=(
  -v "$(pwd)":"$(pwd)"
  -v /media/lucio/easystore:/media/lucio/easystore
  -v /mnt/datawk1:/mnt/datawk1
  -w "$(pwd)"
  -u "$(id -u)":"$(id -g)"
)

# Trova tutti i sample che hanno LP (es: WT_day0_B_LP1 -> WT_day0_B)
mapfile -t SAMPLES_WITH_LP < <(
  ls -1 "$PAIRIX_DIR"/*_LP*.Dd.pairs.gz \
  | sed -E 's#^.*/##' \
  | sed -E 's/\.Dd\.pairs\.gz$//' \
  | sed -E 's/_LP[0-9]+$//' \
  | sort -u
)

if (( ${#SAMPLES_WITH_LP[@]} == 0 )); then
  echo "No *_LP*.Dd.pairs.gz found in $PAIRIX_DIR" >&2
  exit 1
fi

echo "Found ${#SAMPLES_WITH_LP[@]} samples with LP replicates:"
printf '  %s\n' "${SAMPLES_WITH_LP[@]}"
echo

for sample in "${SAMPLES_WITH_LP[@]}"; do
  shopt -s nullglob
  LP_FILES=( "$PAIRIX_DIR/${sample}_LP"*.Dd.pairs.gz )
  shopt -u nullglob

  if (( ${#LP_FILES[@]} < 2 )); then
    echo "[$sample] only ${#LP_FILES[@]} LP file(s) found, skipping merge."
    continue
  fi

  echo "=== [$sample] merging ${#LP_FILES[@]} LP files ==="

  MERGED="$MERGEDIR/${sample}.Dd.mergedLP.pairs.gz"
  SORTED="$MERGEDIR/${sample}.Dd.mergedLP.sorted.pairs.gz"
  TMP="$TMPROOT/$sample"
  mkdir -p "$TMP"

  # merge LP
  sudo docker run --rm "${DOCKER_ARGS[@]}" \
    quay.io/lucidif/microc:0.0.1 \
    bash -lc "pairtools merge -o '$MERGED' $(printf "'%s' " "${LP_FILES[@]}")"

  # sort (importante!)
  sudo docker run --rm "${DOCKER_ARGS[@]}" \
    quay.io/lucidif/microc:0.0.1 \
    bash -lc "pairtools sort --nproc 8 --tmpdir '$TMP' -o '$SORTED' '$MERGED'"

  # cload pairs a tutte le risoluzioni
  for bs in "${BINSIZES[@]}"; do
    COOL="$OUTDIR/${bs}bp_${sample}.cool"
    echo "  cload ${bs} -> $COOL"

    sudo docker run --rm "${DOCKER_ARGS[@]}" \
      quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
      cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        "${CS}:${bs}" "$SORTED" "$COOL"
  done

  echo
done

echo "DONE."
echo "Merged/sorted LP pairs in: $MERGEDIR"
echo "Cool files in: $OUTDIR"
