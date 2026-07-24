#!/usr/bin/env bash
set -euo pipefail

rawsource_folder="/home/lucio/wkdir/data/Lara/RNAseq"   # senza '/' finale va bene uguale
outfolder="out/MLL2_RNAseq_raw"

# ---- liste ----
read -r -a NAMES <<< "D4WT_REP1 D4WT_REP2 D4WT_REP3 D4MLL1KO_REP1 D4MLL1KO_REP2 D4MLL1KO_REP3 D4MLL2KO_REP1 D4MLL2KO_REP2 D4MLL2KO_REP3 D4DoubleKO_REP1 D4DoubleKO_REP2 D4DoubleKO_REP3 D0DoubleKO_REP1 D0DoubleKO_REP2 D0DoubleKO_REP3 D0Mll1KO_REP1 D0Mll1KO_REP2 D0Mll1KO_REP3 D0Mll2KO_REP1 D0Mll2KO_REP2 D0Mll2KO_REP3 D0WTA_REP1 D0WTA_REP2 D0WTA_REP3"

read -r -a FQ1 <<< "rnaseq_lara_june23/1_1.fastq.gz rnaseq_lara_june23/2_1.fastq.gz rnaseq_lara_june23/3_1.fastq.gz rnaseq_lara_june23/4_1.fastq.gz rnaseq_lara_june23/5_1.fastq.gz rnaseq_lara_june23/6_1.fastq.gz rnaseq_lara_june23/7_1.fastq.gz rnaseq_lara_june23/8_1.fastq.gz rnaseq_lara_june23/9_1.fastq.gz rnaseq_lara_june23/10_1.fastq.gz rnaseq_lara_june23/11_1.fastq.gz rnaseq_lara_june23/12_1.fastq.gz rnaseq_lara_fastq_0423/DoubleKOA_1.fastq.gz rnaseq_lara_fastq_0423/DoubleKOB_1.fastq.gz rnaseq_lara_fastq_0423/DoubleKOC_1.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOA_1.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOB_1.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOC_1.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOA_1.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOB_1.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOC_1.fastq.gz rnaseq_lara_fastq_0423/WTA_1.fastq.gz rnaseq_lara_fastq_0423/WTB_1.fastq.gz rnaseq_lara_fastq_0423/WTC_1.fastq.gz"

read -r -a FQ2 <<< "rnaseq_lara_june23/1_2.fastq.gz rnaseq_lara_june23/2_2.fastq.gz rnaseq_lara_june23/3_2.fastq.gz rnaseq_lara_june23/4_2.fastq.gz rnaseq_lara_june23/5_2.fastq.gz rnaseq_lara_june23/6_2.fastq.gz rnaseq_lara_june23/7_2.fastq.gz rnaseq_lara_june23/8_2.fastq.gz rnaseq_lara_june23/9_2.fastq.gz rnaseq_lara_june23/10_2.fastq.gz rnaseq_lara_june23/11_2.fastq.gz rnaseq_lara_june23/12_2.fastq.gz rnaseq_lara_fastq_0423/DoubleKOA_2.fastq.gz rnaseq_lara_fastq_0423/DoubleKOB_2.fastq.gz rnaseq_lara_fastq_0423/DoubleKOC_2.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOA_2.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOB_2.fastq.gz rnaseq_lara_fastq_0423/Mll1_KOC_2.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOA_2.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOB_2.fastq.gz rnaseq_lara_fastq_0423/Mll2_KOC_2.fastq.gz rnaseq_lara_fastq_0423/WTA_2.fastq.gz rnaseq_lara_fastq_0423/WTB_2.fastq.gz rnaseq_lara_fastq_0423/WTC_2.fastq.gz"

# ---- controlli ----
if (( ${#NAMES[@]} != ${#FQ1[@]} || ${#NAMES[@]} != ${#FQ2[@]} )); then
  echo "ERRORE: dimensioni diverse: NAMES=${#NAMES[@]} FQ1=${#FQ1[@]} FQ2=${#FQ2[@]}" >&2
  exit 1
fi

mkdir -p "$outfolder"

ts=$(date +"%Y%m%d-%H%M%S")
checks_file="$outfolder/checksums_$ts.txt"      # MD5 dei DEST (formato per md5sum -c)
report_file="$outfolder/verify_report_$ts.txt"  # confronto SRC vs DEST

echo "# MD5 dei file copiati (destinazione). Generato: $(date)" > "$checks_file"
echo -e "# Verifica copia (src vs dest) — $(date)\n" > "$report_file"

for i in "${!NAMES[@]}"; do
  name="${NAMES[$i]}"

  src1="$rawsource_folder/${FQ1[$i]}"
  src2="$rawsource_folder/${FQ2[$i]}"

  dst1="$outfolder/${name}_R1.fastq.gz"
  dst2="$outfolder/${name}_R2.fastq.gz"

  # Copia (usa -n per NON sovrascrivere se già esiste; togli -n se vuoi forzare)
  cp -v -n "$src1" "$dst1"
  cp -v -n "$src2" "$dst2"

  # Calcola e confronta MD5; scrive checksum dei DEST nel file checksums
  for mate in R1 R2; do
    if [[ "$mate" == "R1" ]]; then src="$src1"; dst="$dst1"; else src="$src2"; dst="$dst2"; fi

    src_md5=$(md5sum "$src" | awk '{print $1}')
    dst_md5=$(md5sum "$dst" | awk '{print $1}')

    printf "%s  %s\n" "$dst_md5" "$dst" >> "$checks_file"

    if [[ "$src_md5" == "$dst_md5" ]]; then
      echo "$name/$mate: OK  $dst" >> "$report_file"
    else
      echo "$name/$mate: **FAILED** MD5 mismatch  src=$src_md5  dst=$dst_md5  ($dst)" >> "$report_file"
    fi
  done
done

echo
echo "Fatto."
echo "MD5 (dest): $checks_file   -> verifica in futuro con:  md5sum -c \"$checks_file\""
echo "Report:     $report_file"
