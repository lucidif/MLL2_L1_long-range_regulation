#!/usr/bin/env bash
#set -euo pipefail
export LC_ALL=C

# === INPUTS ===
D4="nfout/v2_D4_MLL2_broad_0_5/star/mergedLibrary/macs2/broadPeak/D4WTMLL2B_peaks.broadPeak"
D0="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream/nfout/broad_0_5/star/mergedLibrary/macs2/broadPeak/Anti-GFP.mLb.mkD.sorted_peaks.broadPeak"
OUT="otherouts/D0_vs_D4v2_compare"
mkdir -p "$OUT"

[ -s "$D4" ] || { echo "Manca file D4: $D4"; exit 1; }
[ -s "$D0" ] || { echo "Manca file D0: $D0"; exit 1; }

# === SORT ===
sort -k1,1 -k2,2n "$D4" > "$OUT/D4.sorted.bed"
sort -k1,1 -k2,2n "$D0" > "$OUT/D0.sorted.bed"

ND4=$(wc -l < "$OUT/D4.sorted.bed")
ND0=$(wc -l < "$OUT/D0.sorted.bed")
echo "N_D4=$ND4  N_D0=$ND0" | tee "$OUT/summary.txt"

# === OVERLAP (direzionale e reciproco) ===
printf "type\tf\tD4_in_D0\tD4_in_D0_%%\tD0_in_D4\tD0_in_D4_%%\n" > "$OUT/overlap.tsv"
for f in 0.20 0.50 0.80; do
  D4_in_D0=$(bedtools intersect -a "$OUT/D4.sorted.bed" -b "$OUT/D0.sorted.bed" -f $f -u | wc -l)
  D0_in_D4=$(bedtools intersect -a "$OUT/D0.sorted.bed" -b "$OUT/D4.sorted.bed" -f $f -u | wc -l)
  printf "directional\t%s\t%d\t%.2f\t%d\t%.2f\n" "$f" "$D4_in_D0" "$(awk -v x=$D4_in_D0 -v n=$ND4 'BEGIN{print 100*x/n}')" "$D0_in_D4" "$(awk -v x=$D0_in_D4 -v n=$ND0 'BEGIN{print 100*x/n}')" >> "$OUT/overlap.tsv"

  D4_in_D0_r=$(bedtools intersect -a "$OUT/D4.sorted.bed" -b "$OUT/D0.sorted.bed" -f $f -r -u | wc -l)
  D0_in_D4_r=$(bedtools intersect -a "$OUT/D0.sorted.bed" -b "$OUT/D4.sorted.bed" -f $f -r -u | wc -l)
  printf "reciprocal\t%s\t%d\t%.2f\t%d\t%.2f\n" "$f" "$D4_in_D0_r" "$(awk -v x=$D4_in_D0_r -v n=$ND4 'BEGIN{print 100*x/n}')" "$D0_in_D4_r" "$(awk -v x=$D0_in_D4_r -v n=$ND0 'BEGIN{print 100*x/n}')" >> "$OUT/overlap.tsv"
done
echo "-> $OUT/overlap.tsv"

# === SHARED / UNIQUE (reciproco ≥50%) ===
bedtools intersect -a "$OUT/D4.sorted.bed" -b "$OUT/D0.sorted.bed" -f 0.5 -r -u > "$OUT/D4_shared_rec50.bed"
bedtools intersect -a "$OUT/D0.sorted.bed" -b "$OUT/D4.sorted.bed" -f 0.5 -r -u > "$OUT/D0_shared_rec50.bed"  # stesso numero righe atteso
bedtools intersect -a "$OUT/D4.sorted.bed" -b "$OUT/D0.sorted.bed" -f 0.5 -r -v > "$OUT/D4_unique_rec50.bed"
bedtools intersect -a "$OUT/D0.sorted.bed" -b "$OUT/D4.sorted.bed" -f 0.5 -r -v > "$OUT/D0_unique_rec50.bed"

echo "shared_rec50: $(wc -l < "$OUT/D4_shared_rec50.bed")"
echo "D4_unique_rec50: $(wc -l < "$OUT/D4_unique_rec50.bed")"
echo "D0_unique_rec50: $(wc -l < "$OUT/D0_unique_rec50.bed")" | tee -a "$OUT/summary.txt"

# === JACCARD (bp) ===
cut -f1-3 "$OUT/D4.sorted.bed" > "$OUT/D4.3.bed"
cut -f1-3 "$OUT/D0.sorted.bed" > "$OUT/D0.3.bed"
bedtools jaccard -a "$OUT/D4.3.bed" -b "$OUT/D0.3.bed" | tee "$OUT/jaccard_raw.txt"

# versione "merge" (niente doppi conteggi di bp)
bedtools merge -i "$OUT/D4.3.bed" > "$OUT/D4.3m.bed"
bedtools merge -i "$OUT/D0.3.bed" > "$OUT/D0.3m.bed"
bedtools jaccard -a "$OUT/D4.3m.bed" -b "$OUT/D0.3m.bed" | tee "$OUT/jaccard_merged.txt"

# riassunto bp
bpD4=$(awk '{s+=$3-$2} END{print s}' "$OUT/D4.3m.bed")
bpD0=$(awk '{s+=$3-$2} END{print s}' "$OUT/D0.3m.bed")
inter=$(awk 'NR==2{print $1}' "$OUT/jaccard_merged.txt")
unio=$(awk 'NR==2{print $2}' "$OUT/jaccard_merged.txt")
Aw=$(awk -v i=$inter -v a=$bpD4 'BEGIN{printf "%.2f", 100*i/a}')
Bw=$(awk -v i=$inter -v b=$bpD0 'BEGIN{printf "%.2f", 100*i/b}')
printf "bp_D4\t%d\nbp_D0\t%d\nbp_intersection\t%d\nbp_union\t%d\nD4_shared_bp_%%\t%s\nD0_shared_bp_%%\t%s\n" \
  "$bpD4" "$bpD0" "$inter" "$unio" "$Aw" "$Bw" \
  | tee "$OUT/bp_summary.tsv"

echo "FATTO. Risultati in: $OUT/"
