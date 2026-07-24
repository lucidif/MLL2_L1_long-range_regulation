#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  compare_and_write_md5.sh <TSV> <RAW_DEST_DIR> [PROC_DEST_DIR] <OUT_TSV>

Create and compare MD5 of source files (from TSV) vs destination files (already copied/renamed in DEST dirs),
then write a tab-delimited report in the SAME order as the TSV. Each *file* gets one line:
  new_name<TAB>md5_src<TAB>md5_dst<TAB>status<TAB>row<TAB>kind<TAB>src_path<TAB>dst_path

- KIND is FASTQ1, FASTQ2, or PROCESSED (if applicable).
- STATUS is MATCH, MISMATCH, MISSING_SRC_FILE, MISSING_DST_FILE.
- Works with RAW-ONLY TSVs (no processed columns) and FULL TSVs.

Examples:
  # RAW-only (micro-C)
  bash compare_and_write_md5.sh geo_sub_microC.tsv out/MLL2_microC_raw out/microC_md5_compare.tsv

  # FULL (FASTQ + processed)
  bash compare_and_write_md5.sh GEOrename_RNAseq.tsv out/RNA_raw out/RNA_proc out/RNA_md5_compare.tsv
EOF
}

if (( $# < 3 )); then usage; exit 2; fi

TSV="$1"; shift
RAW_DEST="$1"; shift
if (( $# == 1 )); then
  PROC_DEST=""
  OUT="$1"
else
  PROC_DEST="$1"; shift
  OUT="$1"
fi

# prepare out
: > "$OUT"

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then
    md5 -r "$1" | awk '{print $1}'
  else
    echo "ERROR: neither md5sum nor md5 found in PATH" >&2
    exit 3
  fi
}

# print header
printf "new_name\tmd5_src\tmd5_dst\tstatus\trow\tkind\tsrc_path\tdst_path\n" >> "$OUT"

awk -v FS="\t" -v OFS="\t" -v RAW_DEST="$RAW_DEST" -v PROC_DEST="${PROC_DEST:-}" -v OUT="$OUT" '
BEGIN{}
NR==1{
  for(i=1;i<=NF;i++){ gsub(/\r/,"",$i); gsub(/^[ \t]+|[ \t]+$/, "", $i); h[$i]=i }
  # minimal raw columns
  rawreq="GEOname1 GEOname2 originalName1 originalName2 fastq_path"
  n=split(rawreq,r," "); for(i=1;i<=n;i++){ if(!(r[i] in h)) miss=miss r[i] " " }
  if(length(miss)){ print "ERROR: missing raw columns: " miss > "/dev/stderr"; exit 2 }
  has_proc = ("processed_path" in h) && ("GEOprocessed" in h) && ("originalProcessed" in h)
  next
}
{
  # read fields
  fqdir=$(h["fastq_path"]); on1=$(h["originalName1"]); on2=$(h["originalName2"])
  g1=$(h["GEOname1"]); g2=$(h["GEOname2"])

  # clean CR/spaces
  gsub(/\r/,"",fqdir); gsub(/\r/,"",on1); gsub(/\r/,"",on2); gsub(/\r/,"",g1); gsub(/\r/,"",g2);
  gsub(/^[ \t]+|[ \t]+$/, "", fqdir); gsub(/^[ \t]+|[ \t]+$/, "", on1); gsub(/^[ \t]+|[ \t]+$/, "", on2);
  gsub(/^[ \t]+|[ \t]+$/, "", g1); gsub(/^[ \t]+|[ \t]+$/, "", g2);

  src1 = (length(on1)? fqdir "/" on1 : ""); src2 = (length(on2)? fqdir "/" on2 : "")
  dst1 = (length(g1)? RAW_DEST "/" g1 : "");  dst2 = (length(g2)? RAW_DEST "/" g2 : "")

  proc_src=""; geop=""; dstp=""
  if (has_proc){
    pp=$(h["processed_path"]); op=$(h["originalProcessed"]); geop=$(h["GEOprocessed"])
    gsub(/\r/,"",pp); gsub(/\r/,"",op); gsub(/\r/,"",geop);
    gsub(/^[ \t]+|[ \t]+$/, "", pp); gsub(/^[ \t]+|[ \t]+$/, "", op); gsub(/^[ \t]+|[ \t]+$/, "", geop);
    if (op ~ /^\// || op ~ /\//)      proc_src = op;
    else if (pp ~ /\// && length(op)) proc_src = pp "/" op;
    dstp = (length(geop) && length(PROC_DEST)? PROC_DEST "/" geop : "")
  }

  print NR, g1, g2, src1, src2, dst1, dst2, proc_src, dstp, has_proc
}
' "$TSV" | while IFS=$'\t' read -r ROWNO G1 G2 SRC1 SRC2 DST1 DST2 PROCSRC DSTP HASP; do
  # helper to compare and emit one line
  compare_emit() {
    local NEW="$1" SRC="$2" DST="$3" KIND="$4" ROW="$5"
    local smd5="" dmd5="" status=""
    if [[ -n "$SRC" && -e "$SRC" ]]; then smd5="$(md5_of "$SRC")"; else status="MISSING_SRC_FILE"; fi
    if [[ -n "$DST" && -e "$DST" ]]; then dmd5="$(md5_of "$DST")"; else status="${status:+$status|}MISSING_DST_FILE"; fi
    if [[ -z "$status" ]]; then
      if [[ "$smd5" == "$dmd5" ]]; then status="MATCH"; else status="MISMATCH"; fi
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NEW" "${smd5:-}" "${dmd5:-}" "$status" "$ROW" "$KIND" "${SRC:-}" "${DST:-}" >> "$OUT"
  }

  # FASTQ1
  if [[ -n "$G1" ]]; then compare_emit "$G1" "$SRC1" "$DST1" "FASTQ1" "$ROWNO"; fi
  # FASTQ2
  if [[ -n "$G2" ]]; then compare_emit "$G2" "$SRC2" "$DST2" "FASTQ2" "$ROWNO"; fi
  # PROCESSED
  if [[ "$HASP" == "1" && -n "$DSTP" ]]; then compare_emit "$G1" "$PROCSRC" "$DSTP" "PROCESSED" "$ROWNO"; fi
done

echo "MD5 compare written to: $OUT"
