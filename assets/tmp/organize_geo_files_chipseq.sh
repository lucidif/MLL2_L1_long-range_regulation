#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  organize_geo_files.sh <RAW_DEST_DIR> <PROC_DEST_DIR> [TSV] [MD5_OUT]

Arguments:
  RAW_DEST_DIR   Destination directory where FASTQ files will be copied & renamed.
  PROC_DEST_DIR  Destination directory where processed (bigWig) files will be copied & renamed.
  TSV            Input TSV (default: geo_sub_chipSEQ.tsv). Must include columns:
                 GEOname1, GEOname2, originalName1, originalName2, fastq_path,
                 processed_path, GEOprocessed, originalProcessed
  MD5_OUT        Path to write a tab-delimited MD5 report (default: md5_report.tsv).
                 Columns: md5_GEOname1, md5_GEOname2, md5_GEOprocessed (blank if not applicable).

Options (env vars):
  DRY_RUN=1      Show what would be done without copying or hashing.
  OVERWRITE=1    Overwrite destination files if they already exist (default: skip existing).
  ONLY_RAW=1     Process only FASTQ files.
  ONLY_PROC=1    Process only processed files (bigWig).

Examples:
  bash organize_geo_files.sh /data/geo/raw /data/geo/proc
  bash organize_geo_files.sh /data/geo/raw /data/geo/proc geo_sub_chipSEQ.tsv out/md5.tsv
  DRY_RUN=1 bash organize_geo_files.sh raw/ proc/
EOF
}

if (( $# < 2 || $# > 4 )); then usage; exit 2; fi

RAW_DEST="$1"; shift
PROC_DEST="$1"; shift
TSV="${1:-geo_sub_chipSEQ.tsv}"; shift || true
MD5_OUT="${1:-md5_report.tsv}"

mkdir -p "$RAW_DEST" "$PROC_DEST"
: > "$MD5_OUT"

# prefer md5sum, fallback to md5 -r (macOS)
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

copy_file() {
  local src="$1" dst="$2"
  if [[ -z "$src" || -z "$dst" ]]; then return 0; fi
  if [[ ! -e "$src" ]]; then
    echo "WARNING: missing source: $src" >&2
    return 1
  fi
  if [[ -e "$dst" && "${OVERWRITE:-0}" != "1" ]]; then
    echo "INFO: exists, skipping: $dst" >&2
    return 0
  fi
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "DRY-RUN: cp -p \"$src\" \"$dst\"" >&2
  else
    cp -p "$src" "$dst"
  fi
}

# Build a normalized stream of rows using awk (to avoid fragile field indexes)
awk -v FS="\t" -v OFS="\t" '
NR==1{
  for(i=1;i<=NF;i++){ h[$i]=i }
  req="GEOname1 GEOname2 originalName1 originalName2 fastq_path processed_path GEOprocessed originalProcessed"
  n=split(req,r," ")
  for(i=1;i<=n;i++){ if(!(r[i] in h)) miss=miss r[i] " " }
  if(length(miss)){ print "ERROR: missing required columns: " miss > "/dev/stderr"; exit 2 }
  next
}
{
  fq1 = $(h["fastq_path"]) "/" $(h["originalName1"])
  fq2 = $(h["fastq_path"]) "/" $(h["originalName2"])
  bw  = $(h["processed_path"]) "/" $(h["originalProcessed"])
  gsub(/^[ \t]+|[ \t]+$/, "", fq1); gsub(/^[ \t]+|[ \t]+$/, "", fq2); gsub(/^[ \t]+|[ \t]+$/, "", bw)
  print $(h["GEOname1"]),$(h["GEOname2"]),$(h["GEOprocessed"]), fq1, fq2, bw
}
' "$TSV" | while IFS=$'\t' read -r GEO1 GEO2 GEOBW SRC1 SRC2 SRCBW; do
  # Compute destination paths
  DST1="" DST2="" DSTBW=""
  if [[ -n "$GEO1" ]]; then DST1="$RAW_DEST/$GEO1"; fi
  if [[ -n "$GEO2" ]]; then DST2="$RAW_DEST/$GEO2"; fi
  if [[ -n "$GEOBW" ]]; then DSTBW="$PROC_DEST/$GEOBW"; fi

  # Copy files as requested
  if [[ "${ONLY_PROC:-0}" != "1" ]]; then
    copy_file "$SRC1" "$DST1" || true
    copy_file "$SRC2" "$DST2" || true
  fi
  if [[ "${ONLY_RAW:-0}" != "1" ]]; then
    # copy processed only if GEOBW and SRCBW are not empty
    if [[ -n "$GEOBW" && -n "$SRCBW" ]]; then
      copy_file "$SRCBW" "$DSTBW" || true
    fi
  fi

  # MD5 in TSV order (three columns per row; blanks if not applicable)
  md5_1=""; md5_2=""; md5_bw=""
  if [[ "${DRY_RUN:-0}" != "1" ]]; then
    if [[ -n "$DST1" && -e "$DST1" ]]; then md5_1="$(md5_of "$DST1")"; fi
    if [[ -n "$DST2" && -e "$DST2" ]]; then md5_2="$(md5_of "$DST2")"; fi
    if [[ -n "$DSTBW" && -e "$DSTBW" ]]; then md5_bw="$(md5_of "$DSTBW")"; fi
  fi
  printf "%s\t%s\t%s\n" "$md5_1" "$md5_2" "$md5_bw" >> "$MD5_OUT"
done

echo "Done."
echo "  Raw files dest:  $RAW_DEST"
echo "  Proc files dest: $PROC_DEST"
echo "  MD5 report:      $MD5_OUT"
