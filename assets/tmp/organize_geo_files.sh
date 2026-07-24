#!/usr/bin/env bash

usage() {
  cat <<'EOF'
Usage:
  organize_geo_files.sh <RAW_DEST_DIR> <PROC_DEST_DIR> [TSV] [MD5_OUT]

Copia/rinomina i file secondo il TSV e scrive un report MD5 a 6 colonne (DEST):
  name1  md5_1  name2  md5_2  nameP  md5_P
Se il TSV non ha le colonne processed, l’ultima coppia resta vuota.

Env:
  DRY_RUN=1    simula (non copia)
  OVERWRITE=1  sovrascrive file già presenti
  ONLY_RAW=1   copia solo FASTQ
  ONLY_PROC=1  copia solo processed (ignorato se TSV RAW-only)
  RSYNC=1      usa rsync -a (default cp -p)
EOF
}

(( $#<2 || $#>4 )) && { usage; exit 2; }
RAW_DEST="$1"; shift
PROC_DEST="$1"; shift
TSV="${1:-geo_sub.tsv}"; shift || true
MD5_OUT="${1:-md5_report.tsv}"

mkdir -p "$RAW_DEST" "$PROC_DEST"
: > "$MD5_OUT"
# HEADER del report (sempre 6 colonne)
printf "name1\tmd5_1\tname2\tmd5_2\tnameP\tmd5_P\n" >> "$MD5_OUT"

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
  [[ -z "$src" || -z "$dst" ]] && return 0
  if [[ ! -e "$src" ]]; then
    echo "WARNING: missing source: $src" >&2
    return 1
  fi
  if [[ -e "$dst" && "${OVERWRITE:-0}" != "1" ]]; then
    echo "INFO: exists, skipping: $dst" >&2
    return 0
  fi
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    echo "DRY-RUN: ${RSYNC:+rsync -a}${RSYNC:-cp -p} \"$src\" \"$dst\"" >&2
  else
    if [[ "${RSYNC:-0}" == "1" ]]; then rsync -a "$src" "$dst"; else cp -p "$src" "$dst"; fi
  fi
}

# Estraggo e normalizzo i campi dal TSV: GEO1 GEO2 GEOP SRC1 SRC2 PROCSRC HAS_PROC
awk -v FS="\t" -v OFS="\t" '
NR==1{
  for(i=1;i<=NF;i++){ gsub(/\r/,"",$i); gsub(/^[ \t]+|[ \t]+$/, "", $i); h[$i]=i }
  need="GEOname1 GEOname2 originalName1 originalName2 fastq_path"
  n=split(need,r," ")
  for(i=1;i<=n;i++) if(!(r[i] in h)) miss=miss r[i] " "
  if(length(miss)){ print "ERROR: missing raw columns: " miss > "/dev/stderr"; exit 2 }
  has_proc=("processed_path" in h)&&("GEOprocessed" in h)&&("originalProcessed" in h)
  next
}
{
  fq=$(h["fastq_path"]); on1=$(h["originalName1"]); on2=$(h["originalName2"])
  g1=$(h["GEOname1"]);  g2=$(h["GEOname2"])
  gsub(/\r/,"",fq); gsub(/\r/,"",on1); gsub(/\r/,"",on2); gsub(/\r/,"",g1); gsub(/\r/,"",g2)
  gsub(/^[ \t]+|[ \t]+$/, "", fq); gsub(/^[ \t]+|[ \t]+$/, "", on1); gsub(/^[ \t]+|[ \t]+$/, "", on2)
  gsub(/^[ \t]+|[ \t]+$/, "", g1); gsub(/^[ \t]+|[ \t]+$/, "", g2)

  src1=(length(on1)? fq "/" on1 : ""); src2=(length(on2)? fq "/" on2 : "")

  procsrc=""; geop=""
  if (has_proc){
    pp=$(h["processed_path"]); op=$(h["originalProcessed"]); geop=$(h["GEOprocessed"])
    gsub(/\r/,"",pp); gsub(/\r/,"",op); gsub(/\r/,"",geop)
    gsub(/^[ \t]+|[ \t]+$/, "", pp); gsub(/^[ \t]+|[ \t]+$/, "", op); gsub(/^[ \t]+|[ \t]+$/, "", geop)
    if (op ~ /^\// || op ~ /\//)      procsrc = op;
    else if (pp ~ /\// && length(op)) procsrc = pp "/" op;
  }
  # OUTPUT ORDER (fissato): GEO1 GEO2 GEOP SRC1 SRC2 PROCSRC HAS_PROC
  print g1, g2, geop, src1, src2, procsrc, (has_proc?1:0)
}
' "$TSV" | while IFS=$'\t' read -r GEO1 GEO2 GEOP SRC1 SRC2 PROCSRC HASP; do
  # Sanifica nomi → basename senza CR/spazi
  for var in GEO1 GEO2 GEOP; do
    val="${!var}"
    val="${val%$'\r'}"
    val="$(printf "%s" "$val" | sed 's/^[[:space:]]\+//; s/[[:space:]]\+$//')"
    [[ "$val" == */* ]] && val="$(basename -- "$val")"
    printf -v "$var" "%s" "$val"
  done

  DST1=""; DST2=""; DSTP=""
  [[ -n "$GEO1" ]] && DST1="$RAW_DEST/$GEO1"
  [[ -n "$GEO2" ]] && DST2="$RAW_DEST/$GEO2"
  [[ -n "$GEOP" && "${HASP}" = "1" ]] && DSTP="$PROC_DEST/$GEOP"

  # Copia
  if [[ "${ONLY_PROC:-0}" != "1" ]]; then
    copy_file "$SRC1" "$DST1" || true
    copy_file "$SRC2" "$DST2" || true
  fi
  if [[ "${ONLY_RAW:-0}" != "1" && -n "$DSTP" ]]; then
    copy_file "$PROCSRC" "$DSTP" || true
  fi

  # MD5 dei DEST
  m1=""; m2=""; mp=""
  if [[ "${DRY_RUN:-0}" != "1" ]]; then
    [[ -n "$DST1" && -e "$DST1" ]] && m1="$(md5_of "$DST1")" || true
    [[ -n "$DST2" && -e "$DST2" ]] && m2="$(md5_of "$DST2")" || true
    [[ -n "$DSTP" && -e "$DSTP" ]] && mp="$(md5_of "$DSTP")" || true
  fi

  # Riga a 6 colonne (mai array sbagliato)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$GEO1" "$m1" "$GEO2" "$m2" "$GEOP" "$mp" >> "$MD5_OUT"
done

echo "Done."
echo "  Raw files dest:  $RAW_DEST"
echo "  Proc files dest: $PROC_DEST"
echo "  MD5 report:      $MD5_OUT"