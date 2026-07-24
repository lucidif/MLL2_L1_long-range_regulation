#!/usr/bin/env bash
# easy_make_md5.sh
# Uso:
#   bash easy_make_md5.sh -t <tsv_path> -o <out_dir> -n <outname> [-p]

set -euo pipefail

usage() {
  echo "Uso: $0 -t <tsv_path> -o <out_dir> -n <outname> [-p]"
  echo "  -t  Percorso al file TSV (obbligatorio)"
  echo "  -o  Cartella output (obbligatoria)"
  echo "  -n  Prefisso nome file output (obbligatorio)"
  echo "  -p  Abilita il loop 'processed' (opzionale)"
  exit 1
}

TSV=""
OUT=""
OUTNAME=""
PROCESSED="FALSE"

# ---- parsing argomenti (FIX) ----
while getopts ":t:o:n:p" opt; do
  case "$opt" in
    t) TSV="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    n) OUTNAME="$OPTARG" ;;
    p) PROCESSED="TRUE" ;;
    *) usage ;;
  esac
done

# param obbligatori
[[ -z "$TSV" || -z "$OUT" || -z "$OUTNAME" ]] && usage
[[ ! -f "$TSV" ]] && { echo "Errore: TSV non trovato: $TSV" >&2; exit 2; }

mkdir -p "$OUT"
rm -f "${OUT}/${OUTNAME}_1.txt" "${OUT}/${OUTNAME}_2.txt" "${OUT}/${OUTNAME}_processed.txt"

# conta righe (anche senza newline finale)
tot_righe=$(awk 'END{print NR}' "$TSV")
(( tot_righe < 2 )) && { echo "Nessuna riga dati (solo header)."; exit 0; }

# --- loop 1: col5=path, co

# conta righe (anche senza newline finale)
tot_righe=$(awk 'END{print NR}' "$TSV")
if (( tot_righe < 2 )); then
  echo "Nessuna riga dati (solo header o file vuoto)."
  exit 0
fi

echo "TSV: $TSV"
echo "OUT: $OUT"
echo "OUTNAME: $OUTNAME"
echo "PROCESSED: $PROCESSED"
echo "Righe totali: $tot_righe"

# --- loop 1: col5=path, col3=filename -> OUTNAME_1.txt ---
echo "Scrivo: ${OUT}/${OUTNAME}_1.txt"
for (( i=2; i<=tot_righe; i++ )); do
  riga=$(head -n "$i" "$TSV" | tail -n 1)
  path=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f5)
  name1=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f3)
  if [[ -f "$path/$name1" ]]; then
    md5sum -- "$path/$name1" >> "${OUT}/${OUTNAME}_1.txt"
  else
    echo "MISSING  $path/$name1" >> "${OUT}/${OUTNAME}_1.txt"
  fi
done

# --- loop 2: col5=path, col4=filename -> OUTNAME_2.txt ---
echo "Scrivo: ${OUT}/${OUTNAME}_2.txt"
for (( i=2; i<=tot_righe; i++ )); do
  riga=$(head -n "$i" "$TSV" | tail -n 1)
  path=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f5)
  name2=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f4)
  if [[ -f "$path/$name2" ]]; then
    md5sum -- "$path/$name2" >> "${OUT}/${OUTNAME}_2.txt"
  else
    echo "MISSING  $path/$name2" >> "${OUT}/${OUTNAME}_2.txt"
  fi
done

# --- loop processed (opzionale): col6=path2, col8=filename -> OUTNAME_processed.txt ---
if [[ "${PROCESSED}" == "TRUE" ]]; then
  echo "Scrivo: ${OUT}/${OUTNAME}_processed.txt"
  for (( i=2; i<=tot_righe; i++ )); do
    riga=$(head -n "$i" "$TSV" | tail -n 1)
    path2=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f6)
    name3=$(printf '%s' "$riga" | tr -d '\r' | cut -d $'\t' -f8)
    if [[ -f "$path2/$name3" ]]; then
      md5sum -- "$path2/$name3" >> "${OUT}/${OUTNAME}_processed.txt"
    else
      echo "MISSING  $path2/$name3" >> "${OUT}/${OUTNAME}_processed.txt"
    fi
  done
fi

echo "Fatto."

