#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Uso: $0 -i input.tsv -d DEST_DIR -o md5_output.tsv
TSV atteso con header: GEOname1  GEOname2  originalName1  originalName2  fastq_path
- Copia fastq_path/originalName1 -> DEST/GEOname1
- Copia fastq_path/originalName2 -> DEST/GEOname2
- Output: due colonne "<nome>\t<md5>", prima tutte le R1 (GEOname1), poi tutte le R2 (GEOname2)
EOF
}

INPUT="" ; DEST="" ; OUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2;;
    -d|--dest)  DEST="$2";  shift 2;;
    -o|--out)   OUT="$2";   shift 2;;
    -h|--help)  usage; exit 0;;
    *) echo "Argomento sconosciuto: $1" >&2; usage; exit 1;;
  esac
done

if [[ -z "${INPUT}" || -z "${DEST}" || -z "${OUT}" ]]; then
  echo "Errore: -i, -d e -o sono obbligatori." >&2; usage; exit 1
fi
[[ -f "$INPUT" ]] || { echo "Errore: file di input non trovato: $INPUT" >&2; exit 1; }
mkdir -p "$DEST"

# Leggi header e trova gli indici delle colonne richieste
HEADER="$(head -n1 "$INPUT")"
IFS=$'\t' read -r -a HCOLS <<< "$HEADER"

idx_of () {
  local name="$1"
  for i in "${!HCOLS[@]}"; do
    [[ "${HCOLS[$i]}" == "$name" ]] && { echo "$i"; return 0; }
  done
  return 1
}

g1i=$(idx_of "GEOname1" || true)
g2i=$(idx_of "GEOname2" || true)
o1i=$(idx_of "originalName1" || true)
o2i=$(idx_of "originalName2" || true)
ppi=$(idx_of "fastq_path" || true)

if [[ -z "$g1i" || -z "$g2i" || -z "$o1i" || -z "$o2i" || -z "$ppi" ]]; then
  echo "Errore: impossibile trovare colonne richieste nell'header." >&2
  echo "Header: $HEADER" >&2
  echo "Richieste: GEOname1 GEOname2 originalName1 originalName2 fastq_path" >&2
  exit 1
fi

# File temporanei per rispettare l'ordine: prima R1 (GEOname1), poi R2 (GEOname2)
TMP_R1="$(mktemp)"; TMP_R2="$(mktemp)"
trap 'rm -f "$TMP_R1" "$TMP_R2"' EXIT

# Scorri le righe (salta l'header)
tail -n +2 "$INPUT" | while IFS=$'\t' read -r line; do
  # split sicuro su TAB preservando eventuali spazi nei campi
  IFS=$'\t' read -r -a F <<< "$line"

  GEO1="${F[$g1i]}" ; GEO2="${F[$g2i]}"
  ORIG1="${F[$o1i]}" ; ORIG2="${F[$o2i]}"
  PPATH="${F[$ppi]}"

  # R1
  SRC1="${PPATH%/}/${ORIG1}"
  if [[ ! -f "$SRC1" ]]; then
    echo "ATTENZIONE: sorgente non trovata (R1): $SRC1" >&2
  else
    DEST1="${DEST%/}/${GEO1}"
    cp -f -- "$SRC1" "$DEST1"
    MD1="$(md5sum "$DEST1" | awk '{print $1}')"
    printf "%s\t%s\n" "$GEO1" "$MD1" >> "$TMP_R1"
  fi

  # R2
  SRC2="${PPATH%/}/${ORIG2}"
  if [[ ! -f "$SRC2" ]]; then
    echo "ATTENZIONE: sorgente non trovata (R2): $SRC2" >&2
  else
    DEST2="${DEST%/}/${GEO2}"
    cp -f -- "$SRC2" "$DEST2"
    MD2="$(md5sum "$DEST2" | awk '{print $1}')"
    printf "%s\t%s\n" "$GEO2" "$MD2" >> "$TMP_R2"
  fi
done

# Output finale: prima tutte le R1 poi tutte le R2
: > "$OUT"
cat "$TMP_R1" >> "$OUT"
cat "$TMP_R2" >> "$OUT"

echo "Fatto."
echo "File copiati in: $DEST"
echo "MD5 scritto in : $OUT"
