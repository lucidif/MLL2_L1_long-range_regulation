#!/usr/bin/env bash
set -euo pipefail

# Verifica md5 tra sorgenti e file rinominati.
# Schema TSV richiesto (con header):
# GEOname1  GEOname2  originalName1  originalName2  fastq_path
#
# Modalità:
# 1) Con -m md5.tsv: usa md5 attesi dal file (GEOname \t md5).
# 2) Senza -m: calcola md5 dei file rinominati in -d, confronta, e salva l'output in -o
#    con esatto ordine/formato: prima tutte le R1 (GEOname1), poi tutte le R2 (GEOname2).

usage() {
  cat <<EOF
Uso:
  $0 -i input.tsv -d DEST_DIR [-m md5.tsv] [-o out_md5.tsv] [--check-dest]

Opzioni:
  -i, --input       TSV con header: GEOname1 GEOname2 originalName1 originalName2 fastq_path
  -d, --dest        Directory con i file rinominati (GEOname1/GEOname2)
  -m, --md5         (opzionale) md5.tsv da usare come riferimento (GEOname<TAB>md5)
  -o, --out         (richiesto se manca -m) dove scrivere il nuovo md5.tsv generato
      --check-dest  (opz.) mostra anche l'md5 calcolato sui file in -d nel report dei mismatch

Esito:
  - Non stampa nulla se tutto coincide.
  - Stampa righe "MISSING" o "MISMATCH" se ci sono discrepanze.
EOF
}

INPUT="" ; DEST="" ; MD5FILE="" ; OUT="" ; CHECK_DEST=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2;;
    -d|--dest)  DEST="$2";  shift 2;;
    -m|--md5)   MD5FILE="$2"; shift 2;;
    -o|--out)   OUT="$2"; shift 2;;
    --check-dest) CHECK_DEST=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Argomento sconosciuto: $1" >&2; usage; exit 1;;
  esac
done

[[ -n "$INPUT" && -n "$DEST" ]] || { echo "Errore: -i e -d sono obbligatori." >&2; usage; exit 1; }
[[ -f "$INPUT" ]] || { echo "Errore: input TSV non trovato: $INPUT" >&2; exit 1; }
[[ -d "$DEST"  ]] || { echo "Errore: directory di destinazione non trovata: $DEST" >&2; exit 1; }

MODE=""
if [[ -n "$MD5FILE" ]]; then
  [[ -f "$MD5FILE" ]] || { echo "Errore: file md5.tsv non trovato: $MD5FILE" >&2; exit 1; }
  MODE="USE_MD5_FILE"
else
  [[ -n "$OUT" ]] || { echo "Errore: manca -m e non hai specificato -o per salvare il nuovo md5.tsv." >&2; usage; exit 1; }
  MODE="GENERATE_MD5_FROM_DEST"
fi

# -- parsing header -----------------------------------------------------------
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
ppi=$(idx_of "fastq_path"    || true)

if [[ -z "$g1i" || -z "$g2i" || -z "$o1i" || -z "$o2i" || -z "$ppi" ]]; then
  echo "Errore: header TSV non valido. Trovato: $HEADER" >&2
  echo "Richieste: GEOname1 GEOname2 originalName1 originalName2 fastq_path" >&2
  exit 1
fi

# -- carica o genera md5 attesi (per i file rinominati) -----------------------
declare -A MD5_EXPECT

if [[ "$MODE" == "USE_MD5_FILE" ]]; then
  # Carica da file
  while IFS=$'\t' read -r NAME HASH || [[ -n "${NAME-}" ]]; do
    [[ -z "${NAME-}" ]] && continue
    [[ "${NAME:0:1}" == "#" ]] && continue
    MD5_EXPECT["$NAME"]="$HASH"
  done < "$MD5FILE"

else
  # Genera md5 per i file rinominati in -d e scrive in -o con ordine R1 poi R2
  TMP_R1="$(mktemp)"; TMP_R2="$(mktemp)"
  trap 'rm -f "$TMP_R1" "$TMP_R2"' EXIT

  tail -n +2 "$INPUT" | while IFS=$'\t' read -r line; do
    [[ -z "$line" ]] && continue
    IFS=$'\t' read -r -a F <<< "$line"
    GEO1="${F[$g1i]}"; GEO2="${F[$g2i]}"
    DST1="${DEST%/}/${GEO1}"
    DST2="${DEST%/}/${GEO2}"

    if [[ -f "$DST1" ]]; then
      MD1="$(md5sum "$DST1" | awk '{print $1}')"
      printf "%s\t%s\n" "$GEO1" "$MD1" >> "$TMP_R1"
      MD5_EXPECT["$GEO1"]="$MD1"
    else
      echo "MISSING\tR1\tDEST_NON_TROVATO\t$DST1"
    fi

    if [[ -f "$DST2" ]]; then
      MD2="$(md5sum "$DST2" | awk '{print $1}')"
      printf "%s\t%s\n" "$GEO2" "$MD2" >> "$TMP_R2"
      MD5_EXPECT["$GEO2"]="$MD2"
    else
      echo "MISSING\tR2\tDEST_NON_TROVATO\t$DST2"
    fi
  done

  : > "$OUT"
  cat "$TMP_R1" >> "$OUT"
  cat "$TMP_R2" >> "$OUT"
  echo "Nuovo md5.tsv scritto in: $OUT" >&2
fi

# -- confronto con i sorgenti -------------------------------------------------
tail -n +2 "$INPUT" | while IFS=$'\t' read -r line; do
  [[ -z "$line" ]] && continue
  IFS=$'\t' read -r -a F <<< "$line"

  GEO1="${F[$g1i]}"; GEO2="${F[$g2i]}"
  ORIG1="${F[$o1i]}"; ORIG2="${F[$o2i]}"
  PPATH="${F[$ppi]}"

  SRC1="${PPATH%/}/${ORIG1}"
  SRC2="${PPATH%/}/${ORIG2}"
  DST1="${DEST%/}/${GEO1}"
  DST2="${DEST%/}/${GEO2}"

  SRC1_MD5=""; SRC2_MD5=""
  if [[ -f "$SRC1" ]]; then
    SRC1_MD5="$(md5sum "$SRC1" | awk '{print $1}')"
  else
    echo -e "MISSING\tR1\tSORGENTE_NON_TROVATO\t$SRC1"
  fi
  if [[ -f "$SRC2" ]]; then
    SRC2_MD5="$(md5sum "$SRC2" | awk '{print $1}')"
  else
    echo -e "MISSING\tR2\tSORGENTE_NON_TROVATO\t$SRC2"
  fi

  EXP1="${MD5_EXPECT[$GEO1]:-}"; EXP2="${MD5_EXPECT[$GEO2]:-}"
  if [[ -z "$EXP1" ]]; then
    echo -e "MISSING\tR1\tMD5_ATTESO_NON_PRESENTE\t$GEO1"
  fi
  if [[ -z "$EXP2" ]]; then
    echo -e "MISSING\tR2\tMD5_ATTESO_NON_PRESENTE\t$GEO2"
  fi

  if [[ -n "$SRC1_MD5" && -n "$EXP1" && "$SRC1_MD5" != "$EXP1" ]]; then
    if [[ $CHECK_DEST -eq 1 && -f "$DST1" ]]; then
      DST1_MD5="$(md5sum "$DST1" | awk '{print $1}')"
      echo -e "MISMATCH\tR1\tGEO=$GEO1\tSRC=$SRC1_MD5\tEXP=$EXP1\tDST=$DST1_MD5\tSRC_PATH=$SRC1\tDST_PATH=$DST1"
    else
      echo -e "MISMATCH\tR1\tGEO=$GEO1\tSRC=$SRC1_MD5\tEXP=$EXP1\tSRC_PATH=$SRC1"
    fi
  fi

  if [[ -n "$SRC2_MD5" && -n "$EXP2" && "$SRC2_MD5" != "$EXP2" ]]; then
    if [[ $CHECK_DEST -eq 1 && -f "$DST2" ]]; then
      DST2_MD5="$(md5sum "$DST2" | awk '{print $1}')"
      echo -e "MISMATCH\tR2\tGEO=$GEO2\tSRC=$SRC2_MD5\tEXP=$EXP2\tDST=$DST2_MD5\tSRC_PATH=$SRC2\tDST_PATH=$DST2"
    else
      echo -e "MISMATCH\tR2\tGEO=$GEO2\tSRC=$SRC2_MD5\tEXP=$EXP2\tSRC_PATH=$SRC2"
    fi
  fi
done
