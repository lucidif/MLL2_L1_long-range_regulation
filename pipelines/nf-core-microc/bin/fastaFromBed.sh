#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bed_to_many_fastas_docker.sh -b regions.bed -f reference.fa -o outdir [-s]

Opzioni:
  -b  BED file (>=3 colonne: chrom start end). 0-based, end esclusivo (standard BED)
  -f  FASTA di riferimento
  -o  directory output (conterrà un fasta per riga)
  -s  rispetta la strand (se nel BED c'è la 6a colonna +/-) => reverse-complement per '-'

Header FASTA:
  >chrom:start-end|name=...|strand=...
  (name e strand compaiono solo se disponibili/richiesti)

Variabili d'ambiente (opzionali):
  DOCKER=...          (default: "sudo docker")
  BEDTOOLS_IMG=...    (default: quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0)
  EXTRA_MOUNT_SRC=... (default: /media/lucio/easystore; se esiste viene montata)
EOF
}

BED=""
FASTA=""
OUTDIR=""
USE_STRAND=0

while getopts ":b:f:o:s" opt; do
  case "$opt" in
    b) BED="$OPTARG" ;;
    f) FASTA="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) USE_STRAND=1 ;;
    *) usage; exit 1 ;;
  esac
done

[[ -n "$BED" && -n "$FASTA" && -n "$OUTDIR" ]] || { usage; exit 1; }
[[ -f "$BED" ]] || { echo "[ERROR] BED non trovato: $BED" >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "[ERROR] FASTA non trovato: $FASTA" >&2; exit 1; }

DOCKER="${DOCKER:-sudo docker}"
BEDTOOLS_IMG="${BEDTOOLS_IMG:-quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0}"
EXTRA_MOUNT_SRC="${EXTRA_MOUNT_SRC:-/media/lucio/easystore}"

mkdir -p "$OUTDIR"

BED_ABS="$(readlink -f "$BED")"
FASTA_ABS="$(readlink -f "$FASTA")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

# Temp files dentro OUTDIR (visibili al container)
bed_with_names="${OUTDIR_ABS}/.tmp_named.$$.bed"
multi_fa="${OUTDIR_ABS}/.tmp_multi.$$.fa"

cleanup() { rm -f "$bed_with_names" "$multi_fa"; }
trap cleanup EXIT

# Costruisci una 4a colonna che contenga SEMPRE chrom:start-end (+ name/strand se presenti)
# Header finale (col4):
#   chr:start-end|name=NAME|strand=+
# oppure:
#   chr:start-end
awk -v use_strand="$USE_STRAND" 'BEGIN{OFS="\t"}
{
  chr=$1; start=$2; end=$3;
  coord = chr ":" start "-" end;

  hdr = coord;

  if ($4 != "") hdr = hdr "|name=" $4;

  if (use_strand==1 && $6 ~ /^[+-]$/) hdr = hdr "|strand=" $6;

  $4 = hdr;
  print $0
}' "$BED_ABS" > "$bed_with_names"

# Mount dirs
declare -a MOUNTS
BED_DIR="$(dirname "$BED_ABS")"
FASTA_DIR="$(dirname "$FASTA_ABS")"
MOUNTS+=(-v "$BED_DIR:$BED_DIR:ro")
if [[ "$FASTA_DIR" != "$BED_DIR" ]]; then
  MOUNTS+=(-v "$FASTA_DIR:$FASTA_DIR:ro")
fi
MOUNTS+=(-v "$OUTDIR_ABS:$OUTDIR_ABS:rw")
if [[ -d "$EXTRA_MOUNT_SRC" ]]; then
  MOUNTS+=(-v "$EXTRA_MOUNT_SRC:$EXTRA_MOUNT_SRC:rw")
fi

# bedtools getfasta -> multi-fasta
if [[ "$USE_STRAND" -eq 1 ]]; then
  GETFASTA_ARGS=(getfasta -fi "$FASTA_ABS" -bed "$bed_with_names" -name -s -fo "$multi_fa")
else
  GETFASTA_ARGS=(getfasta -fi "$FASTA_ABS" -bed "$bed_with_names" -name -fo "$multi_fa")
fi

"$DOCKER" run --rm -u "$(id -u):$(id -g)" \
  "${MOUNTS[@]}" \
  "$BEDTOOLS_IMG" bedtools "${GETFASTA_ARGS[@]}"

# Split multi-fasta -> un file per entry
awk -v outdir="$OUTDIR_ABS" '
function sanitize(s,    t) {
  t=s
  gsub(/^>/,"",t)
  gsub(/[[:space:]]+/,"_",t)
  gsub(/[^A-Za-z0-9._-]/,"_",t)
  return t
}
BEGIN { file="" }
{
  if ($0 ~ /^>/) {
    if (file != "") close(file)
    fname = sanitize($0)
    file = outdir "/" fname ".fa"
  }
  print $0 >> file
}
' "$multi_fa"

n=$(grep -c '^>' "$multi_fa" || true)
echo "[OK] Creati $n FASTA in: $OUTDIR_ABS" >&2
