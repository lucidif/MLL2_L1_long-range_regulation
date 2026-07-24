#!/usr/bin/env bash
set -euo pipefail

# Estrae sequenze definite in un BED da un FASTA di riferimento usando BEDTOOLS via Docker
# e crea tanti file .fa quanti sono gli intervalli del BED.
#
# Requisiti: Docker (o sudo docker). Non serve bedtools installato localmente.

usage() {
  cat <<'EOF'
Usage:
  bed_to_many_fastas_docker.sh -b regions.bed -f reference.fa -o outdir [-s]

Opzioni:
  -b  BED file (>=3 colonne: chrom start end). 0-based, end esclusivo (standard BED)
  -f  FASTA di riferimento
  -o  directory output (conterrà un fasta per riga)
  -s  rispetta la strand (se nel BED c'è la 6a colonna +/-) => reverse-complement per '-'

Variabili d'ambiente (opzionali):
  DOCKER=...          (default: "sudo docker")
  BEDTOOLS_IMG=...    (default: quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0)
  EXTRA_MOUNT_SRC=... (default: /media/lucio/easystore; se esiste viene montata)

Esempi:
  ./bed_to_many_fastas_docker.sh -b regions.bed -f hg38.fa -o fastas
  ./bed_to_many_fastas_docker.sh -b regions6col.bed -f hg38.fa -o fastas -s
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

#DOCKER="${DOCKER:-sudo docker}"
DOCKER_STR="${DOCKER:-sudo docker}"
read -r -a DOCKER_CMD <<< "$DOCKER_STR"

BEDTOOLS_IMG="${BEDTOOLS_IMG:-quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0}"
EXTRA_MOUNT_SRC="${EXTRA_MOUNT_SRC:-/media/lucio/easystore}"

mkdir -p "$OUTDIR"

# Absolutize paths (serve per montare correttamente in docker)
BED_ABS="$(readlink -f "$BED")"
FASTA_ABS="$(readlink -f "$FASTA")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

# Temp files dentro OUTDIR (così sono sempre visibili al container)
bed_with_names="${OUTDIR_ABS}/.tmp_named.$$.bed"
multi_fa="${OUTDIR_ABS}/.tmp_multi.$$.fa"

cleanup() {
  rm -f "$bed_with_names" "$multi_fa"
}
trap cleanup EXIT

# Prepara BED con una colonna 4 (name) “sicura” per -name:
# - se c'è già la 4a colonna, la usa
# - altrimenti usa chr_start_end
# - se -s e c'è la 6a colonna +/- aggiunge _+ o _-
awk -v use_strand="$USE_STRAND" 'BEGIN{OFS="\t"}
{
  chr=$1; start=$2; end=$3;
  name = ($4 != "" ? $4 : chr"_"start"_"end);
  if (use_strand==1 && $6 ~ /^[+-]$/) name = name"_"$6;
  $4 = name;
  print $0
}' "$BED_ABS" > "$bed_with_names"

# Costruisci opzioni di mount:
# - montiamo le directory di BED e FASTA in-place (stesso path dentro container)
# - OUTDIR montata rw
declare -a MOUNTS
BED_DIR="$(dirname "$BED_ABS")"
FASTA_DIR="$(dirname "$FASTA_ABS")"

MOUNTS+=(-v "$BED_DIR:$BED_DIR:ro")
if [[ "$FASTA_DIR" != "$BED_DIR" ]]; then
  MOUNTS+=(-v "$FASTA_DIR:$FASTA_DIR:ro")
fi
MOUNTS+=(-v "$OUTDIR_ABS:$OUTDIR_ABS:rw")

# Mount extra (se esiste) come nel tuo progetto
if [[ -d "$EXTRA_MOUNT_SRC" ]]; then
  MOUNTS+=(-v "$EXTRA_MOUNT_SRC:$EXTRA_MOUNT_SRC:rw")
fi

# Comando bedtools getfasta dentro container
if [[ "$USE_STRAND" -eq 1 ]]; then
  GETFASTA_ARGS=(getfasta -fi "$FASTA_ABS" -bed "$bed_with_names" -name -s -fo "$multi_fa")
else
  GETFASTA_ARGS=(getfasta -fi "$FASTA_ABS" -bed "$bed_with_names" -name -fo "$multi_fa")
fi

#"$DOCKER" run --rm -u "$(id -u):$(id -g)" \
"${DOCKER_CMD[@]}" run --rm -u "$(id -u):$(id -g)" \
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
