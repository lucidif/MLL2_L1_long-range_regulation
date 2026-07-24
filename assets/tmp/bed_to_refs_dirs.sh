#!/usr/bin/env bash
#set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bed_to_refs_dirs.sh -b regions.bed -f reference.fa -o outdir [-s] [--index]

Cosa fa:
  - Estrae le sequenze dal reference.fa per ogni riga del BED
  - Crea una sottocartella per ogni entry (nome = col4 del BED)
  - Scrive OUTDIR/<name>/genome.fa con header che include chr:start-end
  - Se --index: crea bwa index (genome.fa.*) e samtools faidx (genome.fa.fai)

Opzioni:
  -b  BED (>=3 colonne). Se vuoi cartelle con nome stabile, metti la col4 (name).
  -f  FASTA di riferimento
  -o  output directory
  -s  rispetta la strand (usa col6 +/- e reverse-complement per '-')
  --index  indicizza (bwa index + samtools faidx)

Docker:
  - bedtools è eseguito via Docker.
  - bwa/samtools: se non sono nel PATH, puoi fornire immagini con:
      BWA_IMG="..." SAMTOOLS_IMG="..."

Env:
  DOCKER        default: "sudo docker"
  BEDTOOLS_IMG  default: quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0
  EXTRA_MOUNT_SRC default: /media/lucio/easystore  (se esiste, viene montata)

Note:
  - BED è 0-based start, end esclusivo.
  - Header scritto in genome.fa:
      >NAME|chr:start-end   (e se -s: |strand=+/-)
EOF
}

# ---- parse args
BED=""
FASTA=""
OUTDIR=""
USE_STRAND=0
DO_INDEX=0

# gestiamo anche flag lungo --index
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) BED="$2"; shift 2;;
    -f) FASTA="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -s) USE_STRAND=1; shift 1;;
    --index) DO_INDEX=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Argomento sconosciuto: $1" >&2; usage; exit 1;;
  esac
done

[[ -n "$BED" && -n "$FASTA" && -n "$OUTDIR" ]] || { usage; exit 1; }
[[ -f "$BED" ]] || { echo "[ERROR] BED non trovato: $BED" >&2; exit 1; }
[[ -f "$FASTA" ]] || { echo "[ERROR] FASTA non trovato: $FASTA" >&2; exit 1; }
# ---- normalizza BED (rimuove CRLF) in una copia temporanea, NON modifica l'originale
# BED_CLEAN="${OUTDIR}/.bed_clean.$$.bed"
# sed 's/\r$//' "$BED" > "$BED_CLEAN"
# trap 'rm -f "$tmp_bed" "$tmp_multi" "$BED_CLEAN"' EXIT
# crea OUTDIR subito (serve per scrivere i tmp)
#mkdir -p "$OUTDIR"

# ---- normalizza BED (rimuove CRLF) in una copia temporanea, NON modifica l'originale
#BED_CLEAN="${OUTDIR}/.bed_clean.$$.bed"
BED_CLEAN="$(dirname "$BED")/.bed_clean.$$.bed"
sed 's/\r$//' "$BED" > "$BED_CLEAN"


#DOCKER="${DOCKER:-sudo docker}"

# ---- docker command (supporta anche "sudo docker")
DOCKER_STR="${DOCKER:-}"
if [[ -z "$DOCKER_STR" ]]; then
  # prova docker senza sudo; se non va, ripiega su sudo docker
  if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    DOCKER_CMD=(docker)
  elif command -v sudo >/dev/null 2>&1; then
    DOCKER_CMD=(sudo docker)
  else
    echo "[ERROR] docker non disponibile (né docker né sudo trovati)" >&2
    exit 127
  fi
else
  # se DOCKER="sudo docker" viene splittato in array correttamente
  read -r -a DOCKER_CMD <<<"$DOCKER_STR"
fi


BEDTOOLS_IMG="${BEDTOOLS_IMG:-quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0}"
EXTRA_MOUNT_SRC="${EXTRA_MOUNT_SRC:-/media/lucio/easystore}"

mkdir -p "$OUTDIR"

#BED_ABS="$(readlink -f "$BED")"
BED_ABS="$(readlink -f "$BED_CLEAN")"
FASTA_ABS="$(readlink -f "$FASTA")"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

tmp_bed="${OUTDIR_ABS}/.tmp_named.$$.bed"
tmp_multi="${OUTDIR_ABS}/.tmp_multi.$$.fa"
#cleanup() { rm -f "$tmp_bed" "$tmp_multi"; }
#trap cleanup EXIT
cleanup() { rm -f "$tmp_bed" "$tmp_multi" "$BED_CLEAN"; }
trap cleanup EXIT


# ---- helper: run bwa/samtools (host oppure docker se BWA_IMG / SAMTOOLS_IMG sono settate)
run_bwa() {
  BWA_IMG="${BWA_IMG:-quay.io/biocontainers/bwa:0.7.17--hed695b0_7}"

  "${DOCKER_CMD[@]}" run --rm -u "$(id -u):$(id -g)" \
    -v "$OUTDIR_ABS:$OUTDIR_ABS:rw" \
    -w "$PWD" \
    "$BWA_IMG" bwa "$@"
}

# ---- helper: samtools SEMPRE via Docker (mai dal PATH)
run_samtools() {
  SAMTOOLS_IMG="${SAMTOOLS_IMG:-quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1}"

  "${DOCKER_CMD[@]}" run --rm -u "$(id -u):$(id -g)" \
    -v "$OUTDIR_ABS:$OUTDIR_ABS:rw" \
    -w "$PWD" \
    "$SAMTOOLS_IMG" samtools "$@"
}


# ---- prepara BED con col4 obbligatoria e header informativo
# Se il BED non ha col4, ne genera una: chr_start_end
awk -v use_strand="$USE_STRAND" 'BEGIN{OFS="\t"}
{
  chr=$1; start=$2; end=$3;
  name = ($4 != "" ? $4 : chr"_"start"_"end);
  gsub(/\r/, "", name);              # <-- FIX CRLF
  hdr = name "|" chr ":" start "-" end;
  if (use_strand==1 && $6 ~ /^[+-]$/) hdr = hdr "|strand=" $6;
  $4 = hdr;
  print $0
}' "$BED_ABS" > "$tmp_bed"

# ---- mount per bedtools docker
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

# ---- estrai multi-fasta
if [[ "$USE_STRAND" -eq 1 ]]; then
  "${DOCKER_CMD[@]}" run --rm -u "$(id -u):$(id -g)" \
    "${MOUNTS[@]}" \
    "$BEDTOOLS_IMG" \
    bedtools getfasta -fi "$FASTA_ABS" -bed "$tmp_bed" -name -s -fo "$tmp_multi"
else
  "${DOCKER_CMD[@]}" run --rm -u "$(id -u):$(id -g)" \
    "${MOUNTS[@]}" \
    "$BEDTOOLS_IMG" \
    bedtools getfasta -fi "$FASTA_ABS" -bed "$tmp_bed" -name -fo "$tmp_multi"
fi

# ---- split: OUTDIR/<NAME>/genome.fa
# NAME = parte prima del primo '|'
awk -v outdir="$OUTDIR_ABS" '
function sanitize(s,    t) {
  t=s
  gsub(/[[:space:]]+/,"_",t)
  gsub(/[^A-Za-z0-9._-]/,"_",t)
  return t
}
BEGIN { file="" ; name="" }
{
  if ($0 ~ /^>/) {
    hdr = substr($0,2)
    split(hdr, a, /\|/)
    name = sanitize(a[1])

    dir = outdir "/" name
    # mkdir -p
    cmd = "mkdir -p \"" dir "\""
    system(cmd)

    file = dir "/genome.fa"
    print ">" hdr > file
  } else {
    print $0 >> file
  }
}
' "$tmp_multi"

echo "[OK] Creati genome.fa in sottocartelle sotto: $OUTDIR_ABS" >&2

# ---- indicizzazione (opzionale)
if [[ "$DO_INDEX" -eq 1 ]]; then
  shopt -s nullglob
  for d in "$OUTDIR_ABS"/*; do
    [[ -f "$d/genome.fa" ]] || continue
    (
      cd "$d"
      echo "[INFO] Indexing: $d/genome.fa" >&2
      run_bwa index genome.fa || true
      run_samtools faidx genome.fa || true
    )
  done
  echo "[OK] Indicizzazione completata (dove disponibile)." >&2
fi
