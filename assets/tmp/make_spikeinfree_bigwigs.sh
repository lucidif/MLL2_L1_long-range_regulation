#!/usr/bin/env bash
set -euo pipefail

# =========================
# Inputs:
#   -s  *_SF.txt
#   -b  directory dove cercare i BAM (anche ricorsivo)
#   -g  chrom.sizes
#   -o  output directory
#   -t  target reads (default 15000000)
# =========================

TARGET=15000000
SF_TSV=""
BAM_DIR=""
CHROMSIZES=""
OUTDIR=""

SAMTOOLS_IMG="quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
BEDTOOLS_IMG="quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
BW_IMG="quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0"

EXTRA_MOUNT_SRC="/media/lucio/easystore"
DOCKER="${DOCKER:-sudo docker}"

die(){ echo "[ERROR] $*" >&2; exit 1; }
log(){ echo "[INFO] $*" >&2; }

while getopts ":s:b:g:o:t:" opt; do
  case "$opt" in
    s) SF_TSV="$OPTARG" ;;
    b) BAM_DIR="$OPTARG" ;;
    g) CHROMSIZES="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) TARGET="$OPTARG" ;;
    \?) die "Opzione non valida: -$OPTARG" ;;
    :)  die "Opzione -$OPTARG richiede un argomento" ;;
  esac
done

[[ -f "$SF_TSV" ]] || die "File SF non trovato: $SF_TSV"
[[ -d "$BAM_DIR" ]] || die "Directory BAM non trovata: $BAM_DIR"
[[ -f "$CHROMSIZES" ]] || die "chrom.sizes non trovato: $CHROMSIZES"
[[ -n "$OUTDIR" ]] || die "Specifica -o OUTDIR"
mkdir -p "$OUTDIR"

# docker wrapper (stile tuo)
MOUNTS=(-v "$PWD":"$PWD")
if [[ -d "$EXTRA_MOUNT_SRC" ]]; then
  MOUNTS+=(-v "$EXTRA_MOUNT_SRC":"$EXTRA_MOUNT_SRC")
fi

dkr() {
  local image="$1"; shift
  $DOCKER run --rm \
    -u "$(id -u)":"$(id -g)" \
    "${MOUNTS[@]}" \
    -w "$PWD" \
    "$image" \
    "$@"
}

# indicizza BAM: filename -> full path
log "Indicizzo BAM in: $BAM_DIR"
MAPFILE="$(mktemp)"
trap 'rm -f "$MAPFILE"' EXIT
find "$BAM_DIR" -type f -name "*.bam" -printf "%f\t%p\n" > "$MAPFILE"

find_bam() {
  local id="$1"
  if [[ -f "$id" ]]; then echo "$id"; return 0; fi
  local base; base="$(basename "$id")"
  awk -F'\t' -v b="$base" '$1==b{print $2; exit}' "$MAPFILE"
}

# colonne nel file SF
read -r HEADER < "$SF_TSV"
IFS=$'\t' read -r -a HCOLS <<< "$HEADER"
col_idx() {
  local name="$1"
  for i in "${!HCOLS[@]}"; do
    [[ "${HCOLS[$i]}" == "$name" ]] && { echo $((i+1)); return 0; }
  done
  return 1
}
CID="$(col_idx ID)"
CAB="$(col_idx ANTIBODY)"
CQC="$(col_idx QC)"
CSF="$(col_idx SF)"

log "TARGET=$TARGET | immagini: samtools=$SAMTOOLS_IMG bedtools=$BEDTOOLS_IMG bw=$BW_IMG"

tail -n +2 "$SF_TSV" | while IFS= read -r line; do
  [[ -z "$line" ]] && continue

  id="$(awk -F'\t' -v i="$CID" '{print $i}' <<< "$line")"
  ab="$(awk -F'\t' -v i="$CAB" '{print $i}' <<< "$line")"
  qc="$(awk -F'\t' -v i="$CQC" '{print $i}' <<< "$line")"
  sf="$(awk -F'\t' -v i="$CSF" '{print $i}' <<< "$line")"

  [[ "$qc" == "pass" ]] || continue
  [[ -n "$sf" && "$sf" != "NA" ]] || continue

  bam_path="$(find_bam "$id" || true)"
  if [[ -z "${bam_path:-}" ]]; then
    log "SKIP: BAM non trovato per ID=$id (ANTIBODY=$ab)"
    continue
  fi

  sample="$(basename "$id")"
  sample="${sample%.bam}"

  out_ab_dir="$OUTDIR/$ab"
  mkdir -p "$out_ab_dir"

  bdg_tmp="$out_ab_dir/${sample}.spikeinfree.tmp.bedGraph"
  bdg="$out_ab_dir/${sample}.spikeinfree.bedGraph"
  bw="$out_ab_dir/${sample}.spikeinfree.bw"

  # paired-end: conta reads "properly paired" (primari) e dividi per 2 => frammenti
  reads_pp="$(dkr "$SAMTOOLS_IMG" samtools view -c -f 2 -F 260 "$bam_path")"
  libSize="$(awk -v n="$reads_pp" 'BEGIN{print n/2.0}')"

  scale="$(awk -v t="$TARGET" -v l="$libSize" -v s="$sf" 'BEGIN{print t/(l*s)}')"

  log "ANTIBODY=$ab sample=$sample SF=$sf libSize(frags)=$libSize scale=$scale"
  dkr "$BEDTOOLS_IMG" bedtools genomecov -bg -pc -scale "$scale" -ibam "$bam_path" -g "$CHROMSIZES" > "$bdg_tmp"
  LC_ALL=C sort -k1,1 -k2,2n "$bdg_tmp" > "$bdg"
  rm -f "$bdg_tmp"

  dkr "$BW_IMG" bedGraphToBigWig "$bdg" "$CHROMSIZES" "$bw"
done

log "Fatto."
