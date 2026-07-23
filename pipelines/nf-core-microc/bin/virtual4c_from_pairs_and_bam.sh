#!/usr/bin/env bash
#set -euo pipefail

# ============================================================
# virtual4c_from_pairs_and_bam.sh
#
# Generates a per-sample virtual 4C profile from a Micro-C pairs
# file and its corresponding BAM file, both aligned to a
# locus-specific custom reference.
#
# The procedure extracts all read pairs in which at least one end
# maps within a ±<window> bp interval around the viewpoint coordinate
# (given in mm10 coordinates; converted internally to the contig
# coordinate system by parsing the genomic offset from the contig
# name in the pairs header). Reads corresponding to these pairs are
# then retrieved from the BAM, lifted back to mm10 chromosome
# coordinates (contig renamed + POS/PNEXT shifted), and converted
# to a BigWig coverage track with bamCoverage (deepTools).
#
# Steps:
#   0. Retain only concordantly mapped pairs (both pair_type characters
#      != "N") using pairtools select.
#   1. Select pairs with at least one end in the viewpoint window
#      [viewpoint - window, viewpoint + window] (contig coordinates)
#      using pairtools select.
#   2. Extract the QNAME list from the selected pairs and filter the
#      BAM to retain only matching primary/non-supplementary alignments
#      (samtools view -N -F 0x900).
#   3. Lift the filtered BAM to mm10: rename the contig in the SAM
#      @SQ header to the target chromosome, and add the genomic offset
#      (offset - 1) to all POS and PNEXT fields.
#   4. Generate a BigWig coverage track from the lifted BAM with
#      bamCoverage (--binSize 500; normalization via --scaleFactor or
#      --normalizeUsing, as specified by the caller).
# ============================================================


usage() {
  cat <<'EOF'
Usage:
  virtual4c_from_pairs_and_bam.sh \
    --pairs <in.pairs.gz> \
    --bam <in.bam> \
    --viewpoint <mm10_coord_int> \
    --window <bp_int> \
    --sizes <mm10.sizes> \
    --outdir <outdir> \
    --prefix <name> \
    [--binsize 500] [--smooth 5000] [--threads 8] \
    [--norm none|CPM|RPGC] [--effectiveGenomeSize 2150570000]

Outputs (in outdir):
  <prefix>.mappedOnly.pairs.gz
  <prefix>.inRegion.pairs.gz
  <prefix>.qnames.keep.txt
  <prefix>.viewpoint.filtered.bam (+ .bai)
  <prefix>.viewpoint.filtered.lifted.mm10.bam (+ .bai)
  <prefix>.viewpoint.filtered.lifted.mm10.<norm>.bin<B>.smooth<S>.bw

Notes:
- Requires: bgzip, zcat, samtools on host
- Requires: docker + DOCKER_ARGS from git/Lara_MLL2/bin/docker_env.sh
EOF
}

# ----------------------------
# Defaults
# ----------------------------
BINSIZE=500
SMOOTH=5000
THREADS=8
NORM="none"
EFFECTIVE_GENOME_SIZE=""

# ----------------------------
# Parse args
# ----------------------------
PAIRS=""
BAM=""
VIEWPOINT=""
WINDOW=""
SIZES=""
OUTDIR=""
PREFIX=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pairs) PAIRS="$2"; shift 2;;
    --bam) BAM="$2"; shift 2;;
    --viewpoint) VIEWPOINT="$2"; shift 2;;
    --window) WINDOW="$2"; shift 2;;
    --sizes) SIZES="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --binsize) BINSIZE="$2"; shift 2;;
    --smooth) SMOOTH="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --norm) NORM="$2"; shift 2;;
    --effectiveGenomeSize) EFFECTIVE_GENOME_SIZE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

# ----------------------------
# Validate
# ----------------------------
[[ -n "$PAIRS" && -f "$PAIRS" ]] || { echo "[ERROR] --pairs missing or not found"; exit 1; }
[[ -n "$BAM" && -f "$BAM" ]] || { echo "[ERROR] --bam missing or not found"; exit 1; }
[[ -n "$VIEWPOINT" ]] || { echo "[ERROR] --viewpoint missing"; exit 1; }
[[ -n "$WINDOW" ]] || { echo "[ERROR] --window missing"; exit 1; }
[[ -n "$SIZES" && -f "$SIZES" ]] || { echo "[ERROR] --sizes missing or not found"; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "[ERROR] --outdir missing"; exit 1; }
[[ -n "$PREFIX" ]] || { echo "[ERROR] --prefix missing"; exit 1; }

if [[ "$NORM" == "RPGC" && -z "$EFFECTIVE_GENOME_SIZE" ]]; then
  echo "[ERROR] --norm RPGC requires --effectiveGenomeSize" >&2
  #exit 1
fi
if [[ "$NORM" != "none" && "$NORM" != "CPM" && "$NORM" != "RPGC" ]]; then
  echo "[ERROR] --norm must be one of: none, CPM, RPGC" >&2
  #exit 1
fi

mkdir -p "$OUTDIR"

# ----------------------------
# Load docker volumes/args
# ----------------------------
source "git/Lara_MLL2/bin/docker_env.sh"

# ----------------------------
# Output paths
# ----------------------------
PAIRS_MAPPED_ONLY="${OUTDIR}/${PREFIX}.mappedOnly.pairs.gz"
PAIRS_INREGION="${OUTDIR}/${PREFIX}.inRegion.pairs.gz"
QKEEP="${OUTDIR}/${PREFIX}.qnames.keep.txt"
BAM_FILTERED="${OUTDIR}/${PREFIX}.viewpoint.filtered.bam"
BAM_LIFTED="${OUTDIR}/${PREFIX}.viewpoint.filtered.lifted.mm10.bam"
BW_OUT="${OUTDIR}/${PREFIX}.viewpoint.filtered.lifted.mm10.${NORM}.bin${BINSIZE}.smooth${SMOOTH}.bw"

# ----------------------------
# 0) pairs: both ends mapped
# ----------------------------
TMP_P="${PAIRS_MAPPED_ONLY%.gz}"
sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -i \
  -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
  lucidif/microc:0.0.1 \
  pairtools select 'pair_type[0] != "N" and pair_type[1] != "N"' "$PAIRS" \
  > "$TMP_P"

bgzip -c "$TMP_P" > "$PAIRS_MAPPED_ONLY"
rm -f "$TMP_P"

zcat "$PAIRS_MAPPED_ONLY" | awk '!/^#/ {n++} END{print "[INFO] pairs (both mapped):", n+0}' >&2

# ----------------------------
# 1) pairs: at least one end in VIEWPOINT±WINDOW
#    compute relative window from contig header
# ----------------------------
CONTIG=$(zcat -f "$PAIRS_MAPPED_ONLY" | awk '/^#chromsize:/ {print $2; exit}')
[[ -n "$CONTIG" ]] || { echo "[ERROR] No #chromsize found in $PAIRS_MAPPED_ONLY" >&2; exit 1; }

TARGET_CHR=$(echo "$CONTIG" | sed -E 's/.*\|([^:]+):.*/\1/')
OFFSET=$(echo "$CONTIG" | sed -E 's/.*\|[^:]+:([0-9]+)-.*/\1/')
[[ -n "$TARGET_CHR" && -n "$OFFSET" ]] || { echo "[ERROR] Could not parse TARGET_CHR/OFFSET from CONTIG=$CONTIG" >&2; exit 1; }

VP_REL=$(( VIEWPOINT - OFFSET ))
REL_START=$(( VP_REL - WINDOW ))
REL_END=$(( VP_REL + WINDOW + 1 ))
if (( REL_START < 0 )); then REL_START=0; fi

echo "[INFO] CONTIG=$CONTIG" >&2
echo "[INFO] TARGET_CHR=$TARGET_CHR" >&2
echo "[INFO] OFFSET=$OFFSET" >&2
echo "[INFO] VIEWPOINT(mm10)=$VIEWPOINT -> rel=$VP_REL" >&2
echo "[INFO] rel interval: $REL_START-$REL_END" >&2

sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -i \
  -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
  lucidif/microc:0.0.1 \
  pairtools select \
    "((chrom1==\"$CONTIG\" and pos1>=$REL_START and pos1<$REL_END) or (chrom2==\"$CONTIG\" and pos2>=$REL_START and pos2<$REL_END))" \
    "$PAIRS_MAPPED_ONLY" \
  | bgzip -c > "$PAIRS_INREGION"

zcat "$PAIRS_INREGION" | awk '!/^#/ {n++} END{print "[INFO] pairs in window:", n+0}' >&2

# ----------------------------
# 2) BAM: filter by readID list from pairs subset
# ----------------------------
zcat -f "$PAIRS_INREGION" | awk '!/^#/ {print $1}' | sort -u > "$QKEEP"
echo "[INFO] QNAME kept: $(wc -l < "$QKEEP")" >&2

samtools view -b -N "$QKEEP" -F 0x900 "$BAM" \
| samtools sort -o "$BAM_FILTERED"
samtools index "$BAM_FILTERED"
echo "[INFO] reads in filtered BAM: $(samtools view -c "$BAM_FILTERED")" >&2

# ----------------------------
# 3) Lift BAM to mm10 (rename contig + shift POS/PNEXT) + fix header SQ
# ----------------------------
CHR_LN=$(awk -v c="$TARGET_CHR" '$1==c{print $2}' "$SIZES")
[[ -n "$CHR_LN" ]] || { echo "[ERROR] $TARGET_CHR not found in $SIZES" >&2; exit 1; }

SHIFT=$(( OFFSET - 1 ))
echo "[INFO] SHIFT(for SAM POS/PNEXT)=$SHIFT" >&2
echo "[INFO] $TARGET_CHR LN=$CHR_LN" >&2

HDR_LIFTED=$(mktemp)

samtools view -H "$BAM_FILTERED" \
| awk -v C="$CONTIG" -v T="$TARGET_CHR" -v ln="$CHR_LN" 'BEGIN{OFS="\t"}
    $1=="@SQ" {
      sn=""
      for(i=2;i<=NF;i++) if($i ~ /^SN:/) sn=substr($i,4)
      if(sn==C){ print "@SQ","SN:"T,"LN:"ln; next }
    }
    {print}
  ' > "$HDR_LIFTED"

samtools view "$BAM_FILTERED" \
| awk -v C="$CONTIG" -v T="$TARGET_CHR" -v sh="$SHIFT" 'BEGIN{OFS="\t"}
    {
      orig=$3
      if($3==C){$3=T; $4+=sh}
      if($7=="="){ if(orig==C){$8+=sh} }
      else if($7==C){ $7=T; $8+=sh }
      print
    }' \
| cat "$HDR_LIFTED" - \
| samtools view -bS - \
| samtools sort -o "$BAM_LIFTED"

samtools index "$BAM_LIFTED"
rm -f "$HDR_LIFTED"

# Quick checks
echo "[CHECK] first alignment (QNAME RNAME POS RNEXT PNEXT):" >&2
samtools view "$BAM_LIFTED" | head -n 1 | cut -f1,3,4,7,8 >&2

# ----------------------------
# 4) BigWig
# ----------------------------
DT_ARGS=(--outFileFormat bigwig --binSize "$BINSIZE" --smoothLength "$SMOOTH" --numberOfProcessors "$THREADS")

if [[ "$NORM" == "CPM" ]]; then
  DT_ARGS+=(--normalizeUsing CPM)
elif [[ "$NORM" == "RPGC" ]]; then
  DT_ARGS+=(--normalizeUsing RPGC --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE")
fi

sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -i \
  -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
  quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0 \
  bamCoverage -b "$BAM_LIFTED" -o "$BW_OUT" "${DT_ARGS[@]}"

echo "[OK] BigWig: $BW_OUT" >&2
ls -lh "$BW_OUT" >&2
