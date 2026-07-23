#!/usr/bin/env bash
shopt -s nullglob

# ============================================================
# v4c_recalc_window_from_pairs_and_mm10bam.sh
#
# Purpose:
#   Recalculate virtual4C at a new WINDOW, reusing:
#     - per-sample *.Dd.pairs.gz in NF_OUT/pairix
#     - a pooled lifted mm10 BAM (already available)
#
# Steps:
#   1) Pool per-sample Dd.pairs.gz -> pooled contig pairs (header preserved)
#   2) Convert pooled contig pairs -> mm10 pairs (if needed)
#   3) Run virtual4c_from_pairs_and_bam.sh with --norm none
#   4) Optional RPMgw bigWig using bamCoverage --scaleFactor from total_nodups.tsv
#
# Design choices:
#   - NO set -e/-u/-o pipefail (debug-friendly)
#   - Minimal hard-fail: only on missing required args/files
# ============================================================

usage() {
  cat <<'EOF'
Usage:
  v4c_recalc_window_from_pairs_and_bam.sh \
    --ref-name <REF> \
    --sam-prefix <SAM> \
    --viewpoint <MM10_POS_INT> \
    --window <BP_INT> \
    --mm10-sizes <mm10.sizes> \
    --nf-out <nfout/ref_folder> \
    --bam-mm10 <pooled.lifted.mm10.bam> \
    --outroot <outs/virtual4C> \
    [--v4c-script <path>] \
    [--docker-env <path>] \
    [--pairtools-img <img>] \
    [--bamcov-img <img>] \
    [--binsize 500] [--smooth 0] [--threads 8] \
    [--skip-pool] [--skip-mm10-pairs] [--skip-rpmgw]

Notes:
  - Expects per-sample pairs at: <nf-out>/pairix/*.Dd.pairs.gz
  - Uses total_nodups.tsv at:     <nf-out>/pairix/total_nodups.tsv (for RPMgw)
  - Output dir:
      <outroot>/bw_tracks_<ref-name>/<window>/<sam-prefix>/

EOF
}

# ----------------------------
# Defaults
# ----------------------------
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
V4C_SCRIPT="git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh"
PAIRTOOLS_IMG="lucidif/microc:0.0.1"
BAMCOV_IMG="quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0"

BINSIZE=500
SMOOTH=0
THREADS=8

SKIP_POOL=0
SKIP_MM10_PAIRS=0
SKIP_RPMGW=0

# ----------------------------
# Required args
# ----------------------------
REF_NAME=""
SAM_PREFIX=""
VIEWPOINT=""
WINDOW=""
MM10_SIZES=""
NF_OUT=""
BAM_MM10=""
OUTROOT=""

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref-name) REF_NAME="$2"; shift 2;;
    --sam-prefix) SAM_PREFIX="$2"; shift 2;;
    --viewpoint) VIEWPOINT="$2"; shift 2;;
    --window) WINDOW="$2"; shift 2;;
    --mm10-sizes) MM10_SIZES="$2"; shift 2;;
    --nf-out) NF_OUT="$2"; shift 2;;
    --bam-mm10) BAM_MM10="$2"; shift 2;;
    --outroot) OUTROOT="$2"; shift 2;;

    --v4c-script) V4C_SCRIPT="$2"; shift 2;;
    --docker-env) DOCKER_ENV="$2"; shift 2;;
    --pairtools-img) PAIRTOOLS_IMG="$2"; shift 2;;
    --bamcov-img) BAMCOV_IMG="$2"; shift 2;;

    --binsize) BINSIZE="$2"; shift 2;;
    --smooth) SMOOTH="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;

    --skip-pool) SKIP_POOL=1; shift;;
    --skip-mm10-pairs) SKIP_MM10_PAIRS=1; shift;;
    --skip-rpmgw) SKIP_RPMGW=1; shift;;

    -h|--help) usage; exit 0;;
    *) echo "[WARN] Unknown arg: $1" >&2; shift;;
  esac
done

# ----------------------------
# Helpers (no hard exits except missing required)
# ----------------------------
err() { echo "[ERROR] $*" >&2; }
warn() { echo "[WARN] $*" >&2; }
info() { echo "[INFO] $*" >&2; }

require_nonempty() {
  local name="$1"
  local val="$2"
  if [[ -z "$val" ]]; then
    err "Missing required arg: $name"
    return 1
  fi
  return 0
}

require_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    err "Missing/empty file: $path"
    return 1
  fi
  return 0
}

require_dir() {
  local path="$1"
  if [[ ! -d "$path" ]]; then
    err "Missing dir: $path"
    return 1
  fi
  return 0
}

# ----------------------------
# Validate required args
# ----------------------------
ok=1
require_nonempty --ref-name "$REF_NAME" || ok=0
require_nonempty --sam-prefix "$SAM_PREFIX" || ok=0
require_nonempty --viewpoint "$VIEWPOINT" || ok=0
require_nonempty --window "$WINDOW" || ok=0
require_nonempty --mm10-sizes "$MM10_SIZES" || ok=0
require_nonempty --nf-out "$NF_OUT" || ok=0
require_nonempty --bam-mm10 "$BAM_MM10" || ok=0
require_nonempty --outroot "$OUTROOT" || ok=0
[[ "$ok" -eq 1 ]] || { usage; exit 1; }

PAIRDIR="${NF_OUT}/pairix"
TOTAL_TSV="${PAIRDIR}/total_nodups.tsv"

OUTDIR="${OUTROOT}/bw_tracks_${REF_NAME}/${WINDOW}/${SAM_PREFIX}"
mkdir -p "$OUTDIR"

POOLED_CONTIG="${PAIRDIR}/${SAM_PREFIX}.pooled.Dd.pairs.gz"
POOLED_MM10="${PAIRDIR}/${SAM_PREFIX}.pooled.Dd.mm10.pairs.gz"

PREFIX="${SAM_PREFIX}.pooled"

BW_NONE="${OUTDIR}/${PREFIX}.viewpoint.filtered.lifted.mm10.none.bin${BINSIZE}.smooth${SMOOTH}.bw"
BW_RPMGW="${OUTDIR}/${PREFIX}.RPMgw.bwbin${BINSIZE}.smooth${SMOOTH}.bw"

info "REF_NAME     = $REF_NAME"
info "SAM_PREFIX   = $SAM_PREFIX"
info "VIEWPOINT    = $VIEWPOINT"
info "WINDOW       = $WINDOW"
info "NF_OUT       = $NF_OUT"
info "PAIRDIR      = $PAIRDIR"
info "BAM_MM10     = $BAM_MM10"
info "OUTDIR       = $OUTDIR"
info "BINSIZE/SMOOTH/THREADS = $BINSIZE/$SMOOTH/$THREADS"
info "SKIP_POOL    = $SKIP_POOL"
info "SKIP_MM10_PAIRS = $SKIP_MM10_PAIRS"
info "SKIP_RPMGW   = $SKIP_RPMGW"

# ----------------------------
# Check core files/dirs
# ----------------------------
require_dir "$PAIRDIR" || exit 1
require_file "$MM10_SIZES" || exit 1
require_file "$BAM_MM10" || exit 1
require_file "$V4C_SCRIPT" || exit 1

# ----------------------------
# Docker env
# ----------------------------
if [[ -f "$DOCKER_ENV" ]]; then
  # shellcheck disable=SC1090
  source "$DOCKER_ENV"
  info "Sourced docker env: $DOCKER_ENV"
else
  warn "docker env not found: $DOCKER_ENV (docker runs may fail if volumes not set)"
  DOCKER_ARGS=()
fi
if [[ "${DOCKER_ARGS+set}" != "set" ]]; then DOCKER_ARGS=(); fi
DOCKER_ENVVARS=( -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache )

# ----------------------------
# 1) Pool contig pairs (optional)
# ----------------------------
if [[ "$SKIP_POOL" -eq 1 ]]; then
  info "Skipping pooling (--skip-pool). Expect pooled contig pairs already present: $POOLED_CONTIG"
else
  if [[ -s "$POOLED_CONTIG" ]]; then
    info "Pooled contig pairs exists: $POOLED_CONTIG (reusing)"
  else
    info "Pooling per-sample pairs -> $POOLED_CONTIG"
    FILES=( "${PAIRDIR}"/${SAM_PREFIX}_*_LP*.Dd.pairs.gz )
    if (( ${#FILES[@]} == 0 )); then
      err "No per-sample pairs matched: ${PAIRDIR}/${SAM_PREFIX}_*_LP*.Dd.pairs.gz"
    else
      HDR="$(mktemp)"
      zcat "${FILES[0]}" | awk '/^#/{print} !/^#/{exit}' > "$HDR"

      sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" \
        "${DOCKER_ENVVARS[@]}" \
        "$PAIRTOOLS_IMG" \
        pairtools merge "${FILES[@]}" \
        | grep -v '^#' \
        | cat "$HDR" - \
        | bgzip -c > "$POOLED_CONTIG"

      rm -f "$HDR"

      if zcat "$POOLED_CONTIG" | grep -m1 '^#chromsize:' >/dev/null 2>&1; then
        info "OK pooled contig pairs: $POOLED_CONTIG"
      else
        warn "Pooled contig pairs seems missing #chromsize (check file): $POOLED_CONTIG"
      fi
    fi
  fi
fi

# ----------------------------
# 2) Convert pooled contig pairs -> mm10 pairs (optional)
# ----------------------------
PAIRS_FOR_V4C="$POOLED_CONTIG"

if [[ "$SKIP_MM10_PAIRS" -eq 1 ]]; then
  info "Skipping mm10 conversion (--skip-mm10-pairs). Using pairs as-is: $PAIRS_FOR_V4C"
else
  # If file already exists, reuse; else convert
  if [[ -s "$POOLED_MM10" ]]; then
    info "Pooled mm10 pairs exists: $POOLED_MM10 (reusing)"
    PAIRS_FOR_V4C="$POOLED_MM10"
  else
    # Heuristic: if chromsize contig contains "chr6:" pattern, treat as contig and convert; else assume already mm10
    if [[ -s "$POOLED_CONTIG" ]]; then
      CONTIG=$(zcat "$POOLED_CONTIG" | awk '/^#chromsize:/{print $2; exit}')
      if echo "$CONTIG" | grep -q 'chr6:[0-9]\+-[0-9]\+'; then
        info "Detected contig-style pairs (#chromsize=$CONTIG) -> converting to mm10"
        REGION=$(echo "$CONTIG" | sed -n 's/.*chr6:\([0-9]\+\)-\([0-9]\+\).*/\1 \2/p')
        START=$(echo "$REGION" | awk '{print $1}')
        END=$(echo "$REGION" | awk '{print $2}')
        OFFSET=$((START-1))

        CHR="chr6"
        CHR_SIZE=$(awk -v c="$CHR" '$1==c{print $2; exit}' "$MM10_SIZES")

        info "START=$START END=$END OFFSET=$OFFSET CHR_SIZE=$CHR_SIZE"

        {
          echo "## pairs format v1.0.0"
          echo "#sorted: chr1-chr2-pos1-pos2"
          echo "#shape: upper triangle"
          echo "#genome_assembly: mm10"
          echo "#chromsize: ${CHR} ${CHR_SIZE}"
          echo "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1 sam2"

          zcat "$POOLED_CONTIG" \
            | awk -F'\t' -v OFS='\t' -v contig="$CONTIG" -v chr="$CHR" -v off="$OFFSET" '
                $0 ~ /^#/ {next}
                {
                  if ($2==contig) {$2=chr; $3=$3+off}
                  if ($4==contig) {$4=chr; $5=$5+off}
                  print
                }'
        } | bgzip -c > "$POOLED_MM10"

        PAIRS_FOR_V4C="$POOLED_MM10"
        info "OK pooled mm10 pairs: $POOLED_MM10"
      else
        warn "Pairs do not look contig-style; using pooled contig as-is: $POOLED_CONTIG"
        PAIRS_FOR_V4C="$POOLED_CONTIG"
      fi
    else
      warn "No pooled contig pairs found to convert: $POOLED_CONTIG"
    fi
  fi
fi

# ----------------------------
# 3) Run v4c script (norm none)
# ----------------------------
info "Running v4c_from_pairs_and_bam (norm none)"
chmod +x "$V4C_SCRIPT" 2>/dev/null || true

"$V4C_SCRIPT" \
  --pairs "$PAIRS_FOR_V4C" \
  --bam "$BAM_MM10" \
  --viewpoint "$VIEWPOINT" \
  --window "$WINDOW" \
  --sizes "$MM10_SIZES" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX" \
  --binsize "$BINSIZE" \
  --smooth "$SMOOTH" \
  --threads "$THREADS" \
  --norm none

if [[ -s "$BW_NONE" ]]; then
  info "OK: BW none produced: $BW_NONE"
else
  warn "BW none not found (script may have failed earlier): $BW_NONE"
fi

# ----------------------------
# 4) RPMgw via bamCoverage scaleFactor (optional)
# ----------------------------
if [[ "$SKIP_RPMGW" -eq 1 ]]; then
  info "Skipping RPMgw (--skip-rpmgw)"
else
  if [[ -s "$BW_RPMGW" ]]; then
    info "RPMgw BigWig exists: $BW_RPMGW (reusing)"
  else
    if [[ ! -s "$TOTAL_TSV" ]]; then
      warn "total_nodups.tsv missing/empty: $TOTAL_TSV -> cannot compute RPMgw"
    else
      N_TOTAL=$(awk 'NR>1 && $1!="TOTAL" && $2!="NA"{s+=$2} END{print s+0}' "$TOTAL_TSV")
      if [[ "$N_TOTAL" -le 0 ]]; then
        warn "Computed N_TOTAL=$N_TOTAL from $TOTAL_TSV -> cannot compute RPMgw"
      else
        SCALE=$(LC_ALL=C awk -v n="$N_TOTAL" 'BEGIN{printf "%.12f", 1000000/n}')
        info "RPMgw scaleFactor=$SCALE (N_TOTAL=$N_TOTAL)"

        BAM_OUT="${OUTDIR}/${PREFIX}.viewpoint.filtered.lifted.mm10.bam"
        if [[ ! -s "$BAM_OUT" ]]; then
          warn "Expected BAM not found: $BAM_OUT -> skipping bamCoverage"
        else
          sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -it \
            "${DOCKER_ENVVARS[@]}" \
            "$BAMCOV_IMG" \
            bamCoverage \
              -b "$BAM_OUT" \
              -o "$BW_RPMGW" \
              --outFileFormat bigwig \
              --binSize "$BINSIZE" \
              --smoothLength "$SMOOTH" \
              --numberOfProcessors "$THREADS" \
              --scaleFactor "$SCALE"

          if [[ -s "$BW_RPMGW" ]]; then
            info "OK RPMgw: $BW_RPMGW"
          else
            warn "RPMgw output missing: $BW_RPMGW"
          fi
        fi
      fi
    fi
  fi
fi

info "DONE"
info "PAIRS_FOR_V4C: $PAIRS_FOR_V4C"
info "OUTDIR:        $OUTDIR"
info "BW_NONE:       $BW_NONE"
info "BW_RPMGW:      $BW_RPMGW"