#!/usr/bin/env bash
set -eu

# ============================================================
# virtual4c_reads_counting.sh
#
# Counts the number of valid Micro-C pairs with at least one end
# mapping within a ±<window> bp interval around a specified viewpoint
# coordinate, for all samples in a pairix output directory.
#
# When called with a viewpoint and window (WINDOW mode), this script
# reproduces exactly the pair-selection logic of step 1 in
# virtual4c_from_pairs_and_bam.sh, ensuring that the counts are
# consistent with the pairs actually used to build the virtual 4C
# profiles.
#
# Output: a TSV file (total_nodups.window_vp<VP>_w<W>.tsv) with one
# row per sample, reporting the number of pairs in the viewpoint
# window. A TOTAL row summing all samples is appended. This TSV is
# used downstream by virtual4C_make_bw_tracks_and_pool.sh to compute
# the viewpoint-local normalization scale factor
# (1e6 / sum(pairs_in_window)).
#
# When called without viewpoint/window arguments (STANDARD mode),
# the script instead extracts total_nodups from pairtools stats
# reports and writes total_nodups.tsv. This is used for library-size
# normalization.
# ============================================================

DIR="${1:-}"
ENV_FILE="${2:-}"
VIEWPOINT="${3:-}"   # mm10 coord, es. 51677810
WINDOW="${4:-}"      # bp, es. 1500

IMAGE="${IMAGE:-quay.io/lucidif/microc:0.0.1}"

if [[ -z "$DIR" ]]; then
  echo "Usage: $0 <pairix_dir> [env_to_source.sh] [viewpoint_mm10] [window_bp]" >&2
  echo "If viewpoint or window are omitted: count standard total_nodups (all the reads)" >&2
  return 0 2>/dev/null || true
fi
if [[ ! -d "$DIR" ]]; then
  echo "ERROR: directory not found: $DIR" >&2
  return 0 2>/dev/null || true
fi

if [[ -n "$ENV_FILE" ]]; then
  if [[ -f "$ENV_FILE" ]]; then
    echo "Sourcing: $ENV_FILE"
    source "$ENV_FILE"
  else
    echo "WARN: env file not found, continuing without it: $ENV_FILE" >&2
  fi
fi

if [[ "${DOCKER_ARGS+set}" != "set" ]]; then
  DOCKER_ARGS=()
fi

# ----------------------------
# Modalità: window o standard?
# ----------------------------
USE_WINDOW=0
if [[ -n "$VIEWPOINT" && -n "$WINDOW" ]]; then
  USE_WINDOW=1
  echo "[INFO] WINDOW mode: viewpoint=${VIEWPOINT}, window=${WINDOW}"
  echo "[INFO] pair counts with at least on read (of the couple) in the window"
else
  echo "[INFO] STANDARD mode: count total_nodups from pairtools stats"
fi

shopt -s nullglob

# ----------------------------
# Funzione: conta pairs nella window da un singolo pairs.gz
# Replica la stessa logica di virtual4c_from_pairs_and_bam.sh passo 1
# ----------------------------
count_pairs_in_window() {
  local pairgz="$1"

  # Ricava CONTIG, TARGET_CHR, OFFSET dall'header del pairs file
  local CONTIG TARGET_CHR OFFSET VP_REL REL_START REL_END

  CONTIG=$(zcat -f "$pairgz" | awk '/^#chromsize:/ {print $2; exit}')
  if [[ -z "$CONTIG" ]]; then
    echo "[ERROR] No #chromsize found in $pairgz" >&2
    echo "NA"
    return
  fi

  TARGET_CHR=$(echo "$CONTIG" | sed -E 's/.*\|([^:]+):.*/\1/')
  OFFSET=$(echo "$CONTIG"     | sed -E 's/.*\|[^:]+:([0-9]+)-.*/\1/')

  if [[ -z "$TARGET_CHR" || -z "$OFFSET" ]]; then
    echo "[ERROR] Could not parse TARGET_CHR/OFFSET from CONTIG=$CONTIG" >&2
    echo "NA"
    return
  fi

  VP_REL=$(( VIEWPOINT - OFFSET ))
  REL_START=$(( VP_REL - WINDOW ))
  REL_END=$(( VP_REL + WINDOW + 1 ))
  if (( REL_START < 0 )); then REL_START=0; fi

  echo "[INFO]   CONTIG=${CONTIG} TARGET_CHR=${TARGET_CHR} OFFSET=${OFFSET}" >&2
  echo "[INFO]   rel window: ${REL_START}-${REL_END}" >&2

  # Conta coppie con almeno un'estremità nella window (stesso filtro del passo 1)
  local n
  n=$(sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -i \
    -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
    "$IMAGE" \
    pairtools select \
      "((chrom1==\"${CONTIG}\" and pos1>=${REL_START} and pos1<${REL_END}) or (chrom2==\"${CONTIG}\" and pos2>=${REL_START} and pos2<${REL_END}))" \
      "$pairgz" \
    | awk '!/^#/ {n++} END{print n+0}' 2>/dev/null || echo "NA")

  echo "$n"
}

# ----------------------------
# [1/2] pairtools stats (only standard mode)
# ----------------------------
found_any=0

if [[ "$USE_WINDOW" -eq 0 ]]; then
  echo "[1/2] Generating pairtools stats reports in: $DIR"
  for f in "$DIR"/*.Dd.pairs.gz; do
    found_any=1
    rep="${f%.pairs.gz}.pairs.report.txt"
    if [[ -s "$rep" ]]; then
      echo "  SKIP: $(basename "$rep")"
      continue
    fi
    echo "  RUN : $(basename "$f")"
    sudo docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" \
      -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
      "$IMAGE" \
      pairtools stats "$f" > "$rep" 2> "${rep}.stderr.txt" || true
  done
else
  echo "[1/2] Modalità window: pairtools stats skippato"
  for f in "$DIR"/*.Dd.pairs.gz; do found_any=1; break; done
fi

if [[ "$found_any" -eq 0 ]]; then
  echo "WARN: no *.Dd.pairs.gz found in: $DIR" >&2
  return 0 2>/dev/null || true
fi

# ----------------------------
# [2/2] Scrivi TSV
# ----------------------------
if [[ "$USE_WINDOW" -eq 1 ]]; then
  OUT_TSV="$DIR/total_nodups.window_vp${VIEWPOINT}_w${WINDOW}.tsv"
  COL_LABEL="pairs_in_window"
else
  OUT_TSV="$DIR/total_nodups.tsv"
  COL_LABEL="total_nodups"
fi

echo "[2/2] Writing: $OUT_TSV"
{
  echo -e "sample\t${COL_LABEL}"
  sum=0
  for f in "$DIR"/*.Dd.pairs.gz; do
    sample="$(basename "${f%.Dd.pairs.gz}")"
    echo "  Processing: ${sample}" >&2

    if [[ "$USE_WINDOW" -eq 1 ]]; then
      n=$(count_pairs_in_window "$f")
    else
      rep="${f%.pairs.gz}.pairs.report.txt"
      n="$(awk '$1=="total_nodups"{print $2}' "$rep" 2>/dev/null || true)"
    fi

    if [[ -z "${n:-}" || "$n" == "NA" ]]; then
      echo -e "${sample}\tNA"
      continue
    fi
    echo -e "${sample}\t${n}"
    sum=$((sum + n))
  done
  echo -e "TOTAL\t${sum}"
} > "$OUT_TSV"

echo "Done. Preview:"
tail "$OUT_TSV" || true