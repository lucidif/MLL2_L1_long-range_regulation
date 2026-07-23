#!/usr/bin/env bash
shopt -s nullglob

# ============================================================
# virtual4C_make_bw_tracks_and_pool.sh
#
# A) Run virtual4c_from_pairs_and_bam.sh for each sample (pairs.gz + bam)
# B) (Optional) Pool lifted.mm10 bams for a SAM prefix -> merged bam (+ index)
# C) BigWig from pooled bam (bamCoverage in docker)
#    - If --pairs-counts-tsv provided: uses --scaleFactor = 1e6 / SUM(total_nodups)
#    - Else: fallback to bamCoverage --normalizeUsing (default BPM)
#
# No pipefail. Minimal exits only for --help.
# ============================================================
# ============================================================
# 
# Pools per-sample virtual 4C BAM files for a given condition
# (genotype x differentiation day) and generates a normalized
# BigWig coverage track.
#
# Steps:
#   A. For each sample, run virtual4c_from_pairs_and_bam.sh to
#      extract viewpoint-filtered, mm10-lifted BAM files and
#      per-sample BigWig tracks.
#   B. Merge all per-sample lifted BAM files (*.viewpoint.filtered.
#      lifted.mm10.bam) for the condition into a single pooled BAM
#      using samtools merge; index with samtools index.
#   C. Generate a BigWig from the pooled BAM using bamCoverage
#      (deepTools v3.5.5, --binSize 500).
#      Normalization is controlled by --pairs-counts-tsv:
#        - If a TSV of total_nodups is provided (library-size mode),
#          scaleFactor = 1e6 / sum(total_nodups), normalizing to
#          reads per million valid pairs (RPM).
#        - If a TSV of viewpoint-window pair counts is provided
#          (viewpoint-local mode), scaleFactor = 1e6 / sum(pairs_in_
#          window), normalizing to reads per million viewpoint-
#          captured pairs. Output files are tagged "RPMgw".
#        - If no TSV is provided, bamCoverage --normalizeUsing <NORM>
#          is used as fallback.
# ============================================================


usage() {
  cat <<'EOF'
Usage:
  virtual4C_make_bw_tracks_and_pool.sh \
    --ref-name chr10_Rfx6 \
    --viewpoint 51677810 \
    --window 1500 \
    --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
    --nf-out outs/virtual4C/nfout/chr10_Rfx6_WT_day4 \
    --outroot outs/virtual4C \
    --sam-prefix WT_day4 \
    --docker-env git/Lara_MLL2/bin/docker_env.sh \
    --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
    --pairs-counts-tsv outs/virtual4C/nfout/chr10_Rfx6_WT_day4/pairix/total_nodups.tsv \
    [--no-pool]

Options:
  --no-pool                 Skip pooling + bamCoverage (only per-sample v4c step A)
  --do-pool                 Force pooling + bamCoverage (default)

  --pair-glob '*.Dd.pairs.gz'        (default)
  --bam-suffix '.Lb.bam'             (default)

  --pairs-counts-tsv <file>          TSV with header: sample<TAB>total_nodups and optional TOTAL row.
                                    If provided, BigWig uses:
                                      scaleFactor = 1e6 / SUM(total_nodups excluding NA and TOTAL)

bamCoverage options:
  --binSize 500                      (default)
  --smoothLength 5000                (default)
  --norm BPM                         (default; used only if no --pairs-counts-tsv)
  --threads 8                        (default; used for samtools + bamCoverage)
EOF
}

# ----------------------------
# Defaults
# ----------------------------
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
V4C_SCRIPT="git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh"

REF_NAME=""
VIEWPOINT=""
WINDOW=""
MM10_SIZES=""
NF_OUT=""
OUTROOT="outs/virtual4C"
SAM_PREFIX=""

PAIR_GLOB="*.Dd.pairs.gz"
BAM_SUFFIX=".Lb.bam"

PAIRS_COUNTS_TSV=""   # new

BINSIZE=500
SMOOTH=5000
NORM="BPM"
THREADS=8

DO_POOL=1   # 1=pool merge + bigwig, 0=skip

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --docker-env)  DOCKER_ENV="$2"; shift 2;;
    --v4c-script)  V4C_SCRIPT="$2"; shift 2;;

    --ref-name)    REF_NAME="$2"; shift 2;;
    --viewpoint)   VIEWPOINT="$2"; shift 2;;
    --window)      WINDOW="$2"; shift 2;;
    --mm10-sizes)  MM10_SIZES="$2"; shift 2;;
    --nf-out)      NF_OUT="$2"; shift 2;;
    --outroot)     OUTROOT="$2"; shift 2;;
    --sam-prefix)  SAM_PREFIX="$2"; shift 2;;

    --pair-glob)   PAIR_GLOB="$2"; shift 2;;
    --bam-suffix)  BAM_SUFFIX="$2"; shift 2;;

    --pairs-counts-tsv) PAIRS_COUNTS_TSV="$2"; shift 2;;

    --binSize)       BINSIZE="$2"; shift 2;;
    --smoothLength)  SMOOTH="$2"; shift 2;;
    --norm)          NORM="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;

    --do-pool)     DO_POOL=1; shift;;
    --no-pool)     DO_POOL=0; shift;;

    -h|--help) usage; exit 0;;
    *) echo "[WARN] Unknown arg: $1"; shift;;
  esac
done

echo "[INFO] REF_NAME          = ${REF_NAME}"
echo "[INFO] VIEWPOINT         = ${VIEWPOINT}"
echo "[INFO] WINDOW            = ${WINDOW}"
echo "[INFO] MM10_SIZES         = ${MM10_SIZES}"
echo "[INFO] NF_OUT             = ${NF_OUT}"
echo "[INFO] OUTROOT            = ${OUTROOT}"
echo "[INFO] SAM_PREFIX         = ${SAM_PREFIX}"
echo "[INFO] PAIR_GLOB          = ${PAIR_GLOB}"
echo "[INFO] BAM_SUFFIX         = ${BAM_SUFFIX}"
echo "[INFO] PAIRS_COUNTS_TSV   = ${PAIRS_COUNTS_TSV}"
echo "[INFO] THREADS            = ${THREADS}"
echo "[INFO] DO_POOL            = ${DO_POOL}"
echo "[INFO] bamCoverage NORM   = ${NORM} (fallback only)"
echo "[INFO] bamCoverage bin/sm = ${BINSIZE}/${SMOOTH}"

# ----------------------------
# Docker env (DOCKER_ARGS)
# ----------------------------
if [[ -f "${DOCKER_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${DOCKER_ENV}"
  echo "[INFO] Sourced docker env: ${DOCKER_ENV}"
else
  echo "[WARN] docker env not found: ${DOCKER_ENV} -> docker runs without DOCKER_ARGS"
  DOCKER_ARGS=()
fi

if [[ "${DOCKER_ARGS+set}" != "set" ]]; then
  DOCKER_ARGS=()
fi

# ----------------------------
# Paths
# ----------------------------
PAIRDIR="${NF_OUT}/pairix"
BAMDIR="${NF_OUT}/bwamem"

OUTDIR="${OUTROOT}/bw_tracks_${REF_NAME}/${WINDOW}"
mkdir -p "${OUTDIR}"

# ----------------------------
# A) Run per-sample virtual4C extraction (pairs + bam)
# ----------------------------
echo "[INFO] A) virtual4c_from_pairs_and_bam per sample"

if [[ ! -e "${V4C_SCRIPT}" ]]; then
  echo "[WARN] V4C script missing: ${V4C_SCRIPT}"
  echo "[WARN] Will try to run anyway (may fail)."
else
  chmod +x "${V4C_SCRIPT}" 2>/dev/null || true
fi

if [[ ! -d "${PAIRDIR}" ]]; then
  echo "[WARN] PAIRDIR not found: ${PAIRDIR}"
fi
if [[ ! -d "${BAMDIR}" ]]; then
  echo "[WARN] BAMDIR not found: ${BAMDIR}"
fi
if [[ ! -f "${MM10_SIZES}" ]]; then
  echo "[WARN] mm10 sizes not found: ${MM10_SIZES}"
fi

pairs=( "${PAIRDIR}"/${PAIR_GLOB} )
if (( ${#pairs[@]} == 0 )); then
  echo "[WARN] No pairs found with: ${PAIRDIR}/${PAIR_GLOB}"
else
  for PAIRGZ in "${pairs[@]}"; do
    base="$(basename "${PAIRGZ}")"  # es WT_day4_A_LP1.Dd.pairs.gz

    # ricava sample (adatto al tuo suffisso .Dd.pairs.gz)
    sample="${base%.Dd.pairs.gz}"
    BAM_IN="${BAMDIR}/${sample}${BAM_SUFFIX}"

    if [[ ! -s "${BAM_IN}" ]]; then
      echo "[WARN] BAM missing for ${sample}: ${BAM_IN} (skip)"
      continue
    fi

    PREFIX="${sample}"
    mkdir -p "${OUTDIR}/${PREFIX}"

    echo "[INFO] ==> ${PREFIX}"
    "${V4C_SCRIPT}" \
      --pairs "${PAIRGZ}" \
      --bam "${BAM_IN}" \
      --viewpoint "${VIEWPOINT}" \
      --window "${WINDOW}" \
      --sizes "${MM10_SIZES}" \
      --outdir "${OUTDIR}/${PREFIX}" \
      --prefix "${PREFIX}"
  done
fi

# ----------------------------
# B) Pool lifted.mm10 bams for SAM_PREFIX
# ----------------------------
POOL_DIR="${OUTDIR}"
LIST_TXT="${POOL_DIR}/${SAM_PREFIX}.pooled.bams.txt"
MERGED_DIR="${POOL_DIR}/${SAM_PREFIX}"
MERGED_BAM="${MERGED_DIR}/${SAM_PREFIX}.pooled.viewpoint.filtered.lifted.mm10.bam"

if [[ "${DO_POOL}" -eq 1 ]]; then
  echo "[INFO] B) Pooling lifted.mm10 bams for SAM_PREFIX=${SAM_PREFIX}"

  mkdir -p "${MERGED_DIR}"

  # pattern: <POOL_DIR>/<SAM>_*/<SAM>_*.viewpoint.filtered.lifted.mm10.bam
  BAMS=( "${POOL_DIR}"/${SAM_PREFIX}_*/${SAM_PREFIX}_*.viewpoint.filtered.lifted.mm10.bam )

  if (( ${#BAMS[@]} == 0 )); then
    echo "[WARN] No lifted.mm10 bams found with pattern:"
    echo "       ${POOL_DIR}/${SAM_PREFIX}_*/${SAM_PREFIX}_*.viewpoint.filtered.lifted.mm10.bam"
  else
    printf "%s\n" "${BAMS[@]}" > "${LIST_TXT}"
    n_bams="$(wc -l < "${LIST_TXT}" 2>/dev/null | awk '{print $1}')"
    echo "[INFO] Wrote list: ${LIST_TXT} (${n_bams} bams)"

    echo "[INFO] samtools merge -> ${MERGED_BAM}"
    samtools merge -@ "${THREADS}" -f "${MERGED_BAM}" -b "${LIST_TXT}"
    samtools index "${MERGED_BAM}"
    echo "[INFO] OK: ${MERGED_BAM} (+ .bai)"
  fi
else
  echo "[INFO] B) Pool disabled (--no-pool). Skipping merge/index."
fi

# ----------------------------
# Compute scaleFactor from PAIRS_COUNTS_TSV (optional)
# ----------------------------
SCALE_FACTOR=""
N_TOTAL=0

if [[ -n "${PAIRS_COUNTS_TSV}" && -s "${PAIRS_COUNTS_TSV}" ]]; then
  # sum numeric $2 excluding header and TOTAL row and NA
  N_TOTAL=$(awk 'NR>1 && $1!="TOTAL" && $2!="NA"{s+=$2} END{print s+0}' "${PAIRS_COUNTS_TSV}" 2>/dev/null || echo 0)

  if [[ "${N_TOTAL}" -gt 0 ]]; then
    if [[ "${N_TOTAL}" -gt 0 ]]; then
      SCALE_FACTOR=$(LC_ALL=C awk -v n="${N_TOTAL}" 'BEGIN{printf "%.12f", 1000000/n}')
      SCALE_FACTOR="${SCALE_FACTOR/,/.}"   # safety: comma -> dot
      echo "[INFO] ScaleFactor from TSV: N_TOTAL=${N_TOTAL} -> scaleFactor=${SCALE_FACTOR}"
    fi
    
    echo "[INFO] ScaleFactor from TSV: N_TOTAL=${N_TOTAL} -> scaleFactor=${SCALE_FACTOR}"
  else
    echo "[WARN] TSV present but N_TOTAL computed as 0: ${PAIRS_COUNTS_TSV}"
  fi
else
  if [[ -n "${PAIRS_COUNTS_TSV}" ]]; then
    echo "[ERROR] --pairs-counts-tsv specified but file not found or empty: ${PAIRS_COUNTS_TSV}"
    exit 1
  fi
fi

# ----------------------------
# C) BigWig from pooled bam (bamCoverage)
# ----------------------------
#BW_OUT="${MERGED_DIR}/${SAM_PREFIX}.pooled.${NORM}.bwbin${BINSIZE}.smooth${SMOOTH}.bw"
if [[ -n "${SCALE_FACTOR}" ]]; then
  NORM_TAG="RPMgw"   # o "scaleFactor" / "RPM_pairsSum" come preferisci
else
  NORM_TAG="${NORM}"
fi

BW_OUT="${MERGED_DIR}/${SAM_PREFIX}.pooled.${NORM_TAG}.bwbin${BINSIZE}.smooth${SMOOTH}.bw"


if [[ "${DO_POOL}" -eq 1 ]]; then
  echo "[INFO] C) bamCoverage -> BigWig"
  if [[ ! -s "${MERGED_BAM}" ]]; then
    echo "[WARN] Merged BAM missing: ${MERGED_BAM} -> skipping bamCoverage"
  else
    if [[ -n "${SCALE_FACTOR}" ]]; then
      echo "[INFO] bamCoverage mode: --scaleFactor ${SCALE_FACTOR} (RPM from pairs sum)"
    else
      echo "[INFO] bamCoverage mode: --normalizeUsing ${NORM} (fallback)"
    fi

    sudo -n docker run --rm "${DOCKER_ARGS[@]}" -w "$PWD" -it \
      -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mpl -e XDG_CACHE_HOME=/tmp/.cache \
      quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0 \
      bamCoverage \
        -b "${MERGED_BAM}" \
        -o "${BW_OUT}" \
        --outFileFormat bigwig \
        --binSize "${BINSIZE}" \
        --smoothLength "${SMOOTH}" \
        --numberOfProcessors "${THREADS}" \
        $( [[ -n "${SCALE_FACTOR}" ]] && echo "--scaleFactor ${SCALE_FACTOR}" ) \
        $( [[ -z "${SCALE_FACTOR}" ]] && echo "--normalizeUsing ${NORM}" )

    ls -lh "${BW_OUT}" 2>/dev/null || true  
  fi
else
  echo "[INFO] C) Pool disabled (--no-pool). Skipping bamCoverage."
fi

echo "[DONE] BigWig pipeline finished."
echo "[OUTPUT] OUTDIR: ${OUTDIR}"
echo "[OUTPUT] (if pooled) POOLED_BAM: ${MERGED_BAM}"
echo "[OUTPUT] (if pooled) BIGWIG:     ${BW_OUT}"
