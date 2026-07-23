#!/usr/bin/env bash
#set -euo pipefail  # <-- NON usato (come richiesto)

# ============================================================
# virtual4C_post_microc_6to8.sh
#
# Steps 6-8:
#   6) pairs.gz -> .cool (BIN)
#   7b) .cool -> BEDPE (cooler dump --join) [unlifted]
#   8) lift/rebuild to mm10 + export .hic (lift_and_export_hic.sh)
#
# Inputs expected from nf-core/microc OUTDIR:
#   - ${NF_OUT}/pairix/*.pairs.gz
#   - ${NF_OUT}/chromsizes/genome.fa.sizes
#
# No pipefail, no exits: logs/warnings only.
# ============================================================

usage() {
  cat <<'EOF'
Usage:
  git/nf-core-microc/bin/virtual4C_rebuild_stdref.sh \
    --nf-out outs/virtual4C/nfout/chr10_Rfx6_WT_day4 \
    --bin 5000 \
    --outdir outs/virtual4C/test_virtual4C_chr10_Rfx6_WT_day4 \
    --docker-env git/Lara_MLL2/bin/docker_env.sh \
    --lift-script git/nf-core-microc/bin/lift_and_export_hic.sh \
    --chrom chr10 --offset 1234567 \
    --assembly mm10 \
    --threads 12

Alternative to --chrom/--offset:
  --region-bed outs/virtual4C/customRef/chr10_Rfx6.bed
  (chrom = col1, offset = col2 from first line)

Notes:
- Produces:
  <outdir>/custom_cool/*.cool
  <outdir>/bedpe_from_cool/*.bedpe (unlifted)
  <outdir>/_hic_mm10/<tag>/*.hic + lifted/rebuilt/tmp dirs
EOF
}

# ----------------------------
# Defaults
# ----------------------------
NF_OUT=""
BIN=5000
OUTDIR=""
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
LIFT_SCRIPT="git/nf-core-microc/bin/lift_and_export_hic.sh"
CHROM=""
OFFSET=""
REGION_BED=""
ASSEMBLY="mm10"
THREADS=12

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --nf-out)       NF_OUT="$2"; shift 2;;
    --bin)          BIN="$2"; shift 2;;
    --outdir)       OUTDIR="$2"; shift 2;;
    --docker-env)   DOCKER_ENV="$2"; shift 2;;
    --lift-script)  LIFT_SCRIPT="$2"; shift 2;;
    --chrom)        CHROM="$2"; shift 2;;
    --offset)       OFFSET="$2"; shift 2;;
    --region-bed)   REGION_BED="$2"; shift 2;;
    --assembly)     ASSEMBLY="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    -h|--help)      usage; exit 0;;
    *) echo "[WARN] Unknown arg: $1"; shift;;
  esac
done

echo "[INFO] NF_OUT     = ${NF_OUT}"
echo "[INFO] BIN        = ${BIN}"
echo "[INFO] OUTDIR     = ${OUTDIR}"
echo "[INFO] DOCKER_ENV = ${DOCKER_ENV}"
echo "[INFO] LIFT_SCRIPT= ${LIFT_SCRIPT}"
echo "[INFO] CHROM      = ${CHROM}"
echo "[INFO] OFFSET     = ${OFFSET}"
echo "[INFO] REGION_BED = ${REGION_BED}"
echo "[INFO] ASSEMBLY   = ${ASSEMBLY}"
echo "[INFO] THREADS    = ${THREADS}"

# ----------------------------
# Docker env (DOCKER_ARGS)
# ----------------------------
if [[ -f "${DOCKER_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${DOCKER_ENV}"
  echo "[INFO] Sourced docker env: ${DOCKER_ENV}"
else
  echo "[WARN] docker env not found: ${DOCKER_ENV}"
  echo "[WARN] Will try docker without DOCKER_ARGS."
  DOCKER_ARGS=()
fi

# ----------------------------
# Resolve chrom/offset
# ----------------------------
if [[ -n "${REGION_BED}" && -f "${REGION_BED}" ]]; then
  if [[ -z "${CHROM}" ]]; then
    CHROM="$(awk 'NR==1{print $1}' "${REGION_BED}")"
  fi
  if [[ -z "${OFFSET}" ]]; then
    OFFSET="$(awk 'NR==1{print $2}' "${REGION_BED}")"
  fi
  echo "[INFO] Parsed from REGION_BED -> CHROM=${CHROM} OFFSET=${OFFSET}"
fi

if [[ -z "${CHROM}" || -z "${OFFSET}" ]]; then
  echo "[WARN] Missing CHROM/OFFSET. Step 8 (lift/rebuild/export hic) will likely fail."
fi

# ----------------------------
# Check NF_OUT layout
# ----------------------------
PAIRIX_DIR="${NF_OUT}/pairix"
CHROMSIZES="${NF_OUT}/chromsizes/genome.fa.sizes"

if [[ ! -d "${PAIRIX_DIR}" ]]; then
  echo "[WARN] pairix directory not found: ${PAIRIX_DIR}"
fi
if [[ ! -f "${CHROMSIZES}" ]]; then
  echo "[WARN] chromsizes not found: ${CHROMSIZES}"
fi

# ----------------------------
# Output dirs
# ----------------------------
if [[ -z "${OUTDIR}" ]]; then
  echo "[WARN] Missing --outdir. Defaulting to: ${NF_OUT}/post_6to8"
  OUTDIR="${NF_OUT}/post_6to8"
fi

COOL_DIR="${OUTDIR}/custom_cool"
BEDPE_FROM_COOL_DIR="${OUTDIR}/bedpe_from_cool"

ONEFILE_DIR="${OUTDIR}/_onefile"
LIFTED_DIR="${OUTDIR}/_lifted_${ASSEMBLY}"
REBUILT_DIR="${OUTDIR}/_rebuilt_${ASSEMBLY}"
HIC_DIR="${OUTDIR}/_hic_${ASSEMBLY}"
TMP_DIR="${OUTDIR}/_tmp_rebuild"

mkdir -p "${OUTDIR}" "${COOL_DIR}" "${BEDPE_FROM_COOL_DIR}" \
         "${ONEFILE_DIR}" "${LIFTED_DIR}" "${REBUILT_DIR}" "${HIC_DIR}" "${TMP_DIR}"

# ============================================================
# 6) Build .cool from pairs.gz
# ============================================================
echo "[INFO] (6) Building ${BIN}bp .cool from pairs.gz"

shopt -s nullglob
pairs_list=( "${PAIRIX_DIR}"/*.pairs.gz )
shopt -u nullglob

if [[ "${#pairs_list[@]}" -eq 0 ]]; then
  echo "[WARN] No .pairs.gz found in: ${PAIRIX_DIR}"
else
  for P in "${pairs_list[@]}"; do
    sample="$(basename "${P}" .pairs.gz)"
    OUTCOOL="${COOL_DIR}/${BIN}bp_${sample}.cool"

    echo "[INFO] cooler cload pairs -> ${OUTCOOL}"
    sudo docker run --rm "${DOCKER_ARGS[@]}" \
      quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
      cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        "${CHROMSIZES}:${BIN}" "${P}" "${OUTCOOL}"

    if [[ -f "${OUTCOOL}" ]]; then
      echo "[INFO] Wrote cool: ${OUTCOOL}"
    else
      echo "[WARN] cool not created (check logs above): ${OUTCOOL}"
    fi
  done
fi

# ============================================================
# 7b) Export contacts from each UNLIFTED cool to BEDPE
# ============================================================
echo "[INFO] (7b) cooler dump --join from each unlifted cool -> BEDPE"

shopt -s nullglob
cool_list=( "${COOL_DIR}"/*.cool )
shopt -u nullglob

if [[ "${#cool_list[@]}" -eq 0 ]]; then
  echo "[WARN] No .cool found in: ${COOL_DIR}"
else
  for COOL in "${cool_list[@]}"; do
    base="$(basename "$COOL" .cool)"
    OUTBEDPE="${BEDPE_FROM_COOL_DIR}/${base}.${BIN}bp.cooldump.bedpe"
    ERRLOG="${OUTBEDPE}.err"

    echo "[INFO] cooler dump --join -> ${OUTBEDPE}"

    # non usiamo exit codes per fermare; logghiamo soltanto
    sudo docker run --rm "${DOCKER_ARGS[@]}" \
      quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
      cooler dump --join "$COOL" \
      > "$OUTBEDPE" 2> "$ERRLOG"

    nlines="$(wc -l < "$OUTBEDPE" 2>/dev/null | awk '{print $1}')"
    echo "[INFO] lines=${nlines} cool=$(basename "$COOL")"

    if [[ -s "$ERRLOG" ]]; then
      echo "[WARN] cooler dump stderr (last 20 lines):"
      tail -n 20 "$ERRLOG"
    fi
  done
fi

# ============================================================
# 8) Lift/rebuild/export hic for each cool
# ============================================================
echo "[INFO] (8) Lifting/rebuilding each cool to ${ASSEMBLY} and exporting .hic"

if [[ ! -f "${LIFT_SCRIPT}" ]]; then
  echo "[WARN] lift script not found: ${LIFT_SCRIPT}"
else
  chmod +x "${LIFT_SCRIPT}" 2>/dev/null || true
fi

if [[ -z "${CHROM}" || -z "${OFFSET}" ]]; then
  echo "[WARN] CHROM/OFFSET missing -> skipping lift/rebuild/export hic."
else
  if [[ "${#cool_list[@]}" -eq 0 ]]; then
    echo "[WARN] No cools available for step 8."
  else
    for COOL in "${cool_list[@]}"; do
      base="$(basename "${COOL}")"  # e.g. 5000bp_sample.cool
      tag="${base%.cool}"

      echo "[INFO] Processing lift/rebuild/hic: ${base}"

      # Reset ONEFILE_DIR for clean single-input run
      rm -f "${ONEFILE_DIR}"/*.cool 2>/dev/null || true
      cp -v "${COOL}" "${ONEFILE_DIR}/" 2>/dev/null || true

      # Run lift script
      "${LIFT_SCRIPT}" \
        --indir   "${ONEFILE_DIR}" \
        --lifted  "${LIFTED_DIR}/${tag}" \
        --rebuilt "${REBUILT_DIR}/${tag}" \
        --hicdir  "${HIC_DIR}/${tag}" \
        --tmpdir  "${TMP_DIR}/${tag}" \
        --chrom   "${CHROM}" \
        --offset  "${OFFSET}" \
        --assembly "${ASSEMBLY}" \
        --threads "${THREADS}"

      # Soft checks
      shopt -s nullglob
      hics=( "${HIC_DIR}/${tag}"/*.hic )
      shopt -u nullglob
      if [[ "${#hics[@]}" -gt 0 ]]; then
        echo "[INFO] HIC created: ${hics[0]}"
      else
        echo "[WARN] No .hic found under: ${HIC_DIR}/${tag}"
      fi
    done
  fi
fi

echo "[DONE] Steps 6-8 completed."
echo "[OUTPUT] cools:      ${COOL_DIR}"
echo "[OUTPUT] bedpe:      ${BEDPE_FROM_COOL_DIR}"
echo "[OUTPUT] hic root:   ${HIC_DIR}"
echo "[OUTPUT] lifted:     ${LIFTED_DIR}"
echo "[OUTPUT] rebuilt:    ${REBUILT_DIR}"
echo "[OUTPUT] tmp:        ${TMP_DIR}"
