#!/usr/bin/env bash
#set -euo pipefail  # <-- come richiesto: NON lo usiamo

# ============================================================
# virtual4C_make_reference.sh
#
# Builds a locus-specific reference genome for virtual 4C analysis.
# The reference is a sub-sequence of the mm10 genome spanning the
# genomic region of interest, as defined by coordinates in a TSV
# annotation file. The resulting FASTA is indexed with bwa and
# samtools faidx for use as input to the Micro-C alignment pipeline.
#
# Steps:
#   1. Extract the locus coordinates from the TSV annotation file
#      and write a BED file (chrom, start, end, name).
#   2. Ensure that an mm10 .fai index is present; link it from a
#      source path if necessary.
#   3. Invoke bed_to_refs_dirs.sh to extract the locus sequence from
#      the mm10 FASTA and write it as genome.fa under
#      <custom-ref-root>/<ref-name>/.
#   4. Index the custom FASTA with bwa index and samtools faidx
#      (both run inside Docker).
# ============================================================

usage() {
  cat <<'EOF'
Usage:
  virtual4C_make_reference.sh \
    --ref-name chr10_Rfx6 \
    --tsv sheets/VirtualC_ref_coordinates.tsv \
    --mm10-fa /path/to/mm10.fa \
    --custom-ref-root in/customReferences \
    --outdir-main outs/virtual4C \
    --docker-env git/Lara_MLL2/bin/docker_env.sh \
    [--mm10-fai /path/to/mm10.fa.fai] \
    [--mm10-fai-src /path/to/mm10.fai] \
    [--tsv-name-col 11] [--tsv-chrom-col 4] [--tsv-start-col 9] [--tsv-end-col 10] \
    [--bed-subdir customRef]

Notes:
- Assumes TSV is TAB-separated.
- Default TSV columns (1-based): chrom=4, start=9, end=10, name=11.
- start is written as-is; if you need -1 adjustment, we can add an option later.
EOF
}

# ----------------------------
# Defaults
# ----------------------------
REF_NAME=""
TSV=""
MM10_FA=""
MM10_FAI=""
MM10_FAI_SRC=""
CUSTOM_REF_ROOT="in/customReferences"
OUTDIR_MAIN="outs/virtual4C"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
BED_SUBDIR="customRef"

TSV_NAME_COL=11
TSV_CHROM_COL=4
TSV_START_COL=9
TSV_END_COL=10

# ----------------------------
# Parse args (simple)
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref-name)        REF_NAME="$2"; shift 2;;
    --tsv)             TSV="$2"; shift 2;;
    --mm10-fa)         MM10_FA="$2"; shift 2;;
    --mm10-fai)        MM10_FAI="$2"; shift 2;;
    --mm10-fai-src)    MM10_FAI_SRC="$2"; shift 2;;
    --custom-ref-root) CUSTOM_REF_ROOT="$2"; shift 2;;
    --outdir-main)     OUTDIR_MAIN="$2"; shift 2;;
    --docker-env)      DOCKER_ENV="$2"; shift 2;;
    --bed-subdir)      BED_SUBDIR="$2"; shift 2;;

    --tsv-name-col)    TSV_NAME_COL="$2"; shift 2;;
    --tsv-chrom-col)   TSV_CHROM_COL="$2"; shift 2;;
    --tsv-start-col)   TSV_START_COL="$2"; shift 2;;
    --tsv-end-col)     TSV_END_COL="$2"; shift 2;;

    -h|--help) usage; exit 0;;
    *) echo "[WARN] Unknown arg: $1"; shift;;
  esac
done

# ----------------------------
# Minimal logging of inputs
# ----------------------------
echo "[INFO] REF_NAME=${REF_NAME}"
echo "[INFO] TSV=${TSV}"
echo "[INFO] MM10_FA=${MM10_FA}"
echo "[INFO] MM10_FAI=${MM10_FAI}"
echo "[INFO] MM10_FAI_SRC=${MM10_FAI_SRC}"
echo "[INFO] CUSTOM_REF_ROOT=${CUSTOM_REF_ROOT}"
echo "[INFO] OUTDIR_MAIN=${OUTDIR_MAIN}"
echo "[INFO] DOCKER_ENV=${DOCKER_ENV}"
echo "[INFO] TSV cols: chrom=${TSV_CHROM_COL} start=${TSV_START_COL} end=${TSV_END_COL} name=${TSV_NAME_COL}"

# ----------------------------
# Derived paths
# ----------------------------
mkdir -p "${OUTDIR_MAIN}" "${CUSTOM_REF_ROOT}"

OUTBED="${OUTDIR_MAIN}/${BED_SUBDIR}/${REF_NAME}.bed"
REFDIR="${CUSTOM_REF_ROOT}/${REF_NAME}"
FASTA_PATH="${REFDIR}/genome.fa"

mkdir -p "$(dirname "${OUTBED}")"

# ----------------------------
# Source docker env (DOCKER_ARGS)
# ----------------------------
if [[ -f "${DOCKER_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${DOCKER_ENV}"
  echo "[INFO] Sourced docker env: ${DOCKER_ENV}"
else
  echo "[WARN] docker env not found: ${DOCKER_ENV}"
  echo "[WARN] Will try running docker without DOCKER_ARGS."
  DOCKER_ARGS=()
fi

# ----------------------------
# 1) Extract coordinates from TSV -> BED
# ----------------------------
echo "[INFO] (1) Extracting ${REF_NAME} coordinates -> ${OUTBED}"

if [[ -z "${REF_NAME}" || -z "${TSV}" ]]; then
  echo "[WARN] Missing --ref-name or --tsv. Cannot extract BED."
else
  if [[ ! -f "${TSV}" ]]; then
    echo "[WARN] TSV not found: ${TSV}"
  else
    awk -F'\t' \
      -v NAME="${REF_NAME}" \
      -v NAME_COL="${TSV_NAME_COL}" \
      -v CHR_COL="${TSV_CHROM_COL}" \
      -v S_COL="${TSV_START_COL}" \
      -v E_COL="${TSV_END_COL}" '
    BEGIN { OFS="\t" }
    NR==1 { next }
    {
      # strip CRLF from the name field (common on Windows)
      gsub(/\r$/, "", $NAME_COL)
    }
    $NAME_COL == NAME {
      s = $S_COL
      if (s < 0) s = 0
      print $CHR_COL, s, $E_COL, $NAME_COL
    }
    ' "${TSV}" > "${OUTBED}"

    if [[ ! -s "${OUTBED}" ]]; then
      echo "[WARN] BED is empty: did not find REF_NAME=${REF_NAME} in TSV=${TSV} (name col ${TSV_NAME_COL})."
    else
      echo "[INFO] BED created:"
      cat "${OUTBED}"
    fi
  fi
fi

# ----------------------------
# 2) Ensure mm10 .fai exists
# ----------------------------
echo "[INFO] (2) Ensuring mm10 .fai exists"

# If MM10_FAI not provided, infer as mm10.fa.fai
if [[ -z "${MM10_FAI}" && -n "${MM10_FA}" ]]; then
  MM10_FAI="${MM10_FA}.fai"
  echo "[INFO] MM10_FAI inferred as: ${MM10_FAI}"
fi

if [[ -z "${MM10_FA}" ]]; then
  echo "[WARN] Missing --mm10-fa. Cannot check or build .fai."
else
  if [[ ! -f "${MM10_FA}" ]]; then
    echo "[WARN] mm10 FASTA not found: ${MM10_FA}"
  fi

  if [[ -n "${MM10_FAI}" && ! -f "${MM10_FAI}" ]]; then
    echo "[WARN] Missing mm10 index: ${MM10_FAI}"
    if [[ -n "${MM10_FAI_SRC}" && -f "${MM10_FAI_SRC}" ]]; then
      echo "[INFO] Linking mm10 .fai from source: ${MM10_FAI_SRC} -> ${MM10_FAI}"
      ln -sf "${MM10_FAI_SRC}" "${MM10_FAI}"
    else
      echo "[WARN] No MM10_FAI_SRC provided or file not found."
      echo "[WARN] You can create it with: samtools faidx ${MM10_FA}"
    fi
  else
    echo "[INFO] mm10 .fai OK: ${MM10_FAI}"
  fi
fi

# ----------------------------
# 3) Build custom FASTA + index directory
# ----------------------------
echo "[INFO] (3) Building custom reference dir under ${CUSTOM_REF_ROOT}"

if [[ ! -s "${OUTBED}" ]]; then
  echo "[WARN] BED missing/empty (${OUTBED}). Skipping bed_to_refs_dirs.sh."
else
  if [[ -z "${MM10_FA}" || ! -f "${MM10_FA}" ]]; then
    echo "[WARN] mm10 FASTA missing (${MM10_FA}). Cannot build custom reference."
  else
    if [[ ! -x "git/Lara_MLL2/bin/bed_to_refs_dirs.sh" ]]; then
      echo "[WARN] Helper not executable or not found: git/Lara_MLL2/bin/bed_to_refs_dirs.sh"
      echo "[WARN] Trying to run via bash anyway."
    fi

    bash git/Lara_MLL2/bin/bed_to_refs_dirs.sh \
      -b "${OUTBED}" \
      -f "${MM10_FA}" \
      -o "${CUSTOM_REF_ROOT}/" \
      --index

    if [[ -f "${FASTA_PATH}" ]]; then
      echo "[INFO] Custom FASTA created: ${FASTA_PATH}"
    else
      echo "[WARN] Expected FASTA not found: ${FASTA_PATH}"
    fi
  fi
fi

# ----------------------------
# 4) Index custom FASTA with bwa + samtools (docker)
# ----------------------------
echo "[INFO] (4) Indexing custom FASTA with bwa + samtools (docker)"

if [[ ! -f "${FASTA_PATH}" ]]; then
  echo "[WARN] FASTA not found: ${FASTA_PATH}. Skipping docker indexing."
else
  echo "[INFO] bwa index (docker) -> ${FASTA_PATH}"
  sudo docker run --rm "${DOCKER_ARGS[@]}" \
    quay.io/biocontainers/bwa:0.7.17--hed695b0_7 \
    bwa index "${FASTA_PATH}"

  echo "[INFO] samtools faidx (docker) -> ${FASTA_PATH}"
  sudo docker run --rm "${DOCKER_ARGS[@]}" \
    quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
    samtools faidx "${FASTA_PATH}"
fi

echo "[DONE] Reference generation finished."
echo "[OUTPUT] BED:   ${OUTBED}"
echo "[OUTPUT] REFDIR:${REFDIR}"
echo "[OUTPUT] FASTA: ${FASTA_PATH}"
