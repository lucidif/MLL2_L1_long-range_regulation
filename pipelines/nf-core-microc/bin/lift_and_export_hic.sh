#!/usr/bin/env bash
# SAFE MODE: never exit the shell / terminal.
# - No set -euo pipefail
# - No exit anywhere
# - die() only marks failure and continues (or skips later steps)

FAIL=0

###############################################################################
# lift_and_export_hic_safe.sh
###############################################################################

# ----------------------------- defaults --------------------------------------
INDIR=""
LIFTED_DIR=""
REBUILT_DIR=""
HIC_DIR=""
TMP_ROOT=""

CHROM=""
OFFSET=""
ASSEMBLY="mm10"
THREADS=12

DOCKER_MOUNTS=(
  -v "$(pwd)":"$(pwd)"
  -v /media/lucio/easystore:/media/lucio/easystore
  -v /mnt/datawk1:/mnt/datawk1
)

COOLER_IMAGE="quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0"
HICTK_IMAGE="ghcr.io/paulsengroup/hictk:sha-75a4ba1"

# ----------------------------- helpers ---------------------------------------
die() {
  echo "[ERROR] $*" >&2
  FAIL=1
  return 0
}

info() { echo "[INFO] $*" >&2; }
warn() { echo "[WARN] $*" >&2; }

run() {
  # run CMD... ; never exits; records FAIL if command fails
  echo "+ $*" >&2
  "$@"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    warn "Command failed (rc=$rc): $*"
    FAIL=1
  fi
  return 0
}

docker_run() {
  # docker_run IMAGE CMD...
  local img="$1"; shift
  run sudo docker run --rm \
    "${DOCKER_MOUNTS[@]}" \
    -w "$(pwd)" -u "$(id -u)":"$(id -g)" \
    "$img" \
    "$@"
  return 0
}

usage() {
  sed -n '1,140p' "$0" | sed 's/^# \{0,1\}//'
  return 0
}

# ----------------------------- parse args ------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --indir)     INDIR="$2"; shift 2;;
    --lifted)    LIFTED_DIR="$2"; shift 2;;
    --rebuilt)   REBUILT_DIR="$2"; shift 2;;
    --hicdir)    HIC_DIR="$2"; shift 2;;
    --tmpdir)    TMP_ROOT="$2"; shift 2;;
    --chrom)     CHROM="$2"; shift 2;;
    --offset)    OFFSET="$2"; shift 2;;
    --assembly)  ASSEMBLY="$2"; shift 2;;
    --threads)   THREADS="$2"; shift 2;;
    -h|--help)   usage; return 0;;
    *) die "Unknown argument: $1 (use --help)"; shift;;
  esac
done

# ----------------------------- validate args ---------------------------------
[[ -n "$INDIR" ]]       || die "--indir is required"
[[ -n "$LIFTED_DIR" ]]  || die "--lifted is required"
[[ -n "$REBUILT_DIR" ]] || die "--rebuilt is required"
[[ -n "$HIC_DIR" ]]     || die "--hicdir is required"
[[ -n "$TMP_ROOT" ]]    || die "--tmpdir is required"
[[ -n "$CHROM" ]]       || die "--chrom is required"
[[ -n "${OFFSET:-}" ]]  || die "--offset is required"
[[ "$OFFSET" =~ ^[0-9]+$ ]] || die "--offset must be an integer"

[[ -d "$INDIR" ]] || die "INDIR not found: $INDIR"

# Create output dirs if possible
run mkdir -p "$LIFTED_DIR" "$REBUILT_DIR" "$HIC_DIR" "$TMP_ROOT"

info "Input cools:      $INDIR"
info "Lifted cools:     $LIFTED_DIR"
info "Rebuilt cools:    $REBUILT_DIR"
info "Output hics:      $HIC_DIR"
info "Tmp rebuild dir:  $TMP_ROOT"
info "Target chrom:     $CHROM"
info "Offset:           $OFFSET"
info "Assembly label:   $ASSEMBLY"
info "Threads:          $THREADS"
echo

# ----------------------------- main ------------------------------------------

FAIL=0

# STEP 1: Lift
if [[ $FAIL -eq 0 ]]; then
  info "[STEP 1/4] Lifting cool files..."
  found_any=0
  for IN in "$INDIR"/*.cool; do
    if [[ ! -e "$IN" ]]; then
      continue
    fi
    found_any=1
    BN="$(basename "$IN" .cool)"
    OUT="$LIFTED_DIR/${BN}.lifted.${ASSEMBLY}.cool"
    info "[LIFT] $IN -> $OUT"
    run ./git/Lara_MLL2/bin/lift_cool.py \
      -i "$IN" \
      -o "$OUT" \
      --chrom "$CHROM" \
      --offset "$OFFSET" \
      --assembly "$ASSEMBLY"
  done
  [[ $found_any -eq 1 ]] || die "No .cool files found in $INDIR"
  echo
else
  warn "[SKIP] STEP 1 skipped due to previous errors."
fi

# STEP 2: QC
if [[ $FAIL -eq 0 ]]; then
  info "[STEP 2/4] Quick QC on lifted cool files (python/cooler)..."
  run python3 - <<PY
import glob, cooler, os, sys
indir="${LIFTED_DIR}"
paths=sorted(glob.glob(indir+"/*.cool"))
if not paths:
    print("[QC] No lifted cool files found in", indir, file=sys.stderr)
    sys.exit(1)
for path in paths:
    c=cooler.Cooler(path)
    bins=c.bins()[:]
    print(os.path.basename(path),
          "chroms", c.chromnames,
          "minStart", int(bins["start"].min()),
          "maxEnd", int(bins["end"].max()),
          "nnz", c.info.get("nnz"),
          "sum", c.info.get("sum"))
PY
  echo
else
  warn "[SKIP] STEP 2 skipped due to previous errors."
fi

# STEP 3: Rebuild
# if [[ $FAIL -eq 0 ]]; then
#   info "[STEP 3/4] Rebuilding lifted cool files..."
#   found_any=0
#   for IN in "$LIFTED_DIR"/*.cool; do
#     [[ -e "$IN" ]] || continue
#     found_any=1

#     BN="$(basename "$IN" .cool)"
#     REB="$REBUILT_DIR/${BN}.rebuilt.cool"
#     TMP="$TMP_ROOT/${BN}"
#     run mkdir -p "$TMP"

#     info "[REBUILD] $IN -> $REB"

#     # validate (do not enable set -e ever)
#     docker_run "$HICTK_IMAGE" validate "$IN" >/dev/null 2>&1
#     ok=$?

#     if [[ $ok -eq 0 ]]; then
#       info "[REBUILD] hictk validate PASSED; copying"
#       run cp -f "$IN" "$REB"
#     else
#       info "[REBUILD] validate FAILED -> dump+load"

#       # chrom.sizes
#       run bash -lc "python3 - <<'PY' > '$TMP/chrom.sizes'
# import cooler
# c=cooler.Cooler('$IN')
# bins=c.bins()[:]
# chrom=bins['chrom'].iloc[0]
# chromsize=int(bins['end'].max())
# print(f\"{chrom}\\t{chromsize}\")
# PY"

#       # bin size
#       BINSIZE="$(python3 - <<PY
# import cooler
# c=cooler.Cooler("${IN}")
# print(int(c.info.get("bin-size")))
# PY
# )"
#       if [[ -z "${BINSIZE}" ]]; then
#         die "Could not infer bin-size from $IN"
#         continue
#       fi

#       # dump pixels (IMPORTANT: unbalanced; do not pass '--balanced false')
#       docker_run "$COOLER_IMAGE" bash -lc \
#         "cooler dump -t pixels --join '${IN}' > '${TMP}/pixels.tsv'"

#       # reload
#       docker_run "$COOLER_IMAGE" bash -lc \
#         "cooler load -f bg2 '${TMP}/chrom.sizes:${BINSIZE}' '${TMP}/pixels.tsv' '${REB}'"

#       # validate rebuilt
#       docker_run "$HICTK_IMAGE" validate "$REB"
#     fi
#     echo
#   done
#   [[ $found_any -eq 1 ]] || die "No lifted cool files found in $LIFTED_DIR"
# else
#   warn "[SKIP] STEP 3 skipped due to previous errors."
# fi

# 3) Rebuild each lifted .cool to ensure strict/valid Cooler indexes for hictk
echo "[STEP 3/4] Rebuilding lifted cool files into strict/valid .cool files..."

for IN in "$LIFTED_DIR"/*.cool; do
  [[ -f "$IN" ]] || continue

  BN=$(basename "$IN" .cool)
  REB="$REBUILT_DIR/${BN}.rebuilt.cool"

  TMP="$TMP_ROOT/${BN}"
  mkdir -p "$TMP"

  echo "[REBUILD] $IN -> $REB"

  # bin size from the cool itself
  BINSIZE=$(python3 - <<PY
import cooler
c=cooler.Cooler("${IN}")
print(int(c.info["bin-size"]))
PY
)

  # chrom.sizes (single chrom)
  python3 - <<PY > "$TMP/chrom.sizes"
import cooler
c=cooler.Cooler("${IN}")
bins=c.bins()[:]
chrom=bins["chrom"].iloc[0]
chromsize=int(bins["end"].max())
print(f"{chrom}\t{chromsize}")
PY

  # Dump pixels as bg2-like (7 columns)
  docker_run "$COOLER_IMAGE" bash -lc \
    "cooler dump --join '${IN}' > '${TMP}/pixels.bg2'"

  # Rebuild (regenerates ALL indexes, including /indexes/bin1_offset)
  docker_run "$COOLER_IMAGE" bash -lc \
    "rm -f '${REB}' && cooler load -f bg2 '${TMP}/chrom.sizes:${BINSIZE}' '${TMP}/pixels.bg2' '${REB}'"
done
echo

# STEP 4: Convert
# STEP 4: Convert rebuilt cool files to hic (ALWAYS try)
echo "[STEP 4/4] Converting rebuilt cool files to .hic (hictk convert)..."

shopt -s nullglob
rebuilt_list=( "$REBUILT_DIR"/*.cool )
shopt -u nullglob

if [[ ${#rebuilt_list[@]} -eq 0 ]]; then
  echo "[WARN] No rebuilt .cool files found in $REBUILT_DIR; nothing to convert."
else
  for IN in "${rebuilt_list[@]}"; do
    BN="$(basename "$IN" .cool)"
    OUT="$HIC_DIR/${BN}.hic"

    echo "[CONVERT] $IN -> $OUT"
    docker_run "$HICTK_IMAGE" convert -t "$THREADS" "$IN" "$OUT"
    docker_run "$HICTK_IMAGE" validate "$OUT"
  done
fi


echo
if [[ $FAIL -eq 0 ]]; then
  info "[DONE] Completed without fatal errors."
else
  warn "[DONE] Completed with errors (FAIL=1). Check warnings above."
fi
info "Lifted cools:  $LIFTED_DIR"
info "Rebuilt cools: $REBUILT_DIR"
info "HIC files:     $HIC_DIR"

# Always return success to avoid closing terminal sessions that auto-close on nonzero exit
return 0 2>/dev/null || exit 0
