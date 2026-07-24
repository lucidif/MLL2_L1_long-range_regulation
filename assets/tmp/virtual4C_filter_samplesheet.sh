#!/usr/bin/env bash
#set -euo pipefail  # <-- richiesto: NON lo usiamo

# ============================================================
# virtual4C_filter_samplesheet.sh
#
# - Keep header (first line)
# - Filter lines (from line 2 onward) by one or more regex patterns
# - Case-insensitive by default (like grep -i)
# - Sanitize output:
#     1) CRLF -> LF
#     2) remove ANSI OSC escape sequences + control chars except TAB/LF
#
# No pipefail, no exit codes: only logs/warnings.
# ============================================================

usage() {
  cat <<'EOF'
Usage:
  virtual4C_filter_samplesheet.sh \
    --in sheets/NFSS_Lara_all_microC.csv \
    --out sheets/NFSS_Lara_WT_day4_microC.csv \
    --match 'WT' \
    --match 'day4|d4'

Options:
  --in FILE           Input samplesheet (first line is header)
  --out FILE          Output filtered samplesheet
  --match REGEX       Add a match condition (can be repeated). All must match (AND).
  --no-ignore-case    Make matches case-sensitive (default is ignore-case)
  --dry-run           Print how many lines would pass, and show first 5, but still writes output (lightweight)
EOF
}

IN=""
OUT=""
IGNORE_CASE=1
DRY_RUN=0
MATCHES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in) IN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --match) MATCHES+=("$2"); shift 2;;
    --no-ignore-case) IGNORE_CASE=0; shift;;
    --dry-run) DRY_RUN=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[WARN] Unknown arg: $1"; shift;;
  esac
done

echo "[INFO] Input : ${IN}"
echo "[INFO] Output: ${OUT}"
echo "[INFO] Matches (${#MATCHES[@]}):"
for m in "${MATCHES[@]}"; do
  echo "       - ${m}"
done
echo "[INFO] ignore-case: ${IGNORE_CASE}"
echo "[INFO] dry-run: ${DRY_RUN}"

if [[ -z "${IN}" || -z "${OUT}" ]]; then
  echo "[WARN] Missing --in or --out. Nothing to do."
  exit 0
fi

if [[ ! -f "${IN}" ]]; then
  echo "[WARN] Input file not found: ${IN}"
  # scriviamo comunque un output vuoto? no: non facciamo nulla.
  exit 0
fi

# Build grep pipeline
# We want: header + (tail -n +2 | grep ... | grep ... )
# If no --match is provided: copy all lines (header + body) + sanitize.

mkdir -p "$(dirname "${OUT}")"

tmp_out="${OUT}.tmp.$$"
: > "${tmp_out}"

# 1) copy header
head -n 1 "${IN}" >> "${tmp_out}"

# 2) filter body
if [[ "${#MATCHES[@]}" -eq 0 ]]; then
  echo "[WARN] No --match provided: copying all rows (no filtering)."
  tail -n +2 "${IN}" >> "${tmp_out}"
else
  # Start from body
  if [[ "${IGNORE_CASE}" -eq 1 ]]; then
    body_cmd=(tail -n +2 "${IN}")
    # apply AND matches
    # shellcheck disable=SC2034
    for m in "${MATCHES[@]}"; do
      body_cmd+=( "|" grep -i -E "${m}" )
    done
  else
    body_cmd=(tail -n +2 "${IN}")
    for m in "${MATCHES[@]}"; do
      body_cmd+=( "|" grep -E "${m}" )
    done
  fi

  # Execute pipeline safely without eval gymnastics:
  # We’ll stream through sequential greps.
  stream="$(tail -n +2 "${IN}")"
  for m in "${MATCHES[@]}"; do
    if [[ "${IGNORE_CASE}" -eq 1 ]]; then
      stream="$(printf "%s\n" "${stream}" | grep -i -E "${m}")"
    else
      stream="$(printf "%s\n" "${stream}" | grep -E "${m}")"
    fi
  done
  printf "%s\n" "${stream}" >> "${tmp_out}"
fi

# Sanitize (same logic as your snippet)
# 1) CRLF -> LF
sed -i 's/\r$//' "${tmp_out}"

# 2) remove ESC OSC + control chars except TAB+LF
LC_ALL=C perl -pe 's/\x1b\][^\x07]*(\x07|\x1b\\)//g; s/[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]//g' -i "${tmp_out}"

# Move into place
mv -f "${tmp_out}" "${OUT}"

# Report
nlines="$(wc -l < "${OUT}" 2>/dev/null | awk '{print $1}')"
echo "[INFO] Wrote: ${OUT} (lines=${nlines})"

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "[INFO] Preview (first 5 lines):"
  head -n 5 "${OUT}"
fi

# If only header remains, warn (no exits)
if [[ "${nlines}" -le 1 ]]; then
  echo "[WARN] Filtered sheet has no samples (only header)."
fi
