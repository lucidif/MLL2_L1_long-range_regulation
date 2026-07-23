#!/usr/bin/env bash
set -eu

DIR="${1:-}"
ENV_FILE="${2:-}"
IMAGE="${IMAGE:-quay.io/lucidif/microc:0.0.1}"

if [[ -z "$DIR" ]]; then
  echo "Usage: $0 <pairix_dir> [env_to_source.sh]" >&2
  return 0 2>/dev/null || true
fi

if [[ ! -d "$DIR" ]]; then
  echo "ERROR: directory not found: $DIR" >&2
  return 0 2>/dev/null || true
fi

# Optional source of environment (e.g., defines DOCKER_ARGS)
if [[ -n "$ENV_FILE" ]]; then
  if [[ -f "$ENV_FILE" ]]; then
    echo "Sourcing: $ENV_FILE"
    # shellcheck disable=SC1090
    source "$ENV_FILE"
  else
    echo "WARN: env file not found, continuing without it: $ENV_FILE" >&2
  fi
fi

# Default DOCKER_ARGS to empty if not defined by sourced env
if [[ "${DOCKER_ARGS+set}" != "set" ]]; then
  DOCKER_ARGS=()
fi

shopt -s nullglob

echo "[1/2] Generating pairtools stats reports in: $DIR"

found_any=0
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

if [[ "$found_any" -eq 0 ]]; then
  echo "WARN: no *.Dd.pairs.gz found in: $DIR" >&2
  return 0 2>/dev/null || true
fi

OUT_TSV="$DIR/total_nodups.tsv"
echo "[2/2] Writing: $OUT_TSV"

{
  echo -e "sample\ttotal_nodups"
  sum=0

  for f in "$DIR"/*.Dd.pairs.gz; do
    rep="${f%.pairs.gz}.pairs.report.txt"
    sample="$(basename "${f%.pairs.gz}")"
    n="$(awk '$1=="total_nodups"{print $2}' "$rep" 2>/dev/null || true)"

    if [[ -z "${n:-}" ]]; then
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
