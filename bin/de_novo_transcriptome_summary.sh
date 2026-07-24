#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Summarize outputs from the de novo transcriptome assembly workflow
# =============================================================================

OUTPUT_DIR="${1:-./outs/denovo_transcripts}"

if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "ERROR: Output directory not found: ${OUTPUT_DIR}" >&2
  exit 1
fi

echo "[INFO] Summary for: ${OUTPUT_DIR}"

echo "\n-- Directory tree --"
find "${OUTPUT_DIR}" -maxdepth 3 -type d | sort

echo "\n-- StringTie assembly files --"
find "${OUTPUT_DIR}" -type f \( -iname "*.gtf" -o -iname "*.gtf.gz" -o -iname "*stringtie*" \) | sort | head -n 50

echo "\n-- BigWig coverage files --"
find "${OUTPUT_DIR}" -type f \( -iname "*.bw" -o -iname "*.bigwig" \) | sort | head -n 50

echo "\n-- STAR/Salmon logs and run reports --"
find "${OUTPUT_DIR}" -type f \( -iname "*.log" -o -iname "*.html" -o -iname "*.txt" \) | sort | grep -E "(log|report|trace|timeline|params|README)" | head -n 100 || true

echo "\n-- Top-level stats files --"
for f in "${OUTPUT_DIR}/report.html" "${OUTPUT_DIR}/multiqc.html" "${OUTPUT_DIR}/stringtie/*.gtf"; do
  if [[ -e "${f}" ]]; then
    echo "  ${f}"
  fi
 done

if [[ -f "${OUTPUT_DIR}/DE_NOVO_README.txt" ]]; then
  echo "\n[INFO] Found DE_NOVO_README.txt"
  cat "${OUTPUT_DIR}/DE_NOVO_README.txt"
fi

echo "\n[DONE] De novo transcriptome summary complete"
