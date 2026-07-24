#!/usr/bin/env bash
set -euo pipefail

SHEET="${1:?Usage: $0 GEO_rename_chipSEQ_bamFiles.tsv}"

OLDROOT="${OLDROOT:-/mnt/datawk1/analysis/Lara}"
NEWROOT="${NEWROOT:-/media/lucio/easystore/Lucio/Analysis/Lara}"

ok=0
missing=0

# header report
printf "status\tGEOprocessed\tbamPattern\tbamLocation\tchecked_path\n"

tail -n +2 "$SHEET" | while IFS=$'\t' read -r GEOprocessed bamPattern bamLocation; do
  [[ -z "${GEOprocessed:-}" || -z "${bamPattern:-}" || -z "${bamLocation:-}" ]] && continue

  # rimappa root se serve
  loc="$bamLocation"
  if [[ ! -d "$loc" ]]; then
    loc="${loc/#$OLDROOT/$NEWROOT}"
  fi

  checked="${loc}/${bamPattern}.mLb.mkD.sorted.bam"

  if [[ -f "$checked" ]]; then
    printf "OK\t%s\t%s\t%s\t%s\n" "$GEOprocessed" "$bamPattern" "$bamLocation" "$checked"
    ok=$((ok+1))
  else
    printf "MISSING\t%s\t%s\t%s\t%s\n" "$GEOprocessed" "$bamPattern" "$bamLocation" "$checked"
    missing=$((missing+1))
  fi
done

# Nota: i contatori ok/missing sopra non si aggiornano fuori dal while per via della subshell in pipe.
# Se ti serve il summary finale, dimmelo e te lo faccio in modo che lo stampi correttamente.
