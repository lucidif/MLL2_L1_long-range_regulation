#!/bin/bash

INPUT="sheets/Dsplit_all_samples.tsv"  # <-- aggiusta il path
OUTDIR="outs/quantile_normalization_analysis/tmp"  # <-- aggiusta se serve

# Normalizza i line endings in-memory
tail -n +2 "$INPUT" | tr -d '\r' | while IFS=$'\t' read -r id antibody group; do
    antibody=$(echo "$antibody" | tr -d '"' | xargs)
    group=$(echo "$group" | tr -d '"' | xargs)

    [[ "$antibody" == "INPUT" ]] && continue
    [[ -z "$antibody" || -z "$group" ]] && continue

    outfile="${OUTDIR}/${group}_${antibody}.tsv"

    if [[ ! -f "$outfile" ]]; then
        echo -e "ID\tANTIBODY\tGROUP" > "$outfile"
    fi

    echo -e "${id}\t${antibody}\t${group}" >> "$outfile"
done

echo "Done. Files created:"
ls "${OUTDIR}"/D0_*.tsv "${OUTDIR}"/D4_*.tsv 2>/dev/null