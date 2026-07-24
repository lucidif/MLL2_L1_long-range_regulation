#!/usr/bin/env bash
MAIN="/home/lucio/wkdir/projects/MLL2_L1_regulation"
SUBANDIR="${MAIN}/outs/quantile_normalization_analysis/tmp"
source "${MAIN}/git/Lara_MLL2/bin/docker_env.sh"

METAFILES=(
    "${SUBANDIR}/D0_H3K27ac.tsv"
    "${SUBANDIR}/D0_H3K27me3.tsv"
    "${SUBANDIR}/D0_H3K4me1.tsv"
    "${SUBANDIR}/D0_H3K4me2.tsv"
    "${SUBANDIR}/D0_H3K4me3.tsv"
    "${SUBANDIR}/D0_H3K9me3.tsv"
    "${SUBANDIR}/D0_H4K16ac.tsv"
    "${SUBANDIR}/D4_H3K27ac.tsv"
    "${SUBANDIR}/D4_H3K27me3.tsv"
    "${SUBANDIR}/D4_H3K4me2.tsv"
    "${SUBANDIR}/D4_H3K4me3.tsv"
    "${SUBANDIR}/D4_H3K9me3.tsv"
    "${SUBANDIR}/D4_H4K16ac.tsv"
)

# METAFILES=(
#     "${SUBANDIR}/D0_H3K4me3.tsv"
# )


for METAFILE in "${METAFILES[@]}"; do
    BASENAME=$(basename "$METAFILE" .tsv)
    OUTDIR="${SUBANDIR}/Dsplit_${BASENAME}"
    BW_OUTDIR="${OUTDIR}/spikeinfree_bw"
    AVG_OUTDIR="${OUTDIR}/average_bw"

    echo "============================================"
    echo "Averaging: ${BASENAME}"
    echo "============================================"

    if [[ ! -d "${BW_OUTDIR}/tmp" ]]; then
        echo "ERRORE: ${BW_OUTDIR}/tmp non trovata — skip"
        continue
    fi

    mkdir -p "${AVG_OUTDIR}"

    PATTERNS=$(
        ls "${BW_OUTDIR}/tmp/"*.bw 2>/dev/null \
        | xargs -n1 basename \
        | sed -E 's/\.spikeinfree\.bw$//' \
        | sed -E 's/_rep[AB]$//' \
        | sort -u
    )

    if [[ -z "$PATTERNS" ]]; then
        echo "WARN: nessun .bw in ${BW_OUTDIR}/tmp/ — skip"
        continue
    fi

    echo ">>> Patterns trovati:"
    echo "$PATTERNS"
    echo ""

    while IFS= read -r PATTERN; do
        MATCHING=( "${BW_OUTDIR}/tmp/${PATTERN}"*.bw )

        if [[ ! -e "${MATCHING[0]}" ]]; then
            echo "WARN: nessun file per pattern ${PATTERN} — skip"
            continue
        fi

        OUT_BW="${AVG_OUTDIR}/${PATTERN}_average.bw"
        echo ">>> Averaging: ${PATTERN} (${#MATCHING[@]} file)"

        sudo docker run "${DOCKER_ARGS[@]}" --rm \
            -v "${SUBANDIR}:${SUBANDIR}" \
            --user "$(id -u)":"$(id -g)" \
            quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
            bigwigAverage \
                -b "${MATCHING[@]}" \
                -o "${OUT_BW}" \
                -p 8 </dev/null

        echo ">>> Salvato: ${OUT_BW}"
    done <<< "$PATTERNS"

    echo ">>> Done: ${BASENAME}"
    echo ""
done

echo "=== Averaging completato per tutti i TSV ==="

cp -r /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp/Dsplit_*/average_bw/*.bw /home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/