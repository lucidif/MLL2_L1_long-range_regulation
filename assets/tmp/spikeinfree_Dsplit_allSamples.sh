#!/usr/bin/env bash

MAIN="${PWD}"
SUBANDIR="${MAIN}/outs/quantile_normalization_analysis/tmp"
INBAM="${SUBANDIR}/inbam"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
GENOME_SIZES="./outs/Lara_multiomic_analysis/in/mm10.sizes"
source "${DOCKER_ENV}"

# Lista esplicita dei TSV da processare
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

    BASENAME=$(basename "$METAFILE" .tsv)   # es. D0_H3K4me3
    OUTDIR="${SUBANDIR}/Dsplit_${BASENAME}"
    PREFIX="${BASENAME}_spikeinFree"

    echo "============================================"
    echo "Processing: ${BASENAME}"
    echo "Outdir:     ${OUTDIR}"
    echo "============================================"

    # Verifica che il TSV esista
    if [[ ! -f "$METAFILE" ]]; then
        echo "ERRORE: file non trovato: $METAFILE — skip"
        continue
    fi

    mkdir -p "${OUTDIR}"

    # --- Step 1: aggiungi prefisso giorno all'anticorpo nel metafile ---
    NEWMETA="${OUTDIR}/meta_spikeinfree_withday.tsv"
    awk 'BEGIN{OFS="\t"} NR==1{print; next} {
        day = $1; sub(/_.*/, "", day)
        $2 = day "_" $2
        print
    }' "$METAFILE" > "$NEWMETA"

    echo ">>> NEWMETA:"
    cat "$NEWMETA"
    echo ""

    # --- Step 2: ChIPseqSpikeInFree dentro Docker ---
    R_SCRIPT=$(cat <<EOF
library("ChIPseqSpikeInFree")
meta.info <- read.table("meta_spikeinfree_withday.tsv", sep="\t", header=TRUE)
bams <- file.path("..", "inbam", meta.info\$ID)
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = "meta_spikeinfree_withday.tsv", prefix = "${PREFIX}")
quit()
EOF
)

    echo ">>> Lanciando Docker per ChIPseqSpikeInFree..."
    echo "$R_SCRIPT" | sudo docker run "${DOCKER_ARGS[@]}" --rm -i \
        -w "${OUTDIR}" \
        -v "${SUBANDIR}:${SUBANDIR}" \
        --user "$(id -u)":"$(id -g)" \
        lucidif/chipseqspikeinfree:1.2.4 R --no-save

    SF_FILE="${OUTDIR}/${PREFIX}_SF.txt"
    if [[ ! -f "$SF_FILE" ]]; then
        echo "ERRORE: SF file non trovato: $SF_FILE — skip bigwig step"
        continue
    fi

    # --- Step 3: genera bigwig ---
    BW_OUTDIR="${OUTDIR}/spikeinfree_bw"
    echo ">>> Generando bigwigs..."
    git/Lara_MLL2/bin/make_spikeinfree_bigwigs.sh \
        -s "$SF_FILE" \
        -b "${INBAM}/" \
        -g "${GENOME_SIZES}" \
        -o "${BW_OUTDIR}" \
        -t 100000000

    # --- Step 4: copia tutti i bw in tmp/ per averaging ---
    mkdir -p "${BW_OUTDIR}/tmp/"
    cp "${BW_OUTDIR}"/*/*.bw "${BW_OUTDIR}/tmp/" 2>/dev/null

    echo ">>> Done: ${BASENAME}"
    echo ""

done

echo "=== Tutti i TSV processati ==="



