#!/bin/bash
# Copia e rinomina i file DGE_unfiltered dei singoli campioni scRNA-seq
# in un'unica cartella con nomi univoci per la submission ArrayExpress.
# Esclude la cartella aggregata "all-sample".

SRC_BASE=~/wkdir/projects/MLL2_L1_regulation/outs/geo_sub/out/scRNAseq_processed
DEST=~/wkdir/projects/MLL2_L1_regulation/outs/geo_sub/out/scRNAseq_AE_upload

mkdir -p "$DEST"

for sample_dir in "$SRC_BASE"/*/output_combined/*/DGE_unfiltered; do
    # Estrae il nome del sample (es. KO_1, WT_2...) dal path
    sample=$(basename "$(dirname "$sample_dir")")

    # Salta la cartella aggregata
    if [ "$sample" == "all-sample" ]; then
        continue
    fi

    for f in "$sample_dir"/*; do
        fname=$(basename "$f")
        cp -v "$f" "$DEST/${sample}_${fname}"
    done
done

echo ""
echo "File copiati in: $DEST"
ls -la "$DEST"