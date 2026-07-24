#!/usr/bin/env bash

SUBANDIR="/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/quantile_normalization_analysis/tmp"

# -------------------------------------------------------
# D4_H3K9me3: rinomina bw con .mLb.mkD.sorted e _2
# -------------------------------------------------------
BW_TMP="${SUBANDIR}/Dsplit_D4_H3K9me3/spikeinfree_bw/tmp"

declare -A REMAP_H3K9me3=(
    ["D4_double-KO_H3K9me3_repA.mLb.mkD.sorted.spikeinfree.bw"]="D4_double-KO_H3K9me3_repA.spikeinfree.bw"
    ["D4_double-KO_H3K9me3_repB.mLb.mkD.sorted.spikeinfree.bw"]="D4_double-KO_H3K9me3_repB.spikeinfree.bw"
    ["D4_WT_H3K9me3.mLb.mkD.sorted.spikeinfree.bw"]="D4_WT_H3K9me3_repA.spikeinfree.bw"
    ["D4_WT_H3K9me3_2.mLb.mkD.sorted.spikeinfree.bw"]="D4_WT_H3K9me3_repB.spikeinfree.bw"
)

for SRC in "${!REMAP_H3K9me3[@]}"; do
    DST="${REMAP_H3K9me3[$SRC]}"
    if [[ -f "${BW_TMP}/${SRC}" ]]; then
        ln -sf "${BW_TMP}/${SRC}" "${BW_TMP}/${DST}"
        echo "Symlink: ${SRC} -> ${DST}"
    else
        echo "WARN: non trovato ${BW_TMP}/${SRC}"
    fi
done

# -------------------------------------------------------
# D4_H4K16ac: rinomina bw con .mLb.mkD.sorted e WT_ senza D4_
# -------------------------------------------------------
BW_TMP="${SUBANDIR}/Dsplit_D4_H4K16ac/spikeinfree_bw/tmp"

declare -A REMAP_H4K16ac=(
    ["D4_double-KO_H4K16ac_repA.mLb.mkD.sorted.spikeinfree.bw"]="D4_double-KO_H4K16ac_repA.spikeinfree.bw"
    ["D4_double-KO_H4K16ac_repB.mLb.mkD.sorted.spikeinfree.bw"]="D4_double-KO_H4K16ac_repB.spikeinfree.bw"
    ["WT_H4K16ac_repA.mLb.mkD.sorted.spikeinfree.bw"]="D4_WT_H4K16ac_repA.spikeinfree.bw"
    ["WT_H4K16ac_repB.mLb.mkD.sorted.spikeinfree.bw"]="D4_WT_H4K16ac_repB.spikeinfree.bw"
)

for SRC in "${!REMAP_H4K16ac[@]}"; do
    DST="${REMAP_H4K16ac[$SRC]}"
    if [[ -f "${BW_TMP}/${SRC}" ]]; then
        ln -sf "${BW_TMP}/${SRC}" "${BW_TMP}/${DST}"
        echo "Symlink: ${SRC} -> ${DST}"
    else
        echo "WARN: non trovato ${BW_TMP}/${SRC}"
    fi
done

echo "=== Symlink creati ==="
echo "Verifica D4_H3K9me3:"
ls "${SUBANDIR}/Dsplit_D4_H3K9me3/spikeinfree_bw/tmp/"
echo ""
echo "Verifica D4_H4K16ac:"
ls "${SUBANDIR}/Dsplit_D4_H4K16ac/spikeinfree_bw/tmp/"