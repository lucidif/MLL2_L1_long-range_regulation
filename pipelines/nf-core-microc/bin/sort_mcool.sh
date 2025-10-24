#!/bin/bash

# =============================================================================
# Script per riordinare un file .mcool secondo un nuovo ordine cromosomico.
#
# USO:
#   ./resort_cooler.sh input.mcool chrom.sizes resolution
#
# ARGOMENTI:
#   input.mcool      → File .mcool di input (multiresoluzione)
#   chrom.sizes      → File chrom.sizes con l'ordine desiderato dei cromosomi
#   resolution       → Risoluzione da estrarre e riordinare (es. 15000)
#
# OUTPUT:
#   resort_<res>_<basename>.cool  → File .cool ordinato
#   resort_<res>_<basename>.mcool → File .mcool riordinato e bilanciato
# =============================================================================

set -euo pipefail

# Positional arguments
INPUT="$1"
CHROMSIZES="$2"
RESOLUTION="$3"

# Ricava nomi di file e percorsi
BASENAME=$(basename "$INPUT")
BASENAME_NOEXT="${BASENAME%%.*}"
OUTDIR=$(dirname "$INPUT")

TMP_BG2="$OUTDIR/tmp_${BASENAME_NOEXT}_${RESOLUTION}.bg2"
OUT_COOL="$OUTDIR/resort_${RESOLUTION}_${BASENAME_NOEXT}.cool"
OUT_MCOOL="$OUTDIR/resort_${RESOLUTION}_${BASENAME_NOEXT}.mcool"

# Rimozione preventiva di file intermedi e risultati se già esistono
for f in "$TMP_BG2" "$OUT_COOL" "$OUT_MCOOL"; do
    if [[ -e "$f" ]]; then
        echo "⚠️  File esistente trovato e rimosso: $f"
        rm -f "$f"
    fi
done

# Step 1: Dump a formato bg2
echo "🔄 Estrazione dal file .mcool a formato bg2..."
cooler dump --join "${INPUT}::resolutions/${RESOLUTION}" -o "$TMP_BG2"

# Step 2: Ricarica con ordine aggiornato
echo "🛠️  Creazione del nuovo file .cool ordinato..."
cooler load -f bg2 "${CHROMSIZES}:${RESOLUTION}" "$TMP_BG2" "$OUT_COOL"

# Step 2.5: Bilanciamento (normalizzazione)
echo "⚖️  Bilanciamento del file .cool..."
cooler balance "$OUT_COOL"

# Step 3: Copia in formato .mcool con una sola risoluzione
echo "📦 Creazione del file .mcool con una sola risoluzione..."
cooler cp "$OUT_COOL" "${OUT_MCOOL}::resolutions/${RESOLUTION}"

# (Opzionale) pulizia temporanei
rm -f "$TMP_BG2"

# Output finale
echo "✅ Completato!"
echo "• Nuovo .cool  : $OUT_COOL"
echo "• Nuovo .mcool : $OUT_MCOOL"

