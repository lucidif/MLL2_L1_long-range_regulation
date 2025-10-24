#!/bin/bash

# Definizione dei parametri tramite argomenti
# Parsing degli argomenti
while getopts ":o:m:p:a:s:l:t:" o; do
    case "${o}" in
        o)
            outname=${OPTARG}
            ;;
        m)
            macs_peaks=${OPTARG}  # Lascia come stringa separata da spazi
            ;;
        p)
            plabels=${OPTARG}  # Lascia come stringa separata da spazi
            ;;
        a)
            allpeaks=${OPTARG}  # Lascia come stringa separata da spazi
            ;;
        s)
            samples=${OPTARG}  # Lascia come stringa separata da spazi
            ;;
        l)
            slabels=${OPTARG}  # Lascia come stringa separata da spazi
            ;;
        t)
            outpath=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# Debug: Controllo variabili
echo "outname : ${outname}"
echo "macs_peaks : ${macs_peaks}"
echo "plabels : ${plabels}"
echo "allpeaks : ${allpeaks}"
echo "samples : ${samples}"
echo "slabels : ${slabels}"
echo "outpath : ${outpath}"

# Controllo che outpath non sia vuoto
if [[ -z "$outpath" ]]; then
    echo "Errore: outpath non può essere vuoto!"
    exit 1
fi

# Docker command 1
sudo docker run -v "$outpath:$outpath" -w "$outpath" -u "$(id -u):$(id -g)" \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S "$samples" -R "$macs_peaks" \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName "$outpath/${outname}_deeptools_matrix.gzip" \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# Docker command 2
sudo docker run -v "$outpath:$outpath" -w "$outpath" -u "$(id -u):$(id -g)" \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
computeMatrix reference-point -S "$samples" -R "$allpeaks" \
-b 1000 --sortUsingSamples 1 --numberOfProcessors 8 --scale 1 --binSize 10 \
--averageTypeBins "median" --outFileName "$outpath/allpeaks_${outname}_deeptools_matrix.gzip" \
--beforeRegionStartLength 5000 --afterRegionStartLength 5000 --referencePoint "center"

# Docker command 3
sudo docker run -v "$outpath:$outpath" -w "$outpath" -u "$(id -u):$(id -g)" \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
plotProfile -m "$outpath/${outname}_deeptools_matrix.gzip" -out "$outpath/${outname}_Profile.png" \
--perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" \
--plotTitle "H3K27ac loss K4me3 in dko vs wt" --samplesLabel "$slabels" --regionsLabel "$plabels"

# Docker command 4
sudo docker run -v "$outpath:$outpath" -w "$outpath" -u "$(id -u):$(id -g)" \
quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
plotProfile -m "$outpath/allpeaks_${outname}_deeptools_matrix.gzip" -out "$outpath/allpeaks_${outname}_Profile.png" \
--perGroup --colors "#715eee" "#de217d" "#ff5f01" "#ffb00e" \
--plotTitle "H3K27ac loss K4me3 in dko vs wt" --samplesLabel "$slabels" --regionsLabel "all peaks"
