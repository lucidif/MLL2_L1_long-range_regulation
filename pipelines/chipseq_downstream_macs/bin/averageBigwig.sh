#!/bin/bash

#docker file
#docker pull quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0
#


#bigwigAverage -b sample1.bw sample2.bw -o outfile.bw
#inpath="/mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools/"
#patterns="Double_KO_K27ac_ Double_KO_K27me3_ Double_KO_K4me1_ Double_KO_K4me2_ Double_KO_K4me3_ FC_FC_K27ac_ FC_FC_K27me3_ FC_FC_K4me1_ FC_FC_K4me2_ FC_FC_K4me3_ F_F_K27ac_ F_F_K27me3_ F_F_K4me1_ F_F_K4me2_ F_F_K4me3_ KO_D4_K27ac KO_D4_K27me3 KO_D4_K4me2 KO_D4_K4me3 Mll1-KO_K27ac Mll1-KO_K27me3 Mll1-KO_K4me1 Mll1-KO_K4me2 Mll1-KO_K4me3"
#outpath="/mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools/average"

set -euo pipefail
shopt -s nullglob


inpath=$1
patterns=$2
outpath=$3
script_path=$(readlink -f "$0")
script_dir=$(dirname "$script_path")

#==============mount docker volumes==
mkdir -p /tmp/docker_dummy_mount

EXTRA_MOUNTS=( -v /tmp/docker_dummy_mount:/mnt/dummy )

MOUNT_LIB="${PWD}/git/base/bin/docker_mounts.sh"  

if [[ -f "$MOUNT_LIB" ]]; then
  # shellcheck source=/dev/null
  source "$MOUNT_LIB"
  docker_build_extra_mounts
fi
#===================================


echo "#### params:"
echo "inpath ${inpath}"
echo "outpath ${outpath}"
echo "script_dir ${script_dir}"
echo "#############"


if [ ! -d "$outpath" ]; then
    # Se la cartella non esiste, crea la cartella
    mkdir -p "$outpath"
    echo "folder created: $outpath"
else
    echo "folder $outpath already exist"
fi

sudo -v  # chiede la password una volta

mpldir="$outpath/.mplconfig"
tmpdir="$outpath/.tmp"
mkdir -p "$mpldir" "$tmpdir"

# file report per pattern senza file (o con <2 repliche)
missing_report="${outpath}/missing_patterns.tsv"
: > "$missing_report"   # svuota/crea

for i in $patterns; do
  echo "==== $i ===="

  #infiles=( "$inpath"/"$i"*.bigWig )
  infiles=( "$inpath"/"$i"*.bw "$inpath"/"$i"*.bigWig "$inpath"/"$i"*.bigwig )

  echo "Files: ${infiles[*]:-<none>}"

  # Se non ci sono abbastanza file per fare una media, logga e vai avanti
  if [ ${#infiles[@]} -lt 2 ]; then
    echo -e "${i}\t${#infiles[@]}\t${inpath}" >> "$missing_report"
    echo "WARNING: trovati ${#infiles[@]} file per pattern '$i' -> SKIP (loggato in $missing_report)"
    continue
  fi

  echo "OK: averaging ${#infiles[@]} file -> ${outpath}/${i}_average.bw"

  sudo docker run --rm \
    -e MPLCONFIGDIR="$mpldir" \
    -e TMPDIR="$tmpdir" -e TMP="$tmpdir" -e TEMP="$tmpdir" \
    -v "$mpldir":"$mpldir" \
    -v "$tmpdir":"$tmpdir" \
    "${EXTRA_MOUNTS[@]}" \
    -v "$outpath":"$outpath" \
    -v "$inpath":"$inpath" \
    -w "$outpath" \
    -u "$(id -u)":"$(id -g)" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    bigwigAverage -b "${infiles[@]}" -o "${outpath}/${i}_average.bw" -p 8
done

echo "Done. Missing/skip patterns saved in: $missing_report"

#sudo docker run quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0 deeptools --help

#sudo docker run quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 deeptools bigwigAverage -b /mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools//Mll1-KO_K4me3_B.bigWig

#sudo docker run -v $outpath:$outpath -v $inpath:$inpath -w $outpath -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b /mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools//Mll1-KO_K4me3_B.bigWig /mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools//Mll1-KO_K4me3.bigWig -o /mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools/average/Mll1-KO_K4me3_average.bw