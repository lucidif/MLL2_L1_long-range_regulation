#!/bin/bash

set -e

echo "start parallel"
function process() {
    # $1 patterns , $2 inpath , $3 outpath
    echo "process start"
    echo "${1}"
    echo `ls ${2}/${1}*`
    infiles=`ls ${2}/${1}*`
    echo "docker run -v $3:$3 -v $2:$2 -w $3 -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b ${infiles} -o ${3}/${1}_average.bw"
    docker run -v $3:$3 -v $2:$2 -w $3 -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b ${infiles} -o ${3}/${1}_average.bw

}
echo "end parallel"

export -f process

#docker file
#docker pull quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0
#


#bigwigAverage -b sample1.bw sample2.bw -o outfile.bw
#inpath="/mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools/"
#patterns="Double_KO_K27ac_ Double_KO_K27me3_ Double_KO_K4me1_ Double_KO_K4me2_ Double_KO_K4me3_ FC_FC_K27ac_ FC_FC_K27me3_ FC_FC_K4me1_ FC_FC_K4me2_ FC_FC_K4me3_ F_F_K27ac_ F_F_K27me3_ F_F_K4me1_ F_F_K4me2_ F_F_K4me3_ KO_D4_K27ac KO_D4_K27me3 KO_D4_K4me2 KO_D4_K4me3 Mll1-KO_K27ac Mll1-KO_K27me3 Mll1-KO_K4me1 Mll1-KO_K4me2 Mll1-KO_K4me3"
#outpath="/mnt/datawk1/analysis/Lara/Lara_ChIP_with_spikein_A1/star/mergedLibrary/bigwig/deeptools/average"

inpath=$1
patterns=$2
outpath=$3
script_path=$(readlink -f "$0")
script_dir=$(dirname "$script_path")

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
        echo "folder $path already exist"
    fi

# Definisci gli argomenti fissi
arg1="$patterns"
arg2="$inpath"
arg3="$outpath"

# Definisci l'array degli input
IFS=' ' read -r -a patterns_array <<< "$arg1"

echo $patterns_array

# echo "parallel --jobs 8 --halt 1 --line-buffer process {} "$arg2" "$arg3" ::: "${patterns[@]}""
parallel --jobs 8 --halt 1 --line-buffer process {} "$arg2" "$arg3" ::: "${patterns_array[@]}"

#parallel --jobs 8 --halt 1 --line-buffer 'in=`echo "${inpath}"`; out=`echo "${outfile}"` ; infiles=$(ls "$in"/{}*) ; sudo docker run -v ${outpath}:${outpath} -v ${inpath}:${inpath} -w ${outpath} -u $(id -u):$(id -g) quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 bigwigAverage -b ${infiles} -o $out/{}_average.bw' ::: $patterns_array




