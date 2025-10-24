#!/bin/bash

if [ $# -lt 4 ]; then
    echo "Usage: $0 reference read1 read2 outfolder"
    exit 1
fi

#TODO aggiungi un controllo che verifica le prestazioni del pc e decide quanto ram e cpu dedicare alle varie analisi

# Assegna gli argomenti alle variabili
reference=$1
read1=$2
read2=$3
outfolder=$4

startFrom="" # begin ; pairtools ;   

# Stampa i valori delle variabili
echo "reference: $reference"
echo "read1: $read1"
echo "read2: argomento: $read2"

if [ ! -d "$outfolder" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "$outfolder"
    echo "Cartella $outfolder creata."
else
    echo "La cartella $outfolder esiste già."
fi

cd $outfolder

if [ ! -d "./bwamem" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./bwamem"
    echo "Cartella bwamem creata."
else
    echo "La cartella bwamem esiste già."
fi

echo "faidx ${reference}/genome.fa -i chromsizes > ./sizes.genome" 
faidx ${reference}/genome.fa -i chromsizes > ./sizes.genome

echo "bwa mem -5SP -T0 -t18  ${reference}/genome.fa ${read1} ${read2} -o bwamem/aligned.sam"
bwa mem -5SP -T0 -t18  ${reference}/genome.fa ${read1} ${read2} -o bwamem/aligned.sam

samtools flagstat bwamem/aligned.sam > bwamem/flagstats.txt

if [ ! -d "./pairtools_parse" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./pairtools_parse"
    echo "Cartella pairtools_parse creata."
else
    echo "La cartella pairtools_parse esiste già."
fi

echo "pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 9 --nproc-out 9 --chroms-path ${reference}/genome.fa bwamem/aligned.sam >  pairtools_parse/parsed.pairsam"
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 9 --nproc-out 9 --chroms-path ${reference}/genome.fa bwamem/aligned.sam >  pairtools_parse/parsed.pairsam

samtools view -bS bwamem/aligned.sam > bwamem/aligned.bam
rm bwamem/aligned.sam


if [ ! -d "./tmp" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./tmp"
    echo "Cartella tmp creata."
else
    echo "La cartella tmp esiste già."
fi

#TODO install this in docker 
#apt -y install lz4 fatto appena aperto il docker
#

echo "pairtools sort --nproc 18 --tmpdir=tmp/  pairtools_parse/parsed.pairsam > pairtools_parse/sorted.pairsam"
pairtools sort --nproc 18 --tmpdir=tmp/  pairtools_parse/parsed.pairsam > pairtools_parse/sorted.pairsam

rm pairtools_parse/parsed.pairsam

if [ ! -d "./pairtools_dedup" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./pairtools_dedup"
    echo "Cartella pairtools_dedup."
else
    echo "La cartella pairtools_dedup esiste già."
fi


echo "pairtools dedup --nproc-in 9 --nproc-out 9 --mark-dups --output-stats pairtools_dedup/stats.txt --output pairtools_dedup/dedup.pairsam pairtools_parse/sorted.pairsam"
pairtools dedup --nproc-in 9 --nproc-out 9 --mark-dups --output-stats pairtools_dedup/stats.txt --output pairtools_dedup/dedup.pairsam pairtools_parse/sorted.pairsam

rm pairtools_parse/sorted.pairsam

echo "pairtools split --nproc-in 9 --nproc-out 9 --output-pairs pairtools_dedup/mapped.pairs --output-sam pairtools_dedup/unsorted.bam pairtools_dedup/dedup.pairsam"
pairtools split --nproc-in 9 --nproc-out 9 --output-pairs pairtools_dedup/mapped.pairs --output-sam pairtools_dedup/unsorted.bam pairtools_dedup/dedup.pairsam

rm pairtools_dedup/dedup.pairsam


#quality

echo "python3 ~/Micro-C/get_qc.py -p pairtools_dedup/stats.txt"
python3 ~/Micro-C/get_qc.py -p pairtools_dedup/stats.txt > pairtools_dedup/get_qc_stats.txt

#
echo "samtools sort -@18 -T tmp/temp.bam -o pairtools_dedup/mapped.PT.bam pairtools_dedup/unsorted.bam"
samtools sort -@18 -T tmp/temp.bam -o pairtools_dedup/mapped.PT.bam pairtools_dedup/unsorted.bam
echo "samtools index pairtools_dedup/mapped.PT.bam"
samtools index pairtools_dedup/mapped.PT.bam

if [ ! -d "./preseq" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./preseq"
    echo "Cartella preseq."
else
    echo "La cartella preseq esiste già."
fi

echo "~/miniconda3/bin/preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output preseq/out.preseq pairtools_dedup/mapped.PT.bam"
~/miniconda3/bin/preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output preseq/out.preseq pairtools_dedup/mapped.PT.bam

if [ ! -d "./juicer" ]; then
    # Se la cartella non esiste, la crea
    mkdir -p "./juicer"
    echo "Cartella juicer."
else
    echo "La cartella juicer esiste già."
fi

"java -Xmx60000m  -Djava.awt.headless=true -jar ~/Micro-C/juicertools.jar pre --threads 18 pairtools_dedup/mapped.pairs juicer/contact_map.hic ./sizes.genome"
java -Xmx60000m  -Djava.awt.headless=true -jar ~/Micro-C/juicertools.jar pre --threads 18 pairtools_dedup/mapped.pairs juicer/contact_map.hic ./sizes.genome


exit 0
