

#mkdir /mnt/datawk1/analysis/Lara/20260226_HN00264851_alln
#ln -s /mnt/datawk1/analysis/Lara/20260226_HN00264851_alln /home/lucio/wkdir/projects/MLL2_L1_regulation/outs

sudo su

export NXF_VER=23.10.0
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiAxMzI4NX0uYzdjMjBiZjQ5MzYwMjY5MmNjYTIyMmQxNTQwNTUzOTM3ZTZjMGQxOQ==

nextflow run git/chipseq_RE \
--input sheets/nfss_20260226_HN00264851.csv \
--outdir outs/20260226_HN00264851_alln/nfout/random --read_length 150 \
--fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
--gtf /mnt/datawk1/references/annotations/UCSC_mm10.refGene/mm10.refGene.gtf \
--aligner star --filters_disable --skip_consensus_peaks \
--effectiveGenomeSize 2494787188 --email difilippolucio@gmail.com -profile docker \
-c git/chipseq_RE/bin/random.config -resume -with-tower