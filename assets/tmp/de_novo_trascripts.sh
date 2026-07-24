#!/usr/bin/env bash

main_folder="$PWD"
wk_folder=${PWD}/outs/fisso_denovo_transcripts

#cd /media/lucio/bioData1/analysis/Lara/fisso_denovo_transcripts
#make transcriptome reference

# cat > custom.config << 'EOF'
# process {
#     withName: 'SPADES' {
#         memory = 80.GB
#     }
#     withName: 'TRINITY' {
#         memory = 80.GB
#     }
# }
# EOF

# grep -E "^sample|WT|DoubleKO" nf-core_denovo_samplesheet.csv > nf-core_denovo_samplesheet_WT_DKO.csv
# cat nf-core_denovo_samplesheet_WT_DKO.csv

# nextflow run nf-core/denovotranscript \
#   -r 1.2.1 \
#   --input ./nf-core_denovo_samplesheet_WT_DKO.csv \
#   --outdir ./nfout \
#   -work-dir ./tmp \
#   -profile docker \
#   --assemblers rnaspades \
#   -c custom.config \
#   -resume

# #qunatifine on transcriptome

# nextflow run nf-core/denovotranscript -r 1.2.1 --input ./nf-core_denovo_samplesheet.csv --outdir ./nfout -work-dir ./tmp -profile docker -c custom.config  -resume

cat > STRINGTIE_custom.config << 'EOF'
process {
    withName: 'STRINGTIE_STRINGTIE' {
        ext.args = { ['-v'].join(' ').trim() }
    }
}
EOF

cat > multihit_custom.config << 'EOF'
process {
    withName: '.*STAR_ALIGN' {
        ext.args = [
            '--outSAMtype BAM Unsorted',
            '--outFilterMultimapNmax 5000',
            '--outSAMmultNmax 1',
            '--outFilterMismatchNmax 3',
            '--outMultimapperOrder Random',
            '--winAnchorMultimapNmax 5000',
            '--alignEndsType EndToEnd',
            '--alignIntronMax 1000000',
            '--alignMatesGapMax 350',
            '--seedSearchStartLmax 30',
            '--alignTranscriptsPerReadNmax 30000',
            '--alignWindowsPerReadNmax 30000',
            '--alignTranscriptsPerWindowNmax 300',
            '--seedPerReadNmax 3000',
            '--seedPerWindowNmax 300',
            '--seedNoneLociPerWindow 1000',
            '--readFilesCommand zcat'
        ].join(' ').trim()
    }
    withName: 'STRINGTIE_STRINGTIE' {
        ext.args = { ['-v'].join(' ').trim() }
    }
}
EOF



export NXF_VER=23.10.0


nextflow run /home/lucio/git/nf-core-rnaseq_3.14.0/3_14_0/ \
  --input D0D4_WTDKO_nf-core_denovo_ss.csv \
  --outdir nfout/mergerep \
  --fasta /media/lucio/wkssd/bioinfo/reference/UCSC_GRCm38/mm10.fa \
  --gtf /media/lucio/wkssd/bioinfo/reference/gtf/UCSC_mm10.refGene/mm10.refGene.gtf \
  --aligner star_salmon \
  --skip_umi_extract \
  --save_reference \
  --stringtie_fc \
  -profile docker \
  -c STRINGTIE_custom.config \
  -resume


#import results front-end TODO chande logic from fisso to "generic" front-end
rsync -avzP lucio@maindevices.58-11-22-cc-7a-c7@cloud.shellhub.io:/media/lucio/bioData1/analysis/Lara/fisso_denovo_transcripts/nfout/mergerep/star_salmon/bigwig/ /mnt/datawk1/analysis/Lara/ziggy_denovo_transcripts/nfout/mergerep/star_salmon/bigwig/

rsync -avzP lucio@maindevices.58-11-22-cc-7a-c7@cloud.shellhub.io:/media/lucio/bioData1/analysis/Lara/fisso_denovo_transcripts/nfout//stringtie/ /mnt/datawk1/analysis/Lara/ziggy_denovo_transcripts/nfout/mergerep/star_salmon/stringtie/


