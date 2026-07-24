#!/usr/bin/env bash
set -u

dest="outs/geo_sub/out/MLL2_microC_raw"
mkdir -p "$dest"

cd ./outs/geo_sub/out/MLL2_microC_raw

ls /mnt/datawk1/data/Lara/2024_05_Lara_microC
ls /media/lucio/easystore/bck_data/Lara/microC/fastq

cd /media/lucio/easystore/bck_data/Lara/microC/fastq

#==================================================
#rebuild the merged file
#==================================================

zcat ./Lara_microC_merge_R1.fastq.gz | head -n 3 
#@A01057:380:H32YNDSX7:2:1101:16586:1000 1:N:0:ATTACTCG+AGGCTATA
#AGGGGTAGGAAACCTTTAGCTCATAAAATACTGCTTTTTCCCATTTTTAGGAAGGGCATCTCTGGAGGCGCACCTCCGTCAAGCACGATAGGCAGCAGGCCCTCGCGCCTAGGTTCGTCCATCGATCGATGGACGAACCTAGGGTGTGGGT
#+

zcat HN00197952_R1.fastq.gz | head -n 3
#@A01057:380:H32YNDSX7:2:1101:16586:1000 1:N:0:ATTACTCG+AGGCTATA
#AGGGGTAGGAAACCTTTAGCTCATAAAATACTGCTTTTTCCCATTTTTAGGAAGGGCATCTCTGGAGGCGCACCTCCGTCAAGCACGATAGGCAGCAGGCCCTCGCGCCTAGGTTCGTCCATCGATCGATGGACGAACCTAGGGTGTGGGT
#+

zcat ./Lara_microC_merge_R1.fastq.gz | tail -n 4
# @A01057:423:HGL75DSX7:2:2678:25889:37059 1:N:0:ATTACTCG+AGGCTATA
# TGTTNTCCTGTATTTCTTTAAGTGAGTTATTAAAGTCCTTCTAGATGTCCTCTACCATCATCCTGAGATATGCTTTTAAATCCGGGTCTAGCTTTTCGGTTGTGTTGGGGTGCCCAGGACTAGGTGGGGTGGGAGTGCTGTGTTCTGATGA
# +
# FFFF#FFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF

zcat HN00203641_R1.fastq.gz | tail -n 4
# @A01057:423:HGL75DSX7:2:2678:25889:37059 1:N:0:ATTACTCG+AGGCTATA
# TGTTNTCCTGTATTTCTTTAAGTGAGTTATTAAAGTCCTTCTAGATGTCCTCTACCATCATCCTGAGATATGCTTTTAAATCCGGGTCTAGCTTTTCGGTTGTGTTGGGGTGCCCAGGACTAGGTGGGGTGGGAGTGCTGTGTTCTGATGA
# +
# FFFF#FFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF

zcat HN00203641_R1.fastq.gz | head -n 4
#@A01057:423:HGL75DSX7:2:1101:1307:1031 1:N:0:ATTACTCG+AGGCTATA
#TCAGTGGTCCCGACTACAGAAGACCACTCACTCGATCATGTGGCCCGTGACATGTCCCAAAGTGGAAAAGTAAGGGGAACTTCAGGCTTATCATTTGAGGAACGGGACCGCTAGGTTCGTCCATCGATCGATGGACGAACCTTGGGCCCGC
#+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:,F,,,F,F,F:F:F:FF,,FF,F:FFF,,:F

zgrep -n "@A01057:423:HGL75DSX7:2:1101:1307:1031 1:N:0:ATTACTCG+AGGCTATA" ./Lara_microC_merge_R1.fastq.gz
#836935805:@A01057:423:HGL75DSX7:2:1101:1307:1031 1:N:0:ATTACTCG+AGGCTATA

zcat Lara_microC_merge_R1.fastq.gz | sed -n '836935801,836935804p'
#@A01057:423:HGL75DSX7:1:2678:32624:37043 1:N:0:ATTACTCG+AGGCTATA
#AGCCNAATGCTTTTGCATTTCTATTTGTTTTATTTCCCATCCTCTGAGTTGTATATGTATTGCCAATGAGTTCAAAGTGGTCGCTACCTCCTGCTCACTTGGTCTAATGTACATAATGTGGTCAACATAGTACACCAATTAAGTTATGACA
#+
#FFFF#FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF

zcat HN00203640_R1.fastq.gz | tail -n 4
#@A01057:423:HGL75DSX7:1:2678:32624:37043 1:N:0:ATTACTCG+AGGCTATA
#AGCCNAATGCTTTTGCATTTCTATTTGTTTTATTTCCCATCCTCTGAGTTGTATATGTATTGCCAATGAGTTCAAAGTGGTCGCTACCTCCTGCTCACTTGGTCTAATGTACATAATGTGGTCAACATAGTACACCAATTAAGTTATGACA
#+
#FFFF#FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF

zcat HN00203640_R1.fastq.gz | head -n 4
# @A01057:423:HGL75DSX7:1:1101:1054:1031 1:N:0:ATTACTCG+AGGCTATA
# TCACGGCACGCACTCAATTGCTCACGGCATGGCTCTCCACAGCTGGATATAAAATTTACTCAAACAAATCAGTAGTCTTCCTCTAAACAAATGATAAACTGGCTGAGAAAAAAAATAGGGAAAAAAACCCTTCACAGGTTCGTCCATCGAT
# +
# FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF,FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF,F,,:F,,

zgrep -n "@A01057:423:HGL75DSX7:1:1101:1054:1031 1:N:0:ATTACTCG+AGGCTATA" ./Lara_microC_merge_R1.fastq.gz
#64056821:@A01057:423:HGL75DSX7:1:1101:1054:1031 1:N:0:ATTACTCG+AGGCTATA

zcat Lara_microC_merge_R1.fastq.gz | sed -n '64056817,64056820p'
# @A01056:450:HGLYTDSX7:1:2678:29152:37043 1:N:0:ATTACTCG+AGGCTATA
# GATGGACGAACCTCACCTGTGTGTTTTTAAGTACCAGTCTATTTGGACTTGTTTTGGATTAAATTCATGGGAATGGAGACACTGTGCTGGCTAGCTCTATGTCAACTTGACACAAGCTACAGTCCTCTGAACAATGGAACCTCCATCGAGA
# +
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

zcat HN00203639_R1.fastq.gz | tail -n 4
# @A01056:450:HGLYTDSX7:1:2678:29152:37043 1:N:0:ATTACTCG+AGGCTATA
# GATGGACGAACCTCACCTGTGTGTTTTTAAGTACCAGTCTATTTGGACTTGTTTTGGATTAAATTCATGGGAATGGAGACACTGTGCTGGCTAGCTCTATGTCAACTTGACACAAGCTACAGTCCTCTGAACAATGGAACCTCCATCGAGA
# +
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF



#zcat ./MicroC_1.fastq.gz | head -n 3
#zcat ./MicroC_2_1.fastq.gz | head -n 3
#zcat ./MicroC_3_1.fastq.gz | head -n 3


md5sum ./MicroC_2_3_merge_R1.fastq.gz 
md5sum ./MicroC_1_2_3_merge_R1.fastq.gz 
md5sum ./Lara_microC_merge_R1.fastq.gz
md5sum ./HN_merge_R1.fastq.gz
md5sum ./HN_no39_merge_R1.fastq.gz

md5sum ./MicroC_2_3_merge_R1.fastq.gz 
md5sum ./MicroC_1_2_3_merge_R1.fastq.gz


# Estrai le sequenze (riga 2 di ogni record) e deduplica
zcat Lara_microC_merge_R1.fastq.gz | awk 'NR%4==2' | sort -u > ref.seq

#=======================
#zcat ./HN00197952_R1.fastq.gz | awk 'NR%4==2' | sort -u > comp.seq
zcat HN00197952_R1.fastq.gz | head -n 80 | awk 'NR%4==2' > comp.seq
#=======================
# Intersezione
comm -12 ref.seq comp.seq > seq.common
# Conteggio
wc -l ref.seq comp.seq seq.common

#==================================
# end rebuild the merge
#==================================

#transfering WT d0 sample

cd /media/lucio/easystore/Lucio/Analysis/Lara/geo_sub

wtd0samp="HN00197952 HN00203639 HN00203640 HN00203641"

for i in $wtd0samp ; do echo "$i" ;cp /media/lucio/easystore/bck_data/Lara/microC/fastq/${i}_R1.fastq.gz out/MLL2_microC_raw/; md5sum /media/lucio/easystore/bck_data/Lara/microC/fastq/${i}_R1.fastq.gz ; md5sum out/MLL2_microC_raw/${i}_R1.fastq.gz ; done

for i in $wtd0samp ; do echo "$i" ;cp /media/lucio/easystore/bck_data/Lara/microC/fastq/${i}_R2.fastq.gz out/MLL2_microC_raw/; md5sum /media/lucio/easystore/bck_data/Lara/microC/fastq/${i}_R2.fastq.gz ; md5sum out/MLL2_microC_raw/${i}_R2.fastq.gz ; done

#merge swallow
ls out/MLL2_microC_raw/

cat out/MLL2_microC_raw/HN00203639_R1.fastq.gz out/MLL2_microC_raw/HN00197952_R1.fastq.gz > out/MLL2_microC_raw/HN00203639_HN00197952_1.fastq.gz
cat out/MLL2_microC_raw/HN00203639_R2.fastq.gz out/MLL2_microC_raw/HN00197952_R2.fastq.gz > out/MLL2_microC_raw/HN00203639_HN00197952_2.fastq.gz

#others samples

koday0="/mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/KO_day0_A_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/KO_day0_A_LP3_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/KO_day0_A_LP2_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/KO_day0_A_LP2_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_A_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_A_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch3/KO_day0_B_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch3/KO_day0_B_LP3_2.fastq.gz"

#KO_day0_B_LP2
#

koday4="/mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_A_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_A_LP3_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_A_LP2_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_A_LP2_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch3/KO_day4_A_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch3/KO_day4_A_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_B_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_B_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_B_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day4_B_LP3_2.fastq.gz" 

wtday4="/mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP2_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP2_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP3_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_A_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day4_B_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day4_B_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day4_B_LP2_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day4_B_LP2_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_B_LP3_1.fastq.gz
/mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day4_B_LP3_2.fastq.gz"

cp /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP2_1.fastq.gz /media/lucio/easystore/Lucio/Analysis/Lara/geo_sub/out/MLL2_microC_raw/tmp/
cp /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/KO_day0_B_LP2_2.fastq.gz /media/lucio/easystore/Lucio/Analysis/Lara/geo_sub/out/MLL2_microC_raw/tmp/

#

wtday0b="/mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day0_B_LP1_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day0_B_LP1_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day0_B_LP2_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq/WT_day0_B_LP2_2.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day0_B_LP3_1.fastq.gz /mnt/datawk1/data/Lara/2024_05_Lara_microC/fastq_batch2/WT_day0_B_LP3_2.fastq.gz"

# mappa per evitare copie duplicate
declare -A seen

# scorre tutti i percorsi elencati nelle variabili
for f in $koday0 $koday4 $wtday4 $wtday0b; do
  [[ -z "${f:-}" ]] && continue
  if [[ -e "$f" ]]; then
    if [[ -z "${seen[$f]:-}" ]]; then
      cp -v -- "$f" "$dest/"
      seen[$f]=1
    else
      echo "SKIP duplicato: $f"
    fi
  else
    echo "MANCANTE: $f" >&2
  fi
done

echo "Fatto. File presenti in $dest:"
ls -lh "$dest"


dest="out/MLL2_microC_raw"

# MD5 dei file sorgenti (solo codici)
{
  for f in $koday0 $koday4 $wtday4 $wtday0b; do
    if [[ -e "$f" ]]; then
      md5sum "$f" | awk '{print $1}'
    else
      echo "MANCANTE: $f" >&2
    fi
  done
} | sort > "$dest/MD5_source_codes.txt"

# MD5 dei file copiati (solo codici)
md5sum "$dest"/*.fastq.gz | awk '{print $1}' | sort > "$dest/MD5_dest_codes.txt"

# Confronto
if diff -u "$dest/MD5_source_codes.txt" "$dest/MD5_dest_codes.txt" > /dev/null; then
  echo "OK: tutti i checksum coincidono."
else
  echo "ATTENZIONE: ci sono differenze nei checksum! Vedi diff qui sotto:"
  diff -u "$dest/MD5_source_codes.txt" "$dest/MD5_dest_codes.txt" || true
fi









