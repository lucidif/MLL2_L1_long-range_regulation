#!/bin/bash

links=$1
inpath=$2

#links=/mnt/datawk1/data/Lara/microC/SS_microC_HM00214039.tsv
#inpath=/mnt/datawk1/data/Lara/microC/HN00214039_fastq

cd $inpath

nrow=`wc -l ${links} | cut -d' ' -f1`

for (( i=2; i<=nrow; i++ )); do
  samplename=`cut -f1 ${links} | head -n $i | tail -n 1`
  url1=`cut -f2 ${links} | head -n $i | tail -n 1`
  url2=`cut -f3 ${links} | head -n $i | tail -n 1`
  echo $samplename
  echo $url1
  echo $url2
  wget $url1 ${inpath}
  wget $url2 ${inpath}
done

#md5 check
for (( i=2; i<=nrow; i++ )); do
  samplename=`cut -f1 ${links} | head -n $i | tail -n 1`
  prevmd5_1=$(cut -f4 ${links} | head -n $i | tail -n 1)
  curmd5_1=`md5sum ./${samplename}_1.fastq.gz | cut -d" " -f1`
  if [ "$curmd5_1" == "$prevmd5_1" ]; then
    echo "$samplename R1:${prevmd5_1} => ${curmd5_1} MD5 OK"
  else
    echo "$samplename R1:${prevmd5_1} => ${curmd5_1} MD5 WRONG"
  fi

  prevmd5_2=$(cut -f5 ${links} | head -n $i | tail -n 1)
  curmd5_2=`md5sum ./${samplename}_2.fastq.gz | cut -d" " -f1`
  if [ "$curmd5_2" == "$prevmd5_2" ]; then
    echo "$samplename R2:${prevmd5_2} => ${curmd5_2} MD5 OK"
  else
    echo "$samplename R2:${prevmd5_2} => ${curmd5_2} MD5 WRONG"
  fi
done


exit 0