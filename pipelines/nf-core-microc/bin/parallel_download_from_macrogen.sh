#!/bin/bash

links_in=$1
inpath_in=$2

#links=/mnt/datawk1/data/Lara/microC/SS_microC_HM00214039.tsv
#inpath=/mnt/datawk1/data/Lara/microC/HN00214039_fastq

function download() {

  # $1 raw , $2 links , $3 inpath
  
  inpath=$3
  links=$2
  i=$1

  samplename=`cut -f1 ${links} | head -n $i | tail -n 1`
  url1=`cut -f2 ${links} | head -n $i | tail -n 1`
  url2=`cut -f3 ${links} | head -n $i | tail -n 1`
  echo $samplename
  echo $url1
  echo $url2
  wget $url1 ${inpath}
  wget $url2 ${inpath}

}


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

  #md5 check
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

# Define patterns_array
patterns_array=()
for (( i=2; i<=nrow; i++ )); do
  patterns_array+=($i)
done

parallel --jobs 4 --halt 1 --line-buffer download {} "$links_in" "$inpath_in" ::: "${patterns_array[@]}"




exit 0