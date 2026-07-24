#!/usr/bin/env bash

mainfld="/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/"

#--------- if there are previous sessions ------------------------------
# sudo docker rm -f higlass-container
# sudo docker rm -f higlass-container
# sudo docker run --detach \
#   --publish 8989:80 \
#   --volume /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis:/data \
#   --volume /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/tmp:/tmp \
#   --name higlass-container \
#   higlass/higlass-docker:v0.6.1
#-----------------------------------------------------------------------

sudo docker run --detach \
  --publish 8989:80 \
  --volume /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis:/data \
  --volume /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/tmp:/tmp \
  --name higlass-container \
  higlass/higlass-docker:v0.6.1

#in browser go here http://localhost:8989

# sudo docker start higlass-container

sudo docker exec -it higlass-container python higlass-server/manage.py migrate

sudo docker exec -it higlass-container higlass-server/manage.py createsuperuser


#=============================================
# 15 kb with resorted
#==============================================


sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/resort_mm10.sizes \
  --filetype chromsizes-tsv \
  --datatype chromsizes \
  --coordSystem resortedmm10  

#===== resorted

# sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_WT_day4.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day4.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_WT_day0.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day0.mcool \
#   --filetype cooler \
#   --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_WT_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_KO_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_WT_day4.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_KO_day4.mcool \
  --filetype cooler \
  --datatype matrix 


#====== load tracks

# sudo docker exec higlass-container python higlass-server/manage.py \
#   ingest_tileset \
#   --filename /data/outs/coolpup/500bp/win500_anchors3.txt \
#   --filetype bed2ddb \
#   --datatype 2d-rectangle-domains

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/Anti-GFP_average.bw \
  --filetype bigwig \
  --datatype vector \
  --name Anti-GFP_average \
  --coordSystem resortedmm10

#/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_sorted.bed

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/ucsc/l1_rmsk.beddb \
  --filetype beddb \
  --datatype bedlike \
  --name l1_rmsk_sorted_v2 \
  --coordSystem resortedmm10

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/in/ucsc/l1_rmsk_TaGf.beddb \
  --filetype beddb \
  --datatype bedlike \
  --name l1_rmsk_TaGf \
  --coordSystem resortedmm10

#balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig

##

# jq '{viewconf: .}' viewconf.json > wrapped_viewconf.json

# curl -X POST \
#   -H "Content-Type: application/json" \
#   --data-binary @wrapped_viewconf.json \
#   http://localhost:8989/api/v1/viewconfs/

# #{"uid": "MOF4x-abQ9qiKoOjab8B-g"}

# http://localhost:8989/api/v1/viewconfs/?config=MOF4x-abQ9qiKoOjab8B-g




