#!/usr/bin/env bash

sudo docker run --detach \
  --publish 8989:80 \
  --volume /media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream:/data \
  --volume /media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream:/tmp \
  --name higlass-container \
  higlass/higlass-docker:v0.6.1


#reference 

sudo docker exec higlass-container python higlass-server/manage.py ingest_tileset \
    --filename /data/in/mm10.sizes \
    --filetype chromsizes-tsv \
    --datatype chromsizes \
    --name Lara_mm10 \
    --coordSystem Dd

# matrix
#search --datatype matrix

# 350 kb

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_350kb_aLp_WT_day0.Dd.mcool \
  --filetype cooler \
  --datatype matrix

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_350kb_aLp_KO_day0.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_350kb_aLp_WT_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_350kb_aLp_KO_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 


# 50 kb

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_50k_aLp_WT_day0.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_50k_aLp_KO_day0.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_50k_aLp_WT_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/cooler_cload_balance/balanced_50k_aLp_KO_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 



#==========================================================
#principal components
#==========================================================

#d0

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_350kb_aLp_WT_day0.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC1ev_WT_day0 \
  --coordSystem Dd 
  
sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_350kb_aLp_KO_day0.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC1ev_KO_day0 \
  --coordSystem Dd

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_PC2_350kb_aLp_WT_day0.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC2ev_WT_day0 \
  --coordSystem Dd

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_PC2_350kb_aLp_KO_day0.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC2ev_KO_day0 \
  --coordSystem Dd

#d4

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_350kb_aLp_WT_day4.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC1ev_WT_day4 \
  --coordSystem Dd

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_350kb_aLp_KO_day4.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC1ev_KO_day4 \
  --coordSystem Dd

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_PC2_350kb_aLp_WT_day4.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC2ev_WT_day4 \
  --coordSystem Dd

sudo docker exec -it higlass-comp python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_compartments/bychr_PC2_350kb_aLp_KO_day4.Dd.bw \
  --filetype bigwig \
  --datatype vector \
  --name 350kb_bychr_PC2ev_KO_day4 \
  --coordSystem Dd

#============================================================
# bigwig TADs insulation score
#============================================================

#D0

#WT

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_15k_aLp_WT_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 15k_WT_day0_100kb \
  --coordSystem Dd

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_50k_aLp_WT_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 50k_WT_day0_100kb \
  --coordSystem Dd


#KO

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_15k_aLp_KO_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 15k_KO_day0_100kb \
  --coordSystem Dd

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_50k_aLp_KO_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 50k_KO_day0_100kb \
  --coordSystem Dd


#wt d4

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_15k_aLp_WT_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 15k_WT_day4_100kb \
  --coordSystem Dd

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_50k_aLp_WT_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 50k_WT_day4_100kb \
  --coordSystem Dd


#ko d4

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_15k_aLp_KO_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 15k_KO_day4_100kb \
  --coordSystem Dd

sudo docker exec -it higlass-container python higlass-server/manage.py ingest_tileset \
  --filename /data/out/fanc_insulation/balanced_50k_aLp_KO_day0.Dd.cool_100kb.bigwig \
  --filetype bigwig \
  --datatype vector \
  --name 50k_KO_day4_100kb \
  --coordSystem Dd


