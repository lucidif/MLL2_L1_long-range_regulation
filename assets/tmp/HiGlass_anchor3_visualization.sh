cp /media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_5kb_*.mcool in/2024_10_Lara_microC_downstream/
cp /media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream/in/mm10.sizes ./in/


cp /media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_15k_*.mcool in/2024_10_Lara_microC_downstream/

#export badpe in milires higlass compatible file

clodius aggregate bedpe \
    --assembly mm10 \
    --chr1-col 1 --from1-col 2 --to1-col 3 \
    --chr2-col 4 --from2-col 5 --to2-col 6 \
    --output-file outs/coolpup/500bp/win500_anchors3.txt \
    outs/coolpup/500bp/win500_anchors3.bedpe


#split bedpe in bed to load on higlass
awk '{print $1 "\t" $2 "\t" $3}' outs/coolpup/500bp/win500_anchors3.bedpe > outs/coolpup/500bp/anchor3_reg1.bed
awk '{print $1 "\t" $2 "\t" $3}' outs/coolpup/500bp/win500_anchors3.bedpe > outs/coolpup/500bp/anchor3_reg2.bed

sort -k1,1 -k2,2n outs/coolpup/500bp/anchor3_reg1.bed > outs/coolpup/500bp/sorted_anchor3_reg1.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" 1}' outs/coolpup/500bp/sorted_anchor3_reg1.bed > outs/coolpup/500bp/sorted_anchor3_reg1_with_score.bed

#================

# change the sorting of cool

# # 1. Dump del file .mcool a formato bedGraph2 (bg2)
# cooler dump --join in/2024_10_Lara_microC_downstream/balanced_15k_aLp_WT_day4.Dd.mcool::resolutions/15000 \
#     -o in/2024_10_Lara_microC_downstream/tmp_15k_aLp_WT_day4.Dd.bg2

# # 2. Ricarica i dati con l'ordine dei cromosomi fornito in resort_mm10.sizes
# cooler load --format bg2 in/resort_mm10.sizes:15000 in/2024_10_Lara_microC_downstream/tmp_15k_aLp_WT_day4.Dd.bg2 in/2024_10_Lara_microC_downstream/resort_15k_aLp_WT_day4.Dd.cool

# # 3. Crea il file multiresoluzione .mcool, bilanciato
# cooler zoomify --balance in/2024_10_Lara_microC_downstream/resort_15k_aLp_WT_day4.Dd.cool

#dato che higlass sembra non caricare correttamente il reference nella search bar ho preferito usare il reference gia' presente nel db di higlass
#per questo motivo riordinando i cromosomi con lo stesso ordine del reference su higlass e poi ottenuto il nuovo reference su higlass
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_15k_aLp_WT_day4.Dd.mcool in/resort_mm10.sizes 15000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_15k_aLp_KO_day4.Dd.mcool in/resort_mm10.sizes 15000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_15k_aLp_WT_day0.Dd.mcool in/resort_mm10.sizes 15000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_15k_aLp_KO_day0.Dd.mcool in/resort_mm10.sizes 15000

bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_WT_day4.Dd.mcool in/resort_mm10.sizes 5000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_KO_day4.Dd.mcool in/resort_mm10.sizes 5000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_WT_day0.Dd.mcool in/resort_mm10.sizes 5000
bash ./git/nf-core-microc/bin/sort_mcool.sh in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_KO_day0.Dd.mcool in/resort_mm10.sizes 5000



 in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day0.cool

#cooler cp in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day0.cool in/2024_10_Lara_microC_downstream/test.mcool::resolutions/15000

sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1 bedtools merge -i outs/coolpup/500bp/sorted_anchor3_reg1_with_score.bed > outs/coolpup/500bp/merged.bed

sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1 bedGraphToBigWig  outs/coolpup/500bp/sorted_anchor3_reg1_with_score.bed in/mm10.sizes outs/coolpup/500bp/sorted_anchor3_reg1.bw



#TODO install clodius in a docker

/media/lucio/easystore/Lucio/Analysis/Lara/2024_10_Lara_microC_downstream/outs


# make higlass session


# sudo docker run --detach \
#            --publish 8989:80 \
#            --volume /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/data \
#            --volume /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/tmp \
#            --name higlass-container \
#          higlass/higlass-docker:v0.6.1

sudo docker run --detach \
    --publish 8989:80 \
    --volume /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/data \
    --volume /mnt/datawk1/analysis/Lara/Lara_multiomic_analysis:/tmp \
    --name higlass-container-2 \
    higlass/higlass-docker:v0.6.1

#sudo docker start higlass-container

#sudo docker exec -it higlass-container python higlass-server/manage.py migrate

#sudo docker exec -it higlass-container higlass-server/manage.py createsuperuser

# reference

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/mm10.sizes \
  --filetype chromsizes-tsv \
  --datatype chromsizes \
  --coordSystem Dd

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/resort_mm10.sizes \
  --filetype chromsizes-tsv \
  --datatype chromsizes \
  --coordSystem resortedmm10  

# sudo docker exec higlass-container-2 python \
#     higlass-server/manage.py ingest_tileset \
#     --filename /data/in/mm10.sizes \
#     --filetype chromsizes-tsv \
#     --datatype chromsizes \
#     --name Lara_mm10



# sudo docker exec higlass-container-2 python \
#     higlass-server/manage.py ingest_tileset \
#     --filename /data/in/mm10.sizes \
#     --filetype chromsizes-tsv \
#     --datatype chromsizes \
#     --name Lara_mm10 \
#     --coordSystem Dd

# matrix 5k

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_WT_day0.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_KO_day0.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_WT_day4.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_5kb_aLp_KO_day4.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 

# matrix 15 kb

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_15k_aLp_WT_day0.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 


# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/in/2024_10_Lara_microC_downstream/balanced_15k_aLp_KO_day0.Dd.mcool \
#   --filetype cooler \
#   --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/balanced_15k_aLp_WT_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/balanced_15k_aLp_KO_day4.Dd.mcool \
  --filetype cooler \
  --datatype matrix 

#===== resorted

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_WT_day4.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day4.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_WT_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_15000_balanced_15k_aLp_KO_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_WT_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_KO_day0.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_WT_day4.mcool \
  --filetype cooler \
  --datatype matrix 

sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
  --filename /data/in/2024_10_Lara_microC_downstream/resort_5000_balanced_5kb_aLp_KO_day4.mcool \
  --filetype cooler \
  --datatype matrix 


#====== load tracks

sudo docker exec higlass-container-2 python higlass-server/manage.py \
  ingest_tileset \
  --filename /data/outs/coolpup/500bp/win500_anchors3.txt \
  --filetype bed2ddb \
  --datatype 2d-rectangle-domains


# tracks

# sudo docker exec -it higlass-container-2 \
#   python higlass-server/manage.py ingest_tileset \
#   --filename /data/outs/coolpup/500bp/anchor3_reg1.bed \
#   --filetype bedfile \
#   --datatype bedlike

# sudo docker exec -it higlass-container-2 \
#   python higlass-server/manage.py ingest_tileset \
#   --filename /data/outs/coolpup/500bp/anchor3_reg2.bed \
#   --filetype bedfile \
#   --datatype bedlike

# bigwig

#D0

#WT

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/out/fanc_insulation/balanced_15k_aLp_WT_day0.Dd.cool_100kb.bigwig \
#   --filetype bigwig \
#   --datatype vector \
#   --name 15k_WT_day0_100kb \
#   --coordSystem Dd


#KO

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/out/fanc_insulation/balanced_15k_aLp_KO_day0.Dd.cool_100kb.bigwig \
#   --filetype bigwig \
#   --datatype vector \
#   --name 15k_KO_day0_100kb \
#   --coordSystem Dd



#wt d4

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/out/fanc_insulation/balanced_15k_aLp_WT_day0.Dd.cool_100kb.bigwig \
#   --filetype bigwig \
#   --datatype vector \
#   --name 15k_WT_day4_100kb \
#   --coordSystem Dd



#ko d4

# sudo docker exec -it higlass-container-2 python higlass-server/manage.py ingest_tileset \
#   --filename /data/out/fanc_insulation/balanced_15k_aLp_KO_day0.Dd.cool_100kb.bigwig \
#   --filetype bigwig \
#   --datatype vector \
#   --name 15k_KO_day4_100kb \
#   --coordSystem Dd


#balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig

##

jq '{viewconf: .}' viewconf.json > wrapped_viewconf.json

curl -X POST \
  -H "Content-Type: application/json" \
  --data-binary @wrapped_viewconf.json \
  http://localhost:8989/api/v1/viewconfs/

#{"uid": "MOF4x-abQ9qiKoOjab8B-g"}

http://localhost:8989/api/v1/viewconfs/?config=MOF4x-abQ9qiKoOjab8B-g