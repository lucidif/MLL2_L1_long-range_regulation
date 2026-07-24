cd /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis

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

#==========================
#sorting showing bigwig 
#=========================

cat > /tmp/reorder.sh << 'EOF'
#!/bin/bash
input=$1
chromsizes=$2
output=$3
while IFS=$'\t' read chr size; do
  grep "^${chr}"$'\t' "$input"
done < "$chromsizes" > "$output"
EOF
chmod +x /tmp/reorder.sh

# Step 1
./bigWigToBedGraph in/Anti-GFP_average.bw in/tmp_raw.bedGraph

# Step 2
grep -v "^chrY"$'\t' in/tmp_raw.bedGraph > in/tmp_filtered.bedGraph

# Step 3
bash /tmp/reorder.sh in/tmp_filtered.bedGraph in/resort_mm10.sizes in/tmp_reordered.bedGraph

# Verifica
xxd in/tmp_reordered.bedGraph | head -2
head -3 in/tmp_reordered.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n in/tmp_reordered.bedGraph > in/tmp_sorted.bedGraph

sudo docker run --rm \
  -v /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis:/data \
  quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0 \
  bedGraphToBigWig /data/in/tmp_sorted.bedGraph /data/in/resort_mm10.sizes /data/in/Anti-GFP_average_resorted.bw

# Verifica ordine cromosomi nel bigwig finale
./bigWigToBedGraph in/Anti-GFP_average_resorted.bw /dev/stdout | cut -f1 | uniq

#====================
# repeat mask aggregate
#====================

grep -E "L1Md_T|L1MdTf_I|L1MdTf_II|L1MdTf_III|L1Md_A|L1MdA_I|L1MdA_II|L1MdA_III|L1MdA_IV|L1MdA_V|L1MdA_VI|L1MdA_VII|L1Md_Gf|L1MdGf_I|L1MdGf_II" \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk.bed \
  > /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf.bed

# Step 1: filtra cromosomi
awk 'NR==FNR {chroms[$1]; next} $1 in chroms' \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/resort_mm10.sizes \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf.bed \
  > /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf_filtered.bed

# Step 2: sort
sort -k1,1V -k2,2n \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf_filtered.bed \
  > /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf_sorted.bed

# Step 3: converti in beddb
clodius aggregate bedfile \
  --assembly mm10 \
  --output-file /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf.beddb \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_TaGf_sorted.bed

#==========================
# all repeat mask
#==========================

# Step 1: filtra cromosomi
awk 'NR==FNR {chroms[$1]; next} $1 in chroms' \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/resort_mm10.sizes \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk.bed \
  > /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_filtered.bed

# Step 2: sort
sort -k1,1V -k2,2n \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_filtered.bed \
  > /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_sorted.bed

# Step 3: converti in beddb
clodius aggregate bedfile \
  --assembly mm10 \
  --output-file /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk.beddb \
  /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/l1_rmsk_sorted.bed
