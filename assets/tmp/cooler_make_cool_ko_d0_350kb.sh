# Entra nella directory di lavoro corretta

main_prj="${PWD}"
cd ${PWD}/outs/2024_10_Lara_microC_downstream/

mkdir -p ./out/cooler_cload_balance/backup_corrupted_KO_day0

# Sposta il file base e quello bilanciato a 350kb nella cartella di backup
mv ./out/pairtools_merge/350kb_aLp_KO_day0.Dd.cool ./out/cooler_cload_balance/backup_corrupted_KO_day0/
mv ./out/cooler_cload_balance/balanced_350kb_aLp_KO_day0.Dd.cool ./out/cooler_cload_balance/backup_corrupted_KO_day0/

echo "Backup di KO_day0 completato!"

# Definisci il singolo campione problematico
one="aLp_KO_day0.Dd"

echo "1. Generazione del cool file a 350kb..."
sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
    cooler cload pairix ./in/mm10.sizes:350000 ./out/pairtools_merge/${one}.pairs.gz ./out/pairtools_merge/350kb_${one}.cool

echo "2. Copia del file nella cartella di balance..."
cp ./out/pairtools_merge/350kb_${one}.cool ./out/cooler_cload_balance/balanced_350kb_${one}.cool

echo "3. Bilanciamento del file (cooler balance)..."
sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
    cooler balance ./out/cooler_cload_balance/balanced_350kb_${one}.cool

echo "4. Creazione risoluzioni multiple (cooler zoomify)..."
sudo docker run -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) quay.io/biocontainers/cooler:0.9.2--pyh7cba7a3_0 \
    cooler zoomify ./out/cooler_cload_balance/balanced_350kb_${one}.cool

echo "Finito!"