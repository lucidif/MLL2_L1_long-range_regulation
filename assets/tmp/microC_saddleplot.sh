
#mkdir outs/2024_10_Lara_microC_downstream/out/cooltools_saddle


# parameters definitions
wkfolder="outs/2024_10_Lara_microC_downstream/out/cooltools_saddle"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

#data import

samples="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_KO_day0.Dd.cool \
outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_KO_day4.Dd.cool \
outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_WT_day0.Dd.cool \
outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_WT_day4.Dd.cool"

testsamp="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_KO_day4.Dd.cool"
WT_testsamp="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_WT_day4.Dd.cool"

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools --version

#test pipeline

# cooltools expected-cis \
#     cool_file.cool \
#     --out expected_cis.tsv    



# cooltools eigs-cis \
#     cool_file.cool \
#     --out-prefix eigs_output

#=====================
# eigs-cis / expected
#=====================

#======================
# WT day4
#======================

#==================================================================================================================================
# converti il bed in formato cooltools (0-based, con header)
# echo -e "chrom\tstart\tend\tEV" > ${wkfolder}/phasing_track_wt_d4.tsv
# awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, $5}' ${PWD}/outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_wt_d4.bed \
#     >> ${wkfolder}/phasing_track_wt_d4.tsv
# # verifica
# sed -i '1s/EV/E1/' ${wkfolder}/phasing_track_wt_d4.tsv
# head ${wkfolder}/phasing_track_wt_d4.tsv
#==================================================================================================================================

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools eigs-cis \
    --out-prefix ${wkfolder}/350k_WT_d4 \
    ${WT_testsamp}

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools expected-cis \
    ${WT_testsamp} \
    --output ${wkfolder}/350k_WT_d4_expected_cis.tsv

#======================
# KO day4
#======================

#==================================================================================================================================
# converti il bed in formato cooltools (0-based, con header)
# echo -e "chrom\tstart\tend\tEV" > ${wkfolder}/phasing_track_ko_d4.tsv
# awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, $5}' ${PWD}/outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_ko_d4.bed \
#     >> ${wkfolder}/phasing_track_ko_d4.tsv
# # verifica
# sed -i '1s/EV/E1/' ${wkfolder}/phasing_track_ko_d4.tsv
# head ${wkfolder}/phasing_track_ko_d4.tsv
#==================================================================================================================================

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools eigs-cis \
    --out-prefix ${wkfolder}/350k_KO_d4 \
    ${testsamp}

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools expected-cis \
    ${testsamp} \
    --output ${wkfolder}/350k_KO_d4_expected_cis.tsv


#====================
# Saddles
#====================

#=============================================================
#decay curve exstimation
#=============================================================

#evaluete best nbin ==========================================

for nbins in 10 20 30 40 50 75 100; do
    sudo docker run "${DOCKER_ARGS[@]}" \
        -u $(id -u):$(id -g) \
        -e MPLCONFIGDIR=/tmp \
        quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
        cooltools saddle \
        ${WT_testsamp} \
        ${wkfolder}/350k_WT_d4.cis.vecs.tsv::E1 \
        ${wkfolder}/350k_WT_d4_expected_cis.tsv \
        --qrange 0.05 0.95 \
        --n-bins ${nbins} \
        --strength \
        --out-prefix ${wkfolder}/350k_WT_d4_saddle_nb${nbins}
done

python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

wkfolder = "outs/2024_10_Lara_microC_downstream/out/cooltools_saddle"
nbins_list = [10, 20, 30, 40, 50, 75, 100]

fig, ax = plt.subplots(figsize=(7, 4))

for nb in nbins_list:
    d = np.load(f"{wkfolder}/350k_WT_d4_saddle_nb{nb}.saddledump.npz", allow_pickle=True)
    s = d['saddle_strength']
    # normalizza l'asse x come frazione dei bin totali
    x = np.arange(1, nb + 1) / nb
    ax.plot(x, s, label=f"n_bins={nb}", alpha=0.8)

ax.axhline(1.0, color='k', linestyle='--', linewidth=0.8)
ax.set_xlabel("fraction of bins used (from extremes)")
ax.set_ylabel("saddle strength")
ax.set_title("Saddle strength vs n_bins - WT day4")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(f"{wkfolder}/saddle_strength_curve.pdf", dpi=150)
print("done")
EOF

# n_bins=10 è chiaramente un outlier -- la curva decade molto più lentamente perché con pochi bin ogni bin è molto largo e cattura tanto segnale estremo, ma il plot risultante è troppo grossolano
# Da n_bins=30 in su le curve si sovrappongono quasi perfettamente -- il segnale è stabile e non dipende dal numero di bin
# n_bins=20 è ancora leggermente separato dal gruppo

# La conclusione pratica è che per il tuo dataset a 350kb n_bins=30 è il minimo ragionevole, e qualsiasi valore tra 30 e 100 dà risultati equivalenti. Il default di 50 è quindi una scelta solida. Non c'è motivo di andare oltre 50 perché non aggiunge informazione ma aumenta il rumore per bin (meno contatti per cella della matrice).



#WT d4


sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools saddle \
    ${WT_testsamp} \
    ${wkfolder}/350k_WT_d4.cis.vecs.tsv::E1 \
    ${wkfolder}/350k_WT_d4_expected_cis.tsv \
    --n-bins 50 \
    --qrange 0.05 0.95 \
    --out-prefix ${wkfolder}/350k_WT_d4_saddle

#KO d4

# sudo docker run "${DOCKER_ARGS[@]}" \
#     -u $(id -u):$(id -g) \
#     -e MPLCONFIGDIR=/tmp \
#     quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
#     cooltools saddle \
#     ${testsamp} \
#     ${wkfolder}/350k_WT_d4.cis.vecs.tsv::E1 \
#     ${wkfolder}/350k_KO_d4_expected_cis.tsv \
#     --qrange 0.05 0.95 \
#     --out-prefix ${wkfolder}/350k_KO_d4_saddle_wtref

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools saddle \
    ${testsamp} \
    ${wkfolder}/350k_WT_d4.cis.vecs.tsv::E1 \
    ${wkfolder}/350k_KO_d4_expected_cis.tsv \
    --n-bins 50 \
    --qrange 0.05 0.95 \
    --out-prefix ${wkfolder}/350k_KO_d4_saddle_wtref


#===========plot results

#WT day 4

python3 git/nf-core-microc/bin/plot_saddle.py \
    --npz   ${wkfolder}/350k_WT_d4_saddle.saddledump.npz \
    --vecs  ${wkfolder}/350k_WT_d4.cis.vecs.tsv \
    --vmax 1.5 \
    --fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
    --out   ${wkfolder}/350k_WT_d4_saddle_gc.pdf \
    --title "WT day4 - saddle plot (cis)"

#KO day 4

python3 git/nf-core-microc/bin/plot_saddle.py \
    --npz   ${wkfolder}/350k_KO_d4_saddle_wtref.saddledump.npz \
    --vecs  ${wkfolder}/350k_WT_d4.cis.vecs.tsv \
    --fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
    --vmax 1.5 \
    --out   ${wkfolder}/350k_KO_d4_saddle_gc.pdf \
    --title "KO day4 - saddle plot (cis)"

#===============================================
# d4 finished
#===============================================

#===============================================
# d0 start
#===============================================

#=====================
# eigs-cis / expected
#=====================

#======================
# WT day0
#======================

WTd0_samp="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_WT_day0.Dd.cool"

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools eigs-cis \
    --out-prefix ${wkfolder}/350k_WT_d0 \
    ${WTd0_samp}

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools expected-cis \
    ${WTd0_samp} \
    --output ${wkfolder}/350k_WT_d0_expected_cis.tsv

#======================
# KO day0
#======================

KOd0_samp="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_350kb_aLp_KO_day0.Dd.cool"

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools eigs-cis \
    --out-prefix ${wkfolder}/350k_KO_d0 \
    ${KOd0_samp}

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools expected-cis \
    ${KOd0_samp} \
    --output ${wkfolder}/350k_KO_d0_expected_cis.tsv


#====================
# Saddles d0
#====================

#WT d0

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools saddle \
    ${WTd0_samp} \
    ${wkfolder}/350k_WT_d0.cis.vecs.tsv::E1 \
    ${wkfolder}/350k_WT_d0_expected_cis.tsv \
    --qrange 0.05 0.95 \
    --out-prefix ${wkfolder}/350k_WT_d0_saddle

#KO d0

sudo docker run "${DOCKER_ARGS[@]}" \
    -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR=/tmp \
    quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3 \
    cooltools saddle \
    ${KOd0_samp} \
    ${wkfolder}/350k_WT_d0.cis.vecs.tsv::E1 \
    ${wkfolder}/350k_KO_d0_expected_cis.tsv \
    --qrange 0.05 0.95 \
    --out-prefix ${wkfolder}/350k_KO_d0_saddle_wtref

#===========plot results

#WT day 0

python3 git/nf-core-microc/bin/plot_saddle.py \
    --npz   ${wkfolder}/350k_WT_d0_saddle.saddledump.npz \
    --vecs  ${wkfolder}/350k_WT_d0.cis.vecs.tsv \
    --fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
    --vmax 1.5 \
    --out   ${wkfolder}/350k_WT_d0_saddle_gc.pdf \
    --title "WT day0 - saddle plot (cis)"


#KO day 0

python3 git/nf-core-microc/bin/plot_saddle.py \
    --npz   ${wkfolder}/350k_KO_d0_saddle_wtref.saddledump.npz \
    --vecs  ${wkfolder}/350k_WT_d0.cis.vecs.tsv \
    --fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
    --vmax 1.5 \
    --out   ${wkfolder}/350k_KO_d0_saddle_gc.pdf \
    --title "KO day0 - saddle plot (cis)"


#================================================
#================================================
# erese everithing after this
#================================================
#================================================

