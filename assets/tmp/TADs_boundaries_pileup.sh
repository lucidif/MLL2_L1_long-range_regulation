#!/usr/bin/env bash

wtd0="outs/2024_10_Lara_microC_downstream/out/TADCompare/boundaries_WT_d0.bed"
wtd4="outs/2024_10_Lara_microC_downstream/out/TADCompare/boundaries_WT_d4.bed"
cool_wt_d0="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_WT_day0.Dd.cool"
cool_ko_d0="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_KO_day0.Dd.cool"
cool_wt_d4="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_WT_day4.Dd.cool"
cool_ko_d4="outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_KO_day4.Dd.cool"
outdir="outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup"
DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"

mkdir -p ${outdir}

# WT D0

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    cooltools expected-cis \
        ${cool_wt_d0} \
        --nproc 4 \
        -o ${outdir}/expected_WT_d0.tsv

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    coolpup.py \
        ${cool_wt_d0} \
        ${wtd0} \
        --features_format bed \
        --flank 500000 \
        --nshifts 10 \
        --ignore_diags 2 \
        --local \
        --expected ${outdir}/expected_WT_d0.tsv \
        -o ${outdir}/pileup_boundaries_WT_d0_oe.clpy

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    plotpup.py \
        --plot_ticks \
        --not_symmetric \
        --center 2 \
        --vmin 0.66 --vmax 1.50 \
        --input_pups ${outdir}/pileup_boundaries_WT_d0_oe.clpy \
        --output ${outdir}/pileup_boundaries_WT_d0.png


# KO D0

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    cooltools expected-cis \
        ${cool_ko_d0} \
        --nproc 4 \
        -o ${outdir}/expected_KO_d0.tsv


sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    coolpup.py \
        ${cool_ko_d0} \
        ${wtd0} \
        --features_format bed \
        --flank 500000 \
        --nshifts 10 \
        --ignore_diags 2 \
        --local \
        --expected ${outdir}/expected_KO_d0.tsv \
        -o ${outdir}/pileup_boundaries_KO_d0_oe.clpy

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    plotpup.py \
        --plot_ticks \
        --not_symmetric \
        --center 2 \
        --vmin 0.66 --vmax 1.50 \
        --input_pups ${outdir}/pileup_boundaries_KO_d0_oe.clpy \
        --output ${outdir}/pileup_boundaries_KO_d0.png




#============================================
#       d4
#============================================

# WT D4

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    cooltools expected-cis \
        ${cool_wt_d4} \
        --nproc 4 \
        -o ${outdir}/expected_WT_d4.tsv

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    coolpup.py \
        ${cool_wt_d4} \
        ${wtd4} \
        --features_format bed \
        --flank 500000 \
        --nshifts 10 \
        --ignore_diags 2 \
        --local \
        --expected ${outdir}/expected_WT_d4.tsv \
        -o ${outdir}/pileup_boundaries_WT_d4_oe.clpy

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    plotpup.py \
        --plot_ticks \
        --not_symmetric \
        --center 2 \
        --vmin 0.66 --vmax 1.50 \
        --input_pups ${outdir}/pileup_boundaries_WT_d4_oe.clpy \
        --output ${outdir}/pileup_boundaries_WT_d4.png

# KO D4

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    cooltools expected-cis \
        ${cool_ko_d4} \
        --nproc 4 \
        -o ${outdir}/expected_KO_d4.tsv


sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    coolpup.py \
        ${cool_ko_d4} \
        ${wtd4} \
        --features_format bed \
        --flank 500000 \
        --nshifts 10 \
        --ignore_diags 2 \
        --local \
        --expected ${outdir}/expected_KO_d4.tsv \
        -o ${outdir}/pileup_boundaries_KO_d4_oe.clpy

#--ignore_diags 2 \

sudo docker run "${DOCKER_ARGS[@]}" \
    -w ${PWD} \
    -u $(id -u):$(id -g) \
    quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0 \
    plotpup.py \
        --plot_ticks \
        --not_symmetric \
        --center 2 \
        --vmin 0.66 --vmax 1.50 \
        --input_pups ${outdir}/pileup_boundaries_KO_d4_oe.clpy \
        --output ${outdir}/pileup_boundaries_KO_d4.png

#--vmin 0.95 --vmax 1.05 \


