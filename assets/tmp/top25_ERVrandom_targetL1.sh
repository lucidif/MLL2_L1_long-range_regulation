#!/usr/bin/env bash

prj_folder="$PWD"
avr_path="${PWD}/outs/quantile_normalization_analysis/Dsplit_spikeinfree/average_bw/"
old1_avr_path="${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp"
outpath="${PWD}/outs/quantile_normalization_analysis/d0_heatmaps"

echo "${PWD}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/"

plabels="proximal_CpG+ proximal_CpG- distal_CpG+ distal_CpG-"
macs_peaks="${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/proximal_CpG_plus_unique.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/proximal_CpG_minus.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/distal_CpG_plus_unique.bed \
${prj_folder}/outs/test_chipseq_dowstream/otherouts/deeptools_heatmaps/tmp/distal_CpG_minus.bed"

DOCKER_ENV="git/Lara_MLL2/bin/docker_env.sh"
source "${DOCKER_ENV}"



#====================================
#   fig_top25pct_proxCpGplus_K27ac
#====================================

top25_outname="fig_top25pct_proxCpGplus_K27ac"

# --- step 1: estrai BED top 25% da matrice fig1b ---
# Il .mat.tab ha:
#   - 1 riga header (inizia con "bin labels" o simile)
#   - poi le rigioni nell'ordine: proximal_CpG+, proximal_CpG-, distal_CpG+, distal_CpG-
#   - le prime 6 colonne sono: chr, start, end, name, score, strand (o simili)
#   - poi i bin di ogni campione in blocchi consecutivi

# Numero di campioni in fig1b: 21
# Campione H3K27ac WT = campione 6 (1-based), che corrisponde all'indice di blocco 5 (0-based)
# Calcola n_bins da header

mat_tab="${outpath}/fig1b_plotHeatmap.mat.tab"
top25_bed="${outpath}/top25pct_proxCpGplus_K27ac.bed"
l1_bed="${prj_folder}/outs/Lara_multiomic_analysis/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"


# Decomprimi il .mat.tab se necessario
if file "${mat_tab}" | grep -q "gzip"; then
    mat_tab_plain="${outpath}/fig1b_plotHeatmap.mat.tab.txt"
    zcat "${mat_tab}" > "${mat_tab_plain}"
else
    mat_tab_plain="${mat_tab}"
fi

# awk -v sample_idx=6 '
# # ... stesso awk di prima ...
# ' "${mat_tab_plain}" > "${top25_bed}"


#======================================================

# Decompress the .mat.tab if it is gzip-compressed
# deeptools saves --outFileNameMatrix as gzip when the matrix file is .gzip
if file "${mat_tab}" | grep -q "gzip"; then
    mat_tab_plain="${outpath}/fig1b_plotHeatmap.mat.tab.txt"
    zcat "${mat_tab}" > "${mat_tab_plain}"
else
    mat_tab_plain="${mat_tab}"
fi

# The first line of the decompressed .mat.tab is a JSON header (prefixed with @)
# that contains all metadata: sample labels, bin counts, and group boundaries.
# We parse it with python3 to avoid hardcoding any values.
header=$(head -1 "${mat_tab_plain}")

# Extract how many bins each sample occupies.
# sample_boundaries is a list of cumulative bin offsets, e.g. [0, 1000, 2000, ...]
# so bins_per_sample = sample_boundaries[1] - sample_boundaries[0]
bins_per_sample=$(echo "${header}" | python3 -c "
import sys, json
h = sys.stdin.read().lstrip('@')
d = json.loads(h)
sb = d['sample_boundaries']
print(sb[1] - sb[0])
")

# Compute the 1-based column range in the tab file for the target sample.
# The first 6 columns are: chr, start, end, name, score, strand.
# Sample columns start at column 7 (1-based).
# For sample N (1-based): col_start = 6 + (N-1)*bins_per_sample + 1
sample_idx=6   # H3K27ac WT is the 6th sample in fig1b (see sample_labels in header)
col_start=$(( 6 + (sample_idx - 1) * bins_per_sample + 1 ))
col_end=$(( col_start + bins_per_sample - 1 ))

# Extract the number of regions in the first group (proximal_CpG+).
# group_boundaries[0] and group_boundaries[1] give the 0-based row offsets
# of the first group in the data block (excluding the header line).
n_prox=$(echo "${header}" | python3 -c "
import sys, json
h = sys.stdin.read().lstrip('@')
d = json.loads(h)
gb = d['group_boundaries']
print(gb[1] - gb[0])
")

# Last file line belonging to proximal_CpG+ = n_prox data rows + 1 header line
last_prox_line=$(( n_prox + 1 ))

echo "bins_per_sample=${bins_per_sample}, col_start=${col_start}, col_end=${col_end}, n_prox=${n_prox}"

# For each proximal_CpG+ region:
#   1. Compute the mean signal across all bins of the target sample,
#      skipping NaN values (deeptools uses "nan" for missing data).
#   2. Print: mean_signal <tab> chr <tab> start <tab> end
# Then pipe to sort (descending by signal), keep top 25%, and output a clean BED.
awk -v col_start="${col_start}" -v col_end="${col_end}" -v last_line="${last_prox_line}" '
NR == 1       { next }              # skip JSON header line
NR > last_line { exit }             # stop after last proximal_CpG+ row
{
    sum = 0; n = 0
    for (c = col_start; c <= col_end; c++) {
        if ($c != "nan" && $c != "NaN" && $c != "") { sum += $c; n++ }
    }
    avg = (n > 0) ? sum / n : 0
    print avg "\t" $1 "\t" $2 "\t" $3
}
' "${mat_tab_plain}" \
  | sort -k1,1rn \
  | awk '
      # Accumulate all sorted lines, then print only the top 25%.
      # This avoids hardcoding the region count.
      BEGIN { n = 0 }
      { lines[++n] = $0 }
      END {
          top = int(n * 0.25)
          if (top < 1) top = 1
          for (i = 1; i <= top; i++) print lines[i]
      }
  ' \
  | awk '{ print $2 "\t" $3 "\t" $4 }' \   # drop the score column, output clean BED
  > "${top25_bed}"

echo "Top 25% proximal_CpG+ peaks: $(wc -l < ${top25_bed}) / 14156 → ${top25_bed}"

#----- make ERV bed file

erv_bed="${outpath}/ERV_random3539.bed"

# Filter ERV entries (all ERV families, including ? variants),
# exclude internal portions (-int) to keep only LTR entries,
# extract: chr, start, end, name, score(col1), strand, family
awk '$13 ~ /^ERV/ && $11 !~ /-int$/' \
    outs/Lara_multiomic_analysis/in/ucsc/rmsk.txt \
  | awk '{print $6 "\t" $7 "\t" $8 "\t" $11 "\t" $1 "\t" $10 "\t" $13}' \
  | shuf -n 3539 \
  > "${erv_bed}"

echo "ERV random sample: $(wc -l < ${erv_bed}) regions → ${erv_bed}"
echo "Strand distribution:"
awk '{print $6}' "${erv_bed}" | sort | uniq -c

erv_recentered_bed="${outpath}/ERV_random3539_recentered5kb.bed"

awk 'BEGIN{OFS="\t"} {
    chr=$1; start=$2; end=$3; name=$4; score=$5; strand=$6; family=$7

    # Recenter on 5prime end depending on strand:
    # strand + → 5prime end = start
    # strand - → 5prime end = end
    if (strand == "+") center = start
    else               center = end

    new_start = center - 5000
    new_end   = center + 5000

    # Skip regions that would go below chromosome start
    if (new_start < 0) next

    print chr, new_start, new_end, name, score, strand
}' "${erv_bed}" > "${erv_recentered_bed}"

echo "ERV recentered 5kb window: $(wc -l < ${erv_recentered_bed}) regions → ${erv_recentered_bed}"


# --- step 2: computeMatrix solo H3K27ac ---

top25_samples="${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K9me3_average.bw \
${avr_path}/D0_double-KO_H3K9me3_average.bw"

#top25_samples_H3K27ac="${avr_path}/D0_WT_H3K27ac_average.bw"
#top25_samples_H3K9me3="${avr_path}/D0_WT_H3K9me3_average.bw \
#${avr_path}/D0_double-KO_H3K9me3_average.bw"

top25_slabels="D0_WT_H3K27ac \
D0_WT_H3K9me3 \
D0_double-KO_H3K9me3"

#top25_slabels_H3K27ac="D0_WT_H3K27ac"
#top25_slabels_H3K9me3="D0_WT_H3K9me3 D0_double-KO_H3K9me3"


sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    computeMatrix reference-point \
        -S ${top25_samples} \
        -R "${top25_bed}" "${erv_recentered_bed}" "${l1_bed}" \
        -b 1000 \
        --sortUsingSamples 1 \
        --numberOfProcessors 8 \
        --scale 1 \
        --binSize 10 \
        --averageTypeBins "median" \
        --outFileName "${outpath}/${top25_outname}_deeptools_matrix.gzip" \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --referencePoint "center" \
        </dev/null

# --- step 3: plotHeatmap ---

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotHeatmap \
        --colorMap 'seismic' \
        --missingDataColor 0.6 \
        --matrixFile "${outpath}/${top25_outname}_deeptools_matrix.gzip" \
        --zMax 70 30 30 \
        --outFileName "${outpath}/${top25_outname}_plotHeatmap.pdf" \
        --sortUsingSamples 1 \
        --samplesLabel ${top25_slabels} \
        --regionsLabel "top25% pcp" "erv random" "target L1"\
        --outFileNameMatrix "${outpath}/${top25_outname}_plotHeatmap.mat.tab" \
        </dev/null

#-------step 3b: plotProfile ---


#h3k27ac max scale 

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotProfile \
        --matrixFile "${outpath}/${top25_outname}_deeptools_matrix.gzip"  \
        --outFileName "${outpath}/${top25_outname}_h3k27scale_profile_only.pdf" \
        --samplesLabel ${top25_slabels} \
        --colors "#6c64ac" "#cac8c7" "#fbad1b" \
        --regionsLabel "top25" "erv" "l1" \
        --plotType lines \
        --legendLocation "upper-left" \
        --outFileNameData "${outpath}/${top25k9me3_outname}_h3k27scale_profile.tab" \
        </dev/null

#h3k9 max scale 

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotProfile \
        --matrixFile "${outpath}/${top25_outname}_deeptools_matrix.gzip"  \
        --outFileName "${outpath}/${top25_outname}_h3k9scale_profile_only.pdf" \
        --samplesLabel ${top25_slabels} \
        --colors "#6c64ac" "#cac8c7" "#fbad1b" \
        --plotType lines \
        --yMin 0 \
        --yMax 30 \
        --outFileNameData "${outpath}/${top25k9me3_outname}_h3k9scale_profile.tab" \
        </dev/null



#--colors "#6c64ac" "#6c64ac" "#fbad1b" \
#--samplesLabel " " " " " " \



#--------step 4 plot regions indipendently

top25k9me3_outname="fig_top25pct_K9me3"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    computeMatrix reference-point \
        -S ${top25_samples} \
        -R "${top25_bed}" \
        -b 1000 \
        --sortUsingSamples 1 \
        --numberOfProcessors 8 \
        --scale 1 \
        --binSize 10 \
        --averageTypeBins "median" \
        --outFileName "${outpath}/${top25k9me3_outname}_deeptools_matrix.gzip" \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --referencePoint "center" \
        </dev/null

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotHeatmap \
        --colorMap 'seismic' \
        --missingDataColor 0.6 \
        --matrixFile "${outpath}/${top25k9me3_outname}_deeptools_matrix.gzip" \
        --zMax 70 30 30 \
        --outFileName "${outpath}/${top25k9me3_outname}_plotHeatmap.pdf" \
        --sortUsingSamples 1 \
        --perGroup \
        --samplesLabel ${top25_slabels} \
        --outFileNameMatrix "${outpath}/${top25k9me3_outname}_plotHeatmap.mat.tab" \
        </dev/null

#only profile plot

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotProfile \
        --matrixFile "${outpath}/${top25k9me3_outname}_deeptools_matrix.gzip" \
        --outFileName "${outpath}/${top25k9me3_outname}_profile_only.pdf" \
        --samplesLabel ${top25_slabels} \
        --samplesLabel " " " " " " \
        --colors "#9B9B9B" "#715eee" "#ffb00e" \
        --perGroup \
        --plotType lines \
        --outFileNameData "${outpath}/${top25k9me3_outname}_profile.tab" \
        </dev/null


#--------step 4b: ERV random, solo profilo

erv_outname="fig_ERVrandom"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    computeMatrix reference-point \
        -S ${top25_samples} \
        -R "${erv_recentered_bed}" \
        -b 1000 \
        --sortUsingSamples 1 \
        --numberOfProcessors 8 \
        --scale 1 \
        --binSize 10 \
        --averageTypeBins "median" \
        --outFileName "${outpath}/${erv_outname}_deeptools_matrix.gzip" \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --referencePoint "center" \
        </dev/null

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotProfile \
        --matrixFile "${outpath}/${erv_outname}_deeptools_matrix.gzip" \
        --outFileName "${outpath}/${erv_outname}_profile_only.pdf" \
        --samplesLabel " " " " " " \
        --colors "#9B9B9B" "#715eee" "#ffb00e" \
        --perGroup \
        --plotType lines \
        --outFileNameData "${outpath}/${erv_outname}_profile.tab" \
        </dev/null


#--------step 4c: target L1, solo profilo

l1_outname="fig_targetL1"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    computeMatrix reference-point \
        -S ${top25_samples} \
        -R "${l1_bed}" \
        -b 1000 \
        --sortUsingSamples 1 \
        --numberOfProcessors 8 \
        --scale 1 \
        --binSize 10 \
        --averageTypeBins "median" \
        --outFileName "${outpath}/${l1_outname}_deeptools_matrix.gzip" \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --referencePoint "center" \
        </dev/null

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotProfile \
        --matrixFile "${outpath}/${l1_outname}_deeptools_matrix.gzip" \
        --outFileName "${outpath}/${l1_outname}_profile_only.pdf" \
        --samplesLabel ${top25_slabels} \
        --colors "#9B9B9B" "#715eee" "#ffb00e" \
        --perGroup \
        --plotType lines \
        --legendLocation upper-left \
        --outFileNameData "${outpath}/${l1_outname}_profile.tab" \
        </dev/null




#----------step 5

#===============================================
# top25 fig1b
#===============================================

t251b_outname="top25_fig1b"

t251b_samples="${old1_avr_path}/Anti-GFP_average.bw \
${old1_avr_path}/Anti-Mll1_average.bw \
${old1_avr_path}/Mll2_KO_Mll1_average.bw \
${old1_avr_path}/Anti-RbBp5_average.bw \
${old1_avr_path}/Double_KO_RbBP5_average.bw \
${avr_path}/D0_WT_H3K27ac_average.bw \
${avr_path}/D0_Mll1-KO_H3K27ac_average.bw \
${avr_path}/D0_Mll2-KO_H3K27ac_average.bw \
${avr_path}/D0_double-KO_H3K27ac_average.bw \
${avr_path}/D0_WT_H3K27me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K27me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K27me3_average.bw \
${avr_path}/D0_double-KO_H3K27me3_average.bw \
${avr_path}/D0_WT_H3K4me3_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me3_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me3_average.bw \
${avr_path}/D0_double-KO_H3K4me3_average.bw \
${avr_path}/D0_WT_H3K4me2_average.bw \
${avr_path}/D0_Mll1-KO_H3K4me2_average.bw \
${avr_path}/D0_Mll2-KO_H3K4me2_average.bw \
${avr_path}/D0_double-KO_H3K4me2_average.bw"

t251b_slabels="D0_WT_MLL2 \
D0_WT_MLL1 \
D0_Mll2-KO_MLL1 \
D0_WT_RBBP5 \
D0_double-KO_RBBP5 \
D0_WT_H3K27ac \
D0_Mll1-KO_H3K27ac \
D0_Mll2-KO_H3K27ac \
D0_double-KO_H3K27ac \
D0_WT_H3K27me3 \
D0_Mll1-KO_H3K27me3 \
D0_Mll2-KO_H3K27me3 \
D0_double-KO_H3K27me3 \
D0_WT_H3K4me3 \
D0_Mll1-KO_H3K4me3 \
D0_Mll2-KO_H3K4me3 \
D0_double-KO_H3K4me3 \
D0_WT_H3K4me2 \
D0_Mll1-KO_H3K4me2 \
D0_Mll2-KO_H3K4me2 \
D0_double-KO_H3K4me2"

t251b_thrsholds="10 10 10 15 15 70 70 70 70 80 80 80 80 160 160 160 160 150 150 150 150"

sudo docker run "${DOCKER_ARGS[@]}" -w "${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    computeMatrix reference-point \
        -S ${t251b_samples} \
        -R "${top25_bed}" \
        -b 1000 \
        --sortUsingSamples 1 \
        --numberOfProcessors 8 \
        --scale 1 \
        --binSize 10 \
        --averageTypeBins "median" \
        --outFileName "${outpath}/${t251b_outname}_deeptools_matrix.gzip" \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        --referencePoint "center" \
        </dev/null

# plotHeatmap ---

sudo docker run "${DOCKER_ARGS[@]}" -v "${outpath}:${outpath}" -u $(id -u):$(id -g) \
    -e MPLCONFIGDIR="${outpath}/.matplotlib" \
    quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0 \
    plotHeatmap \
        --colorMap 'seismic' \
        --missingDataColor 0.6 \
        --matrixFile "${outpath}/${t251b_outname}_deeptools_matrix.gzip" \
        --zMax ${t251b_thrsholds} \
        --outFileName "${outpath}/${t251b_outname}_plotHeatmap.pdf" \
        --sortUsingSamples 6 \
        --samplesLabel ${t251b_slabels} \
        --regionsLabel "top25pct_pcp" \
        --outFileNameMatrix "${outpath}/${t251b_outname}_plotHeatmap.mat.tab" \
        </dev/null
