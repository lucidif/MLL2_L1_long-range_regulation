# KMT2B-bound L1 to target gene APA analysis

This guide describes the analysis workflow used to test whether KMT2B-bound L1 elements are in physical proximity with KMT2B target genes, using ChIP-seq and Micro-C integration with the `HicAggR` package.

The repository includes a helper script in `bin/HicAggR_KMT2B_L1_target_gene_APA.R` that imports 5 kb `.mcool` matrices, indexes KMT2B-bound L1 anchors and target gene TSSs, and performs Aggregate Peak Analysis (APA) with insulation score constraints.

## 1. Overview

The analysis follows these major steps:

- define KMT2B-bound L1 anchors and KMT2B target gene pairs from GREAT-derived associations
- convert anchors and gene regions to `GRanges`
- import ICE-balanced Micro-C contact matrices at 5 kb resolution using `ImportHiC`
- build a genome constraint from insulation bigWig scores
- index anchors and baits to Hi-C bins with `IndexFeatures`
- extract oriented submatrices and aggregate them using `APA_analysis`
- visualize the mean contact signal using `ggAPA`

## 2. Required inputs

The script expects the following repository files and directories by default:

- `outs/Lara_multiomic_analysis/outs/coolpup/500bp/nowin_unsorted_anchors3.tsv`
  - anchor/gene pair table used to build the KMT2B-bound L1 to gene GRanges
- `in/mm10.sizes`
  - chromosome sizes for mm10
- `in/geo_micro_mcool/balanced_WT_day0.mcool`
- `in/geo_micro_mcool/balanced_KO_day0.mcool`
- `in/geo_micro_mcool/balanced_WT_day4.mcool`
- `in/geo_micro_mcool/balanced_KO_day4.mcool`
- `in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day0.Dd.cool_150kb.bigwig`
- `in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig`
- `in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig`
- `in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig`

## 3. Running the analysis

From the repository root, run the script with Rscript:

```bash
Rscript bin/HicAggR_KMT2B_L1_target_gene_APA.R
```

The script supports optional command-line arguments if your files are in different locations:

```bash
Rscript bin/HicAggR_KMT2B_L1_target_gene_APA.R \
  --pairs outs/Lara_multiomic_analysis/outs/coolpup/500bp/nowin_unsorted_anchors3.tsv \
  --chrom_sizes in/mm10.sizes \
  --mcool_dir in/geo_micro_mcool \
  --bigwig_dir in/2024_10_Lara_microC_downstream/func_insulation \
  --output_dir outs/Lara_multiomic_analysis/outs/HicAggR \
  --bin_size 5000 \
  --prefix KMT2B_L1_target_gene
```

## 4. Output

The script writes results to the output directory declared by `--output_dir`.

By default it produces:

- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_WT_day0.RData`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_KO_day0.RData`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_WT_day4.RData`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_KO_day4.RData`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_WT_day0_APA.pdf`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_KO_day0_APA.pdf`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_WT_day4_APA.pdf`
- `outs/Lara_multiomic_analysis/outs/HicAggR/KMT2B_L1_target_gene_KO_day4_APA.pdf`

Each PDF contains the aggregated heatmap produced by `ggAPA` for the specified condition.

## 5. Notes on the workflow

- The script uses the helper functions in `pipelines/nf-core-microc/bin/HicAggR_add_fun.R`.
- The `APA_analysis()` function derives TAD-like constraints from the provided insulation bigWig file before indexing features.
- `orientate = TRUE` ensures the directionality of anchors and genes is considered during APA extraction.
- If your anchor/gene pair file uses different column names, update the script or pass a custom formatted TSV file matching the expected columns.

## 6. Customization

If you want to test a different anchor-gene list, prepare a TSV file with these columns:

- `peak.chr`
- `peak.start`
- `peak.end`
- `peak`
- `peak.strand`
- `gene.chr`
- `gene.start`
- `gene.end`
- `gene.strand`
- `gene`

Then pass it with `--pairs path/to/custom_pairs.tsv`.

If your Micro-C matrices or bigwig constraints are stored elsewhere, use `--mcool_dir` and `--bigwig_dir`.

## 7. Troubleshooting

- If `Rscript` fails because a package is missing, install the required packages in your R environment: `HicAggR`, `GenomicRanges`, `InteractionSet`, `S4Vectors`, `GenomeInfoDb`, `rtracklayer`, `IRanges`, and `ggplot2`.
- If the helper file `HicAggR_add_fun.R` cannot be found, confirm that `pipelines/nf-core-microc/bin/HicAggR_add_fun.R` exists relative to the repository root.
- If the `.mcool` or `.bigwig` files are missing, verify the `--mcool_dir` and `--bigwig_dir` paths and the specified filenames.

## 8. Recommended follow-up

- Inspect the generated APA heatmaps for WT vs KO at day0 and day4.
- Compare the mean contact enrichment between conditions.
- Use the saved `.RData` objects to extract numerical APA matrices or to extend the analysis with custom plotting.
