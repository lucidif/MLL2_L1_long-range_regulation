# TE distribution analysis guide

This guide describes the repository workflow for TE distribution analysis using the helper script in `bin/TE_distribution_analysis.R`.

The script uses RepeatMasker annotations and nearest-TE distance calculations to summarize:

- TE macro-class composition near distal CpG-negative peaks
- L1 subfamily proportions near peaks compared with genome-wide L1 distributions

## 1. Required inputs

By default, the script expects the following repository files:

- `in/ucsc/rmsk.txt`
  - RepeatMasker annotation table from UCSC
- `outs/CHiP_postprocessing_line1_dist/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed`
  - distal CpG-negative peak BED file used for the analysis
- `pipelines/downstream_multiomic/bin/distance_functions.R`
  - helper script that provides the `extractdistances()` function used by the analysis

If your data are in different locations, update the variables at the top of `bin/TE_distribution_analysis.R` or modify the script to accept command-line arguments.

## 2. Running the analysis

From the repository root, run:

```bash
Rscript bin/TE_distribution_analysis.R
```

This will source the helper file, read the `rmsk.txt` and BED files, compute nearest TE distances, and write results to:

- `outs/TE_distribution/nearest_TE_within_window.tsv`
- `outs/TE_distribution/TE_macroclass_nearest_500bp.tsv`
- `outs/TE_distribution/TE_macroclass_nearest_500bp_donut.pdf`
- `outs/TE_distribution/L1_subfamily_genome_vs_peaks.tsv`
- `outs/TE_distribution/L1_subfamily_peaks_donut.pdf`
- `outs/TE_distribution/L1_subfamily_genome_donut.pdf`

## 3. Output details

The script produces two major result types:

1. `nearest_TE_within_window.tsv`
   - peak-to-TE nearest neighbor table for peaks within 500 bp of a TE element
   - includes RepeatMasker metadata fields, TE family, element, type, and derived macro-class

2. `TE_macroclass_nearest_500bp.tsv`
   - counts and percentages of the TE macro-classes detected within 500 bp of peaks
   - supported by a donut plot PDF for visualization

3. `L1_subfamily_genome_vs_peaks.tsv`
   - genome-wide L1 subfamily counts versus L1 subfamily counts observed near peaks
   - useful to compare whether peak-proximal L1 elements are enriched for specific subfamilies

## 4. Interpretation notes

- The script filters on the nearest TE element within a fixed `WINDOW_BP` of 500 bp.
- Macro-class calls are derived from RepeatMasker `element` and `type` fields.
- L1 subfamily comparison is based on the RepeatMasker family column.

## 5. Customization

To adapt the analysis for another peak set or window size:

- change `PEAK_BED` to the path of your new BED file
- change `WINDOW_BP` to a different distance threshold
- optionally update `RMSK_FILE` if you are using a different RepeatMasker annotation file

The script is intentionally implemented in base R style without `dplyr` or other tidyverse dependencies, matching the repository’s existing temporary scripts.

## 6. Troubleshooting

- If the script fails because a file is missing, confirm that the path in the top variables is correct and the named file exists.
- If the RepeatMasker file has fewer than 13 columns, verify that `in/ucsc/rmsk.txt` is the full UCSC annotation table and not a truncated subset.
- If no peaks are found within the 500 bp window, the script will stop with a message. In that case, either increase `WINDOW_BP` or verify that the input peak BED and RepeatMasker annotations refer to the same genome build.

## 7. Recommended follow-up

- Inspect the donut plots to see TE macro-class and L1 subfamily distributions.
- Use `nearest_TE_within_window.tsv` to explore the peak-level TE assignments.
- Compare the L1 peak subfamily proportions with genome-wide proportions to identify potential enrichment.
