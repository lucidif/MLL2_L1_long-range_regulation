# Cumulative frequency / distance-based analysis guide

This document describes the cumulative distribution analyses used in the MLL2/L1 long-range regulation project. These analyses compute nearest-feature distances between genomic feature pairs, plot cumulative distributions, and compare groups using Wilcoxon rank-sum tests.

## Overview

The paper-associated cumulative frequency analyses in this repository are implemented with R scripts under `bin/`. They calculate the distance from a primary feature to the nearest secondary feature and summarize the results as empirical cumulative distribution functions (ECDFs). The following comparisons are represented:

- **Genes vs. KMT2B peaks**
  - Distances from the TSS of differentially expressed genes (downregulated, upregulated, and unaffected) to the nearest distal CpG-positive KMT2B peak losing H3K4me3 in double-KO cells.
- **KMT2B peaks vs. transposable elements (TEs)**
  - Distances from four KMT2B peak categories (proximal ± CGI and distal ± CGI) losing H3K4me3 to nearest TE 5'UTRs, including L1 elements, SINEs, and LTRs (ERV and non-ERV subsets).
- **Genes vs. KMT2B-bound L1 elements**
  - Distances from the TSS of differentially expressed genes to the nearest KMT2B-bound L1 element, defined as an L1 5'UTR located within 1 kb of a distal CpG-negative KMT2B peak.

Statistical comparisons between groups are performed using Wilcoxon rank-sum tests.

## Main scripts

The repository includes these primary scripts for cumulative distribution analysis:

- `bin/Cumulative_freq_DEGs_to_dcp_peaks.R`
  - Genes vs. distal CpG-positive KMT2B peaks.
- `bin/cumulative_distribution_te_classes.R`
  - KMT2B peak categories vs. TE classes (LINE/L1, SINE, LTR/ERV).
- `bin/Cumulative_freq_DEGs_L1_MLL2peaks_.R`
  - Genes vs. KMT2B-bound L1 elements.
- `bin/Cumulative_freq_dcp_peaks_by_genes.R`
  - Additional gene vs. KMT2B peak distance distributions.

Each script is designed as a fixed-path analysis script and includes a `USER SETTINGS` section at the top where local input/output paths should be adjusted.

## Dependencies

The scripts rely on these helper sources:

- `pipelines/downstream_multiomic/bin/hypergeometric.R`
- `pipelines/downstream_multiomic/bin/distance_functions.R`
- `pipelines/bioinfo_generics/base/bin/fun_compare_curves.R`

The scripts are also written for the R packages:

- `data.table`
- `ggplot2`

A convenient runtime environment is the Rocker `tidyverse` Docker image.

## Running the analysis

1. Confirm you are in the repository root:

```bash
cd /workspaces/MLL2_L1_long-range_regulation
```

2. Edit the top of the selected script and update the `USER SETTINGS` paths to match your local data layout.

3. Pull the Docker image if needed:

```bash
sudo docker pull rocker/tidyverse:4.5.1
```

4. Run the script with Docker:

```bash
sudo docker run -it -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
  Rscript bin/Cumulative_freq_DEGs_to_dcp_peaks.R
```

For TE class analysis:

```bash
sudo docker run -it -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
  Rscript bin/cumulative_distribution_te_classes.R
```

For KMT2B-bound L1 analysis:

```bash
sudo docker run -it -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
  Rscript bin/Cumulative_freq_DEGs_L1_MLL2peaks_.R
```

For supplementary gene/peak distributions:

```bash
sudo docker run -it -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
  Rscript bin/Cumulative_freq_dcp_peaks_by_genes.R
```

## Output

The scripts write results into local `outs/` directories and generate PDF plots plus summary tables. Each script creates one or more subdirectories under `outs/` according to the analysis type.

## Notes

- The scripts use repository-specific paths by default. Always update the `USER SETTINGS` block before running.
- These analyses are intended as post-processing after ChIP-seq peak calling and RNA-seq differential expression analysis are available.
- The paper methods are reproduced in the repository scripts rather than in a notebook-specific workflow.
