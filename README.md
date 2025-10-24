# MLL2_LINE1_long-range_regulation

## Introduction

**Nextflow Multi-Omic Analysis Workflow** is a modular and reproducible workflow for the integrative analysis of **ChIP-seq**, **RNA-seq**, and **Micro-C** datasets.  
It reproduces and extends the computational procedures described in the *Methods* section of the associated study (https://doi.org/10.1101/2025.08.10.669526), enabling full reproducibility and transparent data processing.

The workflow combines **nf-core community pipelines** with **custom Nextflow modules** and **R/Bash scripts** for multi-omic integration and visualization.

## Pipeline summary

| **Analysis type** | **Workflow / Tools** | **Description** |
|--------------------|----------------------|-----------------|
| **ChIP-seq** | [`nf-core/chipseq v2.0.0`](https://nf-co.re/chipseq) | Read trimming (`Trim Galore`), alignment (`STAR`), duplicate marking (`Picard`), peak calling (`MACS2`), differential analysis (`DiffBind`), and visualization (`deepTools`). |
| **RNA-seq** | [`nf-core/rnaseq v3.14.0`](https://nf-co.re/rnaseq) + [`nf-core/differentialabundance`](https://nf-co.re/differentialabundance) | Alignment (`STAR`), quantification (`RSEM`), differential expression (`DESeq2`). |
| **Micro-C** | Custom Nextflow pipeline | Implements the Dovetail Genomics protocol: alignment (`BWA-MEM`), parsing/deduplication (`pairtools`), contact map generation (`cooler`), and higher-order chromatin analysis (`FAN-C`, `TADCompare`). |

## Integrative downstream analyses

The workflow includes dedicated modules for integrative data analysis:

- **Gene classification** based on H3K27me3 domain width.
- **Cumulative frequency** and distance-based analyses between genomic features.
- **Integration of Micro-C contact maps** and ChIP-seq peaks to identify physical proximity between **MLL2-bound LINE-1 elements** and their putative target genes.

## Software implementation and containers

All processes are executed within **Docker containers** to ensure full reproducibility and portability across environments.

- **Standard containers** are automatically pulled from the [nf-core community registry](https://nf-co.re).  
- **Custom containers** are hosted under [lucidif Docker Hub](https://hub.docker.com/u/lucidif).  
- **Third-party containers** (e.g., `rocker/tidyverse`) are used for specific R-based post-processing and visualization steps.

To use custom or third-party containers within nf-core-compatible workflows, download them and rename with a `quay.io/` prefix if required.

**Example:**
```bash
docker pull lucidif/fanc
# Rename for nf-core compatibility:
# quay.io/lucidif/fanc
```
Custom Docker images used:

- lucidif/fanc
- lucidif/microc
- rocker/tidyverse:4.5.1


## Assets

The `assets/` directory contains supplementary data used by the workflow, including:

- **GREAT analysis results** (gene–regulatory associations).

> **Note:** Due to GitHub storage limits, all assets larger than **100 MB** are **not included** in this repository.  
> These files must be downloaded separately following the instructions provided in the project documentation or upon request from the authors.

### Containers used in this workflow (outside nf-core pipelines)

| **Container** | **Source** | **Purpose** |
|----------------|------------|--------------|
| [`lucidif/fanc`](https://hub.docker.com/r/lucidif/fanc) | Docker Hub | FAN-C environment for Micro-C and chromatin contact map analyses |
| [`lucidif/microc`](https://hub.docker.com/r/lucidif/microc) | Docker Hub | Custom Micro-C processing workflow (pairtools, cooler, FAN-C) |
| [`rocker/tidyverse:4.5.1`](https://hub.docker.com/_/rocker) | Docker Hub | R-based statistical analysis and visualization (ggplot2, dplyr, etc.) |
| [`quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1`](https://quay.io/repository/biocontainers/bedtools) | Biocontainers | Genomic interval operations (intersection, coverage, shuffling) |
| [`quay.io/biocontainers/deeptools:3.5.5--py_0`](https://quay.io/repository/biocontainers/deeptools) | Biocontainers | Signal track generation and heatmap plotting |
| [`quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4`](https://quay.io/repository/biocontainers/macs2) | Biocontainers | Peak calling for ChIP-seq analysis |
| [`quay.io/biocontainers/star:2.7.10a--h43eeafb_0`](https://quay.io/repository/biocontainers/star) | Biocontainers | RNA-seq and ChIP-seq read alignment |
| [`quay.io/biocontainers/rsem:1.3.1--pl526haddd2b5_0`](https://quay.io/repository/biocontainers/rsem) | Biocontainers | Transcript quantification for RNA-seq |
| [`quay.io/biocontainers/diffbind:3.14--r42hdfd78af_0`](https://quay.io/repository/biocontainers/diffbind) | Biocontainers | Differential binding analysis for ChIP-seq |
| [`quay.io/biocontainers/deseq2:1.34.0--r41hc247a5b_0`](https://quay.io/repository/biocontainers/deseq2) | Biocontainers | Differential expression analysis for RNA-seq |

## License

This workflow is distributed under the **GNU General Public License v2.0 (GPL-2.0)**.  
You are free to use, modify, and distribute this software under the terms of the GPL-2.0 license.

See the full text at: [https://www.gnu.org/licenses/old-licenses/gpl-2.0.html](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
