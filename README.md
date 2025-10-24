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

## Assets

The `assets/` directory contains supplementary data used by the workflow, including:

- **GREAT analysis results** (gene–regulatory associations and enrichment outputs).

> **Note:** Due to GitHub storage limits, all assets larger than **100 MB** are **not included** in this repository.  
> These files must be downloaded separately following the instructions provided in the project documentation or upon request from the authors.
