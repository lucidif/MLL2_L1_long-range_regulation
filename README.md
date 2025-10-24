# MLL2_LINE1_long-range_regulation

## Introduction

**Nextflow Multi-Omic Analysis Workflow** is a modular and reproducible pipeline for the integrative analysis of **ChIP-seq**, **RNA-seq**, and **Micro-C** datasets.  
It reproduces and extends the computational procedures described in the *Methods* section of the associated study, enabling full reproducibility and transparent data processing.

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




## Docker and Pipelines

All processes are analyzed using Docker containers, or are already included within the nf-core pipelines (see the documentation of each individual pipeline).

For Docker images not included in nf-core pipelines, the corresponding containers can be downloaded from my Docker Hub:
https://hub.docker.com/u/lucidif
To use my custom Docker images inside a Nextflow pipeline, it is required to download the image from the indicated Docker Hub and rename it by prefixing quay.io/ to the image name.

Example:

```bash
docker pull lucidif/fanc
# Rename for nf-core compatibility:
# quay.io/lucidif/fanc
```

List of custom Docker images used outside standard nf-core community pipelines:

- lucidif/fanc
- lucidif/microc
- rocker/tidyverse:4.5.1
- 
