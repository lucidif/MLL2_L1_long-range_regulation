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
### Containers used in this workflow (outside nf-core pipelines)

Below is the partial list of Docker images that are used outside of the nf-core standard workflows.

#### Core and supplementary workflow containers

| **Container** | **Source / Reference** | **Purpose** |
|----------------|------------------------|--------------|
| [`lucidif/fanc`](https://hub.docker.com/r/lucidif/fanc) | Docker Hub | Supplementary image for Micro-C and FAN-C analyses |
| [`lucidif/microc`](https://hub.docker.com/r/lucidif/microc) | Docker Hub | Custom container used across Micro-C workflow processes |
| [`rocker/tidyverse:4.5.1`](https://hub.docker.com/_/rocker) | Docker Hub | R environment for statistical analysis and visualization |
| [`lucidif/edger:0.0.1`](https://hub.docker.com/r/lucidif/edger) | Docker Hub | Custom image used for standalone edgeR-based RNA-seq analysis |

####  Standalone scripts (`bin/`, `pipelines/bioinfo_generics`)

| **Container** | **Source** | **Purpose** |
|----------------|-------------|--------------|
| [`rocker/tidyverse:4.5.1`](https://hub.docker.com/_/rocker) | Docker Hub | Used in R scripts for cumulative analysis and plotting |
| [`quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1`](https://quay.io/repository/biocontainers/bedtools) | Biocontainers | Invoked by `peaks_classification.sh` for genomic interval operations |
| [`quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0`](https://quay.io/repository/biocontainers/deeptools) | Biocontainers | Used in `fun_deeptools_heatmap_and_profile.sh` for matrix and profile generation |

> For the complete list of containers (including all nf-core process-specific images), see the [full container overview](docs/containers_full_list.md).


## Assets

The `assets/` directory contains supplementary data used by the workflow, including:

- **GREAT analysis results** (gene–regulatory associations).

> **Note:** Due to GitHub storage limits, all assets larger than **100 MB** are **not included** in this repository.  
> These files must be downloaded separately following the instructions provided in the project documentation or upon request from the authors.

## System requirements

The workflow is implemented using **Nextflow v23.04.4** and tested on **Linux (Ubuntu 22.04)** systems.  
It is compatible with **local** environments that support dockers. 

### Software dependencies

All dependencies are containerized. No local installations are required, but the following must be available on the host system:

- **Nextflow ≥ 23.04.4**  
- **Docker ≥ 24.0**  
- **Git ≥ 2.34**

### Tested environments

| System name | CPU / RAM | Storage | OS | Status |
|--------------|------------|----------|-----|---------|
| **device1** | AMD Ryzen 9 5900X (12 cores / 24 threads), 96 GiB RAM | SSD 1 TB + HDD 2 TB × 3 | Ubuntu 22.04 LTS | ✅ Tested |
| **device2** | AMD Ryzen 9 7900 (12 cores / 24 threads), 126 GiB RAM | SSD 2 TB + HDD 4 TB × 2 | Ubuntu 22.04 LTS | ✅ Tested |

### Hardware requirements

- **Minimum:** 8 CPU cores, 32 GB RAM, ~100 GB disk space  
  *(Minimum memory increased to meet STAR genome indexing requirements for mouse mm10)*  
- **Recommended:** ≥ 32 CPU cores, ≥ 64 GB RAM, ≥ 200 GB disk space for Micro-C workflows  
- **Tested:** 12 cores / 24 threads with 96–126 GB RAM on AMD Ryzen 9 systems

## License

This workflow is distributed under the **GNU General Public License v2.0 (GPL-2.0)**.  
You are free to use, modify, and distribute this software under the terms of the GPL-2.0 license.

See the full text at: [https://www.gnu.org/licenses/old-licenses/gpl-2.0.html](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
