# KMT2B_LINE1_long-range_regulation

## Introduction

**Nextflow Multi-Omic Analysis Workflow** is a modular and reproducible workflow for the integrative analysis of **ChIP-seq**, **RNA-seq**, **Micro-C**, and **scRNA-seq** datasets.
It reproduces and extends the computational procedures described in the *Methods* section of the associated study (https://doi.org/10.1101/2025.08.10.669526).

The workflow combines **nf-core community pipelines** with **custom Nextflow modules** and **R/Bash scripts** for multi-omic integration and visualization.

For installation and usage instructions, please refer to the [Installation Guide](docs/INSTALL.md).


## Pipeline summary

| **Analysis type** | **Workflow / Tools** | **Description** |
|--------------------|----------------------|-----------------|
| **ChIP-seq** | [`nf-core/chipseq v2.0.0`](https://nf-co.re/chipseq) | Read trimming (`Trim Galore` v0.6.7 / `cutadapt` v3.4), alignment (`STAR` v2.6.1d, "mouse random mode" for TE-compatible mapping), sorting/dedup (`samtools` v1.15.1, `Picard` v2.27.4), spike-in-free normalization (`ChIPseqSpikeInFree` v1.2.4), peak calling (`MACS2` v2.2.7.1, `--broad`), differential binding (`DiffBind` v3.14), coverage/visualization (`deepTools` v3.5.5). |
| **RNA-seq (quantification & DEG)** | [`nf-core/rnaseq v3.14.0`](https://nf-co.re/rnaseq) (`--aligner star_rsem`) + [`nf-core/differentialabundance`](https://nf-co.re/differentialabundance) | Alignment (`STAR` 2.7.10a), quantification (`RSEM` 1.3.1), differential expression (`DESeq2` v1.34.0), functional enrichment (`WebGestalt` GO "biological process", `Enrichr` "ChEA 2022"). |
| **RNA-seq (de novo assembly)** | `nf-core/rnaseq v3.14.0` (`--aligner star_salmon`) + `StringTie` | De novo transcript assembly used to test whether KMT2B-bound L1 elements act as alternative promoters, by checking for continuous transcriptional signal between L1 5'UTR and neighboring gene exons (visual inspection in genome browser). |
| **TE/L1 subfamily analysis** | `RepeatMasker` (UCSC mm10 annotations) + `TEtranscripts` | Classification of TE macro-classes/L1 subfamilies near KMT2B peaks; differential expression of L1 subfamilies between genotypes (DESeq2 within `TEtranscripts`). |
| **Micro-C** | Custom Nextflow pipeline (Dovetail Genomics protocol) | Alignment (`BWA-MEM`), parsing/deduplication (`pairtools`, `--min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30`), contact matrix generation + ICE balancing (`cooler`), A/B compartment analysis (`FAN-C`), saddle plots (`cooltools saddle`), TAD calling (`SpectralTAD`). |
| **Virtual 4C** | Custom Nextflow pipeline | Locus-specific reference extraction (TSS–L1 region, ±50 Kb padding), alignment (`BWA-MEM`), viewpoint-based read selection (`pairtools select`, viewpoint window ±1–1.5 Kb), custom liftover of filtered reads to mm10 coordinates, viewpoint-local normalization (scale factor = 1e6 / Σ total_nodups), track generation (`bamCoverage`, deepTools). |
| **Micro-C / ChIP-seq integration (APA)** | `HicAggR` v1.0.2 | Aggregate Peak Analysis (`HicAggR`) to quantify physical proximity between KMT2B-bound L1 elements and the promoters of putative target genes. |


## Integrative downstream analyses

The workflow includes dedicated modules for integrative data analysis:

- **Gene classification** based on H3K27me3 domain width at the TSS (broad ≥6 Kb, narrow <6 Kb, negative = 0).
- **Cumulative frequency / distance-based analyses** (Wilcoxon rank-sum test) between KMT2B ChIP-seq peaks, TE/L1 5'UTRs, and differentially expressed genes.
- **TE/L1 composition analysis** (donut plots) of the TE macro-classes/subfamilies nearest to distal CGI-negative KMT2B peaks.
- **Aggregate Peak Analysis** (`HicAggR`) integrating Micro-C contact maps and ChIP-seq peaks to test physical proximity between KMT2B-bound L1 elements and target gene promoters.
- **Locus-specific Virtual 4C profiles** for individual validation of L1–gene physical interactions.
- **Functional enrichment** of DEG sets (`WebGestalt` GO; `Enrichr` ChEA 2022 TF targets).

## Software implementation and containers

All processes are executed within **Docker containers** to ensure full reproducibility and portability across environments.

- **Standard containers** are automatically pulled from the [nf-core community registry](https://nf-co.re).
- **Custom containers** are hosted under [lucidif Docker Hub](https://hub.docker.com/u/lucidif).
- **Third-party containers** (e.g., `rocker/tidyverse`, Biocontainers) are used for specific R-based post-processing, statistics, and visualization steps.

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
| [`lucidif/fanc`](https://hub.docker.com/r/lucidif/fanc) | Docker Hub | A/B compartment analysis (FAN-C) on Micro-C data |
| [`lucidif/microc`](https://hub.docker.com/r/lucidif/microc) | Docker Hub | Custom container used across Micro-C workflow processes (BWA-MEM, pairtools, cooler) |
| [`lucidif/chipseqspikeinfree:1.2.4`](https://hub.docker.com/r/lucidif/chipseqspikeinfree) | Docker Hub | Spike-in-free ChIP-seq normalization (`ChIPseqSpikeInFree` R package) |
| [`rocker/tidyverse:4.5.1`](https://hub.docker.com/_/rocker) | Docker Hub | R environment for statistical analysis and visualization |
| [`quay.io/biocontainers/cooltools:0.7.1--py311h93dcfea_3`](https://quay.io/repository/biocontainers/cooltools) | Biocontainers | Saddle plot generation (`cooltools saddle`) for compartment strength analysis |


####  Standalone scripts (`bin/`, `pipelines/bioinfo_generics`)

| **Container** | **Source** | **Purpose** |
|----------------|-------------|--------------|
| [`rocker/tidyverse:4.5.1`](https://hub.docker.com/_/rocker) | Docker Hub | Used in R scripts for cumulative distance analysis, gene/H3K27me3-domain classification, and plotting |
| [`quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1`](https://quay.io/repository/biocontainers/bedtools) | Biocontainers | Invoked by `peaks_classification.sh` for genomic interval operations (proximal/distal, CGI +/- classification) |
| [`quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0`](https://quay.io/repository/biocontainers/deeptools) | Biocontainers | Used in `fun_deeptools_heatmap_and_profile.sh` for matrix/heatmap/profile generation and BAM→bigwig coverage (incl. Virtual 4C tracks) |

> For the complete list of containers (including all nf-core process-specific images), see the [full container overview](docs/containers_full_list.md).


## Assets

The `assets/` directory contains supplementary data used by the workflow, including:

- **GREAT analysis results** (gene–regulatory associations).
- **Locus-specific reference coordinates** for Virtual 4C analysis (`VirtualC_ref_coordinates.tsv`).

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
