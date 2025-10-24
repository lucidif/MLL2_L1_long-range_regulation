### Containers used in this workflow (outside nf-core pipelines)

Below is the complete list of Docker images used in this repository, grouped by usage area.  
These containers are referenced in the Nextflow processes, standalone scripts, or documentation files.

---

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


---

####  nf-core/microc pipeline

| **Container** | **Source** | **Purpose** |
|----------------|-------------|--------------|
| [`lucidif/microc:0.0.1`](https://hub.docker.com/r/lucidif/microc) | Docker Hub | Custom Micro-C workflow image (pairtools, cooler, FAN-C) |
| [`quay.io/biocontainers/samtools:1.12--hd5e65b6_0`](https://quay.io/repository/biocontainers/samtools) | Biocontainers | Used for chromosome size calculation (`CHROMSIZES`) |

---

####  nf-core/chipseq_RE pipeline

| **Container** | **Source** | **Purpose** |
|----------------|-------------|--------------|
| `mulled-v2-*` images (various) | Biocontainers | Used for MACS2 plotting, Homer annotation, DESeq2 QC, BAM filtering, Bowtie2 alignment, FRiP, etc. |
| [`quay.io/biocontainers/star:2.6.1d--0`](https://quay.io/repository/biocontainers/star) | Biocontainers | STAR aligner for ChIP-seq reads |
| [`quay.io/biocontainers/r-base:3.5.1`](https://quay.io/repository/biocontainers/r-base) | Biocontainers | Base R environment for PhantomPeakQualTools reports |
| [`quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0`](https://quay.io/repository/biocontainers/bedtools) | Biocontainers | Blacklist filtering and coverage computation |
| [`quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1`](https://quay.io/repository/biocontainers/multiqc) | Biocontainers | Aggregation of QC results |
| [`quay.io/biocontainers/perl:5.26.2`](https://quay.io/repository/biocontainers/perl) | Biocontainers | GTF → BED conversions |
| [`quay.io/biocontainers/samtools:1.15.1--h1170115_0`](https://quay.io/repository/biocontainers/samtools) | Biocontainers | Sorting, indexing, and general BAM operations |
| [`quay.io/biocontainers/samtools:1.18--h50ea8bc_1`](https://quay.io/repository/biocontainers/samtools) | Biocontainers | FASTA indexing (`SAMTOOLS_FAIDX`) |
| [`quay.io/biocontainers/subread:2.0.1--hed695b0_0`](https://quay.io/repository/biocontainers/subread) | Biocontainers | Feature quantification with `featureCounts` |
| [`quay.io/biocontainers/phantompeakqualtools:1.2.2--0`](https://quay.io/repository/biocontainers/phantompeakqualtools) | Biocontainers | SPP quality assessment |
| [`quay.io/biocontainers/homer:4.11--pl526hc9558a2_3`](https://quay.io/repository/biocontainers/homer) | Biocontainers | Peak annotation and motif discovery |
| [`quay.io/biocontainers/fastqc:0.11.9--0`](https://quay.io/repository/biocontainers/fastqc) | Biocontainers | Quality control of raw FASTQ files |
| [`quay.io/biocontainers/preseq:3.1.2--h445547b_2`](https://quay.io/repository/biocontainers/preseq) | Biocontainers | Library complexity estimation |
| [`quay.io/biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3`](https://quay.io/repository/biocontainers/macs2) | Biocontainers | MACS2 peak calling |
| [`quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2`](https://quay.io/repository/biocontainers/khmer) | Biocontainers | Unique k-mer estimation |
| [`quay.io/biocontainers/picard:2.27.4--hdfd78af_0`](https://quay.io/repository/biocontainers/picard) | Biocontainers | Picard utilities (e.g. MarkDuplicates) |
| [`quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0`](https://quay.io/repository/biocontainers/bowtie2) | Biocontainers | Index building and alignment |
| [`quay.io/biocontainers/python:3.8.3`](https://quay.io/repository/biocontainers/python) | Biocontainers | Samplesheet validation and helper scripts |

---

####  chipseq_downstream_macs pipeline

| **Container** | **Source** | **Purpose** |
|----------------|-------------|--------------|
| [`quay.io/biocontainers/bwa:0.7.17--hed695b0_7`](https://quay.io/repository/biocontainers/bwa) | Biocontainers | BWA indexing |
| [`quay.io/biocontainers/chromap:0.2.1--hd03093a_0`](https://quay.io/repository/biocontainers/chromap) | Biocontainers | Chromap index construction and alignment |
| [`quay.io/biocontainers/deeptools:3.5.1--py_0`](https://quay.io/repository/biocontainers/deeptools) | Biocontainers | Heatmap and profile plotting (`computeMatrix`, `plotHeatmap`, `plotProfile`) |
| [`quay.io/biocontainers/gffread:0.12.1--h8b12597_0`](https://quay.io/repository/biocontainers/gffread) | Biocontainers | GFF → GTF conversion |
| [`quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h446ed27_1`](https://quay.io/repository/biocontainers/ucsc-bedgraphtobigwig) | Biocontainers | BedGraph → BigWig conversion |
| [`quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0`](https://quay.io/repository/biocontainers/trim-galore) | Biocontainers | FASTQ trimming prior to alignment |

---

####  Documentation examples

| **Container** | **Source / Example** | **Purpose** |
|----------------|----------------------|--------------|
| [`quay.io/biocontainers/fastqc`](https://quay.io/repository/biocontainers/fastqc) | Mentioned in usage examples for container profiles | FASTQ quality check example |
| [`quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0`](https://quay.io/repository/biocontainers/pangolin) | Example of module container override | Demonstration of custom container substitution in Nextflow modules |
