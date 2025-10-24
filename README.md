# MLL2_LINE1_long-range_regulation

Docker and Pipelines

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
