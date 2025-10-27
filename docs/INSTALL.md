# Installation Guide (Docker-based setup)

This guide explains how to install and run the **MLL2_LINE1_long-range_regulation** workflow using Docker containers.  
All dependencies are containerized to guarantee reproducibility and avoid local package conflicts.

---

## 1. Requirements

Before installation, ensure that the following software is installed and accessible in your environment:

- **Docker** ≥ 24.0  
- **Nextflow** ≥ 23.04.4  
- **Git** ≥ 2.34  

Optional: **Apptainer/Singularity** ≥ 1.2 (for HPC systems without Docker support)

---

## 2. Clone the repository

Download the repository and enter the project folder:

```bash
git clone https://github.com/lucidif/MLL2_L1_long-range_regulation.git
cd MLL2_L1_long-range_regulation
```

---

## 3. Install Pre-pull extra images (outside nf-core pipelines)

The nf-core-managed images are pulled automatically by Nextflow at runtime.  
However, this repository also relies on a set of **custom** and **third-party** images used outside the standard nf-core pipelines (standalone scripts, supplementary modules).  
Pre-pull them with the following command to avoid delays during the first run.

```bash
# Pull all extra images used outside nf-core pipelines
# You can paste this block as-is in your shell.

images=(
  "lucidif/fanc@sha256:f9737f9538068033854c473fcfe4ef416e8732e1dd488c3b07e8976dbf1599a2"
  "lucidif/microc:0.0.1"
  "lucidif/edger:0.0.1"
  "rocker/tidyverse:4.5.1"
  "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
  "quay.io/biocontainers/deeptools:3.5.5--py_0"
  "quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4"
  "quay.io/biocontainers/star:2.7.10a--h43eeafb_0"
  "quay.io/biocontainers/rsem:1.3.1--pl526haddd2b5_0"
  "quay.io/biocontainers/diffbind:3.14--r42hdfd78af_0"
  "quay.io/biocontainers/deseq2:1.34.0--r41hc247a5b_0"
)

for img in "${images[@]}"; do
  echo "Pulling $img"
  docker pull "$img"
done
```

**Notes:**

If your execution profile requires tagging custom images with a `quay.io` prefix, you can tag them after pulling, for example:

```bash
docker pull lucidif/fanc@sha256:f9737f9538068033854c473fcfe4ef416e8732e1dd488c3b07e8976dbf1599a2
docker tag $(docker images --no-trunc -q lucidif/fanc | head -n1) quay.io/lucidif/fanc:pin-f9737f9
```

---

## 4. Run the container

To start an interactive container session with access to your current directory:

```bash
docker run -it --rm   -v "$(pwd)":"$(pwd)"   lucidif/microc:0.0.1 /bin/bash
```

This mounts your current working directory into the container under the same path.  
You can now inspect tools or manually test pipeline commands inside this environment.  
Additional folders can be mounted using extra `-v` options (see Docker documentation for details).

---

## 5. Launch the Nextflow pipeline

Once the environment is ready, you can run the main Nextflow pipeline, for example:

```bash
nextflow run main.nf -profile docker   --input samplesheet.csv   --outdir results/
```

Replace `samplesheet.csv` with your own input sample sheet (see the examples in the `assets/` folder).  
All workflow dependencies are automatically pulled as Docker containers when the pipeline runs with the `-profile docker` flag.

---

## 6. Verify the installation

To confirm that all required tools are available within the container, run:

```bash
nextflow -version
docker --version
```

If both commands return valid version numbers, your setup is ready.

---

## 7. System and resource recommendations

- Tested on **Ubuntu 22.04 LTS**  
- Compatible with both local and HPC environments  

**Minimum (for STAR alignment on mouse mm10):**
- 8 CPU cores  
- 32 GB RAM  
- ~100 GB disk space  

**Recommended (for full Micro-C analysis):**
- ≥ 32 CPU cores  
- ≥ 64 GB RAM  
- ≥ 200 GB disk space  

---

## 8. Troubleshooting

If Docker cannot access your files due to permissions, try:

```bash
sudo usermod -aG docker $USER
newgrp docker
```

or run the container with elevated privileges (`sudo docker run ...`).

If Nextflow fails to download a container, check your network proxy or use:

```bash
nextflow pull lucidif/MLL2_L1_long-range_regulation
nextflow run lucidif/MLL2_L1_long-range_regulation -profile docker
```

---

## 9. References

- Docker Hub: https://hub.docker.com/u/lucidif  
- nf-core community pipelines: https://nf-co.re  
- Nextflow documentation: https://www.nextflow.io/docs/latest/  

---

## Author and License

**Author:** Lucio Di Francesco  
**License:** GNU General Public License v2.0 (GPL-2.0)  
**Repository:** https://github.com/lucidif/MLL2_L1_long-range_regulation