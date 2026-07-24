sudo docker pull rocker/r-ubuntu:24.04
sudo docker run rocker/r-ubuntu:24.04 R --version

# R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
# Copyright (C) 2025 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu

# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under the terms of the
# GNU General Public License versions 2 or 3.
# For more information about these matters see
# https://www.gnu.org/licenses/.



sudo docker run --rm -it \
  -v "$PWD":"$PWD" -w "$PWD" \
  rocker/r-ubuntu:24.04

#root@ed08bee73eed:/home/lucio/wkdir/projects/MLL2_L1_regulation
apt-get update

apt-get install -y \
  build-essential pkg-config \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype6-dev \
  libharfbuzz-dev libfribidi-dev \
  libcairo2-dev \
  libpng-dev libjpeg-dev libtiff5-dev libwebp-dev \
  libglpk-dev

R


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("preprocessCore")
install.packages("data.table")
BiocManager::install("edgeR")

#in un altro teminare:
sudo docker ps 

sudo docker commit ba38a89abba1 quantilenormalization:edger4.8.1

sudo docker tag  quantilenormalization:edger4.8.1 lucidif/quantilenormalization:edger4.8.1
sudo docker push lucidif/quantilenormalization:edger4.8.1
