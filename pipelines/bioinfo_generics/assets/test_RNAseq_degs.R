sudo docker run -it -v /home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals:/home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals lucidif/edger:0.0.1

R

source("/home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals/RNAseq/bin/fun_statistical.R")

differential(prjFolder="/home/lucio/",
             prjName = "Lara_RNAseq_day0_4",
             rawCounts = "/home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals/assets/rawcounts.tsv.csv",
             designTablePath ="/home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals/assets/SS_RNAseq_Lara_edgeR_SS_D0.tsv",
             annoPath = "/home/lucio/MEGA/linux/ziggy/git/bioinfoGenerals/assets/gencode_GRCm38_p4_annotations.tsv",
             previous = FALSE,
             accurateFiltering = FALSE
             )