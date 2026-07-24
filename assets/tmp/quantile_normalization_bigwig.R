install_packages<-FALSE

#installation 
if(install_packages==TRUE) {

    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("epigenomix")

}



#load packages

library(rtracklayer)     # per import/export bedGraph/bigWig
library(GenomicRanges)
library(epigenomix)  # per normalize


#import files

bw_files <- list.files(
  path = "outs/geo_sub/out/MLL2_chipSEQ_processed",
  pattern = "\\.bigWig$",
  full.names = TRUE
)
repA_only<-bw_files[grep("_repA", bw_files)]
repB_only<-bw_files[grep("_repB", bw_files)]

markers<-c("H3K27ac","H3K27me3","H3K4me1","H3K4me2","H3K4me3","RBBP5","RING1B","MLL1", "MLL2")

tarmar<-1
tar_files<-repA_only[grep(paste0("_",markers[tarmar],"_"),repA_only)]

# Import dei 6 bigWig
sample_names <- sub("\\.bigWig$", "", basename(tar_files))
gr_list <- lapply(tar_files, import)
names(gr_list) <- sample_names

#check bin uniformity

ref <- gr_list[[1]]

for (nm in names(gr_list)) {
  gr <- gr_list[[nm]]
  message("Checking ", nm)
  stopifnot(
    identical(seqnames(ref), seqnames(gr)),
    identical(start(ref),    start(gr)),
    identical(end(ref),      end(gr))
  )
}


#extract bin coordinates
rowRanges <- ref
names(rowRanges) <- paste0(
  as.character(seqnames(rowRanges)), ":",
  start(rowRanges), "-",
  end(rowRanges)
)




# matrix of bins × samples 

chip_mat <- do.call(
  cbind,
  lapply(gr_list, function(gr) mcols(gr)$score)
)

colnames(chip_mat) <- sample_names
rownames(chip_mat) <- names(rowRanges)

dim(chip_mat)
