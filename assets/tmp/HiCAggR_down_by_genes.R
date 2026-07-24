suppressPackageStartupMessages({
  library(HicAggR)
  library(GenomicRanges)
  library(InteractionSet)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(IRanges)
  library(zoo)
  library(cowplot)
})

#==============support functions
load.anchors<-function(){
  anchor2 <- read.table("outs/coolpup/500bp/nowin_unsorted_d4_anchors2.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
anchor2$strand.peak[anchor2$strand.peak == "."] <- "*"

anc2_peaks.gr<-GRanges(
  seqnames = anchor2$chr.peak,
  ranges   = IRanges(start = as.integer(anchor2$start.peak),
                     end   = as.integer(anchor2$end.peak)),
  strand   = ifelse(anchor2$strand.peak %in% c("+","-"), anchor2$strand.peak, "*"),
  score= anchor2$score
)
mcols(anc2_peaks.gr)$name  <- anchor2$peak

anc2_genes.gr <- GRanges(
  seqnames = anchor2$chr.gene,
  ranges   = IRanges(start = as.integer(anchor2$start.gene),
                     end   = as.integer(anchor2$end.gene)),
  strand   = ifelse(anchor2$gene.strand %in% c("+","-"), anchor2$gene.strand, "*")
)
mcols(anc2_genes.gr)$name  <- anchor2$gene

return(
  list(
    anc2_peaks.gr,
    anc2_genes.gr
  )
)

}




#=============load helper functions, files and parameters

source("git/nf-core-microc/bin/HicAggR_add_fun.R")

chromSizes <- read.table("in/mm10.sizes", sep = "\t", header = FALSE,
                         col.names = c("seqnames","seqlengths"),
                         stringsAsFactors = FALSE)

binSize=5000

anchor2 <- read.table("outs/coolpup/500bp/nowin_unsorted_d4_anchors2.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
anchor2$strand.peak[anchor2$strand.peak == "."] <- "*"

anc2_peaks.gr<-GRanges(
  seqnames = anchor2$chr.peak,
  ranges   = IRanges(start = as.integer(anchor2$start.peak),
                     end   = as.integer(anchor2$end.peak)),
  strand   = ifelse(anchor2$strand.peak %in% c("+","-"), anchor2$strand.peak, "*"),
  score= anchor2$score
)
mcols(anc2_peaks.gr)$name  <- anchor2$peak

anc2_genes.gr <- GRanges(
  seqnames = anchor2$chr.gene,
  ranges   = IRanges(start = as.integer(anchor2$start.gene),
                     end   = as.integer(anchor2$end.gene)),
  strand   = ifelse(anchor2$gene.strand %in% c("+","-"), anchor2$gene.strand, "*")
)
mcols(anc2_genes.gr)$name  <- anchor2$gene






stopifnot(
  is.data.frame(chromSizes),
  all(c("seqnames","seqlengths") %in% names(chromSizes)),
  is.numeric(chromSizes$seqlengths),
  is.numeric(binSize),
  "score" %in% colnames(mcols(anc2_peaks.gr))
)

# anchor indexing (works; you verified it)

#import contact matrix
WT_day4 <- ImportHiC(
        file    = "in/geo_micro_mcool/balanced_WT_day4.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = TRUE
)



KO_day4 <- ImportHiC(
        file    = "in/geo_micro_mcool/balanced_KO_day4.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = TRUE
)


#=========
# D0
#=========


WT_day0 <- ImportHiC(
        file    = "in/geo_micro_mcool/balanced_WT_day0.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = FALSE
)


#import contact matrix
KO_day0 <- ImportHiC(
        file    = "in/geo_micro_mcool/balanced_KO_day0.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = FALSE
)


attributes(KO_day0)$mtx
summary(as.matrix(KO_day0$chr1_chr1@matrix))

#get_contact_matrix(KO_day0$chr1_chr1)


#====================
#   APA
#====================


APA_WT_day4<-APA_analysis(
  resolution = 5000,
  anchor_grange = anc2_peaks.gr,
  bait_grange = anc2_genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate=TRUE
)

WT_day4_p<-ggAPA(APA_WT_day4$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.8
)
WT_day4_p


APA_KO_day4<-APA_analysis(
  resolution = 5000,
  anchor_grange = anc2_peaks.gr,
  bait_grange = anc2_genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate=TRUE
)

KO_day4_p<-ggAPA(APA_KO_day4$APA, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.8
)

KO_day4_p

APA_WT_day0<-APA_analysis(
  resolution = 5000,
  anchor_grange = anc2_peaks.gr,
  bait_grange = anc2_genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day0.Dd.cool_150kb.bigwig",
  hic_in=WT_day0,
  orientate=TRUE
)

WT_day0_p<-ggAPA(APA_WT_day0$APA, title = "5kb_APA_WT_day0 APA"
  ,colMin = 0.1
  ,colMax = 0.8
)
WT_day0_p


APA_KO_day0<-APA_analysis(
  resolution = 5000,
  anchor_grange = anc2_peaks.gr,
  bait_grange = anc2_genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig",
  hic_in=KO_day0,
  orientate=TRUE
)

KO_day0_p<-ggAPA(APA_KO_day0$APA, title = "5kb_APA_KO_day0 APA"
  ,colMin = 0.1
  ,colMax = 0.8
)

KO_day0_p





