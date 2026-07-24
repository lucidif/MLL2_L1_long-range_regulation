## =============================================================================
## HicAggR analysis workflow — cleaned up to IndexFeatures()
## TODO make docker 
## =============================================================================

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

pkgs <- c("HicAggR", "GenomicRanges", "InteractionSet", "S4Vectors",
          "GenomeInfoDb", "rtracklayer", "IRanges", "zoo", "cowplot")

sapply(pkgs, function(p) as.character(packageVersion(p)))

#=============load helper functions

source("git/nf-core-microc/bin/HicAggR_add_fun.R")


# =============================================================================
# 1) Load parametes and features (peaks/genes) and chromSizes 
#    - your features file is anchor+gene pairs; we create two GRanges
#    - ensure a numeric 'score' metadata column (required by IndexFeatures)
# =============================================================================

binSize=5000


anchor3 <- read.table("outs/Lara_multiomic_analysis/outs/coolpup/500bp/nowin_unsorted_anchors3.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
anchor3$peak.strand[anchor3$peak.strand == "."] <- "*"

l1_strand<-read.table("outs/Lara_multiomic_analysis/outs/coolpup/500bp/anchors1_line_associated_peaks.tsv",
header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

anchor3$..idx <- seq_len(nrow(anchor3))

anchor3_l1_strand<-merge(anchor3, l1_strand, by.x="peak", by.y="associated_peak", sort=TRUE, all.x=TRUE)

anchor3_l1_strand<-anchor3_l1_strand[order(anchor3_l1_strand$..idx), ]

head(anchor3_l1_strand)

peaks.gr <- GRanges(
  seqnames = anchor3$peak.chr,
  ranges   = IRanges(start = as.integer(anchor3$peak.start),
                     end   = as.integer(anchor3$peak.end)),
  strand   = ifelse(anchor3$peak.strand %in% c("+","-"), anchor3$peak.strand, "*")
)
mcols(peaks.gr)$name  <- anchor3$peak
mcols(peaks.gr)$score <- 1


genes.gr <- GRanges(
  seqnames = anchor3$gene.chr,
  ranges   = IRanges(start = as.integer(anchor3$gene.start),
                     end   = as.integer(anchor3$gene.end)),
  strand   = ifelse(anchor3$gene.strand %in% c("+","-"), anchor3$gene.strand, "*")
)
mcols(genes.gr)$name  <- anchor3$gene
mcols(genes.gr)$score <- 1

NROW(peaks.gr)
NROW(genes.gr)
# stopifnot(all(c("peak.chr","peak.start","peak.end","peak","line_strand") %in% names(anchor3_l1_strand)))

# anchor_l1.gr <- GenomicRanges::GRanges(
#   seqnames = anchor3_l1_strand$peak.chr,
#   ranges   = IRanges::IRanges(start = as.integer(anchor3_l1_strand$peak.start),
#                               end   = as.integer(anchor3_l1_strand$peak.end)),
#   strand   = ifelse(anchor3_l1_strand$line_strand %in% c("+","-"),
#                     anchor3_l1_strand$line_strand, "*")
# )
# S4Vectors::mcols(anchor_l1.gr)$name  <- anchor3_l1_strand$peak
# S4Vectors::mcols(anchor_l1.gr)$score <- 1



chromSizes <- read.table("outs/Lara_multiomic_analysis/in/mm10.sizes", sep = "\t", header = FALSE,
                         col.names = c("seqnames","seqlengths"),
                         stringsAsFactors = FALSE)


#broad_genes.gr <-rtracklayer::import("outs/Lara_multiomic_analysis/outs/broad_all_genes_with_strand_filtered_tabs.bed", format = "BED")
#mcols(broad_genes.gr)$score <- as.numeric(mcols(broad_genes.gr)$score)

# down regulated genes interactions

anchor2 <- read.table("outs/Lara_multiomic_analysis/outs/coolpup/500bp/nowin_unsorted_d4_anchors2.tsv",
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



split_strand<-make_strand_pairs(anchor3_l1_strand, chromSizes = chromSizes)

split_strand$anch_plus_bait_plus[[1]]
split_strand$anch_plus_bait_plus[[2]]

# =============================================================================
# 5) Final checks and IndexFeatures() (anchor side)
#    - binSize is the Hi-C resolution you’ll use downstream (e.g., 5000)
#    - metadataColName must exist and be numeric (we set 'score' = 1)
# =============================================================================



stopifnot(
  is.data.frame(chromSizes),
  all(c("seqnames","seqlengths") %in% names(chromSizes)),
  is.numeric(chromSizes$seqlengths),
  is.numeric(binSize),
  "score" %in% colnames(mcols(peaks.gr))
)

# anchor indexing (works; you verified it)

#import contact matrix

attributes(KO_day0)$mtx
summary(as.matrix(KO_day0$chr1_chr1@matrix))

#get_contact_matrix(KO_day0$chr1_chr1)

WT_day0 <- ImportHiC(
        file    = "outs/Lara_multiomic_analysis/in/geo_micro_mcool/balanced_WT_day0.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = FALSE
)

attributes(WT_day0)$mtx

KO_day0 <- ImportHiC(
        file    = "outs/Lara_multiomic_analysis/in/geo_micro_mcool/balanced_KO_day0.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = FALSE
)


WT_day4 <- ImportHiC(
        file    = "outs/Lara_multiomic_analysis/in/geo_micro_mcool/balanced_WT_day4.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = TRUE
)

attributes(WT_day4)$mtx
try(S4Vectors::metadata(WT_day4[[1]])$mtx, silent = TRUE)

KO_day4 <- ImportHiC(
        file    = "outs/Lara_multiomic_analysis/in/geo_micro_mcool/balanced_KO_day4.mcool",
        hicResolution     = 5000,
        chrom_1 = chromSizes$seqnames,
        chrom_2 = chromSizes$seqnames
        #,cool_balanced = TRUE
        #,cool_weight_name = "weight"
        #,cool_divisive_weights = TRUE
)



#==================================================
# all regions
#==================================================

APA_KO_day4<-APA_analysis(
  resolution = 5000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = TRUE
)

KO_day4_p<-ggAPA(APA_KO_day4$APA, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0
  ,colMax = 0.8
)

KO_day4_p

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_KO_day4_APA.pdf", plot = KO_day4_p, width = 7, height = 7, units = "in")


#===========================
# WT DAY4 APA
#===========================


APA_WT_day4<-APA_analysis(
  resolution = 5000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate=TRUE
)

save(APA_WT_day4, file = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/APA_WT_day4.RData")
attributes(APA_WT_day4$H_used)$mtx

#APA_WT_day4$M_or
#APA_WT_day4$H_used
#attr(APA_WT_day4$APA, "referencePoint") <- "pf"

WT_day4_p<-ggAPA(APA_WT_day4$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.8
)
WT_day4_p

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_WT_day4_APA.pdf", plot = WT_day4_p, width = 7, height = 7, units = "in")


#===========================================
# scale optimizzed for low covarage regions D4
#===========================================

WT_day4_p_rescale<-ggAPA(APA_WT_day4$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0.15
  ,colMax = 0.6
)

#WT_day4_p_rescale

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/rescale_5kb_WT_day4_APA.pdf", plot = WT_day4_p_rescale, width = 7, height = 7, units = "in")

KO_day4_p_rescale<-ggAPA(APA_KO_day4$APA, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0.15
  ,colMax = 0.6
)

#KO_day4_p_rescale

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/rescale_5kb_KO_day4_APA.pdf", plot = KO_day4_p_rescale, width = 7, height = 7, units = "in")




#=========================================================
#                       run  D0 apa analysis
#=========================================================

APA_WT_day0<-APA_analysis(
  resolution = 5000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day0.Dd.cool_150kb.bigwig",
  hic_in=WT_day0
)

save(APA_WT_day0, file = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/APA_WT_day0.RData")

WT_day0_p<-ggAPA(APA_WT_day0$APA, title = "5kb_APA_WT_day0 APA"
,colMin = 0
  ,colMax = 0.8
)

WT_day0_p

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_WT_day0_APA.pdf", plot = WT_day0_p, width = 7, height = 7, units = "in")


APA_KO_day0<-APA_analysis(
  resolution = 5000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig",
  hic_in=KO_day0
)

save(APA_KO_day0, file = "outs/Lara_multiomic_analysis/in/2024_10_Lara_microC_downstream/func_insulation/APA_KO_day0.RData")

KO_day0_p<-ggAPA(APA_KO_day0$APA, title = "APA_KO_day0 APA"
,colMin = 0
  ,colMax = 0.8
)

KO_day0_p
ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_KO_day0_APA.pdf", plot = KO_day0_p, width = 7, height = 7, units = "in")



#===========================================
# scale optimizzed for low covarage regions D0
#===========================================

WT_day0_p_rescale<-ggAPA(APA_WT_day0$APA, title = "5kb_APA_WT_day0 APA"
  ,colMin = 0.15
  ,colMax = 0.6
)

#WT_day0_p_rescale

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/rescale_5kb_WT_day0_APA.pdf", plot = WT_day0_p_rescale, width = 7, height = 7, units = "in")

KO_day0_p_rescale<-ggAPA(APA_KO_day0$APA, title = "5kb_APA_KO_day0 APA"
  ,colMin = 0.15
  ,colMax = 0.6
)

#KO_day0_p_rescale

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/rescale_5kb_KO_day0_APA.pdf", plot = KO_day0_p_rescale, width = 7, height = 7, units = "in")















#===========================
# TODEL 
#===========================



#======
# D4
#======

APA_WT_day4<-APA_analysis(
  resolution = 5000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate=TRUE
)

#attr(APA_WT_day4$APA, "referencePoint") <- "pf"

WT_day4_p<-ggAPA(APA_WT_day4$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.8
)
WT_day4_p

ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_WT_day4_APA.pdf", plot = WT_day4_p, width = 7, height = 7, units = "in")

APA_WT_day4_pp<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_plus_bait_plus[[1]] ,
  bait_grange = split_strand$anch_plus_bait_plus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4
)

WT_day4_pp_p<-ggAPA(APA_WT_day4_pp$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.8
)

WT_day4_pp_p


ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_WT_day4_APA_pp.pdf", plot = WT_day4_pp_p, width = 7, height = 7, units = "in")

APA_WT_day4_mm<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_minus_bait_minus[[1]] ,
  bait_grange = split_strand$anch_minus_bait_minus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4
  )
WT_day4_mm_p<-ggAPA(APA_WT_day4_mm, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.8
)
WT_day4_mm_p
ggplot2::ggsave("outs/Lara_multiomic_analysis/outs/HicAggR/temp/5kb_WT_day4_APA_mm.pdf", plot = WT_day4_mm_p, width = 7, height = 7, units = "in")



APA_KO_day4_pp<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_plus_bait_plus[[1]],
  bait_grange = split_strand$anch_plus_bait_plus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4
)

KO_day4_pp_p<-ggAPA(APA_KO_day4_pp, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0
  ,colMax = 0.6
)
KO_day4_pp_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_KO_day4_APA_pp.pdf", plot = KO_day4_pp_p, width = 7, height = 7, units = "in")


APA_KO_day4_mm<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_minus_bait_minus[[1]],
  bait_grange = split_strand$anch_minus_bait_minus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4
)

KO_day4_mm_p<-ggAPA(APA_KO_day4_mm, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0
  ,colMax = 0.6
)
KO_day4_mm_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_KO_day4_APA.pdf", plot = KO_day4_p, width = 7, height = 7, units = "in")



#15kb 


APA_WT_day4_15kb<-APA_analysis(
  resolution = 15000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4_15kb
)

WT_day4_15kb_p<-ggAPA(APA_WT_day4_15kb$APA, title = "15kb_APA_WT_day4 APA"
  #,colMin = 0
  #,colMax = 20000
)
WT_day4_15kb_p

APA_KO_day4_15kb<-APA_analysis(
  resolution = 15000,
  anchor_grange = peaks.gr,
  bait_grange = genes.gr,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4_15kb
)

KO_day4_15kb_p<-ggAPA(APA_KO_day4_15kb, title = "15kb_APA_KO_day4 APA"
  #,colMin = 0
  #,colMax = 20000
)
KO_day4_15kb_p



