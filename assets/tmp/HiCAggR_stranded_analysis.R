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

#=============load helper functions, files and parameters

source("git/nf-core-microc/bin/HicAggR_add_fun.R")


binSize=5000

degs_d4<-read.table("in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv",
                  #sep="\t",
                  header=TRUE
)

degs_d4<-degs_d4[which(degs_d4$padj<=0.05),]

anchor3 <- read.table("outs/coolpup/500bp/nowin_unsorted_anchors3.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
anchor3$peak.strand[anchor3$peak.strand == "."] <- "*"

l1_strand<-read.table("outs/coolpup/500bp/anchors1_line_associated_peaks.tsv",
header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

chromSizes <- read.table("in/mm10.sizes", sep = "\t", header = FALSE,
                         col.names = c("seqnames","seqlengths"),
                         stringsAsFactors = FALSE)

broad_genes.gr <-rtracklayer::import("outs/broad_all_genes_with_strand_filtered_tabs.bed", format = "BED")
mcols(broad_genes.gr)$score <- as.numeric(mcols(broad_genes.gr)$score)

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

#=======================================



anchor3$..idx <- seq_len(nrow(anchor3))

anchor3_l1_strand<-merge(anchor3, l1_strand, by.x="peak", by.y="associated_peak", sort=TRUE, all.x=TRUE)

anchor3_l1_std_plus<-anchor3_l1_strand[which(anchor3_l1_strand$line_strand=="+"),]

paste0("anchor3_l1_strand rows:",nrow(anchor3_l1_strand))
paste0("anchor3_l1_std_plus rows:",nrow(anchor3_l1_std_plus))

#=======ALL combination out
anchor3_l1_strand_all_comb<-make_strand_pairs(anchor3_l1_strand,chromSizes = chromSizes, use_gene_strand=TRUE)

anchor3_l1_strand_all_comb$anch_plus_bait_plus[[1]]
anchor3_l1_strand_all_comb$anch_plus_bait_plus[[2]]
anchor3_l1_strand_all_comb$anch_minus_bait_minus[[1]]
anchor3_l1_strand_all_comb$anch_plus_bait_minus[[1]]
anchor3_l1_strand_all_comb$anch_minus_bait_plus[[1]]


#take all genes WITH L1 +
anchor3_l1_std_plus_split<-make_strand_pairs(anchor3_l1_std_plus, chromSizes = chromSizes, use_gene_strand=FALSE)
NROW(anchor3_l1_std_plus_split$anch_plus_bait_plus[[1]])

#take only the genes that came after the lines
anchor3_l1_std_plus<-anchor3_l1_std_plus[which(anchor3_l1_std_plus$peak.end<=anchor3_l1_std_plus$gene.start),]
paste0("anchor3_l1_std_plus (upstream genes removed) rows:",nrow(anchor3_l1_std_plus))

regions_distances<-(anchor3_l1_std_plus$gene.start+(anchor3_l1_std_plus$gene.end - anchor3_l1_std_plus$gene.start))/2 - (anchor3_l1_std_plus$peak.start+(anchor3_l1_std_plus$peak.end-anchor3_l1_std_plus$peak.start))/2 

plot(density(regions_distances))
summary(regions_distances)

colnames(anchor3_l1_strand)

#anchor3_l1_strand<-anchor3_l1_strand[order(anchor3_l1_strand$..idx), ]

split_strand<-make_strand_pairs(anchor3_l1_std_plus, chromSizes = chromSizes, use_gene_strand=FALSE)

split_strand$anch_plus_bait_plus[[1]]
NROW(split_strand$anch_plus_bait_plus[[1]])
split_strand$anch_plus_bait_plus[[2]]

split_strand$anch_minus_bait_plus[[1]]
split_strand$anch_minus_bait_plus[[2]]


anchor3_l1_std_minus<-anchor3_l1_strand[which(anchor3_l1_strand$line_strand=="-"),]

#take only the genes that came after the lines
anchor3_l1_std_minus<-anchor3_l1_std_minus[which(anchor3_l1_std_minus$peak.end>=anchor3_l1_std_minus$gene.start),]

split_strand_minus<-make_strand_pairs(anchor3_l1_std_minus, chromSizes = chromSizes, use_gene_strand=FALSE)

split_strand_minus$anch_plus_bait_plus[[1]]
split_strand_minus$anch_plus_bait_plus[[2]]

split_strand_minus$anch_minus_bait_minus[[1]]
split_strand_minus$anch_minus_bait_minus[[2]]
NROW(split_strand_minus$anch_minus_bait_minus[[2]])

gr <- split_strand_minus$anch_minus_bait_minus[[2]]
table(strand(gr))


#TODO test overlap sgnificance
#intgenes_minus<-intersect(degs_d4[,1],genes)
#intgenes_minus<-merge(intgenes_minus,degs_d4, by=1)

gr_plus<- split_strand$anch_plus_bait_plus[[2]]
table(strand(gr_plus))
plus_genes <- unique(S4Vectors::mcols(gr_plus)$name)

#TODO test overlap sgnificance
intgenes_plus<-intersect(degs_d4[,1],plus_genes)
intgenes_plus<-merge(intgenes_plus,degs_d4, by=1)

nrow(intgenes_plus)

#length(which(intgenes_plus$log2FoldChange<0))/nrow(intgenes_plus)
#length(which(intgenes_minus$log2FoldChange<0))/nrow(intgenes_minus)



#=========================================================
#                       run all apa analysis
#=========================================================

#======
# D4
#======

#===========================================
# plus all combination
#===========================================

NROW(anchor3_l1_std_plus_split$anch_plus_bait_plus[[1]])
NROW(anchor3_l1_std_plus_split$anch_plus_bait_plus[[2]])

#FILTER BY gene position

pp_int<-cbind(as.data.frame(anchor3_l1_strand_all_comb$anch_plus_bait_plus[[1]]), 
as.data.frame(anchor3_l1_strand_all_comb$anch_plus_bait_plus[[2]]))

filtered_ids<-which(pp_int[,9]+(pp_int[,10]-pp_int[,9])/2>= (pp_int[,2]+ round((pp_int[,3]-pp_int[,2])/2)) )

pp_anchor<-anchor3_l1_strand_all_comb$anch_plus_bait_plus[[1]][filtered_ids,]
pp_bait<-anchor3_l1_strand_all_comb$anch_plus_bait_plus[[2]][filtered_ids,]

#++
APA_WT_day4_p<-APA_analysis(
  resolution = 5000,
  anchor_grange = pp_anchor ,
  bait_grange = pp_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)
#297 matrices are oriented
WT_day4_p_p<-ggAPA(APA_WT_day4_p$APA, title = "++ 5kb_APA_WT_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

WT_day4_p_p

#anchor3_l1_strand_all_comb$anch_minus_bait_minus

APA_KO_day4_p<-APA_analysis(
  resolution = 5000,
  anchor_grange = pp_anchor ,
  bait_grange = pp_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)
#297 matrices are oriented
KO_day4_p_p<-ggAPA(APA_KO_day4_p$APA, title = "++ 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

KO_day4_p_p





# - - 

mm_int <- cbind(as.data.frame(anchor3_l1_strand_all_comb$anch_minus_bait_minus[[1]]), 
as.data.frame(anchor3_l1_strand_all_comb$anch_minus_bait_minus[[2]]))
filtered_ids<-which(mm_int[,9]+(mm_int[,10]-mm_int[,9])/2<= (mm_int[,2]+ round((mm_int[,3]-mm_int[,2])/2)) )

mm_anchor <- anchor3_l1_strand_all_comb$anch_minus_bait_minus[[1]][filtered_ids,]
mm_bait <- anchor3_l1_strand_all_comb$anch_minus_bait_minus[[2]][filtered_ids,]

APA_WT_day4_m<-APA_analysis(
  resolution = 5000,
  anchor_grange = mm_anchor ,
  bait_grange = mm_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)
#297 matrices are oriented
WT_day4_m_m<-ggAPA(APA_WT_day4_m$APA, title = "-- 5kb_APA_WT_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

WT_day4_m_m

APA_KO_day4_m<-APA_analysis(
  resolution = 5000,
  anchor_grange = mm_anchor ,
  bait_grange = mm_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)
#297 matrices are oriented
KO_day4_m_m<-ggAPA(APA_KO_day4_m$APA, title = "-- 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

KO_day4_m_m




# +-

pm_int <- cbind(as.data.frame(anchor3_l1_strand_all_comb$anch_plus_bait_minus[[1]]), 
as.data.frame(anchor3_l1_strand_all_comb$anch_plus_bait_minus[[2]]))


filtered_ids<-which(pm_int[,9]+(pm_int[,10]-pm_int[,9])/2<= (pm_int[,2]+ round((pm_int[,3]-pm_int[,2])/2)) )

pm_anchor <- anchor3_l1_strand_all_comb$anch_plus_bait_minus[[1]][filtered_ids,]
pm_bait <- anchor3_l1_strand_all_comb$anch_plus_bait_minus[[2]][filtered_ids,]

APA_WT_day4_pm<-APA_analysis(
  resolution = 5000,
  anchor_grange = pm_anchor ,
  bait_grange = pm_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)
WT_day4_p_m<-ggAPA(APA_WT_day4_pm$APA, title = "+- 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

WT_day4_p_m


APA_KO_day4_pm<-APA_analysis(
  resolution = 5000,
  anchor_grange = pm_anchor ,
  bait_grange = pm_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)

KO_day4_p_m<-ggAPA(APA_KO_day4_pm$APA, title = "+- 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

KO_day4_p_m


# - + 

mp_int <- cbind(as.data.frame(anchor3_l1_strand_all_comb$anch_minus_bait_plus[[1]]), 
as.data.frame(anchor3_l1_strand_all_comb$anch_minus_bait_plus[[2]]))


filtered_ids<-which(mp_int[,9]+(mp_int[,10]-mp_int[,9])/2>= (mp_int[,2]+ round((mp_int[,3]-mp_int[,2])/2)) )

mp_anchor <- anchor3_l1_strand_all_comb$anch_minus_bait_plus[[1]][filtered_ids,]
mp_bait <- anchor3_l1_strand_all_comb$anch_minus_bait_plus[[2]][filtered_ids,]

APA_WT_day4_mp<-APA_analysis(
  resolution = 5000,
  anchor_grange = mp_anchor ,
  bait_grange = mp_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)
WT_day4_p_m<-ggAPA(APA_WT_day4_mp$APA, title = "-+ 5kb_APA_WT_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)
WT_day4_p_m



APA_KO_day4_mp<-APA_analysis(
  resolution = 5000,
  anchor_grange = mp_anchor ,
  bait_grange = mp_bait,
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)

KO_day4_p_m<-ggAPA(APA_KO_day4_mp$APA, title = "-+ 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

KO_day4_p_m










#ci_WT <- bootstrap_APA(APA_WT_day4_p$M_or, B = 500)
#ci_KO <- bootstrap_APA(APA_WT_day4_p$M_or, B = 500)
#ci_WT; ci_KO


ggplot2::ggsave("./outs/HicAggR/temp/L1p_5kb_WT_day4_APA_p.pdf", plot = WT_day4_p_p, width = 7, height = 7, units = "in")

APA_KO_day4_p<-APA_analysis(
  resolution = 5000,
  anchor_grange = anchor3_l1_std_plus_split$anch_plus_bait_plus[[1]] ,
  bait_grange = anchor3_l1_std_plus_split$anch_plus_bait_plus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = TRUE
)
#297 matrices are oriented
KO_day4_p_p<-ggAPA(APA_KO_day4_p$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.7
)

KO_day4_p_p

ggplot2::ggsave("./outs/HicAggR/temp/L1p_5kb_KO_day4_APA_p.pdf", plot = KO_day4_p_p, width = 7, height = 7, units = "in")

#ci_WT <- bootstrap_APA(APA_WT_day4_p$M_or, B = 500)
#ci_KO <- bootstrap_APA(APA_WT_day4_p$M_or, B = 500)
#ci_WT; ci_KO

ggplot2::ggsave("./outs/HicAggR/temp/L1p_5kb_WT_day4_APA_p.pdf", plot = WT_day4_p_p, width = 7, height = 7, units = "in")


#===========================================
#plus l1 with genes on the right
#===========================================

NROW(split_strand$anch_plus_bait_plus[[1]])
NROW(split_strand$anch_plus_bait_plus[[2]])


APA_WT_day4_p<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_plus_bait_plus[[1]] ,
  bait_grange = split_strand$anch_plus_bait_plus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)

WT_day4_p_p<-ggAPA(APA_WT_day4_p$APA, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.7
)

WT_day4_p_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_WT_day4_APA_p.pdf", plot = WT_day4_p_p, width = 7, height = 7, units = "in")


APA_KO_day4_p<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand$anch_plus_bait_plus[[1]],
  bait_grange = split_strand$anch_plus_bait_plus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)

KO_day4_p_p<-ggAPA(APA_KO_day4_p, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0
  ,colMax = 0.7
)
KO_day4_p_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_KO_day4_APA_p.pdf", plot = KO_day4_p_p, width = 7, height = 7, units = "in")


#========================================
# minus l1 with the genes on the left
#========================================

NROW(split_strand_minus$anch_minus_bait_minus[[1]])
NROW(split_strand_minus$anch_minus_bait_minus[[2]])

split_strand_minus$anch_minus_bait_minus[[1]]


APA_WT_day4_m<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand_minus$anch_minus_bait_minus[[1]] ,
  bait_grange = split_strand_minus$anch_minus_bait_minus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig",
  hic_in=WT_day4,
  orientate = FALSE
)

WT_day4_m_p<-ggAPA(APA_WT_day4_m, title = "5kb_APA_WT_day4 APA"
  ,colMin = 0
  ,colMax = 0.7
)

WT_day4_m_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_WT_day4_APA_m.pdf", plot = WT_day4_m_p, width = 7, height = 7, units = "in")


APA_KO_day4_m<-APA_analysis(
  resolution = 5000,
  anchor_grange = split_strand_minus$anch_minus_bait_minus[[1]] ,
  bait_grange = split_strand_minus$anch_minus_bait_minus[[2]],
  chromSizes=chromSizes,
  path_constraint_bigwig = "in/2024_10_Lara_microC_downstream/func_insulation/balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig",
  hic_in=KO_day4,
  orientate = FALSE
)

KO_day4_m_p<-ggAPA(APA_KO_day4_m, title = "5kb_APA_KO_day4 APA"
  ,colMin = 0
  ,colMax = 0.7
)

KO_day4_m_p

ggplot2::ggsave("./outs/HicAggR/temp/5kb_KO_day4_APA_m.pdf", plot = KO_day4_m_p, width = 7, height = 7, units = "in")

KO_day4_m_p<-ggAPA(APA_KO_day4_m$APA, title = "-- 5kb_APA_KO_day4 APA"
  ,colMin = 0.1
  ,colMax = 0.85
)

KO_day4_m_m

