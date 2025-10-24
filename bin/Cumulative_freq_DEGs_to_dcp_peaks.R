#!/usr/bin/env Rscript
# ======================================================================
# Script: cumulative_distance_analysis_fixed_paths.R
# Purpose: Same as the parametric version, but all paths are fixed
#          within the script. Edit the "USER SETTINGS" section below
#          to change input/output paths.
#
# HOW TO RUN THIS SCRIPT IN DOCKER (rocker/tidyverse:4.5.1)
# ----------------------------------------------------------------------
# 1) Install Docker: https://docs.docker.com/get-docker/
# 2) Pull the Rocker image with tidyverse (includes ggplot2 + data.table):
#      docker pull rocker/tidyverse:4.5.1
# 3) From your project folder (where this script and your data live), run:
#      sudo docker run -it -v /media/lucio/easystore:/media/lucio/easystore -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
#        Rscript Cumulative_freq_DEGs_to_dcp_peaks.R
#
# Notes:
#   - Requirements: data.table, ggplot2
#   - Local sources required:
#       ./git/downstream_multiomic/bin/hypergeometric.R
#       ./git/downstream_multiomic/bin/distance_functions.R
#       ./git/bioinfo_generics/base/bin/fun_compare_curves.R
# ======================================================================


# ------------------------------ USER SETTINGS ---------------------------
IN_RNASEQ_DIR <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/build38_DEseq2_RNAseq"
IN_GREAT_DIR  <- "git/MLL2_LINE1_long-range_regulation/assets/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp"
IN_ANNO_DIR   <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/build38_DEseq2_RNAseq"
REFSEQ_PATH   <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/build38_mm10_ncbiRefSeqCurated.txt.gz"

# Peaks (near LINE-1) files
PEAKS_DCM_NEAR_FILE <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm.target.peaks.bed"
PEAKS_DCP_NEAR_FILE <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcp.target.peaks.bed"

# Output directory
OUTDIR <- "outs"

# Function sources
SRC_HYPER  <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/hypergeometric.R"
SRC_DIST   <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/distance_functions.R"
SRC_CURVES <- "./git/MLL2_LINE1_long-range_regulation/pipelines/bioinfo_generics/base/bin/fun_compare_curves.R"

# ------------------------------ PREP I/O --------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "cumulative_dist_degsVSdcp"), showWarnings = FALSE, recursive = TRUE)

#==============================
# Environment
#==============================

#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

library(ggplot2)

allpreaks=FALSE

#=======================================
#  functions
#=======================================

source(SRC_HYPER)
source(SRC_DIST)
source(SRC_CURVES)

exportBed<-function(distancePlotResult, file){
  
  D0_down_pro<-distancePlotResult
  
  D0_down_pro
  chr<-c()
  start<-c()
  end<-c()
  direction<-c()
  
  for (i in 1:nrow(D0_down_pro$dis.CpG.plus.degs)){
    D0_down_pro$dis.CpG.plus.degs[i,]
    chr[i]<-D0_down_pro$dis.CpG.plus.degs[i,1]
    
    if(D0_down_pro$dis.CpG.plus.degs[i,"from_position"]<=D0_down_pro$dis.CpG.plus.degs[i,"to_position"]){
      start[i]<-round(D0_down_pro$dis.CpG.plus.degs[i,"from_position"],0)
      end[i]<-round(D0_down_pro$dis.CpG.plus.degs[i,"to_position"],0)
      direction[i]<-"peak.right"
    }else{
      start[i]<-round(D0_down_pro$dis.CpG.plus.degs[i,"to_position"],0)
      end[i]<-round(D0_down_pro$dis.CpG.plus.degs[i,"from_position"],0)
      direction[i]<-"peak.left"
    }
  }
  
  D0_down_pcm_bed<-cbind(chr, start, end, direction)
  
  write.table(D0_down_pcm_bed, file=file, 
              sep="\t",
              quote=FALSE,
              col.names = FALSE,
              row.names = FALSE
  )
}


#================================
# inputs
#================================

degs_D4<-read.table(file.path(IN_RNASEQ_DIR, "D4_DKO_VS_WT.deseq2.results.tsv"), header=TRUE, sep="\t")

degs_D0<-read.table(file.path(IN_RNASEQ_DIR, "D0_DKO_vs_WT.deseq2.results.tsv"), header=TRUE, sep="\t")

degs_wtD4D0<-read.table(file.path(IN_RNASEQ_DIR, "D4WT_VS_D0WT.deseq2.results.tsv"), header=TRUE, sep="\t")

allanno.genes<-read.table("git/MLL2_LINE1_long-range_regulation/assets/GREATv4.genes.mm10.tsv",sep="\t",header=FALSE)

great.dcm<-read.table(file.path(IN_GREAT_DIR, "dcm", "20241126-public-4.0.4-hpxvmj-mm10-all-gene.txt"), sep="\t", header=FALSE)

great.dcp<-read.table(file.path(IN_GREAT_DIR, "dcp", "20241128-public-4.0.4-qlnior-mm10-all-gene.txt"), sep="\t", header=FALSE)

#nfcore RNAseq annotation extracted from gtf
anno<-read.table(file.path(IN_ANNO_DIR, "mm10.anno.tsv"), header=TRUE, sep="\t")

#gtf of annotation used as input also for RNAseq analysis
#anno.gtf<-read.table("in/build38_DEseq2_RNAseq/mm10.refGene.gtf.gz", sep="\t")
#anno.gtf<-anno.gtf[which(anno.gtf$V3=="CDS"),]

#TSS annotated reference
refseq.anno<-read.table(REFSEQ_PATH, sep="\t")
refseq.anno[,"V5"]<-refseq.anno[,"V5"]+1
refseq.anno<-refseq.anno[,c("V13","V3","V5","V6","V4")]

#REMOVE the genes strand "+" that have same gene,start and the genes strand "-" gene, end
refseq.anno.stplus<-refseq.anno[which(refseq.anno$V4=="+"),]
refseq.anno.stplus<-refseq.anno.stplus[which(duplicated(refseq.anno.stplus[,c("V13","V3", "V5")])==FALSE),]

refseq.anno.stmin<-refseq.anno[which(refseq.anno$V4=="-"),]
refseq.anno.stmin<-refseq.anno.stmin[which(duplicated(refseq.anno.stmin[,c("V13","V3", "V6")])==FALSE),]

refseq.anno<-rbind(refseq.anno.stplus, refseq.anno.stmin)

#add ensembl to refseq

refseq.anno<-merge(anno[,c("gene_name","gene_id")], refseq.anno, all.y=TRUE, by=1)
refseq.anno<-refseq.anno[which(duplicated(refseq.anno)==FALSE),]

refseq.anno<-refseq.anno[which(is.na(refseq.anno$gene_id)==FALSE),]

#unid<-rep(NA,nrow(refseq.anno))
refseq.anno<-refseq.anno[grep("_.", refseq.anno$V3, invert = TRUE),]
#unid[which(duplicated(refseq.anno$gene_name)==FALSE)] <-refseq.anno[which(duplicated(refseq.anno$gene_name)==FALSE),""] 

refseq.anno.reform<-data.frame(
  GeneID=refseq.anno$gene_id,
  unid=paste0(refseq.anno$gene_id,"_",row.names(refseq.anno)),
  GeneName=refseq.anno$gene_name,
  Chromosome=refseq.anno$V3,
  start=refseq.anno$V5,
  end=refseq.anno$V6,
  strand=refseq.anno$V4
)


norm.counts<-read.table(file.path(IN_RNASEQ_DIR, "all.normalised_counts.tsv"),
                        sep="\t",
                        header=TRUE)

#gene.annWithPosition<-read.table("in/ucsc/martexport.txt", sep="\t", header=TRUE)

GeneID=c()
GeneIDversion=c()
GeneName=c()
Chromosome=c()
start=c()
end=c()
strand=c()


gene.annWithPosition <- refseq.anno.reform
colnames(gene.annWithPosition)<-c( 
                                "GeneID",
                                "uniq",
                                "GeneName",
                                "Chromosome",
                                 "start",
                                 "end",
                                 "strand"
                                 )


annobed<-data.frame(chr=gene.annWithPosition$Chromosome , 
                    start=gene.annWithPosition$start , 
                    end=gene.annWithPosition$end, 
                    gene=paste(gene.annWithPosition$GeneName,gene.annWithPosition$strand),
                    strand=gene.annWithPosition$strand
                    )

write.table(annobed, file.path(OUTDIR, "build38_mm10_ncbiRefSeqCurated.bed"),
            sep="\t",
            col.names = FALSE,
            row.names = FALSE,
            quote= FALSE
)


unique(gene.annWithPosition$Chromosome)

peaks.near.lines.file <- PEAKS_DCM_NEAR_FILE

dcp_peaks.near.lines.file <- PEAKS_DCP_NEAR_FILE


#===========================================================

dcm_peak.near.lines<-read.table(peaks.near.lines.file,
                                sep="\t",
                                header=FALSE
)

dcp_peaks.near.lines<-read.table(dcp_peaks.near.lines.file,
                                 sep="\t",
                                 header=FALSE
)

if (allpreaks==TRUE){


  pcm_peaks.near.lines<-read.table(pcm_peaks.near.lines.file,
                                   sep="\t",
                                   header=FALSE
  )

  pcp_peaks.near.lines<-read.table(pcp_peaks.near.lines.file,
                                   sep="\t",
                                   header=FALSE
  )
}


#==========add annotation to DEGs

degs_D4_ann_diff_down<-addAnnotationToDegs(degs_D4, gene.annWithPosition, direction = "down", log2fc = -1)

degs_D0_ann_diff_down<-addAnnotationToDegs(degs_D0, gene.annWithPosition, direction = "down", log2fc = -1)

degs_wtD4D0_ann_diff_down<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, direction = "down", log2fc = -1  )

degs_D4_ann_diff_up<-addAnnotationToDegs(degs_D4, gene.annWithPosition, 
                                         direction = "up",
                                         log2fc = 1,
                                         )

degs_D0_ann_diff_up<-addAnnotationToDegs(degs_D0, gene.annWithPosition, 
                                         direction = "up",
                                         log2fc = 1
                                         )

degs_wtD4D0_ann_diff_up<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, direction = "up", log2fc = 1  )

degs_wtD4D0_ann_nodiff<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, log2fc = 1  )

degs_D4_ann_diff_all<-addAnnotationToDegs(degs_D4, gene.annWithPosition, 
                                         direction = "all"
)

degs_D0_ann_diff_all<-addAnnotationToDegs(degs_D0, gene.annWithPosition, 
                                         direction = "all"
)

allgenes_D0_ann_diff_all<-addAnnotationToDegs(degs_D0, gene.annWithPosition,
                                              padj = 1,
                                              log2fc = 0,
                                              direction = "all"
)

allgenes_D4_ann_diff_all<-addAnnotationToDegs(degs_D4, gene.annWithPosition,
                                              padj = 1,
                                              log2fc = 0,
                                              direction = "all"
)

degs_wtD4D0_ann_diff<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, 
                                           direction = "all",
                                           log2fc = 1  )

degs_wtD4D0_ann_allgenes<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, 
                                           direction = "all",
                                           padj = 1,
                                           log2fc = 0)

degs_wtD4D0_ann_allgenes<-addAnnotationToDegs( degs_wtD4D0, gene.annWithPosition, 
                                               direction = "all",
                                               padj = 1,
                                               log2fc = 0)


nodiff_D0 <- allgenes_D0_ann_diff_all[which(allgenes_D0_ann_diff_all$padj > 0.05),]
nodiff_D4 <- allgenes_D4_ann_diff_all[which(allgenes_D4_ann_diff_all$padj > 0.05),]

degs_D4_ann_diff_down$ensembl
degs_wtD4D0_ann_diff$ensembl

#[which( degs_wtD4D0_ann_diff$ensembl %in% degs_D4_ann_diff_down$ensembl),]

overlap_down.torem<-intersect(degs_wtD4D0_ann_diff$ensembl, degs_D4_ann_diff_down$ensembl)
D4_no_wtD4_D0_down<-degs_D4_ann_diff_down[which(!degs_D4_ann_diff_down$ensembl %in% overlap_down.torem),]

overlap_up.torem<-intersect(degs_wtD4D0_ann_diff$ensembl, degs_D4_ann_diff_up$ensembl)
D4_no_wtD4_D0_up<-degs_D4_ann_diff_up[which(!degs_D4_ann_diff_up$ensembl %in% overlap_up.torem),]


#overlap_nodiff.torem<-intersect(degs_wtD4D0_ann_diff$ensembl, nodiff_D4$ensembl)

#==========================
#        D0
#==========================

#=====distal density

#  


dcm_peak.near.lines

if (allpreaks==TRUE){
  
  dcp_peaks.near.lines
  pcm_peaks.near.lines
  pcp_peaks.near.lines
  
  D0_down_pro<-distance_plot(minus_peaks=pcm_peaks.near.lines, 
                             plus_peaks = pcp_peaks.near.lines  ,
                             degs = degs_D0_ann_diff_down,
                             title="proximal peaks near line D0 down",
                             window = 5000,
                             window.type = "inside"
  )
  
  exportBed(D0_down_pro, file=file.path(OUTDIR, "D0_down_pcp_L1near.bed"))  
  
  D0_up_pro<-distance_plot(minus_peaks=pcm_peaks.near.lines, 
                           plus_peaks = pcp_peaks.near.lines,
                           degs = degs_D0_ann_diff_up,
                           title="proximal peaks near line D0 up",
                           window = 5000,
                           window.type = "inside"
  )
  
}


D0_down_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                               plus_peaks = dcp_peaks.near.lines ,
                               degs=degs_D0_ann_diff_down,
                               title="distal peaks near line D0 down",
                               window = 1000000,
                               window.type = "inside"
)

D0_nodiff_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                                 plus_peaks = dcp_peaks.near.lines ,
                                 degs=nodiff_D0,
                                 title="distal peaks near line D0 all annotated genes",
                                 window = 1000000,
                                 window.type = "inside"
)

D0_up_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                         plus_peaks = dcp_peaks.near.lines ,
                         degs=degs_D0_ann_diff_up,
                         title="distal peaks near line D0 up",
                         window = 1000000,
                         window.type = "inside"
)

df_windowed<-df_all<-as.data.frame(
  rbind(

    cbind(distances=D0_down_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_DOWN")
    ),

    cbind(distances=D0_up_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_UP")
    ),
    
    cbind(distances=D0_nodiff_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_NODIFF")
    )
  )
)





#cumulative distribution plot

df_windowed$distances<-as.numeric(df_windowed$distances)

ggp1<-ggplot(df_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D0 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.minus_DOWN" = "#3f3092",
                                "dis.CpG.minus_UP" ="#608dd1",
                                "dis.CpG.minus_NODIFF" = "#bdbbba"
  ))

ggsave(file.path(OUTDIR, "cumulative_dist_degsVSdcp", "D0_dis_downDEGs_dcm.png"), ggp1,
       width = 8, height = 6, dpi = 300, units = "in")


grps<-c("dis.CpG.minus_DOWN",   "dis.CpG.minus_UP", "dis.CpG.minus_NODIFF")

kmatrix<-matrix(nrow=length(grps), ncol=length(grps))
colnames(kmatrix)<-grps
rownames(kmatrix)<-grps

# dcm DOWN vs dcm UP
a<- grps[1] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[2] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

# dcm DOWN vs no diff
a<- grps[1] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[3] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

# dcm UP vs no diff
a<- grps[2] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[3] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

a<- grps[1] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[1] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

a<- grps[2] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[2] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

a<- grps[3] #df_windowed[which( df_windowed$group== grps[1] ), ]
b<- grps[3] #df_windowed[which( df_windowed$group== grps[2] ), ]
kmatrix<-compair.curve(df=df_windowed,a,b,kmatrix,grps)

write.table(kmatrix, file=file.path(OUTDIR, "cumulative_dist_degsVSdcp", "D0_dis_downDEGs_dcm_pval_wilcox.tsv"), col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")


#====================================================================================================================
#====================================================================================================================
#
#           cumulative distribution differentially expressed gene to nearest MLL2-bound L1 element
#
#====================================================================================================================
#====================================================================================================================



D4_down_dis_win<-distance_plot(minus_peaks = dcm_peak.near.lines, 
                               plus_peaks = dcp_peaks.near.lines ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

D4_nodiff_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                                 plus_peaks = dcp_peaks.near.lines ,
                                 degs=nodiff_D4,
                                 title="distal peaks near line D4 all annotated genes",
                                 window = 1000000,
                                 window.type = "inside"
)


# tt_D4_nodiff_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
#                                  plus_peaks = dcp_peaks.near.lines ,
#                                  degs=nodiff_D4,
#                                  title="distal peaks near line D4 all annotated genes",
#                                  window = NULL,
#                                  window.type = "inside"
# )

nrow(nodiff_D4)
#nrow(tt_D4_nodiff_dis_win$dis.CpG.plus.degs)
#nrow(tt_D4_nodiff_dis_win$dis.CpG.minus.degs)

#setdiff(tt_D4_nodiff_dis_win$dis.CpG.minus.degs$from_identifier, nodiff_D4$ensembl)


D4_up_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                             plus_peaks = dcp_peaks.near.lines ,
                             degs=degs_D4_ann_diff_up,
                             title="distal peaks near line D4 up",
                             window = 1000000,
                             window.type = "inside"
)


#D4_down_induced
#D4_up_induced
#d4_nodiff_induced


df_d4_windowed<-as.data.frame(
  rbind(

    cbind(distances=D4_down_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_DOWN")
    ),
    
    cbind(distances=D4_up_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_UP")
    ),

    cbind(distances=D4_nodiff_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_NODIFF")
    )
  )
)


ggplot(df_d4_windowed, aes(as.numeric(as.character(distances)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic()


nodiff_D4 <- allgenes_D4_ann_diff_all[which(allgenes_D4_ann_diff_all$padj > 0.05),]

overlap_downdko_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, degs_D4_ann_diff_down$ensembl)


#overlap_down.torem <- intersect(degs_wtD4D0_ann_diff$ensembl, degs_D4_ann_diff_down$ensembl)

#D4_no_wtD4_D0_down <- degs_D4_ann_diff_down[which(!degs_D4_ann_diff_down$ensembl %in% overlap_down.torem),]
#D4_in_wtD4_D0_down <- degs_D4_ann_diff_down[which( degs_D4_ann_diff_down$ensembl %in% overlap_down.torem ),]

D4_downdko_upwt <- degs_D4_ann_diff_down[which( degs_D4_ann_diff_down$ensembl %in% overlap_downdko_upwt ),]

D4_down_wtD4D0_up_dis<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                                    plus_peaks = dcp_peaks.near.lines ,
                                    degs=D4_downdko_upwt,
                                    title="distal peaks near line",
                                    window = 1000000,
                                    window.type = "inside"
)

overlap_updko_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, degs_D4_ann_diff_up$ensembl)
D4_updko_upwt <- degs_D4_ann_diff_up[which( degs_D4_ann_diff_up$ensembl %in% overlap_updko_upwt ),]

D4_up_wtD4D0_up_dis<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                                    plus_peaks = dcp_peaks.near.lines ,
                                    degs= D4_updko_upwt,
                                    title="distal peaks near line",
                                    window = 1000000,
                                    window.type = "inside"
)

overlap_undiff_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, nodiff_D4$ensembl)
D4_undiff_upwt <- nodiff_D4[which( nodiff_D4$ensembl %in% overlap_undiff_upwt ),]

D4_undiff_wtD4D0_up_dis<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                                    plus_peaks = dcp_peaks.near.lines ,
                                    degs= D4_undiff_upwt,
                                    title="distal peaks near line",
                                    window = 1000000,
                                    window.type = "inside"
)

df_windowed<-as.data.frame(
  rbind(

    cbind(distances=D4_down_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_DOWN")
    ),

    cbind(distances=D4_up_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_UP")
    ),
    
    cbind(distances=D4_nodiff_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_NODIFF")
    ),

    cbind(distances=D4_down_wtD4D0_up_dis$dis.CpG.plus.degs$distances,
          group=rep("induced.dis.CpG.plus_DOWN")
    )

    # cbind(distances=D4_up_wtD4D0_up_dis$dis.CpG.plus.degs$distances ,
    #       group=rep("induced.dis.CpG.plus_UP")
    # ),

    # cbind(distances=D4_undiff_wtD4D0_up_dis$dis.CpG.plus.degs$distances ,
    #       group=rep("induced.dis.CpG.plus_NODIFF")
    # )

    
  )
)


#cumulative distribution plot

df_windowed$distances<-as.numeric(df_windowed$distances)

ggp1<-ggplot(df_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.plus_DOWN" = "#3f3092",
                                "dis.CpG.plus_UP" ="#608dd1",
                                "dis.CpG.plus_NODIFF" = "#bdbbba",
                                "induced.dis.CpG.plus_DOWN" = "red"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"



  ))                                       

ggsave(file.path(OUTDIR, "cumulative_dist_degsVSdcp", "D4_dis_downDEGs_dcm.png"), ggp1,
       width = 8, height = 6, dpi = 300, units = "in")



df_windowed<-as.data.frame(
  rbind(

    cbind(distances=D4_down_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_DOWN")
    ),

    cbind(distances=D4_up_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_UP")
    ),
    
    cbind(distances=D4_nodiff_dis_win$dis.CpG.minus.degs$distances,
          group=rep("dis.CpG.minus_NODIFF")
    ),

    cbind(distances=D4_down_wtD4D0_up_dis$dis.CpG.minus.degs$distances,
          group=rep("induced.dis.CpG.minus_DOWN")
    )

    # cbind(distances=D4_up_wtD4D0_up_dis$dis.CpG.plus.degs$distances ,
    #       group=rep("induced.dis.CpG.plus_UP")
    # ),

    # cbind(distances=D4_undiff_wtD4D0_up_dis$dis.CpG.plus.degs$distances ,
    #       group=rep("induced.dis.CpG.plus_NODIFF")
    # )

    
  )
)

df_windowed$distances<-as.numeric(df_windowed$distances)

table(df_windowed$group)

ggp1<-ggplot(df_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D0 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.minus_DOWN" = "#3f3092",
                                "dis.CpG.minus_UP" ="#608dd1",
                                "dis.CpG.minus_NODIFF" = "#bdbbba",
                                "induced.dis.CpG.minus_DOWN" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"



  )) 

