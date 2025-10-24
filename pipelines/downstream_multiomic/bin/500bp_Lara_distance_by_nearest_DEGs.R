#==============================
# Environment
#==============================

#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

library(ggplot2)

allpreaks=FALSE

#=======================================
#  functions
#=======================================

source("./git/downstream_multiomic/bin/hypergeometric.R")
source("./git/downstream_multiomic/bin/distance_functions.R")
source ("./git/bioinfoGenerals/base/bin/fun_compare_curves.R")

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

degs_D4<-read.table("in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv", header=TRUE, sep="\t")

degs_D0<-read.table("in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv", header=TRUE, sep="\t")

degs_wtD4D0<-read.table("in/build38_DEseq2_RNAseq/D4WT_VS_D0WT.deseq2.results.tsv", header=TRUE, sep="\t")

allanno.genes<-read.table("in/GREATv4/GREATv4.genes.mm10.tsv",sep="\t",header=FALSE)

great.dcm<-read.table("outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcm/20241126-public-4.0.4-hpxvmj-mm10-all-gene.txt", sep="\t", header=FALSE)

great.dcp<-read.table("outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcp/20241128-public-4.0.4-qlnior-mm10-all-gene.txt", sep="\t", header=FALSE)

#nfcore RNAseq annotation extracted from gtf
anno<-read.table("in/build38_DEseq2_RNAseq/mm10.anno.tsv", header=TRUE, sep="\t")

#gtf of annotation used as input also for RNAseq analysis
#anno.gtf<-read.table("in/build38_DEseq2_RNAseq/mm10.refGene.gtf.gz", sep="\t")
#anno.gtf<-anno.gtf[which(anno.gtf$V3=="CDS"),]

#TSS annotated reference
refseq.anno<-read.table("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/in/ucsc/build38_mm10_ncbiRefSeqCurated.txt.gz", sep="\t")
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


norm.counts<-read.table("in/build38_DEseq2_RNAseq/all.normalised_counts.tsv",
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

#print(nrow(anno.gtf))
# for(i in 1:nrow(anno.gtf)){
#   
#   if(i %% 10000 == 0){
#     cat("i: ",i,"\n")
#   }
#   
#   meta<-strsplit(anno.gtf$V9[i],"; ")[[1]]
#   
#   GeneID[i]<-gsub("gene_id ", "" ,meta[grep("gene_id", meta)])
#   GeneIDversion[i]<-paste(
#     GeneID[i],
#     gsub("gene_version ", "" ,meta[grep("^gene_version", meta)]),
#     sep="."
#   )
#   GeneName[i]<-gsub("gene_name ", "" ,meta[grep("gene_name", meta)])
#   
#   Chromosome[i]<-anno.gtf$V1[i]
#   start[i]<-anno.gtf$V4[i]
#   end[i]<-anno.gtf$V5[i]
#   strand[i]<-anno.gtf$V7[i]
#   
# }
# 
# gene.annWithPosition <- data.frame(
#   GeneID,
#   GeneIDversion,
#   GeneName,
#   Chromosome,
#   start,
#   end,
#   strand
# )



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

write.table(annobed, "in/ucsc/build38_mm10_ncbiRefSeqCurated.bed",
            sep="\t",
            col.names = FALSE,
            row.names = FALSE,
            quote= FALSE
)


unique(gene.annWithPosition$Chromosome)

peaks.near.lines.file<-"outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm.target.peaks.bed"

dcp_peaks.near.lines.file<-"outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcp.target.peaks.bed"

#pcm_peaks.near.lines.file<-"outs/distfilter_DKO_K4me3_pcm.target.peaks.bed"

#pcp_peaks.near.lines.file<-"outs/distfilter_DKO_K4me3_pcp.target.peaks.bed"

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

#TODO modify this and change variables dko.dis.CpG.minus.peak instead of change name to lina1.peaks  

#dko.dis.CpG.minus.peak<-peaks.near.lines
#dko.dis.CpG.plus.peak
#dko.pro.CpG.plus.peak
#dko.pro.CpG.minus.peak

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
  
  exportBed(D0_down_pro, file="./outs/D0_down_pcp_L1near.bed")  
  
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


df_windowed$distances<-as.numeric(df_windowed$distances)


#cumulative distribution plot

ggp1<-ggplot(df_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D0 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.minus_DOWN" = "#3f3092",
                                "dis.CpG.minus_UP" ="#608dd1",
                                "dis.CpG.minus_NODIFF" = "#bdbbba"
  ))

ggsave("outs/distribution_distances_downDEGS_dcm/dis_downDEGs_dcm.png", ggp1,
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

write.table(kmatrix, file="outs/distribution_distances_downDEGS_dcm/dis_downDEGs_dcm_pval_wilcox.tsv", col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")

#======================================
# D4
#======================================

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

D4_up_dis_win<-distance_plot(minus_peaks=dcm_peak.near.lines, 
                             plus_peaks = dcp_peaks.near.lines ,
                             degs=degs_D4_ann_diff_up,
                             title="distal peaks near line D4 up",
                             window = 1000000,
                             window.type = "inside"
)

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



# wtD4D0_undiff<-distance_plot(minus_peaks=dcm_peak.near.lines, 
#                                            plus_peaks = dcp_peaks.near.lines ,
#                                            degs=degs_wtD4D0_ann_nodiff,
#                                            title="distan peaks near line no diff genes",
#                                            window = 1000000,
#                                            window.type = "inside"
# )
# 
# df_wt_up_dko_diff$distances<-as.numeric(df_wt_up_dko_diff$distances)
# 
# sample_size = df_wt_up_dko_diff %>% group_by(group) %>% summarize(num=n())
# 
# ggplot(df_wt_up_dko_diff, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
#   labs(title="Empirical Cumulative \n Density Function",
#        y = "perc", x="distance")+
#   theme_classic()


#======================================
#  genes that overlap with great
#======================================

# degs_D0_ann_diff_down
# 
# dcm_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/dcm/20240719-public-4.0.4-nxvS8V-mm10-all-gene.txt", sep="\t",header=FALSE)
# 
# dcp_nearest_ggenes<-read.table("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/outs/great/filtered_by_distance/nearest/dcp/20240719-public-4.0.4-qoHtEr-mm10-all-gene.txt", 
#                                sep="\t",
#                                header=FALSE
# )
# 
# pcp_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/pcp/20240719-public-4.0.4-9RkkYc-mm10-all-gene.txt",
#                                sep="\t",
#                                header=FALSE                             
# )
# 
# 
# pcm_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/pcm/20240719-public-4.0.4-Zd7s4R-mm10-all-gene.txt",
#                                sep="\t",
#                                header=FALSE                             
# )
# 
# degs.anno<-read.table("in/DEseq2_RNAseq/fixed_Mus_musculus.anno.tsv", header=TRUE, sep="\t")
# degs.anno<-cbind(ensembl=degs.anno$gene_id, geneName=degs.anno$gene_name)
# degs.anno<-degs.anno[which(duplicated(degs.anno)==FALSE),]
# 
# degs_D0<-read.table("in/DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_D0<-merge(degs.anno,degs_D0,  by=1)
# 
# degs_D4<-read.table("in/DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_D4<-merge(degs.anno,degs_D4,  by=1)
# 
# degs_WT<-read.table("in/DEseq2_RNAseq/D4WT_VS_D0WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_WT<-merge(degs.anno,degs_WT,  by=1)
# 
# D0.up<-extract_degs(degs_D0, direction="up" )
# D0.down<-extract_degs(degs_D0, direction="down" )
# D0.unch<-extract_degs(degs_D0, pval = NULL )
# 
# D4.up<-extract_degs(degs_D4, direction="up" )
# D4.down<-extract_degs(degs_D4, direction="down" )
# D4.unch<-extract_degs(degs_D4, pval = NULL )
# 
# WT.sign.down<-extract_degs(degs_WT, direction="down" )
# WT.sign.up<-extract_degs(degs_WT, direction="up" )
# 
# D0.down.dcp.overlap<-hyper_degs_tgenes(degs=D0.down, 
#                   targetGenes=dcp_nearest_ggenes, 
#                   anno=allanno.genes, 
#                   extract="overlap")
# 
# D0.down.dcm.overlap<-hyper_degs_tgenes(degs=D0.down, 
#                                        targetGenes=dcm_nearest_ggenes, 
#                                        anno=allanno.genes, 
#                                        extract="overlap")
# 
# D0.down.pcm.overlap<-hyper_degs_tgenes(degs=D0.down, 
#                                        targetGenes=pcm_nearest_ggenes, 
#                                        anno=allanno.genes, 
#                                        extract="overlap")
# 
# D0.down.pcp.overlap<-hyper_degs_tgenes(degs=D0.down, 
#                                        targetGenes=pcp_nearest_ggenes, 
#                                        anno=allanno.genes, 
#                                        extract="overlap")
# 
# D0.down.dcp.overlap<-degs_D0_ann_diff_down[which (degs_D0_ann_diff_down$gene %in% D0.down.dcp.overlap),]
# D0.down.dcm.overlap<-degs_D0_ann_diff_down[which (degs_D0_ann_diff_down$gene %in% D0.down.dcm.overlap),]
# 
# D0.down.pcm.overlap<-degs_D0_ann_diff_down[which (degs_D0_ann_diff_down$gene %in% D0.down.pcm.overlap),]
# D0.down.pcp.overlap<-degs_D0_ann_diff_down[which (degs_D0_ann_diff_down$gene %in% D0.down.pcp.overlap),]
# 
# pro.CpG.minus.degs_D0<-extract_distances(peak.coords =  pcm_peaks.near.lines ,
#                                          startend.coords = D0.down.pcm.overlap[,c(3,4,5,2,7,6)],
#                                          from="startend",
#                                          to="peak",
#                                          peak.coordtype = "centered",
#                                          startend.coordtype = "fprime",
#                                          to_name="to",
#                                          from_name="from"
# )
# 
# pro.CpG.plus.degs_D0<-extract_distances(peak.coords =  pcp_peaks.near.lines ,
#                                          startend.coords = D0.down.pcp.overlap[,c(3,4,5,2,7,6)],
#                                          from="startend",
#                                          to="peak",
#                                          peak.coordtype = "centered",
#                                          startend.coordtype = "fprime",
#                                          to_name="to",
#                                          from_name="from"
# )
# 
# 
# dis.CpG.minus.degs_D0<-extract_distances(peak.coords =  dcm_peak.near.lines ,
#                                          startend.coords = D0.down.dcm.overlap[,c(3,4,5,2,7,6)],
#                                          from="startend",
#                                          to="peak",
#                                          peak.coordtype = "centered",
#                                          startend.coordtype = "fprime",
#                                          to_name="to",
#                                          from_name="from"
# )
# 
# dis.CpG.plus.degs_D0<-extract_distances(peak.coords =  dcp_peaks.near.lines ,
#                                          startend.coords = D0.down.dcp.overlap[,c(3,4,5,2,7,6)],
#                                          from="startend",
#                                          to="peak",
#                                          peak.coordtype = "centered",
#                                          startend.coordtype = "fprime",
#                                          to_name="to",
#                                          from_name="from"
# )
# 
# 
# D0_down_dis.overlap<-distance_plot(minus_peaks=dcm_peak.near.lines,
#                            plus_peaks = dcp_peaks.near.lines ,
#                            degs=D0.down.dcm.overlap,
#                            title="distal peaks near line D0 down"
# )
# 
# 

#======================================
#
#======================================

# library(plotly)
# library(ggplot2)
# library(ggdendro)
# 
# model <- hclust(dist(similarity.matrix.d), "ave")
# dhc <- as.dendrogram(model)
# 
# data <- dendro_data(dhc, type = "rectangle")
# 
# p <- ggdendrogram(dhc, rotate = FALSE, size = 2)
# 
# ggplotly(p)
# 
# shapiro.test(df[which(df$group=="dis.CpG.minus_DOWN"),"distances"])
# shapiro.test(df[which(df$group=="dis.CpG.plus_DOWN"),"distances"])
# shapiro.test(df[which(df$group=="pro.CpG.minus_DOWN"),"distances"])
# shapiro.test(df[which(df$group=="pro.CpG.plus_DOWN"),"distances"])
# shapiro.test(df[which(df$group=="dis.CpG.minus_NOTDIFF"),"distances"])
# shapiro.test(df[which(df$group=="dis.CpG.plus_NOTDIFF"),"distances"])
# shapiro.test(df[which(df$group=="pro.CpG.minus_NOTDIFF"),"distances"])
# shapiro.test(df[which(df$group=="pro.CpG.plus_NOTDIFF"),"distances"])
# 
# krut<-kruskal.test(distances ~ group, data = df)
# 
# pwilcox<-pairwise.wilcox.test(df$distances, df$group,
#                      p.adjust.method = "BH",
#                      paired=FALSE
#                      )
# 
# 
# model <- hclust(dist(pwilcox$p.value), "ave")
# dhc <- as.dendrogram(model)
# 
# data <- dendro_data(dhc, type = "rectangle")
# 
# p <- ggdendrogram(dhc, rotate = FALSE, size = 2)
# 
# ggplotly(p)
# 
# #==================================================
# #                valuate GREAT genes
# #==================================================
# 
# #instead of degs i want use genes predicted by GREAT 
# 
# 
# #allanno.genes
# #great.dcm
# #great.dcp
# 
# great.dcm.anno<-merge(allanno.genes, great.dcm, by.x=5,by.y=1,all.y==TRUE)
# great.dcp.anno<-merge(allanno.genes, great.dcp, by.x=5,by.y=1,all.y==TRUE)
# 
# #TODO make this better 
# 
# great.dcm.anno<-data.frame(ensembl=great.dcm.anno$V1, 
#                            gene=great.dcm.anno$V5, 
#                            chromosome = great.dcm.anno$V2.x, 
#                            start = great.dcm.anno$V3, 
#                            end = great.dcm.anno$V3+1, 
#                            strand=great.dcm.anno$V4,
#                            val1=great.dcm.anno$V4,
#                            val2=great.dcm.anno$V4
#                            )
# 
# great.dcp.anno<-data.frame(ensembl = great.dcp.anno$V1, 
#                            gene = great.dcp.anno$V5, 
#                            chromosome = great.dcp.anno$V2.x, 
#                            start = great.dcp.anno$V3, 
#                            end = great.dcp.anno$V3+1, 
#                            strand = great.dcp.anno$V4,
#                            val1 = great.dcp.anno$V4,
#                            val2 = great.dcp.anno$V4
#                            )
# 
# allgenes_D0_ann_diff_all
# 
# #nrow(great.dcm.anno)
# #nrow(great.dcm.anno[great.dcm.anno$V1 %in% allgenes_D0_ann_diff_all$ensembl,])
# 
# #head(D4_no_wtD4_D0_down[,c(3,4,5,2,7,6)])
# 
# # D4_no_wtD4_D0_down_pro<-distance_plot(minus_peaks=dcm_peak.near.lines , 
# #                                       plus_peaks =dcp_peak.near.lines  ,
# #                                       degs=,
# #                                       title="distal peaks near line D0 all annotated genes"
# # )
# 
# 
# #great.dcm.anno<-merge(great.dcm.anno, anno, by.x=2, by.y=1, fill.x=TRUE, fill.y=FALSE)
# 
# #great.dcm.anno<-great.dcm.anno[,c(1,18,7,8,9,5,10,27)]
# 
# 
# #
# dis.CpG.minus.great_D0<-extract_distances(peak.coords = dcm_peak.near.lines,
#                                          startend.coords = great.dcm.anno[,c(3,4,5,2,7,6)],
#                                          from="startend",
#                                          to="peak",
#                                          peak.coordtype = "centered",
#                                          startend.coordtype = "fprime",
#                                          to_name="to",
#                                          from_name="from"
#                                         )
# 
# dis.CpG.plus.great_D0<-extract_distances(peak.coords = dcp_peaks.near.lines,
#                                           startend.coords = great.dcp.anno[,c(3,4,5,2,7,6)],
#                                           from="startend",
#                                           to="peak",
#                                           peak.coordtype = "centered",
#                                           startend.coordtype = "fprime",
#                                           to_name="to",
#                                           from_name="from"
# )
# 
# sign<-t.test( dis.CpG.plus.great_D0$distances ,  dis.CpG.minus.great_D0$distances, 
#               alternative = "two.sided", var.equal = FALSE)
# 
# plot(density(dis.CpG.plus.great_D0$distances), 
#      main=paste0(
#                  " KStest pval = ",
#                  sign$p.value,
#                  " D = ",
#                  round(sign$statistic,3)
#       )
#     ) 
# lines(density(dis.CpG.minus.great_D0$distances), col = "red") 
# legend("topright", c("Plus", "Minus"), 
#        col =c("black","red"), lty=1)
# 
