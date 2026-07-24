#==============================
# Environment
#==============================

#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

library(ggplot2)

allpreaks=FALSE

#=======================================
#  functions
#=======================================

wkfolder=paste0(getwd(),"/outs/Lara_multiomic_analysis")

source("./git/downstream_multiomic/bin/hypergeometric.R")
source("./git/downstream_multiomic/bin/distance_functions.R")
source("git/bioinfoGenerals/base/bin/fun_compare_curves.R")

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

degs_D4<-read.table(paste0(wkfolder,"/in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv"), header=TRUE, sep="\t")

degs_D0<-read.table(paste0(wkfolder,"/in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv"), header=TRUE, sep="\t")

degs_wtD4D0<-read.table(paste0(wkfolder,"/in/build38_DEseq2_RNAseq/D4WT_VS_D0WT.deseq2.results.tsv"), header=TRUE, sep="\t")

allanno.genes<-read.table(paste0(wkfolder,"/in/GREATv4/GREATv4.genes.mm10.tsv"), sep="\t", header=FALSE)

broad_narrow_distribution<-read.table(paste0(wkfolder,"/outs/broad_narrow_distribution.tsv"), sep="\t", header=TRUE)

#great.dcm<-read.table("outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcm/20241126-public-4.0.4-hpxvmj-mm10-all-gene.txt", sep="\t", header=FALSE)

#great.dcp<-read.table("outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcp/20241128-public-4.0.4-qlnior-mm10-all-gene.txt", sep="\t", header=FALSE)

#nfcore RNAseq annotation extracted from gtf
anno<-read.table(paste0(wkfolder,"/in/build38_DEseq2_RNAseq/mm10.anno.tsv"), header=TRUE, sep="\t")


#gtf of annotation used as input also for RNAseq analysis
#anno.gtf<-read.table("in/build38_DEseq2_RNAseq/mm10.refGene.gtf.gz", sep="\t")
#anno.gtf<-anno.gtf[which(anno.gtf$V3=="CDS"),]

#TSS annotated reference
refseq.anno<-read.table(paste0(wkfolder,"/in/ucsc/build38_mm10_ncbiRefSeqCurated.txt.gz"), sep="\t")
refseq.anno[,"V5"]<-refseq.anno[,"V5"]+1
refseq.anno<-refseq.anno[,c("V13","V3","V5","V6","V4")]

#REMOVE the genes strand "+" that have same gene,start and the genes strand "-" gene, end
refseq.anno.stplus<-refseq.anno[which(refseq.anno$V4=="+"),]
refseq.anno.stplus<-refseq.anno.stplus[which(duplicated(refseq.anno.stplus[,c("V13","V3", "V5")])==FALSE),]

refseq.anno.stmin<-refseq.anno[which(refseq.anno$V4=="-"),]
refseq.anno.stmin<-refseq.anno.stmin[which(duplicated(refseq.anno.stmin[,c("V13","V3", "V6")])==FALSE),]

refseq.anno<-rbind(refseq.anno.stplus, refseq.anno.stmin)


#==============EXTRACT BROAD AND NARROW PEAKS ASSOCIATED TO DOWN DEGs

set1 <- broad_narrow_distribution[which(broad_narrow_distribution$V2 == "proximal_dko_wt_loose"), "intg"]
set2 <- broad_narrow_distribution[which(broad_narrow_distribution$V2 == "down"), "intg"]

length(setdiff(set1,set2))

nodown<-broad_narrow_distribution[which(broad_narrow_distribution$intg %in% setdiff(set1,set2)),]

nrow(nodown)

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

# write.table(annobed, "in/ucsc/build38_mm10_ncbiRefSeqCurated.bed",
#             sep="\t",
#             col.names = FALSE,
#             row.names = FALSE,
#             quote= FALSE
# )


unique(gene.annWithPosition$Chromosome)


#===========================================
# peaks
#===========================================


peaks.near.ref.file<-"in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
dcp_peaks.near.ref.file<-"in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed"

dcm_peak.near.ref<-read.table(peaks.near.ref.file,
                                sep="\t",
                                header=FALSE
)

dcp_peaks.near.ref<-read.table(dcp_peaks.near.ref.file,
                                 sep="\t",
                                 header=FALSE
)

dcp_peaks.near.ref<-data.frame(dcp_peaks.near.ref$V1,
           dcp_peaks.near.ref$V2,
           dcp_peaks.near.ref$V3,
           paste0("dcp_peak_",row.names(dcp_peaks.near.ref)),
           rep(NA,nrow(dcp_peaks.near.ref)),
           rep(".",nrow(dcp_peaks.near.ref))
           )

colnames(dcp_peaks.near.ref)<-c("V1","V2","V3","V4","V5","V6")

if (allpreaks==TRUE){


  pcm_peaks.near.ref<-read.table(pcm_peaks.near.ref.file,
                                   sep="\t",
                                   header=FALSE
  )

  pcp_peaks.near.ref<-read.table(pcp_peaks.near.ref.file,
                                   sep="\t",
                                   header=FALSE
  )
}


dcp_near_k27ac.file<-"outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed"

dcp_near_k27ac<-read.table(dcp_near_k27ac.file,
                                sep="\t",
                                header=FALSE
)

dcp_near_k27ac<-cbind(dcp_near_k27ac,
                  paste0("dcp_near_k27ac_",row.names(dcp_near_k27ac)),
                  rep(NA,nrow(dcp_near_k27ac)),
                  rep(".",nrow(dcp_near_k27ac))
                )

colnames(dcp_near_k27ac)<-c("V1","V2","V3","V4","V5","V6")



dcm_near_k27ac.file<-"outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcm_lossk4me3_near_D4_k27ac.bed"

dcm_near_k27ac<-read.table(dcm_near_k27ac.file,
                                sep="\t",
                                header=FALSE
)

dcm_near_k27ac<-cbind(dcm_near_k27ac,
                  paste0("dcm_near_k27ac_",row.names(dcm_near_k27ac)),
                  rep(NA,nrow(dcm_near_k27ac)),
                  rep(".",nrow(dcm_near_k27ac))
                )

colnames(dcm_near_k27ac)<-c("V1","V2","V3","V4","V5","V6")


#==============line peaks




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

#dko.dis.CpG.minus.peak<-peaks.near.ref
#dko.dis.CpG.plus.peak
#dko.pro.CpG.plus.peak
#dko.pro.CpG.minus.peak

#dcm_peak.near.ref

if (allpreaks==TRUE){

  #dcp_peaks.near.ref
  pcm_peaks.near.ref
  pcp_peaks.near.ref

  D0_down_pro<-distance_plot(minus_peaks=pcm_peaks.near.ref,
                             plus_peaks = pcp_peaks.near.ref  ,
                             degs = degs_D0_ann_diff_down,
                             title="proximal peaks near line D0 down",
                             window = 5000,
                             window.type = "inside"
  )

  #exportBed(D0_down_pro, file="./outs/D0_down_pcp_L1near.bed")

  # D0_up_pro<-distance_plot(minus_peaks=pcm_peaks.near.ref,
  #                          plus_peaks = pcp_peaks.near.ref,
  #                          degs = degs_D0_ann_diff_up,
  #                          title="proximal peaks near line D0 up",
  #                          window = 5000,
  #                          window.type = "inside"
  # )

}


D0_down_dis_win<-distance_plot(minus_peaks=dcm_peak.near.ref,
                               plus_peaks = dcp_peaks.near.ref ,
                               degs=degs_D0_ann_diff_down,
                               title="distal peaks near line D0 down",
                               window = 1000000,
                               window.type = "inside"
)

D0_nodiff_dis_win<-distance_plot(minus_peaks=dcm_peak.near.ref,
                                 plus_peaks = dcp_peaks.near.ref ,
                                 degs=nodiff_D0,
                                 title="distal peaks near line D0 all annotated genes",
                                 window = 1000000,
                                 window.type = "inside"
)

D0_up_dis_win<-distance_plot(minus_peaks=dcm_peak.near.ref,
                             plus_peaks = dcp_peaks.near.ref ,
                             degs=degs_D0_ann_diff_up,
                             title="distal peaks near line D0 up",
                             window = 1000000,
                             window.type = "inside"
)

df_windowed<-as.data.frame(
  rbind(

    cbind(distances=D0_down_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_DOWN")
    ),

    cbind(distances=D0_up_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_UP")
    ),

    cbind(distances=D0_nodiff_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_NODIFF")
    )
  )
)


df_windowed$distances<-as.numeric(df_windowed$distances)


#cumulative distribution plot

ggplot(df_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D0 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
    scale_color_manual(values = c("dis.CpG.plus_DOWN" = "#3f3092",
                                "dis.CpG.plus_UP" ="#608dd1",
                                "dis.CpG.plus_NODIFF" = "#bdbbba"
                                ))


ggsave("outs/distribution_distances_downDEGs_dcp/dis_downDEGs_dcp.png", dpi = 300, width = 8, height = 6, units = "in")

#grps<-c("dis.CpG.plus_DOWN",   "dis.CpG.plus_UP", "dis.CpG.plus_NODIFF")


wilcox.matrix(df_windowed = df_windowed, file="outs/distribution_distances_downDEGs_dcp/dis_downDEGs_dcp_pval_wilcox.tsv")



#====================================================================================================================
#====================================================================================================================
#
#                 D4   
#
#====================================================================================================================
#====================================================================================================================



#===========================================
# distal cpg plus (with the addiction of induced down in h27ac)
#===========================================


D4_down_dis_win<-distance_plot(minus_peaks = dcm_peak.near.ref,
                               plus_peaks = dcp_peaks.near.ref ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)


D4_nodiff_dis_win<-distance_plot(minus_peaks=dcm_peak.near.ref,
                                 plus_peaks = dcp_peaks.near.ref ,
                                 degs=nodiff_D4,
                                 title="distal peaks near line D4 all annotated genes",
                                 window = 1000000,
                                 window.type = "inside"
)

D4_up_dis_win<-distance_plot(minus_peaks=dcm_peak.near.ref,
                             plus_peaks = dcp_peaks.near.ref ,
                             degs=degs_D4_ann_diff_up,
                             title="distal peaks near line D4 up",
                             window = 1000000,
                             window.type = "inside"
)

#==========================add curve of subsetted peaks by h3k27ac

overlap_downdko_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, degs_D4_ann_diff_down$ensembl)


#overlap_down.torem <- intersect(degs_wtD4D0_ann_diff$ensembl, degs_D4_ann_diff_down$ensembl)

#D4_no_wtD4_D0_down <- degs_D4_ann_diff_down[which(!degs_D4_ann_diff_down$ensembl %in% overlap_down.torem),]
#D4_in_wtD4_D0_down <- degs_D4_ann_diff_down[which( degs_D4_ann_diff_down$ensembl %in% overlap_down.torem ),]

D4_downdko_upwt <- degs_D4_ann_diff_down[which( degs_D4_ann_diff_down$ensembl %in% overlap_downdko_upwt ),]

#overlap_downdko_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, degs_D4_ann_diff_down$ensembl)
#D4_downdko_upwt <- degs_D4_ann_diff_down[which( degs_D4_ann_diff_down$ensembl %in% overlap_downdko_upwt ),]

D4k27ac_downDKO_wtD4D0up<-distance_plot(minus_peaks=dcm_near_k27ac,
                                     plus_peaks = dcp_near_k27ac ,
                                     degs=D4_downdko_upwt,
                                     title="distal peaks near line",
                                     window = 1000000,
                                     window.type = "inside"
)

df_d4_plus_windowed<-as.data.frame(
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

    cbind(distances=D4k27ac_downDKO_wtD4D0up$dis.CpG.plus.degs$distances,
          group=rep("induced.dis.CpG.plus_DOWN")
    )

  )
)


df_d4_plus_windowed$distances<-as.numeric(df_d4_plus_windowed$distances)

ggp1<-ggplot(df_d4_plus_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.plus_DOWN" = "#3f3092",
                                "dis.CpG.plus_UP" ="#608dd1",
                                "dis.CpG.plus_NODIFF" = "#bdbbba",
                                "induced.dis.CpG.plus_DOWN" = "#ee2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  ))                                       

ggsave("outs/distribution_distances_downDEGs_dcp/D4_distribution_distances_downDEGs_dcp.png", ggp1,
       width = 8, height = 6, dpi = 300, units = "in")

wilcox.matrix(df_windowed = df_d4_plus_windowed, file="outs/distribution_distances_downDEGs_dcp/D4_dis_downDEGs_dcp_pval_wilcox.tsv")


#===========================================
#===========================================
# distal cpg plus only subset down in h27ac
#===========================================
#===========================================

D4k27ac_down_dis_win<-distance_plot(minus_peaks = dcm_near_k27ac ,
                               plus_peaks = dcp_near_k27ac ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

D4k27ac_nodiff_dis_win<-distance_plot(minus_peaks=dcm_near_k27ac,
                                 plus_peaks = dcp_near_k27ac ,
                                 degs=nodiff_D4,
                                 title="distal peaks near line D4 all annotated genes",
                                 window = 1000000,
                                 window.type = "inside"
)

D4k27ac_up_dis_win<-distance_plot(minus_peaks=dcm_near_k27ac,
                             plus_peaks = dcp_near_k27ac ,
                             degs=degs_D4_ann_diff_up,
                             title="distal peaks near line D4 up",
                             window = 1000000,
                             window.type = "inside"
)

df_d4_near_k27ac_windowed<-as.data.frame(
  rbind(

    cbind(distances=D4k27ac_down_dis_win$dis.CpG.plus.degs$distances ,
          group=rep("dis.CpG.plus_lossk4me3_near_k27ac_DOWN")
    ),
    
    cbind(distances=D4k27ac_up_dis_win$dis.CpG.plus.degs$distances ,
          group=rep("dis.CpG.plus_lossk4me3_near_k27ac_UP")
    ),

    cbind(distances=D4k27ac_nodiff_dis_win$dis.CpG.plus.degs$distances,
          group=rep("dis.CpG.plus_lossk4me3_near_k27ac_NODIFF")
    ),

    cbind(distances=D4k27ac_downDKO_wtD4D0up$dis.CpG.plus.degs$distances ,
          group=rep("dis.CpG.plus_lossk4me3_near_k27ac_induced_genes_DOWN")
    )

  )
)

df_d4_near_k27ac_windowed$distances<-as.numeric(df_d4_near_k27ac_windowed$distances)

ggp1<-ggplot(df_d4_near_k27ac_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.plus_lossk4me3_near_k27ac_DOWN" = "#3f3092",
                                "dis.CpG.plus_lossk4me3_near_k27ac_UP" ="#608dd1",
                                "dis.CpG.plus_lossk4me3_near_k27ac_NODIFF" = "#bdbbba",
                                "dis.CpG.plus_lossk4me3_near_k27ac_induced_genes_DOWN" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"



  )) 

ggsave("outs/distribution_distances_downDEGs_dcp/D4_distri_distan_downDEGs_dcp_nearK27ac.png", ggp1,
       width = 8, height = 6, dpi = 300, units = "in")


wilcox.matrix(df_windowed = df_d4_near_k27ac_windowed, file="outs/distribution_distances_downDEGs_dcp/D4_distri_distan_downDEGs_dcp_nearK27ac_pval_wilcox.tsv")


#===========================================
#===========================================
# distal cpg minus 
#===========================================
#===========================================



df_d4_minus_windowed<-as.data.frame(
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


df_d4_minus_windowed$distances<-as.numeric(df_d4_minus_windowed$distances)


ggp1<-ggplot(df_d4_minus_windowed, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("dis.CpG.minus_DOWN" = "#3f3092",
                                "dis.CpG.minus_UP" ="#608dd1",
                                "dis.CpG.minus_NODIFF" = "#bdbbba"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  )) 

ggsave("outs/distribution_distances_downDEGS_dcm/D4_distribution_distances_downDEGs_dcm.png", ggp1,
       width = 8, height = 6, dpi = 300, units = "in")

#==========================================
#==========================================
#distances from line elements
#==========================================
#==========================================


#================================
#genes distance by nearest dcp line
#================================

#dcm.l1<-read.table("outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed",sep="\t",header=FALSE)
#dcp.l1<-read.table("./outs/CHiP_postprocessing_line1_dist/dcp.l1.bed",sep="\t",header=FALSE)

l1_coords<-read.table("outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed",sep="\t",header=FALSE)


#head(dcp.l1)


l1near_down<-distance_plot(minus_peaks = l1_coords ,
                               plus_peaks = l1_coords ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_up<-distance_plot(minus_peaks = l1_coords ,
                               plus_peaks = l1_coords ,
                               degs=degs_D4_ann_diff_up,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)


l1near_undiff<-distance_plot(minus_peaks = l1_coords ,
                               plus_peaks = l1_coords ,
                               degs=nodiff_D4,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_downDKO_wtD4D0up<-distance_plot(minus_peaks= l1_coords,
                                     plus_peaks = l1_coords ,
                                     degs=D4_downdko_upwt,
                                     title="distal peaks near line",
                                     window = 1000000,
                                     window.type = "inside"
)

l1near_df<-as.data.frame(
  rbind(

    cbind(distances=l1near_down$dis.CpG.minus.degs$distances,
          group=rep("l1near_down")
    ),

    cbind(distances=l1near_up$dis.CpG.minus.degs$distances,
          group=rep("l1near_up")
    ),

    cbind(distances=l1near_undiff$dis.CpG.minus.degs$distances,
          group=rep("l1near_undiff")
    ),

    cbind(distances=l1near_downDKO_wtD4D0up$dis.CpG.minus.degs$distances,
          group=rep("l1near_induced_down")
    )

  )
)


l1near_df$distances<-as.numeric(l1near_df$distances)


ggp1<-ggplot(l1near_df, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("l1near_down" = "#3f3092",
                                "l1near_up" ="#608dd1",
                                "l1near_undiff" = "#bdbbba",
                                "l1near_induced_down" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  )) 

ggsave("outs/distribution_distances_downDEGS_dcm/D4_distribution_distances_downDEGs_dcm_l1_peaks.png", ggp1,
       width = 8, height = 6, dpi = 300, units = "in")

wilcox.matrix(df_windowed = l1near_df, file="outs/distribution_distances_downDEGs_dcp/D4_distribution_distances_downDEGs_dcm_l1_peaks_pval_wilcox.tsv")


#===============high dist
l1_coords_high<-read.table("outs/l1_d4_heatmap/high_sort2_l1_d4mll2_plotHeatmap.bed",sep="\t",header=FALSE)

l1near_down<-distance_plot(minus_peaks = l1_coords_high ,
                               plus_peaks = l1_coords_high ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_up<-distance_plot(minus_peaks = l1_coords_high ,
                               plus_peaks = l1_coords_high ,
                               degs=degs_D4_ann_diff_up,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)


l1near_undiff<-distance_plot(minus_peaks = l1_coords_high ,
                               plus_peaks = l1_coords_high ,
                               degs=nodiff_D4,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_downDKO_wtD4D0up<-distance_plot(minus_peaks= l1_coords_high,
                                     plus_peaks = l1_coords_high ,
                                     degs=D4_downdko_upwt,
                                     title="distal peaks near line",
                                     window = 1000000,
                                     window.type = "inside"
)

l1near_df<-as.data.frame(
  rbind(

    cbind(distances=l1near_down$dis.CpG.minus.degs$distances,
          group=rep("l1near_down")
    ),

    cbind(distances=l1near_up$dis.CpG.minus.degs$distances,
          group=rep("l1near_up")
    ),

    cbind(distances=l1near_undiff$dis.CpG.minus.degs$distances,
          group=rep("l1near_undiff")
    ),

    cbind(distances=l1near_downDKO_wtD4D0up$dis.CpG.minus.degs$distances,
          group=rep("l1near_induced_down")
    )

  )
)


l1near_df$distances<-as.numeric(l1near_df$distances)


ggp1_high<-ggplot(l1near_df, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("l1near_down" = "#3f3092",
                                "l1near_up" ="#608dd1",
                                "l1near_undiff" = "#bdbbba",
                                "l1near_induced_down" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  )) 

ggsave("outs/l1_d4_heatmap/high_D4_dist_downDEGs_dcm_l1_peaks.png", ggp1_high,
       width = 8, height = 6, dpi = 300, units = "in")


#=======low only dist

l1_coords_mid<-read.table("outs/l1_d4_heatmap/mid_sort2_l1_d4mll2_plotHeatmap.bed",sep="\t",header=FALSE)

l1near_down<-distance_plot(minus_peaks = l1_coords_mid ,
                               plus_peaks = l1_coords_mid ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_up<-distance_plot(minus_peaks = l1_coords_mid ,
                               plus_peaks = l1_coords_mid ,
                               degs=degs_D4_ann_diff_up,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)


l1near_undiff<-distance_plot(minus_peaks = l1_coords_mid ,
                               plus_peaks = l1_coords_mid ,
                               degs=nodiff_D4,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_downDKO_wtD4D0up<-distance_plot(minus_peaks= l1_coords_mid,
                                     plus_peaks = l1_coords_mid ,
                                     degs=D4_downdko_upwt,
                                     title="distal peaks near line",
                                     window = 1000000,
                                     window.type = "inside"
)

l1near_df<-as.data.frame(
  rbind(

    cbind(distances=l1near_down$dis.CpG.minus.degs$distances,
          group=rep("l1near_down")
    ),

    cbind(distances=l1near_up$dis.CpG.minus.degs$distances,
          group=rep("l1near_up")
    ),

    cbind(distances=l1near_undiff$dis.CpG.minus.degs$distances,
          group=rep("l1near_undiff")
    ),

    cbind(distances=l1near_downDKO_wtD4D0up$dis.CpG.minus.degs$distances,
          group=rep("l1near_induced_down")
    )

  )
)


l1near_df$distances<-as.numeric(l1near_df$distances)


ggp1_mid<-ggplot(l1near_df, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("l1near_down" = "#3f3092",
                                "l1near_up" ="#608dd1",
                                "l1near_undiff" = "#bdbbba",
                                "l1near_induced_down" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  )) 

ggsave("outs/l1_d4_heatmap/mid_D4_dist_downDEGs_dcm_l1_peaks.png", ggp1_mid,
       width = 8, height = 6, dpi = 300, units = "in")



#======= low only dist

l1_coords_low<-read.table("outs/l1_d4_heatmap/low_sort2_l1_d4mll2_plotHeatmap.bed",sep="\t",header=FALSE)

l1near_down<-distance_plot(minus_peaks = l1_coords_low ,
                               plus_peaks = l1_coords_low ,
                               degs=degs_D4_ann_diff_down,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_up<-distance_plot(minus_peaks = l1_coords_low ,
                               plus_peaks = l1_coords_low ,
                               degs=degs_D4_ann_diff_up,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)


l1near_undiff<-distance_plot(minus_peaks = l1_coords_low ,
                               plus_peaks = l1_coords_low ,
                               degs=nodiff_D4,
                               title="distal peaks near line D4 down",
                               window = 1000000,
                               window.type = "inside"
)

l1near_downDKO_wtD4D0up<-distance_plot(minus_peaks= l1_coords_low,
                                     plus_peaks = l1_coords_low ,
                                     degs=D4_downdko_upwt,
                                     title="distal peaks near line",
                                     window = 1000000,
                                     window.type = "inside"
)

l1near_df<-as.data.frame(
  rbind(

    cbind(distances=l1near_down$dis.CpG.minus.degs$distances,
          group=rep("l1near_down")
    ),

    cbind(distances=l1near_up$dis.CpG.minus.degs$distances,
          group=rep("l1near_up")
    ),

    cbind(distances=l1near_undiff$dis.CpG.minus.degs$distances,
          group=rep("l1near_undiff")
    ),

    cbind(distances=l1near_downDKO_wtD4D0up$dis.CpG.minus.degs$distances,
          group=rep("l1near_induced_down")
    )

  )
)


l1near_df$distances<-as.numeric(l1near_df$distances)


ggp1_low<-ggplot(l1near_df, aes(distances, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="D4 Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = c("l1near_down" = "#3f3092",
                                "l1near_up" ="#608dd1",
                                "l1near_undiff" = "#bdbbba",
                                "l1near_induced_down" = "#ed2c85"
                                #"induced.dis.CpG.plus_UP" = "green",
                                #"induced.dis.CpG.plus_NODIFF" ="yellow"
  )) 

ggsave("outs/l1_d4_heatmap/low_D4_dist_downDEGs_dcm_l1_peaks.png", ggp1_low,
       width = 8, height = 6, dpi = 300, units = "in")


