#1)
setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

library(ggplot2)
library(RColorBrewer)

#=======functions

source("./git/downstream_multiomic/bin/distance_functions.R")

#========input files

refpath.fasta<-"/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa"

antimml2.annotated.peaksFile="in/test_chipseq_downstream/macs_broadpeaks/Anti-GFP.mLb.mkD.sorted_peaks.annotatePeaks.txt"

antimml2.peaksfile="in/test_chipseq_downstream/bedtools_window/coordinate.bed"

#dko.pro.peaksFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_Double_KO_vs_F_F.bed"
#dko.dis.peaksFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_Double_KO_vs_F_F.bed"

#dko.pro.CpG.plus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_CpG_plus_Double_KO_vs_F_F.bed"
#dko.pro.CpG.minus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_CpG_minus_Double_KO_vs_F_F.bed"

#dko.dis.CpG.plus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_CpG_plus_Double_KO_vs_F_F.bed"
#dko.dis.CpG.minus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_CpG_minus_Double_KO_vs_F_F.bed"

#==================parameters
infolder.diffpeaks="in/test_chipseq_downstream/diffbind/"
infolder.mll2peaks="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/"
outfolder="outs/"
infolder.ucsc="in/ucsc/"

#==================number of peaks by subclassification
antimml2.annotated.peaks<-read.table(antimml2.annotated.peaksFile, sep="\t", header=TRUE)

peaks.strands<-cbind(id=paste0(antimml2.annotated.peaks$Chr,"_",
                               antimml2.annotated.peaks$Start,"_",
                               antimml2.annotated.peaks$End
),antimml2.annotated.peaks$Strand)

antimml2.peaks<-read.table(antimml2.peaksfile,sep="\t", header=FALSE)
total.peaks<-nrow(antimml2.peaks)

#dko.pro.CpG.plus.peak<-read.table(dko.pro.CpG.plus.peakFile, sep="\t", header=FALSE)
#tot.dko.pro.CpG.plus<-nrow(dko.pro.CpG.plus.peak)

#dko.pro.CpG.minus.peak<-read.table(dko.pro.CpG.minus.peakFile, sep="\t", header=FALSE)
#tot.dko.pro.CpG.minus<-nrow(dko.pro.CpG.minus.peak)

#dko.dis.CpG.plus.peak<-read.table(dko.dis.CpG.plus.peakFile, sep="\t", header=FALSE)
#tot.dko.dis.CpG.plus<-nrow(dko.dis.CpG.plus.peak)

#dko.dis.CpG.minus.peak<-read.table(dko.dis.CpG.minus.peakFile, sep="\t", header=FALSE)
#tot.dko.dis.CpG.minus<-nrow(dko.dis.CpG.minus.peak)

#==========================================
# extract line annotations from repeat mask
#==========================================
repeat_mask<-read.table("in/ucsc/rmsk.txt",sep="\t",header=FALSE)
nrow(repeat_mask)

print(unique(repeat_mask$V12))

#  [1] "LINE"           "Simple_repeat"  "LTR"            "SINE"           "Low_complexity"
#  [6] "DNA"            "snRNA"          "Other"          "Satellite"      "Unknown"       
# [11] "SINE?"          "srpRNA"         "tRNA"           "LTR?"           "RNA"           
# [16] "scRNA"          "RC"             "rRNA"           "DNA?"           "RC?"           
# [21] "LINE?"          "Retroposon"

repeat_mask_line<-repeat_mask[grep("\\b(LINE|LINE?)\\b",repeat_mask$V12),]

refline_prop<-table(repeat_mask_line$V11)[order(table(repeat_mask_line$V11),decreasing = TRUE)]

df_refline_prop<-data.frame(refline_prop)
colnames(df_refline_prop)<-c("family","n")

df_refline_prop$family<-as.character(df_refline_prop$family)

#top12+others

df_refline_prop[13:nrow(df_refline_prop),1]<-"Others"

sum(df_refline_prop[which(df_refline_prop$family=="Others"),2])

df_refline_prop<-rbind(df_refline_prop[1:12,],
                      data.frame(family="Others",
                      n=sum(df_refline_prop[which(df_refline_prop$family=="Others"),2])
                        )
                      )

df_refline_prop

ggplot(data=df_refline_prop, aes(x=family, y=as.numeric(n), fill=family)) +
  geom_bar(stat="identity")



bedrm<-repeat_mask[,c("V6","V7","V8","V11")]


write.table(bedrm ,file="./in/ucsc/l1_rmsk.bed", sep="\t", 
            col.names = FALSE, row.names = FALSE,
            quote=FALSE
            )



#==========================================
#calc distribution of distances between  
#==========================================


#atb<-c("K4me3","K4me2")

atb<-c("K4me3")
distances_list<-NULL
for (i in 1:length(atb)){
  
  # if(i==0){
  #   file_pattern="_proximal_Double_KO_vs_F_F.bed"
  # }else{
    file_pattern="_distal_CpG_minus_Double_KO_vs_F_F.bed"  
  # }
  
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_line,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
                               )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  file_pattern="_distal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_line,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  #=============================================================================
  
  file_pattern="_proximal_CpG_minus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_line,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  #==========================================================================
  
  file_pattern="_proximal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_line,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  #================================
  
  # file_pattern="_proximal_Double_KO_vs_F_F.bed"
  
  # d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  # distances<-extractdistances( peaks.coords = d,  
  #                              repeat_mask_line=repeat_mask_line,
  #                              addfamily = TRUE,
  #                              addpeaks = TRUE,
  #                              addlines = TRUE,
  #                              add5prime = TRUE
  # )
  
  #distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  
  
  
  
  
}

#cumulative distribution plot

for (i in 1:length(distances_list)){
  pre<-data.frame(distance=as.numeric(as.character(get(distances_list[i])[,1])),
                  group=distances_list[i])
  if(i ==1 ){
    final<-pre
  }else{
    final<-rbind(final,pre)
  }
}

df<-final

unique(df$group)

npeaks_dcm_under1kb<-length(which(df$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed" &df$distance <= 1000))
n_peaks_dcm_tot<-length(which(df$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"))

# Basic ECDF plot

df<-df[which(df$group!="K4me3_proximal_Double_KO_vs_F_F.bed"),]

cum_dist<-ggplot(df, aes(distance, colour = group)) + 
  stat_ecdf(geom = "step") +
  labs(title = "Cumulative Distance Distribution",
       y = "perc", x = "distance") +
  theme_classic() +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "#6a6a6a", size = 0.8) +
  annotate("text", x = 1000 + 2500, y = 1, label = npeaks_dcm_under1kb/n_peaks_dcm_tot, 
           vjust = -0.5, color = "#6a6a6a", size = 4) +
  scale_color_manual(values = c("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed" = "#5b92cb",
                                "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed" ="#50c7ef",
                                "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed" = "#ed2c85",
                                "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed" ="#ad2f93"
  )) +
  theme(axis.line = element_line(color = "black"),  
        axis.text = element_text(color = "black"),  
        axis.title = element_text(color = "black"))

ggsave("outs/CHiP_postprocessing_line1_dist/Cumulative_Density_plot.png", plot = cum_dist, width = 13, height = 8, dpi = 300)

# grps<-c("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed",   "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed", "K4me3_proximal_Double_KO_vs_F_F.bed" )

grps<-c("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed",   "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed")

df[which(df$group==grps[1]),"distance"]

compair.curve<-function(a,b,kmatrix,grps, test=c("Wilcox","KS")){

  # a_VS_b<-t.test( df[which(df$group==a),"distance"] ,  
  #                     df[which(df$group==b),"distance"], 
  #                     alternative = "two.sided", var.equal = FALSE)
  
  if(test[1]=="Wilcox"){
    a_VS_b<-wilcox.test(df$distance[df$group == a], df$distance[df$group == b])
  }

  if(test[1]=="KS"){
    a_VS_b<-ks.test(df$distance[df$group == a], df$distance[df$group == b])
  }
  
  kmatrix[a,b]<-a_VS_b$p.value
  kmatrix[b,a]<-a_VS_b$p.value
  
  return(kmatrix)
  
}


kmatrix<-matrix(nrow=length(grps), ncol=length(grps))
colnames(kmatrix)<-grps
rownames(kmatrix)<-grps

# dcm vs dcp
a<-grps[1]
b<-grps[2]
kmatrix<-compair.curve(a,b,kmatrix,grps)


# dcm vs pcm
a<-grps[1]
b<-grps[3]
kmatrix<-compair.curve(a,b,kmatrix,grps)

#dcm vs pcp
a<-grps[1]
b<-grps[4]
kmatrix<-compair.curve(a,b,kmatrix,grps)

#dcm vs p
# a<-grps[1]
# b<-grps[5]
# kmatrix<-compair.curve(a,b,kmatrix,grps)

#dcp vs pcm
a<-grps[2]
b<-grps[3]
kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[2]
b<-grps[4]
kmatrix<-compair.curve(a,b,kmatrix,grps)

# a<-grps[2]
# b<-grps[5]
# kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[3]
b<-grps[4]
kmatrix<-compair.curve(a,b,kmatrix,grps)

# a<-grps[3]
# b<-grps[5]
# kmatrix<-compair.curve(a,b,kmatrix,grps)

# a<-grps[4]
# b<-grps[5]
# kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[1]
b<-grps[1]
kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[2]
b<-grps[2]
kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[3]
b<-grps[3]
kmatrix<-compair.curve(a,b,kmatrix,grps)

a<-grps[4]
b<-grps[4]
kmatrix<-compair.curve(a,b,kmatrix,grps)

# a<-grps[5]
# b<-grps[5]
# kmatrix<-compair.curve(a,b,kmatrix,grps)


kr.matrix<-kmatrix
colnames(kr.matrix)<-c("dcm",   "dcp", "pcm", "pcp")
rownames(kr.matrix)<-c("dcm",   "dcp", "pcm", "pcp")

write.table(kr.matrix, file="outs/CHiP_postprocessing_line1_dist/Cumulative_Density_pval_wilcox.tsv", col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")

boxdf<-df


#====================================
#  line1 distribution around subgroups
#====================================

#1 pieDonut on repeat_mask.coordinates_extended

 repeat_mask_line<-repeat_mask_line[,c(6,7,8,10,11,12,13)]

 colnames(repeat_mask_line)<-c("chr","start","end","strand","family","element","type")

 repeat_mask_line_plus<-repeat_mask_line[which(repeat_mask_line$strand=="+"),]
 repeat_mask_line_minus<-repeat_mask_line[which(repeat_mask_line$strand=="-"),]

 repeat_mask.5prime<-rbind(cbind(repeat_mask_line_plus, fprime=repeat_mask_line_plus$start),
                           cbind(repeat_mask_line_minus, fprime=repeat_mask_line_minus$end))

 repeat_mask.coordinates_extended=cbind(repeat_mask.5prime$chr,
                                        repeat_mask.5prime$fprime,
                                        repeat_mask.5prime$family
 )

line.freqs<-as.data.frame(table(repeat_mask.coordinates_extended[,3]))
summary(line.freqs$Freq)

#remove under first quantile lines families (familyes not rapresented in annotation file)
su.line.freqs<-line.freqs[which(line.freqs$Freq>=as.numeric(summary(line.freqs$Freq)[2])),]

top.line.freqs<-line.freqs[order(line.freqs$Freq, decreasing = TRUE)[1:7],]

other.line.freqs<-line.freqs[order(line.freqs$Freq, decreasing = TRUE),]
levels(other.line.freqs[,"Var1"])<-c(levels(other.line.freqs[,"Var1"]),"Others")
other.line.freqs[8:nrow(other.line.freqs),"Var1"]<-as.factor("Others")

#su.line.freqs<-cbind(su.line.freqs,perc=su.line.freqs$Freq/sum(su.line.freqs$Freq))

#1 see rapresantated family in distal cpg minus group

dgm.line.freqs<-as.data.frame(table(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"family"]))
su.dgm.line.fraqs<-dgm.line.freqs[which(dgm.line.freqs$Freq>=
                                          as.numeric(summary(dgm.line.freqs$Freq)[2])),]
top.dgm.line.fraqs<-dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE)[1:7],]

other.dgm.line.freqs<-dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE),]
levels(other.dgm.line.freqs[,"Var1"])<-c(levels(other.dgm.line.freqs[,"Var1"]),"Others")
other.dgm.line.freqs[8:nrow(other.dgm.line.freqs),"Var1"]<-as.factor("Others")

#1 generate comparative plot

laneFreq.dataplot<-as.data.frame(rbind(cbind(su.line.freqs,group=rep("all")),
      cbind(su.dgm.line.fraqs,group=rep("only.dcm"))
      ))

#ggplot(laneFreq.dataplot, aes(fill=group, y=Freq, x=Var1)) + 
#  geom_bar(position="fill", stat="identity")

#ggplot(laneFreq.dataplot, aes(fill=Var1, y=Freq, x=group)) + 
#  geom_bar(position="fill", stat="identity")

top.laneFreq.dataplot<-as.data.frame(rbind(cbind(top.line.freqs,group=rep("all")),
                                           cbind(top.dgm.line.fraqs,group=rep("only.dcm"))
))

# ggplot(top.laneFreq.dataplot, aes(fill=Var1, y=Freq, x=group)) + 
#   geom_bar(position="fill", stat="identity")

other.laneFreq.dataplot<-as.data.frame(rbind(cbind(other.line.freqs,group=rep("all")),
                                             cbind(other.dgm.line.freqs,group=rep("only.dcm"))
))

length(unique(other.laneFreq.dataplot$Var1))

ggplot(other.laneFreq.dataplot, aes(fill=Var1, y=Freq, x=group)) + 
  geom_bar(position="fill", stat="identity")

#========================================
#====extract interesting line coordinates
#========================================

tarfam<-as.character(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group=="only.dcm" 
                              & other.laneFreq.dataplot$Var1 != "Others"),"Var1"])

tarfam_anno<-K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,c("line_chr","line_start","line_end","family","distances","line_strand")]
tarfam_anno[,"family"]<-paste0(tarfam_anno[,"family"],"_",rep(1:nrow(tarfam_anno)))

tarpeak<-K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,c("peaks_chr","peaks_start","peaks_end")]
tarpeak<-cbind(tarpeak,
              name=paste0("peak_",rep(1:nrow(tarpeak))),
              score=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"distances"],
              strand=rep(".",nrow(tarpeak))
              )
        


write.table(tarfam_anno, file="./outs/CHiP_postprocessing_line1_dist/dcm.l1.bed", 
            col.names = FALSE, 
            row.names = FALSE,
            quote= FALSE,
            sep="\t"
            )

write.table(tarpeak, file="./outs/CHiP_postprocessing_line1_dist/dcm.target.peaks.bed",
            col.names = FALSE, 
            row.names = FALSE,
            quote= FALSE,
            sep="\t"
            )

#================================
# extract other line1 coordinate
#================================

dcp_anno<-K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[,c("line_chr","line_start","line_end","family","distances","line_strand")]
dcp_anno[,"family"]<-paste0(dcp_anno[,"family"],"_",rep(1:nrow(dcp_anno)))

dcp_tarpeak<-K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[,c("peaks_chr","peaks_start","peaks_end")]
dcp_tarpeak<-cbind(dcp_tarpeak,
               name=paste0("peak_",rep(1:nrow(dcp_tarpeak))),
               score=K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[,"distances"],
               strand=rep(".",nrow(dcp_tarpeak))
)



write.table(dcp_anno, file="./outs/CHiP_postprocessing_line1_dist/dcp.l1.bed", 
            col.names = FALSE, 
            row.names = FALSE,
            quote= FALSE,
            sep="\t"
)

write.table(dcp_tarpeak, file="./outs/CHiP_postprocessing_line1_dist/dcp.target.peaks.bed",
            col.names = FALSE, 
            row.names = FALSE,
            quote= FALSE,
            sep="\t"
)


#filtered by distance peaks

filterbydistance(extractdistances.out=K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=2500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter_DKO_K4me3_dcp")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=2500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=2500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter_DKO_K4me3_pcm")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=2500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter_DKO_K4me3_pcp")

filterbydistance(extractdistances.out=K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=1000, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter1000_DKO_K4me3_dcp")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=1000, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter1000_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=1000, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter1000_DKO_K4me3_pcm")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=1000, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter1000_DKO_K4me3_pcp")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcp")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist=500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_pcp")

filterbydistance(extractdistances.out=K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=500, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_pcm")

filterbydistance(extractdistances.out=K4me3_proximal_Double_KO_vs_F_F.bed,
                 fildist=100000,
                 out="outs/CHiP_postprocessing_line1_dist/DKO_K4me3_proximal")

filterbydistance(extractdistances.out=K4me3_proximal_Double_KO_vs_F_F.bed,
                 fildist=1000,
                 out="outs/CHiP_postprocessing_line1_dist/disfilter1000_DKO_K4me3_proximal")

filterbydistance(extractdistances.out=K4me3_proximal_Double_KO_vs_F_F.bed,
                 fildist=1000,
                 out="outs/CHiP_postprocessing_line1_dist/disfilter500_DKO_K4me3_proximal")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=400, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter400_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=300, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter300_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=200, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter200_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=100, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter100_DKO_K4me3_dcm")

filterbydistance(extractdistances.out=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist=0, 
                 out="outs/CHiP_postprocessing_line1_dist/distfilter0_DKO_K4me3_dcm")


#distance plot after filtering
#distfilter500_DKO_K4me3_dcm.l1.bed

#center on 5prime
disf_dcm_complete<-read.table("outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm.l1.bed",
sep="\t",
header=FALSE
)
center_win=5000

for(i in 1:nrow(disf_dcm_complete)){
  if(disf_dcm_complete[i,6]=="+"){
    fprime<-disf_dcm_complete[i,2]
  }else{
    if(disf_dcm_complete[i,6]=="-"){
      fprime<-disf_dcm_complete[i,3]
    }
  }
  disf_dcm_complete[i,2]<-fprime-center_win
  disf_dcm_complete[i,3]<-fprime+center_win
}

write.table(disf_dcm_complete, file="outs/CHiP_postprocessing_line1_dist/recentered_distfil500_DKO_K4me3_dcm.l1.bed",
            sep="\t",
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE
)


distances_files<-c(
  "distfilter500_DKO_K4me3_dcm.l1.bed",
  "distfilter500_DKO_K4me3_dcp.l1.bed",
  "distfilter500_DKO_K4me3_pcm.l1.bed",
  "distfilter500_DKO_K4me3_pcp.l1.bed"
)

distances_list<-distances_list[which(distances_list!="K4me3_proximal_Double_KO_vs_F_F.bed")]

for (i in 1:length(distances_list)){
  
  distvar<-read.table(paste0("./outs/CHiP_postprocessing_line1_dist/",distances_files[i]),
              sep="\t",
              header =FALSE
              )
  

             
  pre<-data.frame(distance=as.numeric(as.character(distvar[,5])),
                  group=distances_files[i])
  if(i ==1 ){
    final<-pre
  }else{
    final<-rbind(final,pre)
  }
}

df<-final

ggplot(df, aes(distance, colour = group)) + stat_ecdf(geom = "step")+
  labs(title="Empirical Cumulative \n Density Function",
       y = "perc", x="distance")+
  theme_classic()

#======================================
#families by cutoff
#======================================

#test line families inside peaks bed by distance

files<-c("distfilter_DKO_K4me3_dcm.l1.bed",
         "distfilter1000_DKO_K4me3_dcm.l1.bed",
         "distfilter500_DKO_K4me3_dcm.l1.bed",
         "distfilter400_DKO_K4me3_dcm.l1.bed",
         "distfilter300_DKO_K4me3_dcm.l1.bed",
         "distfilter200_DKO_K4me3_dcm.l1.bed",
         "distfilter100_DKO_K4me3_dcm.l1.bed"
)

famout<-paste0("fm_rap_",files)

fam1val<-c()
fam2val<-c()
fam3val<-c()
fam4val<-c()
nm1<-c()
nm2<-c()
nm3<-c()
nm4<-c()

for (i in 1:length(files)){
  
  print("#=======================================")
  
  bd<-read.table(paste0("./outs/CHiP_postprocessing_line1_dist/",files[i]))
  
  arr<-bd[,4]
  
  transformed_array <- gsub("_\\d+$", "", arr)
  
  sror<-table(transformed_array)
  
  tb<-sror[order(sror,decreasing = TRUE)]
  
  assign(famout[i],tb)
  
  print(paste0(files[i]))
  print(tb[1:4])
  print(paste0("other lines",sum(tb[5:length(tb)])))
  
  print("ratio top4 vs others")
  print(sum(tb[1:4])/sum(tb[5:length(tb)]))
  
  fam1val[i]<-tb[1]
  fam2val[i]<-tb[2]
  fam3val[i]<-tb[3]
  fam4val[i]<-tb[4]
  nm1[i]<-names(tb[1])
  nm2[i]<-names(tb[2])
  nm3[i]<-names(tb[3])
  nm4[i]<-names(tb[4])
  
  if(i!=1){
    
    nmpre<-c(nm1[i-1],nm2[i-1],nm3[i-1],nm4[i-1])
    valpre<-c(fam1val[i-1],fam2val[i-1],fam3val[i-1],fam4val[i-1])
    
    preval1<-valpre[grep(nm1[i],nmpre)]
    preval2<-valpre[grep(nm2[i],nmpre)]
    preval3<-valpre[grep(nm3[i],nmpre)]
    preval4<-valpre[grep(nm4[i],nmpre)]
    
    print(paste0(nm1[i]," ratio with previous value ",fam1val[i]/preval1))
    print(paste0(nm2[i]," ratio with previous value ",fam2val[i]/preval2))
    print(paste0(nm3[i]," ratio with previous value ",fam3val[i]/preval3))
    print(paste0(nm4[i]," ratio with previous value ",fam4val[i]/preval4))
    
  }
  
  print("#=======================================")
  
  
}

for (i in 1:length(famout)){
    famfreq<-get(famout[i])  
    famfreq<-as.data.frame(cbind(famfreq,family=names(famfreq)))
    famfreq$famfreq<-as.numeric(famfreq$famfreq)
    famfreq[11:nrow(famfreq),2]<-"Others"
    
    ggplot(famfreq, aes(fill=family, y=famfreq, x=famout[i])) + 
      geom_bar(position="stack", stat="identity")
    
}




