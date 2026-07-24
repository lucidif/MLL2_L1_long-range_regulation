#1)
#

library(ggplot2)
library(RColorBrewer)

#=============parametes

options(scipen = 999)
te_group_colors <- c(
  "K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"   = "#5b92cb",
  "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed"    = "#50c7ef",
  "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed" = "#ed2c85",
  "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed"  = "#ad2f93"
)


#=======functions

source("./git/downstream_multiomic/bin/distance_functions.R")


find_distances<-function(class="LINE", infolder.mll2peaks, outfolder,repeat_mask, add_info = FALSE){

pattern <- paste0("\\b(", class, "|", class, "?)\\b")

repeat_mask_su <- repeat_mask[grep(pattern, repeat_mask$V12),]

#repeat_mask_su<-repeat_mask[grep("\\b(LINE|LINE?)\\b",repeat_mask$V12),]
#invece di prendere solo le LINE fai la stessa cosa con tutti gli altri timpi di elementi ripetuti e vedi cehe succede

#subset only L1 

refline_prop<-table(repeat_mask_su$V11)[order(table(repeat_mask_su$V11),decreasing = TRUE)]

df_refline_prop<-data.frame(refline_prop)
colnames(df_refline_prop)<-c("family","n")

#df_refline_prop$family<-as.character(df_refline_prop$family)

#df_refline_prop[13:nrow(df_refline_prop),1]<-"Others"

# sum(df_refline_prop[which(df_refline_prop$family=="Others"),2])

# df_refline_prop<-rbind(df_refline_prop[1:12,],
#                       data.frame(family="Others",
#                       n=sum(df_refline_prop[which(df_refline_prop$family=="Others"),2])
#                         )
#                       )

# df_refline_prop

bedrm<-repeat_mask[,c("V6","V7","V8","V11")]


#==========================================
#calc distribution of distances between  
#==========================================

atb<-c("K4me3")
distances_list<-NULL
for (i in 1:length(atb)){
  
    file_pattern="_distal_CpG_minus_Double_KO_vs_F_F.bed"  
  
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_su,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
                               )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],"_",class,"_",file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest ", class),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  file_pattern="_distal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_su,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],"_",class,"_",file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest ", class),
       sub=paste0("median=",median(densval)))
      abline(v=median(densval))
  dev.off()
  
  #=============================================================================
  
  file_pattern="_proximal_CpG_minus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_su,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],"_",class,"_",file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest ", class),
       sub=paste0("median=",median(densval)))
  abline(v=median(densval))
  dev.off()
  
  #==========================================================================
  
  file_pattern="_proximal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  distances<-extractdistances( peaks.coords = d,  
                               repeat_mask_line=repeat_mask_su,
                               addfamily = TRUE,
                               addpeaks = TRUE,
                               addlines = TRUE,
                               add5prime = TRUE
  )
  
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  densval<-as.numeric(as.character(distances[,1]))
  pdf(paste0(outfolder,atb[i],"_",class,"_",file_pattern,".pdf"))
  plot(density(densval),
       main=paste0(atb[i],file_pattern," nearest ", class),
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

##add coordinates 

for (i in 1:length(distances_list)){
  pre<-data.frame(distance=as.numeric(as.character(get(distances_list[i])[,1])),
                  group=distances_list[i], 
                  class= rep(class, nrow(get(distances_list[i]))) 
                  )

  pre_add_info<- data.frame(peak_chr=as.character(get(distances_list[i])[,2]),
                  peak_start=as.numeric(as.character(get(distances_list[i])[,3])), 
                  peak_end=as.numeric(as.character(get(distances_list[i])[,4])),
                  line_chr=as.character(get(distances_list[i])[,6]),
                  line_start=as.numeric(as.character(get(distances_list[i])[,7])),
                  line_end=as.numeric(as.character(get(distances_list[i])[,8])),
                  line_strand=as.character(get(distances_list[i])[,9]),
                  family=as.character(get(distances_list[i])[,5])
                  )

  if(i ==1 ){
    final<-pre
    final_add_info<-pre_add_info
  }else{
    final<-rbind(final,pre)
    final_add_info<-rbind(final_add_info,pre_add_info)
  }
}

if(add_info==TRUE){
  df<-cbind(final, final_add_info)
}else{
  df<-final
}



return(df)

}


#========input files

refpath.fasta<-"in/references/fasta/UCSC_GRCm38/mm10.fa"

peaks_root_folder="outs/Lara_multiomic_analysis/outs/CHiP_postprocessing_line1_dist"

antimml2.annotated.peaksFile="outs/Lara_multiomic_analysis/in/test_chipseq_downstream/macs_broadpeaks/Anti-GFP.mLb.mkD.sorted_peaks.annotatePeaks.txt"

antimml2.peaksfile="outs/Lara_multiomic_analysis/in/test_chipseq_downstream/bedtools_window/coordinate.bed"


#==================parameters
infolder.diffpeaks="outs/Lara_multiomic_analysis/in/test_chipseq_downstream/diffbind/"
infolder.mll2peaks="outs/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/"
outfolder="outs/Lara_multiomic_analysis/outs/"
infolder.ucsc="in/ucsc/"

#==================number of peaks by subclassification
antimml2.annotated.peaks<-read.table(antimml2.annotated.peaksFile, sep="\t", header=TRUE)

peaks.strands<-cbind(id=paste0(antimml2.annotated.peaks$Chr,"_",
                               antimml2.annotated.peaks$Start,"_",
                               antimml2.annotated.peaks$End
),antimml2.annotated.peaks$Strand)

antimml2.peaks<-read.table(antimml2.peaksfile,sep="\t", header=FALSE)
total.peaks<-nrow(antimml2.peaks)

#==========================================
# extract line annotations from repeat mask
#==========================================

repeat_mask<-read.table("outs/Lara_multiomic_analysis/in/ucsc/rmsk.txt",sep="\t",header=FALSE)
#repeat_mask<-read.table("in/ucsc/rmsk.txt",sep="\t",header=FALSE)
nrow(repeat_mask)

print(unique(repeat_mask$V12))

#  [1] "LINE"           "Simple_repeat"  "LTR"            "SINE"           "Low_complexity"
#  [6] "DNA"            "snRNA"          "Other"          "Satellite"      "Unknown"       
# [11] "SINE?"          "srpRNA"         "tRNA"           "LTR?"           "RNA"           
# [16] "scRNA"          "RC"             "rRNA"           "DNA?"           "RC?"           
# [21] "LINE?"          "Retroposon"

print(unique(repeat_mask$V12))

theme_classic <- function(base_size = 11, base_family = "") {
  ggplot2::theme_classic(base_size = 7, base_family = base_family) +
    theme(
      legend.text  = element_text(size = 5),
      legend.title = element_blank(),
      axis.text    = element_text(size = 6),
      axis.title   = element_text(size = 7),
      plot.title   = element_text(size = 7)
    )
}

df<-find_distances(class="LINE",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask, add_info=TRUE )

#df<-find_distances(class="LINE", infolder.mll2peaks="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/", outfolder="outs/", repeat_mask = repeat_mask )


LINE_plot<-ggplot(df, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="LINE Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))

ggsave(paste0(outfolder,"LINE_cumulative_distribution.pdf"), plot = LINE_plot, width = 14, height = 10, units="cm", device = "pdf")


df_SINE<-find_distances(class="SINE",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask, add_info=TRUE)

SINE_plot<-ggplot(df_SINE, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="SINE Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))

ggsave(paste0(outfolder,"SINE_cumulative_distribution.pdf"), plot = SINE_plot, width = 14, height = 10, units="cm", device = "pdf")


df_LTR<-find_distances(class="LTR",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask, add_info=TRUE)

LTR_plot<-ggplot(df_LTR, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="LTR Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))


ggsave(paste0(outfolder,"LTR_cumulative_distribution.pdf"), plot = LTR_plot, width = 14, height = 10, units="cm", device = "pdf")


df_RC<-find_distances(class="RC",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask, add_info=TRUE)

RC_plot<-ggplot(df_RC, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="RC Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic()  +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))


ggsave(paste0(outfolder,"RC_cumulative_distribution.pdf"), plot = RC_plot, width = 14, height = 10, units="cm", device = "pdf")


df_DNA<-find_distances(class="DNA",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask, add_info=TRUE)

DNA_plot<-ggplot(df_DNA, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="DNA Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))


ggsave(paste0(outfolder,"DNA_cumulative_distribution.pdf"), plot = DNA_plot, width = 14, height = 10, units="cm", device = "pdf")


#df_retroposon<-find_distances(class="Retroposon",
#infolder.mll2peaks=infolder.mll2peaks,
#outfolder=outfolder)

#ggplot(df_retroposon, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
#  labs(title="D4 Distance Cumulative Distribution",
#       y = "perc", x="distance")+
#  theme_classic()

#ERV only

#repeat_mask

class<-"LTR"
pattern <- paste0("\\b(", class, "|", class, "?)\\b")

repeat_mask_ltr <- repeat_mask[grep(pattern, repeat_mask$V12),]

repeat_mask_erv<-repeat_mask_ltr[grepl("ERV", repeat_mask_ltr$V13),]

df_ERV<-find_distances(class="LTR",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = repeat_mask_erv, add_info=TRUE)

ERV_plot<-ggplot(df_ERV, aes(as.numeric(as.character(distance)), colour = group)) + stat_ecdf(geom = "step")+
  labs(title="ERV Distance Cumulative Distribution",
       y = "perc", x="distance")+
  theme_classic() +
  scale_color_manual(values = te_group_colors) +
  theme(axis.line  = element_line(color = "black"),
        axis.text  = element_text(color = "black"),
        axis.title = element_text(color = "black"))


ggsave(paste0(outfolder,"ERV_cumulative_distribution.pdf"), plot = ERV_plot, width = 14, height = 10, units="cm", device = "pdf")

#merged cumulatice distribution

 LINE_dcm<-df[which(df$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),] 
 SINE_dcm<-df_SINE[which(df_SINE$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]
 
 LTR_dcm<-df_LTR[which(df_LTR$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]
 ERV_dcm<-df_ERV[which(df_ERV$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]
 head(ERV_dcm)
 ERV_dcm$class<-"ERV (LTR subset)"

 RC_dcm<-df_RC[which(df_RC$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]
 DNA_dcm<-df_DNA[which(df_DNA$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]

TE_dcm<-rbind(
  LINE_dcm,
  SINE_dcm, 
LTR_dcm, 
ERV_dcm#, 
#RC_dcm, 
#DNA_dcm
)

te_colors <- c(
  "LTR"              = "#ed4b8e",
  "SINE"             = "#5d93cc",
  "LINE"             = "#50caf1",
  "ERV (LTR subset)" = "#b2469a"
)

# TE_plot <- ggplot(TE_dcm, aes(x = as.numeric(as.character(distance)), colour = class)) + 
#   stat_ecdf(geom = "step", linewidth = 1) + # Aggiunto linewidth = 1 per rendere la linea un po' più visibile (opzionale)
#   scale_color_manual(values = te_colors) +  # <--- Ecco il comando magico
#   labs(
#     title = "TE Distance Cumulative Distribution from distal CpG minus",
#     y = "perc", 
#     x = "distance"
#   ) +
#   theme_classic() +
#   theme(
#     axis.line  = element_line(color = "black"),
#     axis.text  = element_text(color = "black"),
#     axis.title = element_text(color = "black")
#   )

TE_plot <- ggplot(TE_dcm, aes(x = as.numeric(as.character(distance)), colour = class)) + 
  stat_ecdf(geom = "step", linewidth = 0.4) +
  scale_color_manual(values = te_colors) +
  labs(
    title = "TE - Distance cumulative distribution",
    y = "Fraction of peaks", 
    x = "Distance from nearest L1 (bp)"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0, size = 11),
    axis.line     = element_line(color = "black"),
    axis.text     = element_text(color = "black"),
    axis.title    = element_text(color = "black"),
    legend.title  = element_blank(),
    legend.position = "right"
  )

# TE_plot<-ggplot(TE_dcm, aes(as.numeric(as.character(distance)), colour = class)) + stat_ecdf(geom = "step")+
#   labs(title="TE Distance Cumulative Distribution from distal CpG minus",
#        y = "perc", x="distance")+
#   theme_classic() +
#   theme(axis.line  = element_line(color = "black"),
#         axis.text  = element_text(color = "black"),
#         axis.title = element_text(color = "black"))


ggsave(paste0(outfolder,"TE_cumulative_distribution.pdf"), plot = TE_plot, width = 14, height = 10, units="cm", device = "pdf")

#===================================================================
#extimate pvalue with wilcoxon test for each pair of TE classes
#===================================================================

compair.curve.TE<-function(a,b,kmatrix,data,group_col="class",value_col="distance", test=c("Wilcox","KS")){

  va<-as.numeric(as.character(data[data[[group_col]] == a, value_col]))
  vb<-as.numeric(as.character(data[data[[group_col]] == b, value_col]))

  if(test[1]=="Wilcox"){
    a_VS_b<-wilcox.test(va, vb)
  }

  if(test[1]=="KS"){
    a_VS_b<-ks.test(va, vb)
  }

  kmatrix[a,b]<-a_VS_b$p.value
  kmatrix[b,a]<-a_VS_b$p.value

  return(kmatrix)

}

TE_classes<-unique(TE_dcm$class)

kmatrix.TE<-matrix(nrow=length(TE_classes), ncol=length(TE_classes))
colnames(kmatrix.TE)<-TE_classes
rownames(kmatrix.TE)<-TE_classes

for(i in 1:length(TE_classes)){
  for(j in i:length(TE_classes)){
    a<-TE_classes[i]
    b<-TE_classes[j]
    kmatrix.TE<-compair.curve.TE(a,b,kmatrix.TE,data=TE_dcm,group_col="class",value_col="distance",test="Wilcox")
  }
}

kmatrix.TE

#outfolder /media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs

write.table(kmatrix.TE, file=paste0(outfolder,"TE_cumulative_distribution_pval_wilcox.tsv"),
            col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")


##piechart

re_labeled_rm <- repeat_mask 

unique(re_labeled_rm$V12)
#  [1] "LINE"           "Simple_repeat"  "LTR"            "SINE"           "Low_complexity"
#  [6] "DNA"            "snRNA"          "Other"          "Satellite"      "Unknown"       
# [11] "SINE?"          "srpRNA"         "tRNA"           "LTR?"           "RNA"           
# [16] "scRNA"          "RC"             "rRNA"           "DNA?"           "RC?"           
# [21] "LINE?"          "Retroposon"


library(dplyr)
re_labeled_rm <- re_labeled_rm %>%
  filter(V12 %in% c("LINE", "SINE", "LTR"))

unique(re_labeled_rm$V12)

re_labeled_rm$V11<-paste(re_labeled_rm$V12,re_labeled_rm$V13,sep="_")
re_labeled_rm$V12<-"ALL"

all_TE_distances<-find_distances(class="ALL",
infolder.mll2peaks=infolder.mll2peaks,
outfolder=outfolder, repeat_mask = re_labeled_rm, add_info=TRUE)

TE_dcm_1kb<-all_TE_distances[which(all_TE_distances[,"distance"]<=1000 & all_TE_distances$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]

TE_dcm_500b<-all_TE_distances[which(all_TE_distances[,"distance"]<=500 & all_TE_distances$group=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"),]

#TE_dcm_1kb$family <- sub("_.*", "", TE_dcm_1kb$family)

#====================================================================
#========================================================================================

TE_dcm_1kb$class <- ifelse(
  grepl("^LTR_ERV", TE_dcm_1kb$family),   # family è tipo "LTR_ERV1", "LTR_ERVK"
  "ERV",
  sub("_.*", "", TE_dcm_1kb$family)        # LINE, SINE, LTR (non-ERV)
)

TE_dcm_500b$class <- ifelse(
  grepl("^LTR_ERV", TE_dcm_500b$family),   # family è tipo "LTR_ERV1", "LTR_ERVK"
  "ERV",
  sub("_.*", "", TE_dcm_500b$family)        # LINE, SINE, LTR (non-ERV)
)

#=====================================
TE_dcm_1kb <- TE_dcm_1kb %>%
  mutate(macro_class = case_when(
    grepl("ERV", family) ~ "ERV",   # Se contiene ERV, classificalo come ERV (ha la precedenza)
    grepl("LTR", family) ~ "LTR",   # Se contiene LTR, classificalo come LTR
    grepl("LINE", family) ~ "LINE", # Se contiene LINE, classificalo come LINE
    grepl("SINE", family) ~ "SINE", # Se contiene SINE, classificalo come SINE
    TRUE ~ "OTHER"                  # Per eventuali altri elementi
  ))

TE_dcm_500b <- TE_dcm_500b %>%
  mutate(macro_class = case_when(
    grepl("ERV", family) ~ "ERV",
    grepl("LTR", family) ~ "LTR",
    grepl("LINE", family) ~ "LINE",
    grepl("SINE", family) ~ "SINE",
    TRUE ~ "OTHER"
  ))
#==========================================

TE_dcm_1kb_piedf <- as.data.frame(table(TE_dcm_1kb$macro_class))
TE_dcm_500b_piedf <- as.data.frame(table(TE_dcm_500b$macro_class))

TE_dcm_1kb_piedf$pct <- round(TE_dcm_1kb_piedf$Freq / sum(TE_dcm_1kb_piedf$Freq) * 100, 1)
TE_dcm_500b_piedf$pct <- round(TE_dcm_500b_piedf$Freq / sum(TE_dcm_500b_piedf$Freq) * 100, 1)

TE_dcm_1kb_piedf$label <- paste0(TE_dcm_1kb_piedf$Var1, "\n", TE_dcm_1kb_piedf$pct, "%")
TE_dcm_500b_piedf$label <- paste0(TE_dcm_500b_piedf$Var1, "\n", TE_dcm_500b_piedf$pct, "%")

# poi usa $class invece di $family per il piechart
# TE_dcm_1kb_piedf <- as.data.frame(table(TE_dcm_1kb$class))
# TE_dcm_500b_piedf <- as.data.frame(table(TE_dcm_500b$class))

# TE_dcm_1kb_piedf$pct <- round(TE_dcm_1kb_piedf$Freq / sum(TE_dcm_1kb_piedf$Freq) * 100, 1)
# TE_dcm_500b_piedf$pct <- round(TE_dcm_500b_piedf$Freq / sum(TE_dcm_500b_piedf$Freq) * 100, 1)

# TE_dcm_1kb_piedf$label <- paste0(TE_dcm_1kb_piedf$Var1, "\n", TE_dcm_1kb_piedf$pct, "%")
# TE_dcm_500b_piedf$label <- paste0(TE_dcm_500b_piedf$Var1, "\n", TE_dcm_500b_piedf$pct, "%")

cols <- c(
  "LTR"  = "#ed4b8e",
  "ERV"  = "#b2469a",   # viola più saturo per distinguere da LTR
  "LINE" = "#50caf1",
  "SINE" = "#5d93cc"
)

ggplot(TE_dcm_1kb_piedf, aes(x = 2, y = Freq, fill = Var1)) +
  geom_col(width = 1, color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(10, 10, 10, 10),
    plot.title       = element_text(size = 21, face = "bold", hjust = 0.5)
  ) +
  geom_text(
    aes(label = paste0(pct, "%")),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  labs(fill = "TE class", title = "TE within 1kb of distance from dcm")

write.table()
ggsave(paste0(outfolder, "TE_cumulative_piechart.pdf"), plot = last_plot(), width = 6, height = 5)


ggplot(TE_dcm_500b_piedf, aes(x = 2, y = Freq, fill = Var1)) +
  geom_col(width = 1, color = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  xlim(0.5, 3.5) + # <--- aumentato ulteriormente per lasciare spazio alle etichette
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(10, 10, 10, 10),
    plot.title       = element_text(size = 21, face = "bold", hjust = 0.5)
  ) +
  geom_text(
    aes(x = 3.1, label = paste0(pct, "%")), # <--- spostato più all'esterno (era 2.7)
    position = position_stack(vjust = 0.5),
    size = 3.5, 
    color = "black", 
    fontface = "bold"
  ) +
  labs(fill = "TE class", title = "TE within 500 b of distance from dcm")
  
ggsave(paste0(outfolder, "500b_TE_cumulative_piechart.pdf"), plot = last_plot(), width = 6, height = 5)

unique(TE_dcm_1kb_piedf$Var1)
