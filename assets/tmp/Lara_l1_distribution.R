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

#subset only L1 

refline_prop<-table(repeat_mask_line$V11)[order(table(repeat_mask_line$V11),decreasing = TRUE)]

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

grps<-c("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed",   "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed", "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed")

#boxdf<-df


#====================================
#  line1 distribution around subgroups
#====================================

#1 pieDonut on repeat_mask.coordinates_extended

#subset only the l1 subfemilies

 bck_repeat_mask_line<-repeat_mask_line
 
 repeat_mask_line<-repeat_mask_line[which(repeat_mask_line$V13=="L1"),]
 unique(repeat_mask_line$V13)

 repeat_mask_line<-repeat_mask_line[,c(6,7,8,10,11,12,13)]
tt$V13
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

#all L1 families in the references
all.families<-line.freqs$Var1


# su.line.freqs<-line.freqs[which(line.freqs$Freq>=as.numeric(summary(line.freqs$Freq)[2])),]

# top.line.freqs<-line.freqs[order(line.freqs$Freq, decreasing = TRUE)[1:7],]

# other.line.freqs<-line.freqs[order(line.freqs$Freq, decreasing = TRUE),]
# levels(other.line.freqs[,"Var1"])<-c(levels(other.line.freqs[,"Var1"]),"Others")
# other.line.freqs[8:nrow(other.line.freqs),"Var1"]<-as.factor("Others")

#1 see rapresantated family in distal cpg minus group

dgm.line.freqs<-as.data.frame(table(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"family"]))

#subset only l1 family
dgm.line.freqs<-dgm.line.freqs[which( dgm.line.freqs$Var1 %in% line.freqs$Var1),]


full.table.lines.abundancy<-merge(line.freqs, dgm.line.freqs, by=1 ,all.x=TRUE)

# Ordina il dataframe in base a Freq.y decrescente, escludendo gli NA
full.table.lines.abundancy.sorted <- full.table.lines.abundancy[order(-full.table.lines.abundancy$Freq.y, na.last = TRUE), ]
to_export<-full.table.lines.abundancy.sorted
colnames(to_export)<-c("family","overall","dcm_loss_k4me3")

full.table.lines.abundancy.sorted$Freq.y[is.na(full.table.lines.abundancy.sorted$Freq.y)] <- 0


# Prendi i nomi Var1 delle prime 5 righe con Freq.y più alto
top5 <- full.table.lines.abundancy.sorted$Var1[order(-full.table.lines.abundancy.sorted$Freq.y)][1:5]


# Converti Var1 in carattere, poi assegna "Others" a tutti i valori che non sono nei top5
full.table.lines.abundancy.sorted$Var1 <- as.character(full.table.lines.abundancy.sorted$Var1)
full.table.lines.abundancy.sorted$Var1[!(full.table.lines.abundancy.sorted$Var1 %in% top5)] <- "Others"

# Riconverti Var1 in fattore
full.table.lines.abundancy.sorted$Var1 <- factor(full.table.lines.abundancy.sorted$Var1)



library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Calcola percentuali per Freq.x e Freq.y
df_plot <- full.table.lines.abundancy.sorted %>%
  group_by(Var1) %>%
  summarise(
    Freq.x = sum(Freq.x),
    Freq.y = sum(Freq.y)
  ) %>%
  mutate(
    perc_x = Freq.x / sum(Freq.x) * 100,
    perc_y = Freq.y / sum(Freq.y) * 100
  ) %>%
  select(Var1, perc_x, perc_y)

# 2. Converti in formato lungo (long format) per ggplot
df_long <- df_plot %>%
  pivot_longer(
    cols = c(perc_x, perc_y),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = recode(Category,
                      perc_x = "Overall",
                      perc_y = "DCM_only")
  )

df_long$Category <- factor(df_long$Category, levels = c("Overall", "DCM_only"))
df_long$Var1 <- factor(df_long$Var1, levels = c("Others", setdiff(unique(df_long$Var1), "Others")))


# 3. Barplot affiancato
p_famfreq<-ggplot(df_long, aes(x = Var1, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Family abundancy",
    x = "Family",
    y = "Percentage",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outs/CHiP_postprocessing_line1_dist/families_frequency_plot.pdf",
       plot = p_famfreq, width = 10, height = 6)

write.table(to_export, file="outs/CHiP_postprocessing_line1_dist/families_frequency.tsv",sep="\t", 
col.names=TRUE, 
row.names = FALSE,
quote=FALSE,
)

























# #su.dgm.line.fraqs<-dgm.line.freqs[which(dgm.line.freqs$Freq>=
# #                                          as.numeric(summary(dgm.line.freqs$Freq)[2])),]

# #all.rapresented<-

# #subset distal cpg minus with L1 only family
# #line.freqs$Var1

# #bck_su.dgm.line.fraqs<-su.dgm.line.fraqs

# su.dgm.line.fraqs<-su.dgm.line.fraqs[which( su.dgm.line.fraqs$Var1 %in% line.freqs$Var1),]

# top.dgm.line.fraqs<-dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE)[1:7],]

# other.dgm.line.freqs<-dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE)[8:nrow(dgm.line.freqs)],]

# levels(other.dgm.line.freqs[,"Var1"])<-c(levels(other.dgm.line.freqs[,"Var1"]),"Others")

# other.dgm.line.freqs[1:nrow(other.dgm.line.freqs),"Var1"]<-as.factor("Others")

# #1 generate comparative plot

# laneFreq.dataplot<-as.data.frame(rbind(cbind(su.line.freqs,group=rep("all")),
#       cbind(su.dgm.line.fraqs,group=rep("only.dcm"))
#       ))


# top.laneFreq.dataplot<-as.data.frame(rbind(cbind(top.line.freqs,group=rep("all")),
#                                            cbind(top.dgm.line.fraqs,group=rep("only.dcm"))
# ))

# #prendi le top in dcm
# #top.laneFreq.dataplot<-cbind(top.dgm.line.fraqs,group=rep("only.dcm"))

# # ggplot(top.laneFreq.dataplot, aes(fill=Var1, y=Freq, x=group)) + 
# #   geom_bar(position="fill", stat="identity")

# other.laneFreq.dataplot<-as.data.frame(rbind(cbind(other.line.freqs,group=rep("all")),
#                                              cbind(other.dgm.line.freqs,group=rep("only.dcm"))
# ))

# length(unique(other.laneFreq.dataplot$Var1))

# ggplot(other.laneFreq.dataplot, aes(fill=Var1, y=Freq, x=group)) + 
#   geom_bar(position="fill", stat="identity")

# #====extract family percentage on two condition

# tot.all<-sum(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group=="all"),"Freq"])

# tot.only.dcm<-sum(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group=="only.dcm"),"Freq"])

# tar_families<-unique(other.laneFreq.dataplot$Var1)
# # [1] L1Md_F2 Lx8     Lx9     L1_Mus3 Lx7     L1_Mur3 L1_Mus1 Others  L1Md_T 
# #[10] L1Md_A  L2a     L2b     L1Md_Gf

# onlydcm_perc<-c()
# all_perc<-c()

# for(i in 1:length(tar_families)){

#   tarfm_num_only.dcm<-sum(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group=="only.dcm" & other.laneFreq.dataplot$Var1 == tar_families[i] ),"Freq"])

#   onlydcm_perc[i]<-tarfm_num_only.dcm/tot.only.dcm

#   tarfm_num_all<-sum(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group=="all" & other.laneFreq.dataplot$Var1 == tar_families[i] ),"Freq"])

#   all_perc[i]<-tarfm_num_all/tot.all

# }

# famfreq<-data.frame(family=as.character(tar_families),overall=all_perc,dcm_loss_k4me3=onlydcm_perc)

# write.table(famfreq, file="outs/CHiP_postprocessing_line1_dist/families_frequency.tsv",sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)

# gg_famfreq<-rbind(data.frame(family=as.character(tar_families),freq=all_perc,type=rep("overall",length(all_perc))),data.frame(family=as.character(tar_families),freq=onlydcm_perc,type=rep("dcm_loss_k4me3",length(onlydcm_perc))))

# library(dplyr)

# # Crei l’ordinamento desiderato, escludendo "Others"
# ordered_families <- gg_famfreq %>%
#   filter(type == "dcm_loss_k4me3", family != "Others") %>%
#   arrange(desc(freq)) %>%
#   pull(family)

# # Poi metti "Others" in testa
# gg_famfreq$family <- factor(gg_famfreq$family,
#   levels = c("Others", ordered_families)
# )

# # Ordini anche i type, così "overall" viene prima
# gg_famfreq$type <- factor(gg_famfreq$type, levels = c("overall", "dcm_loss_k4me3"))

# p_famfreq<-ggplot(data = gg_famfreq, aes(x = family, y = freq, fill = type)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylab("Frequency")

# ggsave("outs/CHiP_postprocessing_line1_dist/families_frequency_plot.pdf",
#        plot = p_famfreq, width = 10, height = 6)



