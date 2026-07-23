#============
#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")


#======= generate stats for chipseq results 





#========input files

refpath.fasta<-"/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa"

antimml2.annotated.peaksFile="/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/macs_broadpeaks/Anti-GFP.mLb.mkD.sorted_peaks.annotatePeaks.txt"

antimml2.peaksfile="in/test_chipseq_downstream/bedtools_window/coordinate.bed"
distal.peaksfile="in/test_chipseq_downstream/bedtools_window/distal_peaks.bed"
proximal.peaksfile="in/test_chipseq_downstream/bedtools_window/proximal_preaks.bed"

dis.CpG.plus.peaksFile="in/test_chipseq_downstream/bedtools_window/distal_CpG_plus_unique.bed"
dis.CpG.min.peaksFile="in/test_chipseq_downstream/bedtools_window/distal_CpG_minus.bed"

pro.CpG.plus.peaksFile="in/test_chipseq_downstream/bedtools_window/proximal_CpG_plus_unique.bed"
pro.CpG.min.peaksFile="in/test_chipseq_downstream/bedtools_window/proximal_CpG_minus.bed"

dko.pro.peaksFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_Double_KO_vs_F_F.bed"
dko.dis.peaksFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_Double_KO_vs_F_F.bed"

dko.pro.CpG.plus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_CpG_plus_Double_KO_vs_F_F.bed"
dko.pro.CpG.minus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/proximal_CpG_minus_Double_KO_vs_F_F.bed"

dko.dis.CpG.plus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_CpG_plus_Double_KO_vs_F_F.bed"
dko.dis.CpG.minus.peakFile="in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/distal_CpG_minus_Double_KO_vs_F_F.bed"

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

distal.peaks<-read.table(distal.peaksfile,sep="\t", header=FALSE)
total.distal<-nrow(distal.peaks)

proximal.peaks<-read.table(proximal.peaksfile,sep="\t", header=FALSE)
total.proximal<-nrow(proximal.peaks)

dis.CpG.plus.peaks<-read.table(dis.CpG.plus.peaksFile,sep="\t", header=FALSE)
total.dis.CpG.plus<-nrow(dis.CpG.plus.peaks)

dis.CpG.min.peaks<-read.table(dis.CpG.min.peaksFile, sep="\t", header=FALSE)
total.dis.CpG.min<-nrow(dis.CpG.min.peaks)

pro.CpG.plus.peaks<-read.table(pro.CpG.plus.peaksFile,sep="\t", header=FALSE)
total.pro.CpG.plus<-nrow(pro.CpG.plus.peaks)

pro.CpG.min.peaks<-read.table(pro.CpG.min.peaksFile, sep="\t", header=FALSE)
total.pro.CpG.min<-nrow(pro.CpG.min.peaks)

# Load webr
library(webr)
library(ggplot2)

#create dataset

stats <- structure(
  list(
    TSSdist = structure(
      c(1L, 1L, 2L, 2L),
      .Label = c("distal", "proximal"),
      class = "factor"
    ),
    NMIsdist = structure(
      c(1L, 2L, 1L, 2L),
      .Label = c("NMIs_plus", "NMIs_minus"),
      class = "factor"
    ),
    peaks = c(total.dis.CpG.plus,
              total.dis.CpG.min,
              total.pro.CpG.plus,
              total.pro.CpG.min)
  ),
  .Names = c("TSSdist", "NMIsdist", "peaks"),
  row.names = c(NA, -4L),
  class = "data.frame"
)


# Create a Pie-Donut chart using webr library
PieDonut(stats, aes(x = TSSdist, y = NMIsdist, count = peaks),
         ratioByGroup = FALSE)

#==================

dko.pro.peaks<-read.table(dko.pro.peaksFile,sep="\t", header=FALSE)
tot.dko.pro<-nrow(dko.pro.peaks)

dko.dis.peaks<-read.table(dko.dis.peaksFile,sep="\t", header=FALSE)
tot.dko.dis<-nrow(dko.dis.peaks)

dko.pro.CpG.plus.peak<-read.table(dko.pro.CpG.plus.peakFile, sep="\t", header=FALSE)
tot.dko.pro.CpG.plus<-nrow(dko.pro.CpG.plus.peak)

dko.pro.CpG.minus.peak<-read.table(dko.pro.CpG.minus.peakFile, sep="\t", header=FALSE)
tot.dko.pro.CpG.minus<-nrow(dko.pro.CpG.minus.peak)

dko.dis.CpG.plus.peak<-read.table(dko.dis.CpG.plus.peakFile, sep="\t", header=FALSE)
tot.dko.dis.CpG.plus<-nrow(dko.dis.CpG.plus.peak)

dko.dis.CpG.minus.peak<-read.table(dko.dis.CpG.minus.peakFile, sep="\t", header=FALSE)
tot.dko.dis.CpG.minus<-nrow(dko.dis.CpG.minus.peak)


stats2 <- structure(
  list(
    subgroups = structure(
      c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L),
      .Label = c("distal_GpG_plus", "distal_GpG_minus", "proximal_CpG_plus","proximal_CpG_minus"
                 ),
      class = "factor"
    ),
    diffpeaks = structure(
      c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
      .Label = c("anti_mll2_only", "dko_diffpeaks"),
      class = "factor"
    ),
    peaks = c(total.dis.CpG.plus-tot.dko.dis.CpG.plus,
              total.dis.CpG.min-tot.dko.dis.CpG.minus,
              total.pro.CpG.plus-tot.dko.pro.CpG.plus,
              total.pro.CpG.min-tot.dko.pro.CpG.minus,
              tot.dko.dis.CpG.plus,
              tot.dko.dis.CpG.minus,
              tot.dko.pro.CpG.plus,
              tot.dko.pro.CpG.minus)
  ),
  .Names = c("subgroups", "diffpeaks", "peaks"),
  row.names = c(NA, -8L),
  class = "data.frame"
)

PieDonut(stats2, aes(x = subgroups, y = diffpeaks, count = peaks),
         ratioByGroup = FALSE)

#==================

stats3 <- structure(
  list(
    TSSdist = structure(
      c(1L, 1L, 2L, 2L),
      .Label = c("distal", "proximal"),
      class = "factor"
    ),
    NMIsdist = structure(
      c(1L, 2L, 1L, 2L),
      .Label = c("NMIs_plus", "NMIs_minus"),
      class = "factor"
    ),
    peaks = c(tot.dko.dis.CpG.plus,
              tot.dko.dis.CpG.minus,
              tot.dko.pro.CpG.plus,
              tot.dko.pro.CpG.minus)
  ),
  .Names = c("TSSdist", "NMIsdist", "peaks"),
  row.names = c(NA, -4L),
  class = "data.frame"
)


# Create a Pie-Donut chart using webr library
PieDonut(stats3, aes(x = TSSdist, y = NMIsdist, count = peaks),
         ratioByGroup = FALSE)

repeat_mask<-read.table("in/ucsc/rmsk.txt",sep="\t",header=FALSE)
nrow(repeat_mask)

print(unique(repeat_mask$V12))

#  [1] "LINE"           "Simple_repeat"  "LTR"            "SINE"           "Low_complexity"
#  [6] "DNA"            "snRNA"          "Other"          "Satellite"      "Unknown"       
# [11] "SINE?"          "srpRNA"         "tRNA"           "LTR?"           "RNA"           
# [16] "scRNA"          "RC"             "rRNA"           "DNA?"           "RC?"           
# [21] "LINE?"          "Retroposon"

repeat_mask_line<-repeat_mask[grep("\\b(LINE|LINE?)\\b",repeat_mask$V12),]
repeat_mask_line<-repeat_mask_line[,c(6,7,8,10,11,12,13)]

colnames(repeat_mask_line)<-c("chr","start","end","strand","family","element","type")


#GC content evaluation

library(seqinr)

ref.fasta<-seqinr::read.fasta(file=refpath.fasta)

allfml<-unique(repeat_mask_line$family)
fml.content<-c()
fprime.content<-c()
fprime.window<-300

for (i in 1:length(allfml)){
  print(paste0(i,"/",length(allfml)))
  print(allfml[i])
  tarfml<-allfml[i]
  tar_repeat_mask_line<-repeat_mask_line[which(repeat_mask_line$family==tarfml),]
  
  infml.content<-c()
  infprime.content<-c()
  
  for(j in 1:nrow(tar_repeat_mask_line)){
    
    seq<-ref.fasta[[tar_repeat_mask_line[j,1]]][tar_repeat_mask_line[j,"start"]:tar_repeat_mask_line[j,"end"]]
    
    if(length(seq)>fprime.window){
      if(tar_repeat_mask_line[j,"strand"]=="-"){
        fprimseq<-seq[(length(seq)-fprime.window+1):length(seq)]
      }else{
        if(tar_repeat_mask_line[j,"strand"]=="+"){
          fprimseq<-seq[1:fprime.window+1]
        }else{
          fprimseq<-NA
        }
      }      
    }else{
      fprimseq<-NA
    }
    
  if(is.null(seq)==TRUE){
    infml.content[j]<-NA
    infprime.content[j]<-NA
  }else{
    a.content<-length(which(seq=="a" | seq=="A"))
    t.content<-length(which(seq=="t" | seq=="T"))
    g.content<-length(which(seq=="g" | seq=="G"))
    c.content<-length(which(seq=="c" | seq=="C"))
    
    gprime.content<-length(which(fprimseq=="g" | fprimseq=="G"))
    cprime.content<-length(which(fprimseq=="c" | fprimseq=="C"))
    
    #gc.content<-(g.content+c.content)/(a.content+t.content+g.content+c.content)  
    gc.content<-(g.content+c.content)/(tar_repeat_mask_line[j,"end"] - tar_repeat_mask_line[j,"start"])
    gc.prime.content<-(gprime.content+cprime.content)/fprime.window
    
    infml.content[j]<-gc.content
    infprime.content[j]<-gc.prime.content
  }  
    
  }
  
  infml.content<-infml.content[which(!is.na(infml.content))]
  fml.content[i]<-median(infml.content)
  infprime.content<-infprime.content[which(!is.na(infprime.content))]
  fprime.content[i]<-median(infprime.content)
}


fml.cnt<-data.frame(family=allfml, gc=fml.content, gcprime=fprime.content)
fml.cnt<-fml.cnt[which(!is.na(fml.cnt$gc)),]

fml.cnt<-fml.cnt[order(fml.cnt$gc),]

plot(density(fml.cnt$gc))

# gene_id<-c()
# for (i in 1:nrow(repeat_mask)){
#   if(i %% 10000 == 0){
#     cat("rowa: ",i,"\n")
#    }
#   
#   splitarray<-strsplit(repeat_mask$V9[i],";")
#   gene_id[i]<-splitarray[[1]][grep("gene_id",splitarray)]
#   
# }



#gene_id2<-gsub("gene_id ","",gene_id)
#gene_id2<-gsub("\\(","",gene_id2)
#gene_id2<-gsub("\\)","_",gene_id2)

#repeat_mask<-cbind(repeat_mask,gene_id2)


#find distance vs repeat mask elements
#distance by 5' 
#divide by strand

repeat_mask_line_plus<-repeat_mask_line[which(repeat_mask_line$strand=="+"),]
repeat_mask_line_minus<-repeat_mask_line[which(repeat_mask_line$strand=="-"),]

print("nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line)")
nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line)

repeat_mask.5prime<-rbind(cbind(repeat_mask_line_plus, fprime=repeat_mask_line_plus$start),
cbind(repeat_mask_line_minus, fprime=repeat_mask_line_minus$end))

repeat_mask.coordinates<-cbind(repeat_mask.5prime$chr,repeat_mask.5prime$fprime)
repeat_mask.coordinates_extended=cbind(repeat_mask.5prime$chr,
                                       repeat_mask.5prime$fprime,
                                       repeat_mask.5prime$family
                                       )

extractdistances<-function(mean.coords, l1mask.meancoords, addfamily=FALSE){
  distances<-c()
  family<-NULL
  for(j in 1:nrow(mean.coords)){
    if(j %% 1000 == 0){
      cat("j: ",j,"\n")
    }
    tarreg<-mean.coords[j,]
    l1mask.mc.sub<-l1mask.meancoords[which(l1mask.meancoords[,1]==tarreg[1]),]
    
    allvals<-abs(as.numeric(as.character(l1mask.mc.sub[,2]))
        -as.numeric(as.character(tarreg[2])))
    
    distances[j]<-min(allvals)
    
    if(addfamily==TRUE){
      family[j]<-l1mask.mc.sub[which(allvals==distances[j])[1],3]
    }
    
    
    
  }
  
  
  if(addfamily==TRUE){
    return(cbind(distances, family)) 
  }else{
    return(distances) 
  }
  
  
}

#calc distribution of distances between  

atb<-c("K4me3","K4me2")
distances_list<-NULL
for (i in 1:length(atb)){
  
  #distal
  file_pattern="_distal_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks, atb[i],file_pattern),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder,atb[i],"_distal_Double_KO_vs_F_F.pdf"))
  plot(density(distances),
       main=paste0(atb[i],"_distal_Double_KO_vs_F_F nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()
  
  file_pattern="_distal_CpG_minus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(distances),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()
  
  file_pattern="_distal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(distances),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()
  
  #proximal
  d<-read.table(paste0(infolder.mll2peaks,atb[i],"_proximal_Double_KO_vs_F_F.bed"),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder,atb[i],"_proximal_Double_KO_vs_F_F.pdf"))
  plot(density(distances),
       main=paste0(atb[i],"_proximal_Double_KO_vs_F_F nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()
  
  file_pattern="_proximal_CpG_minus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder,atb[i],file_pattern,".pdf"))
  plot(density(distances),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()
  
  file_pattern="_proximal_CpG_plus_Double_KO_vs_F_F.bed"
  d<-read.table(paste0(infolder.mll2peaks,atb[i],file_pattern),sep="\t")
  mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))
  
  distances<-extractdistances( mean.coords=mean.coords,  l1mask.meancoords=repeat_mask.coordinates)
  distances_list[length(distances_list)+1]<-paste0(atb[i],file_pattern)
  assign(paste0(atb[i],file_pattern),distances)
  
  pdf(paste0(outfolder, atb[i],file_pattern,".pdf"))
  plot(density(distances),
       main=paste0(atb[i],file_pattern," nearest LINE"),
       sub=paste0("median=",median(distances)))
  abline(v=median(distances))
  dev.off()

}

for( i in 1:length(distances_list)){
  
  if ( i == 1 ){
    final<-get(distances_list[i])
    final<-cbind(final,rep(distances_list[i]))
  }else{
    add<-get(distances_list[i])
    add<-cbind(add,rep(distances_list[i]))
    final<-rbind(final,add)
    
  }
  
}

densities<-data.frame(val=as.numeric(final[,1]),subgroup=final[,2])

densities<-densities[which(densities$subgroup=="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed" | densities$subgroup=="K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed" | densities$subgroup=="K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed" |  densities$subgroup=="K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed" ),]

densities$subgroup<-gsub("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed","distal_CpG_minus",densities$subgroup)
densities$subgroup<-gsub("K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed","distal_CpG_plus",densities$subgroup)
densities$subgroup<-gsub("K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed","proximal_CpG_minus",densities$subgroup)
densities$subgroup<-gsub("K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed","proximal_CpG_plus",densities$subgroup)

library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

ggplot( data=densities, aes(x=log(val), color=subgroup)) +
  geom_density(alpha=0.6) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  #geom_text( data=annot, aes(x=x, y=y, label=text, color=text), hjust=0, size=4.5) +
  theme_ipsum() +
  geom_vline(xintercept=5.48)+
  geom_vline(xintercept=4.80) +
  #theme(
  #  legend.position="none"
  #) +
  ylab("") +
  xlab("Log distance by LINE")

file_pattern="K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
d<-read.table(paste0(infolder.mll2peaks,file_pattern),sep="\t")
mean.coords<-cbind(d$V1,d$V2+round((d$V3-d$V2)/2))

targetlines<-which(repeat_mask.coordinates_extended[,3]=="L1Md_T" 
                     )

distances<-extractdistances( mean.coords = mean.coords,  
                             l1mask.meancoords = repeat_mask.coordinates_extended,
                             addfamily = TRUE
                             )
family_dist<-data.frame(val=as.numeric(distances[,1]),family=distances[,2])

family_dist<-family_dist[order(family_dist$val, decreasing = TRUE),]
family_dist_su<-family_dist[1:20,]
#color=family
ggplot( data=family_dist_su, aes(x=log(val))) +
  geom_density(alpha=0.6) +
  #scale_fill_viridis(discrete=TRUE) +
  #scale_color_viridis(discrete=TRUE) +
  #geom_text( data=family_dist, aes(x=log(val), y= ,label=family), hjust=0, size=1.5) +
  theme_ipsum() +
  geom_vline(xintercept=5.48)+
  theme(
   legend.position="none"
  ) +
  ylab("") +
  xlab("Log distance by LINE family")

fmls<-unique(family_dist$family)

ggplot( family_dist,  aes(x=val, y=family, fill=family)) + 
  geom_boxplot() +
  xlab("class") +
  theme(legend.position="none",
        axis.text.y = element_text(size=6.5)
        ) +
  xlab("")

nsamples<-c()
median<-c()
for (i in 1:length(fmls)){
  idx<-which(family_dist$family==fmls[i])
  nsamples[i]<-length(idx)
  median[i]<-median(family_dist$val[idx])
}

family_summary<-data.frame(family=fmls,
                           nsamples=as.numeric(as.character(nsamples)),
                           median=as.numeric(as.character(median)))

summary(family_summary$nsamples)

target_family<-family_summary[which(family_summary$nsamples>=10),]
target_family<-target_family[which(target_family$median<=1000),]

family_sub<-family_dist[which(family_dist$family %in% target_family$family),]

ggplot( data=family_sub, aes(x=log2(val), color=family )) +
  geom_density(alpha=0.6) +
  #scale_fill_viridis(discrete=TRUE) +
  #scale_color_viridis(discrete=TRUE) +
  #geom_text( data=family_dist, aes(x=log(val), y= ,label=family), hjust=0, size=1.5) +
  theme_ipsum() +
  # theme(
  #   legend.position="none"
  # ) +
  ylab("") +
  xlab("Log distance by LINE family")

#selected family gc content

inte_family<-unique(family_sub$family)
fml.cnt[fml.cnt$family %in% inte_family,]



#target_family<-c("L1Md_T", "L1Md_A" ,"L1Md_Gf")
l1mask<-repeat_mask[grep("\\b(L1Md_T|L1Md_A|L1Md_Gf)\\b",repeat_mask$gene_id),]

l1mask.meancoords<-cbind(l1mask$V1,l1mask$V4+round((l1mask$V5-l1mask$V4)/2))


write.table(cbind(l1mask[,c(1,4,5)],paste0(l1mask[,10],"_",l1mask[,7])),
            file="l1md_a_t_gf.bed",sep="\t",
            row.names =FALSE,
            quote=FALSE,
            col.names=FALSE
)





pcalled.m2<-read.table("in/test_chipseq_downstream/diffbind/fold2_th0_05_K4me2_Double_KO_vs_F_F_DBA_DESEQ2_results.txt",
                       sep="\t",
                       header=TRUE
)

pcalled.m3<-read.table("in/test_chipseq_downstream/diffbind/fold2_th0_05_K4me3_Double_KO_vs_F_F_DBA_DESEQ2_results.txt",
                       sep="\t",
                       header=TRUE
)

nmi<-read.table("in/test_chipseq_downstream/bedtools_window/NMIs_mm10_mESC_afterLiftover_from_mm9_to_mm10.bed", 
                sep="\t", header=FALSE)

tss<-read.table("in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/TSS.bed",
                sep="\t", header=FALSE
)

pcall.m2.down<-pcalled.m2[which(pcalled.m2$Fold<0),]
meancoord_pcall.m2.down<-cbind(pcall.m2.down, 
                               meancoord=pcall.m2.down$start+(pcall.m2.down$end-pcall.m2.down$start)/2)

final.mnisel.mean<-NULL
final.mnisel.chr<-NULL

for (i in 1:nrow(meancoord_pcall.m2.down)){
  
  print(i)
  
  check.coord<-meancoord_pcall.m2.down[i,]
  
  nmisel<-nmi[which(nmi$V1 == check.coord$seqnames[1] ),]
  
  #nmisel[min(which(nmisel$V2 >= check.coord$meancoord[1])),]
  
  #nmisel<-nmisel[which(nmisel$V2 >= check.coord$meancoord[1] & nmisel$V3 < check.coord$meancoord[1]),]
  
  mnisel.mean<-nmisel$V2+((nmisel$V3-nmisel$V2)/2)
  mnisel.chr<-nmisel$V1
  
  nmisel<-nmi[which(nmi$V1==check.coord$seqnames[1]),]
  
  if( i !=1 ){
    final.mnisel.chr<-c(pre.chr,final.mnisel.chr)
    final.mnisel.mean<-c(pre,final.mnisel.mean)
  }else{
    pre<-mnisel.mean
    pre.chr<-mnisel.chr
  }
  
  
  
  # if(nrow(nmisel)!=0){
  #   # distances[length(distances)+1]
  #   # 
  #   # for(j in 1:nrow(nmisel)){
  #   #   
  #   # }
  #   print(i)
  # }
  
  
  
}


