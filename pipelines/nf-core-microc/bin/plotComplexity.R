pattern=".lc_extrap.txt"
folder="/mnt/datawk1/analysis/Lara/2024_03_Lara_microC/nfout/preseq_lcextrap/"
in.files=list.files(folder,pattern)

samples=c()
treads=c()
expected=c()
for (i in 1:length(in.files)){
  add<-read.delim(paste0(folder,"/",in.files[i]))
  add<-cbind(add, sample=rep(gsub(".psq.lc_extrap.txt","",in.files[i]),nrow(add)))
  
  if(i==1){
    final<-add
  }else{
    final<-rbind(final,add)
  }
  
}

reference.path="/mnt/datawk1/analysis/Lara/microC_MLL/dovetail_pipe/preseq/out.preseq"
refsample=read.delim(reference.path)
refsample<-cbind(refsample, sample="ref_sample")

final<-rbind(final,refsample)

library(ggplot2)
library(gridExtra)

ggplot(final, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT, z=sample)) + 
  geom_path(aes(color = factor(sample))) +
  geom_vline(xintercept=300000000, linetype="dashed", 
             color = "red", size=1) +
  geom_segment(aes(xend = TOTAL_READS + 1, yend = EXPECTED_DISTINCT, color = factor(sample)))


samples=c("KO_day0","KO_day4","WT_day0","WT_day4")

for (i in 1:length(samples)){
  
  subfinal<-final[grep(samples[i],final$sample),]
  reffinal<-final[which(final$sample=="ref_sample"),]
  tarfinal<-rbind(reffinal, subfinal)
  
  gg<-ggplot(tarfinal, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT, z=sample)) + 
    geom_path(aes(color = factor(sample))) +
    geom_vline(xintercept=300000000, linetype="dashed", 
               color = "red", size=1) +
    geom_segment(aes(xend = TOTAL_READS + 1, yend = EXPECTED_DISTINCT, color = factor(sample))) +
    ggtitle(samples[i])
  
  assign(paste0("gg",i),gg)
  
}

grid.arrange(gg1,gg2,gg3,gg4,nrow=2)


splitrep=c("KO_day4_A","KO_day4_B")

for (i in 1:length(splitrep)){
  subfinal<-final[grep(splitrep[i],final$sample),]
  reffinal<-final[which(final$sample=="ref_sample"),]
  tarfinal<-rbind(reffinal, subfinal)
  
  gg<-ggplot(tarfinal, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT, z=sample)) + 
    geom_path(aes(color = factor(sample))) +
    geom_vline(xintercept=300000000, linetype="dashed", 
               color = "red", size=1) +
    geom_segment(aes(xend = TOTAL_READS + 1, yend = EXPECTED_DISTINCT, color = factor(sample))) +
    ggtitle(splitrep[i])
  
  assign(paste0("gg_split",i),gg)
}

grid.arrange(gg1,gg2,gg3,gg4,gg_split1,gg_split2,nrow=3)


