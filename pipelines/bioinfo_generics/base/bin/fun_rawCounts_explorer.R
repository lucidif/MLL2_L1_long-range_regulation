countingStats<-function(countsFile,
                        cpmThr=0,
                        ggplotForm=FALSE,
                        topThs=c(1,5,20),
                        rownames=TRUE, #If there are row names in input files (like output of nf-core mir)
                        gsutil=FALSE
                      ){
  #countsFile<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2020/5_HTFT7DRXX/Extracellular_vesicle_RNA_Extracellular_vesicle_RNA_rawCounts.txt"

  #countingStats("gs://ngdx-runs/Analysis/Other/PDCD10/mirna_seq_Analysis/raw_counts/miRBase_mature/mature_raw_read_counts.txt", gsutil=TRUE)
  
  #Es. countingStats(countsFile,0)

  if(gsutil==TRUE){
    
    id<-sample(1000000:9999999)[1]
    
    if (file.exists(paste0(id,"_mirnaRwctCounts.txt"))){
      file.remove(paste0(id,"_mirnaRwctCounts.txt"))
    }
    
    name<-paste0(id,"_mirnaRwctCounts.txt")
    system(paste("gsutil cp",countsFile,name,sep=" "))
    
    countsTable<-read.table(name, sep="\t", header=TRUE, check.names = FALSE)
    
    file.remove(paste0(id,"_mirnaRwctCounts.txt"))
    
  }else{
    countsTable<-read.table(countsFile, sep="\t", header=TRUE, check.names = FALSE)
  }
  
  
  

  countsTable2<-countsTable

  if(rownames==FALSE){
    row.names(countsTable2)<-countsTable2[,1]
    countsTable2<-countsTable2[,-1]

  }

  totalCnt<-c()
  sNames<-c()
  for(i in 1:length(countsTable2[1,])){
    sNames[i]<-colnames(countsTable2)[i]
    #detectedMirna[i]<-length(subset(contestabileCnt2[,i],contestabileCnt2[,i]>0))
    totalCnt[i]<-round(sum(countsTable2[,i]),0)
  }

  total<-data.frame(sNames,totalCnt)
  colnames(total)<-c("sample","value")

  normCnt<-edgeR::cpm(countsTable2)

  #cpmThr<-0
  sNames<-c()
  detectedMirna<-c()
  print(length(normCnt[,1]))
  #1966
  totalCnt<-c()
  for(i in 1:length(normCnt[1,])){
    sNames[i]<-colnames(normCnt)[i]
    filtered<-subset(normCnt[,i],as.numeric(normCnt[,i])>cpmThr)
    #filtered<-subset(filtered,as.numeric(filtered)<100)
    #detectedMirna[i]<-length(subset(conormCnt[,i],as.numeric(conormCnt[,i])>cpmThr))
    detectedMirna[i]<-length(filtered)


  }
  detected<-data.frame(sNames,detectedMirna)
  colnames(detected)<-c("sample","value")

  #perc Expression

  sNames<-c()
  t1<-c()
  t5<-c()
  t20<-c()
  for (i in 1:length(normCnt[1,])){
    sNames[i]<-colnames(normCnt)[i]
    sortvalues<-normCnt[order(normCnt[,i],decreasing =TRUE),i]

    tot<-sum(normCnt[,i])
    #first gene
    first<-sortvalues[1:topThs[1]]
    ftop<-(first/tot)*100
    t1[i]<-ftop
    #top 5
    top5<-sum(sortvalues[1:topThs[2]])
    top5<-(top5/tot)*100
    t5[i]<-top5
    #top 20
    top20<-sum(sortvalues[1:topThs[3]])
    top20<-(top20/tot)*100
    t20[i]<-top20
  }

  top1table<-data.frame(sNames,t1)
  colnames(top1table)<-c("sample","value")

  top5table<-data.frame(sNames,t5)
  colnames(top5table)<-c("sample","value")

  top20table<-data.frame(sNames,t20)
  colnames(top20table)<-c("sample","value")

  #print(total)
  #print(detected)

  if(ggplotForm==FALSE){
    final<-merge(total,detected,by="sample")
    final<-merge(final,top1table,by="sample")
    final<-merge(final,top5table,by="sample")
    final<-merge(final,top20table,by="sample")
    colnames(final)<-c("sample","totalCounts","detectedGene","top1","top5","top20")
  }else{
    measure<-c(rep("totalReads",length(normCnt[1,])),
               rep("detected",length(normCnt[1,])),
               rep("1_top1",length(normCnt[1,])),
               rep("2_top5",length(normCnt[1,])),
               rep("3_top20",length(normCnt[1,]))
               )
    final<-rbind(total,detected,top1table,top5table,top20table)
    final<-cbind(final,measure)
  }

  return(final)

}


topPlot<-function(countingFunResults, dodge=TRUE, order=FALSE){

  #countsFile<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2020/5_HTFT7DRXX/Extracellular_vesicle_RNA_Extracellular_vesicle_RNA_rawCounts.txt"
  #countingFunResults<-countingStats(countsFile,0,ggplotForm=TRUE)
  varpalette<-rev(viridis::viridis(3))
  cntStats<-countingFunResults

  topStats<-cntStats[grep("top",cntStats$measure),]

  if(order==TRUE){
    topStats_sorted<-data.frame(
        reorder(topStats$sample,topStats$value),
        reorder(topStats$value,topStats$value),
        reorder(topStats$measure,topStats$value)

    )

    colnames(topStats_sorted)<-colnames(topStats)
    topStats_sorted[,2]<-as.numeric(as.character(topStats_sorted[,2]))
    topStats<-topStats_sorted
  }



  library(ggplot2)
  ggplot(data=topStats, aes(x=sample, y=value, fill=measure)) +

    if(dodge==TRUE){
      geom_bar(stat="identity",position=position_dodge())
    }else{
      geom_bar(stat="identity",position="identity", alpha=1)
    }



    #scale_fill_hue(l=10, c=40)
    #scale_color_manual(values=c("#FDE725FF", "#21908CFF", "#440154FF"))
    #scale_fill_manual(breaks =topStats$measure,
    #                  values=varpalette)


}





#Francesca Vanni

#exaRes<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2020/5_HTFT7DRXX/Extracellular_vesicle_RNA_Extracellular_vesicle_RNA_rawCounts.txt"


#Pietro Carotenuto

#exaRes<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2020/5_HTFT7DRXX/smallRNA_CCA_16_10_2020_rawCounts.txt"

#Contestabile

#exaRes<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2019/1/miR_mmu_1/oldContestabileMirna_rawCounts.txt"


#Tarallo

#exaRes<-"/home/lucio/MEGA/bioinformatics/data/miRNA/2020/3/nfcorepipe/miRNA_Tarallo2_newNGDpipe_rawCounts.txt"


#View(countingStats(exaRes))
