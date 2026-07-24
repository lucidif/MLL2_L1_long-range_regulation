#setwd("/media/lucio/external.wk/bioinfo/wkdir/Lara/Lara_multiomic_analysis")
#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis")

make.anch<-function(outfile,
                    anch1.path,
                    peaks.path){

  anch1<-read.table(anch1.path,
                    sep="\t",
                    header=FALSE
  )

  proximal.peaks<-read.table(peaks.path,
                             sep="\t",
                             header=FALSE)

  anch3_unsort<-merge(anch1,
                      proximal.peaks,
                      by.x=8,
                      by.y=1)

  #anch3<-cbind(anch3, anch3$V10)
  anch3_unsort<-anch3_unsort[,c(-1,-8,-9)]
  #anch3<-anch3[,c(1,2,3,5,6,7,4,8)]

  colnames(anch3_unsort)<-c("chrA","startA","endA","chrB","startB","endB")

  anch3<-anch3_unsort

  for (i in 1:nrow(anch3_unsort)){

    meanA <- round((anch3_unsort$startA[i] + ((anch3_unsort$endA[i] - anch3_unsort$startA [i] )/2)),0)
    meanB <- round((anch3_unsort$startB[i] + ((anch3_unsort$endB[i] - anch3_unsort$startB [i] )/2)),0)

    if(meanB<meanA){
      anch3[i,1]<-anch3_unsort[i,4]
      anch3[i,2]<-anch3_unsort[i,5]
      anch3[i,3]<-anch3_unsort[i,6]
      anch3[i,4]<-anch3_unsort[i,1]
      anch3[i,5]<-anch3_unsort[i,2]
      anch3[i,6]<-anch3_unsort[i,3]
    }


  }

  write.table(anch3,
              outfile,
              quote=FALSE,
              sep="\t",
              col.names = FALSE,
              row.names = FALSE
  )


}


make.anch.fromTsv<-function(outfile,
                    anch1.path,
                    peaks.path){

  anch1<-read.table(anch1.path,
                    sep="\t",
                    header=FALSE
  )

  proximal.peaks<-read.table(peaks.path,
                             sep="\t",
                             header=FALSE)

  anch3_unsort<-merge(anch1,
                      proximal.peaks,
                      by.x=10,
                      by.y=1)

  #anch3<-cbind(anch3, anch3$V10)
  #anch3_unsort<-anch3_unsort[,c(-1,-8,-9)]
  anch3_unsort<-cbind(anch3_unsort[,-1],anch3_unsort[,1])
  anch3_unsort<-anch3_unsort[,-10]
  #anch3<-anch3[,c(1,2,3,5,6,7,4,8)]

  colnames(anch3_unsort)<-c("peak.chr","peak.start","peak.end","peak.strand","peak",
  "gene.chr","gene.start","gene.end","gene.strand","gene")

#  anch3<-anch3_unsort

#   for (i in 1:nrow(anch3_unsort)){

#     meanA <- round((anch3_unsort$startA[i] + ((anch3_unsort$endA[i] - anch3_unsort$startA [i] )/2)),0)
#     meanB <- round((anch3_unsort$startB[i] + ((anch3_unsort$endB[i] - anch3_unsort$startB [i] )/2)),0)

#     if(meanB<meanA){
#       anch3[i,1]<-anch3_unsort[i,4]
#       anch3[i,2]<-anch3_unsort[i,5]
#       anch3[i,3]<-anch3_unsort[i,6]
#       anch3[i,4]<-anch3_unsort[i,1]
#       anch3[i,5]<-anch3_unsort[i,2]
#       anch3[i,6]<-anch3_unsort[i,3]
#     }


#   }

  write.table(anch3_unsort,
              outfile,
              quote=FALSE,
              sep="\t",
              col.names = TRUE,
              row.names = FALSE
  )


}




make.anch(outfile = "outs/coolpup/500bp/win500_anchors3.bedpe",
          anch1.path = "outs/coolpup/500bp/win500_anchors1.bedpe",
          peaks.path = "outs/great/Double_KO_vs_F_F/basal/proximal/20241126-public-4.0.4-xCT3F1-mm10-all-gene.txt"
          )

make.anch(outfile = "outs/coolpup/500bp/win1_anchors3.bedpe",
          anch1.path = "outs/coolpup/500bp/win1_anchors1.bedpe",
          peaks.path = "outs/great/Double_KO_vs_F_F/basal/proximal/20241126-public-4.0.4-xCT3F1-mm10-all-gene.txt"
)

make.anch(outfile = "outs/coolpup/500bp/win1000_anchors3.bedpe",
          anch1.path = "outs/coolpup/500bp/win1000_anchors1.bedpe",
          peaks.path = "outs/great/Double_KO_vs_F_F/basal/proximal/20241126-public-4.0.4-xCT3F1-mm10-all-gene.txt"
)

#TODO add nowin ad bedpe (reorder with strand at the end) TODO devi riordinare nowin secondo dei criteri che make.anch possa usare

make.anch.fromTsv(outfile = "outs/coolpup/500bp/nowin_unsorted_anchors3.tsv",
                anch1.path = "outs/coolpup/500bp/nowin_anchors1_unsorted_withStrand.tsv",
                peaks.path = "outs/great/Double_KO_vs_F_F/basal/proximal/20241126-public-4.0.4-xCT3F1-mm10-all-gene.txt"
)

