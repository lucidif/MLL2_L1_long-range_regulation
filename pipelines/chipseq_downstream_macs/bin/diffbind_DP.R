diffbind_DP<-function(inssfile, antibodies, comparisons, 
                      out="./" ,th=0.05, fold=NULL,
                      method="DBA_DESEQ2"
                      ){
  

  # inssfile<-"Diffbind_SS.csv"
  # antibodies<-c("K4me2","K4me3")
  # 
  # inss<-inss[which(inss$Antibody%in%antibodies),]
  # 
  # comparisons<-data.frame(
  #   c1=c(numerator="FC_FC",denominator="F_F"),
  #   c2=c(numerator="Mll1-KO",denominator="F_F"),
  #   c3=c(numerator="Double_KO",denominator="F_F")
  # )
  #comparisons<-data.frame( c1=c(numerator="FC_FC",denominator="F_F"), c2=c(numerator="Mll1-KO",denominator="F_F"), c3=c(numerator="Double_KO",denominator="F_F") )
  #fold=0.01

  library(DiffBind)
  
  inss<-read.table(inssfile,sep=",",header=TRUE)
  #complete.dba=dba(sampleSheet=inss)
  #dba.plotPCA(complete.dba,DBA_CONDITION,label=DBA_ID,main="all_samples")
  inss<-inss[which(inss$Antibody%in%antibodies),]
  
  antb<-unique(inss$Antibody)
  
  for (i in 1:length(antb)){
    tar.ant.smp<-inss[which(inss$Antibody==antb[i]),] 
    tar.ant.dba<-dba(sampleSheet=tar.ant.smp)
    
    tar.ant.dba <- dba.count(tar.ant.dba)

    pdf(file = paste0(antb[i],"_pca.pdf"))
    dba.plotPCA(tar.ant.dba,DBA_CONDITION,label=DBA_REPLICATE)
    dev.off()
    
    tar.ant.dba <- dba.normalize(tar.ant.dba)
    # tar.ant.dba <- dba.contrast(tar.ant.dba,
    #                             reorderMeta=list(Condition=),
    #                             minMembers=2)
    # tar.ant.dba<-dba.analyze(tar.ant.dba, bGreylist=FALSE)
    #tar.ant.report <- dba.report(tar.ant.dba)
    
    
  }
  
  
  #divide by antibody (histone modification)
  
  atb=unique(inss$Antibody)
  models.ls<-c()
  analysis.ls<-c()
  peakfiles=read.table(inss$Peaks[1],sep="\t",header=FALSE)
  for (i in 1:length(atb)){
    tarss<-inss[which(inss$Antibody==atb[i]),]
    
    for (j in 1:length(comparisons)){
      
      tarnum=comparisons[1,j]
      tarden=comparisons[2,j]
      
      compname=paste0("th",gsub("\\.","_",th),"_",atb[i],"_",tarnum, "_vs_", tarden,"_",method)
      
      if(is.null(fold)==FALSE){
        compname<-paste0("fold",gsub("\\.","_",fold),"_",compname)
      }
      
      compss<-tarss[c(which(tarss$Condition==tarnum),which(tarss$Condition==tarden)),]
      
      tardba<-dba(sampleSheet=compss)
      tardba <- dba.count(tardba, bUseSummarizeOverlaps=TRUE)
      #tardba <- dba.peakset(tardba, peaks=)
      pdf(file = paste0(compname,"_pca.pdf"))
      dba.plotPCA(tardba,DBA_CONDITION,label=DBA_REPLICATE)
      dev.off()
      
      tardba <- dba.normalize(tardba)
      tardba<-dba.contrast(tardba,
                           reorderMeta=list(Condition=tarden),
                           minMembers=2
      )
      

      
      tardba<-dba.analyze(tardba, bGreylist=FALSE, method=get(method))
      #tarba.report <- dba.report(tardba, method=DBA_DESEQ2, th=0.1, fold=log2(1.5))

      #tardba<-dba.analyze(tardba, bGreylist=FALSE)

      tarba.report <- dba.report(tardba, th=th,method=get(method) ,bNormalized=TRUE)
      df=as.data.frame(tarba.report)
      
      if (is.null(fold)==FALSE){
        df<-df[which(abs(df$Fold)>=log2(fold)),]
        df.up<-df[which(df$Fold>0),]
        df.down<-df[which(df$Fold<0),]  
        
      }

      # deseqresults<-as.data.frame(tardba[["contrasts"]][[1]][["DESeq2"]])
      # coords<-c()
      # for(i in 1:nrow(deseqresults)){
      #   
      #   if(length(which(df$Fold==deseqresults$de.fold[i] & df$FDR==deseqresults$de.padj[i] ))){
      #     z=which(df$Fold==deseqresults$de.fold[i])
      #     print("fold match")
      #     coords[i+1]<-paste0(df$seqnames[z],":",df$start[z],"-",df$end[z],"_",df$Fold[z],"_",df$FDR[z])
      #   } else{
      #     print("fold nomatch")
      #   }
      #   
      # } 
      
      write.table(df, paste0(out, "/", compname,"_results.txt"),
                  sep="\t",
                  col.names = TRUE,
                  quote=FALSE,
                  row.names = FALSE
      )

      write.table(df[,c(1,2,3)], paste0(out, "/", compname,".bed"),
                  sep="\t",
                  col.names = FALSE,
                  quote=FALSE,
                  row.names = FALSE
      )
      
      write.table(df.up[,c(1,2,3)], paste0(out, "/", compname,"_UP.bed"),
                  sep="\t",
                  col.names = FALSE,
                  quote=FALSE,
                  row.names = FALSE
      )
      
      write.table(df.down[,c(1,2,3)], paste0(out, "/", compname,"_DOWN.bed"),
                  sep="\t",
                  col.names = FALSE,
                  quote=FALSE,
                  row.names = FALSE
      )

      
      print (tardba)
      
      assign(paste0(compname,"_dba"),tardba)
      models.ls[length(models.ls)+1]<-compname
      analysis.ls[length(analysis.ls)+1]<-paste0(compname,"_dba")
      
    }
    
    #selraws=c(which(inss$Condition==),)
    
  }
  
  
  
  for (i in 1:length(analysis.ls)){
    
    pdf(file = paste0(out, "/",analysis.ls[i],"_venn.pdf"))
    dba.plotVenn(get(analysis.ls[i]), contrast=1, bDB=TRUE,
                 bGain=TRUE, bLoss=TRUE, bAll=FALSE, main=analysis.ls[i])
    dev.off()
    
    pdf(file = paste0(out, "/",analysis.ls[i],"_MA.pdf"))
    dba.plotMA(get(analysis.ls[i]))
    dev.off()
    pdf(file = paste0(out, "/",analysis.ls[i],"_volcano.pdf"))
    dba.plotVolcano(get(analysis.ls[i]), fold=log2(fold)
                    , method = get(method)
                    , th=th
                    )
    dev.off()
    
  }
  
}


