

#=======================
#funs list

#----diffExp_unpaired   # unpaired differential expression case
#----diffExp_paired     # paired differential expression case
#----pca_plot           # make pca plot
#----differential       #execute interly diffenretial expression workflow
#----statRawTable
#----filedgeR()         #flter gene by DEseq like strategy

#=======================

#samtools mpileup -uf /home/lucio/localData/references/references_index_BWA_GRCm39_GRCm39.primary_assembly.genome.fa /home/lucio/MEGA/bioinformatics/data/OtherAnalysis/2021/committers/Auricchio/210429_offTarget/bwa-mem/pool_1_47os.bam | bcftools call -mv > /home/lucio/MEGA/bioinformatics/data/OtherAnalysis/2021/committers/Auricchio/210429_offTarget/bcftools_call/var.raw.vcf

# pca_plot <- function(data, groups, gplabs = groups, pca_f, pcaloc = 'topright', legendpage = FALSE, monoLayout=TRUE, sampleNames=TRUE){
#   names(gplabs) <- groups
#   cc <- as.factor(groups)
#   names(cc) <- gplabs
#   pca_cl <- prcomp(t(log2(data + 1)))
#   pdf(pca_f)
#     # Add extra space to right of plot area; change clipping to figure
#     if (monoLayout){par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, pin=c(4.4,4))}
#     plot(pca_cl$x[, c(1,2)], pch = unclass(cc), col = unclass(cc))
#     if (sampleNames){
#       text(pca_cl$x[, c(1,2)], row.names(pca_cl$x), cex = 0.5, pos = 3)
#     }else{
#         text(pca_cl$x[, c(1,2)], cex = 0.5, pos = 3)
#     }
#     if (legendpage) {
#       plot.new()
#     }
#     if (monoLayout){
#     legend("topright", inset=c(-0.2,0) , legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA', cex=0.8)
#     }else{
#     legend(pcaloc, legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA')
#     }
#
#   dev.off()
#   # norm_data <- merge(ann, round(data, 2), by.x = 1, by.y = 0)
#   # return(norm_data)
# }

diffExp_unpaired <- function(data, groups,
                             gplabs = groups,
                             ann,
                             ngth = 4,
                             advFilter=FALSE) {

  print("parameters")
  print(paste("ngth=",ngth))
  library(edgeR)
  #print("start analysis")

  # data = rawc; groups = pjsamples$condition
  # if(gplabs == NULL) {
  #   gplabs <- groups
  # }
  ngth<-as.numeric(as.character(ngth))
  names(gplabs) <- groups
  # make design matrix
  th <- min(table(groups)) - 1
  #data <- as.matrix(data[rowSums(cpm(data) > 1) > th, ])

  #advFilter

  if(advFilter==TRUE){
    passGenes<-filedgeR(data=data,groups=groups,ann=ann)
    #data2<-data[passGenes,]
    #data2<-data
    rownames(data)<-ann[,1]
    data<-data[passGenes,]
    ann<-as.matrix(ann[match(passGenes,ann[,1]),])
    rownames(data)<-NULL
  }


  #========

  design <- model.matrix(~0 + groups)
  rownames(design) <- colnames(data)
  colnames(design) <- gsub("groups", "", colnames(design))
  y <- DGEList(counts=data,genes=ann[,1])
  y <- calcNormFactors(y)
  y <- estimateGLMRobustDisp(y,design)   #robust tagwise dispersions
  #gn<-  y$genes
  fit <- glmFit(y,design)
  ng <- length(unique(groups))



  if(ng <= ngth){
    vv <- seq(ng - 1)
  }else{
    vv <- grep('^CTR', colnames(design))
  }
  colnames(design) <- gplabs[colnames(design)]
  # vv
  #ann<-as.matrix(ann[intersect(rownames(ann),gn),1])
  #ann2<-as.matrix(ann[intersect(rownames(ann),gn),1])
  #firstcol<-colnames(ann)[1]

  for(i in vv){
    for(j in c((i+1):ng)){
      cont <- rep(0, ng)
      cont[i] = -1
      cont[j] = 1
      print(cont)
      lrt <- glmLRT(fit, contrast = cont)
      nr <- dim(data)[1]
      #res <- topTags(lrt, n = nr)$table[, c(1,4,5)]
      res <- topTags(lrt, n = nr)$table[, c(2,5,6)]
      row.names(res)<-topTags(lrt, n = nr)$table[,1]
      colnames(res) <- paste(colnames(res), '_', colnames(design)[j], 'vs', colnames(design)[i], sep = '')

      # ann <- merge(ann, round(res, 3), by.x = 1, by.y = 0)
      #ann <- merge(ann, res, by.x = 1, by.y = 0)
      #geneSyn<-as.matrix(ann[intersect(rownames(ann),rownames(res)),1])
      ann<-merge(ann,res,by.x=1,by.y=0)
      rownames(ann)<-ann[,1]
      #ann<-ann[,-1]
    }
  }
  colnames(ann)[1]<-"Genes"
  ann
}

diffExp_paired <- function(data, groups, gplabs = groups, subject, ann, ngth = 4,
                      force=FALSE # ignore design controls
) {
  # data = rawc; groups = pjsamples$condition; subject = pjsamples$paired
  # if(gplabs == NULL) {
  #   gplabs <- groups
  # }
  print("start paired analysis")
  ngth<-as.numeric(as.character(ngth))
  grpCtrl<-data.frame(groups, as.numeric(subject))
  colnames(grpCtrl)<-c("groups","subject")
  allvals<-unique(as.numeric(subject))
  for (i in 1:length(allvals)){
    findgr<-grpCtrl$groups[which(grpCtrl$subject==allvals[i])]
    if (length(unique(grpCtrl$groups))!=length(unique(findgr))&&force==FALSE){ #desing control 1
      stop("conditions are not paired in correct way")
    }
  }
  print("parameters")
  print(paste("ngth=",ngth))

  nctrls <- length(unique(groups[grep('^[C,c]', groups)]))
  nexp <- length(unique(groups[grep('^[E,e]', groups)]))

  names(gplabs) <- groups
  subject <- factor(subject)
  # make design matrix
  th <- min(table(groups)) - 1
  data <- as.matrix(data[rowSums(cpm(data) > 1) > th, ])
  nr <- dim(data)[1]
  design <- model.matrix(~subject + groups)
  rownames(design) <- colnames(data)
  colnames(design) <- gsub("groups", "", colnames(design))
  y <- edgeR::DGEList(counts=data)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateGLMRobustDisp(y,design)
  fit <- edgeR::glmFit(y,design)
  ng <- length(unique(groups))
  ns <- length(unique(subject))
  gn<-  rownames(y$counts)
  #CTR1 vs all EXP
  vv <- c(grep('^CTR', colnames(design)), grep('^EXP', colnames(design)))
  ann<-as.matrix(ann[intersect(rownames(ann),gn),1])
  for(i in vv){
    #print(i)
        lrt<- edgeR::glmLRT(fit, coef=i)
        res <- edgeR::topTags(lrt, n = nr)$table[, c(1,4,5)]

        print(i)
        print(head(res))

        colnames(res) <- paste(colnames(res), '_', colnames(design)[i], 'vsCTR', sep = '')
        rownames(res)<-gn
        #ann <- merge(ann, round(res, 3), by.x = 1, by.y = 0)
        #ann <- merge(ann, res, by.x = 1, by.y = 0)

        #geneSyn<-as.matrix(ann[intersect(rownames(ann),rownames(res)),1])
        #ann<-cbind(ann,res)
        ann<-merge(ann,res,by.x=0,by.y=0)
  }

  #es. 1CTRvs5EXP
  #colnames(ann)[grep("vs",colnames(ann))]
 #  [1] "logFC_EXP1vsCTR"  "PValue_EXP1vsCTR" "FDR_EXP1vsCTR"    "logFC_EXP2vsCTR"
 #  [5] "PValue_EXP2vsCTR" "FDR_EXP2vsCTR"    "logFC_EXP3vsCTR"  "PValue_EXP3vsCTR"
 #  [9] "FDR_EXP3vsCTR"    "logFC_EXP4vsCTR"  "PValue_EXP4vsCTR" "FDR_EXP4vsCTR"
 # [13] "logFC_EXP5vsCTR"  "PValue_EXP5vsCTR" "FDR_EXP5vsCTR"

  #write.table(ann,file="/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/internalPrecision/Rosalind_PrecisionMed_paired/stAnalysis/DEGs.csv",sep="\t",col.names = TRUE, row.names = FALSE)

  if(nctrls > 1){
     vv <- c(grep('^CTR', colnames(design)))
     for(i in vv){
      for(j in ((i+1):(ns + ng - 1))) {
        cont <- rep(0, ns+ng-2+1)
        cont[i] = -1
        cont[j] = 1
        print(cont)
        cat(colnames(design)[j], 'vs', colnames(design)[i], '\n')
        lrt<-glmLRT(fit, contrast=cont)
        res <- topTags(lrt, n = nr)$table[, c(1,4,5)]
        colnames(res) <- paste(colnames(res), '_', colnames(design)[j], 'vs', colnames(design)[i], sep = '')

        #ann <- merge(ann, res, by.x = 1, by.y = 0)
        #ann <- merge(ann, round(res, 3), by.x = 1, by.y = 0)

        #geneSyn<-as.matrix(ann[intersect(rownames(ann),rownames(res)),1])
        #ann<-cbind(ann,res)
        ann<-merge(ann,res,by.x=0,by.y=0)
      }
    }
  }

  #all exp vs all exp
  #if(nexp < ngth && nctrls > 1) {
if( ng <= ngth && nexp>1) {
  print ("EXP all vs all start")
   vv <- c(grep('^EXP', colnames(design)))
   for(i in (vv[1]:(ns + ng - 2))){
      for(j in ((i+1):(ns + ng - 1))) {
        cont <- rep(0, ns+ng-2+1)
        cont[i] = -1
        cont[j] = 1
        print(cont)
        cat(colnames(design)[j], 'vs', colnames(design)[i], '\n')
        lrt<-edgeR::glmLRT(fit, contrast=cont)
        res <- edgeR::topTags(lrt, n = nr)$table[, c(1,4,5)]
        colnames(res) <- paste(colnames(res), '_', colnames(design)[j], 'vs', colnames(design)[i], sep = '')
        #ann <- merge(ann, res, by.x = 1, by.y = 0)
        #ann <- merge(ann, round(res, 3), by.x = 1, by.y = 0)
        #ann<-cbind(ann,res)
        ann<-merge(ann,res,by.x=0,by.y=0)
      }
    }
}
  colnames(ann)[1]<-"Genes"
  return(ann)
}


pca_plot <- function(data, groups, gplabs = groups, pca_f, pcaloc = 'topright', legendpage = FALSE, monoLayout=TRUE, sampleNames=TRUE,legendSize=0.6,p1=1,p2=2){

  #pca_f --> folder and pca file name
  #data --> cpm raw counts


  #==================old one
  # names(gplabs) <- groups
  # cc <- as.factor(groups)
  # names(cc) <- gplabs
  # pca_cl <- prcomp(t(log2(data + 1)))
  # pdf(pca_f)
  #   # Add extra space to right of plot area; change clipping to figure
  #   if (monoLayout){par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, pin=c(4.4,4))}
  #   plot(pca_cl$x[, c(1,2)], pch = unclass(cc), col = unclass(cc))
  #   if (sampleNames){
  #     text(pca_cl$x[, c(1,2)], row.names(pca_cl$x), cex = 0.5, pos = 3)
  #   }else{
  #       text(pca_cl$x[, c(1,2)], cex = 0.5, pos = 3)
  #   }
  #   if (legendpage) {
  #     plot.new()
  #   }
  #   if (monoLayout){
  #   legend("topright", inset=c(-0.2,0) , legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA', cex=0.8)
  #   }else{
  #   legend(pcaloc, legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA')
  #   }
  #
  # dev.off()
  # # norm_data <- merge(ann, round(data, 2), by.x = 1, by.y = 0)
  # # return(norm_data)
  #======================

  names(gplabs) <- groups
  cc <- as.factor(groups)
  names(cc) <- gplabs
  pca_cl <- prcomp(t(log2(data + 1)))
  sumrepo<-as.data.frame(as.matrix(summary(pca_cl)$importance))
  sumrepo<-cbind(row.names(sumrepo),sumrepo)
  colnames(sumrepo)[1]<-" "
  write.table(sumrepo,file="pca_summary.txt",sep="\t",row.names=FALSE,quote=FALSE)
  pdf(pca_f)
    # Add extra space to right of plot area; change clipping to figure
    if (monoLayout){par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, pin=c(4.4,4))}
    plot(pca_cl$x[, c(p1,p2)], pch = unclass(cc), col = unclass(cc), xlab=paste0(colnames(pca_cl$x[, c(p1,p2)])[1],"(",as.character(round(sumrepo[2,which(colnames(sumrepo)==paste0("PC",p1))],2)),")"),
    ylab=paste0(colnames(pca_cl$x[, c(p1,p2)])[2],"(",as.character(round(sumrepo[2,which(colnames(sumrepo)==paste0("PC",p2))],2)),")")
  )
    if (sampleNames){
      text(pca_cl$x[, c(p1,p2)], row.names(pca_cl$x), cex = 0.5, pos = 3)
    }else{
        #text(pca_cl$x[, c(1,2)], cex = 0.5, pos = 3)
    }
    if (legendpage) {
      plot.new()
    }
    if (monoLayout){
    legend("topright", inset=c(-0.2,0) , legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA', cex=legendSize)
    }else{
    legend(pcaloc, legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA')
    }

  dev.off()


   return (pca_cl)



}


differential<-function(prjFolder,
                       prjName,
                       rawCounts=NULL, #must have only ENSEMBL id col before "by gene expression" samples col
                       designTablePath=NULL,
                       annoPath=NULL,
                       previous=FALSE, #if TRUE then the fun load rawcount table & design table saved in project stStatical folder
                       accurateFiltering=FALSE #ancora non è stato integrato
                       ){

  library(edgeR)

  #if annoPath is NULL , ENSEMBL id aren't used and all

  #Example parameters==========
  # prjFolder="/media/lucio/FastData/3.SuperGrive/Mega/bioinformatics/data/base/statFiltering/"
  # prjName="statFiltering"
  # rawCounts="/media/lucio/FastData/3.SuperGrive/Mega/bioinformatics/data/base/statFiltering/raw_counts_ECM_gel_Small_Intestinal_Organoid_Characterisation.csv"
  # designTablePath<-"/media/lucio/FastData/3.SuperGrive/Mega/bioinformatics/data/base/statFiltering/statisticsSS.csv"
  # annoPath<-"/media/lucio/FastData/3.SuperGrive/Mega/bioinformatics/data/base/statFiltering/quant_hg38cr_annotation.txt"
  # source("/mnt/groups/data/gbioinfo/difilippo/Pipelines/4.RNAseq_allTra/2.wkSpace/fun_statistical.R")
  #============================

  outpath=paste0(prjFolder,"/stAnalysis")

      options(stringsAsFactors=FALSE)



    if (dir.exists(outpath)==FALSE){
      dir.create(outpath)
    }

    if(previous==TRUE){
      rawCounts=paste0(outpath,"/rawcounts.csv")
      designTablePath=paste0(outpath,"/designTable.csv")

    }else{
      #copy file in wk folder
      file.copy(from=designTablePath,to=paste0(outpath,"/designTable.csv"),overwrite = TRUE)
      file.copy(from=rawCounts,to=paste0(outpath,"/rawcounts.csv"),overwrite = TRUE)
      #=====================
    }


  rwct<-read.table(rawCounts,sep="\t",header=TRUE,check.names=FALSE)
  rwctNoID<-rwct[,-1]
  geneID<-as.character(rwct[,1])
  clnamrwct<-colnames(rwct)
  #rm(rwct)
  rwctcl=length(rwctNoID[1,])
  rwctrw=length(rwctNoID[,1])
  rwctnum<-matrix(mapply(rwctNoID,FUN=as.numeric),ncol=rwctcl,nrow=rwctrw)
  colnames(rwctnum)<-clnamrwct[-1]
  print("raw counts imported as numeric matrix")
  print(head(rwctnum))
  # rwctnum<-edgeR::cpm(rwctnum)
  # rwctnum<-round(rwctnum,2)
  # rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  # colnames(rwctnum)<-clnamrwct

  #rwct<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/RNASeqAnalysis/Project_NDUF_mir181/stAnalysis/rawcounts.csv",sep="\t",header=TRUE)

  #designTable<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/RNASeqAnalysis/Project_NDUF_mir181/stAnalysis/statisticsSS_NDUF_mir181.csv", sep="\t",header=TRUE)
  designTable<-read.table(designTablePath, sep="\t", header=TRUE, check.names=FALSE)

  anno<-as.matrix(rwct[,1])
  rwct<-as.matrix(rwct[,2:length(rwct[1,])])

  mins <- min(table(as.character(designTable$Condition))) - 1

  remid=which(designTable$Removed=="t")
  print("removed samples")
  print(designTable$Original_Sample_Name[remid])
  if(length(remid)>0){
    removedSamples=as.character(designTable$Original_Sample_Name[remid])
    rwct<-rwct[ , -which(colnames(rwct) %in% removedSamples)]
    rwctnum<-rwctnum[,-which(colnames(rwctnum) %in% removedSamples)]
    clnamrwct <-c(clnamrwct[1],colnames(rwctnum))
    print("samples removed from raw matrix")
  }

  rwctnum<-edgeR::cpm(rwctnum)
  rwctnum<-round(rwctnum,2)
  rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  print("edgeR cpm raw normalizated")

  outCPM<-cbind(anno,rwctnum)
  colnames(outCPM)[1]<-"ENSEMBL"

  write.table(outCPM,file=paste0(outpath,"/cpm_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)

  rm(outCPM)

    if (accurateFiltering==FALSE){
        rwctnum<-rwctnum[rowSums(rwctnum > 1) > mins,]
        #subset rwct
        rownames(rwct)<-geneID
        rwct<-rwct[c(rwctnum$geneID),]
        anno<-as.matrix(rownames(rwct))
        print("anno")
        print(head(anno))
        #subset geneId
        geneID<-rownames(rwct)

        print(head(rwct))
        print(head(rwctnum))
        colnames(rwctnum)<-clnamrwct
    }





  #anno<-read.table("/home/ngsworkspace/references/RNASeq/quant_mm10_annotation.txt",header=FALSE)
  #anno<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/references/RNASeq/quant_mm10_annotation.txt",header=TRUE)

  grp<-c()
  sbj<-c()
  for(i in 1:length(rwct[1,])){
    sample<-colnames(rwct)[i]
    grp[i]<-as.character(designTable[which(designTable$Original_Sample_Name==sample),"Condition"])
    sbj[i]<-as.numeric(designTable[which(designTable$Original_Sample_Name==sample),"Paired"])
  }


  if ( is.null(annoPath)==FALSE){
      rostone<-read.delim(annoPath,header=TRUE,check.names=FALSE,sep="\t")
      geneNames<-c()
      leno<-length(anno[,1])
      for( i in 1:leno){
        print(paste0(i,"/",leno))
        findNam<-as.character(rostone[which(as.character(rostone[,1])==as.character(anno[i,1])),2])
        if(length(findNam)>0){

             geneNames[i]<-paste0(i,"_",findNam)


          #geneNames[i]<-as.character(rostone[which(as.character(rostone[,1])==as.character(anno[i,1])),2])
        }else{
          print("symbol not find in annotation file")
          geneNames[i]<-as.character(anno[i,1])
        }

      }

      if (any(is.na(geneNames))==TRUE){print("NA value in geneNames Array")}
      #"
      geneNames<-as.matrix(geneNames)
      rownames(geneNames)<-anno[,1]
      rownames(rwct)<-anno[,1]
  } else {
    #rownames(rwct)<-geneID
    geneNames<-geneID

    geneNames<-as.matrix(geneNames)

    rownames(geneNames)<-geneNames[,1]
    rownames(geneNames)<-geneNames
  }

  pca_plot(data=rwct,group=grp,pca_f=paste0(outpath,"/pca_",prjName,".pdf"))

print(paste0("unique(sbj): ",unique(sbj)))
print(paste0("sbj: ",sbj))

  if (length(unique(sbj))==1 && unique(sbj)[1]==0){
    resu<-diffExp_unpaired(data=rwct,groups=grp,ann=geneNames)
  } else {
   resu<-diffExp_paired(data=rwct,groups=grp,subject=sbj,ann=geneNames)
  }

  resu<-cbind(row.names(resu),resu)
  rwNames<-cbind(row.names(geneNames),geneNames)
  resu<-merge(rwNames,resu,by.x=2,by.y=1)

  #resu<-de_paired(data=rwct,groups=grp,subject=sbj,ann=geneNames)
  #pvcols <- grep('^PValue_', colnames(resu))
  #results[,pvcols] <- round(results[,pvcols], 5)
  #fccols <- grep('^logFC_', colnames(resu))
  #results[,fccols] <- round(results[,fccols], 3)
  #fdrcols <- grep('^FDR_', colnames(resu))

  if(is.null(annoPath)==FALSE){
    #resu<-cbind(rownames(resu),resu)

    colnames(resu)[1]<-"Annotation"
    colnames(resu)[2]<-"ENSEMBL"
    colnames(resu)[3:length(colnames(resu))]
    resu<-resu[,c("ENSEMBL","Annotation",colnames(resu)[3:length(colnames(resu))])]


    mergcol_resu=1
  } else {
    colnames(resu)[1]<-"ENSEMBL"
    resu<-resu[,-2]
    mergcol_resu=1
  }


  #"
  write.table(resu, file=paste0(outpath,"/edgeR_de_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)

  # raw counts instead cpm
  #rwct<-read.table(rawCounts,sep="\t",header=TRUE,check.names=FALSE)


  # rwctNoID<-rwct[,-1]
  # geneID<-as.character(rwct[,1])
  # clnamrwct<-colnames(rwct)
  # #rm(rwct)
  # rwctcl=length(rwctNoID[1,])
  # rwctrw=length(rwctNoID[,1])
  # rwctnum<-matrix(mapply(rwctNoID,FUN=as.numeric),ncol=rwctcl,nrow=rwctrw)
  # rwctnum<-edgeR::cpm(rwctnum)
  # rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  # colnames(rwctnum)<-clnamrwct
  print("saving fulltable")

  mergcol_rwctnum=1

  #prefinal<-merge(geneNames,rwctnum,by.x=0,by.y=1)

  #colnames(prefinal)[1]<-"ENSEMBL"
  #colnames(prefinal)[2]<-"Annotation"

  #rwct ha i ensembl mentre la resu
  fullResults<-merge(rwctnum,resu, by.x=mergcol_rwctnum,by.y=mergcol_resu,all=TRUE)
  colnames(fullResults)[1]<-"ENSEMBL"

  clnames<-colnames(fullResults)[which(!colnames(fullResults)%in%c("ENSEMBL","Annotation"))]

  fullResults<-fullResults[,c("ENSEMBL","Annotation",clnames)]

  #fullResults<-merge(rwctnum,resu, by.x="Genes",by.y="Genes",all=TRUE)


  print("rwctnum")
  print(head(rwctnum))
  print("resu")
  print(head(resu))

  print("fullResults")
  print(head(fullResults))

  colnam_rwctnum=colnames(rwctnum)[mergcol_rwctnum]
  colnam_resu=colnames(resu)[mergcol_resu]
  print(paste0("colnam_rwctnum=",colnam_rwctnum))

  print(paste0("colnam_resu=",colnam_resu))

  #setdiff(colnames(C),c("geneId","geneSyn"))
  #fullResults<-fullResults[,c("geneId","Genes",setdiff(colnames(fullResults),c("geneId","Genes")))]

  #fullResults<-fullResults[,c(colnam_rwctnum,colnam_resu,setdiff(colnames(fullResults),c(colnam_rwctnum,colnam_resu)))]

  write.table(fullResults, file=paste0(outpath,"/fullResults_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)

}


#=================
#
#=================


statRawTable<-function(prjName, rawCounts, condition, paired, outpath=getwd(), removeColumns=FALSE, HTseqMM_Filtering=NULL, filtering=TRUE, fixmin=NULL){

  #statistics from raw table

  # HTseqMM_Filtering  path multimap raw table (if isn't defined, skip hitc multimap filtering strategy)



  #Example
  #statRawTable(prjName="geoGSE110051",rawCounts="/home/users/difilippo/mywd/otherAnalysis/Rossella/GSE110051_count_table_RAW.txt.gz",condition=c("CTR1","CTR1","CTR1","EXP1","EXP1","EXP1","EXP2","EXP2","EXP2"),paired=c(0,0,0,0,0,0,0,0,0))
  # HTseqMM_Filtering="/home/users/difilippo/mywd/quantseq/Committers/Luni/multi_counts_mm_Naive2.txt.csv"

  #removeColumns apply for quantseq results for remove MM and gene id column before apply statistics

  library(edgeR)

  #outpath=paste0(prjFolder,"/stAnalysis")

  options(stringsAsFactors=FALSE)

  # if (dir.exists(outpath)==FALSE){
  #   dir.create(outpath)
  # }

  #if(previous==TRUE){
    #rawCounts=paste0(outpath,"/rawcounts.csv")
    #designTablePath=paste0(outpath,"/designTable.csv")

  #}else{
    #copy file in wk folder
    #file.copy(from=designTablePath,to=paste0(outpath,"/designTable.csv"),overwrite = TRUE)
    #file.copy(from=rawCounts,to=paste0(outpath,"/rawcounts.csv"),overwrite = TRUE)
    #=====================
  #}

  # alpha threshold for FDR
  ath <- 0.05
  #%MM reads threshold
  mmth <- 0.20


  sbj=paired
  grp=condition

  rwct<-read.table(rawCounts,sep="\t",header=TRUE,check.names=FALSE)

  if (removeColumns==TRUE){
    rwct<-rwct[,-3]
    rwct<-rwct[,-2]
  }

  rwctNoID<-rwct[,-1]
  rownames(rwctNoID)<-rwct[,1]
  geneID<-as.character(rwct[,1])
  clnamrwct<-colnames(rwct)
  #rm(rwct)



  designMatrix<-data.frame(colnames(rwct)[-1],condition,paired)
  colnames(designMatrix)<-c("sampleName","condition","paired")
  write.table(designMatrix, file=paste0(outpath,"/design_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)


  rwctcl=length(rwctNoID[1,])
  rwctrw=length(rwctNoID[,1])

  rwctnum<-matrix(mapply(rwctNoID,FUN=as.numeric),ncol=rwctcl,nrow=rwctrw)
  colnames(rwctnum)<-clnamrwct[-1]
  rownames(rwctnum)<-rownames(rwctNoID)
  print("raw counts imported as numeric matrix")
  print(head(rwctnum))

  # rwctnum<-edgeR::cpm(rwctnum)
  # rwctnum<-round(rwctnum,2)
  # rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  # colnames(rwctnum)<-clnamrwct

  #rwct<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/RNASeqAnalysis/Project_NDUF_mir181/stAnalysis/rawcounts.csv",sep="\t",header=TRUE)

  #designTable<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/RNASeqAnalysis/Project_NDUF_mir181/stAnalysis/statisticsSS_NDUF_mir181.csv", sep="\t",header=TRUE)


  # designTable<-read.table(designTablePath, sep="\t", header=TRUE)
  #
  # rostone<-read.table(annoPath,header=TRUE)
  #
  # anno<-as.matrix(rwct[,1])
  # rwct<-as.matrix(rwct[,2:length(rwct[1,])])
  #
  # remid=which(designTable$Removed=="t")
  # print("removed samples")
  # print(designTable$Original_Sample_Name[remid])
  # if(length(remid)>0){
  #   removedSamples=as.character(designTable$Original_Sample_Name[remid])
  #   rwct<-rwct[ , -which(colnames(rwct) %in% removedSamples)]
  #   rwctnum<-rwctnum[,-which(colnames(rwctnum) %in% removedSamples)]
  #   clnamrwct <-c(clnamrwct[1],colnames(rwctnum))
  #   print("samples removed from raw matrix")
  # }
  #

  #filtering  HERE

  cpmctnum<-edgeR::cpm(rwctnum)
  cpmctnum<-round(cpmctnum,2)

  if (filtering==TRUE){

    if(is.null(fixmin)==TRUE){
      mins <- min(table(as.character(condition))) - 1
      print(paste0("mins fixed to:", mins))

    } else {
      mins <- fixmin
      print(paste0("mins fixed to:", mins))
    }

      keep_id2 <- names(rowSums(rwctnum > 1) > mins)[rowSums(rwctnum > 1) > mins]

    #exctly same result of this:
    #row.names(subset(rwctnum, rowSums(rwctnum > 1) > mins ))
    cat(length(keep_id2), 'genes passing low reads filter with threshold', mmth, '\n')
    if(is.null(HTseqMM_Filtering)==FALSE){
      #TODO calculate mean percent
      mm2<-read.table(HTseqMM_Filtering,sep="\t",header=TRUE)
      #colnames(mm2)[1]<-"EnsemblGeneID"
      rownames(mm2)<-mm2[,1]
      mm2<-mm2[,-1]

      # pmcunts = % non-unique reads
      #if(length(unique(which(colnames(mm2)==colnames(cpmctnum))))) TODO column need to have same order you need make a check of this

      pmcounts<- 1 - mm2/rwctnum

      mean_perc <- rowMeans(pmcounts, na.rm = TRUE)

      #names(mean_perc)

      #MM<-cbind(names(mean_perc),mean_perc)
      #colnames(MM)<-c("Row.names","Perc_MM_reads")

      #resultMM<-merge(mean_perc,rwctnum,by="row.names")

      keep_id1 <- names(mean_perc)[mean_perc <= mmth & !is.na(mean_perc)]
      mean_perc<-as.data.frame(mean_perc)
      MM_perc<-cbind(rownames(mean_perc),mean_perc[,1])
      colnames(MM_perc)<-c("GeneName","MM_perc")

      filteredc <- cpmctnum[intersect(keep_id1, keep_id2),]

    }else{
      filteredc <- cpmctnum[keep_id1,]
    }

  Removed_gene <- as.data.frame(!(row.names(cpmctnum) %in% intersect(keep_id1, keep_id2)))
  Removed_gene<-cbind(row.names(cpmctnum),Removed_gene)
  colnames(Removed_gene)<-c("GeneName","Removed_gene")

  #cpmctnum<-filteredc
  "GeneName"
  filtab<-merge(Removed_gene, cpmctnum, by.x="GeneName", by.y=0)
  write.table(filtab, file=paste0(outpath,"/filtered_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)

} else {
    write.table(cpmctnum, file=paste0(outpath,"/cpmCounts_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)
}
  #END




  #cpmctnum<-cbind(rwct[,1],cpmctnum)
  #colnames(cpmctnum)[1]<-"GeneName"

  #rownames(rwctnum)<-rwct[,1]
  #rownames(cpmctnum)<-rwct[,1]


  #rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  print("edgeR cpm raw normalizated")
  print(head(rwct))
  print(head(rwctnum))
  print(head(cpmctnum))
  #colnames(rwctnum)<-clnamrwct[-1]



  #anno<-read.table("/home/ngsworkspace/references/RNASeq/quant_mm10_annotation.txt",header=FALSE)
  # #anno<-read.table("/home/lucio/sshfsFolder/mntNgsWorkspace/references/RNASeq/quant_mm10_annotation.txt",header=TRUE)
  #
  # grp<-c()
  # sbj<-c()
  # for(i in 1:length(rwct[1,])){
  #   sample<-colnames(rwct)[i]
  #   grp[i]<-as.character(designTable[which(designTable$Original_Sample_Name==sample),"Condition"])
  #   sbj[i]<-as.numeric(designTable[which(designTable$Original_Sample_Name==sample),"Paired"])
  # }
  #
  #
  # geneNames<-c()
  # leno<-length(anno[,1])
  # for( i in 1:leno){
  #   print(paste0(i,"/",leno))
  #   geneNames[i]<-as.character(rostone[which(as.character(rostone[,1])==as.character(anno[i,1])),2])
  # }
  #
  # if (any(is.na(geneNames))==TRUE){print("NA value in geneNames Array")}
  # #"
  # geneNames<-as.matrix(geneNames)
  # rownames(geneNames)<-anno[,1]
  # rownames(rwct)<-anno[,1]
  #
  # pca_plot(data=rwct,group=grp,pca_f=paste0(outpath,"/pca_",prjName,".pdf"))

  #rwctnum<-matrix(mapply(rwctNoID,FUN=as.numeric),ncol=rwctcl,nrow=rwctrw)
  geneID<-as.matrix(geneID)
  rownames(geneID)<-geneID[,1]
  colnames(geneID)[1]<-"GeneName"

  #data2<<-rwctnum
  #grp2<<-grp
  #ann<<-geneID

  if (length(unique(sbj))==1 && sbj==0){
    if(filtering==TRUE){
        geneID<-as.matrix(rownames(filteredc))
        rownames(geneID)<-geneID[,1]
        colnames(geneID)[1]<-"GeneName"
        resu<-diffExp_unpaired(data=filteredc,groups=grp,ann=geneID)
    }else{
        resu<-diffExp_unpaired(data=rwctnum,groups=grp,ann=geneID)
    }

  } else {
    if(filtered==TRUE){
      geneID<-as.matrix(rownames(filteredc))
      rownames(geneID)<-geneID[,1]
      colnames(geneID)[1]<-"GeneName"
      resu<-diffExp_paired(data=filteredc,groups=grp,subject=sbj,
        ann=geneID)
    }else{
      resu<-diffExp_paired(data=rwctnum,groups=grp,subject=sbj,ann=geneID)
    }

  }


  #resu<-de_paired(data=rwct,groups=grp,subject=sbj,ann=geneNames)
  #pvcols <- grep('^PValue_', colnames(resu))
  #results[,pvcols] <- round(results[,pvcols], 5)
  #fccols <- grep('^logFC_', colnames(resu))
  #results[,fccols] <- round(results[,fccols], 3)
  #fdrcols <- grep('^FDR_', colnames(resu))
  #resu<-cbind(rownames(resu),resu)
  colnames(resu)[1]<-"GeneName"
  #"
  write.table(resu, file=paste0(outpath,"/edgeR_de_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)

  #rwct<-read.table(rawCounts,sep="\t",header=TRUE)

  # rwctNoID<-rwct[,-1]
  # geneID<-as.character(rwct[,1])
  # clnamrwct<-colnames(rwct)
  # #rm(rwct)
  # rwctcl=length(rwctNoID[1,])
  # rwctrw=length(rwctNoID[,1])
  # rwctnum<-matrix(mapply(rwctNoID,FUN=as.numeric),ncol=rwctcl,nrow=rwctrw)
  # rwctnum<-edgeR::cpm(rwctnum)
  # rwctnum<-cbind(geneID,as.data.frame(rwctnum))
  # colnames(rwctnum)<-clnamrwct



  fullResults<-merge(cpmctnum,resu, by.x=0,by.y="GeneName",all=TRUE)
  colnames(fullResults)[1]<-"GeneName"

  if(filtering==TRUE){
    if(is.null(HTseqMM_Filtering)==FALSE){
        adding<-merge(MM_perc,Removed_gene,by="GeneName")
        fullResults<-merge(adding,fullResults, by="GeneName",all=TRUE)
    }else{
        fullResults<-merge(Removed_gene,fullResults, by="GeneName",all=TRUE)
    }
  }else{
    #fullResults<-merge(cpmctnum,resu, by.x=0,by.y="GeneName",all=TRUE)
  }


  #setdiff(colnames(C),c("geneId","geneSyn"))
  #fullResults<-fullResults[,c("geneId","geneSyn",setdiff(colnames(fullResults),c("geneId","geneSyn")))]
  write.table(fullResults, file=paste0(outpath,"/fullResults_",prjName,".csv"), sep = '\t', row.names = FALSE, col.names=TRUE, quote=FALSE)
}




filedgeR<-function(data,
                   ann,
                   groups
                   ){
    print("parameters")

    library(edgeR)
    #print("start analysis")

    #names(gplabs) <- groups
    # make design matrix
    th <- min(table(groups)) - 1
    #data <- as.matrix(data[rowSums(cpm(data) > 1) > th, ])
    gnNames<-rownames(data)
    design <- model.matrix(~0 + groups)
    rownames(design) <- colnames(data)
    colnames(design) <- gsub("groups", "", colnames(design))
    y <- DGEList(counts=data,genes=ann[,1])
    y <- calcNormFactors(y)
    #if(advFilter==TRUE){
        y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
        y <- estimateGLMTrendedDisp(y, design)
    #}
    y <- estimateGLMRobustDisp(y,design)   #robust tagwise dispersions
    gn<-  rownames(y$counts)
    fit <- glmFit(y,design)
    ng <- length(unique(groups))

    #advFilter
    y.lrt<-glmLRT(fit)
    reads.cpm <- cpm(y)
    unfiltered.results <- data.frame(id=y$genes, reads.cpm, y.lrt$table)
    filter <- apply(X=reads.cpm, MARGIN=1,
                    FUN=function(datain) {datain[order(rank(datain), decreasing=TRUE)[2]]})

    lowerQuantile <- mean(filter == 0)
    if (lowerQuantile < .95){
        upperQuantile <- .95
    }else{
        upperQuantile <- 1}
    theta <- seq(lowerQuantile, upperQuantile, length=50)
    filtPadj <- genefilter::filtered_p(filter=filter, test=unfiltered.results$PValue,
                                       theta=theta, method="BH")
    min.fdr <- 0.05
    numRej <- colSums(filtPadj < min.fdr, na.rm = TRUE)

    filter.quantiles <- quantile(filter, probs=theta)

    lo.fit.theta <- lowess(numRej ~ theta, f=1/5)

    if (max(numRej) <= 10) {
        j <- 1
    } else {
        residual <- if (all(numRej==0)) {
            0
        } else {
            numRej[numRej > 0] - lo.fit.theta$y[numRej > 0]
        }
        thresh <- max(lo.fit.theta$y) - sqrt(mean(residual^2))
        j <- if (any(numRej > thresh)) {
            which(numRej > thresh)[1]
        } else {
            1
        }
    }

    plot(theta, numRej, type="b", xlab="", ylab="", lwd=3,
         frame.plot=FALSE, col="black")



    filtered.results <- unfiltered.results
    filtered.results$FDR <- filtPadj[, j, drop=TRUE]

    filtered.results$de <- sign(filtered.results$logFC)*(filtered.results$FDR < 0.05)
    filtered.results$sig <- abs(filtered.results$de)

    if(length(which(is.na(filtered.results$sig)==TRUE))>0){
      filtered.results[is.na(filtered.results$sig),]$sig <- 0
    }


    ###
    ### This command does a "head" on just the significantly DE genes. This let's
    ### us, at a glance, make sure things worked as expected.
    ###
    head(filtered.results[filtered.results$sig==1,])
    filteredGenes<-filtered.results[filtered.results$sig==1,1]
    return(filteredGenes)

    #========

}
