gitwk="C:/gitwk/bioinfoGenerals/RNAseq/bin"
source(paste0(gitwk, "/fun_statistical.R"))
wkdir="C:/wkdir/Lara/RNAseq"

#unite unique counts
rawfolder="C:/wkdir/Lara/RNAseq/star_sem"
files<-list.files(rawfolder, pattern=".genes")


for (i in 1:length(files)){
  intab<-read.table(paste0(rawfolder,"/",files[i]),header=TRUE)
  midtab<-cbind(intab$gene_id,intab$expected_count)
  colnames(midtab)<-c("gene_id",files[i])
  
  if(i ==1 ){
    outtab <- midtab
  }else{
    outtab <- merge(outtab, midtab, by=1)
  }
  
}

colnames(outtab)<-gsub(".genes.results","",colnames(outtab))

write.table(outtab, file="C:/wkdir/Lara/RNAseq/star_sem/rawcounts.tsv.csv", 
            sep="\t" , col.names=TRUE, row.names=FALSE, quote=FALSE)

# colnames(outtab) <- gsub(".genes.results", "",colnames(outtab))
  annotations.path <- "C:/wkdir/Lara/RNAseq/gencode.vM10.annotation.gtf.gz"
# 
   annotations<-read.table(annotations.path, sep="\t", header=FALSE)
# 
# annotations<-read.table("C:/wkdir/Lara/RNAseq/gencode.vM10.annotation.gtf.gz", sep="\t", header=FALSE)
# 
  annotations<-annotations[which(annotations[,3]=="transcript"),]
# 
  gene.id<-grep("gene_id", unlist(strsplit(annotations[,9], "; ")),value=TRUE)
  gene.name<-grep("gene_name", unlist(strsplit(annotations[,9], "; ")),value=TRUE)
# 
  anno.2<-data.frame(gene.id,gene.name)
  #colnames(anno.2)<-c("gene_id","gene_name")
  anno.2$gene.id<-gsub("gene_id ", "", anno.2$gene.id)
  anno.2$gene.id<-gsub("\\..*","",anno.2$gene.id)
  anno.2$gene.name<-gsub("gene_name ", "", anno.2$gene.name)
  anno.3<-anno.2[which(duplicated(anno.2)==FALSE),]
  
  length(unique(anno.3$gene.id))==length(anno.3$gene.id)
  length(unique(anno.3$gene.name))==length(anno.3$gene.name)
  anno.3$gene.name[which(duplicated(anno.3$gene.name)==TRUE)]<-
    paste0("dupgene",which(duplicated(anno.3$gene.name)==TRUE),"_",
         anno.3$gene.name[which(duplicated(anno.3$gene.name)==TRUE)])

  write.table(anno.3, file="C:/wkdir/Lara/RNAseq/gencode_GRCm38_p4_annotations.tsv", row.names = FALSE, col.names = TRUE, quote=FALSE, sep="\t")