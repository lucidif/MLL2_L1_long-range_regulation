#setwd("outs/Lara_multiomic_analysis/")

D0degs<-read.table("in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv",sep="\t", header=TRUE)

annotations<-read.table("/mnt/datawk1/references/annotations/ncbi/mm10_Build38/mm10.refGene.gtf.gz", sep="\t",header=FALSE)

broad_narrow_dist<-read.table("./outs/H3K27me3_broad_narrow/downGenes_peaksType_distribution.tsv",sep="\t",header=TRUE)

window<-5000

tranno<-annotations[which(annotations$V3=="transcript"),]

gene_id<-c()
#chr<-c()
#transcript_start<-c()
#transcript_end<-c()
for(i in 1:nrow(tranno)){
  print(i)
  attr <- strsplit (tranno[i,9], "; ")
  attr <- gsub (";","", attr[[1]])
  gn <- attr[grep("gene_id", attr)]
  gn <- gsub( "gene_id ", "", gn )
  gene_id[i]<-gn
}

tranno<-cbind(gene_id=gene_id,tranno)

tss<-c()
for(i in 1:nrow(tranno)){
  if(tranno[i,"V7"]=="+"){
    tss[i]<-tranno[i,"V4"]
  }else{
    if(tranno[i,"V7"]=="-"){
      tss[i]<-tranno[i,"V5"]
    }else{ #no orientation
      tss[i]<-NA
    }
  }
}

tranno<-cbind(tranno,tss)

bed<-cbind(chr=tranno[,"V1"],
           start=as.numeric(tranno[,"tss"]) - window,
           end=as.numeric(tranno[,"tss"]) + window,
           gene=tranno[,"gene_id"],
           score=tranno[,"V5"]-tranno[,"V4"],
           strand=tranno[,"V7"]
)


down<-D0degs[which(D0degs$padj<=0.05 & D0degs$log2FoldChange <= -1 ),]


subset.bed <- merge (bed, down, by.x=4, by.y=1)

subed<-subset.bed[,c(2,3,4,1,5,6)]
subed<-subed[which(duplicated(subed)==FALSE),]

write.table(subed,
            file="outs/downGenesHeatmap/D0dkoVSwt.downgenes.bed",
            col.names = FALSE,
            row.names = FALSE,
            sep="\t",
            quote = FALSE
)

broad_genes<-broad_narrow_dist[which(broad_narrow_dist$category=="broad"),"gene"]
broad_subed<-subed[which(subed$gene %in% broad_genes),]
broad_subed<-broad_subed[!duplicated(broad_subed$gene),]

narrow_genes<-broad_narrow_dist[which(broad_narrow_dist$category=="narrow"),"gene"]
narrow_subed<-subed[which(subed$gene %in% narrow_genes),]
narrow_subed<-narrow_subed[!duplicated(narrow_subed$gene),]

negative_genes<-broad_narrow_dist[which(broad_narrow_dist$category=="negative"),"gene"]
negative_subed<-subed[which(subed$gene %in% negative_genes),]
negative_subed<-negative_subed[!duplicated(negative_subed$gene),]

#TODO cosa fare con i geni che si overlappano ? (+ TSS  associati allo stesso gene)

write.table(broad_subed,
            file="outs/downGenesHeatmap/D0dkoVSwt.broad.downgenes.bed",
            col.names = FALSE,
            row.names = FALSE,
            sep="\t",
            quote = FALSE
)

write.table(narrow_subed,
            file="outs/downGenesHeatmap/D0dkoVSwt.narrow.downgenes.bed",
            col.names = FALSE,
            row.names = FALSE,
            sep="\t",
            quote = FALSE
)

write.table(negative_subed,
            file="outs/downGenesHeatmap/D0dkoVSwt.negative.downgenes.bed",
            col.names = FALSE,
            row.names = FALSE,
            sep="\t",
            quote = FALSE
)
