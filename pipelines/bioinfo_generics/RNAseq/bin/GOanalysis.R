goTerms<-function(compose,dbfile="/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/Imparato/mart_export.txt",geneLengthFile="/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/Imparato/maxlength_byGene.txt", termReduction=FALSE, filterBy="padj"){

  #dbfile<-"/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/Imparato/mart_export.txt"

  library(goseq)
  library(rrvgo)

  db<-read.delim(dbfile,sep="\t",header=TRUE)
  colnames(db)<-c("Ensembl.Gene.ID","GO.Term.Accession")
  gl<-read.delim(geneLengthFile,sep="\t",header=TRUE)

  selgenes<-intersect(names(compose),gl$Ensembl.Gene.ID)

  #names(compose)[which(!names(compose)%in%gl$Ensembl.Gene.ID)]

  #grep("ENSG00000257999",gl$Ensembl.Gene.ID)

  rownames(gl)<-gl$Ensembl.Gene.ID

  gl<-gl[selgenes,]
  gl<-gl[,-1]

  names(gl)<-selgenes

  compose<-compose[selgenes]



  #Esempio edgeR
  #tested=exactTest(disp) tabella (cols: ENSEMBL logFC logCPM PValue FDR)
  #genes=as.integer(p.adjust(tested$table$PValue[tested$table$logFC!=0],
  #+ method="BH")<.05)
  #names(genes)=row.names(tested$table[tested$table$logFC!=0,])
  #> table(genes)
  #genes
  #  0     1
  #19535 3208

  #1 extract gene length from bash:
  #remove any comment lines from the GTF file and use awk to convert column 1 to "1"
  #grep -v '#' $HSA_GTF | awk '{OFS="\t"} $1=1' >> $HSA_GTF.tmp\

  #now run GTF tools
  #gtftools.py  -l Homo_sapiens.GRCh38.90.gtf.genelength  Homo_sapiens.GRCh38.90.gtf.tmp



  #Determining Genome and Gene ID

  #supportedOrganisms()[grep("h38",supportedOrganisms()$Genome),]
  #supportedOrganisms()[supportedOrganisms()$Genome=="hg38",]

  #doesn't supported

  #The nullp will fit a model to account for gene length biases in our data. In order to use the human genome version hg38 we need the TxDb.Hsapiens.UCSC.hg38.knownGene package to be installed.

  #BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

  #GO analysis

  #1  Fitting the Probability Weighting Function (PWF)

  par(mar=c(1,1,1,1))

  pwf <- nullp(compose,bias.data=gl)

  #gene set enrichment analysis with the goseq function.

  goterms<-goseq(pwf,gene2cat = db) #use_genes_without_cat=TRUE

  goterms.allinfo<-cbind(goterms,p.adjust(goterms$over_represented_pvalue,method = "BH"))

  colnames(goterms.allinfo)[ncol(goterms.allinfo)]<-"padj"

  #length(which(goterms$padj<0.05))

  #goterms.randomaSampling<-goseq(pwf, "hg38","geneSymbol",use_genes_without_cat=TRUE,method = "Sampling",repcnt=1000)

  #Biological Process
  #goterms.bp <- goseq(pwf, "hg38","geneSymbol",test.cats="GO:BP")

  # Over-represented means that there are more DE genes in the
  # category than we would expect given the size of the category and the
  # gene length distribution so that would be enriched for DE genes.
  # Under-represented means that there are fewer DE genes in the category
  # than we would expect by chance. The p-value relates to the probability
  # of observing this number of DE genes in the category by chance.

  #Molecular Function

  #goterms.mf <- goseq(pwf, "hg38","geneSymbol",test.cats="GO:MF")

  #Cell Component
  #goterms.cc <- goseq(pwf, "hg38","geneSymbol",test.cats="GO:CC")

  #goterms.allinfo<-goseq(pwf, "hg38","geneSymbol",use_genes_without_cat=TRUE)


  #goterms.allinfo<-rbind(goterms.bp,goterms.mf,goterms.cc)

  goterms.sign<-goterms.allinfo[which(goterms.allinfo[,filterBy]<=0.05),]

  goterms.sign<-goterms.sign[order(goterms.sign$over_represented_pvalue,by=goterms.sign$padj),]

  goterm.final<-goterms.sign[,c("category", "padj","numDEInCat","numInCat","term","ontology")]

  goterm.final<-cbind(goterm.final,goterm.final$numDEInCat/goterm.final$numInCat)

  colnames(goterm.final)[ncol(goterm.final)]<-"ratio"

  if(termReduction==TRUE){
    simMatrix <- calculateSimMatrix(goterm.final$category,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")

    reducedTerms <- reduceSimMatrix(simMatrix,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")

    reducedTerms_2<- cbind(reducedTerms$go,reducedTerms$cluster)

    goterm.results<-merge(reducedTerms_2,goterm.final,by=1)
  }else{
    goterm.results<-goterm.final
  }


  if(length(which(!is.na(goterm.results$term)))>1){
    goterm.results<-goterm.results[which(!is.na(goterm.results$term)),]
  }




  return(goterm.results)



}

godot<-function(goterms,fileout=NULL,type="barplot", w=1200, h=1000,plotTitle="",xtitle="",ytitle=""){

  library(ggplot2)

  #type: barplot ; dotplot
  #fileout="/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/Imparato/GOenrichmentAnalysis/test.png"
  #goterms=padj.terms

  #if(!is.null(colorCat)){
    #customPal<-c("red","blue","green")

  #}




  if(type=="barplot"){
    #png(fileout, width = w, height = h, res = 300)

    if(xtitle==""){
      xtitle="GOterms"
    }

    if(ytitle==""){
      ytitle="-log(padj)"
    }

#-log(padj)
    p<-ggplot(data=goterms, aes(x= reorder(term, -padj) , y=-log(padj),fill=ontology)) +
      geom_bar(stat="identity", position=position_dodge()) +
      facet_grid(ontology ~., scale="free" )+
      ggtitle(plotTitle)+
      xlab(xtitle)+ #il plot è flippato
      ylab(ytitle)+
      coord_flip()+
      scale_fill_manual(values=c("BP"="firebrick","CC"="deepskyblue1","MF"="aquamarine4"))+
      geom_text(data=goterms, aes(x=reorder(term, -padj), y=-log(padj)/2, label=round(ratio,2)), col='black', size=5)
    #dev.off()
  }

  if(type=="dotplot"){


    if(xtitle==""){
      xtitle="-log(padj)"
    }

    if(ytitle==""){
      ytitle="GOterms"
    }

    #png(fileout, width = w, height = h, res = 300)
    p<-ggplot(goterms, aes(x=-log(padj), y=term, group=ontology )) +
      facet_grid(ontology ~., scale="free" )+
      geom_point(mapping = aes(size=numDEInCat,color=ratio)

                 #,aes(color=ratio)
      )+
      scale_size(name   = "numInCat",
                 breaks = fivenum(goterms$numDEInCat),
                 labels = fivenum(goterms$numDEInCat))+
      scale_color_gradientn(colours = rainbow(5))+
      xlab(xtitle)+
      ylab(ytitle)+
      ggtitle(plotTitle)
    #dev.off()
  }

  if(!is.null(fileout)){
    ggsave(p,file=fileout,width = w , height = h, units = "px", device = "png")
  }



  return(p)


}

inTermDEGs<-function(compose, goTerms.results, mart_export.file="/home/lucio/MEGA/bioinformatics/data/QuantSeq/2021/commitments/Imparato/mart_export.txt",gensyn=TRUE,database="hsapiens_gene_ensembl"){


  martexp<-read.table(mart_export.file,sep="\t",header=TRUE)

  mart.tar<-martexp[which(martexp$Gene.stable.ID%in%names(which(compose>0))),]

  if(gensyn==TRUE){
    library('biomaRt')
    mart <- useDataset(database, useMart("ensembl"))
    #df$id <- NA
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                              "entrezgene_id",
                                                              "external_gene_name",
                                                              "description"),
                    values=mart.tar$Gene.stable.ID,mart= mart)

    ensembl<-c()
    gene<-c()
    for(i in 1:length(G_list[,1])){

      ensembl[i]<-G_list$ensembl_gene_id[i]
      gene[i]<-G_list$external_gene_name[i]

    }
    targen<-data.frame(ensembl,gene)
  }else{
    targen<-data.frame(names(compose),names(compose)) # è un accrocchio
    colnames(targen)<-c("ensembl","gene")
  }




  termGenes<-c()
  for(j in 1:nrow(goTerms.results)){
    cat<-goTerms.results$category[j]
    selmart<-mart.tar[which(mart.tar$GO.term.accession==cat),]
    selmart<-selmart[!duplicated(selmart),]
    allg<-merge(targen,selmart,by=1)
    termGenes[j]<-paste(allg[!duplicated(allg),"gene"],collapse=",")
  }


  out<-cbind(goTerms.results,termGenes)

  return(out)

}
