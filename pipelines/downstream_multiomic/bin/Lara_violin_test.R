
#====================set environment

library(networkD3)
library(target)
library(GenomicRanges)

setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

source("git/downstream_multiomic/bin/hypergeometric.R")

#====================input
int_basal<-"outs/distfiltered_great/basal/20240704-public-4.0.4-HYyfkD-mm10-all-gene.txt"
int_nearest<-"outs/distfiltered_great/nearest/20240704-public-4.0.4-LTZU3b-mm10-all-gene.txt"

gene.annWithPosition<-read.table("in/ucsc/martexport.txt", sep="\t", header=TRUE)

int<-read.table(int_basal, sep="\t",header=FALSE)

degs_D4<-read.table("in/DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv", header=TRUE, sep="\t")

degs_D0<-read.table("in/DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv", header=TRUE, sep="\t")

allanno.genes<-read.table("in/GREATv4/GREATv4.genes.mm10.tsv",sep="\t",header=FALSE)

anno<-read.table("in/DEseq2_RNAseq/fixed_Mus_musculus.anno.tsv", header=TRUE, sep="\t")
norm.counts<-read.table("in/DEseq2_RNAseq/all.normalised_counts.tsv",sep="\t",header=TRUE)


degs <- degs_D4

#which(int.ensembl$V2.y=="ENSMUSG00000000078")
#==========================

violin.expression<-function(ctable, intgroups, logfc=FALSE){
  
  #colnames(tarpeaks_degs_counts) 
  # [1] "x"               "D4WT_REP1"       "D4WT_REP2"       "D4WT_REP3"      
  # [5] "D4MLL1KO_REP1"   "D4MLL1KO_REP2"   "D4MLL1KO_REP3"   "D4MLL2KO_REP1"  
  # [9] "D4MLL2KO_REP2"   "D4MLL2KO_REP3"   "D4DoubleKO_REP1" "D4DoubleKO_REP2"
  # [13] "D4DoubleKO_REP3" "D0DoubleKO_REP1" "D0DoubleKO_REP2" "D0DoubleKO_REP3"
  # [17] "D0Mll1KO_REP1"   "D0Mll1KO_REP2"   "D0Mll1KO_REP3"   "D0Mll2KO_REP1"  
  # [21] "D0Mll2KO_REP2"   "D0Mll2KO_REP3"   "D0WTA_REP1"      "D0WTA_REP2"     
  # [25] "D0WTA_REP3" 
  #note : x = ensembl
  
  #intgroups<-c("D0DoubleKO","D0WTA")
  
  tarpeaks_degs_counts<-ctable
  
  id1<-grep(intgroups[1],colnames(tarpeaks_degs_counts))
  id2<-grep(intgroups[2],colnames(tarpeaks_degs_counts))
  
  id1.tarpeaks_degs_counts<-tarpeaks_degs_counts[,c(1,id1)]
  id2.tarpeaks_degs_counts<-tarpeaks_degs_counts[,c(1,id2)]
  
  #signgenes.down<-degs[which(degs$padj<=0.05 & degs$log2FoldChange<=-1),1]
  #signgenes.up<-degs[which(degs$padj<=0.05 & degs$log2FoldChange>=1),1]
  #signgenes<-c(signgenes.down,signgenes.up)
  
  #id1.tarpeaks_degs_counts<-merge(signgenes, id1.tarpeaks_degs_counts, by=1)
  #id2.tarpeaks_degs_counts<-merge(signgenes, id2.tarpeaks_degs_counts, by=1)
  
  row.names(id1.tarpeaks_degs_counts)<-id1.tarpeaks_degs_counts[,1]
  id1.tarpeaks_degs_counts<-id1.tarpeaks_degs_counts[,-1]
  
  row.names(id2.tarpeaks_degs_counts)<-id2.tarpeaks_degs_counts[,1]
  id2.tarpeaks_degs_counts<-id2.tarpeaks_degs_counts[,-1]
  
  id1.tarpeaks_degs_means<-as.data.frame(rowMeans(id1.tarpeaks_degs_counts))
  id2.tarpeaks_degs_means<-as.data.frame(rowMeans(id2.tarpeaks_degs_counts))
  
  id1.tarpeaks_degs_means<-cbind(ensembl=row.names(id1.tarpeaks_degs_means),
                                 id1.tarpeaks_degs_means,rep(intgroups[1]))
  id2.tarpeaks_degs_means<-cbind(ensembl=row.names(id2.tarpeaks_degs_means),
                                 id2.tarpeaks_degs_means,rep(intgroups[2]))
  
  colnames(id1.tarpeaks_degs_means)<-c("ensembl","meanExpr","condition")
  colnames(id2.tarpeaks_degs_means)<-c("ensembl","meanExpr","condition")
  
  row.names(id1.tarpeaks_degs_means)<-NULL
  row.names(id2.tarpeaks_degs_means)<-NULL
  
  
  tarpeaks_degs_exprdf<-as.data.frame(rbind(id1.tarpeaks_degs_means,id2.tarpeaks_degs_means))

  if(logfc==FALSE){
    return(tarpeaks_degs_exprdf)
  }else{
    logfc<-log2(rowMeans(id1.tarpeaks_degs_counts+1)) - log2(rowMeans(id2.tarpeaks_degs_counts+1))

    return(logfc)
  }  
  

  
}

violin.analysis<-function(annotation, degs, counts,identified_genes){
  
  #annotation=anno
  #degs=degs_D4
  #counts=
  
  int<-identified_genes
  
  nw<-as.data.frame(int)
  
  anno<-annotation
  degs_gene<-degs
  norm.counts<-counts
  
  genes.anno<-cbind(ensembl=anno$gene_id, gene=anno$gene_name)
  genes.anno<-genes.anno[which(!duplicated(genes.anno)),]
  
  degs_gene<-merge(genes.anno,degs, by=1)
  
  colnames(degs_gene)[1:2]<-c("ensembl","genename")
  
  tarpeaks_degs<-merge(nw,degs_gene, by.x=1, by.y=2)
  
  print(paste0("tarpeaks_degs:",nrow(tarpeaks_degs)))
  
  tarpeaks_degs_counts<-merge(tarpeaks_degs[,"ensembl"],norm.counts, by=1)
  
  colnames(tarpeaks_degs_counts)[1]<-"ensembl"
  
  tarpeaks_degs_exprdf<-violin.expression(tarpeaks_degs_counts, intgroups = c("D4DoubleKO","D4WT"))
  
  #control pairing
  head(tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4DoubleKO"),])
  head(tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition!="D4DoubleKO"),])
  
  library(ggplot2)
  library(dplyr)
  library(hrbrthemes)
  library(viridis)
  
  sample_size = tarpeaks_degs_exprdf %>% group_by(condition) %>% summarize(num=n())
  
  tarpeaks_degs_exprdf %>%
    left_join(sample_size) %>%
    mutate(myaxis = paste0(condition, "\n", "n=", num)) %>%
    ggplot( aes(x=myaxis, y=log10(meanExpr+1), fill=condition)) +
    geom_violin(width=0.8) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    scale_fill_viridis(discrete = TRUE) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Expression level of genes identified by GREAT") +
    xlab("")
  
  x<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4DoubleKO"),"meanExpr"]
  y<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4WT"),"meanExpr"]
  
  #t.test(x, y, alternative = "two.sided", var.equal = FALSE)
  
  wilcox.test(tarpeaks_degs_exprdf$meanExpr ~ tarpeaks_degs_exprdf$condition)
  
  wilcox.test(x, y, paired=TRUE)
  
  
}


#==========================
#compare degs

degs_D0_sign<- degs_D0[which(degs_D0$padj <= 0.05),]
degs_D0_sign<-degs_D0_sign[which(degs_D0_sign$log2FoldChange <= -1 | 
                                   degs_D0_sign$log2FoldChange >= 1),]

degs_D4_sign<- degs_D4[which(degs_D4$padj <= 0.05 ),]
degs_D4_sign<-degs_D4_sign[which(degs_D4_sign$log2FoldChange <= -1 | 
                                   degs_D4_sign$log2FoldChange >= 1),]


ensembl.d0d4.intersect<-intersect(degs_D4_sign$gene_id, 
                                  degs_D0_sign$gene_id)

# Generate 3 sets of 200 words

# Chart

VennDiagram::venn.diagram(
  x = list(degs_D0_sign$gene_id, degs_D4_sign$gene_id),
  category.names = c("D0" , "D4"),
  filename = './outs/D0_D4_venn_diagramm.png',
  output=TRUE
)



#==========================

print(paste0("int:",nrow(int)))

#int<-read.table("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/outs/great/20240703-public-4.0.4-vOagjL-mm10-all-region.txt",sep="\t",header=FALSE)

nw<-as.data.frame(int)

#nw[,2]<-gsub("\\s*\\([+-]?\\d+\\)","",nw[,2])


#p <- simpleNetwork(nw, height="100px", width="100px")


#plot(density(table(nw[,2])))

#unique(nw$V2)

genes.anno<-cbind(ensembl=anno$gene_id, gene=anno$gene_name)
genes.anno<-genes.anno[which(!duplicated(genes.anno)),]

degs_gene<-merge(genes.anno,degs, by=1)

colnames(degs_gene)[1:2]<-c("ensembl","genename")
#nw$V2 %in% degs_gene$V2

tarpeaks_degs<-merge(nw,degs_gene, by.x=1, by.y=2)

print(paste0("tarpeaks_degs:",nrow(tarpeaks_degs)))

#tarpeaks_degs_su<-cbind(gene=tarpeaks_degs$V2, ensembl=tarpeaks_degs$V1.y, padj=tarpeaks_degs$padj)
#tarpeaks_degs_counts<-merge(tarpeaks_degs_su[,"ensembl"],norm.counts, by=1)

tarpeaks_degs_counts<-merge(tarpeaks_degs[,"ensembl"],norm.counts, by=1)

colnames(tarpeaks_degs_counts)[1]<-"ensembl"

tarpeaks_degs_exprdf<-violin.expression(tarpeaks_degs_counts, intgroups = c("D4DoubleKO","D4WT"))

#control pairing
head(tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4DoubleKO"),])
head(tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition!="D4DoubleKO"),])

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)


#tarpeaks_degs_exprdf<-tarpeaks_degs_exprdf[which(signgenes %in% row.names(tarpeaks_degs_exprdf)),]
#tarpeaks_degs_exprdf<-tarpeaks_degs_exprdf[which(row.names(tarpeaks_degs_exprdf) %in% signgenes),]


sample_size = tarpeaks_degs_exprdf %>% group_by(condition) %>% summarize(num=n())

tarpeaks_degs_exprdf %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(condition, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=log10(meanExpr+1), fill=condition)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Expression level of genes identified by GREAT") +
  xlab("")


#tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$ensembl==tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$meanExpr>=10000),"ensembl"]),]

#test distribution

#install.packages("ggpubr")

#x<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D0DoubleKO"),"meanExpr"]
#y<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D0WTA"),"meanExpr"]

x<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4DoubleKO"),"meanExpr"]
y<-tarpeaks_degs_exprdf[which(tarpeaks_degs_exprdf$condition=="D4WT"),"meanExpr"]

#t.test(x, y, alternative = "two.sided", var.equal = FALSE)

wilcox.test(tarpeaks_degs_exprdf$meanExpr ~ tarpeaks_degs_exprdf$condition)

wilcox.test(x, y, paired=TRUE)

#TRY to see the distribution on entire dataset

allgenes.counts<-merge(allanno.genes[,1],norm.counts,by=1)

allgenes.df<-violin.expression(ctable=allgenes.counts, intgroups=c("D4WT","D4DoubleKO"))

sample_size = allgenes.df %>% group_by(condition) %>% summarize(num=n())

allgenes.df %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(condition, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=log10(meanExpr+1), fill=condition)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Expression level of all genes") +
  xlab("")

x<-allgenes.df[which(allgenes.df$condition=="D0DoubleKO"),"meanExpr"]
y<-allgenes.df[which(allgenes.df$condition=="D0WTA"),"meanExpr"]

wilcox.test(allgenes.df$meanExpr ~ allgenes.df$condition)

wilcox.test(x, y, paired=TRUE)

# print(paste0("wilcox.allgenes$p.value>wilcox.tarpeaks$p.value ",
#        wilcox.allgenes$p.value>wilcox.tarpeaks$p.value       
#        ))
# print(paste0("wilcox.allgenes$p.value is ",
#        wilcox.allgenes$p.value/wilcox.tarpeaks$p.value,
#        " more then wilcox.tarpeaks$p.value"
# ))

#box plot with logfc file

tarpeaks_degs<-merge(nw,degs_gene, by.x=1, by.y=2)
tarpeaks_degs<-cbind(type=rep("all",nrow(tarpeaks_degs)),tarpeaks_degs)

#sample_size = allgenes.df %>% group_by(condition) %>% summarize(num=n())

allgenes.fc<-violin.expression(ctable=allgenes.counts, intgroups=c("D4DoubleKO","D4WT"),logfc = TRUE)

tarpeaks.fc<-violin.expression(ctable=tarpeaks_degs_counts, intgroups=c("D4DoubleKO","D4WT"),logfc = TRUE)

allgenes.fc.df<-data.frame(type="all_genes",log2FoldChange=allgenes.fc)

tarpeaks.fc.df<-data.frame(type="intpeaks",log2FoldChange=tarpeaks.fc)

fc.df<-rbind(allgenes.fc.df, tarpeaks.fc.df)


allgenes.logfc<-ggplot(tarpeaks_degs, aes(x=type, y=log2FoldChange)) + 
  geom_violin()

sample_size = fc.df %>% group_by(type) %>% summarize(num=n())

fc.df %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(type, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=log2FoldChange, fill=type)) +
  geom_violin(width=0.8) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Expression level of all genes") +
  xlab("")

x<-fc.df[which(fc.df$type=="all_genes"),"log2FoldChange"]
y<-fc.df[which(fc.df$type=="intpeaks"),"log2FoldChange"]

x<-x[which(!is.na(x))]
x<-x[which(!is.infinite(x))] 

y<-y[which(!is.na(y))]
y<-y[which(!is.infinite(y))]

t.test(x, y, var.equal = FALSE)

#==============================================================================
#distance by near peak a degs down / up / not with or without distances filter
#==============================================================================







#==============================================================
#explore significant gene results
#==============================================================


#all genes thet are inside great ref are mantained , the genes of rnaseq ref that are not presente inside great are removed from rnaseq results degs 
#allanno.genes$V1 norm.counts$gene_id

pres_rnaseq<-setdiff( norm.counts$gene_id, allanno.genes$V1 )

#which(allanno.genes$V1=="ENSMUSG00000074357")
#which(norm.counts$gene_id=="ENSMUSG00000074357")

int.ensembl <- merge(int,cbind(allanno.genes$V5,allanno.genes$V1),by=1)



degs_D4_sign<-degs_D4[which( degs_D4$padj<=0.05 ),]

degs_D4_sign.down<-degs_D4[which( degs_D4$padj<=0.05 & degs_D4$log2FoldChange <= -0.5 ),]
degs_D4_sign.up<-degs_D4[which( degs_D4$padj<=0.05 & degs_D4$log2FoldChange >=0.5 ),]

sign_int_d4 <- degs_D4_sign[which( degs_D4_sign$gene_id %in% int.ensembl$V2.y ),]

sign_int_d4.up <- degs_D4_sign[which( degs_D4_sign.up$gene_id %in% int.ensembl$V2.y ),]
sign_int_d4.down <- degs_D4_sign[which( degs_D4_sign.down$gene_id %in% int.ensembl$V2.y ),]

#q = length(overlap) - 1 overla:all differential expressed gene UP or DOWN present also in great_genes ; vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
#m = length(great_genes) the number of white balls in the urn.
#n = length(universe) - length(great_genes)  universe : all genes in the dataset the number of black balls in the urn.
#k = length(DEGs) degs:all differential expressed gene UP or DOWN the number of balls drawn from the urn, hence must be in 
#lower.tail = FALSE probability, it must be between 0 and 1.
#phyper number of observations. If length(nn) > 1, the length is taken to be the number required.

length(sign_int_d4.down$gene_id)
fil.sign_int_d4.down<-sign_int_d4.down[which(!sign_int_d4.down$gene_id %in% pres_rnaseq),]
fil.norm.counts<-norm.counts[which(!norm.counts$gene_id %in% pres_rnaseq),]
fil.degs_D4_sign.down<-degs_D4_sign.down[which(!degs_D4_sign.down$gene_id %in% pres_rnaseq),]

phyper(
  q = nrow(sign_int_d4.down) - 1,
  m = nrow(int),
  n = nrow(allanno.genes) - nrow(int),
  k = nrow(fil.degs_D4_sign.down),
  lower.tail = FALSE
)

d4_interesting_genes <-merge(allanno.genes, sign_int_d4.down, by=1)

paste(d4_interesting_genes$V5, collapse = " ")

#==========================================
#intersect with common genes of d0 and d4
#==========================================

degs_D0_sign<-degs_D0[which( degs_D0$padj <= 0.05 ),]

degs_D0_sign.down<-degs_D0[which( degs_D0$padj <= 0.05 & degs_D0$log2FoldChange <= -0.5 ),]
degs_D0_sign.up<-degs_D0[which( degs_D0$padj <= 0.05 & degs_D0$log2FoldChange >= 0.5 ),]

sign_int_d0 <- degs_D0_sign[which( degs_D0_sign$gene_id %in% int.ensembl$V2.y ),]

sign_int_d0.up <- degs_D0_sign[which( degs_D0_sign.up$gene_id %in% int.ensembl$V2.y ),]
sign_int_d0.down <- degs_D0_sign[which( degs_D0_sign.down$gene_id %in% int.ensembl$V2.y ),]

#d0d4.down<-intersect( sign_int_d0.down$gene_id , sign_int_d4.down$gene_id)

sign_int_d4d0.down<-merge(sign_int_d0.down, sign_int_d4.down, by=1)

length(sign_int_d4d0.down$gene_id)
fil.sign_int_d4d0.down<-sign_int_d4d0.down[which(!sign_int_d4.down$gene_id %in% pres_rnaseq),]

#fil.norm.counts<-norm.counts[which(!norm.counts$gene_id %in% pres_rnaseq),]
degs_D0D4_sign.down<-merge(degs_D0_sign.down, degs_D4_sign.down, by = 1 )
fil.degs_D0D4_sign.down<-degs_D0D4_sign.down[which(!degs_D0D4_sign.down$gene_id %in% pres_rnaseq),]

phyper(
  q = nrow(sign_int_d4d0.down) - 1,
  m = nrow(int),
  n = nrow(allanno.genes) - nrow(int),
  k = nrow(fil.degs_D0D4_sign.down),
  lower.tail = FALSE
)

#day 0 hyper

fil.degs_D0_sign.down<-degs_D0_sign.down[which(!degs_D0_sign.down$gene_id %in% pres_rnaseq),]
phyper(
  q = nrow(sign_int_d0.down) - 1,
  m = nrow(int),
  n = nrow(allanno.genes) - nrow(int),
  k = nrow(fil.degs_D0_sign.down),
  lower.tail = FALSE
)



#test target

library(target)
library(GenomicRanges)

data("real_peaks")
data("real_transcripts")

data("sim_peaks")
data("sim_transcripts")

summary(real_transcripts$adj.P.Val)

dt <- direct_targets(real_peaks, real_transcripts, 'ID', 't')

# with my data

head(genes.anno)


degs_D0_anno<-merge(gene.annWithPosition, degs_D0, by=1)
colnames(degs_D0_anno)<-c("ensembl","version",
                          "gene_name","scaffold",
                          "start", "end",
                          "strand", "baseMean",
                          "log2FC", "lfcSE",
                          "pvalue", "padj"
                          )

df <- data.frame(chr="chr1", start=11:15, end=12:16,
                 strand=c("+","-","+","*","."), score=1:5)
df
makeGRangesFromDataFrame(df)

degs_D0_anno$strand[which(degs_D0_anno$strand=="-1")]<-"-"
degs_D0_anno$strand[which(degs_D0_anno$strand=="1")]<-"+"

degs_D0_anno_gr<-makeGRangesFromDataFrame(degs_D0_anno,
                        keep.extra.columns=TRUE,
                         start.field="start",
                         seqinfo=degs_D0_anno$ensembl,
                         end.field="end",
                         strand.field = "Strand",
                         seqnames.field = "ensembl"
                         )

# gr <- GRanges(
#   seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#   ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
#   strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#   score = 1:10,
#   GC = seq(1, 0, length=10))



