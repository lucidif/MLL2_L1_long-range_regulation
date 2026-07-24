#load libraries

library(GenomicFeatures)
library(edgeR)
library(ggplot2)

#TODO transform fpkm in function

#parameters

setwd("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4")
annotations<-"/mnt/datawk1/references/annotations/ncbi/mm10_Build38/mm10.refGene.gtf.gz"

# load data

D0.DKOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D0_DKO_vs_WT.deseq2.results.tsv", sep="\t", header = TRUE)
D0.MLL1KOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D0Mll1KO_vs_WT.deseq2.results.tsv", sep="\t", header = TRUE)
D0.MLL2KOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D0Mll2KO_vs_WT.deseq2.results.tsv", sep="\t", header = TRUE)

D4.DKOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D4_DKO_VS_WT.deseq2.results.tsv", sep="\t", header = TRUE)
D4.MLL1KOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D4MLL1KO_VS_WT.deseq2.results.tsv", sep="\t", header = TRUE)
D4.MLL2KOvsWT<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D4MLL2KO_vs_WT.deseq2.results.tsv", sep="\t", header = TRUE)

D4WTvsD0wt<-read.table("nfout/build38/differentialabundance_ALLcomp/tables/differential/D4WT_VS_D0WT.deseq2.results.tsv", sep="\t", header=TRUE)

rawcounts<-read.table("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4/in/rawcounts.tsv",
                      sep="\t",
                      header=TRUE
)

normcounts<-read.table("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4/nfout/build38/differentialabundance_ALLcomp/tables/processed_abundance/all.normalised_counts.tsv", sep="\t", header=TRUE)



#start analysis

txdb <- makeTxDbFromGFF(annotations)

tx_lengths <- transcriptLengths(txdb, with.cds_len = TRUE)
longest_tx_per_gene <- aggregate(tx_lengths$tx_len, by=list(tx_lengths$gene_id), max)

colnames(longest_tx_per_gene) <- c("GeneID", "Length")

write.csv(longest_tx_per_gene, "/mnt/datawk1/references/annotations/ncbi/mm10_Build38/longest_transcript_lengths.csv", row.names = FALSE)

sign_dodko<-D0.DKOvsWT[which(D0.DKOvsWT$padj<=0.05),]
undiff<-D0.DKOvsWT[which(D0.DKOvsWT$padj>0.05 | is.na(D0.DKOvsWT$padj) ),]

sign_dodko<-sign_dodko[which(sign_dodko$log2FoldChange<=-1 ),]

rawcounts<-read.table("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4/in/rawcounts.tsv",
                      sep="\t",
                      header=TRUE
)

raw.glength<-merge(longest_tx_per_gene,rawcounts, by.x="GeneID", by.y="gene_id")


samples<-colnames(raw.glength)[3:length(colnames(raw.glength))]

for (i in 1:length(samples)){
  print(i)
  sam<-samples[i]
  counts <- raw.glength[,sam]
  genelengths<-raw.glength[,"Length"]
  dge <- DGEList(counts = counts)
  fpkm <- rpkm(dge, gene.length = genelengths)

  if(i==1){
    fpkm.tab<-data.frame(raw.glength[,1],fpkm)
    colnames(fpkm.tab)<-c("geneID",sam)
  }else{
    fpkm.tab<-cbind(fpkm.tab,fpkm)
    colnames(fpkm.tab)[length(colnames(fpkm.tab))]<-sam
  }

}

sign_dodko

#which(sign_dodko$gene_id %in% fpkm.tab$geneID)
down.genes<-fpkm.tab[which( fpkm.tab$geneID %in% sign_dodko$gene_id),]
undiff.fpkm<-fpkm.tab[which(fpkm.tab$geneID %in% undiff$gene_id),]
all.fpkm<-fpkm.tab

#groups<-c("D0DoubleKO", "D0Mll1KO","D0Mll2KO", "D0WTA")
groups<-c("D0DoubleKO", "D0WTA")

for (i in 1:length(groups)){

  tarval<-data.frame(down.genes[,grep(paste0(groups[i],"_"),colnames(down.genes))])
  tarmean<-data.frame(down.genes[,1],rowMeans(tarval))

  undiff.val<-data.frame(undiff.fpkm[,grep(paste0(groups[i],"_"),colnames(undiff.fpkm))])
  undiff.means<-data.frame(undiff.fpkm[,1],rowMeans(undiff.val))

  all.val <- data.frame(all.fpkm[,grep(paste0(groups[i],"_"),colnames(all.fpkm))])
  all.mean <- data.frame(all.fpkm[,1],rowMeans(all.val))

  out<-rbind(
    data.frame(gene=tarmean[,1],
               fpkm.mean=tarmean[,2],
               type=rep("down",nrow(tarmean))
    ),
    data.frame(gene=all.mean[,1],
               fpkm.mean=all.mean[,2],
               type=rep("all",nrow(all.mean))
    )
  )

  assign(groups[i],out)

}

# ggplot(D0DoubleKO, aes(x=type, y=log(fpkm.mean))) +
#   geom_boxplot(outlier.colour="red", outlier.shape=8,
#                outlier.size=4,
#                )

# ggplot(D0WTA, aes(x = type, y = log2(fpkm.mean), fill=type)) +
#    geom_boxplot(outlier.shape = NA) +  # Nasconde gli outlier
#    ggtitle("D0WTA gene fpkm")


head(D0WTA)

all<-rbind(
  cbind(D0WTA,line=rep("D0WTA",nrow(D0WTA))),
  cbind(D0DoubleKO,line=rep("D0DoubleKO",nrow(D0DoubleKO)))
  )

all$line <- factor(all$line, levels = c("D0WTA", "D0DoubleKO"))

#c064ad
#7879b9

ggplot(all, aes(x = type, y = fpkm.mean, colour=line)) +
  geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +  # Nasconde gli outlier
  ggtitle("fpkm boxplot") +
  coord_cartesian(ylim = c(0, 25)) +
  theme_minimal() +
  scale_colour_manual(values = c("D0WTA" = "#c064ad", "D0DoubleKO" = "#7879b9"))+
  theme(legend.title = element_blank())


# ggplot(df, aes(x = variabile_x, y = variabile_y)) +
#   geom_boxplot() +
#   scale_x_continuous(limits = c(min, max))


#
# ggplot(D0DoubleKO, aes(x = type, y = log2(fpkm.mean+1), fill=type)) +
#   geom_boxplot(outlier.shape = NA) +  # Nasconde gli outlier
#   ggtitle("D0DoubleKO")
#
# ggplot(D0DoubleKO, aes(x = type, y = log2(fpkm.mean), fill=type)) +
#   geom_boxplot(outlier.shape = NA) +  # Nasconde gli outlier
#   ggtitle("D0DoubleKO")
#

ggsave("othouts/fpkm_distribution/plot_fpkm.png", dpi = 300, width = 8, height = 6)

head(D0DoubleKO)

#format and export outs

colnames(rawcounts)
#  [1] "gene_id"         "D0DoubleKO_REP1" "D0DoubleKO_REP2" "D0DoubleKO_REP3"
#  [5] "D0Mll1KO_REP1"   "D0Mll1KO_REP2"   "D0Mll1KO_REP3"   "D0Mll2KO_REP1"  
#  [9] "D0Mll2KO_REP2"   "D0Mll2KO_REP3"   "D0WTA_REP1"      "D0WTA_REP2"     
# [13] "D0WTA_REP3"      "D4DoubleKO_REP1" "D4DoubleKO_REP2" "D4DoubleKO_REP3"
# [17] "D4MLL1KO_REP1"   "D4MLL1KO_REP2"   "D4MLL1KO_REP3"   "D4MLL2KO_REP1"  
# [21] "D4MLL2KO_REP2"   "D4MLL2KO_REP3"   "D4WT_REP1"       "D4WT_REP2"      
# [25] "D4WT_REP3"


fulltab.D0.DKOvsWT<-merge(normcounts[,c("gene_id", "D0DoubleKO_REP1", 
                                                  "D0DoubleKO_REP2", 
                                                  "D0DoubleKO_REP3",
                                                  "D0WTA_REP1",
                                                  "D0WTA_REP2",
                                                  "D0WTA_REP3")], D0.DKOvsWT, by=1)

fulltab.D0.MLL1KOvsWT<-merge(normcounts[,c("gene_id", "D0Mll1KO_REP1", 
                                                  "D0Mll1KO_REP2", 
                                                  "D0Mll1KO_REP3",
                                                  "D0WTA_REP1",
                                                  "D0WTA_REP2",
                                                  "D0WTA_REP3")], D0.MLL1KOvsWT, by=1)

fulltab.D0.MLL2KOvsWT<-merge(normcounts[,c("gene_id", "D0Mll2KO_REP1", 
                                                  "D0Mll2KO_REP2", 
                                                  "D0Mll2KO_REP3",
                                                  "D0WTA_REP1",
                                                  "D0WTA_REP2",
                                                  "D0WTA_REP3")], D0.MLL2KOvsWT, by=1)


# fulltab.D0.vsWT<-merge(normcounts[,c("gene_id", "D0DoubleKO_REP1", 
#                                                   "D0DoubleKO_REP2", 
#                                                   "D0DoubleKO_REP3",
#                                                   "D0WTA_REP1",
#                                                   "D0WTA_REP2",
#                                                   "D0WTA_REP3")], D0.DKOvsWT, by=1)


fulltab.D4.DKOvsWT<-merge(normcounts[,c("gene_id", "D4DoubleKO_REP1", 
                                                  "D4DoubleKO_REP2", 
                                                  "D4DoubleKO_REP3",
                                                  "D4WT_REP1",
                                                  "D4WT_REP2",
                                                  "D4WT_REP3")], D4.DKOvsWT, by=1)

fulltab.D4.MLL1KOvsWT<-merge(normcounts[,c("gene_id", "D4MLL1KO_REP1", 
                                                  "D4MLL1KO_REP2", 
                                                  "D4MLL1KO_REP3",
                                                  "D4WT_REP1",
                                                  "D4WT_REP2",
                                                  "D4WT_REP3")], D4.MLL1KOvsWT, by=1)

fulltab.D4.MLL2KOvsWT<-merge(normcounts[,c("gene_id", "D4MLL2KO_REP1", 
                                                  "D4MLL2KO_REP2", 
                                                  "D4MLL2KO_REP3",
                                                  "D4WT_REP1",
                                                  "D4WT_REP2",
                                                  "D4WT_REP3")], D4.MLL2KOvsWT, by=1)

fulltab.D4WTvsD0WT<-merge(normcounts[,c("gene_id",   "D4WT_REP1",
                                            "D4WT_REP2",
                                            "D4WT_REP3",
                                            "D0WTA_REP1",
                                            "D0WTA_REP2",
                                            "D0WTA_REP3"
                                            )], D4WTvsD0wt, by=1)

write.table(fulltab.D0.DKOvsWT, file="othouts/expression_bar_plots/D0_DKOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)

write.table(fulltab.D0.MLL1KOvsWT, file="othouts/expression_bar_plots/D0_MLL1KOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)

write.table(fulltab.D0.MLL2KOvsWT, file="othouts/expression_bar_plots/D0_MLL2KOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)


write.table(fulltab.D4.DKOvsWT, file="othouts/expression_bar_plots/D4_DKOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)

write.table(fulltab.D4.MLL1KOvsWT, file="othouts/expression_bar_plots/D4_MLL1KOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)

write.table(fulltab.D4.MLL2KOvsWT, file="othouts/expression_bar_plots/D4_MLL2KOvsWT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)


write.table(fulltab.D4WTvsD0WT, file="othouts/expression_bar_plots/D4WTvsD0WT.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
)



quit()