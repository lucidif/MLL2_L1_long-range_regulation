setwd("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4")

library("ggplot2")

infolder="./nfout/build38/differentialabundance_ALLcomp/tables/differential/"

target_comparisons=c("D4_MLL1KO_VS_WT","D4MLL2KO_vs_WT","D4_DKO_VS_WT", "D4WT_VS_D0WT")
target_files=c("D4MLL1KO_VS_WT.deseq2.results.tsv", "D4MLL2KO_vs_WT.deseq2.results.tsv","D4_DKO_VS_WT.deseq2.results.tsv","D4WT_VS_D0WT.deseq2.results.tsv")

for (i in 1:length(target_comparisons)){
  
  tmpvar<-read.table(file = paste0(infolder,target_files[i]),
                     sep="\t",
                     header=TRUE
                     )
  assign(target_comparisons[i],tmpvar)
  
  
}

D4MLL2KO_vs_WT_degs <- as.data.frame(D4MLL2KO_vs_WT[ , c("gene_id", "log2FoldChange", "padj") ] )

length(which(D4MLL2KO_vs_WT_degs$padj <= 0.05 & D4MLL2KO_vs_WT_degs$log2FoldChange <= -1 ))
length(which(D4MLL2KO_vs_WT_degs$padj <= 0.05 & D4MLL2KO_vs_WT_degs$log2FoldChange >= 1 ))
nrow(D4MLL2KO_vs_WT_degs[which(is.na(D4MLL2KO_vs_WT_degs$padj)),])
length(which(D4MLL2KO_vs_WT_degs$padj > 0.05 )
#degs_dko_wt_toplot<-degs_dko_wt_toplot[which(!is.na(degs_dko_wt_toplot$padj)),]

D4MLL2KO_vs_WT_gn_down<-D4MLL2KO_vs_WT_degs[which(D4MLL2KO_vs_WT_degs$padj <= 0.05 & D4MLL2KO_vs_WT_degs$log2FoldChange <= -1 ),"gene_id"]
D4MLL2KO_vs_WT_gn_up<-D4MLL2KO_vs_WT_degs[which(D4MLL2KO_vs_WT_degs$padj <= 0.05 & D4MLL2KO_vs_WT_degs$log2FoldChange >= 1 ),"gene_id"]

D4MLL2KO_vs_WT_gdown <-paste(D4MLL2KO_vs_WT_gn_down, collapse =" ")
D4MLL2KO_vs_WT_gup<-paste(D4MLL2KO_vs_WT_gn_up, collapse =" ")

D4_MLL1KO_VS_WT_degs <- as.data.frame(D4_MLL1KO_VS_WT[ , c("gene_id", "log2FoldChange", "padj") ] )

length(which(D4_MLL1KO_VS_WT_degs$padj <= 0.05 & D4_MLL1KO_VS_WT_degs$log2FoldChange <= -1 ))
length(which(D4_MLL1KO_VS_WT_degs$padj <= 0.05 & D4_MLL1KO_VS_WT_degs$log2FoldChange >= 1 ))
nrow(D4_MLL1KO_VS_WT_degs[which(is.na(D4_MLL1KO_VS_WT_degs$padj)),])
length(which(D4_MLL1KO_VS_WT_degs$padj > 0.05 ))

D4_MLL1KO_VS_WT_gn_down<-D4_MLL1KO_VS_WT_degs[which(D4_MLL1KO_VS_WT_degs$padj <= 0.05 & D4_MLL1KO_VS_WT$log2FoldChange <= -1 ),"gene_id"]
D4_MLL1KO_VS_WT_gn_up<-D4_MLL1KO_VS_WT_degs[which(D4_MLL1KO_VS_WT$padj <= 0.05 & D4_MLL1KO_VS_WT_degs$log2FoldChange >= 1 ),"gene_id"]

D4_MLL1KO_VS_WT_gdown <-paste(D4_MLL1KO_VS_WT_gn_down, collapse =" ")
D4_MLL1KO_VS_WT_gup<-paste(D4_MLL1KO_VS_WT_gn_up, collapse =" ")

D4_DKO_VS_WT_degs <- as.data.frame(D4_DKO_VS_WT[ , c("gene_id", "log2FoldChange", "padj") ] )

length(which(D4_DKO_VS_WT_degs$padj <= 0.05 & D4_DKO_VS_WT_degs$log2FoldChange <= -1 ))
length(which(D4_DKO_VS_WT_degs$padj <= 0.05 & D4_DKO_VS_WT_degs$log2FoldChange >= 1 ))
nrow(D4_DKO_VS_WT_degs[which(is.na(D4_DKO_VS_WT_degs$padj)),])
length(which(D4_DKO_VS_WT_degs$padj > 0.05 ))

D4_DKO_VS_WT_gn_down<-D4_DKO_VS_WT_degs[which(D4_DKO_VS_WT_degs$padj <= 0.05 & D4_DKO_VS_WT_degs$log2FoldChange <= -1 ),"gene_id"]
D4_DKO_VS_WT_gn_up<-D4_DKO_VS_WT_degs[which(D4_DKO_VS_WT_degs$padj <= 0.05 & D4_DKO_VS_WT_degs$log2FoldChange >= 1 ),"gene_id"]

D4_DKO_VS_WT_gdown <-paste(D4_DKO_VS_WT_gn_down, collapse =" ")
D4_DKO_VS_WT_gup<-paste(D4_DKO_VS_WT_gn_up, collapse =" ")

D4WT_VS_D0WT_degs <- as.data.frame(D4WT_VS_D0WT[ , c("gene_id", "log2FoldChange", "padj") ] )

length(which(D4WT_VS_D0WT_degs$padj <= 0.05 & D4WT_VS_D0WT_degs$log2FoldChange <= -1 ))
length(which(D4WT_VS_D0WT_degs$padj <= 0.05 & D4WT_VS_D0WT_degs$log2FoldChange >= 1 ))
nrow(D4WT_VS_D0WT_degs[which(is.na(D4WT_VS_D0WT_degs$padj)),])
length(which(D4WT_VS_D0WT_degs$padj > 0.05 ))

D4WT_VS_D0WT_gn_down<-D4WT_VS_D0WT_degs[which(D4WT_VS_D0WT_degs$padj <= 0.05 & D4WT_VS_D0WT_degs$log2FoldChange <= -1 ),"gene_id"]
D4WT_VS_D0WT_gn_up<-D4WT_VS_D0WT_degs[which(D4WT_VS_D0WT_degs$padj <= 0.05 & D4WT_VS_D0WT_degs$log2FoldChange >= 1 ),"gene_id"]
D4WT_VS_D0WT_gdown <-paste(D4WT_VS_D0WT_gn_down, collapse =" ")
D4WT_VS_D0WT_gup<-paste(D4WT_VS_D0WT_gn_up, collapse =" ")

degs_genes<-rbind(D4MLL2KO_vs_WT_down=D4MLL2KO_vs_WT_gdown, 
D4MLL2KO_vs_WT_up=D4MLL2KO_vs_WT_gup, 
D4_MLL1KO_VS_WT_down=D4_MLL1KO_VS_WT_gdown, 
D4_MLL1KO_VS_WT_up=D4_MLL1KO_VS_WT_gup,
D4_DKO_VS_WT_down=D4_DKO_VS_WT_gdown,
D4_DKO_VS_WT_up=D4_DKO_VS_WT_gup,
D4WT_VS_D0WT_down=D4WT_VS_D0WT_gdown,
D4WT_VS_D0WT_up=D4WT_VS_D0WT_gup
)

write.table(degs_genes, file="D4_diff_genes_list.tsv", row.names = TRUE, col.names = FALSE, sep="\t")


induced_d4_genes<-D4WT_VS_D0WT_degs[which( D4WT_VS_D0WT_degs$padj <= 0.05 & D4WT_VS_D0WT_degs$log2FoldChange >= 1 ),"gene_id"]

dowregulated_dko<-D4_DKO_VS_WT_degs[which(D4_DKO_VS_WT_degs$padj <= 0.05 & D4_DKO_VS_WT_degs$log2FoldChange <= -1 ),"gene_id"]

overlap_downDKO_inducedWT<-intersect(dowregulated_dko, induced_d4_genes)
length(overlap_downDKO_inducedWT)

q<-433
m<-1075
k<-3795
tot<-20981

phyper(q - 1, m, tot - m , k, lower.tail = FALSE)






