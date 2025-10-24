
setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

#=====================================================
#             functions
#=====================================================

source("git/downstream_multiomic/bin/hypergeometric.R")
source("git/downstream_multiomic/bin/distance_functions.R")

#=====================================================
#
#=====================================================

outpath=paste0(getwd(),"/outs/overlap/")

#extract degs test

#degs,targetGenes,anno
allanno.genes<-read.table("in/GREATv4/GREATv4.genes.mm10.tsv",sep="\t",header=FALSE)

gtf.file<-read.table("in/build38_DEseq2_RNAseq/mm10.refGene.gtf.gz", sep="\t")

degs.anno<-read.table("in/build38_DEseq2_RNAseq/mm10.anno.tsv", header=TRUE, sep="\t")
degs.anno<-cbind(ensembl=degs.anno$gene_id, geneName=degs.anno$gene_name)
degs.anno<-degs.anno[which(duplicated(degs.anno)==FALSE),]


length(unique(degs.anno[,1]))
length(unique(degs.anno[,2]))

#int_basal<-read.table("outs/distfiltered_great/basal/20240704-public-4.0.4-HYyfkD-mm10-all-gene.txt", sep="\t",header=FALSE)

#int_nearest <- read.table("outs/distfiltered_great/nearest/20240704-public-4.0.4-LTZU3b-mm10-all-gene.txt", sep="\t",header=FALSE)

dcm_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/dcm/20240719-public-4.0.4-nxvS8V-mm10-all-gene.txt", sep="\t",header=FALSE)

dcm_basal_ggenes<-read.table("outs/great/filtered_by_distance/basal/dcm/20240719-public-4.0.4-nxvS8V-mm10-all-gene.txt",
                             sep="\t",header=FALSE
                             )

dcp_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/dcp/20240719-public-4.0.4-qoHtEr-mm10-all-gene.txt", 
                               sep="\t",
                               header=FALSE
                               )
dcp_basal_ggenes<-read.table("outs/great/filtered_by_distance/basal/dcp/20240719-public-4.0.4-qoHtEr-mm10-all-gene.txt", 
                               sep="\t",
                               header=FALSE
)

pcm_basal_ggenes<-read.table("outs/great/filtered_by_distance/basal/pcm/20240719-public-4.0.4-Zd7s4R-mm10-all-gene.txt",
                             sep="\t",
                             header=FALSE                             
                             )

pcp_basal_ggenes<-read.table("outs/great/filtered_by_distance/basal/pcp/20240719-public-4.0.4-9RkkYc-mm10-all-gene.txt",
                             sep="\t",
                             header=FALSE                             
)

pcp_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/pcp/20240719-public-4.0.4-9RkkYc-mm10-all-gene.txt",
                             sep="\t",
                             header=FALSE                             
)


pcm_nearest_ggenes<-read.table("outs/great/filtered_by_distance/nearest/pcm/20240719-public-4.0.4-Zd7s4R-mm10-all-gene.txt",
                               sep="\t",
                               header=FALSE                             
)



degs_D0<-read.table("in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv", header=TRUE, sep="\t")
degs_D0<-merge(degs.anno,degs_D0,  by=1)

degs_D4<-read.table("in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv", header=TRUE, sep="\t")
degs_D4<-merge(degs.anno,degs_D4,  by=1)

degs_WT<-read.table("in/build38_DEseq2_RNAseq/D4WT_VS_D0WT.deseq2.results.tsv", header=TRUE, sep="\t")
degs_WT<-merge(degs.anno,degs_WT,  by=1)

D0.up<-extract_degs(degs_D0, direction="up" )
D0.down<-extract_degs(degs_D0, direction="down" )
D0.unch<-extract_degs(degs_D0, pval = NULL )

D4.up<-extract_degs(degs_D4, direction="up" )
D4.down<-extract_degs(degs_D4, direction="down" )
D4.unch<-extract_degs(degs_D4, pval = NULL )

WT.sign.down<-extract_degs(degs_WT, direction="down" )
WT.sign.up<-extract_degs(degs_WT, direction="up" )

D4degs_down_lessWT <- setdiff(D4.down$ensembl, WT.sign.down$ensembl)
D4degs_down_lessWT<-D4.down[which(D4degs_down_lessWT %in% D4.down$ensembl),]

D4degs_up_lessWT <- setdiff(D4.up$ensembl, WT.sign.up$ensembl)
D4degs_up_lessWT<-D4.up[which(D4degs_up_lessWT %in% D4.up$ensembl),]


#library(VennDiagram)

# venn.diagram(
#   x = list(WT.sign.down$ensembl, D4.down$ensembl),
#   category.names = c("WTsign" , "D4down" ),
#   filename = './output/degs_WTvsD4down',
#   output=TRUE
# )

#allanno.genes<-allanno.genes[,c("V5","V1","V2","V3","V4")]
#tt.down<-merge(allanno.genes,tt.down, by.x=1,by.y=2)

#D0 hypers

pcm.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
pcm.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
pcm.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

pcp.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcp_basal_ggenes, anno=allanno.genes)
pcp.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
dcm.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcm.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcp.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

#===nearest

pcm.D0.down.nearest<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcm_nearest_ggenes , anno=allanno.genes, anno.col = 5)
pcm.D0.up.nearest<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)
  pcm.D0.undiff.nearest<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)

pcp.D0.down.nearest<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D0.up.nearest<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D0.undiff.nearest<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D0.down.nearest<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcm_nearest_ggenes , anno=allanno.genes, anno.col = 5)
dcm.D0.up.nearest<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)
dcm.D0.undiff.nearest<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)

dcp.D0.down.nearest<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D0.up.nearest<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D0.undiff.nearest<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)


# D4
pcm.D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=pcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
pcm.D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
pcm.D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

pcp.D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=dcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
dcm.D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcm.D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcp.D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcp.D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

#===nearest

pcm.D4.down.nearest<-hyper_degs_tgenes(degs=D4.down, targetGenes=pcm_nearest_ggenes , anno=allanno.genes, anno.col = 5)
pcm.D4.up.nearest<-hyper_degs_tgenes(degs=D4.up, targetGenes=pcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)
pcm.D4.undiff.nearest<-hyper_degs_tgenes(degs=D4.unch, targetGenes=pcm_nearest_ggenes, anno=allanno.genes, anno.col = 5)

pcp.D4.down.nearest<-hyper_degs_tgenes(degs=D4.down, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D4.up.nearest<-hyper_degs_tgenes(degs=D4.up, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)
pcp.D4.undiff.nearest<-hyper_degs_tgenes(degs=D4.unch, targetGenes=pcp_nearest_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D4.down.nearest<-hyper_degs_tgenes(degs=D4.down, targetGenes=dcm_nearest_ggenes , anno=allanno.genes)
dcm.D4.up.nearest<-hyper_degs_tgenes(degs=D4.up, targetGenes=dcm_nearest_ggenes, anno=allanno.genes)
dcm.D4.undiff.nearest<-hyper_degs_tgenes(degs=D4.unch, targetGenes=dcm_nearest_ggenes, anno=allanno.genes)

dcp.D4.down.nearest<-hyper_degs_tgenes(degs=D4.down, targetGenes=dcp_nearest_ggenes, anno=allanno.genes)
dcp.D4.up.nearest<-hyper_degs_tgenes(degs=D4.up, targetGenes=dcp_nearest_ggenes, anno=allanno.genes)
dcp.D4.undiff.nearest<-hyper_degs_tgenes(degs=D4.unch, targetGenes=dcp_nearest_ggenes, anno=allanno.genes)

#===============D4 wtD4D0less down / up

pcm.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcm_basal_ggenes , anno=allanno.genes)
pcm.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcm_basal_ggenes, anno=allanno.genes)

pcp.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcp_basal_ggenes, anno=allanno.genes)
pcp.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcp_basal_ggenes, anno=allanno.genes)

dcm.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcm_basal_ggenes , anno=allanno.genes)
dcm.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcm_basal_ggenes, anno=allanno.genes)

dcp.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcp_basal_ggenes, anno=allanno.genes)
dcp.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcp_basal_ggenes, anno=allanno.genes)

#===nearest

pcm.D4lessD4D0.down.nearest<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcm_nearest_ggenes , anno=allanno.genes)
pcm.D4lessD4D0.up.nearest<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcm_nearest_ggenes, anno=allanno.genes)

pcp.D4lessD4D0.down.nearest<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcp_nearest_ggenes, anno=allanno.genes)
pcp.D4lessD4D0.up.nearest<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcp_nearest_ggenes, anno=allanno.genes)

dcm.D4lessD4D0.down.nearest<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcm_nearest_ggenes , anno=allanno.genes)
dcm.D4lessD4D0.up.nearest<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcm_nearest_ggenes, anno=allanno.genes)

dcp.D4lessD4D0.down.nearest<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcp_nearest_ggenes, anno=allanno.genes)
dcp.D4lessD4D0.up.nearest<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcp_nearest_ggenes, anno=allanno.genes)





#=================
#
#=================




# print(paste0("D0 basal hyper test: down_pval=", D0.down.basal[6],
#              " up_pval=",D0.up.basal[6],
#              " undiff_pval=",D0.undiff.basal[6]
#             )
#       )
# 
# print(paste0("D0 nearest hyper test: down_pval=", D0.down.nearest[6],
#              " up_pval=",D0.up.nearest[6],
#              " undiff_pval=",D0.undiff.nearest[6]
#             )
# )


#D4 hypers

# D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=int_basal, anno=allanno.genes)
# D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=int_basal, anno=allanno.genes)
# D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=int_basal, anno=allanno.genes)
# 
# D4.down.nearest<-hyper_degs_tgenes(degs=D4.down, targetGenes=int_nearest, anno=allanno.genes)
# D4.up.nearest<-hyper_degs_tgenes(degs=D4.up, targetGenes=int_nearest, anno=allanno.genes)
# D4.undiff.nearest<-hyper_degs_tgenes(degs=D4.unch, targetGenes=int_nearest, anno=allanno.genes)
# 
# print(paste0("D4 basal hyper test: down_pval=", D4.down.basal[6],
#              " up_pval=",D4.up.basal[6],
#              " undiff_pval=",D4.undiff.basal[6]
# )
# )
# 
# print(paste0("D4 nearest hyper test: down_pval=", D4.down.nearest[6],
#              " up_pval=",D4.up.nearest[6],
#              " undiff_pval=",D4.undiff.nearest[6]
# )
# )

#====== intersect WT d4/d0 with 

WT.all<-extract_degs(degs_WT, direction="all" )

WT_D4<-intersect(WT.all$ensembl,D4.down$ensembl)
noWT_D4<-setdiff(D4.down$ensembl,WT_D4) 

#extract nearest peaks at this genes (GREAT ?)

noWT_D4.down<-D4.down[which( noWT_D4 %in% D4.down$ensembl ),]
# noWT_D4.down.nearest <- hyper_degs_tgenes(degs=noWT_D4.down, targetGenes=int_nearest, anno=allanno.genes)
# 
# WT_D4.down<-D4.down[which( WT_D4 %in% D4.down$ensembl ),]
# WT_D4.down.nearest <- hyper_degs_tgenes(degs=WT_D4.down, targetGenes=int_nearest, anno=allanno.genes)
# 
# 
#=====================================
#distance by significant genes
#=====================================

dcp.D0.down.nearest.genes<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcp_nearest_ggenes, anno=allanno.genes, extract="overlap")

targen<-intersect(allanno.genes$V5, dcp.D0.down.nearest.genes)
length(dcp.D0.down.nearest.genes)
length(targen)

targen.coordinates<-allanno.genes[which(targen %in% allanno.genes$V5 ),]
targen.coordinates<-data.frame(syn=targen.coordinates$V5,
                               start=targen.coordinates[,"V3"]-1,
                               end=targen.coordinates[,"V3"]+1,
                               strand=targen.coordinates$V4
                               )


#extractdistances()


#====================================
# write outputs
#====================================

#/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/outs/overlap/

#outpath=paste0(getwd(),"/outs/overlap/")

#write.table(file=paste0(outpath,"D0.down.tsv"),sep="\t", col.names=TRUE)

# degs_D0<-read.table("in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_D0<-merge(degs.anno,degs_D0,  by=1)

# degs_D4<-read.table("in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_D4<-merge(degs.anno,degs_D4,  by=1)

# degs_WT<-read.table("in/build38_DEseq2_RNAseq/D4WT_VS_D0WT.deseq2.results.tsv", header=TRUE, sep="\t")
# degs_WT<-merge(degs.anno,degs_WT,  by=1)

#D0.up<-extract_degs(degs_D0, direction="up" )
#D0.down
write.table(D0.down, file=paste0(outpath,"D0.down.tsv"),sep="\t", 
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
           )
#D0.unch<-extract_degs(degs_D0, pval = NULL )

#D4.up<-extract_degs(degs_D4, direction="up" )
#D4.down
write.table(D4.down, file=paste0(outpath,"D4.down.tsv"),sep="\t", 
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
           )
#D4.unch<-extract_degs(degs_D4, pval = NULL )

#WT.sign.down<-extract_degs(degs_WT, direction="down" )
#WT.sign.up<-extract_degs(degs_WT, direction="up" )



# file.copy()
