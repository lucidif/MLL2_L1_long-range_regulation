
setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

#=====================================================
#             functions
#=====================================================

source("git/downstream_multiomic/bin/hypergeometric.R")
source("git/downstream_multiomic/bin/distance_functions.R")

#=====================================================
#
#=====================================================

outpath=paste0(getwd(),"/outs/overlap_cutoff500bp/")

#extract degs test

#degs,targetGenes,anno
allanno.genes<-read.table("in/GREATv4/GREATv4.genes.mm10.tsv",sep="\t",header=FALSE)

gtf.file<-read.table("in/build38_DEseq2_RNAseq/mm10.refGene.gtf.gz", sep="\t")

degs.anno<-read.table("in/build38_DEseq2_RNAseq/mm10.anno.tsv", header=TRUE, sep="\t")
degs.anno<-cbind(ensembl=degs.anno$gene_id, geneName=degs.anno$gene_name)
degs.anno<-degs.anno[which(duplicated(degs.anno)==FALSE),]


length(unique(degs.anno[,1]))
length(unique(degs.anno[,2]))


dcm_basal_ggenes<-read.table("outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcm/20241126-public-4.0.4-hpxvmj-mm10-all-gene.txt",
                             sep="\t",header=FALSE
                             )

#====================================
#    DEGs
#====================================


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


#D0 hypers

#pcm.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
#pcm.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
#pcm.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

#pcp.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
#pcp.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=pcp_basal_ggenes, anno=allanno.genes)
#pcp.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=pcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
dcm.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcm.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D0.hyper<-rbind(dcm.D0.down.basal,dcm.D0.up.basal,dcm.D0.undiff.basal)

dcm.D0.hyper <- cbind(row.names(dcm.D0.hyper), dcm.D0.hyper)

write.table(dcm.D0.hyper, 
            file=paste0(outpath,"/dcm.D0.hyper.tsv"),
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE
            )

#dcp.D0.down.basal<-hyper_degs_tgenes(degs=D0.down, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
#dcp.D0.up.basal<-hyper_degs_tgenes(degs=D0.up, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)
#dcp.D0.undiff.basal<-hyper_degs_tgenes(degs=D0.unch, targetGenes=dcp_basal_ggenes, anno=allanno.genes, anno.col = 5)

# D4

dcm.D4.down.basal<-hyper_degs_tgenes(degs=D4.down, targetGenes=dcm_basal_ggenes , anno=allanno.genes, anno.col = 5)
dcm.D4.up.basal<-hyper_degs_tgenes(degs=D4.up, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)
dcm.D4.undiff.basal<-hyper_degs_tgenes(degs=D4.unch, targetGenes=dcm_basal_ggenes, anno=allanno.genes, anno.col = 5)

dcm.D4.hyper<-rbind(dcm.D4.down.basal,dcm.D4.up.basal,dcm.D4.undiff.basal)

dcm.D4.hyper <- cbind(row.names(dcm.D4.hyper), dcm.D4.hyper)

write.table(dcm.D4.hyper, 
            file=paste0(outpath,"/dcm.D4.hyper.tsv"),
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE
)

#===============D4 wtD4D0less down / up

#pcm.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcm_basal_ggenes , anno=allanno.genes)
#pcm.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcm_basal_ggenes, anno=allanno.genes)

#pcp.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=pcp_basal_ggenes, anno=allanno.genes)
#pcp.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=pcp_basal_ggenes, anno=allanno.genes)

dcm.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcm_basal_ggenes , anno=allanno.genes)
dcm.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcm_basal_ggenes, anno=allanno.genes)

#dcp.D4lessD4D0.down.basal<-hyper_degs_tgenes(degs=D4degs_down_lessWT, targetGenes=dcp_basal_ggenes, anno=allanno.genes)
#dcp.D4lessD4D0.up.basal<-hyper_degs_tgenes(degs=D4degs_up_lessWT, targetGenes=dcp_basal_ggenes, anno=allanno.genes)

#=================
#
#=================

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
