
tegtf<-read.table("/mnt/datawk1/references/annotations/TEtranscripts/mm10_rmsk_TE.gtf", 
sep="\t",
header=FALSE)

#d0_sigdiff_tb<-read.table("othouts/TEtranscripts/TEexpr_WTAd0_VS_DKOd0_sigdiff_gene_TE.txt")
d0_te_tb<-read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd0_VS_DKOd0_gene_TE_analysis.txt")
d0_cnt<-read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd0_VS_DKOd0.cntTable", header=TRUE)

d4_te_tb<-read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd4_VS_DKOd4_gene_TE_analysis.txt")
d4_cnt<-read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd4_VS_DKOd4.cntTable", header=TRUE)

# Esegui lo split
split_list <- strsplit(tegtf$V9, ";")

# Prendi solo il primo elemento di ogni lista
gene_id <- sapply(split_list, function(x) x[1])
gene_id <- gsub("gene_id ", "", gene_id)

class_id <- sapply(split_list, function(x) x[4])
class_id <- gsub("class_id ", "", class_id)

family_id <- sapply(split_list, function(x) x[3])
family_id <- gsub("family_id ", "", family_id)

TEdf<-data.frame(gene_id, class_id, family_id)

unique(TEdf$class_id)

unique(TEdf[grep("SINE",TEdf$class_id),"gene_id"])
#  [1] "B3"       "ID_B1"    "ID4_"     "B4A"      "B2_Mm2"   "B1_Mus1" 
#  [7] "PB1D10"   "PB1D11"   "B2_Mm1t"  "RSINE1"   "B1F"      "B1_Mm"   
# [13] "B1F2"     "B3A"      "B1F1"     "MIRb"     "PB1D7"    "MIR3"    
# [19] "MIR"      "MIRc"     "B1_Mur4"  "B1_Mus2"  "MER131"   "ID4"     
# [25] "PB1"      "B1_Mur2"  "B1_Mur1"  "B4"       "B1_Mur3"  "PB1D9"   
# [31] "ID2"      "B2_Mm1a"  "ID"       "AmnSINE2" "AmnSINE1" "FAM"     

unique(TEdf[grep("LTR",TEdf$class_id),"gene_id"])


allTEfam<-unique(gene_id)

allTEfam[grep("L1Md",allTEfam)]
allTEfam[grep("L1_Mur3",allTEfam)]
allTEfam[grep("L1_Mur2",allTEfam)]
allTEfam[grep("L1_Mur1",allTEfam)]
allTEfam[grep("L2a",allTEfam)]
allTEfam[grep("L2b",allTEfam)]

intfam<-c("L1Md_T", 
"L1Md_F2", 
"L1Md_A", 
"L1Md_Gf", 
"L1_Mur3",
"L1_Mur2",
"L1_Mur1",
"L2a",
"L2b",
"Lx7",
"Lx8",
"Lx9"
)

intfam <- c(
  # L1
  "L1Md_T", 
  "L1Md_F2", 
  "L1Md_A", 
  "L1Md_Gf", 
  "L1_Mur3",
  "L1_Mur2",
  "L1_Mur1",
  "L2a",
  "L2b",
  "Lx7",
  "Lx8",  # (controlla se il ":" è voluto)
  "Lx9",
  
  # SINE (ID/B1/B2/B4/MIR)
  "ID",
  "ID2",
  "ID4",
  "ID_B1",
  "B1_Mm",
  "B1_Mus1",
  "B1_Mus2",
  "B1_Mur1", 
  "B1_Mur2", 
  "B1_Mur3", 
  "B1_Mur4",
  "B1F",
  "B1F1",
  "B1F2",
  "B2_Mm1a",
  "B2_Mm1a",
  "B2_Mm2",
  "B4",
  "B4A",
  "MIR",
  "MIRb",
  "MIR3",
  "MIRc",
  
  ## LTR – ERV vari
  "ERVB7_1-LTR_MM",
  "ERVB7_3-LTR_MM",
  "ERVB7_4-LTR_MM",
  "ERVB7_2-LTR_MM",
  "ERVB4_1B-LTR_MM",
  "ERVB4_2-LTR_MM",
  "ERVB3_1-LTR_MM",
  "ERVB5_1-LTR_MM",
  "ERVB5_2-LTR_MM",
  "ERVL-int",
  "ERVL-E-int",
  "ERVL-B4-int",
  "MERVL-int",
  "MERVL_2A-int",
  "MERV1_I-int",
  "MERV1_LTR",
  "RodERV21-int",
  "MurERV4-int",
  "MurERV4_19-int",
  "MURVY-int",
  "MURVY-LTR",
  "MYSERV-int",
  "MYSERV6-int",
  "MYSERV16_I-int",
  "MMERVK9E_I-int",
  "MMERVK9C_I-int",
  "MMERVK10C-int",
  "MMERVK10D3_I-int",
  "MMERVK10D3_LTR",
  "RNERVK23-int",
  "MMVL30-int",
  "HERV16-int",
  "HERVL40-int",
  "HERVL74-int",
  "MuLV-int",
  "MMTV-int",
  "SRV_MM-int",
  "MuRRS-int",
  "MuRRS4-int",
  "MamGyp-int",
  "MamGypLTR1a",
  "MamGypLTR1b",
  "MamGypLTR1c",
  "MamGypLTR1d",
  "MamGypLTR2b",
  "MamGypLTR2c",
  "MamGypLTR3",
  "MamGypLTR3a",
  
  ## LTR – IAP
  "IAPLTR1a_Mm",
  "IAPLTR1_Mm",
  "IAPLTR2a",
  "IAPLTR2a2_Mm",
  "IAPLTR2b",
  "IAPLTR2_Mm",
  "IAPLTR3",
  "IAPLTR3-int",
  "IAPLTR4",
  "IAPLTR4_I",
  "IAP1-MM_LTR",
  "IAP1-MM_I-int",
  "IAP-d-int",
  "IAPEy-int",
  "IAPEz-int",
  "IAPEY2_LTR",
  "IAPEY3_LTR",
  "IAPEY3-int",
  "IAPEY3C_LTR",
  "IAPEY4_LTR",
  "IAPEY4_I-int",
  "IAPEY5_LTR",
  "IAPEY5_I-int",
  "IAPA_MM-int",
  
  ## LTR – ETn / MMETn / ETnERV
  "MMETn-int",
  "ETnERV-int",
  "ETnERV2-int",
  "ETnERV3-int"
)



#merge(intfam, d0_sigdiff_tb, by.x=1, by.y=1)

for (i in 1:length(intfam)){
    print(d0_cnt[grep(intfam[i], d0_cnt[,1]),])
}


for (i in 1:length(intfam)){
    if(i==1){
        tmptt<-d0_te_tb[grep(intfam[i], rownames(d0_te_tb)),]
    }else{
        adder<-d0_te_tb[grep(intfam[i], rownames(d0_te_tb)),]
        tmptt<-rbind(tmptt, adder)
    }
}

tmptt



for (i in 1:length(intfam)){
    if(i==1){
        tmpttd4<-d4_te_tb[grep(intfam[i], rownames(d4_te_tb)),]
    }else{
        adderd4<-d4_te_tb[grep(intfam[i], rownames(d4_te_tb)),]
        tmpttd4<-rbind(tmpttd4, adderd4)
    }
}

tmpttd4

View(tmpttd4)

for (i in 1:length(intfam)){
    if(i==1){
        cntd4<-d4_cnt[grep(intfam[i], d4_cnt[,1]),]
    }else{
        adderd4<-d4_cnt[grep(intfam[i], d4_cnt[,1]),]
        cntd4<-rbind(cntd4, adderd4)
    }
}


for (i in 1:length(intfam)){
    if(i==1){
        cntd4<-d4_cnt[grep(intfam[i], d4_cnt[,1]),]
    }else{
        adderd4<-d4_cnt[grep(intfam[i], d4_cnt[,1]),]
        cntd4<-rbind(cntd4, adderd4)
    }
}



#================================== d4
#setwd("/mnt/datawk1/analysis/Lara/DE_RNAseq_lara_day0_4/othouts/TEtranscripts")
data <- read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd4_VS_DKOd4.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",3)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)

dds<-estimateSizeFactors(dds)

norm.cnt<-counts(dds, normalized = TRUE)
colnames(norm.cnt) <- sub(".*RNAseq_lara_day0_4_random\\.(.*)\\.markdup\\.sorted\\.bam\\..*", "\\1", colnames(norm.cnt))


dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)
res <- results(dds)
# write.table(res, file="TEexpr_WTAd4_VS_DKOd4_gene_TE_analysis.txt", sep="\t",quote=F)
# resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
# write.table(resSig, file="TEexpr_WTAd4_VS_DKOd4_sigdiff_gene_TE.txt",sep="\t", quote=F)

#================================== d0
data <- read.table("outs/DE_RNAseq_lara_day0_4/othouts/TEtranscripts/TEexpr_WTAd0_VS_DKOd0.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",3)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds<-estimateSizeFactors(dds)

norm.cnt.d0<-counts(dds, normalized = TRUE)
colnames(norm.cnt.d0) <- sub(".*RNAseq_lara_day0_4_random\\.(.*)\\.markdup\\.sorted\\.bam\\..*", "\\1", colnames(norm.cnt.d0))

#dds$groups = relevel(dds$groups,ref="CGroup")
#dds <- DESeq(dds)
#res <- results(dds)
#write.table(res, file="TEexpr_WTAd0_VS_DKOd0_gene_TE_analysis.txt", sep="\t",quote=F)
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
#write.table(resSig, file="TEexpr_WTAd0_VS_DKOd0_sigdiff_gene_TE.txt",sep="\t", quote=F)

#=====================================

#==================================

value<-c()
sample<-c()
family<-c()
time<-c()

for (i in 1:length(intfam)){
    value <- c(value,as.matrix(norm.cnt[grep(intfam[i], rownames(norm.cnt)),]))
    print("==")
    print(i)
    print(length(as.matrix(norm.cnt[grep(intfam[i], rownames(norm.cnt)),])))
    sample <- c(sample,colnames(norm.cnt))#c(sample,as.matrix(norm.cnt[grep(intfam[i], rownames(norm.cnt)),]))
    print("==")
    family<-c(family,rep(intfam[i],ncol(norm.cnt)))
    time<- c(time,rep("d4",ncol(norm.cnt)))

}

df<-data.frame(sample=sample,
           family=family,
           value=value,
           time=time, 
           value=value,
           name=sub("_REP\\d+", "", sample)        
)

value.d0<-c()
sample.d0<-c()
family.d0<-c()
time.d0<-c()
for (i in 1:length(intfam)){ 
    fmtar<-as.matrix(norm.cnt.d0[grep(paste0(intfam[i],":"), rownames(norm.cnt.d0)),]) 
    #print(fmtar)
    length(fmtar)

    if(length(fmtar)==0){
      print(intfam[i])
    }
    
    value.d0 <- c(value.d0,as.matrix(norm.cnt.d0[grep(paste0(intfam[i],":"), rownames(norm.cnt.d0)),]))
    
    #print("==")
    #print(i)
    #print(length(as.matrix(norm.cnt[grep(intfam[i], rownames(norm.cnt.d0)),])))
    sample.d0 <- c(sample.d0,colnames(norm.cnt.d0))#c(sample,as.matrix(norm.cnt[grep(intfam[i], rownames(norm.cnt)),]))
    #print("==")

    if(length(rep(intfam[i],ncol(norm.cnt.d0)))>6){
      print (intfam)
    }

    family.d0<-c(family.d0,rep(intfam[i],ncol(norm.cnt.d0)))



    time.d0<- c(time.d0,rep("d0",ncol(norm.cnt.d0)))

}

length(value.d0)
length(family.d0)


df.d0<-data.frame(sample=sample.d0,
           family=family.d0,
           value=value.d0,
           time=time.d0, 
           value=value.d0,
           name=sub("_REP\\d+", "", sample.d0)        
)

length(sample.d0)
length(family.d0)
length(value.d0)
length(time.d0)

value.d0

df$grouped_family <- paste(df$name, df$family, sep = "_")
df.d0$grouped_family <- paste(df.d0$name, df.d0$family, sep = "_")

df.final<-rbind(df,df.d0)

ordered_families <- c(
  "L1Md_T",
  "L1Md_A",
  "L1Md_Gf", 
  "L1Md_F2", 
  "L1_Mur3",
  "L1_Mur2",
  "L1_Mur1",
  "L2a",
  "L2b",
  "Lx7",
  "Lx8:",
  "Lx9"
)

df$family <- factor(df$family, levels = ordered_families)
df.d0$family <- factor(df.d0$family, levels = ordered_families)

df$name <- factor(df$name, levels = c("D4WT", "D4DoubleKO"))
df.d0$name <- factor(df.d0$name, levels = c("D0WTA", "D0DoubleKO"))

library(ggplot2)

dot.d4<-ggplot(df, aes(x = family, y = log(value), color = name)) +
  geom_point(position = position_dodge(width = 0.7), size = 4, alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(x = "Family",
       y = "Expression",
       color = "Condition",
       title = "D4 WT vs KO") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

dot.d0<-ggplot(df.d0, aes(x = family, y = log(value), color = name)) +
  geom_point(position = position_dodge(width = 0.7), size = 4, alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(x = "Family",
       y = "Expression",
       color = "Condition",
       title = "D0 WT vs KO") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave("./othouts/TEtranscripts/dot.d4.png", plot = dot.d4, width = 10, height = 6, dpi = 300)
ggsave("./othouts/TEtranscripts/dot.d0.png", plot = dot.d0, width = 10, height = 6, dpi = 300)

# ggplot(df, aes(x = grouped_family, y = value, color = name)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
#              size = 2, alpha = 0.7) +
#   theme_minimal(base_size = 12) +
#   labs(x = "Condition + Family",
#        y = "Expression",
#        color = "Condition",
#        title = "Dot plot per family e condizione") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(hjust = 0.5))

df$log2_value <- log2(df$value)
summary_log2 <- df %>%
  mutate(log2_value = log2(value)) %>%
  group_by(family, name) %>%
  summarise(mean = mean(log2_value), sd = sd(log2_value), .groups = "drop")

library(dplyr)
summary_df <- df %>%
  group_by(family, name) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )

sd.dot.plot.d4<-ggplot(df, aes(x = family, y = log2_value, color = name)) +
  # Barre di errore (SD)
  geom_errorbar(data = summary_log2,
                aes(x = family, y = mean, ymin = mean - sd, ymax = mean + sd, color = name),
                position = position_dodge(width = 0.6),
                width = 0.2,
                linewidth = 1.3,
                inherit.aes = FALSE) +

  # Punto medio
  geom_point(data = summary_log2,
             aes(x = family, y = mean, color = name),
             position = position_dodge(width = 0.6),
             shape = 18, size = 3,
             inherit.aes = FALSE) +
  scale_color_manual(values = c("D4WT" = "#715eee", "D4DoubleKO" = "#ffb00e")) +
  theme_minimal(base_size = 12) +
  labs(x = "Family", y = "log2(Expression)", color = "Condition",
       title = "Dot plot log2 con SD - d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(5, 15))
#geom_pointrange(aes(ymin=value-sd, ymax=value+sd))

df.d0$log2_value <- log2(df.d0$value)

summary_log2_d0 <- df.d0 %>%
  group_by(family, name) %>%
  summarise(
    mean = mean(log2_value),
    sd = sd(log2_value),
    .groups = "drop"
  )

sd.dot.plot.d0<-ggplot(df.d0, aes(x = family, y = log2_value, color = name)) +
  # Punti individuali
  #geom_point(position = position_dodge(width = 0.6), size = 6, alpha = 0.7) +

  # Barre di errore (SD)
  geom_errorbar(data = summary_log2_d0,
                aes(x = family, y = mean, ymin = mean - sd, ymax = mean + sd, color = name),
                position = position_dodge(width = 0.6),
                width = 0.2,
                linewidth = 1.3,
                inherit.aes = FALSE) +

  # Punto medio
  geom_point(data = summary_log2_d0,
             aes(x = family, y = mean, color = name),
             position = position_dodge(width = 0.6),
             shape = 18, size = 3,
             inherit.aes = FALSE) +
  scale_color_manual(values = c("D0WTA" = "#715eee", "D0DoubleKO" = "#ffb00e")) +
  theme_minimal(base_size = 12) +
  labs(x = "Family", y = "log2(Expression)", color = "Condition",
       title = "Dot plot log2 con SD - d0") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
         coord_cartesian(ylim = c(5, 15))

ggsave("./othouts/TEtranscripts/sd_dot_plot_d4.png", plot = sd.dot.plot.d4, width = 10, height = 6, dpi = 300)
ggsave("./othouts/TEtranscripts/sd_dot_plot_d0.png", plot = sd.dot.plot.d0, width = 10, height = 6, dpi = 300)

sd.dot.plot.d4.nogrid <- ggplot(df, aes(x = family, y = log2_value, color = name)) +
  geom_errorbar(data = summary_log2,
                aes(x = family, y = mean, ymin = mean - sd, ymax = mean + sd, color = name),
                position = position_dodge(width = 0.6),
                width = 0.2, linewidth = 1.3, inherit.aes = FALSE) +
  geom_point(data = summary_log2,
             aes(x = family, y = mean, color = name),
             position = position_dodge(width = 0.6),
             shape = 18, size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = c("D4WT" = "#715eee", "D4DoubleKO" = "#ffb00e")) +
  theme_minimal(base_size = 12) +
  labs(x = "Family", y = "log2(Expression)", color = "Condition",
       title = "Dot plot log2 con SD - d4") +
  coord_cartesian(ylim = c(5, 15)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

sd.dot.plot.d0.nogrid <- ggplot(df.d0, aes(x = family, y = log2_value, color = name)) +
  geom_errorbar(data = summary_log2_d0,
                aes(x = family, y = mean, ymin = mean - sd, ymax = mean + sd, color = name),
                position = position_dodge(width = 0.6),
                width = 0.2, linewidth = 1.3, inherit.aes = FALSE) +
  geom_point(data = summary_log2_d0,
             aes(x = family, y = mean, color = name),
             position = position_dodge(width = 0.6),
             shape = 18, size = 3, inherit.aes = FALSE) +
  scale_color_manual(values = c("D0WTA" = "#715eee", "D0DoubleKO" = "#ffb00e")) +
  theme_minimal(base_size = 12) +
  labs(x = "Family", y = "log2(Expression)", color = "Condition",
       title = "Dot plot log2 con SD - d0") +
  coord_cartesian(ylim = c(5, 15)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("./othouts/TEtranscripts/nogrid_sd_dot_plot_d4.png", plot = sd.dot.plot.d4.nogrid, width = 10, height = 6, dpi = 300)
ggsave("./othouts/TEtranscripts/nogrid_sd_dot_plot_d0.png", plot = sd.dot.plot.d0.nogrid, width = 10, height = 6, dpi = 300)







df.final<-rbind(df.d0, df)
#df.final$group <- paste(df.final$name, df.final$time, sep = "_")
df.final$group <- df.final$name

#df.final$group <- paste(df.final$name, df.final$time, sep = "_")
# df.final$group <- factor(df.final$group,
#   levels = c("D4WT_d0", "D4DoubleKO_d0", "D4WT_d4", "D4DoubleKO_d4")
# )

df.final$group <- factor(df.final$group,
  levels = c("D0WTA", "D0DoubleKO", "D4WT", "D4DoubleKO")
)


boxplot<-ggplot(df.final, aes(x = family, y = log2(value), fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.9) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.5, alpha = 0.6) +
  theme_minimal(base_size = 13) +
  labs(x = "Family", y = "log2(Expression)", fill = "Condition_Time",
       title = "Expression per family - WT/KO a d0 e d4") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("./othouts/TEtranscripts/boxplot.png", plot = boxplot, width = 10, height = 6, dpi = 300)

