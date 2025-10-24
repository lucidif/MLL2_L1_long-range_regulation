
#use rocker/tidyverse:4.5.1 to run this script

#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")
#setwd("/media/lucio/external.wk/bioinfo/wkdir/Lara/Lara_multiomic_analysis/")

#setwd("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/")

#================ define functions
replace_outliers <- function(x) {
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  x[(x < Q1 - 1.5 * IQR) | (x > Q3 + 1.5 * IQR)] <- NA
  x
}

category_rapresentation<-function(intg_preaktss6,cat.x,cat.y) {

  intg_preaktss_all<-intg_preaktss6[which(intg_preaktss6$V2==cat.x),]

  all_broad_gn<-intersect(intg_preaktss6[which(intg_preaktss6$V2==cat.x),"intg"]
                          ,intg_preaktss6[which(intg_preaktss6$V2==cat.y),"intg"])

  all_broad_gn_peaks<-intg_preaktss_all[which(intg_preaktss_all$intg %in% all_broad_gn ),-2]
  all_broad_gn_peaks<-all_broad_gn_peaks[which(duplicated(all_broad_gn_peaks)==FALSE),]

  return(all_broad_gn_peaks)

}

#gtf<-

#================= define main params

broad_thresholds<-6000


#================  load files

load("assets/mESC_H3K27me3_levels_linked_toEachTSS.RData")
d0dkovswt<-read.table("in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv",
                      sep="\t",
                      header=TRUE
)

proximal_dko_vs_wt<-read.table("/mnt/datawk1/analysis/Lara/Lara_multiomic_analysis/outs/great/Double_KO_vs_F_F/nearest/K4me3_proximal_Double_KO_vs_FF/20250123-public-4.0.4-dq6pJp-mm10-all-gene.txt",
                               sep="\t",
                               header=FALSE
)

#=====================

#source("git/downstream_multiomic/bin/distance_functions.R")
source("git/downstream_multiomic/bin/hypergeometric.R")

library("dplyr")
library("ggplot2")

#=====================


D0.down<-extract_degs(d0dkovswt, direction="down" )

tss_h3k27me3Level<-tss_h3k27me3Level[order(tss_h3k27me3Level$longitudPicos, decreasing = TRUE),]

tss_h3k27me3Level<-data.frame(id=c(1:nrow(tss_h3k27me3Level)),tss_h3k27me3Level)

tss_h3k27me3Level <- tss_h3k27me3Level[!duplicated(tss_h3k27me3Level[, -c(1,5)]),]

which(duplicated(tss_h3k27me3Level$gene))
p<-380
tss_h3k27me3Level[which(tss_h3k27me3Level$gene==tss_h3k27me3Level$gene[p]),]

tss_h3k27me3Level <- tss_h3k27me3Level[which(!duplicated(tss_h3k27me3Level$gene)),]

#=====merge proximal peaks great and tss_h3k27me3 anno

peakFind=rep(NA,nrow(tss_h3k27me3Level))
proxiNoGene<-c()

near_proxi_genes<-unique(proximal_dko_vs_wt$V1)
near_proxi_genes<-cbind(near_proxi_genes, peakFind="yes")

tss_inside_proximal<-merge(tss_h3k27me3Level,near_proxi_genes,by.x=2,by.y=1, all.x=TRUE)

tss_inside_proximal<-tss_inside_proximal[which(tss_inside_proximal$peakFind=="yes"),]


intg<-D0.down$gene_id[D0.down$gene_id %in% tss_h3k27me3Level$gene]

intg<-cbind(intg, rep("down",length(intg)))

intg_preaktss<-merge(intg,tss_h3k27me3Level, by.x=1, by.y=2, all.y=TRUE, sort=TRUE)

intg_preaktss[which(is.na(intg_preaktss$V2)==TRUE),"V2"]<-"undiff"

length(which(intg_preaktss$longitudPicos==0))

length(which(intg_preaktss$longitudPicos<=broad_thresholds))


length(which(intg_preaktss$longitudPicos>broad_thresholds))



intg_preaktss <- intg_preaktss %>%
  arrange(ifelse(V2 == "down", 1, 0))



ggplot(tss_h3k27me3Level, aes(x=id, y=longitudPicos)) +
  geom_point() +
  geom_hline(yintercept = broad_thresholds, linetype = "dashed", color = "red") +
  geom_vline(xintercept = length(which(intg_preaktss$longitudPicos>=broad_thresholds)), linetype = "dashed", color = "red")+
  geom_vline(xintercept = length(intg_preaktss$longitudPicos)-length(which(intg_preaktss$longitudPicos<=0)), linetype = "dashed", color = "red") +
  annotate("text", x = length(which(intg_preaktss$longitudPicos>=broad_thresholds))/2 , y = 50000, label = "broad", color = "red", size = 5, angle = 90)+
  annotate("text", x = length(which(intg_preaktss$longitudPicos>=broad_thresholds))*2 , y = 50000, label = "narrow", color = "red", size = 5, angle = 90)+
  annotate("text", x = (length(intg_preaktss$longitudPicos)-length(which(intg_preaktss$longitudPicos<=0)))*1.5 , y = 50000, label = "negative", color = "red", size = 5, angle = 90)



summary(intg_preaktss[which(intg_preaktss$V2=="down"),"longitudPicos"])

intg_preaktss.down<-intg_preaktss[which(intg_preaktss$V2=="down"),"longitudPicos"]

length(which(intg_preaktss.down==0))

length(which(intg_preaktss.down<=broad_thresholds))


length(which(intg_preaktss.down>broad_thresholds))

intg_preaktss_all<-intg_preaktss

intg_preaktss_all$V2<-"all"

intg_preaktss2<-rbind(intg_preaktss,
                      intg_preaktss_all)

undiff.peaks <- intg_preaktss2[which(intg_preaktss2$V2=="undiff"),"longitudPicos"]
down.peaks <- intg_preaktss2[which(intg_preaktss2$V2=="down"),"longitudPicos"]
all.peaks <- intg_preaktss2[which(intg_preaktss2$V2=="all"),"longitudPicos"]

summary(undiff.peaks)
summary(down.peaks)
summary(all.peaks)

plot(density(undiff.peaks))
plot(density(down.peaks))
plot(density(all.peaks))

intg_preaktss2<-intg_preaktss2[which(intg_preaktss2$V2!="undiff"),]

intg_preaktss_broad<-intg_preaktss2[which(intg_preaktss2$V2=="all" & intg_preaktss2$longitudPicos >= broad_thresholds),]

intg_preaktss_broad$V2<-"broad"

intg_preaktss_narrow<-intg_preaktss2[which(intg_preaktss2$V2=="all" & intg_preaktss2$longitudPicos < broad_thresholds & intg_preaktss2$longitudPicos > 0),]

intg_preaktss_narrow$V2<-"narrow"

intg_preaktss_negative<-intg_preaktss2[which(intg_preaktss2$V2=="all" & intg_preaktss2$longitudPicos == 0),]
intg_preaktss_negative$V2<-"negative"


intg_preaktss3<-rbind(intg_preaktss2,intg_preaktss_broad)
intg_preaktss5<-rbind(intg_preaktss2,
                      intg_preaktss_broad,
                      intg_preaktss_narrow,
                      intg_preaktss_negative
)

colnames(tss_inside_proximal)
colnames(intg_preaktss5)

tss_inside_proximal<-tss_inside_proximal[,c("gene","peakFind","id","peakStart","peakEnd","TSS","chr","longitudPicos")]

colnames(tss_inside_proximal)<-colnames(intg_preaktss5)

tss_inside_proximal$V2<-"proximal_dko_wt_loose"

intg_preaktss6<-rbind(intg_preaktss5,
                      tss_inside_proximal
)


outlayer_rem_intg_preaktss3<-intg_preaktss3
outlayer_rem_intg_preaktss3$longitudPicos<-replace_outliers(intg_preaktss3$longitudPicos)

plot(density(intg_preaktss3[which(intg_preaktss3$V2=="broad"),"longitudPicos"]))
summary(intg_preaktss3[which(intg_preaktss3$V2=="broad"),"longitudPicos"])

broad.peaks <- intg_preaktss3[which(intg_preaktss3$V2=="broad"),"longitudPicos"]

plot(density(broad.peaks))

ggplot(intg_preaktss3, aes(x = longitudPicos, color = V2)) +
  geom_density(size = 1) +
  labs(title = "H3K27me3 peaks TSS dist",
       x = "peak length",
       y = "Density",
       color = "") +
  theme_minimal()

#bin barplot

intg_peaktss4bis<-intg_preaktss6[which(intg_preaktss6$V2=="all" |
                                         intg_preaktss6$V2=="broad" |
                                         intg_preaktss6$V2=="down" |
                                         intg_preaktss6$V2=="proximal_dko_wt_loose"
),]
intg_peaktss_tar<-intg_peaktss4bis

intg_peaktss_tar$longitudPicos

bin_name=c("bin1","bin2","bin3", "bin4", "bin5", "bin6")
binDiv_start<-c(NA,1,3000,broad_thresholds,10000,20000)
binDiv_end<-c(0,3000,broad_thresholds,10000,20000,NA)

assigned_bin<-c()

for (i in 1:length(bin_name)){


  if(!is.na(binDiv_start[i])){
    inbin_start <- which(intg_peaktss_tar$longitudPicos > binDiv_start[i])
  }

  if(!is.na(binDiv_end[i])){
    inbin_end <- which(intg_peaktss_tar$longitudPicos <= binDiv_end[i])
  }

  if(exists("inbin_end") & exists("inbin_start")){
    inside_bin<-intersect(inbin_end, inbin_start)
    rm(inbin_start)
    rm(inbin_end)
  }

  if(!exists("inbin_end") & exists("inbin_start")){
    inside_bin<-inbin_start
    rm(inbin_start)
  }

  if(exists("inbin_end") & !exists("inbin_start")){
    inside_bin<-inbin_end
    rm(inbin_end)
  }

  assigned_bin[inside_bin]<-bin_name[i]



}


histo_preaktss3<-cbind(intg_peaktss_tar,assigned_bin)
histo_preaktss3<-cbind(intg_peaktss_tar,assigned_bin)

nbin=6

down.counts<-t(rbind(table(histo_preaktss3[which(histo_preaktss3$V2=="down"),"assigned_bin"]),
                     rep("down",nbin)
))

down.counts=cbind(bin=row.names(down.counts),down.counts)

all.counts<-t(rbind(
  table(histo_preaktss3[which(histo_preaktss3$V2=="all"),"assigned_bin"]),
  rep("all",nbin)
))

all.counts=cbind(bin=row.names(all.counts),all.counts)

broad.counts<-t(rbind(
  table(histo_preaktss3[which(histo_preaktss3$V2=="broad"),"assigned_bin"]),
  rep("broad",3)
))

broad.counts=cbind(bin=row.names(broad.counts),broad.counts)

proximal.counts<-t(rbind(
  table(histo_preaktss3[which(histo_preaktss3$V2=="proximal_dko_wt_loose"),"assigned_bin"]),
  rep("proximal_dko_wt_loose",nbin)
))

proximal.counts=cbind(bin=row.names(proximal.counts),proximal.counts)

histo.df<-as.data.frame(rbind(down.counts,all.counts,broad.counts, proximal.counts))
colnames(histo.df)<-c("bin","counts","group")

histo.df <- data.frame(
  bin=as.factor(histo.df$bin)
  ,counts=as.numeric(histo.df$counts),
  group=as.character(histo.df$group)
)

perc<-c()

all.total<-sum(histo.df[which(histo.df$group=="all"),"counts"])
down.total<-sum(histo.df[which(histo.df$group=="down"),"counts"])
broad.total<-sum(histo.df[which(histo.df$group=="broad"),"counts"])
proximal.total<-sum(histo.df[which(histo.df$group=="proximal_dko_wt_loose"),"counts"])

perc[which(histo.df$group=="all")]<-histo.df[which(histo.df$group=="all"),"counts"]/all.total
perc[which(histo.df$group=="down")]<-histo.df[which(histo.df$group=="down"),"counts"]/down.total
perc[which(histo.df$group=="broad")]<-histo.df[which(histo.df$group=="broad"),"counts"]/broad.total
perc[which(histo.df$group=="proximal_dko_wt_loose")]<-histo.df[which(histo.df$group=="proximal_dko_wt_loose"),"counts"]/proximal.total

histo.df<-cbind(histo.df,perc=perc)

histo.df$bin <- factor(histo.df$bin, levels = bin_name)

histo.df$bin <- factor(histo.df$bin, levels = c("bin6", "bin5", "bin4", "bin3", "bin2", "bin1"))

ggplot(histo.df, aes(x = bin, y = perc, fill = group)) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  labs(title = "",
       x = "Peaks length (bin)",
       y = "% of Genes within each category") +
  scale_x_discrete(labels = c("bin1" = "0",
                              "bin2" = "1-3000",
                              "bin3" = paste0("3000-",broad_thresholds),
                              "bin4" = paste0(broad_thresholds,"-10000"),
                              "bin5" = "10000-20000",
                              "bin6" = ">20000"
  )) +  
  theme_minimal() +
  coord_flip ()

m<-length(which(intg_preaktss5$V2=="narrow"))+length(which(intg_preaktss5$V2=="broad"))

n<-length(which(intg_preaktss5$V2=="all"))-m



k<-length(which(intg_preaktss5$V2=="down"))

q <- sum(
  intg_preaktss5[which(intg_preaktss5$V2 == "down"), "intg"] %in% unique(
    c(
      intg_preaktss5[which(intg_preaktss5$V2 == "narrow"), "intg"],
      intg_preaktss5[which(intg_preaktss5$V2 == "broad"), "intg"]
    )
  )
)

pval<-phyper(
  q = q - 1,
  m = m,
  n = n,
  k = k,
  lower.tail = FALSE
)


#broad vs negative + narrow

# M (Successi nella popolazione)
m<-length(which(intg_preaktss5$V2=="broad"))

# N (total population). tutti i geni possibili.
n<-length(which(intg_preaktss5$V2=="all"))-m
#is the same that length(which(intg_preaktss5$V2 %in% c("broad", "narrow", "negative")))

# k La dimensione del campione ( DEGs totali).
k<-length(which(intg_preaktss5$V2=="down"))

# q
#q: Il numero di successi osservati nel campione

q <- sum(
  intg_preaktss5[which(intg_preaktss5$V2 == "down"), "intg"] %in% unique(
    intg_preaktss5[which(intg_preaktss5$V2 == "broad"), "intg"]
  )
)

pval<-phyper(
  #q = length(D0.down$gene_id) - 1,
  q = q - 1,
  m = m,
  n = n,
  k = k,
  lower.tail = FALSE
)



#negative vs broad + narrow

# M (Successi nella popolazione)

m<-length(which(intg_preaktss5$V2=="negative"))

# N (total population). tutti i geni possibili.
n<-length(which(intg_preaktss5$V2=="all"))-m
#is the same that length(which(intg_preaktss5$V2 %in% c("broad", "narrow", "negative")))

# k La dimensione del campione ( DEGs totali).
k<-length(which(intg_preaktss5$V2=="down"))

# q
#q: Il numero di successi osservati nel campione

q <- sum(
  intg_preaktss5[which(intg_preaktss5$V2 == "down"), "intg"] %in% unique(
    intg_preaktss5[which(intg_preaktss5$V2 == "negative"), "intg"]
  )
)

pval<-phyper(
  #q = length(D0.down$gene_id) - 1,
  q = q - 1,
  m = m,
  n = n,
  k = k,
  lower.tail = FALSE
)


# narrow vs negative + broad

# M (Successi nella popolazione)

m<-length(which(intg_preaktss5$V2=="narrow"))

# N (total population). tutti i geni possibili.
n<-length(which(intg_preaktss5$V2=="all"))-m
#is the same that length(which(intg_preaktss5$V2 %in% c("broad", "narrow", "negative")))

# k La dimensione del campione ( DEGs totali).
k<-length(which(intg_preaktss5$V2=="down"))

# q
#q: Il numero di successi osservati nel campione

q <- sum(
  intg_preaktss5[which(intg_preaktss5$V2 == "down"), "intg"] %in% unique(
    intg_preaktss5[which(intg_preaktss5$V2 == "narrow"), "intg"]
  )
)

pval<-phyper(
  #q = length(D0.down$gene_id) - 1,
  q = q - 1,
  m = m,
  n = n,
  k = k,
  lower.tail = FALSE
)


dwgn<-intg_peaktss4bis$intg[which(intg_peaktss4bis$V2=="down")]
progn<-intg_peaktss4bis$intg[which(intg_peaktss4bis$V2=="proximal_dko_wt_loose")]

broad_gn_list <- intg_peaktss4bis[which(intg_peaktss4bis$longitudPicos>broad_thresholds & intg_peaktss4bis$V2=="all" ),]

broad_gn_list<-broad_gn_list[,c("chr","TSS","TSS","intg")]
broad_gn_list$TSS<-broad_gn_list$TSS-5000
broad_gn_list$TSS.1<-broad_gn_list$TSS+5000

write.table(broad_gn_list,
            "./outs/broad_all_genes.bed",
            col.names = FALSE,
            row.names = FALSE,
            sep="\t", quote=FALSE)


write.table(intg_peaktss4bis,
            "/home/lucio/Pictures/broad_narrow_distribution.tsv",
            col.names = TRUE,
            row.names = FALSE,
            sep="\t", quote=FALSE)

length(intersect(dwgn,progn))


length(intg_peaktss4bis$intg[which(intg_peaktss4bis$V2=="down")])


length(progn)

# q number of success draws
q<-length(intersect(dwgn,progn))


# M (Successi nella popolazione)
m<-length(intg_peaktss4bis$intg[which(intg_peaktss4bis$V2=="down")])

# N (total population). tutti i geni possibili.
#n the number of black balls in the urn.
n<-length(which(intg_peaktss4bis$V2=="all"))-m

# k numbers of draws
k<-length(progn)


# Calcolo del valore p (probabilità cumulativa ipergeometrica)

p_value <- phyper(q - 1, m, n, k, lower.tail = FALSE)

#fraction_of_catalogated_peaks

#all down proximal_k4me3loss

#from intg_preaktss6
#all intersect narrow and count
#all intersect broad and count
#all intersect negative and count
#all counts sums need to be the number od total all

#down intersect narrow and count
#down intersect broad and count
#down intersect negative and count
#sums of all counts must be the total numebr of down

#same with proximal
#make bar plot in percentage

#fraction_of_catalogated_peaks

#all down proximal_k4me3loss

#from intg_preaktss6
#all intersect narrow and count

intg_preaktss_all

all_narrow_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "all" , cat.y= "narrow")

all_broad_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "all" , cat.y= "broad")

all_negative_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "all" , cat.y= "negative")

all_narrow<-nrow(all_narrow_gn_peaks)
all_broad<-nrow(all_broad_gn_peaks)
all_neg<-nrow(all_negative_gn_peaks)

nrow(all_narrow_gn_peaks) + nrow(all_broad_gn_peaks) + nrow(all_negative_gn_peaks)


#all down

down_narrow_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "down" , cat.y= "narrow")

down_broad_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "down" , cat.y= "broad")

down_negative_gn_peaks<-category_rapresentation(intg_preaktss6, cat.x = "down" , cat.y= "negative")

down_neg<-nrow(down_negative_gn_peaks)
down_broad<-nrow(down_broad_gn_peaks)
down_narrow<-nrow(down_narrow_gn_peaks)
tdown<-down_neg + down_broad + down_narrow

#c(down_negative_gn_peaks,)

down_distr<-rbind(
  cbind(gene=down_negative_gn_peaks$intg,
        category=rep("negative",nrow(down_negative_gn_peaks))
  ),
  cbind(gene=down_broad_gn_peaks$intg,
        category=rep("broad",nrow(down_broad_gn_peaks))
  ),
  cbind(gene=down_narrow_gn_peaks$intg,
        category=rep("narrow",nrow(down_narrow_gn_peaks))
  )

)

write.table(down_distr,"outs/H3K27me3_broad_narrow/downGenes_peaksType_distribution.tsv",
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE
)

#============================= proximal

proximal_narrow_gn_peaks <- category_rapresentation(intg_preaktss6, cat.x = "proximal_dko_wt_loose" , cat.y= "narrow")
proximal_broad_gn_peaks <- category_rapresentation(intg_preaktss6, cat.x = "proximal_dko_wt_loose" , cat.y= "broad")
proximal_negative_gn_peaks <- category_rapresentation(intg_preaktss6, cat.x = "proximal_dko_wt_loose" , cat.y= "negative")

pro_nar<-nrow(proximal_narrow_gn_peaks)
pro_road<-nrow(proximal_broad_gn_peaks)
pro_neg<-nrow(proximal_negative_gn_peaks)

length(which(intg_preaktss6$V2=="proximal_dko_wt_loose"))

tprox=nrow(proximal_narrow_gn_peaks) + nrow(proximal_broad_gn_peaks) + nrow(proximal_negative_gn_peaks)

df<-rbind(
  data.frame(peak.type=c("narrow","broad","negative"),
             regions=c("all",
                       "all",
                       "all"),
             val=c(all_narrow/(all_narrow+all_broad+all_neg),all_broad/(all_narrow+all_broad+all_neg),all_neg/(all_narrow+all_broad+all_neg))
  ),
  data.frame(peak.type=c("narrow","broad","negative"),
             regions=c("down",
                       "down",
                       "down"),
             val=c(down_narrow/(down_narrow+down_broad+down_neg),down_broad/(down_narrow+down_broad+down_neg),down_neg/(down_narrow+down_broad+down_neg))
  ),
  data.frame(peak.type=c("narrow","broad","negative"),
             regions=c("proximal_dko_wt_loss",
                       "proximal_dko_wt_loss",
                       "proximal_dko_wt_loss"),
             val=c(pro_nar/(pro_nar+pro_road+pro_neg),pro_road/(pro_nar+pro_road+pro_neg),pro_neg/(pro_nar+pro_road+pro_neg))
  )
)

write.table(df, file="outs/H3K27me3_broad_narrow/D0_df_barplot.tsv", sep="\t")

ggplot(data=df, aes(x=regions, y=val, fill=peak.type)) +
  geom_bar(stat="identity") +
  xlab("genes") +
  ylab("Percentage per peak type")

ggsave("./outs/H3K27me3_broad_narrow/peak_distribution.png", width = 10, height = 6, dpi = 300)
