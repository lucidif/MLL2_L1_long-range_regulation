#==========================================
#              load functions
#==========================================
source ("git/nf-core-microc/bin/microC_MAplot_fun.R")

#==========================================
#               set env
#==========================================

hic_path="/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/2024_10_Lara_microC_downstream/out/pairtools_merge/hic_by_replicate"

hic_name<-c("aLp_KO_day0_A.Dd.hic",
"aLp_KO_day0_B.Dd.hic", 
"aLp_KO_day4_A.Dd.hic", 
"aLp_KO_day4_B.Dd.hic", 
"aLp_WT_day0_A.Dd.hic", 
"aLp_WT_day0_B.Dd.hic",
"aLp_WT_day4_A.Dd.hic",
"aLp_WT_day4_B.Dd.hic"
)

#aLp_WT_day4_A.Dd.hic


bedpe.path<-"outs/Lara_multiomic_analysis/outs/coolpup/500bp/win500_anchors3.bedpe"
ach3<-read.table(bedpe.path, sep="\t")
chrsize.path<-"/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/2024_10_Lara_microC_downstream/in/mm10.sizes"

chromsize<-read.table(chrsize.path)
bin.size=5000


library(strawr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)

#==========================================
#          Execute analysis 
#==========================================


for (i in 1:length(hic_name)){

    print("extract target contacts from")
    print(hic_name[i])

    bedpe<-read.table(bedpe.path)

    res_list <- vector("list", nrow(ach3))  # pre-alloc
    validr<-0
    out <- data.frame(
        anchor = character(nrow(bedpe)),
        bait   = character(nrow(bedpe)),
        value  = rep(NA_real_, nrow(bedpe)),
        stringsAsFactors = FALSE
    )
    for (j in 1:nrow(bedpe)){
        
        anch <- bedpe[j, c(1,2,3)]
        bait <- bedpe[j, c(4,5,6)]

        delta.centered.anch <- ((anch[,3]-anch[,2])/2)
        delta.centered.bait <- ((bait[,3]-bait[,2])/2)

        anch.centered<- cbind(anch[,1], anch[,2]+delta.centered.anch, anch[,2]+delta.centered.anch+1)
        bait.centered<- cbind(bait[,1], bait[,2]+delta.centered.bait, bait[,2]+delta.centered.bait+1)

        imported_hic <- straw(
            norm    = "KR",
            fname   = paste0(hic_path, "/", hic_name[i]),
            chr1loc = paste(anch, collapse=":"),
            chr2loc = paste(bait, collapse=":"),
            unit    = "BP",
            binsize = bin.size
        )
        out$anchor[j] <- paste(anch.centered, collapse=":")
        out$bait[j]   <- paste(bait.centered, collapse=":")
        out$value[j]  <- if (nrow(imported_hic) > 0) sum(imported_hic$counts) else NA_real_

        if (nrow(imported_hic) > 0) {
            validr <- validr + 1
            res_list[[j]] <- list(
            i = j,
            anch = anch,
            bait = bait,
            contacts = imported_hic
            )
        } else {
            res_list[[j]] <- NULL
        }
    }
    assign(paste0(hic_name[i], "_central_bin"),out)

    sxcontact<-extract_contacts(hicpath=paste0(hic_path,"/",hic_name[i]), 
                      bedpe.path=bedpe.path,
                      window_bins=5,
                      viewpoint="center",
                      bin.size=5000,
                      chromsize.file="/home/lucio/wkdir/projects/MLL2_L1_regulation/outs/2024_10_Lara_microC_downstream/in/mm10.sizes",
                      return_mode="sum_df"
     )

    assign(hic_name[i],sxcontact)

}


allsamples_tb<-cbind(get(hic_name[1]),
get(hic_name[2])[,3],
get(hic_name[3])[,3],
get(hic_name[4])[,3],
get(hic_name[5])[,3],
get(hic_name[6])[,3],
get(hic_name[7])[,3],
get(hic_name[8])[,3]
)

colnames(allsamples_tb)<-c("anchor","bait",hic_name[1],
hic_name[2],
hic_name[3],
hic_name[4],
hic_name[5],
hic_name[6],
hic_name[7],
hic_name[8]
)

write.table(allsamples_tb, 
file="outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/taget_regions_contact_matrix_sum.tsv", 
col.names = TRUE,
sep="\t",
quote=FALSE,
row.names=FALSE
)


allsamples_tb_central<-cbind(get(paste0(hic_name[1], "_central_bin")),
get(paste0(hic_name[2], "_central_bin"))[,3],
get(paste0(hic_name[3], "_central_bin"))[,3],
get(paste0(hic_name[4], "_central_bin"))[,3],
get(paste0(hic_name[5], "_central_bin"))[,3],
get(paste0(hic_name[6], "_central_bin"))[,3],
get(paste0(hic_name[7], "_central_bin"))[,3],
get(paste0(hic_name[8], "_central_bin"))[,3]
)


colnames(allsamples_tb_central)<-c("anchor","bait",hic_name[1],
hic_name[2],
hic_name[3],
hic_name[4],
hic_name[5],
hic_name[6],
hic_name[7],
hic_name[8]
)

#come ti spieghi tutti questu NA ?

write.table(allsamples_tb_central, 
file="outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/taget_regions_central_bin.tsv", 
col.names = TRUE,
sep="\t",
quote=FALSE,
row.names=FALSE
)

long_df <- allsamples_tb |>
  pivot_longer(
    cols = -c(anchor, bait),
    names_to = "sample",
    values_to = "value"
  )

long_df2 <- long_df %>%
  mutate(
    geno = if_else(str_detect(sample, "_KO_"), "KO", "WT"),
    day  = if_else(str_detect(sample, "day0"), "day0", "day4"),
    group = paste(geno, day, sep = "_")
  )

p<-ggplot(long_df2, aes(x = sample, y = value, fill = group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "KO_day0" = "#d73027",  # rosso
    "KO_day4" = "#fc8d59",  # rosso più chiaro/arancio
    "WT_day0" = "#1f78b4",  # blu
    "WT_day4" = "#a6cee3"   # blu più chiaro
  )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = NULL,
    y = "Sum(KR counts) [log10]",
    fill = "Group",
    title = "Contacts per sample (log scale)"
  )

ggsave("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/contacts_per_sample_log_boxplot.pdf", plot = p, width = 10, height = 6, units = "in")

save.image(file = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/microC_pileupPlot_fromhic.RData")
