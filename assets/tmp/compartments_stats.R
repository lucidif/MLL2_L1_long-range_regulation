#parameters

#getwd()

source("git/nf-core-microc/bin/compartments_stats.R")

ref<-read.table("outs/2024_10_Lara_microC_downstream/in/mm10.refGene.gtf", sep="\t")

tarcomp<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/PC_selection_micro_Lara.tsv",sep="\t", header=TRUE)

unique(ref$V3)

futr<-ref[which(ref$V3=="5UTR"),]

head(futr)

tss<-c()

tss[which(futr$V7=="+")]<-futr$V4[which(futr$V7=="+")]
tss[which(futr$V7=="-")]<-futr$V5[which(futr$V7=="-")]

futr<-cbind(futr,tss=tss)

unifutr<-futr[,c("V1","tss")]

unifutr<-unifutr[which(!duplicated(unifutr)),]

chromo_list<-read.table("outs/2024_10_Lara_microC_downstream/in/mm10.sizes")
chromo_list<- chromo_list[which(chromo_list[,1] != "chrM"),]
chromo_list<- chromo_list[which(chromo_list[,1] != "chrY"),]

chromo_list<-cbind(chromo_list, chr_comp_pc=rep(NA,nrow(chromo_list)))


pc1_inbed_wt_d0<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_350kb_aLp_WT_day0.Dd.bed")
pc2_inbed_wt_d0<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_PC2_350kb_aLp_WT_day0.Dd.bed")
pc1_inbed_wt_d4<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_350kb_aLp_WT_day4.Dd.bed")
pc2_inbed_wt_d4<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_PC2_350kb_aLp_WT_day4.Dd.bed")

pc1_inbed_ko_d0<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_350kb_aLp_KO_day0.Dd.bed")
pc2_inbed_ko_d0<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_PC2_350kb_aLp_KO_day0.Dd.bed")
pc1_inbed_ko_d4<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_350kb_aLp_KO_day4.Dd.bed")
pc2_inbed_ko_d4<-read.table("outs/2024_10_Lara_microC_downstream/out/fanc_compartments/bychr_PC2_350kb_aLp_KO_day4.Dd.bed")



#make inbed with preferencial pc

bestpc_inbed_wt_d0<-c()

sampl<-c("inbed_wt_d0", 
        "inbed_ko_d0",
        "inbed_wt_d4",
        "inbed_ko_d4"
        )

colnns<-c("PC_d0_wt",
"PC_d0_ko",
"PC_d4_wt",
"PC_d4_ko"
)        

for ( i in 1:length(sampl)){

    pc1_tarbed<-get(paste0("pc1_", sampl[i]))
    pc2_tarbed<-get(paste0("pc2_", sampl[i]))

    pc1chr<-tarcomp[which(tarcomp[,colnns[i]]==1),1]

    pc2chr<-tarcomp[which(tarcomp[,colnns[i]]==2),1]

    #pc1_tarbed[which( pc1_tarbed$V1 %in% pc1chr ), ]   
    #pc2_tarbed[which( pc2_tarbed$V1 %in% pc2chr ), ]

    #pcmerge<-matrix(NA ,ncol=ncol(pc1_tarbed),nrow=nrow(pc2_tarbed))
    pcmerge<-data.frame(V1=rep(NA, nrow(pc1_tarbed)), 
    V2=rep(NA, nrow(pc1_tarbed)),
    V3=rep(NA, nrow(pc1_tarbed)),
    V4=rep(NA, nrow(pc1_tarbed)),
    V5=rep(NA, nrow(pc1_tarbed)),
    V6=rep(NA, nrow(pc1_tarbed))
    )
    #pcmerge<-matrix(NA ,ncol=ncol(pc1_tarbed),nrow=nrow(pc2_tarbed))
    pcmerge[which( pc1_tarbed$V1 %in% pc1chr ),]<-pc1_tarbed[which( pc1_tarbed$V1 %in% pc1chr ), ]
    pcmerge[which( pc2_tarbed$V1 %in% pc2chr ),]<-pc2_tarbed[which( pc2_tarbed$V1 %in% pc2chr ), ]


    #remove NA columns (columns with chromosomes not present in pc list)
    pcmerge<-pcmerge[which(!is.na(pcmerge$V1)),]

    assign(paste0("selectedPC_", sampl[i]),pcmerge)

    

}



# ABgene_wt_d0<-genes.per.compartments(pc1_inbed_wt_d0, unifutr, chromo_list)
# ttABgene_wt_d0<-genes.per.compartments(pc1_inbed_wt_d0, 
# unifutr, 
# chromo_list,
# chromo_directions=c()
# )


ABgene_wt_d0_chr_dir=rep(1,nrow(chromo_list))

ABgene_ko_d0_chr_dir=rep(1,nrow(chromo_list))
ABgene_ko_d0_chr_dir[which(chromo_list$V1=="chr18")]<- -1

ABgene_wt_d4_chr_dir=rep(1,nrow(chromo_list))
ABgene_wt_d4_chr_dir[which(chromo_list$V1=="chr18")]<- -1
ABgene_wt_d4_chr_dir[which(chromo_list$V1=="chr19")]<- -1

ABgene_ko_d4_chr_dir=rep(1,nrow(chromo_list))
ABgene_ko_d4_chr_dir[which(chromo_list$V1=="chr18")]<- -1
ABgene_ko_d4_chr_dir[which(chromo_list$V1=="chr19")]<- -1


ABgene_wt_d0<-genes.per.compartments(pc1_inbed_wt_d0, unifutr, chromo_list)
ABgene_wt_d4<-genes.per.compartments(pc1_inbed_wt_d4, unifutr, chromo_list, ABgene_wt_d4_chr_dir)
ABgene_ko_d0<-genes.per.compartments(pc1_inbed_ko_d0, unifutr, chromo_list, ABgene_ko_d0_chr_dir)
ABgene_ko_d4<-genes.per.compartments(pc1_inbed_ko_d4, unifutr, chromo_list, ABgene_ko_d4_chr_dir)

pc2_ABgene_wt_d0<-genes.per.compartments(pc2_inbed_wt_d0, unifutr, chromo_list)
pc2_ABgene_wt_d4<-genes.per.compartments(pc2_inbed_wt_d4, unifutr, chromo_list, ABgene_wt_d4_chr_dir)
pc2_ABgene_ko_d0<-genes.per.compartments(pc2_inbed_ko_d0, unifutr, chromo_list, ABgene_ko_d0_chr_dir)
pc2_ABgene_ko_d4<-genes.per.compartments(pc2_inbed_ko_d4, unifutr, chromo_list, ABgene_ko_d4_chr_dir)


Agenes_PC_d0_wt<-c()
Agenes_PC_d4_wt<-c()
Agenes_PC_d0_ko<-c()
Agenes_PC_d4_ko<-c()

Bgenes_PC_d0_wt<-c()
Bgenes_PC_d4_wt<-c()
Bgenes_PC_d0_ko<-c()
Bgenes_PC_d4_ko<-c()

d0_ko_selected_pc<-c()
d0_wt_selected_pc<-c()
d4_ko_selected_pc<-c()
d4_wt_selected_pc<-c()


for (i in 1:nrow(tarcomp)){

    icomp<-tarcomp[i,1]

##PC_d0_wt
    if(tarcomp[i,2]==2){ #prendi la pc2 
     d0_wt_selected_pc[i] <- "pc2"
     Agenes_PC_d0_wt[i]<-pc2_ABgene_wt_d0[which(pc2_ABgene_wt_d0$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d0_wt[i]<-pc2_ABgene_wt_d0[which(pc2_ABgene_wt_d0$chr==tarcomp[i,1]),"Bgenes"]  
    }else{ #altrimenti prendi sempre la 1
    d0_wt_selected_pc[i] <- "pc1"
     Agenes_PC_d0_wt[i]<-ABgene_wt_d0[which(ABgene_wt_d0$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d0_wt[i]<-ABgene_wt_d0[which(ABgene_wt_d0$chr==tarcomp[i,1]),"Bgenes"]
    }

##PC_d0_ko
    if(tarcomp[i,2]==2){ #prendi la pc2 
    d0_ko_selected_pc[i] <- "pc2"
    Agenes_PC_d0_ko[i]<-pc2_ABgene_ko_d0[which(pc2_ABgene_ko_d0$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d0_ko[i]<-pc2_ABgene_ko_d0[which(pc2_ABgene_ko_d0$chr==tarcomp[i,1]),"Bgenes"] 
    }else{ #altrimenti prendi sempre la 1
    d0_ko_selected_pc[i] <- "pc1"
    Agenes_PC_d0_ko[i]<-ABgene_ko_d0[which(ABgene_ko_d0$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d0_ko[i]<-ABgene_ko_d0[which(ABgene_ko_d0$chr==tarcomp[i,1]),"Bgenes"]
    }
##PC_d4_wt
    if(tarcomp[i,2]==2){ #prendi la pc2 
    d4_wt_selected_pc[i] <- "pc2"
    Agenes_PC_d4_wt[i]<-pc2_ABgene_wt_d4[which(pc2_ABgene_wt_d4$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d4_wt[i]<-pc2_ABgene_wt_d4[which(pc2_ABgene_wt_d4$chr==tarcomp[i,1]),"Bgenes"]  
    }else{ #altrimenti prendi sempre la 1
    d4_wt_selected_pc[i] <- "pc1"
     Agenes_PC_d4_wt[i]<-ABgene_wt_d4[which(ABgene_wt_d4$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d4_wt[i]<-ABgene_wt_d4[which(ABgene_wt_d4$chr==tarcomp[i,1]),"Bgenes"]
    }
##PC_d4_ko
    if(tarcomp[i,2]==2){ #prendi la pc2 
    d4_ko_selected_pc[i] <- "pc2"
    Agenes_PC_d4_ko[i]<-pc2_ABgene_ko_d4[which(pc2_ABgene_ko_d4$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d4_ko[i]<-pc2_ABgene_ko_d4[which(pc2_ABgene_ko_d4$chr==tarcomp[i,1]),"Bgenes"]
    }else{ #altrimenti prendi sempre la 1
    d4_ko_selected_pc[i] <- "pc1"
    Agenes_PC_d4_ko[i]<-ABgene_ko_d4[which(ABgene_ko_d4$chr==tarcomp[i,1]),"Agenes"]
     Bgenes_PC_d4_ko[i]<-ABgene_ko_d4[which(ABgene_ko_d4$chr==tarcomp[i,1]),"Bgenes"]
    }

}

all_vals_tab<-cbind(tarcomp,
                    Agenes_PC_d0_wt,
                    Bgenes_PC_d0_wt,
                    Agenes_PC_d0_ko,
                    Bgenes_PC_d0_ko,
                    Agenes_PC_d4_wt,
                    Bgenes_PC_d4_wt,
                    Agenes_PC_d4_ko,
                    Bgenes_PC_d4_ko
                    )

selected_pc<-cbind(tarcomp$chr, 
                    d0_wt_selected_pc,
                    d0_ko_selected_pc,
                    d4_wt_selected_pc,
                    d4_ko_selected_pc
                    )

write.table(all_vals_tab,
file="outs/2024_10_Lara_microC_downstream/out/fanc_compartments/ABgenes.tsv",
sep="\t", 
col.names = TRUE, 
row.names = FALSE,
quote=FALSE
)

#invert the direction of eighrnvector

# ABgene_ko_d0_chr_dir=rep(1,nrow(chromo_list))
# ABgene_ko_d0_chr_dir[which(chromo_list$V1=="chr18")]<- -1

# ABgene_wt_d4_chr_dir=rep(1,nrow(chromo_list))
# ABgene_wt_d4_chr_dir[which(chromo_list$V1=="chr18")]<- -1
# ABgene_wt_d4_chr_dir[which(chromo_list$V1=="chr19")]<- -1

# ABgene_ko_d4_chr_dir=rep(1,nrow(chromo_list))
# ABgene_ko_d4_chr_dir[which(chromo_list$V1=="chr18")]<- -1
# ABgene_ko_d4_chr_dir[which(chromo_list$V1=="chr19")]<- -1


dchange_selpc_inbed_ko_d0<- invertCompartments(as.data.frame(selectedPC_inbed_ko_d0),"chr18") 

dchange_selpc_inbed_wt_d4 <- invertCompartments(as.data.frame(selectedPC_inbed_wt_d4),"chr18") 
dchange_selpc_inbed_wt_d4 <- invertCompartments(as.data.frame(dchange_selpc_inbed_wt_d4),"chr19") 

dchange_selpc_inbed_ko_d4 <- invertCompartments(as.data.frame(selectedPC_inbed_ko_d4),"chr18") 
dchange_selpc_inbed_ko_d4 <- invertCompartments(as.data.frame(dchange_selpc_inbed_ko_d4),"chr19") 

compartments.transitions(selectedPC_inbed_wt_d0, 
                        dchange_selpc_inbed_ko_d0,
                        "d0",
                        "outs/2024_10_Lara_microC_downstream/out/fanc_compartments/"
                        )

compartments.transitions(dchange_selpc_inbed_wt_d4, 
                        dchange_selpc_inbed_ko_d4,
                        "d4",
                        "outs/2024_10_Lara_microC_downstream/out/fanc_compartments/"
                        )

#export bed files

write.table(dchange_selpc_inbed_wt_d4, 
file="outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_wt_d4.bed", quote=FALSE, sep="\t",col.names = FALSE, row.names = FALSE)

write.table(dchange_selpc_inbed_ko_d4, 
file="outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_ko_d4.bed", quote=FALSE, sep="\t",col.names = FALSE, row.names = FALSE)

write.table(selectedPC_inbed_wt_d0, 
file="outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_wt_d0.bed", quote=FALSE, sep="\t",col.names = FALSE, row.names = FALSE)

write.table(dchange_selpc_inbed_ko_d4, 
file="outs/2024_10_Lara_microC_downstream/out/fanc_compartments/dchange_ko_d0.bed", quote=FALSE, sep="\t",col.names = FALSE, row.names = FALSE)

