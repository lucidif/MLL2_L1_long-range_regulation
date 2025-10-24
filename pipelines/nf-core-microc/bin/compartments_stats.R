#parameters

genes.per.compartments<-function( pc1_inbed_wt_d0,
unifutr,
chromo_list, chromo_direction=NULL){
  
  
  chr<-c()
  Agenes<-c()
  Bgenes<-c()

if(is.null(chromo_direction)==TRUE){
  chromo_direction<-rep(1,nrow(chromo_list))
}


for (i in 1:nrow(chromo_list)){
   tar_chr<-chromo_list[i,"V1"]
   print(tar_chr) 
   tarfutr<- unifutr[which(unifutr$V1==tar_chr),]
    head(tarfutr)
    comp<-c()
   for (j in 1:nrow(tarfutr)){
      comp[j]<-pc1_inbed_wt_d0[which ( pc1_inbed_wt_d0$V2 <= as.numeric(tarfutr$tss[j]) & 
      pc1_inbed_wt_d0$V3 > as.numeric(tarfutr$tss[j]) & 
      pc1_inbed_wt_d0$V1 == tar_chr ),"V4"]  
   }

   tarfutr<-cbind(tarfutr, comp=comp)
   
    chr[i]<-tar_chr

    if(chromo_direction[i]==1){
    Agenes[i]<-length(which(tarfutr$comp=="A"))
    Bgenes[i]<-length(which(tarfutr$comp=="B"))
    }else{ #chromo_direction[i]!=1

      if(chromo_direction[i]==-1){
        print("inverted direction")
        Agenes[i]<-length(which(tarfutr$comp=="B"))
        Bgenes[i]<-length(which(tarfutr$comp=="A"))}
    }




}

return(data.frame(chr,Agenes,Bgenes))

}


compartments.transitions<-function(inbed1, 
                                    inbed2, 
                                    name="", 
                                    outfolder=""
                                    )
                                    {

getwd()
  #inbed1<-read.table("out/fanc_compartments/1M_aLp_WT_day0.Dd.bed")
  #inbed2<-read.table("out/fanc_compartments/1M_aLp_KO_day0.Dd.bed")
  #name="1M_WTd0_to_KOd0"
  #outfolder="out/fanc_compartments/"

  #libraries

  library(ggplot2)


  #code

  if(nrow(inbed1)==nrow(inbed2)){
      total_bins=nrow(inbed1)
  }

  AAbin<-length(which(inbed1$V4=="A" & inbed2$V4=="A"))
  BBbin<-length(which(inbed1$V4=="B" & inbed2$V4=="B"))

  ABbin<-length(which(inbed1$V4=="A" & inbed2$V4=="B"))
  BAbin<-length(which(inbed1$V4=="B" & inbed2$V4=="A"))


  paste0("bin A->A = ", AAbin)
  paste0("bin B->B = ", BBbin)
  paste0("bin A->B = ", ABbin)
  paste0("bin B->A = ", BAbin)

  paste0("AAbin + BBbin + ABbin + BAbin = ", AAbin + BBbin + ABbin + BAbin)
  paste0("total_bins = ", total_bins)
  paste0("AAbin + BBbin + ABbin + BAbin == total_bins -> ",AAbin + BBbin + ABbin + BAbin == total_bins)

  df<-as.data.frame(cbind(comp=c("AA","BB","AB","BA"),rbind(AAbin,BBbin,ABbin,BAbin)))

  rownames(df)=NULL

  df_withtot<-rbind(df,cbind(comp="tot",V2=total_bins))
  write.table(df_withtot, file=paste0(outfolder,"/bin_AB_changes_",name,".tsv"), sep="\t",quote=FALSE)


  #make plot



  df_perc<-df
  df_perc$V2<-as.numeric(as.character(df$V2))/total_bins

  df_perc$comp <- factor(df_perc$comp,
                        levels = c("AA", "BB", "AB", "BA"),
                        labels = c("A -> A", "B -> B", "A -> B", "B -> A"))

  p <- ggplot(df_perc, aes(x=comp, y=V2, fill=comp)) +
    geom_bar(stat="identity", width=0.6) +
    geom_text(aes(label=sprintf("%.1f%%", V2*100)), vjust=-0.5, size=5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    labs(title=paste("changes of compartments :", name),
        x="Transition type",
        y="Percentage of total bins",
        fill="Transition",
        face = "bold"
        ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "right"
    )


  p  

  ggsave(
    filename = paste0(outfolder, "plot_bin_AB_changes_", name, ".png"),
    plot = p,
    width = 6,
    height = 4,
    dpi = 300
  )

}


invertCompartments<-function(inbed, chr){

invbed<-inbed
invbed[which(invbed$V1==chr & invbed$V4=="A"),4]<-1
invbed[which(invbed$V1==chr & invbed$V4=="B"),4]<-2
invbed[which(invbed$V1==chr & invbed$V4==1),4]<-"B"
invbed[which(invbed$V1==chr & invbed$V4==2),4]<-"A"

#invert sign 
invbed[which(invbed$V1==chr ),5]<-invbed[which(invbed$V1==chr ) ,5]*(-1)


return (invbed)

}


