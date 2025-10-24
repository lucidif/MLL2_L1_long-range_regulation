extract_distances<-function(peak.coords, startend.coords,
                            from=c("peak","startend"),
                            to=c("startend","peak"),
                            peak.coordtype=c("centered","fprime"),
                            startend.coordtype=c("fprime","centered"),
                            to_name="to",
                            from_name="from",
                            window=NULL,
                            window.type=c("inside","outside")
                            #only.nearest.to=FALSE #if FALSE the argorithm calculate the distance from each nearest from coordniate by to coordinate; but in this case if to coordinates is near more from coordinates , it mantains all distance from the same to coordinate.  
){
  
  #startend.coords=degs_D0_ann_diff_down[,c(3,4,5,2,7,6)]
  
  #from mean coords to startend.coords
  #define what is the from (center of calculation) and witch is to (the element that you what know the sitances)
  
  find.fprime<-function(coordinates){
    fprime<-c()
    for (i in 1:nrow(coordinates)){
      if(coordinates[i,6]=="-"){
        fprime[i]<-coordinates[i,3]
      }else{
        
        if(coordinates[i,6]=="+"){
          fprime[i]<-coordinates[i,2]
        }else{
          #fprime[i]<-coordinates[i,4]+((coordinates[i,5]-coordinates[i,4])/2)
          fprime[i]<-NA
        }
        
        
      }
      
    }
    
    num.fprime <- as.numeric(as.character(fprime))
    
    return(data.frame(chr=coordinates[,1],
                      fprime=num.fprime, 
                      identifier=coordinates[,4])
    )
    #num.mean.tarreg<- as.numeric(as.character(mean_tarreg_coordinates))
  }
  
  
  
  centered.peak<-data.frame(chr=peak.coords[,1],
                            mean.crd=((peak.coords[,3]-peak.coords[,2])/2)+peak.coords[,2],
                            identifier=peak.coords[,4]
  )
  
  centered.starend<-data.frame(chr=startend.coords[,1],
                               mean.crd=((startend.coords[,3]-startend.coords[,2])/2)+startend.coords[,2],
                               identifier=startend.coords[,4]
  )                             
  
  fprime.peak <- find.fprime(coordinates=peak.coords) 
  fprime.startend <- find.fprime(coordinates=startend.coords)
  
  #set the from
  
  if(from[1]=="peak" & peak.coordtype[1] == "centered"){
    from.coordinates<-centered.peak
    #to.coordinates<-startend.coords
    
    #tarreg<-data.frame(chr=from.coordinates[,1],mean.crd=(from.coordinates[,3]-from.coordinates[,2]/2)+from.coordinates[,2])
  }
  
  if(from[1]=="peak" & peak.coordtype[1] == "fprime"){
    from.coordinates<-fprime.peak
  }
  
  if(from[1]=="startend" & startend.coordtype[1] == "centered"){
    from.coordinates<-centered.starend
  }
  
  if(from[1]=="startend" & startend.coordtype[1] == "fprime"){
    from.coordinates<-fprime.startend
  }
  
  if(to[1]=="peak" & peak.coordtype[1] == "centered"){
    to.coordinates<-centered.peak
  }
  
  if(to[1]=="peak" & peak.coordtype[1] == "fprime"){
    to.coordinates<-fprime.peak
  }
  
  if(to[1]=="startend" & startend.coordtype[1] == "fprime"){
    to.coordinates<-fprime.startend
  }
  
  if(to[1]=="startend" & startend.coordtype[1] == "centered"){
    to.coordinates<-centered.starend
  }
  
  
  #set the to
  
  #calculate new coordinates
  #tarreg<-peak.coords[j,] 
  #l1mask.mc.sub<-startend.coords[which(startend.coords[,3]==as.character(tarreg[1])),]
  #mean_tarreg_coordinates<-((tarreg[1,3]-tarreg[1,2])/2)+tarreg[1,2]
  distances<-c()
  to_chr<-c()
  to_position<-c()
  to_identifier<-c()
  
  from_chr<-c()
  from_position<-c()
  from_identifier<-c()
  
  
  for(j in 1:nrow(from.coordinates)){
    if(j %% 1000 == 0){
      cat("j: ",j,"\n")
    }
    
    tarreg<-from.coordinates[j,]
    #l1mask.mc.sub<-startend.coords[which(startend.coords[,3]==as.character(tarreg[1])),]
    to.reg<-to.coordinates[which(to.coordinates[,1]==as.character(tarreg[1])),]
    
    if (nrow(to.reg) != 0 ) {
      allvals<-abs(to.reg[,2] - tarreg[,2])
      
      #to_ensembl[j]<-to.reg[which(allvals==min(allvals)),1]
      #to_gene[j]<-to.reg[which(allvals==min(allvals)),2]
      to_chr[j]<-to.reg[which(allvals==min(allvals)),1]
      to_position[j]<-to.reg[which(allvals==min(allvals)),2]
      to_identifier[j]<-to.reg[which(allvals==min(allvals)),3]
      #to_end[j]<-to.reg[which(allvals==min(allvals)),5]
      #to_strand[j]<-to.reg[which(allvals==min(allvals)),6]
      #to_fprime[j]<-
      
      distances[j]<-min(allvals)
      
      from_chr[j]<-tarreg[1,1]
      from_position[j]<-tarreg[1,2]
      from_identifier[j]<-tarreg[1,3]
      #from_end[j]<-tarreg[1,3] 
    } 
    
    
  }
  
  
  from.nearest<-data.frame(
    from_chr,
    from_position,
    #from_end,
    to_chr,
    to_position,
    #to_end,
    #to_strand,
    distances,
    from_identifier,
    to_identifier
  )
  
  #remove the raw that are FALSE in if (nrow(to.reg) != 0 )
  from.nearest<-na.omit(from.nearest)
  
  #bring only the minimum distance region interaction when from identifier is the same
  dup.from<-names(which(table(from.nearest$from_identifier)>1))
  if(length(dup.from)>0){
    
    for ( z in 1:length(dup.from) ){
      dup.nearest <- from.nearest[which(from.nearest$from_identifier==dup.from[z]),]
      dup.selected.add <- dup.nearest[which(dup.nearest$distances==min(dup.nearest$distances)),]
      if(z ==1 ){
        dup.selected<-dup.selected.add
      }else{
        dup.selected<-rbind(dup.selected, dup.selected.add)
      }
    }
    
    less.dup<-from.nearest[!from.nearest$from_identifier %in% dup.from  , ]
    from.nearest<-rbind(less.dup,dup.selected)
    
  }
  
  
  
  
  # if(only.nearest.to==TRUE){
  #   
  #   chrtr<-unique(from.nearest$from_chr)
  #   
  #   for(i in 1:length(chrtr)){
  #     tar.chr<-chrtr[i]
  #     tar.nearest<-from.nearest[which(from.nearest$from_chr==tar.chr & from.nearest$to_chr== tar.chr ),]
  #     
  #     # 1 order by distance
  #     head(tar.nearest[order(tar.nearest$distances),])
  #     
  #     
  #     
  #   }
  #   
  #   
  # }
  
  if(is.null(window)==FALSE){
    #IF window is not null do :
    if(window.type[1]=="inside"){
      from.nearest <-from.nearest[which(from.nearest$distances<=window),]
    }
    if(window.type[1]=="outside"){
      from.nearest <-from.nearest[which(from.nearest$distances>=window),]
    }    
    
  }
  
  return(
    from.nearest
  )
  
  
}




addAnnotationToDegs<-function(degs, gene.annWithPosition, padj=0.05, log2fc=1, 
                              direction=c("all","up","down")
){
  
  degs_annotation <- merge(gene.annWithPosition ,degs, by=1)
  
  
  #degs_annotation$Chromosome.scaffold.name<-paste0("chr",degs_annotation$Chromosome.scaffold.name)
  
  #degs_annotation[which(degs_annotation$Strand=="-1"),"Strand"]<-"-"
  #degs_annotation[which(degs_annotation$Strand=="1"),"Strand"]<-"+"
  #=============================================================================
  
  degs_ann_diff<-degs_annotation[which(degs_annotation$padj<=padj),]
  
  if(direction[1]=="all"){
    degs_ann_diff_down<-degs_ann_diff[which(degs_ann_diff$log2FoldChange <= - abs(log2fc) | degs_ann_diff$log2FoldChange >= abs(log2fc)
    ),]
  }
  
  if(direction[1]=="down"){
    degs_ann_diff_down<-degs_ann_diff[which(degs_ann_diff$log2FoldChange <= log2fc
    ),]
  }
  
  if(direction[1]=="up"){
    degs_ann_diff_down<-degs_ann_diff[which(degs_ann_diff$log2FoldChange >= log2fc
    ),]
  } 
  
  degs_ann_diff_down<-degs_ann_diff_down[,c("GeneID",
                                            "GeneName",
                                            "Chromosome",
                                            "start",
                                            "end",
                                            "strand",
                                            "log2FoldChange",
                                            "padj"
  )]
  
  colnames(degs_ann_diff_down)<-c("ensembl", "gene", "chromosome","start","end","strand","log2FoldChange", "padj")
  
  return(degs_ann_diff_down)
  
}


distance_plot<-function(minus_peaks, 
                        plus_peaks,
                        degs, 
                        title="distance plot",
                        window=NULL,
                        window.type=c("inside","outside")
){
  
  
  # head(minus_peaks)
  #    V1      V2      V3     V4  V5 V6
  # 1 chr1 3851928 3853587 peak_1  75  .
  # 2 chr1 5076717 5077203 peak_2 162  .
  # 3 chr1 5188428 5190006 peak_3  74  .
  # 4 chr1 5349541 5351436 peak_4  96  .
  # 5 chr1 5992727 5994069 peak_5 303  .
  # 6 chr1 6150526 6152734 peak_6 196  .
  
  
  #dko.dis.CpG.minus.peak<-peaks.near.lines
  
  # minus_peaks<-minus_peaks
  # plus_peaks<-plus_peaks
  # degs<-degs
  degs_reform<-degs[,c(3,4,5,2,7,6)]
  
  
  # down
  dis.CpG.minus.degs_D0<-extract_distances(peak.coords = minus_peaks,
                                           startend.coords = degs_reform,
                                           from="startend",
                                           to="peak",
                                           peak.coordtype = "centered",
                                           startend.coordtype = "fprime",
                                           to_name="to",
                                           from_name="from",
                                           window=window,
                                           window.type=window.type
  )
  
  #plot(density(dis.CpG.minus.degs_D0$distances))
  summary(dis.CpG.minus.degs_D0$distances)
  
  dis.CpG.plus.degs_D0<-extract_distances(peak.coords = plus_peaks,
                                          startend.coords = degs_reform,
                                          from="startend",
                                          to="peak",
                                          peak.coordtype = "centered",
                                          startend.coordtype = "fprime",
                                          to_name="to",
                                          from_name="from",
                                          window=window,
                                          window.type=window.type
  )
  
  summary(dis.CpG.plus.degs_D0$distances)
  
  sign<-t.test( dis.CpG.plus.degs_D0$distances ,  dis.CpG.minus.degs_D0$distances, 
                alternative = "two.sided", var.equal = FALSE)
  
  #sign<-wilcox.test(dis.CpG.plus.degs_D0$distances ~ dis.CpG.minus.degs_D0$distances)
  # sign<-wilcox.test(x = dis.CpG.plus.degs_D0$distances, y = dis.CpG.minus.degs_D0$distances,
  #   mu=0, alt="two.sided", conf.int=T, conf.level=0.99,
  #   paired=FALSE, exact=T, correct=T)
  # distance between cdfs
  
  dis.CpG.plus.degs_D0<-dis.CpG.plus.degs_D0[which(is.na(dis.CpG.plus.degs_D0$distances)==FALSE),]
  dis.CpG.minus.degs_D0<-dis.CpG.minus.degs_D0[which(is.na(dis.CpG.minus.degs_D0$distances)==FALSE),]
  #sign<-ks.test(x = dis.CpG.plus.degs_D0$distances, y = dis.CpG.minus.degs_D0$distances)
  
  print(sign)
  
  plot(density(dis.CpG.plus.degs_D0$distances), main=paste0(title,
                                                            " KStest pval = ",
                                                            sign$p.value,
                                                            " D = ",
                                                            round(sign$statistic,3)
  )
  ) 
  lines(density(dis.CpG.minus.degs_D0$distances), col = "red") 
  legend("topright", c("Plus", "Minus"), 
         col =c("black","red"), lty=1)
  
  return(list(dis.CpG.plus.degs=dis.CpG.plus.degs_D0,dis.CpG.minus.degs=dis.CpG.minus.degs_D0))
  
}

distance_plot.one.sample<-function(minus_peaks, 
                        plus_peaks,
                        degs, 
                        title="distance plot",
                        window=NULL,
                        window.type=c("inside","outside")
){
  
  
  # head(minus_peaks)
  #    V1      V2      V3     V4  V5 V6
  # 1 chr1 3851928 3853587 peak_1  75  .
  # 2 chr1 5076717 5077203 peak_2 162  .
  # 3 chr1 5188428 5190006 peak_3  74  .
  # 4 chr1 5349541 5351436 peak_4  96  .
  # 5 chr1 5992727 5994069 peak_5 303  .
  # 6 chr1 6150526 6152734 peak_6 196  .
  
  
  #dko.dis.CpG.minus.peak<-peaks.near.lines
  
  # minus_peaks<-minus_peaks
  # plus_peaks<-plus_peaks
  # degs<-degs
  degs_reform<-degs[,c(3,4,5,2,7,6)]
  
  
  # down
  dis.CpG.minus.degs_D0<-extract_distances(peak.coords = minus_peaks,
                                           startend.coords = degs_reform,
                                           from="startend",
                                           to="peak",
                                           peak.coordtype = "centered",
                                           startend.coordtype = "fprime",
                                           to_name="to",
                                           from_name="from",
                                           window=window,
                                           window.type=window.type
  )
  
  #plot(density(dis.CpG.minus.degs_D0$distances))
  summary(dis.CpG.minus.degs_D0$distances)
  
  # dis.CpG.plus.degs_D0<-extract_distances(peak.coords = plus_peaks,
  #                                         startend.coords = degs_reform,
  #                                         from="startend",
  #                                         to="peak",
  #                                         peak.coordtype = "centered",
  #                                         startend.coordtype = "fprime",
  #                                         to_name="to",
  #                                         from_name="from",
  #                                         window=window,
  #                                         window.type=window.type
  # )
  
  # summary(dis.CpG.plus.degs_D0$distances)
  
  # sign<-t.test( dis.CpG.plus.degs_D0$distances ,  dis.CpG.minus.degs_D0$distances, 
  #               alternative = "two.sided", var.equal = FALSE)
  
  #sign<-wilcox.test(dis.CpG.plus.degs_D0$distances ~ dis.CpG.minus.degs_D0$distances)
  # sign<-wilcox.test(x = dis.CpG.plus.degs_D0$distances, y = dis.CpG.minus.degs_D0$distances,
  #   mu=0, alt="two.sided", conf.int=T, conf.level=0.99,
  #   paired=FALSE, exact=T, correct=T)
  # distance between cdfs
  
  # dis.CpG.plus.degs_D0<-dis.CpG.plus.degs_D0[which(is.na(dis.CpG.plus.degs_D0$distances)==FALSE),]
  dis.CpG.minus.degs_D0<-dis.CpG.minus.degs_D0[which(is.na(dis.CpG.minus.degs_D0$distances)==FALSE),]
  #sign<-ks.test(x = dis.CpG.plus.degs_D0$distances, y = dis.CpG.minus.degs_D0$distances)
  
  # print(sign)
  
  plot(density(dis.CpG.minus.degs_D0$distances), main=paste0(title,
                                                            " KStest pval = ",
                                                            sign$p.value,
                                                            " D = ",
                                                            round(sign$statistic,3)
  )
  ) 
  
  return(dis.CpG.minus.degs_D0)
  
}

extractdistances<-function(peaks.coords, repeat_mask_line, 
                           addfamily=FALSE, 
                           addpeaks=FALSE, 
                           add5prime=FALSE, 
                           addlines=FALSE){
  distances<-c()
  family<-NULL
  peaks_chr<-c()
  peaks_start<-c()
  peaks_end<-c()
  line_chr<-c()
  line_start<-c()
  line_end<-c()
  line_strand<-c()
  fprime<-c()
  mean.coords<-cbind(peaks.coords[,1],peaks.coords[,2]+round((peaks.coords[,3]-peaks.coords[,2])/2))
  
  repeat_mask_line<-repeat_mask_line[,c(6,7,8,10,11,12,13)]
  
  colnames(repeat_mask_line)<-c("chr","start","end","strand","family","element","type")
  
  repeat_mask_line_plus<-repeat_mask_line[which(repeat_mask_line$strand=="+"),]
  repeat_mask_line_minus<-repeat_mask_line[which(repeat_mask_line$strand=="-"),]
  
  print("nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line)")
  print(nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line))
  
  print("nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line)")
  print(nrow(repeat_mask_line_plus)+nrow(repeat_mask_line_minus)==nrow(repeat_mask_line))
  
  repeat_mask.5prime<-rbind(cbind(repeat_mask_line_plus, fprime=repeat_mask_line_plus$start),
                            cbind(repeat_mask_line_minus, fprime=repeat_mask_line_minus$end))
  
  #repeat_mask.coordinates<-cbind(repeat_mask.5prime$chr,repeat_mask.5prime$fprime)
  #repeat_mask.coordinates_extended=cbind(repeat_mask.5prime$chr,
  #                                       repeat_mask.5prime$fprime,
  #                                       repeat_mask.5prime$family
  #)
  
  l1mask.meancoords<-repeat_mask.5prime
  removed_lines<-0
  for(j in 1:nrow(mean.coords)){
    if(j %% 1000 == 0){
      cat("j: ",j,"\n")
    }
    tarreg<-mean.coords[j,]
    l1mask.mc.sub<-l1mask.meancoords[which(l1mask.meancoords$chr==tarreg[1]),]
    allvals<-abs(as.numeric(as.character(l1mask.mc.sub$fprime))
                 -as.numeric(as.character(tarreg[2])))
    
    distances[j]<-min(allvals)
    
    if(addfamily==TRUE){
      family[j]<-l1mask.mc.sub[which(allvals==distances[j])[1],"family"]
    }
    
    if(add5prime==TRUE){

      
      if(length(allvals)!=nrow(l1mask.mc.sub)){
        print(j)
        print(length(allvals))
        print(dim(l1mask.mc.sub))
      }
      
      if(length(which(allvals==min(allvals)))>1){
        print(j)
        print(l1mask.mc.sub[which(allvals==min(allvals)),])
        print("more annotated line with same five prime coordinate (same minimal distance by peak), only the first annotated is considered:")
        l1mask.mc.sub[1,]<-l1mask.mc.sub[which(allvals==min(allvals))[1],]
        print(l1mask.mc.sub[1,])
        removed_lines<-removed_lines+1
      }

      fprime[j]<-l1mask.mc.sub[which(allvals==min(allvals)),"fprime"]
    }
    
    if(addpeaks==TRUE){
      peaks_chr[j]<-peaks.coords[j,1]
      peaks_start[j]<-peaks.coords[j,2]
      peaks_end[j]<-peaks.coords[j,3]
    }  
    
    if(addlines==TRUE){
      line_chr[j]<-l1mask.mc.sub[which(allvals==min(allvals)),"chr"]
      line_start[j]<-l1mask.mc.sub[which(allvals==min(allvals)),"start"]
      line_end[j]<-l1mask.mc.sub[which(allvals==min(allvals)),"end"]
      line_strand[j]<-l1mask.mc.sub[which(allvals==min(allvals)),"strand"]
    }
  }
  
  print(paste0("removed lines: ",
               removed_lines,
               " of ",  
               nrow(mean.coords) , 
               " totals")
        )
  
  if(addpeaks==TRUE){
    distances<-cbind(distances,peaks_chr,peaks_start,peaks_end)
  }
  
  if(addfamily==TRUE){
    distances<-cbind(distances, family) 
  }
  
  if(addlines==TRUE){
    distances<-cbind(distances, line_chr, line_start, line_end, line_strand)
  }
  
  if(add5prime==TRUE){
    distances<-cbind(distances,fprime=fprime)
  }
  
  
  
  return(distances) 
  
}

filterbydistance<-function(extractdistances.out, fildist=2500, out="/outs/distfilter"){
  
  K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed=extractdistances.out
  
  peaks_cdr<- data.frame(chr=K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"peaks_chr"],
                         start=as.numeric(as.character(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"peaks_start"])),
                         end=as.numeric(as.character(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"peaks_end"]))
  )
  
  #plot(density(peaks_cdr$end - peaks_cdr$start))
  
  #fildist<-2500
  
  subdist<-K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[which(as.numeric(as.character(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[,"distances"]))<=fildist),]
  
  subdist.tarfam_anno<-subdist[,c("line_chr","line_start","line_end","family","distances","line_strand")]
  
  familyfreq<-as.data.frame(table(subdist.tarfam_anno[,"family"]))
  
  subdist.tarfam_anno[,"family"]<-paste0(subdist.tarfam_anno[,"family"],"_",
                                         rep(1:nrow(subdist.tarfam_anno)))
  
  
  print(paste0("at distance <=" ,
               fildist ,
               " retained ",
               nrow(subdist),
               "/",
               nrow(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed),
               " of total peaks"
  ))
  
  
  subdist.tarpeak<-subdist[,c("peaks_chr","peaks_start","peaks_end")]
  subdist.tarpeak<-cbind(subdist.tarpeak,
                         name=paste0("peak_",rep(1:nrow(subdist.tarpeak))),
                         score=subdist[,"distances"],
                         strand=rep(".",nrow(subdist.tarpeak))
  )
  
  
  familyfreq<-cbind(familyfreq,group=rep("distance_filtered"))
  
  familyfreq<-familyfreq[order(familyfreq$Freq, decreasing = TRUE),]
  familyfreq2<-familyfreq
  
  levels(familyfreq2$Var1) <- c(levels(familyfreq2$Var1), "Other")
  familyfreq2[15:nrow(familyfreq),"Var1"]<-"Other"
  
  ggplot(familyfreq2, aes(fill=Var1, y=Freq, x=group)) + 
    geom_bar(position="fill", stat="identity")
  
  
  
  write.table(subdist.tarfam_anno, file=paste0(out,".l1.bed"), 
              col.names = FALSE, 
              row.names = FALSE,
              quote= FALSE,
              sep="\t"
  )
  
  write.table(subdist.tarpeak, file=paste0(out,".target.peaks.bed"),
              col.names = FALSE, 
              row.names = FALSE,
              quote= FALSE,
              sep="\t"
  )
}

