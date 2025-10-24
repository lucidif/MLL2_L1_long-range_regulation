extract_degs<-function(degs, log2FC=1, direction=c("all","up","down"),pval=0.05){
  
  #set pval to null to obtain only unchanged genes
  
  if(is.null(pval)==TRUE){
    degs_sign<-degs[which( degs$padj > 0.05 ),]
  } else {
    degs_sign<-degs[which( degs$padj <= pval ),]
    if(direction[1] == "all"){
      degs_sign<-degs_sign[which( 
        degs_sign$log2FoldChange <= - abs(log2FC)
        | degs_sign$log2FoldChange >= log2FC ),]
    }
    
    if(direction[1] =="up"){
      degs_sign<-degs_sign[which( degs_sign$log2FoldChange >= log2FC ),]
    }
    
    if(direction[1] == "down"){
      degs_sign<-degs_sign[which( degs_sign$log2FoldChange <= - abs(log2FC) ),]
    }
  }
  
    
      
    
    return(degs_sign)
    
  #sign_int <- degs_sign[which( degs_sign$gene_id %in% int.ensembl$V2.y ),]
  
  #sign_int.up <- degs_sign[which( degs_sign.up$gene_id %in% int.ensembl$V2.y ),]
  #sign_int.down <- degs_sign[which( degs_sign.down$gene_id %in% int.ensembl$V2.y ),]
  
  #d0d4.down<-intersect( sign_int_d0.down$gene_id , sign_int_d4.down$gene_id)
  
  #sign_int_d4d0.down<-merge(sign_int_d0.down, sign_int_d4.down, by=1)
}



hyper_degs_tgenes<-function(degs, targetGenes, anno, lower.tail = FALSE, 
                            extract=c("hypergeom", "overlap"), 
                            degs.col=1, 
                            anno.col=1,
                            target.col=1 )
                            {
  
  out.universe<-degs[which(! degs[,degs.col] %in% anno[,anno.col]),degs.col]
  uni.degs<-degs[which( degs[,degs.col] %in% anno[,anno.col]),]
  #uni.degs<-unique(degs[, degs.col],degs[, anno.col])
  
  overlap<-intersect(targetGenes[,target.col], uni.degs[,2])
  #overlap<-intersect(targetGenes[,target.col], uni.degs)
  
  #n.overlap<-length(intersect(targetGenes[,1], uni.degs[,2]))
  #remove outuniverse overlapping genes
  n.overlap<-length(setdiff(overlap,out.universe))
  
  pval<-phyper(
    q = n.overlap - 1,
    m = nrow(targetGenes),
    n = nrow(anno) - nrow(targetGenes),
    k = nrow(uni.degs),
    lower.tail = FALSE
  )
  
  if(extract[1]=="hypergeom"){
    return(
      c(q=as.character(n.overlap - 1),
        m=as.character(nrow(targetGenes)),
        n=as.character(nrow(anno) - nrow(targetGenes)),
        k = as.character(nrow(uni.degs)),
        lower.tail = as.character(lower.tail) ,
        pval = as.character(pval),
        all.anno=nrow(anno),
        targetGenes=nrow(targetGenes),
        out.universe=length(out.universe)
      )
    )
  }else{
    return(overlap)
  }
  

  
}

