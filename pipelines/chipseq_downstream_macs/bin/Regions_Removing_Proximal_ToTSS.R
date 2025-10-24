######################################################################################
## Function to remove genome regions if they are closer than a threshold to a TSS 
## and returning the rows selected by distance (those regions further away from TSS)
######################################################################################

##creator : Victor Sanchez Gaya
##edited : Lucio Di Filippo

regions_removingProximal<-function(regions_bedfile, TSS_bedfile, kbDistance, startend_mean = TRUE){
  print("Check that chr in both files present the same namings. That the characters are not factor codified and that positions are in numeric class.")
  print("It is required a bed file with the enhancer coordinates. With column names: chr, start, end")
  print("The TSS anotation with column names: chr TSS gene")
  
  #regions_bedfile = "./otherouts/deeptools_heatmaps/proximal_NMIs_preaks.bed"
  #TSS_bedfile = "./otherouts/deeptools_heatmaps/NMIs_mm10_mESC_afterLiftover_from_mm9_to_mm10.bed"
  #kbDistance=0
  
  regionsBed<-read.table(regions_bedfile)
  TSS_in<-read.table(TSS_bedfile)

  TSS_loc<-data.frame(chrom=TSS_in[,1],
                      TSS=rowMeans(TSS_in[,c(2,3)]),
                      gene=rep("null",nrow(TSS_in))
                      )
  
  
  
  # head(regionsBed)
  # V1      V2      V3    V4 V5 V6
  # chr1 3062338 3063675  .  0  .
  # chr1 3343542 3344454  .  0  .

  # head(TSS_loc)
  # chrom     TSS   gene 
  # chr1  3671498   Xkr4   
  # chr1  4360303    Rp1  

  ###
  kbDistance<-kbDistance*1000

  selectedRegions<-integer() ##To track the positions of the maintained regions
  proximalRegions<-integer()
  
  for(nRegion in 1:nrow(regionsBed)){
    
    ##nRegion will also be used to select and track the mantained rows
    
    #to know the process status
    #we print each thousand multiple
    if(nRegion %% 1000 == 0){
      cat("nStudiedREgions: ",nRegion,"\n")
    }
    
    #######
    targetRegionStart<-regionsBed[nRegion,2]
    targetRegionEnd<-regionsBed[nRegion,3]
    targetRegionChr<-regionsBed[nRegion,1]
    
    
    ###Check if any TSS included in this region
    ## ----10kb-----startEnhancer----endEnhancer-------10KB---
    
    TSS_loc$isChr<-TSS_loc$chrom==targetRegionChr
    ##Start comprehends enhancer starts - Thresholdkb
    
    if(startend_mean == TRUE){
      TSS_loc$isAfterStart<-(TSS_loc$TSS >= (targetRegionStart - kbDistance))  
      TSS_loc$isBeforeEnd<-(TSS_loc$TSS <= (targetRegionEnd + kbDistance))      
    }
    

    
    logicalInfo<-TSS_loc[,c("isChr","isAfterStart","isBeforeEnd")]
    TSSmapped<-any(rowSums(logicalInfo)==3) ##if TRUE, at least a TSS located in the region enh+-10kb
    #print(paste0("TSSmapped=",TSSmapped,",nRegion=",nRegion))
    
    if(TSSmapped == FALSE){
      ##So, for the current studied region, no overlap with any TSS
      ##Remember this row, will be selected afterwards
      selectedRegions<-c(selectedRegions, nRegion)
    } else {
      proximalRegions<-c(proximalRegions, nRegion)
    }
    
    
  }
  
  ##Let's mantain the rows distal at windows
  resultMatrix<-regionsBed[selectedRegions,]
  
  ##Let's mantain the rows proximal at window
  proximalMatrix<-regionsBed[proximalRegions,]
  
  results<-list(distal=resultMatrix, proximal=proximalMatrix)

  return(results)
}


