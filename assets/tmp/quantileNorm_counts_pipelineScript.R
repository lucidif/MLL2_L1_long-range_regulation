#=====inside docker example
# sudo docker run --rm -it \
#   --user "$(id -u)":"$(id -g)" \
#   -v  /media/lucio/easystore:/media/lucio/easystore \
#   -v "$PWD":"$PWD" \
#   -w "$PWD" \
#   lucidif/quantilenormalization:edger4.8.1 \
#   R


##################################################################
## R SCRIPT for quantile normalizing counts from FeatureCounts
##################################################################
library(preprocessCore)
library(data.table)
library(edgeR)

#!/usr/bin/env Rscript
inputArguments = commandArgs(trailingOnly=TRUE) #Capturing arguments passed after the script
fileName <- inputArguments[1] #bins.countmatrix
print(fileName)

currentDirectory<- inputArguments[2] ##Directory from where the bash script has been executed
print(currentDirectory)

fullPathInputCounts<-paste0(currentDirectory,"/",fileName)
print(fullPathInputCounts)
## Load count matrix, removing all columns that we are not interested in:
raw.counts <- fread(fullPathInputCounts, data.table=FALSE, header=TRUE, skip=c(1))

#Sorting by Chr & Start in case required downstream
raw.counts<-raw.counts[order(raw.counts$Chr, raw.counts$Start),]

metadataPeaks<-raw.counts
rownames(metadataPeaks)<-rownames(raw.counts)<-raw.counts$Geneid
raw.counts <- raw.counts[,7:ncol(raw.counts), drop = FALSE] ##7 first columns are always metadata
##drop = FALSE, added to preserve dimensionality when one sample analysed (applied to get just  RPKMs) https://adv-r.hadley.nz/subsetting.html#simplify-preserve
##drop = FALSE added 15 april 2024

#Improving colnames (quitar el folder path, mantener solo el baseName, si genero outputName usando el path  va a ser problematico)
colnames(raw.counts)
colnames(raw.counts)<-basename(colnames(raw.counts))
colnames(raw.counts)

# class(raw.counts) ##DATA.FRAME, i have checked

#######################################################################################################################################
# Computing RPKMs
# First, computing CPMs (counts per million). Then, multiplying by 10.so... getting RPKMs
# Bec, we are working with 100bp bins, which in kb is 0.1, if we divide by 0.1 same as x10
# And as stated in bamCoverage: RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb))
# RPKM calculation will be done regardless of whether quantile norm applied or not downstream
#######################################################################################################################################
##I have checked locally that CPM function only divides reads by total read sum (in millions)
##"By default, the normalized library sizes are used in the computation for DGEList objects but simple column sums for matrices."
##But simple column sums for matrices (as it is the case since not using DGElist object)
##Tested locally in TestingCPM.R (/home/victor/Dropbox/Cantabria/PhD_Project/ScriptsPhd/ScriptsParaUsoLocal/Colaboraciones/maria/ATAC_Seq/DiffAccessibility_Analyses/NormalizationBWs_strategy_ATpoint_TMM/TestingCPM.R)
##alternativa via Sweep disponible tambien. ##Using cpm function in case it is more efficient in terms of performance (time)
rpkm.counts<-as.data.frame(cpm(raw.counts,log = FALSE)*10)##*10 scaling factor to make number larger, 1-more user friendly. 2-to avoid that 0.0 gets rounded by 0 in any step

##############################
## QUANTILE NORMALIZING DATA #
##############################
##Passing the rpkm.counts to quantile normalization (if only one sample analysed the result is the same)
quantileNorm.counts<-as.data.frame(normalize.quantiles(as.matrix(rpkm.counts)))

colnames(quantileNorm.counts)<-colnames(raw.counts)
rownames(quantileNorm.counts)<-rownames(raw.counts)

if(ncol(raw.counts)>1){
  #Consideration added 15 april 2024
  ##>1 Sample analysed (for only one no change, and also lapply function does not handle properly)

  print("Pre quantile norm summary")
  preNorm_summary<-lapply(raw.counts,summary)
  print(preNorm_summary)
  
  print("After quantile norm summary")
  postNorm_summary<-lapply(quantileNorm.counts,summary)  
  print(postNorm_summary)
  
}

##########################################
## ASSEMBLING AND WRITING BEDGRAPHS
##########################################

#If quantile norm performed working with qnorm df otherwise select the rpkm dataframe
masterNormalizedDataFrame<-quantileNorm.counts

for(targetFile in colnames(raw.counts)){
  print(targetFile)
  
  targetNormDF<-data.frame("chr"=metadataPeaks$Chr,
                           "start"=metadataPeaks$Start,
                           "end"=metadataPeaks$End,
                           "counts"=masterNormalizedDataFrame[[targetFile]])
  
  ##Remove file extension, i.e. .bam rename file ending
  noExtensionTargetFile<-unlist(strsplit(x=targetFile, split=".bam",fixed=TRUE))[1]##The script deals with bam files, so the ending has to be a .bam extension otherwise ERROR
  
  outputBedgraphName<-paste0(currentDirectory,"/",noExtensionTargetFile,".bdg")
  
  write.table(x=targetNormDF,
              file = outputBedgraphName,
              sep="\t",
              quote = F,
              row.names = F,
              col.names = F)
}


