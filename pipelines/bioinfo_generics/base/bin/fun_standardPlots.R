#fun index
#----1)plot.correlation
#----2)plot.pca

#====Requirements
#----a)ggplot2

#============

#============startFuns

plot.correlation<-function(x,y,xname="x",yname="y",title="CorrelationPlot",limit=FALSE,fixLimit=NULL, symmetric=FALSE,cor.method="pearson", labelPoints=NULL, xlim=NULL, ylim=NULL, debug=FALSE, norm.x=NULL,norm.y=NULL,color=NULL){

  #norm.x ; norm.y = c(method,steps)
  #Accepted methods exsamples:
  #c("quantile",10)
  #c("intervals",10)

  library(ggplot2)

  x.val<-x
  y.val<-y


  # x.val<-as.numeric(as.character(x))
  # y.val<-as.numeric(as.character(y))

  xna<-which(is.na(x.val))


  if(length(xna)>0){
    x.val<-x.val[-xna]
    y.val<-y.val[-xna]
  }

  yna<-which(is.na(y.val))

  if(length(yna)>0){
    x.val<-x.val[-yna]
    y.val<-y.val[-yna]
  }

  #normalize data

  int.quantileNorm<-function(val,normStrat){
    steps<-as.numeric(as.character(normStrat[2]))
    num<-as.numeric((as.character(val)))

    st<-1/steps
    bins<-seq(st,1,st)

    quant<-quantile(num,bins)

    #TODO controlla se più quantile hanno lo stesso valore (in questo caso sono troppi quantili)

    switch<-c()

    for(i in 1:length(quant)){
      quant[i]
      if(i == 1){
        switch[which(as.numeric(as.character(val))<=as.numeric(quant[i]))]<-gsub("%","",names(quant[i]))
      }else{
        switch[which(as.numeric(as.character(val))<=as.numeric(quant[i]) & as.numeric(as.character(val)) > as.numeric(quant[i-1]))]<-gsub("%","",names(quant[i]))
      }


    }

    return(switch)

  }

  if(!is.null(norm.x)){

    if(norm.x[1]=="quantile"){

      x.val<-as.factor(as.numeric(int.quantileNorm(val=x.val,normStrat=norm.x)))

    }

    if(norm.x[1]=="intervals"){

    }

  }

  if(!is.null(norm.y)){

    if(norm.y[1]=="quantile"){

      y.val<-as.factor(as.numeric(int.quantileNorm(val=y.val,normStrat=norm.y)))

    }

    if(norm.y[1]=="intervals"){

    }

  }



  corr<-cor(as.numeric(x.val),as.numeric(y.val),method=cor.method)

  corrframe<-data.frame(x.val,y.val)
  colnames(corrframe)<-c(xname,yname)

  if(debug==TRUE){
    print("tag1")
  }

  if(limit==TRUE){
    if(is.null(fixLimit)){ # se non viene selezionato un fixLimit prende uppertquartile

      #if(symmetric==TRUE){
      #  limit.val<-max(quantile(x.val)[4],quantile(y.val)[4])
      #  limit.xval<-limit.val
      #  limit.yval<-limit.val
      #}else{
        limit.xval<-max(quantile(as.numeric(x.val))[4])
        limit.yval<-max(quantile(as.numeric(y.val))[4])
      #}

    }else{
      limit.xval<-fixLimit[1]
      limit.yval<-fixLimit[2]
    }
  }else{
    #limit.val<-max(x.val,y.val)
    limit.xval<-max(as.numeric(x.val))
    limit.yval<-max(as.numeric(y.val))
  }

  if(debug==TRUE){
    print("tag2")
  }


  if(symmetric==TRUE){
    limit.val<-max(limit.xval,limit.yval)
    limit.xval<-limit.val
    limit.yval<-limit.val
  }

  if(debug==TRUE){
    print("tag3")
  }



  if(is.null(labelPoints)){

    if(debug==TRUE){
      print("tag4")
    }

    ggplot(corrframe, aes(x=x.val, y=y.val) ) + geom_point(aes(color=color))+ geom_smooth(method=glm,fullrange=FALSE,se=FALSE)+ xlab(label=xname)+ ylab(label=yname)+ ggtitle(paste0(title,"(corr=",round(corr,3),")"))+ theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(face = "bold"))+ theme(plot.title = element_text(lineheight = 3.2)) #+ xlim(0,limit.xval)+ ylim(0,limit.yval)

    #TODO risolvi il problema dei limit val x e y con i factors

  }else{

    if(debug==TRUE){
      print("tag5")
    }

    corrframe<-as.data.frame(cbind(labelPoints,corrframe))
    ggplot(corrframe, aes(x=x.val, y=y.val ,label=labelPoints) ) + geom_point(aes(color=color))+ geom_text(aes(label=labelPoints),hjust=0, vjust=0)+ geom_smooth(method=glm,fullrange=FALSE,se=FALSE)+ xlab(label=xname)+ ylab(label=yname)+ ggtitle(paste0(title,"(corr=",round(corr,3),")"))+ theme(plot.title = element_text(hjust = 0.5))+ theme(plot.title = element_text(face = "bold"))+ theme(plot.title = element_text(lineheight = 3.2))#+ xlim(0,limit.xval)+ ylim(0,limit.yval)
  }





}





plot.pca <- function(data, groups, gplabs = groups, pca_f, pcaloc = 'topright', legendpage = FALSE, monoLayout=TRUE, sampleNames=TRUE,legendSize=0.6,px=1,py=2){

  #pca_f --> folder and pca file name
  #data --> cpm raw counts


  #==================old one
  # names(gplabs) <- groups
  # cc <- as.factor(groups)
  # names(cc) <- gplabs
  # pca_cl <- prcomp(t(log2(data + 1)))
  # pdf(pca_f)
  #   # Add extra space to right of plot area; change clipping to figure
  #   if (monoLayout){par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, pin=c(4.4,4))}
  #   plot(pca_cl$x[, c(1,2)], pch = unclass(cc), col = unclass(cc))
  #   if (sampleNames){
  #     text(pca_cl$x[, c(1,2)], row.names(pca_cl$x), cex = 0.5, pos = 3)
  #   }else{
  #       text(pca_cl$x[, c(1,2)], cex = 0.5, pos = 3)
  #   }
  #   if (legendpage) {
  #     plot.new()
  #   }
  #   if (monoLayout){
  #   legend("topright", inset=c(-0.2,0) , legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA', cex=0.8)
  #   }else{
  #   legend(pcaloc, legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA')
  #   }
  #
  # dev.off()
  # # norm_data <- merge(ann, round(data, 2), by.x = 1, by.y = 0)
  # # return(norm_data)
  #======================

  names(gplabs) <- groups
  cc <- as.factor(groups)
  names(cc) <- gplabs
  pca_cl <- prcomp(t(log2(data + 1)))
  sumrepo<-as.data.frame(as.matrix(summary(pca_cl)$importance))
  sumrepo<-cbind(row.names(sumrepo),sumrepo)
  colnames(sumrepo)[1]<-" "
  write.table(sumrepo,file="pca_summary.txt",sep="\t",row.names=FALSE,quote=FALSE)
  pdf(pca_f)
  # Add extra space to right of plot area; change clipping to figure
  if (monoLayout){par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, pin=c(4.4,4))}
  plot(pca_cl$x[, c(px,py)], pch = unclass(cc), col = unclass(cc), xlab=paste0(colnames(pca_cl$x[, c(px,py)])[1],"(",as.character(round(sumrepo[2,which(colnames(sumrepo)==paste0("PC",px))],2)),")"),
       ylab=paste0(colnames(pca_cl$x[, c(px,py)])[2],"(",as.character(round(sumrepo[2,which(colnames(sumrepo)==paste0("PC",py))],2)),")")
  )
  if (sampleNames){
    text(pca_cl$x[, c(px,py)], row.names(pca_cl$x), cex = 0.5, pos = 3)
  }else{
    #text(pca_cl$x[, c(1,2)], cex = 0.5, pos = 3)
  }
  if (legendpage) {
    plot.new()
  }
  if (monoLayout){
    legend("topright", inset=c(-0.2,0) , legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA', cex=legendSize)
  }else{
    legend(pcaloc, legend = unique(names(cc)), pch = unclass(as.factor(unique(groups))), col = unique(unclass(cc)), fill = 'transparent', border = 'NA')
  }

  dev.off()


return(pca_cl)




}
