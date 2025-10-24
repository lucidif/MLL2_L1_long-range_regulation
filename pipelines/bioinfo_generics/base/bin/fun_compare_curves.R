compair.curve<-function(df,a,b,kmatrix,grps, test=c("Wilcox","KS")){

  # a_VS_b<-t.test( df[which(df$group==a),"distance"] ,  
  #                     df[which(df$group==b),"distance"], 
  #                     alternative = "two.sided", var.equal = FALSE)
  
  if(test[1]=="Wilcox"){
    a_VS_b<-wilcox.test(df$distance[df$group == a], df$distance[df$group == b])
  }

  if(test[1]=="KS"){
    a_VS_b<-ks.test(df$distance[df$group == a], df$distance[df$group == b])
  }
  
  kmatrix[a,b]<-a_VS_b$p.value
  kmatrix[b,a]<-a_VS_b$p.value
  
  return(kmatrix)
  
}


wilcox.matrix<-function(df_windowed , file=file){

                #         head(df_windowed)
                #   distances             group
                # 1  443993.5 dis.CpG.plus_DOWN
                # 2  147220.5 dis.CpG.plus_DOWN
                # 3   75409.5 dis.CpG.plus_DOWN
                # 4  102460.0 dis.CpG.plus_DOWN
                # 5    9588.0 dis.CpG.plus_DOWN
                # 6   40921.0 dis.CpG.plus_DOWN    

                #grps<-c("dis.CpG.plus_DOWN",   "dis.CpG.plus_UP", "dis.CpG.plus_NODIFF")

                grps<-unique(df_windowed$group)

                kmatrix<-matrix(nrow=length(grps), ncol=length(grps)) # nolint
                colnames(kmatrix)<-grps
                rownames(kmatrix)<-grps

                # dcp DOWN vs dcp UP # nolint

                for (i in 1:length(grps)){ # nolint

                    a<- grps[i]

                    for(j in 1:(length(grps)-1))
                    {
                        b <- grps[j+1]
                        kmatrix<-compair.curve(df=df_windowed, a, b, kmatrix, grps)
                    }


                }

              kmatrix<-cbind(samples=row.names(kmatrix),kmatrix)

                write.table(kmatrix, file=file, 
                col.names = TRUE, 
                row.names = FALSE, 
                quote=FALSE, 
                sep="\t")

                return(kmatrix)


        }

