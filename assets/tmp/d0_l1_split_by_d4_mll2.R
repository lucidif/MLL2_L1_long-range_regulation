
heatmap_tab2<-read.delim("outs/l1_d4_heatmap/sort2_l1_d4mll2_plotHeatmap.mat.tab", skip = 1, header = FALSE, sep = "\t")


trackbins<-ncol(heatmap_tab2[,7:ncol(heatmap_tab2)])/25

track2_colstart<-6+trackbins+1
track2_colend<-6+trackbins+trackbins

#median.trk2 <- apply(heatmap_tab2[, track2_colstart:track2_colend], 1, function(row) {
#  mean(as.numeric(row), na.rm = TRUE)
#})

delta<-round(nrow(heatmap_tab2)/3,0)

heatmap_tab2<-heatmap_tab2[order(median.trk2, decreasing = TRUE),]

high<-heatmap_tab2[1:delta,]
mid<-heatmap_tab2[(delta+1):(delta*2),]
low<-heatmap_tab2[(delta*2):nrow(heatmap_tab2),]

write.table(high[,1:6],
file="outs/l1_d4_heatmap/high_sort2_l1_d4mll2_plotHeatmap.bed",
sep="\t",
quote=FALSE,
row.names = FALSE,
col.names = FALSE
)

write.table(mid[,1:6],
file="outs/l1_d4_heatmap/mid_sort2_l1_d4mll2_plotHeatmap.bed",
sep="\t",
quote=FALSE,
row.names = FALSE,
col.names = FALSE
)

write.table(low[,1:6],
file="outs/l1_d4_heatmap/low_sort2_l1_d4mll2_plotHeatmap.bed",
sep="\t",
quote=FALSE,
row.names = FALSE,
col.names = FALSE
)

mean(as.numeric(mid[nrow(mid),7:ncol(mid)]))
mean(as.numeric(low[1,7:ncol(mid)]))


length(7:ncol(mid))

heatmap_tab3<-read.delim("outs/l1_d4_heatmap/div_l1_d4mll2_deeptools_matrix.gzip", skip = 1, header = FALSE, sep = "\t")




mean(as.numeric(mid[nrow(mid),track2_colstart:track2_colend]))
mean(as.numeric(low[1,track2_colstart:track2_colend]))
median(as.numeric(mid[nrow(mid),track2_colstart:track2_colend]))
median(as.numeric(low[1,track2_colstart:track2_colend]))


