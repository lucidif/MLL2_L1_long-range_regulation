virtual4C_scatter<-function(bedpe, viewpoint, target=NULL, outfile, plot="bars"){

    #bedpe=""
    #viewpoint="chr10:51744068" #L1
    #target="chr10:51677810"
    #plot="bars"
    
    #viewpoint="chr10:51677810"
    #target="chr10:51744068"
    
    #outfile="outs/Lara_multiomic_analysis/outs/test_virtual4C/_hic_mm10_onefile2/virtual4C_vpgene_KO_day0_B.5kb.pdf"

    window="FULL" # add this as parameter, for now is fixed 

    intbin<-read.table(file=bedpe, sep="\t",header=FALSE)
    vp_chr=strsplit(viewpoint, ":")[[1]][1]
    vp_coord=as.numeric(strsplit(viewpoint, ":")[[1]][2])

    # target
    t_chr   <- strsplit(target, ":")[[1]][1]
    t_coord <- as.numeric(strsplit(target, ":")[[1]][2])
    
    idx <- which(
                (intbin$V1 == vp_chr & intbin$V2 <= vp_coord & intbin$V3 >= vp_coord) |
                (intbin$V4 == vp_chr & intbin$V5 <= vp_coord & intbin$V6 >= vp_coord)
            )

    keep <- rep(FALSE, nrow(intbin))
    keep[idx] <- TRUE

    vp_bin <- intbin
    vp_bin[!keep, 7] <- 0

    tt<-intbin[idx,]

    # intbin = bedpe letto, vp_bin = bedpe con V7 azzerato dove non interessa
    #vp_chr   <- strsplit(viewpoint, ":")[[1]][1]
    #vp_coord <- as.numeric(strsplit(viewpoint, ":")[[1]][2])

    # X = midpoint dell’anchor partner, Y = score (colonna 7)
    mid1 <- (vp_bin$V2 + vp_bin$V3) / 2
    mid2 <- (vp_bin$V5 + vp_bin$V6) / 2

    vp_on_1 <- (vp_bin$V1 == vp_chr & vp_bin$V2 <= vp_coord & vp_bin$V3 >= vp_coord)
    x_pos   <- ifelse(vp_on_1, mid2, mid1)

    ok <- vp_bin$V7 > 0 & !is.na(x_pos)

    # plot(x_pos[ok], vp_bin$V7[ok],
    #     pch=16, cex=0.5,
    #     xlab="posizione (midpoint anchor partner)", ylab="score (V7)")

    ord <- order(x_pos[ok])

ord_all <- order(x_pos)
x <- x_pos[ord_all]
y <- vp_bin$V7[ord_all]

# stima una larghezza barra (mediana dei delta tra punti)
w <- stats::median(diff(x), na.rm=TRUE)


    #dist <- x_pos - vp_coord
    
    


pdf(outfile, width = 7, height = 5)
    #==============================================
    # vp_chr, vp_coord già definiti
    # vp_on_1 <- (vp_bin$V1 == vp_chr & vp_bin$V2 <= vp_coord & vp_bin$V3 >= vp_coord)

    # # partner intervals (start/end) e score
    # p_chr   <- ifelse(vp_on_1, vp_bin$V4, vp_bin$V1)
    # p_start <- ifelse(vp_on_1, vp_bin$V5, vp_bin$V2)
    # p_end   <- ifelse(vp_on_1, vp_bin$V6, vp_bin$V3)
    # score   <- vp_bin$V7

    # # (opzionale ma consigliato) tieni solo partner sullo stesso chr del viewpoint
    # sel <- (p_chr == vp_chr) & (score > 0)

    # # ordina per start
    # ord <- order(p_start[sel])

    # x1 <- p_start[sel][ord]
    # x2 <- p_end[sel][ord]
    # y  <- score[sel][ord]

    # plot(c(min(x1), max(x2)), c(0, max(y)),
    #     type="n",
    #     xlab="coordinates (bins)", ylab="score (V7)")

    # rect(x1, 0, x2, y)

    # abline(v = t_coord, col="red", lty=2)
    # abline(v = vp_coord, col="blue", lty=2)
xmin <- min(x_pos[ok], vp_coord, t_coord, na.rm=TRUE)
xmax <- max(x_pos[ok], vp_coord, t_coord, na.rm=TRUE)

plot(x_pos[ok][ord], vp_bin$V7[ok][ord],
     type="b", pch=16, cex=0.4,
     xlab="positionn", ylab="score (V7)",
     xlim=c(xmin, xmax))

abline(v=t_coord, col="red", lty=2)
abline(v=vp_coord, col="blue", lty=2)

dev.off()



#===============================================
    

    

    # ord <- order(dist[ok])
    # plot(dist[ok][ord], vp_bin$V7[ok][ord],
    #     type="l",
    #     xlab="distanza dal viewpoint (bp)", ylab="score (V7)")
    # abline(v=0, lty=2)    

}


