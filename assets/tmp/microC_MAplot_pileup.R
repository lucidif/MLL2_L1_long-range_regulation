# ============================================================
# Read coolpup extracted matrices (avg_num_over_n) for two conditions,
# reshape them into long format (binX, binY, ave), merge by pixel,
# and build MA-plots + density plots.
#
# Notes:
# - Input TSVs are matrices exported from the .clpy (avg = num/n per pixel).
# - M = log2(KO/WT) is the log-fold-change per pixel.
# - A = 0.5*log2(KO*WT) is the log-intensity (geometric mean, in log space).
# - nologA = sqrt(KO*WT) is the geometric mean WITHOUT log (linear X axis).
# ============================================================

# ----------------------------
# Helper: read a TSV matrix and return a long dataframe
# ----------------------------
coolpop_MA_plot_in <- function(avema1path, avema2path){
  # Read condition 1 matrix (TSV -> matrix -> table -> long dataframe)
  # as.table() keeps row/col names (bin indices) as factors Var1/Var2
  ma1 <- as.table(as.matrix(read.table(avema1path, sep="\t")))
  ma1 <- as.data.frame(ma1)
  colnames(ma1) <- c("binX", "binY", "WT_ave")  # rename columns (cond1 assumed WT)

  # Read condition 2 matrix (same transformation)
  ma2 <- as.table(as.matrix(read.table(avema2path, sep="\t")))
  ma2 <- as.data.frame(ma2)
  colnames(ma2) <- c("binX", "binY", "KO_ave")  # rename columns (cond2 assumed KO)

  # Quick sanity peek (prints to console)
  head(ma1)
  head(ma2)

  # Create a unique pixel ID per (binX, binY) so we can merge the two conditions
  ma1 <- cbind(unid = paste0(ma1$binX, "_", ma1$binY), ma1)
  ma2 <- cbind(unid = paste0(ma2$binX, "_", ma2$binY), ma2)

  # Merge the two long tables by pixel ID
  # Output contains WT_ave and KO_ave for the same pixel
  ma_all <- merge(ma1, ma2, by="unid")

  return(ma_all)
}

# chromsize: data.frame con colonne V1 (chr) e V2 (len)
get_chr_len <- function(chr, chromsize) {
  chromsize$V2[match(chr, chromsize$V1)]
}

pick_viewpoint <- function(start, end, method=c("center","start","end")) {
  method <- match.arg(method)
  if (method == "center") return(as.integer(floor((start + end) / 2)))
  if (method == "start")  return(as.integer(start))
  as.integer(end)
}

# dato chr + coord (bp), ritorna bin index (0-based) e start/end del bin
coord_to_bin <- function(coord, binsize) {
  bin_id <- as.integer(floor(coord / binsize))  # 0-based
  bin_start <- bin_id * binsize
  bin_end   <- bin_start + binsize
  list(bin_id=bin_id, bin_start=bin_start, bin_end=bin_end)
}

# costruisce [start,end] esteso di window_bins prima/dopo del bin del viewpoint
make_window_locus <- function(chr, start, end, chromsize, binsize,
                              viewpoint=c("center","start","end"),
                              window_bins=10) {
  viewpoint <- match.arg(viewpoint)
  chr_len <- get_chr_len(chr, chromsize)
  if (is.na(chr_len)) stop("Chrom non trovato in chromsize: ", chr)

  vp <- pick_viewpoint(start, end, viewpoint)
  b <- coord_to_bin(vp, binsize)

  win_start <- b$bin_start - window_bins * binsize
  win_end   <- b$bin_end   + window_bins * binsize

   # IMPORTANTISSIMO: make the coordinate of the bin unique
   win_end <- win_end - 1

  # clamp ai limiti del cromosoma (assumendo coordinate 0-based come Juicer)
  win_start <- max(0, win_start)
  win_end   <- min(chr_len, win_end)

  list(
    chr = chr,
    vp = vp,
    bin_id = b$bin_id,
    locus_vec = c(chr, win_start, win_end) # per paste(..., collapse=":")
  )
}

#===============================
# main functions
#===============================

extract_contacts<-function(hicpath, bedpe.path,  
window_bins=11, 
viewpoint="center", 
bin.size=5000,
chromsize.file,
return_mode = c("list", "sum_df"),
verbose=FALSE
){

# Read inputs robustly
ach3 <- read.table(
  bedpe.path,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c("character","integer","integer","character","integer","integer")
)


# Pre-alloc outputs
  res_list <- vector("list", nrow(ach3))
  out_df <- data.frame(
    anchor = character(nrow(ach3)),
    bait   = character(nrow(ach3)),
    value  = rep(NA_real_, nrow(ach3)),
    stringsAsFactors = FALSE
  )


chromsize <- read.table(
  chromsize.file,
  header = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c("character", "integer")
)
#colnames(chromsize) <- c("chr", "len")

validr <- 0

for(i in 1:nrow(ach3)){

  anch <- ach3[i, c(1,2,3)]  # chr, start, end
  bait <- ach3[i, c(4,5,6)]

  # costruisci loci estesi basati sul bin del viewpoint
  anch_win <- make_window_locus(
    chr = as.character(anch[1]),
    start = as.integer(anch[2]),
    end   = as.integer(anch[3]),
    chromsize = chromsize,
    binsize = bin.size,
    viewpoint = viewpoint,
    window_bins = window_bins
  )

  bait_win <- make_window_locus(
    chr = as.character(bait[1]),
    start = as.integer(bait[2]),
    end   = as.integer(bait[3]),
    chromsize = chromsize,
    binsize = bin.size,
    viewpoint = viewpoint,
    window_bins = window_bins
  )


  anch_loc <- paste(anch_win$locus_vec, collapse=":")
  bait_loc <- paste(bait_win$locus_vec, collapse=":")

  imported_hic <- straw(
    norm    = "KR",
    fname   = hicpath,
    chr1loc = anch_loc,
    chr2loc = bait_loc,
    unit    = "BP",
    binsize = bin.size
  )

  # Fill dataframe output always (NA if no contacts)
    out_df$anchor[i] <- anch_loc
    out_df$bait[i]   <- bait_loc
    out_df$value[i]  <- if (nrow(imported_hic) > 0) sum(imported_hic$counts) else NA_real_


 # Fill list output (NULL if no contacts)
    if(nrow(imported_hic) > 0){
      validr <- validr + 1
      res_list[[i]] <- list(
        i = i,
        anch = anch,
        bait = bait,
        anch_locus = anch_win$locus_vec,
        bait_locus = bait_win$locus_vec,
        anch_bin_id = anch_win$bin_id,
        bait_bin_id = bait_win$bin_id,
        contacts = imported_hic
      )

      if (verbose) {
        cat("i=", i,
            " anch_bin=", anch_win$bin_id,
            " bait_bin=", bait_win$bin_id,
            " anch_locus=", anch_loc,
            " bait_locus=", bait_loc,
            "\n")
      }
    } else {
      res_list[[i]] <- NULL
    }
  }

  if (return_mode == "sum_df") {
    return(out_df)
  } else {
    return(res_list)
  }

  

}


# ----------------------------
# Load D0 and D4 matrices (WT vs KO)
# ----------------------------
d0_ave <- coolpop_MA_plot_in(
  avema1path = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/WT_day0.Dd_extracted/avg_num_over_n_row0.tsv",
  avema2path = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/D0/5kb/KO_day0.Dd_extracted/avg_num_over_n_row0.tsv"
)

d4_ave <- coolpop_MA_plot_in(
  avema1path = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/WT_day4.Dd_extracted/avg_num_over_n_row0.tsv",
  avema2path = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/D4/5kb/KO_day4.Dd_extracted/avg_num_over_n_row0.tsv"
)

# ----------------------------
# Compute MA quantities (D0)
# ----------------------------

# If you ever have zeros, you can add a small eps to avoid log(0).
# eps <- 1e-6

# A (log-intensity) = 0.5 * log2(KO*WT)  == log2( sqrt(KO*WT) )
# M (log-fold-change) = log2(KO/WT)
d0_ave$A <- 0.5 * log2(d0_ave$KO_ave * d0_ave$WT_ave)
d0_ave$M <- log2(d0_ave$KO_ave / d0_ave$WT_ave)

# Same as A but without log on X: geometric mean in linear scale
d0_ave$nologA <- sqrt(d0_ave$WT_ave * d0_ave$KO_ave)

# Summaries of per-pixel averages (useful to sanity-check ranges)
summary(d0_ave$KO_ave)
summary(d0_ave$WT_ave)

# ----------------------------
# MA-plot D0 (log X)
# ----------------------------
pdf("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/MAplot_D0_KO_vs_WT.pdf",
    width = 6, height = 5)

plot(d0_ave$A, d0_ave$M,
     pch = 16, cex = 0.4,
     xlab = "A = 0.5*log2(KO*WT)",
     ylab = "M = log2(KO/WT)")
abline(h = 0, lty = 2)  # M=0 line means KO==WT

dev.off()

# ----------------------------
# MA-plot D0 (linear X via geometric mean)
# ----------------------------
pdf("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/nologX_MAplot_D0_KO_vs_WT.pdf",
    width = 6, height = 5)

plot(d0_ave$nologA, d0_ave$M,
     pch = 16, cex = 0.4,
     xlab = "A = sqrt(WT * KO)",
     ylab = "M = log2(KO/WT)")
abline(h = 0, lty = 2)

dev.off()

# ----------------------------
# Extra sanity checks on sign/range of logFC
# ----------------------------
range(log2(d0_ave$KO_ave / d0_ave$WT_ave))
range(log2(d0_ave$WT_ave / d0_ave$KO_ave))

# Count how many pixels have KO > WT (should match sign expectation for M)
table(d0_ave$KO_ave > d0_ave$WT_ave)

# ----------------------------
# Compute MA quantities (D4)
# ----------------------------
d4_ave$A <- 0.5 * log2(d4_ave$KO_ave * d4_ave$WT_ave)
d4_ave$M <- log2(d4_ave$KO_ave / d4_ave$WT_ave)

d4_ave$nologA <- sqrt(d4_ave$WT_ave * d4_ave$KO_ave)

# ----------------------------
# MA-plot D4 (log X)
# ----------------------------
pdf("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/MAplot_D4_KO_vs_WT.pdf",
    width = 6, height = 5)

plot(d4_ave$A, d4_ave$M,
     pch = 16, cex = 0.4,
     xlab = "A = 0.5*log2(KO*WT)",
     ylab = "M = log2(KO/WT)")
abline(h = 0, lty = 2)

dev.off()

# ----------------------------
# MA-plot D4 (linear X via geometric mean)
# ----------------------------
pdf("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/nologX_MAplot_D4_KO_vs_WT.pdf",
    width = 6, height = 5)

plot(d4_ave$nologA, d4_ave$M,
     pch = 16, cex = 0.4,
     xlab = "A = sqrt(WT * KO)",
     ylab = "M = log2(KO/WT)")
abline(h = 0, lty = 2)

dev.off()

# ----------------------------
# Compare discretization / uniqueness of values across timepoints
# (this can hint at different smoothing/quantization in matrices)
# ----------------------------
length(unique(round(d0_ave$WT_ave, 4)))
length(unique(round(d4_ave$WT_ave, 4)))

# ============================================================
# Density plots: compare distributions of per-pixel averages
# ============================================================

pdf("outs/Lara_multiomic_analysis/outs/coolpup/500bp/MAplot/tmp/density_WT_KO_D0_D4.pdf",
    width = 7, height = 5)

# Kernel density estimates (overall distribution of pixel values)
d0_wt <- density(d0_ave$WT_ave, na.rm=TRUE)
d0_ko <- density(d0_ave$KO_ave, na.rm=TRUE)
d4_wt <- density(d4_ave$WT_ave, na.rm=TRUE)
d4_ko <- density(d4_ave$KO_ave, na.rm=TRUE)

plot(d0_wt, main="Density ave: WT/KO (D0 vs D4)", xlab="ave", lwd=2)
lines(d0_ko, lwd=2, lty=2)
lines(d4_wt, lwd=2, lty=3)
lines(d4_ko, lwd=2, lty=4)

legend("topright",
       legend=c("D0 WT","D0 KO","D4 WT","D4 KO"),
       lwd=2, lty=c(1,2,3,4), bty="n")

dev.off()

# ============================================================
# Highlight central 16x16 block in MA-plot (center pixels of pileup)
# ============================================================

# Define 16 central bins along X (your custom naming scheme)
center16x <- c("M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","A1","B1")

# Define 16 central bins along Y (V13..V28)
center16y <- paste0("V", 13:28)

# NOTE: tar_unid below creates only 16 paired IDs (diagonal-like),
# NOT the full 16x16 cartesian product. It's not used later.
tar_unid <- paste0(center16x, "_", center16y)

# Select pixels that fall inside the 16x16 central block.
# After merge(), binX/binY columns are duplicated as binX.x/binY.x and binX.y/binY.y.
# Use consistently the .x version (WT side) to define pixel coordinates.
is_center_d0 <- d0_ave$binX.x %in% center16x & d0_ave$binY.x %in% center16y
table(is_center_d0)  # expect TRUE = 256 if you indeed select a 16x16 block

is_center_d4 <- d4_ave$binX.x %in% center16x & d4_ave$binY.x %in% center16y
table(is_center_d4)

# Quick interactive plots (not saved) to check highlighting
plot(d0_ave$A, d0_ave$M, pch=16, cex=0.5, col="grey70")
points(d0_ave$A[is_center_d0], d0_ave$M[is_center_d0], pch=16, cex=0.6, col="red")
abline(h=0, lty=2)

plot(d4_ave$A, d4_ave$M, pch=16, cex=0.5, col="grey70")
points(d4_ave$A[is_center_d4], d4_ave$M[is_center_d4], pch=16, cex=0.6, col="red")
abline(h=0, lty=2)

# Attach center flag to d4_ave for inspection of selected rows
test_d4_ave <- cbind(d4_ave, is_center_d4 = is_center_d4)

# Print only the central-block pixels (256 rows) for inspection
test_d4_ave[which(test_d4_ave$is_center_d4 == TRUE),]




