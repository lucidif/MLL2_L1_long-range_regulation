# =========================
# R script to build anchors1, anchors2, and anchors3
# (merges your 3 original scripts, fixes minor bugs, adds safety checks)
# =========================

suppressWarnings({
  options(stringsAsFactors = FALSE, scipen = 999)
})

# -------------------------
# Parameters (edit as needed)
# -------------------------
outpath <- "outs/coolpup/500bp"             # base output folder
windowing <- TRUE                            # whether to generate windowed anchors1
windows <- c(1, 500, 1000)                   # windows (bp half-size) for anchors1 and anchors3

# Input paths
lines_file_path <- "./outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm.l1.bed" 
peaks_file_path <- "./outs/CHiP_postprocessing_line1_dist/distfilter500_DKO_K4me3_dcm.target.peaks.bed"
great_genes_path <- "./outs/great/Double_KO_vs_F_F/basal/filtered_by_distance_500bp/dcm/20241126-public-4.0.4-hpxvmj-mm10-all-gene.txt"
anno_path <- "./in/GREATv4/GREATv4.genes.mm10.tsv"  # GREAT annotation; must contain gene symbol in column V5

# Optional inputs for script 2 and 3
# Script 2 (DEGs to filter anchors1)
deg_path <- "outs/overlap/D0.down.tsv"      # header=TRUE, first column gene symbols

# Script 3 (proximal peaks list for make.anch)
proximal_peaks_path <- "outs/great/Double_KO_vs_F_F/basal/proximal/20241126-public-4.0.4-xCT3F1-mm10-all-gene.txt"

# -------------------------
# Helpers
# -------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_tsv_any <- function(path, header = FALSE) {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path))
  read.table(path, sep = "\t", header = header, check.names = FALSE, quote = "")
}

centered_window <- function(start, end, half_win) {
  mid <- round(start + (end - start)/2, 0)
  data.frame(start = mid - half_win, end = mid + half_win)
}

orient_bedpe <- function(df) {
  # df must have columns: chrA,startA,endA,chrB,startB,endB
  stopifnot(all(c("chrA","startA","endA","chrB","startB","endB") %in% colnames(df)))
  out <- df
  meanA <- round(df$startA + (df$endA - df$startA)/2, 0)
  meanB <- round(df$startB + (df$endB - df$startB)/2, 0)
  swap <- which(meanB < meanA)
  if (length(swap)) {
    out[swap, c("chrA","startA","endA","chrB","startB","endB")] <-
      df[swap, c("chrB","startB","endB","chrA","startA","endA")]
  }
  out
}

write_bedpe <- function(df, file) {
  ensure_dir(dirname(file))
  write.table(df, file = file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}

# -------------------------
# Step 1: Build anchors1 (from GREAT genes + peaks)
# -------------------------
build_anchors1 <- function(out_dir = outpath, windows = windows, windowing = windowing) {
  message("[anchors1] Reading inputs …")
  # Optional file (not used), kept for parity with original script
  if (file.exists(lines_file_path)) invisible(read_tsv_any(lines_file_path, header = FALSE))

  peaks.file <- read_tsv_any(peaks_file_path, header = FALSE)
  rownames(peaks.file) <- peaks.file$V4
  great.basal.dcm.genes <- read_tsv_any(great_genes_path, header = FALSE)
  anno <- read_tsv_any(anno_path, header = FALSE)

  message("[anchors1] Expanding GREAT peak lists per gene …")
  # original: V1 = gene, V2 = comma-separated lists containing strings with "peak" tokens
  tt <- lapply(seq_len(nrow(great.basal.dcm.genes)), function(i){
    gene <- great.basal.dcm.genes[i,1]
    s <- as.character(great.basal.dcm.genes$V2[i])
    toks <- unlist(strsplit(s, ","))
    toks <- unlist(strsplit(toks, " "))
    peaks <- toks[grepl("peak", toks)]
    if (length(peaks)) data.frame(peak = peaks, gene = gene, stringsAsFactors = FALSE) else NULL
  })
  ttfinal <- do.call(rbind, tt)
  if (is.null(ttfinal) || !nrow(ttfinal)) stop("[anchors1] No peak↔gene associations parsed from GREAT file.")

  message("[anchors1] Joining with gene annotation and peak coordinates …")
  # anno: expects gene symbol in V5
  gene.coords <- merge(ttfinal, anno, by.x = "gene", by.y = "V5")
  colnames(gene.coords)[2:6] <- c("peaks","ensembl","gene.chr","gene.tss","gene.strand")

  gene.peaks.coords <- merge(gene.coords, peaks.file, by.x = "peaks", by.y = "V4")
  colnames(gene.peaks.coords)[7:11] <- c("peaks.chr","peak.start","peak.end","peak.value","peak.strand")

  bedpe1_unsort <- data.frame(
    chrA   = gene.peaks.coords$peaks.chr,
    startA = gene.peaks.coords$peak.start,
    endA   = gene.peaks.coords$peak.end,
    chrB   = gene.peaks.coords$gene.chr,
    startB = gene.peaks.coords$gene.tss - 1,
    endB   = gene.peaks.coords$gene.tss + 1,
    peak   = gene.peaks.coords$peaks,
    gene   = gene.peaks.coords$gene,
    check.names = FALSE
  )

  bedpe1 <- orient_bedpe(bedpe1_unsort)

  if (isTRUE(windowing)) {
    for (tarwin in windows) {
      message(sprintf("[anchors1] Writing win%d_anchors1.bedpe", tarwin))
      meanA <- round(bedpe1$startA + (bedpe1$endA - bedpe1$startA)/2, 0)
      meanB <- round(bedpe1$startB + (bedpe1$endB - bedpe1$startB)/2, 0)
      windowed <- data.frame(
        chrA = bedpe1$chrA,
        startA = meanA - tarwin,
        endA = meanA + tarwin,
        chrB = bedpe1$chrB,              # FIX: use chrB (not chrA)
        startB = meanB - tarwin,
        endB = meanB + tarwin,
        peak = bedpe1$peak,
        gene = bedpe1$gene,
        check.names = FALSE
      )
      write_bedpe(windowed, file.path(out_dir, sprintf("win%d_anchors1.bedpe", tarwin)))
    }
  } else {
    message("[anchors1] Writing 500_anchors1.bedpe")
    write_bedpe(bedpe1, file.path(out_dir, "500_anchors1.bedpe"))
  }

  invisible(bedpe1)
}

# -------------------------
# Step 2: Build anchors2 (filter anchors1 by DEGs and re-window)
# -------------------------
build_anchors2 <- function(out_dir = outpath, window = 500, anchors1_file = NULL, deg_file = deg_path) {
  if (is.null(anchors1_file)) anchors1_file <- file.path(out_dir, sprintf("win%d_anchors1.bedpe", window))

  message("[anchors2] Reading anchors1 and DEG list …")
  anch1 <- read_tsv_any(anchors1_file, header = FALSE)
  colnames(anch1) <- c("chrA","startA","endA","chrB","startB","endB","peak","gene")

  d0.down.degs <- read_tsv_any(deg_file, header = TRUE)
  if (!ncol(d0.down.degs)) stop("[anchors2] DEG file seems empty or malformed.")
  gene_col <- colnames(d0.down.degs)[1]

  anch1.d0.down <- merge(anch1, d0.down.degs, by.x = "gene", by.y = gene_col)
  if (!nrow(anch1.d0.down)) stop("[anchors2] No overlap between anchors1 genes and DEG list.")

  anch2_unsort <- anch1.d0.down[, c("chrA","startA","endA","chrB","startB","endB")]
  colnames(anch2_unsort) <- c("chrA","startA","endA","chrB","startB","endB")

  anch2 <- orient_bedpe(anch2_unsort)

  message(sprintf("[anchors2] Re-windowing at ±%d bp and writing win%d_anchors2.bedpe", window, window))
  meansA <- round(anch2$startA + (anch2$endA - anch2$startA)/2, 0)
  meansB <- round(anch2$startB + (anch2$endB - anch2$startB)/2, 0)
  anch2$startA <- meansA - window
  anch2$endA   <- meansA + window
  anch2$startB <- meansB - window
  anch2$endB   <- meansB + window

  write_bedpe(anch2, file.path(out_dir, sprintf("win%d_anchors2.bedpe", window)))
  invisible(anch2)
}

# -------------------------
# Step 3: Build anchors3 (make.anch) — merge anchors1 with proximal peaks list
# -------------------------
make_anch <- function(outfile, anch1.path, peaks.path) {
  anch1 <- read_tsv_any(anch1.path, header = FALSE)
  colnames(anch1) <- c("chrA","startA","endA","chrB","startB","endB","peak","gene")

  proximal.peaks <- read_tsv_any(peaks.path, header = FALSE)
  # proximal.peaks: first column must be gene symbol (or the key used at anch1$gene)

  anch3_unsort <- merge(anch1, proximal.peaks, by.x = "gene", by.y = 1)
  # Drop merge key and the extra V2 (if present) like original intent
  # Keep only 6 bedpe columns
  anch3_unsort <- anch3_unsort[, c("chrA","startA","endA","chrB","startB","endB")]
  colnames(anch3_unsort) <- c("chrA","startA","endA","chrB","startB","endB")

  anch3 <- orient_bedpe(anch3_unsort)
  write_bedpe(anch3, outfile)
  invisible(anch3)
}

# -------------------------
# Main: run all steps
# -------------------------
run_all <- function() {
  ensure_dir(outpath)

  # 1) anchors1 for all windows
  build_anchors1(out_dir = outpath, windows = windows, windowing = windowing)

  # 2) anchors2 for window 500 only (as in your original script 2)
  if (file.exists(deg_path)) {
    try({ build_anchors2(out_dir = outpath, window = 500, anchors1_file = file.path(outpath, "win500_anchors1.bedpe"), deg_file = deg_path) }, silent = FALSE)
  } else {
    message("[anchors2] DEG file not found; skipping anchors2.")
  }

  # 3) anchors3 for windows 1, 500, 1000 (as in your original script 3)
  for (w in windows) {
    anch1_file <- file.path(outpath, sprintf("win%d_anchors1.bedpe", w))
    outfile <- file.path(outpath, sprintf("win%d_anchors3.bedpe", w))
    if (file.exists(anch1_file) && file.exists(proximal_peaks_path)) {
      message(sprintf("[anchors3] Building %s", basename(outfile)))
      make_anch(outfile = outfile, anch1.path = anch1_file, peaks.path = proximal_peaks_path)
    } else {
      message(sprintf("[anchors3] Skipping win%d: missing %s or %s", w, anch1_file, proximal_peaks_path))
    }
  }

  message("All done.")
}


