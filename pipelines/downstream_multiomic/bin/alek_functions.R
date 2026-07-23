# =============================================================================
# SCRIPT INDEX / OVERVIEW
# =============================================================================
# This script is organized into functional blocks for typical Hi-C / Micro-C
# analysis workflows. Below is a map of what is inside and what each block does.
#
# 1) Preprocessing: genomic ranges & binning
#    - getGR(chr, start, end, offset)
#        Create GRanges intervals (optionally extended by 'offset') to represent
#        anchors/regions and enable overlap-based operations.
#    - binAnno(chrSizes, chroms, binSize) # nolint: commented_code_linter.
#        Build a genome bin annotation table (chr/start/end) that tiles each
#        chromosome with fixed-width bins (used as the coordinate system for matrices).
#
# 2) I/O + data structures: dumped interactions -> sparse matrices
#    - transformToSparseMatrix(df, row.annotation, col.annotation)
#        Convert a 3-column "triplet" interaction table (rowBin, colBin, value)
#        into a sparse Matrix aligned to the provided bin annotation.
#        (Includes +1 correction to match 0-based dumped coordinates to 1-based R bins.)
#    - read.hic_files(directory, prefix, suffix, binAnnotation, chroms)
#        Read per-chromosome dumped interaction files from disk and return a list
#        of sparse intra-chromosome matrices (one per chromosome).
#
# 3) Genome-wide matrix assembly (block -> whole genome)
#    - getGWmatrix(list.of.sm, ga)
#        Assemble a genome-wide sparse matrix from named blocks representing
#        chromosome pairs (e.g. "chr1 chr1", "chr1 chr2") by shifting block indices
#        into genome-wide coordinates defined by 'ga'.
#    - getGWmatrix_dedicated(list.of.sm, ga, which.object)
#        Same goal as getGWmatrix(), but for objects returned by IPF() per chromosome:
#        select either raw ("LFM") or balanced ("balanced") matrix and assemble genome-wide.
#    - getGWmatrix_intra(list.of.sm, ga, which.object)
#        Intra-chromosome version of the above (currently identical logic, kept separate
#        for clarity / future customization).
#
# 4) Hi-C normalization (matrix balancing / ICE-like)
#    - loss.function(mat)
#        Quantify how close a matrix is to being balanced (row/col sums ~ 1).
#    - Create_complete_matrix(x)
#        Force symmetry by rebuilding a symmetric sparse matrix from stored entries.
#    - IPF_alg(lfm, numberOfIterations)
#        Iterative Proportional Fitting: iteratively rescales rows/cols to balance the matrix.
#    - normalize_LFM_iteratively(x, numberOfIterations)
#        Wrapper that symmetrizes input, runs IPF, collapses scaling factors into biases,
#        and applies final balancing.
#    - IPF(x, numberOfIterations)
#        High-level interface returning a 3-part object: raw (LFM), balanced, and biases.
#
# 5) Loop-centric utilities (APA-style extraction and loop scoring)
#    - getPixelsAround(ROW, COL, SIZE, ID)
#        Generate a square grid of pixel coordinates around a focal pixel (loop center).
#    - GetCentroidSignal(loops, hic, size, CHROMS, which_file)
#        Extract a (2*size+1)^2 pixel window around each loop from per-chromosome matrices
#        (raw or balanced), returning pixel-level signal.
#    - ProcessLoops(lo)
#        Reshape the pixel-level output into one small square matrix per loop (for APA plots).
#    - GetCentroidSignalGW(loops, hic, size, CHROMS, ga)
#        Same as GetCentroidSignal(), but reads from a pre-assembled genome-wide matrix.
#    - compute_putative_loop_signal(ctcf_ctcf_loop_signal, central_sq)
#        Compute a central enrichment metric: mean core signal vs sampled background windows.
#    - getLoopSignal(ll, D)
#        Apply compute_putative_loop_signal() across many loops and return log2FC values.
#
# 6) Regions / visualization helpers
#    - getMatrixForARegion(lowertri, uppertri, CHROM, START, END, ga, upperLimit)
#        Create a “hybrid” matrix for plotting comparisons (one condition in lower triangle,
#        the other in upper triangle), optionally capping extreme values for visualization.
#
# 7) BED / narrowPeak parsing + motif annotation
#    - readBed_filterChroms(filePath, chroms, scoreCol)
#        Read BED-like file -> GRanges, add "chr" prefix, filter chromosomes, attach score.
#    - readBed_filterChromsStraded(filePath, chroms, scoreCol)
#        As above but also keep strand (from column 6).
#    - readNarrowPeak2getSummit(filePath, chroms, scoreCol)
#        Read narrowPeak -> 1bp summit GRanges (peakStart + summitOffset).
#    - appendMotifInformation(peakFile, motifFile, summitFile)
#        For each peak, count motifs near summit (±25bp) and store best motif score/strand/pos.
#
# 8) Feature signal in bins (ChIP/ATAC/RNA windows)
#    - getSignalInBins(bins, modification, featureColumn)
#        Overlap features onto bins and assign a per-bin signal (max over overlaps),
#        then reshape into a matrix based on your "bin groups" encoding.
#    - getSignalInBins400(bins, modification, featureColumn)
#        Same idea but fixed-width windows of 401 bins and using mean over overlaps.
#    - reorderMatrixBasedOnStrand(m, peaks)
#        Strand-normalize matrices around motifs by reversing rows for '-' strand peaks.
#    - GetTPM(m, cols2take, gn)
#        Compute TPM from raw counts and gene lengths.
#
# 9) Insulation score
#    - getSUMMEDsignal4bins(B, D, A, M)
#        Sum Hi-C signal in square windows around offset positions; core helper.
#    - InsulationScore(bin, mat, GAGR, distance, Area, HOWFAR)
#        Compute left/middle/right insulation components from balanced per-chrom matrices.
#    - InsulationScorePublishedKR(bin, mat, GAGR, distance, Area, HOWFAR)
#        Same but assumes matrices are directly stored (already KR-normalized).
#    - processIS(IS, GAGR)
#        Convert left/middle/right components into a single insulation score per locus.
#    - cleanIS(x)
#        Replace NA/Inf with 0 for robust downstream plotting/statistics.
#
# Notes / assumptions:
# - Coordinate conventions: many dumped Hi-C tables are 0-based; bins in R are 1-based.
#   This is why you often see "+1" when matching dumped coordinates to bin starts.
# - Ordering: genome-wide assembly assumes 'ga' is ordered by chromosome then start.
# - Duplicates: if the same function is defined twice later in the script, the last
#   definition overwrites earlier ones (consider removing duplicates to avoid confusion).
# =============================================================================


# ==========================
# Preprocessing functions
# ==========================
# These helpers are used to:
#  - represent genomic intervals as GRanges
#  - create fixed-size genome bins (bin annotation)
#  - convert sparse "triplet" interaction tables into sparse matrices
#  - read per-chromosome matrices dumped from .hic (or similar) into R objects
#  - assemble genome-wide matrices from per-chromosome sparse blocks
#
# NOTE: Most functions assume:
#  - 1-based coordinates for bins/annotations (typical in R/GenomicRanges)
#  - dumped interaction files store bin starts in 0-based (common for Juicer dumps),
#    hence the "+1" adjustments in transformToSparseMatrix().
#

#' getGR
#'
#' Build a GRanges object from vectors of chromosome, start, end.
#' Optionally expands each interval by an "offset" on both sides.
#'
#' @param chr    character vector of chromosome names (e.g. "chr1")
#' @param start  integer/numeric vector of start positions (same length as chr)
#' @param end    integer/numeric vector of end positions (same length as chr)
#' @param offset integer scalar; bases to extend left and right (can be 0)
#'
#' @return A GRanges with:
#'   - seqnames = chr
#'   - ranges   = [start-offset, end+offset]
#'   - strand   = "*" for all entries
#'   - metadata column "loop" as 1..N (useful as an ID)
#'
#' @usage
#'   gr <- getGR(chr = loops$chr1, start = loops$x1, end = loops$x2, offset = 500)
getGR <- function(chr, start, end, offset) {

  # Create GRanges representing the genomic intervals.
  # Rle() is used for memory efficiency when values repeat a lot.
  tp <- GRanges(
    seqnames = Rle(chr),
    ranges   = IRanges(start - offset, end = end + offset),
    strand   = Rle(rep("*", length(start))),
    loop     = seq_len(length(start))
  )

  # Force UCSC chromosome naming style (e.g., "chr1" not "1"),
  # so downstream overlap calls are consistent.
  seqlevelsStyle(tp) <- "ucsc"

  return(tp)
}


#' binAnno
#'
#' Create a data.frame of fixed-width bins covering a set of chromosomes.
#' Bins are 1-based and inclusive (start/end), matching typical GRanges usage.
#'
#' @param chrSizes  object that supports seqlengths(), e.g. a Seqinfo or BSgenome
#' @param chroms    character vector of chromosomes to bin (e.g. c("chr1","chr2"))
#' @param binSize   integer; bin size in bp (e.g. 10000 for 10 kb)
#'
#' @return data.frame with columns:
#'   - chr   chromosome
#'   - start bin start (1-based)
#'   - end   bin end (1-based)
#'
#' @details
#' If chromosome length is not divisible by binSize, the function pads the
#' effective length so that the last bin fully covers the tail.
#'
#' @usage
#'   ga <- binAnno(chrSizes = seqinfo_object, chroms = paste0("chr",1:19), binSize = 10000)
binAnno <- function(chrSizes, chroms, binSize) {

  do.call("rbind", lapply(chroms, function(chr) {

    # Extract chromosome length from a Seqinfo/BSgenome-like object
    chr.length <- seqlengths(chrSizes[chr])

    # If length is not multiple of binSize, pad it so bins tile evenly
    dongling <- chr.length %% binSize
    if (dongling > 0) {
      # + (binSize - dongling) makes it divisible; +1 avoids edge issues
      # when later generating sequences with "to" endpoints.
      chr.length <- chr.length + (binSize - dongling) + 1
    }

    # Produce consecutive bins:
    # start: 1, 1+binSize, ...
    # end:   binSize, 2*binSize, ...
    data.frame(
      chr   = chr,
      start = seq(1, chr.length - binSize, by = binSize),
      end   = seq(binSize, chr.length,       by = binSize),
      stringsAsFactors = FALSE
    )
  }))
}


#' transformToSparseMatrix
#'
#' Convert a dumped interaction table (triplet format) into a sparse Matrix.
#'
#' @param df             data.frame with at least 3 columns:
#'                        - df[,1] = row-bin coordinate (often 0-based start)
#'                        - df[,2] = col-bin coordinate (often 0-based start)
#'                        - df[,3] = contact count / weight
#' @param row.annotation data.frame for row bins (typically binAnnotation filtered to one chr)
#' @param col.annotation data.frame for col bins (typically same as row.annotation for intra-chr)
#'
#' @return A sparse Matrix with dimensions:
#'   nrow(row.annotation) x nrow(col.annotation)
#'
#' @details
#' The "+1" is intentional:
#'  - Many Hi-C dump formats store bin starts as 0-based.
#'  - Your binAnno uses 1-based starts (1, 1+binSize, ...).
#' So df[,1]+1 is matched against annotation start positions (annotation[,2] == start).
#'
#' @usage
#'   sm <- transformToSparseMatrix(df, row.annotation = ga_chr1, col.annotation = ga_chr1)
transformToSparseMatrix <- function(df, row.annotation, col.annotation) {

  # Pre-allocate an all-zero sparse matrix.
  # NOTE: Matrix(0, ...) creates a dense Matrix class; for large objects you may
  # prefer sparseMatrix(i=..., j=..., x=..., dims=...) directly.
  sm <- Matrix(0, nrow = nrow(row.annotation), ncol = nrow(col.annotation))

  # Map dumped coordinates to bin indices in the annotation:
  # annotation[,2] is assumed to be "start" (see binAnno output).
  tp <- cbind(
    row = match(as.numeric(df[, 1] + 1), row.annotation[, 2]),
    col = match(as.numeric(df[, 2] + 1), col.annotation[, 2])
  )

  # Fill matrix at matched (row, col) positions with the contact value.
  sm[tp] <- df[, 3]

  return(sm)
}

transformToSparseMatrix_bpStarts <- function(df, row.annotation, col.annotation) {
  # df must have columns: x, y, counts where x/y are bin starts in bp (as returned by strawr)
  # row.annotation/col.annotation must have a 'start' column (binAnno output)

  # Standardize column names if needed
  if (all(c("x","y","counts") %in% colnames(df))) {
    x <- as.numeric(df$x)
    y <- as.numeric(df$y)
    v <- as.numeric(df$counts)
  } else {
    x <- as.numeric(df[,1])
    y <- as.numeric(df[,2])
    v <- as.numeric(df[,3])
  }

  # Match directly on bin starts (no +1)
  i <- match(x, row.annotation$start)
  j <- match(y, col.annotation$start)

  # Drop rows that failed to match (safety)
  ok <- !(is.na(i) | is.na(j) | is.na(v))
  i <- i[ok]; j <- j[ok]; v <- v[ok]

  # Build sparse matrix directly (more efficient than filling a zero matrix)
  sparseMatrix(
    i = i,
    j = j,
    x = v,
    dims = c(nrow(row.annotation), nrow(col.annotation))
  )
}



#' read.hic_files
#'
#' Read per-chromosome interaction matrices that were dumped to text files
#' (e.g., from .hic using juicer tools "dump" or a similar exporter).
#'
#' @param directory      path to dumped files (string, should end with "/" if needed)
#' @param prefix         filename prefix before chromosome (string)
#' @param suffix         filename suffix after chromosome (string)
#' @param binAnnotation  data.frame created by binAnno() for all chromosomes
#' @param chroms         character vector of chromosomes to load
#'
#' @return A list of sparse matrices, one per chromosome, in the same order as chroms.
#'         Each element is an intra-chromosome sparse contact matrix.
#'
#' @usage
#'   mats <- read.hic_files(
#'     directory = dumped_dir,
#'     prefix    = "ES1_30_",
#'     suffix    = ".txt",
#'     binAnnotation = ga,
#'     chroms    = paste0("chr", 1:19)
#'   )
read.hic_files <- function(directory, prefix, suffix, binAnnotation, chroms) {

  # Make a named list of chromosomes so lapply keeps names
  chrlist <- as.list(chroms)
  names(chrlist) <- chroms

  lapply(as.list(chrlist), function(chromosome) {

    # Build path to the dumped file for this chromosome
    file_path <- paste0(directory, prefix, chromosome, suffix)

    # Read dump: expected 3 columns (rowStart, colStart, value)
    df <- read.delim(file_path, as.is = TRUE, header = FALSE)

    # Remove rows with missing interaction values (safety)
    df <- df[!is.na(df$V3), ]

    # Ensure numeric type for indexing
    df[, 1] <- as.numeric(df[, 1])

    # Subset bin annotation to this chromosome
    rowAnnotation <- binAnnotation[binAnnotation$chr == chromosome, ]
    colAnnotation <- binAnnotation[binAnnotation$chr == chromosome, ]

    print(chromosome)

    # Convert triplet table to sparse matrix
    transformToSparseMatrix(
      df,
      row.annotation = rowAnnotation,
      col.annotation = colAnnotation
    )
  })
}


# ==========================
# Genome-wide assembly helpers
# ==========================

#' getGWmatrix
#'
#' Assemble a genome-wide sparse matrix from a list of sparse blocks.
#'
#' @param list.of.sm named list of sparse matrices (blocks).
#'        Names are expected to encode chromosome pairs like "chr1 chr1", "chr1 chr2", etc.
#' @param ga data.frame genome bin annotation (from binAnno), concatenated across chromosomes
#'
#' @return A sparse Matrix of size nrow(ga) x nrow(ga) containing all blocks.
#'
#' @details
#' This function offsets per-chromosome block indices into the global coordinate
#' system defined by 'ga' (which stacks chromosomes one after another).
#'
#' IMPORTANT ASSUMPTION:
#'  - 'ga' is ordered by chromosome then by start coordinate.
#'  - list.of.sm names use the exact chromosome strings present in ga$chr.
getGWmatrix <- function(list.of.sm, ga) {

  # Pre-allocate an empty genome-wide sparse matrix
  nummSM <- Matrix(0, nrow = nrow(ga), ncol = nrow(ga))

  # Convert each block into triplets (i,j,x) and shift indices by chromosome offsets
  ids <- do.call("rbind", lapply(as.list(names(list.of.sm)), function(chrom_comb) {

    print(paste0("including ", chrom_comb))

    chrom <- unlist(strsplit(chrom_comb, " "))

    # Compute the global offsets for the start of each chromosome in 'ga'
    # start.row and start.col are 0-based offsets that will be added to the block i/j.
    if (chrom[1] == "chr1" & chrom[2] == "chr1") {
      start.row <- 0
      start.col <- 0
    }
    if (chrom[1] == "chr1" & chrom[2] != "chr1") {
      start.row <- 0
      start.col <- min(which(ga$chr == chrom[2])) - 1
    }
    if (chrom[1] != "chr1") {
      start.row <- min(which(ga$chr == chrom[1])) - 1
      start.col <- min(which(ga$chr == chrom[2])) - 1
    }

    # Extract non-zero entries as triplets from sparse matrix
    tp <- as.data.frame(summary(list.of.sm[[chrom_comb]]), stringsAsFactors = FALSE)

    # Shift local block indices into genome-wide indices
    tp$i <- tp$i + start.row
    tp$j <- tp$j + start.col

    return(tp)
  }))

  print("assembling the matrix, depending on the resolution this might take a while")

  # Fill genome-wide matrix with block entries
  nummSM[cbind(ids[, 1], ids[, 2])] <- ids[, 3]

  return(nummSM)
}


#' getGWmatrix_dedicated
#'
#' Genome-wide assembly for data structures returned by IPF().
#' Each chromosome element is expected to be a list with:
#'   [[1]] raw LFM, [[2]] balanced matrix, [[3]] biases
#'
#' @param list.of.sm named list where each element is itself a list (IPF output)
#' @param ga bin annotation (binAnno output for the whole genome)
#' @param which.object character, one of c("LFM", "balanced")
#'
#' @return genome-wide sparse Matrix for the chosen object.
getGWmatrix_dedicated <- function(list.of.sm, ga, which.object) {

  nummSM <- Matrix(0, nrow = nrow(ga), ncol = nrow(ga))

  # Select which object from IPF output:
  # [[1]] = LFM (raw), [[2]] = balanced
  if (which.object == "LFM") i <- 1
  if (which.object == "balanced") i <- 2

  ids <- do.call("rbind", lapply(as.list(names(list.of.sm)), function(chrom) {

    print(paste0("including ", chrom))

    # Offsets for this chromosome in the genome-wide annotation
    if (chrom == "chr1") {
      start.row <- 0; start.col <- 0
    }
    if (chrom != "chr1") {
      start.row <- min(which(ga$chr == chrom)) - 1
      start.col <- min(which(ga$chr == chrom)) - 1
    }

    # Extract non-zero triplets from the selected matrix
    tp <- as.data.frame(summary(list.of.sm[[chrom]][[i]]), stringsAsFactors = FALSE)

    # Shift indices into genome-wide coordinates
    tp$i <- tp$i + start.row
    tp$j <- tp$j + start.col

    return(tp)
  }))

  print("assembling the matrix, depending on the resolution this might take a while")

  nummSM[cbind(ids[, 1], ids[, 2])] <- ids[, 3]
  return(nummSM)
}

ICE_normalize <- function(M,
                          n_iter = 50,
                          tol = 1e-6,
                          exclude_diagonals = 0L,   # metti 1 o 2 se vuoi togliere diag/diag+1
                          verbose = TRUE) {
  stopifnot(inherits(M, "Matrix"), nrow(M) == ncol(M))

  # forza simmetria (ICE assume tipicamente matrice simmetrica)
  # se hai già Create_complete_matrix() nel tuo script, usa quello
  if (!isSymmetric(M)) {
    M <- (M + t(M)) / 2
  }

  A <- M

  # opzionale: azzera diagonali vicine (spesso si fa in Hi-C)
  if (exclude_diagonals > 0L) {
    n <- nrow(A)
    for (k in 0:exclude_diagonals) {
      idx <- seq_len(n - k)
      A[cbind(idx, idx + k)] <- 0
      A[cbind(idx + k, idx)] <- 0
    }
  }

  rs <- Matrix::rowSums(A)
  valid <- rs > 0
  if (sum(valid) < 2) stop("Troppi pochi bin validi (rowSums==0).")

  # lavora solo sul sottoinsieme valido
  A2 <- A[valid, valid, drop = FALSE]
  m <- nrow(A2)

  # bias (moltiplicativo)
  b <- rep(1, m)

  for (it in seq_len(n_iter)) {
    s <- as.numeric(A2 %*% b)
    s[s == 0] <- 1
    b_new <- b / s

    # rinormalizza scala per evitare drift numerico
    b_new <- b_new / exp(mean(log(b_new)))

    delta <- max(abs(b_new - b))
    b <- b_new

    if (verbose) message("ICE iter ", it, "  delta=", signif(delta, 3))
    if (delta < tol) break
  }

  # matrice bilanciata nel sottoinsieme: diag(b) A2 diag(b)
  A2_bal <- (b * A2) * rep(b, each = m)

  # rimappa a dimensione piena
  A_bal <- Matrix::Matrix(0, nrow = nrow(A), ncol = ncol(A), sparse = TRUE)
  A_bal[valid, valid] <- A2_bal

  # bias full-length (NA sui bin invalidi)
  b_full <- rep(NA_real_, nrow(A))
  b_full[valid] <- b

  list(balanced = A_bal, bias = b_full, valid = valid)
}



#' getGWmatrix_intra
#'
#' Same idea as getGWmatrix_dedicated(), but explicitly intended for intra-chromosomal blocks.
#' (In practice it is identical in your code; keeping it separate is fine if you later
#' want different behavior for intra-only matrices.)
#'
#' @param list.of.sm list like in getGWmatrix_dedicated(): per chromosome IPF outputs
#' @param ga bin annotation
#' @param which.object "LFM" or "balanced"
getGWmatrix_intra <- function(list.of.sm, ga, which.object) {

  if (which.object == "LFM") i <- 1
  if (which.object == "balanced") i <- 2

  nummSM <- Matrix(0, nrow = nrow(ga), ncol = nrow(ga))

  ids <- do.call("rbind", lapply(as.list(names(list.of.sm)), function(chrom) {

    print(paste0("including ", chrom))

    if (chrom == "chr1") {
      start.row <- 0; start.col <- 0
    }
    if (chrom != "chr1") {
      start.row <- min(which(ga$chr == chrom)) - 1
      start.col <- min(which(ga$chr == chrom)) - 1
    }

    tp <- as.data.frame(summary(list.of.sm[[chrom]][[i]]), stringsAsFactors = FALSE)
    tp$i <- tp$i + start.row
    tp$j <- tp$j + start.col

    return(tp)
  }))

  print("assembling the matrix, depending on the resolution this might take a while")

  nummSM[cbind(ids[, 1], ids[, 2])] <- ids[, 3]
  return(nummSM)
}


# ==========================
# Hi-C normalisation
# ==========================

#' loss.function
#'
#' Compute a simple convergence measure for matrix balancing:
#' sums of absolute deviations of row/col sums from 1, ignoring empty rows/cols.
#'
#' @param mat numeric matrix (dense or sparse)
#' @return numeric scalar; smaller means closer to balanced.
loss.function <- function(mat) {

  tmp <- rowSums(mat)
  tmp <- tmp[tmp > 0]
  rs <- abs(tmp - 1)

  tmp <- colSums(mat)
  tmp <- tmp[tmp > 0]
  cs <- abs(tmp - 1)

  0.5 * sum(c(rs, cs))
}


#' Create_complete_matrix
#'
#' Ensure a symmetric sparse matrix by mirroring the stored triangle.
#'
#' @param x sparse Matrix, typically storing only one triangle (or incomplete symmetry)
#' @return symmetric sparseMatrix with row/col names preserved
#'
#' @details
#' This uses summary(x) triplets and rebuilds a symmetric matrix.
Create_complete_matrix <- function(x) {

  t.sum <- as.data.frame(summary(x), stringsAsFactors = FALSE)

  test <- sparseMatrix(
    i = t.sum$i,
    j = t.sum$j,
    x = t.sum$x,
    dims = dim(x),
    symmetric = TRUE
  )

  rownames(test) <- rownames(x)
  colnames(test) <- colnames(x)

  return(test)
}


#' IPF_alg
#'
#' Iterative Proportional Fitting (IPF) to balance a matrix so that
#' row sums and column sums approach 1.
#'
#' @param lfm numeric matrix (usually the full symmetric contact matrix)
#' @param numberOfIterations integer
#'
#' @return list with:
#'  - lfm balanced matrix
#'  - x row scaling factors per iteration
#'  - y column scaling factors per iteration
#'  - lf vector of loss values over iterations
#'
#' @usage
#'  out <- IPF_alg(cm, numberOfIterations = 50)
IPF_alg <- function(lfm, numberOfIterations) {

  d <- dim(lfm)[2]

  # Store scaling factors for each iteration
  x <- matrix(1, nrow = d, ncol = numberOfIterations)
  y <- matrix(1, nrow = d, ncol = numberOfIterations)

  lf <- numeric(numberOfIterations)
  counter <- 1

  while (counter <= numberOfIterations) {

    # Row-normalize (avoid dividing by 0 for empty rows)
    s.temp <- rowSums(lfm)
    x[, counter] <- ifelse(s.temp > 0, s.temp, 1)
    lfm <- lfm / x[, counter]

    # Column-normalize
    s.temp <- colSums(lfm)
    y[, counter] <- ifelse(s.temp > 0, s.temp, 1)
    lfm <- t(t(lfm) / y[, counter])

    # Track convergence
    lf <- c(lf, loss.function(lfm))
    print(paste("loss function: ", lf[length(lf)], "counter: ", counter))

    counter <- counter + 1
  }

  return(list(lfm = lfm, x = x, y = y, lf = lf))
}


#' normalize_LFM_iteratively
#'
#' Wrapper around IPF_alg() that:
#'  - symmetrizes the input sparse matrix
#'  - runs IPF
#'  - collapses per-iteration scaling factors into a single bias per bin
#'  - applies final balancing using these biases
#'
#' @param x sparse Matrix (typically lower triangle or raw LFM sparse matrix)
#' @param numberOfIterations integer number of IPF iterations
#'
#' @return list with:
#'  - lfm balanced matrix
#'  - x row scalings per iteration
#'  - y col scalings per iteration
normalize_LFM_iteratively <- function(x, numberOfIterations) {

  # Make sure we balance a complete symmetric matrix
  cm <- Create_complete_matrix(x)

  # Run iterative proportional fitting
  tp <- IPF_alg(cm, numberOfIterations = numberOfIterations)

  x <- tp$x
  y <- tp$y

  # Convert scaling matrices to per-bin cumulative product across iterations:
  # exp(sum(log(...))) == product(...) but numerically safer.
  tp.x <- data.frame(x)
  test.x <- exp(Reduce("+", log(tp.x)))

  tp.y <- data.frame(y)
  test.y <- exp(Reduce("+", log(tp.y)))

  # Estimate a robust central ratio between x and y cumulative biases
  # (requires smhuber() to exist in your environment)
  fac <- smhuber(test.x / test.y)$mu

  # Handle bins with all-zero rows (unmappable / filtered bins)
  van <- rowSums(cm) == 0
  test.x[van] <- fac

  # Normalize bias vector by sqrt(fac) so scale is centered
  test.x <- test.x / sqrt(fac)

  # Apply balancing: divide rows and columns by bias
  tp <- cm / test.x
  tp <- t(t(tp) / test.x)

  return(list(lfm = tp, x = x, y = y))
}


#' IPF
#'
#' High-level API to produce:
#'  - raw matrix (LFM)
#'  - balanced matrix
#'  - biases (row/col scalings)
#'
#' @param x sparse Matrix (raw contacts)
#' @param numberOfIterations integer
#'
#' @return named list:
#'   - "LFM"       raw input matrix
#'   - "balanced"  balanced/symmetrized matrix
#'   - "IC_biases" list of scaling factors (x and y across iterations)
IPF <- function(x, numberOfIterations) {

  twm <- normalize_LFM_iteratively(x, numberOfIterations = numberOfIterations)

  biases <- twm[2:3]
  corrected <- twm[[1]]

  res <- vector("list", 3)
  res[[1]] <- x
  res[[2]] <- corrected
  res[[3]] <- biases

  names(res) <- c("LFM", "balanced", "IC_biases")
  return(res)
}


# ==========================
# Hi-C / loop-centric utilities
# ==========================

#' getPixelsAround
#'
#' Build a "pixel coordinate grid" around a focal pixel (ROW, COL).
#' Used to extract a square window around loop anchors from a Hi-C matrix.
#'
#' @param ROW  integer; center row index (bin index)
#' @param COL  integer; center col index (bin index)
#' @param SIZE integer; half-window size in pixels (SIZE=1 -> 3x3 window)
#' @param ID   identifier for the loop (e.g. rowname of loop table)
#'
#' @return data.frame with columns: rows, col, id
getPixelsAround <- function(ROW, COL, SIZE, ID) {

  idx <- seq(ROW - SIZE, ROW + SIZE)

  data.frame(
    rows = rep(idx, length(idx)),
    col  = rep(seq(COL - SIZE, COL + SIZE), each = length(idx)),
    id   = ID
  )
}


#' GetCentroidSignal
#'
#' For each chromosome:
#'  - subset loops on that chromosome
#'  - for each loop, build a (2*size+1)^2 set of pixels around the loop center
#'  - extract Hi-C signal for those pixels from raw or balanced matrix
#'
#' @param loops data.frame with at least:
#'   - X.chr1 (chromosome)
#'   - left_binID, right_binID (bin indices in that chromosome matrix)
#' @param hic   list of per-chromosome objects; each must contain $balanced and $LFM
#' @param size  integer half-window size in pixels
#' @param CHROMS character vector of chromosomes to process
#' @param which_file "balanced" or "raw"
#'
#' @return list (one per chromosome) of data.frames with pixel coords + hic_signal
GetCentroidSignal <- function(loops, hic, size, CHROMS, which_file) {

  lapply(as.list(CHROMS), function(chr) {

    print(chr)

    theseloops <- loops[loops$X.chr1 == chr, ]

    # Choose matrix type
    if (which_file == "balanced") this_hic <- hic[[chr]]$balanced
    if (which_file == "raw")      this_hic <- hic[[chr]]$LFM

    # Avoid window spilling outside matrix bounds
    theseloops <- theseloops[theseloops$left_binID - size > 1, ]
    theseloops <- theseloops[theseloops$right_binID + size < nrow(this_hic), ]

    # Create pixel windows for each loop and concatenate
    coords <- do.call("rbind", lapply(split(theseloops, rownames(theseloops)), function(l) {
      getPixelsAround(
        ROW  = l$left_binID,
        COL  = l$right_binID,
        SIZE = size,
        ID   = rownames(l)
      )
    }))

    # Extract signal from the matrix for each pixel coordinate
    coords$hic_signal <- this_hic[cbind(coords$rows, coords$col)]

    return(coords)
  })
}


#' ProcessLoops
#'
#' Convert the output of GetCentroidSignal() into a list of small matrices,
#' one per loop ID, suitable for aggregate plots (APA-style analyses).
#'
#' @param lo list of per-chromosome data.frames returned by GetCentroidSignal()
#' @return list where each element is a square numeric matrix (window around loop)
ProcessLoops <- function(lo) {

  tp <- do.call("rbind", lo)

  # Split by loop id
  tp <- split(tp, tp$id)

  # For each loop, reshape the linear signal vector into a square matrix
  tp <- lapply(tp, function(x) {
    matrix(x$hic_signal, nrow = sqrt(nrow(x)), ncol = sqrt(nrow(x)))
  })

  return(tp)
}


#' GetCentroidSignalGW
#'
#' Same concept as GetCentroidSignal(), but reads signal from a *genome-wide*
#' matrix (already assembled), extracting the per-chromosome submatrix via 'ga'.
#'
#' @param loops data.frame with X.chr1, left_binID, right_binID
#' @param hic   genome-wide matrix (nrow(ga) x nrow(ga))
#' @param size  half-window size (pixels)
#' @param CHROMS chromosomes to process
#' @param ga    bin annotation used to assemble the genome-wide matrix
GetCentroidSignalGW <- function(loops, hic, size, CHROMS, ga) {

  lapply(as.list(CHROMS), function(chr) {

    print(chr)

    # Extract chromosome block from genome-wide matrix
    chr_idx <- min(which(ga$chr == chr)):max(which(ga$chr == chr))
    this_hic <- hic[chr_idx, chr_idx]

    theseloops <- loops[loops$X.chr1 == chr, ]

    # Avoid boundary issues
    theseloops <- theseloops[theseloops$left_binID - size > 1, ]
    theseloops <- theseloops[theseloops$right_binID + size < nrow(this_hic), ]

    coords <- do.call("rbind", lapply(split(theseloops, rownames(theseloops)), function(l) {
      getPixelsAround(
        ROW  = l$left_binID,
        COL  = l$right_binID,
        SIZE = size,
        ID   = rownames(l)
      )
    }))

    coords$hic_signal <- this_hic[cbind(coords$rows, coords$col)]
    return(coords)
  })
}


# ==========================
# BED / NarrowPeak parsing helpers
# ==========================

#' readBed_filterChromsStraded
#'
#' Read a BED-like file, add "chr" prefix, filter to selected chromosomes,
#' and return a GRanges carrying a score column and strand.
#'
#' @param filePath path to file
#' @param chroms   chromosomes to keep (must include "chr" prefix)
#' @param scoreCol integer index of the column used as score
#'
#' @return GRanges with seqnames, ranges, strand (from column 6), score
readBed_filterChromsStraded <- function(filePath, chroms, scoreCol) {

  x <- read.delim(filePath, as.is = TRUE, header = FALSE)

  # Force UCSC style "chrN"
  x[, 1] <- paste0("chr", x[, 1])

  # Filter chromosomes
  x <- x[x[, 1] %in% chroms, ]

  GRanges(
    seqnames = Rle(x[, 1]),
    ranges   = IRanges(as.numeric(x[, 2]), end = as.numeric(x[, 3]), names = rownames(x)),
    strand   = x[, 6],
    score    = x[, scoreCol]
  )
}


#' readNarrowPeak2getSummit
#'
#' Read narrowPeak and return GRanges at the *summit* position (1bp wide).
#'
#' @param filePath path to narrowPeak file
#' @param chroms   chromosomes to keep (must include "chr" prefix)
#' @param scoreCol column index used as score
#'
#' @return GRanges where start=end=(peakStart + summitOffset)
#'
#' @details
#' narrowPeak column 10 is "peak" (summit offset from start) in ENCODE format.
readNarrowPeak2getSummit <- function(filePath, chroms, scoreCol) {

  x <- read.delim(filePath, as.is = TRUE, header = FALSE)
  x[, 1] <- paste0("chr", x[, 1])
  x <- x[x[, 1] %in% chroms, ]

  summit <- as.numeric(x[, 2]) + as.numeric(x[, 10])

  GRanges(
    seqnames = Rle(x[, 1]),
    ranges   = IRanges(summit, end = summit, names = rownames(x)),
    score    = x[, scoreCol]
  )
}


#' appendMotifInformation
#'
#' Attach motif overlap information to a set of peaks:
#'  - count motifs within ±25bp of the summit
#'  - store the best motif per peak (max motif score), plus strand and position
#'
#' @param peakFile   GRanges of peaks (will be returned with extra metadata columns)
#' @param motifFile  GRanges of motif hits with a 'score' metadata column and strand
#' @param summitFile GRanges of peak summits; will be expanded ±25bp internally
#'
#' @return peakFile with new columns:
#'   number_of_motifs, motif_score, motif_strand, motif_start, motifID
appendMotifInformation <- function(peakFile, motifFile, summitFile) {

  # Restrict motif search to a small window around the summit (±25 bp)
  start(summitFile) <- start(summitFile) - 25
  end(summitFile)   <- end(summitFile) + 25

  # Count how many motifs overlap each summit window
  peakFile$number_of_motifs <- countOverlaps(summitFile, motifFile)

  # Initialize motif annotations
  peakFile$motif_score  <- 0
  peakFile$motif_strand <- factor(rep("*", length(peakFile)), levels = c("*", "+", "-"))

  # For motifs on '-' strand, you may want a consistent "start" coordinate:
  # this sets realStart to motif end for '-' (often used as a strand-aware anchor).
  motifFile$realStart <- start(motifFile)
  motifFile$realStart[which(strand(motifFile) == "-")] <- end(motifFile[which(strand(motifFile) == "-")])

  motifFile$motifID <- 0

  # Find overlaps between summit windows and motifs
  peak_motif <- findOverlaps(summitFile, motifFile)

  # Build a table of motif attributes for each overlap
  tp <- data.frame(
    score       = motifFile$score,
    strand      = strand(motifFile),
    motif_start = motifFile$realStart,
    motif_id    = seq_along(motifFile),
    stringsAsFactors = FALSE
  )

  # Split by peak (query) and keep only motifs overlapping each peak
  tp <- split(tp[subjectHits(peak_motif), ], queryHits(peak_motif))

  # For each peak, pick the motif hit with maximum score
  th <- lapply(tp, function(x) x[which.max(x$score), ])

  # Store best-motif annotations back into peakFile
  peakFile$motif_score[as.numeric(names(th))]  <- unlist(lapply(th, function(f) f$score))
  peakFile$motif_strand[as.numeric(names(th))] <- unlist(lapply(th, function(f) f$strand))
  peakFile$motif_start[as.numeric(names(th))]  <- unlist(lapply(th, function(f) f$motif_start))
  peakFile$motifID[as.numeric(names(th))]      <- unlist(lapply(th, function(f) f$motif_id))

  return(peakFile)
}


# ==========================
# Loop signal statistics
# ==========================

#' compute_putative_loop_signal
#'
#' Compute a simple loop enrichment metric from a small square matrix centered on the loop:
#'  - core signal: mean of a central square of size (2*central_sq+1)
#'  - background: mean of 15 randomly sampled 5x5 windows elsewhere in the matrix
#'  - outputs log2 fold-change and a t-test p-value (background vs core as mean)
#'
#' @param ctcf_ctcf_loop_signal numeric matrix (e.g. output of ProcessLoops() per loop)
#' @param central_sq integer; half-size of the core window around the center
#'
#' @return named vector with:
#'   - interact_signal : core mean
#'   - l2fc_v.bg       : log2(core / mean(background))
#'   - p.value         : p-value from t.test(background, mu=core)
compute_putative_loop_signal <- function(ctcf_ctcf_loop_signal, central_sq) {

  signal_samples <- c()

  # Central window indices (around matrix midpoint)
  mid <- ceiling(ncol(ctcf_ctcf_loop_signal) / 2)
  core_range <- (mid - central_sq):(mid + central_sq)

  # Core mean signal
  core_signal <- mean(ctcf_ctcf_loop_signal[core_range, core_range])

  # Avoid sampling background squares overlapping core region
  no_go_range <- (min(core_range) - 4):max(core_range)

  check_indexs <- 1:(ncol(ctcf_ctcf_loop_signal) - 4)
  check_indexs <- check_indexs[!(check_indexs %in% no_go_range)]

  # Sample 15 background windows (top-left corners for 5x5 squares)
  sample_rows <- sample(check_indexs, 15, replace = TRUE)
  sample_cols <- sample(check_indexs, 15, replace = TRUE)

  for (i in 1:15) {
    signal_samples <- c(
      signal_samples,
      mean(ctcf_ctcf_loop_signal[
        sample_rows[i]:(sample_rows[i] + 4),
        sample_cols[i]:(sample_cols[i] + 4)
      ])
    )
  }

  return(c(
    interact_signal = core_signal,
    l2fc_v.bg       = log2(core_signal / mean(signal_samples)),
    t.test(signal_samples, mu = core_signal)[3]
  ))
}


#' getMatrixForARegion
#'
#' Build a "hybrid" matrix for plotting comparisons between conditions:
#'  - lower triangle from one matrix (lowertri)
#'  - upper triangle from another matrix (uppertri)
#'
#' Commonly used to visualize WT vs KO in one heatmap:
#'  - WT on one triangle, KO on the other.
#'
#' @param lowertri list of per-chromosome matrices with $balanced
#' @param uppertri list of per-chromosome matrices with $balanced
#' @param CHROM chromosome to extract
#' @param START genomic start coordinate
#' @param END genomic end coordinate
#' @param ga bin annotation (binAnno output)
#' @param upperLimit numeric; cap values for visualization
#'
#' @return numeric matrix for the selected region, flipped for plotting
getMatrixForARegion <- function(lowertri, uppertri, CHROM, START, END, ga, upperLimit) {

  # Find bins fully contained in [START, END]
  bins <- ga[ga$chr == CHROM, ]
  bins <- which(bins$start > START & bins$end < END)

  # Extract region from each matrix
  a <- lowertri[[CHROM]]$balanced
  b <- uppertri[[CHROM]]$balanced
  a <- a[min(bins):max(bins), min(bins):max(bins)]
  b <- b[min(bins):max(bins), min(bins):max(bins)]

  # Combine: lower triangle from a, upper triangle from b
  ab <- a
  ab[upper.tri(ab)] <- b[upper.tri(b)]
  diag(ab) <- 0

  # Flip vertically (common for image-like plotting in R)
  ab <- t(matrix(apply(ab, 2, rev), nrow(ab), ncol(ab)))

  # Cap extreme values for nicer contrast in plots
  ab[ab > upperLimit] <- upperLimit

  return(ab)
}


#' getSignalInBins
#'
#' Assign a genomic signal (e.g., ChIP/ATAC score) to bins by overlap.
#' For each bin, keeps the maximum overlapping feature score.
#'
#' @param bins GRanges of bins (often a tiled window around peaks/anchors)
#' @param modification GRanges of features with a numeric metadata column
#' @param featureColumn name or index of metadata column in 'modification'
#'
#' @return numeric matrix:
#'   rows correspond to "groups" (e.g. one row per peak)
#'   columns correspond to bins per group (sum(bins$peak==1) in your encoding)
#'
#' @details
#' This assumes your 'bins' object encodes repeated groups and includes a
#' bins$peak==1 marker to define group width. Keep as-is because it matches
#' your downstream code.
getSignalInBins <- function(bins, modification, featureColumn) {

  PsDFvsX <- findOverlaps(bins, modification)

  # Collect overlapping feature scores per bin
  PsDFvsX_score <- split(
    as.numeric(elementMetadata(modification)[, featureColumn][subjectHits(PsDFvsX)]),
    queryHits(PsDFvsX)
  )

  # Use max score per bin (peak-like signal)
  PsDFvsX_score <- unlist(lapply(PsDFvsX_score, function(x) max(x, na.rm = TRUE)))

  bm <- rep(0, length(bins))
  bm[as.numeric(names(PsDFvsX_score))] <- PsDFvsX_score

  # Reshape into a matrix: one row per group, columns per bin within group
  bm <- matrix(
    bm,
    nrow = length(bins) / sum(bins$peak == 1),
    ncol = sum(bins$peak == 1),
    byrow = TRUE
  )

  rownames(bm) <- names(bins)[seq(1, length(bins), by = sum(bins$peak == 1))]

  return(bm)
}


#' getLoopSignal
#'
#' Convenience wrapper: compute loop enrichment (log2 FC vs background) for many loops.
#'
#' @param ll list of per-loop matrices (e.g. output of ProcessLoops())
#' @param D integer; central half-window size passed to compute_putative_loop_signal()
#'
#' @return numeric vector of l2fc values, one per loop
getLoopSignal <- function(ll, D) {

  tp <- lapply(ll, function(x) compute_putative_loop_signal(x, D))

  # Return the second element ("l2fc_v.bg") for each loop
  unlist(lapply(tp, function(x) x[[2]]))
}


#' importBEDPE_Filter_Loops
#'
#' Read HiCCUPS (or similar) loop calls in BEDPE, filter by max span,
#' and optionally extend anchors by a fixed number of bp.
#'
#' @param loopdir directory containing subfolders with merged_loops.bedpe
#' @param max_size numeric; keep loops with span < max_size
#' @param extend_anchors_by numeric; extend x1/x2/y1/y2 by this many bp
#'
#' @return data.frame with loops concatenated across subfolders
importBEDPE_Filter_Loops <- function(loopdir, max_size, extend_anchors_by) {

  do.call("rbind", lapply(as.list(list.files(loopdir)), function(x) {

    tp <- read.delim(paste0(loopdir, x, "/merged_loops.bedpe"), header = TRUE)

    # Drop the first row if it is a duplicated header or metadata row (as in your original code)
    tp <- tp[-1, ]

    tp$loop_span <- tp$y2 - tp$x1
    tp <- tp[tp$loop_span < max_size, ]

    # Extend anchors
    tp$x1 <- tp$x1 - extend_anchors_by
    tp$x2 <- tp$x2 + extend_anchors_by
    tp$y1 <- tp$y1 - extend_anchors_by
    tp$y2 <- tp$y2 + extend_anchors_by

    tp
  }))
}


#' compareLoops
#'
#' Compare two sets of loops and split them into common vs specific.
#' A loop is considered overlapping if BOTH anchors overlap (left and right).
#'
#' @param loops1 data.frame with columns: chr1, x1, x2, chr2, y1, y2 (in this order)
#' @param loops2 same format as loops1
#' @param offset bp extension used for overlap (tolerance)
#' @param maxLoopSize maximum allowed span (filter applied to both sets)
#'
#' @return named list with 4 data.frames:
#'   - First_set_common
#'   - First_set_specific
#'   - Second_set_common
#'   - Second_set_specific
compareLoops <- function(loops1, loops2, offset, maxLoopSize) {

  loops1 <- loops1[loops1$loop_span < maxLoopSize, ]
  loops2 <- loops2[loops2$loop_span < maxLoopSize, ]

  # Build GRanges for anchors with optional tolerance (offset)
  leftGR1  <- getGR(loops1[, 1], loops1[, 2], loops1[, 3], offset = offset)
  leftGR2  <- getGR(loops2[, 1], loops2[, 2], loops2[, 3], offset = offset)
  rightGR1 <- getGR(loops1[, 4], loops1[, 5], loops1[, 6], offset = offset)
  rightGR2 <- getGR(loops2[, 4], loops2[, 5], loops2[, 6], offset = offset)

  # Overlaps for left and right anchors
  leftOV  <- findOverlaps(leftGR1, leftGR2)
  rightOV <- findOverlaps(rightGR1, rightGR2)

  # A loop pair is "common" if both left and right overlaps exist for same (query,subject)
  overlapping1 <- table(c(
    paste(queryHits(leftOV),  subjectHits(leftOV),  sep = "_"),
    paste(queryHits(rightOV), subjectHits(rightOV), sep = "_")
  ))

  overlapping_loops1 <- unique(as.numeric(
    unlist(strsplit(names(overlapping1[overlapping1 == 2]), "_"))[
      seq(1, 2 * length(overlapping1[overlapping1 == 2]), by = 2)
    ]
  ))

  overlapping_loops1 <- overlapping_loops1[order(overlapping_loops1)]

  loops1_common   <- loops1[overlapping_loops1, ]
  loops1_specific <- loops1[-overlapping_loops1, ]

  # Repeat in the opposite direction to get loops2 common/specific
  leftOV  <- findOverlaps(leftGR2, leftGR1)
  rightOV <- findOverlaps(rightGR2, rightGR1)

  overlapping2 <- table(c(
    paste(queryHits(leftOV),  subjectHits(leftOV),  sep = "_"),
    paste(queryHits(rightOV), subjectHits(rightOV), sep = "_")
  ))

  overlapping_loops2 <- unique(as.numeric(
    unlist(strsplit(names(overlapping2[overlapping2 == 2]), "_"))[
      seq(1, 2 * length(overlapping2[overlapping2 == 2]), by = 2)
    ]
  ))

  loops2_common   <- loops2[overlapping_loops2, ]
  loops2_specific <- loops2[-overlapping_loops2, ]

  res <- vector("list", length = 4)
  res[[1]] <- loops1_common
  res[[2]] <- loops1_specific
  res[[3]] <- loops2_common
  res[[4]] <- loops2_specific

  names(res) <- c(
    "First_set_common", "First_set_specific",
    "Second_set_common", "Second_set_specific"
  )

  return(res)
}


# ==========================
# Insulation analysis
# ==========================

#' getSUMMEDsignal4bins
#'
#' Helper for insulation score: for each bin index B, sum signal in a square region
#' offset by distance D and half-width A.
#'
#' @param B integer vector of bin indices (named)
#' @param D integer offset distance (in bins)
#' @param A integer half-size of the square window
#' @param M matrix from which to extract values
#'
#' @return numeric vector; sum per bin (names propagated from B)
getSUMMEDsignal4bins <- function(B, D, A, M) {

  # Build row/col indices for (2A+1)x(2A+1) squares per B
  rows <- unlist(lapply(as.list(B), function(b) {
    rep(((b - D) + seq(-A, A)), (1 + 2 * A))
  }))

  cols <- unlist(lapply(as.list(B), function(b) {
    rep(((b + D) + seq(-A, A)), each = (1 + 2 * A))
  }))

  ids <- cbind(rows, cols)

  # Extract matrix values (vectorized)
  m <- M[ids]

  IDS <- rep(names(B), each = (1 + 2 * A)^2)
  matrices <- split(m, IDS)

  unlist(lapply(matrices, function(m) sum(m, na.rm = TRUE)))
}


#' InsulationScore
#'
#' Compute insulation components (left/middle/right) around bin centers for each chromosome,
#' using the *balanced* matrix stored in mat[[chr]]$balanced.
#'
#' @param bin GRanges of positions where you want insulation (often bins or boundary candidates)
#' @param mat list of per-chromosome objects with $balanced matrix
#' @param GAGR GRanges of all bins with metadata column 'binid'
#' @param distance integer in bins
#' @param Area integer half-size of square window
#' @param HOWFAR integer shift for left/right windows (in bins)
#'
#' @return data.frame with columns left/middle/right for all processed bins (stacked)
InsulationScore <- function(bin, mat, GAGR, distance, Area, HOWFAR) {

  do.call("rbind", lapply(as.list(unique(as.character(chrom(bin)))), function(chr) {

    print(paste0("processing ", chr))

    thism <- mat[[chr]]$balanced

    # Use 1bp centers for bin positions
    theseb <- GenomicRanges::resize(bin[which(chrom(bin) == chr)], 1, fix = "center")

    # Map GRanges positions to bin IDs using overlaps with the global bin GRanges
    bins <- GAGR$binid[queryHits(findOverlaps(GAGR, theseb))]
    names(bins) <- names(theseb[subjectHits(findOverlaps(GAGR, theseb))])

    # Filter out bins too close to matrix edges
    bins <- bins[
      bins - distance - HOWFAR - Area > 0 &
      bins + distance + HOWFAR + Area < (length(GAGR[which(chrom(GAGR) == chr)]))
    ]

    data.frame(
      left   = getSUMMEDsignal4bins(B = bins - HOWFAR, D = distance, A = Area, M = thism),
      middle = getSUMMEDsignal4bins(B = bins,          D = distance, A = Area, M = thism),
      right  = getSUMMEDsignal4bins(B = bins + HOWFAR, D = distance, A = Area, M = thism)
    )
  }))
}


#' InsulationScorePublishedKR
#'
#' Variant of InsulationScore() that assumes mat[[chr]] is already a matrix
#' (e.g., published KR-normalized matrices), not a list containing $balanced.
InsulationScorePublishedKR <- function(bin, mat, GAGR, distance, Area, HOWFAR) {

  do.call("rbind", lapply(as.list(unique(as.character(chrom(bin)))), function(chr) {

    print(paste0("processing ", chr))

    thism <- mat[[chr]]
    theseb <- GenomicRanges::resize(bin[which(chrom(bin) == chr)], 1, fix = "center")

    bins <- GAGR$binid[queryHits(findOverlaps(GAGR, theseb))]
    names(bins) <- names(theseb[subjectHits(findOverlaps(GAGR, theseb))])

    bins <- bins[
      bins - distance - HOWFAR - Area > 0 &
      bins + distance + HOWFAR + Area < (length(GAGR[which(chrom(GAGR) == chr)]))
    ]

    data.frame(
      left   = getSUMMEDsignal4bins(B = bins - HOWFAR, D = distance, A = Area, M = thism),
      middle = getSUMMEDsignal4bins(B = bins,          D = distance, A = Area, M = thism),
      right  = getSUMMEDsignal4bins(B = bins + HOWFAR, D = distance, A = Area, M = thism)
    )
  }))
}


#' processIS
#'
#' Convert left/middle/right insulation components into a single insulation score:
#'   score = log2( mean(left,right) / middle )
#'
#' @param IS data.frame/matrix with columns 1=left, 2=middle, 3=right
#' @param GAGR GRanges used as a template for output positions (must have names matching rownames(IS))
#'
#' @return GRanges with 'score' metadata set for bins found in IS
processIS <- function(IS, GAGR) {

  res <- GAGR
  res$binid <- NULL

  res$score <- 0
  res$score[match(rownames(IS), names(res))] <-
    log2(rowMeans((0.001 + IS[, c(1, 3)])) / (0.001 + IS[, 2]))

  return(res)
}


#' cleanIS
#'
#' Replace NA / Inf values with 0 in insulation score vectors/matrices.
cleanIS <- function(x) {
  x[is.na(x)] <- 0
  x[!is.finite(x)] <- 0
  x
}


# ==========================
# RNA-seq and ChIP-seq functions
# ==========================

#' getSignalInBins400
#'
#' Similar to getSignalInBins(), but uses MEAN overlap score per bin and reshapes
#' into a matrix with exactly 401 columns (your fixed window encoding).
#'
#' @param bins GRanges of bins (length must be multiple of 401)
#' @param modification GRanges of features with numeric metadata column
#' @param featureColumn metadata column name/index to extract
#'
#' @return numeric matrix (nrow = length(bins)/401, ncol = 401)
getSignalInBins400 <- function(bins, modification, featureColumn) {

  PsDFvsX <- findOverlaps(bins, modification)

  PsDFvsX_score <- split(
    as.numeric(elementMetadata(modification)[, featureColumn][subjectHits(PsDFvsX)]),
    queryHits(PsDFvsX)
  )

  # Use mean score per bin (coverage-like signal)
  PsDFvsX_score <- unlist(lapply(PsDFvsX_score, function(x) mean(x, na.rm = TRUE)))

  bm <- rep(0, length(bins))
  bm[as.numeric(names(PsDFvsX_score))] <- PsDFvsX_score

  bm <- matrix(bm, nrow = length(bins) / 401, ncol = 401, byrow = TRUE)
  rownames(bm) <- (bins$peak)[seq(1, length(bins), by = 401)]

  return(bm)
}


#' reorderMatrixBasedOnStrand
#'
#' Reorder and strand-correct a signal matrix around motifs:
#'  - '+' motifs keep original orientation
#'  - '-' motifs are reversed so motif direction is consistent
#'
#' @param m numeric matrix (rows = peaks, cols = window bins)
#' @param peaks GRanges/data.frame with column motif_strand in {"+","-"}
#'
#' @return numeric matrix with '-' rows reversed and appended after '+' rows
reorderMatrixBasedOnStrand <- function(m, peaks) {

  m_fwd <- m[which(as.character(peaks$motif_strand) == "+"), ]
  m_rev <- t(apply(m[which(as.character(peaks$motif_strand) == "-"), ], 1, rev))

  rbind(m_fwd, m_rev)
}


#' GetTPM
#'
#' Compute TPM (Transcripts Per Million) from raw counts and gene lengths.
#'
#' @param m data.frame with at least:
#'   - m$Length gene length (bp)
#'   - count columns specified in cols2take
#' @param cols2take integer vector of columns to use as raw counts
#' @param gn character vector of gene names for rownames of output
#'
#' @return numeric matrix of TPM values with genes as rows and samples as columns
GetTPM <- function(m, cols2take, gn) {

  mb <- m[, cols2take]
  cs <- colSums(mb) # not used downstream, but kept from original code

  f1 <- do.call("cbind", lapply(seq_len(ncol(mb)), function(x) {

    # Reads Per Kilobase (RPK)
    rpk <- 1000 * (mb[, x] / m$Length)

    # Scaling factor: sum(RPK) / 1e6
    SF <- sum(rpk) / 1000000

    # TPM
    rpk / SF
  }))

  colnames(f1) <- colnames(mb)
  rownames(f1) <- gn

  return(f1)
}


#' readBed_filterChroms
#'
#' Read a BED-like file, add "chr" prefix, filter chromosomes,
#' and return GRanges with a score column (no strand).
readBed_filterChroms <- function(filePath, chroms, scoreCol) {

  x <- read.delim(filePath, as.is = TRUE, header = FALSE)
  x[, 1] <- paste0("chr", x[, 1])
  x <- x[x[, 1] %in% chroms, ]

  GRanges(
    seqnames = Rle(x[, 1]),
    ranges   = IRanges(as.numeric(x[, 2]), end = as.numeric(x[, 3]), names = rownames(x)),
    score    = x[, scoreCol]
  )
}

# NOTE: Your snippet ends with duplicate definitions of:
#  - readBed_filterChromsStraded
#  - readNarrowPeak2getSummit
# If these are in the same script, the later definitions will overwrite earlier ones.
# Consider keeping only one copy to avoid confusion.
