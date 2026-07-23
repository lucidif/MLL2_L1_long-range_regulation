# ============================================================
# Utilities for:
# 1) Reading coolpup exported matrices and reshaping to long format
# 2) Handling genomic loci and binning (0-based, fixed bin size)
# 3) Extracting Hi-C contacts from .hic (Juicer) via strawr::straw()
#
# Coordinate conventions:
# - We assume BEDPE-like intervals are 0-based half-open [start, end)
#   (typical BED convention). The "viewpoint" coordinate is derived from
#   these inputs.
# - Juicer .hic queries via straw() use chr:start:end strings. In practice,
#   using an "end inclusive" convention avoids pulling the next bin when
#   end falls exactly on a bin boundary; this is why we set win_end <- win_end - 1.
# ============================================================


# ----------------------------
# coolpop_MA_plot_in()
# ----------------------------
# Read two TSV matrices (e.g., coolpup avg matrices) and convert them into a
# merged long-format dataframe keyed by pixel (binX, binY).
#
# Inputs:
# - avema1path: path to TSV matrix for condition 1 (here assumed WT)
# - avema2path: path to TSV matrix for condition 2 (here assumed KO)
#
# Output:
# - A dataframe where each row is a matrix pixel with:
#   unid, binX, binY, WT_ave, KO_ave
#
# Notes:
# - read.table() -> matrix -> as.table() -> as.data.frame() converts a dense
#   matrix into a long table (binX/binY/value).
# - binX/binY come out as factors; merge() will still work, but you can
#   convert to integer if needed after loading.
coolpop_MA_plot_in <- function(avema1path, avema2path) {

  # Read condition 1 matrix (TSV -> matrix -> long dataframe)
  ma1 <- as.table(as.matrix(read.table(avema1path, sep = "\t")))
  ma1 <- as.data.frame(ma1)
  colnames(ma1) <- c("binX", "binY", "WT_ave")

  # Read condition 2 matrix
  ma2 <- as.table(as.matrix(read.table(avema2path, sep = "\t")))
  ma2 <- as.data.frame(ma2)
  colnames(ma2) <- c("binX", "binY", "KO_ave")

  # Optional sanity check (prints first rows to console)
  head(ma1)
  head(ma2)

  # Create a unique pixel ID so we can merge the two conditions safely
  ma1 <- cbind(unid = paste0(ma1$binX, "_", ma1$binY), ma1)
  ma2 <- cbind(unid = paste0(ma2$binX, "_", ma2$binY), ma2)

  # Merge by pixel ID
  ma_all <- merge(ma1, ma2, by = "unid")

  return(ma_all)
}


# ----------------------------
# get_chr_len()
# ----------------------------
# Utility: given a chromosome name and a chromsize table (V1=chr, V2=len),
# return the chromosome length.
#
# Inputs:
# - chr: chromosome string, e.g. "chr1"
# - chromsize: data.frame with columns V1 (chr) and V2 (len)
#
# Output:
# - integer chromosome length (or NA if chr not found)
get_chr_len <- function(chr, chromsize) {
  chromsize$V2[match(chr, chromsize$V1)]
}


# ----------------------------
# pick_viewpoint()
# ----------------------------
# Choose a single "viewpoint" coordinate (bp) from an interval [start, end).
#
# Inputs:
# - start, end: interval boundaries (expected 0-based BED-style)
# - method: "center", "start", or "end"
#
# Output:
# - integer coordinate (bp) used as the viewpoint
#
# Notes:
# - "center" uses floor((start + end) / 2), which is stable for even/odd spans.
pick_viewpoint <- function(start, end, method = c("center", "start", "end")) {
  method <- match.arg(method)
  if (method == "center") return(as.integer(floor((start + end) / 2)))
  if (method == "start")  return(as.integer(start))
  as.integer(end)
}


# ----------------------------
# coord_to_bin()
# ----------------------------
# Convert a genomic coordinate (bp) into a fixed-size bin index and bin boundaries.
#
# Inputs:
# - coord: genomic coordinate in bp (integer)
# - binsize: bin size in bp (e.g. 5000)
#
# Output:
# - list(bin_id, bin_start, bin_end)
#
# Conventions:
# - bin_id is 0-based: bin_id = floor(coord / binsize)
# - bin_start = bin_id * binsize
# - bin_end   = bin_start + binsize  (end boundary, half-open)
coord_to_bin <- function(coord, binsize) {
  bin_id <- as.integer(floor(coord / binsize))  # 0-based
  bin_start <- bin_id * binsize
  bin_end   <- bin_start + binsize
  list(bin_id = bin_id, bin_start = bin_start, bin_end = bin_end)
}


# ----------------------------
# make_window_locus()
# ----------------------------
# Build a juicer-style locus chr:start:end around the viewpoint bin, expanded by
# a symmetric window of +/- window_bins bins.
#
# Example:
# - binsize=5000, window_bins=5 => total bins per axis = 2*5+1 = 11
#
# Inputs:
# - chr, start, end: interval defining the region (we use it only to get viewpoint)
# - chromsize: table with chr lengths (V1=chr, V2=len)
# - binsize: bin size in bp
# - viewpoint: "center" | "start" | "end"
# - window_bins: number of bins to include on each side of the viewpoint bin
#
# Output:
# - list with:
#   chr, vp (viewpoint bp), bin_id (viewpoint bin index),
#   locus_vec = c(chr, win_start, win_end)
#
# Notes:
# - We clamp win_start/win_end to chromosome boundaries.
# - IMPORTANT: we set win_end <- win_end - 1 to avoid including the next bin
#   when win_end lands exactly on a bin boundary (common off-by-one in .hic queries).
make_window_locus <- function(chr, start, end, chromsize, binsize,
                              viewpoint = c("center", "start", "end"),
                              window_bins = 10) {

  viewpoint <- match.arg(viewpoint)

  # Chromosome length (integer)
  chr_len <- get_chr_len(chr, chromsize)
  if (is.na(chr_len)) stop("Chrom not found in chromsize: ", chr)

  # Choose viewpoint coordinate and map it to a bin
  vp <- pick_viewpoint(start, end, viewpoint)
  b  <- coord_to_bin(vp, binsize)

  # Window around the viewpoint bin
  win_start <- b$bin_start - window_bins * binsize
  win_end   <- b$bin_end   + window_bins * binsize

  # Make end boundary inclusive (avoid accidental extra bin at boundary)
  win_end <- win_end - 1

  # Clamp to chromosome boundaries (0-based)
  win_start <- max(0, win_start)
  win_end   <- min(chr_len, win_end)

  list(
    chr = chr,
    vp = vp,
    bin_id = b$bin_id,
    locus_vec = c(chr, win_start, win_end)  # ready for paste(..., collapse=":")
  )
}


# ============================================================
# extract_contacts()
# ============================================================
# Main function: for each BEDPE pair, query a .hic file using strawr::straw()
# on a window around the viewpoint bin. Returns either:
# - a list of per-pair sparse contact tables (straw output), OR
# - a dataframe with per-pair summed contact strength (sum of counts), with NA
#   when no contacts are returned for that query.
#
# Inputs:
# - hicpath: full path to a .hic file
# - bedpe.path: path to a 6-column BEDPE-like file:
#     chrA startA endA chrB startB endB
# - window_bins: +/- bins around viewpoint bin (per axis)
# - viewpoint: "center" | "start" | "end" for viewpoint selection
# - bin.size: binsize passed to straw(); must exist in the .hic resolutions
# - chromsize.file: path to chromosome sizes file (chr, length)
# - return_mode: "list" or "sum_df"
# - verbose: if TRUE, print per-entry debug info to console
#
# Output:
# - if return_mode == "sum_df":
#     data.frame(anchor, bait, value)
#   where anchor/bait are the exact loci queried (chr:start:end),
#   and value is sum(imported_hic$counts) or NA if no rows returned.
#
# - if return_mode == "list":
#     list of entries, each with metadata + contacts (sparse table).
#     Entries with no contacts are stored as NULL.
#
# Notes:
# - straw() returns a sparse table with columns x, y, counts where x/y are
#   bin start coordinates (multiples of binsize) and counts are normalized
#   if norm="KR" (can be non-integers).
extract_contacts <- function(hicpath, bedpe.path,
                             window_bins = 11,
                             viewpoint = "center",
                             bin.size = 5000,
                             chromsize.file,
                             return_mode = c("list", "sum_df"),
                             verbose = FALSE) {

  return_mode <- match.arg(return_mode)

  # Ensure strawr is available (straw() comes from strawr)
  if (!requireNamespace("strawr", quietly = TRUE)) {
    stop("Package 'strawr' not installed. Install with install.packages('strawr').")
  }

  # Read BEDPE robustly: chr columns as character, coordinates as integer
  ach3 <- read.table(
    bedpe.path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    colClasses = c("character","integer","integer","character","integer","integer")
  )

  # Read chrom sizes robustly
  chromsize <- read.table(
    chromsize.file,
    header = FALSE,
    stringsAsFactors = FALSE,
    colClasses = c("character", "integer")
  )

  # Pre-allocate outputs
  res_list <- vector("list", nrow(ach3))

  out_df <- data.frame(
    anchor = character(nrow(ach3)),
    bait   = character(nrow(ach3)),
    value  = rep(NA_real_, nrow(ach3)),
    stringsAsFactors = FALSE
  )

  validr <- 0

  for (i in 1:nrow(ach3)) {

    # Anchor and bait intervals
    anch <- ach3[i, c(1,2,3)]
    bait <- ach3[i, c(4,5,6)]

    # Build window loci around viewpoint bin
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

    # Convert to juicer locus strings
    anch_loc <- paste(anch_win$locus_vec, collapse = ":")
    bait_loc <- paste(bait_win$locus_vec, collapse = ":")

    # Query .hic (KR = balanced/ICE-like normalization)
    imported_hic <- strawr::straw(
      norm    = "KR",
      fname   = hicpath,
      chr1loc = anch_loc,
      chr2loc = bait_loc,
      unit    = "BP",
      binsize = bin.size
    )

    # Always fill dataframe output (NA if no contacts returned)
    out_df$anchor[i] <- anch_loc
    out_df$bait[i]   <- bait_loc
    out_df$value[i]  <- if (nrow(imported_hic) > 0) sum(imported_hic$counts) else NA_real_

    # Fill list output only if there are contacts
    if (nrow(imported_hic) > 0) {

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

  # Return requested output type
  if (return_mode == "sum_df") {
    return(out_df)
  } else {
    return(res_list)
  }
}
