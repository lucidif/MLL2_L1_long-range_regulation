suppressPackageStartupMessages({
  library(HicAggR)
  library(GenomicRanges)
  library(InteractionSet)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(IRanges)
  library(zoo)
  library(cowplot)
})

# -----------------------------------------------------------------------------
# Helper 1: call borders from an insulation bigWig
#   - imports insulation track
#   - optional NA interpolation + light smoothing
#   - calls local minima
#   - returns ~bin-wide windows centered on borders, with name/class columns
# -----------------------------------------------------------------------------
hicaggr_borders_from_bigwig <- function(
  bw_path,
  k = 5,                 # smoothing window for zoo::rollmedian
  binSize = NULL,        # width of returned windows; defaults to median bin width
  interpolate = TRUE,    # fill internal NA gaps only (ends remain NA)
  min_separation = 0L    # merge boundaries closer than this distance (bp)
){
  # 1) Import insulation track and sanitize scores (NaN -> NA)
  gr <- rtracklayer::import(bw_path, as = "GRanges")
  sc <- S4Vectors::mcols(gr)$score
  sc[is.nan(sc)] <- NA_real_
  S4Vectors::mcols(gr)$score <- sc
  gr <- GenomicRanges::sort(gr)  # ensure increasing genomic order

  # Infer bin size from the track if not provided
  if (is.null(binSize)) {
    binSize <- as.integer(stats::median(IRanges::width(gr), na.rm = TRUE))
  } else {
    binSize <- as.integer(binSize)
  }

  # Helper: interpolate only internal NA gaps (leave ends as NA)
  interp_internal_na <- function(x) {
    if (!interpolate) return(x)
    if (all(!is.finite(x))) return(x)
    idx <- which(is.finite(x))
    if (length(idx) < 2) return(x)               # not enough points to interpolate
    xx <- seq_along(x)
    y  <- x
    yhat <- approx(x = xx[idx], y = x[idx], xout = xx, rule = 2)$y
    left  <- min(idx); right <- max(idx)
    fill  <- (xx >= left & xx <= right & !is.finite(x))
    y[fill] <- yhat[fill]
    y[xx < left | xx > right] <- NA_real_        # keep leading/trailing NAs
    y
  }

  # 2) Call local minima per chromosome
  borders_list <- vector("list", length = length(GenomeInfoDb::seqlevels(gr)))
  names(borders_list) <- GenomeInfoDb::seqlevels(gr)

  for (chr in GenomeInfoDb::seqlevels(gr)) {
    g <- gr[GenomicRanges::seqnames(gr) == chr]
    if (length(g) == 0L) next

    s  <- as.numeric(S4Vectors::mcols(g)$score)
    s2 <- interp_internal_na(s)
    if (sum(is.finite(s2)) < 3L) next             # need at least 3 finite points

    # Light smoothing; requires zoo
    s_sm <- zoo::rollmedian(s2, k = k, fill = NA)

    # Find strict local minima relative to immediate neighbors
    lag  <- c(NA_real_, s_sm[-length(s_sm)])
    lead <- c(s_sm[-1], NA_real_)
    idx  <- which(is.finite(s_sm) & s_sm < lag & s_sm < lead)
    if (!length(idx)) next

    # Convert bins to integer midpoints
    centers <- as.integer(IRanges::start(g)[idx] + (IRanges::width(g)[idx] %/% 2L))

    # Optionally merge boundaries that are too close (cluster & take midpoints)
    if (min_separation > 0L && length(centers) > 1L) {
      centers <- sort(centers)
      blocks  <- IRanges::reduce(IRanges::IRanges(start = centers, width = 1L),
                                 min.gapwidth = as.integer(min_separation))
      centers <- as.integer(IRanges::start(blocks) + (IRanges::width(blocks) %/% 2L))
    }

    borders_list[[chr]] <- GenomicRanges::GRanges(seqnames = chr,
                                                  ranges   = IRanges::IRanges(centers, width = 1L))
  }

  # 3) Concatenate all chromosomes and sort
  nonempty <- Filter(function(x) methods::is(x, "GRanges") && length(x) > 0L, borders_list)
  borders  <- if (length(nonempty) == 0L) {
    GenomicRanges::GRanges()
  } else {
    unlist(GenomicRanges::GRangesList(nonempty), use.names = FALSE)
  }
  borders <- GenomicRanges::sort(borders)

  if (length(borders) == 0L) {
    stop("No boundaries called. Try lowering 'k', enabling 'interpolate', or reducing 'min_separation'.")
  }

  # 4) Return ~bin-sized windows centered at boundaries, with minimal annotation
  windows <- GenomicRanges::resize(borders, width = binSize, fix = "center")
  S4Vectors::mcols(windows)$name  <- paste0("BORDER_", seq_along(windows))
  S4Vectors::mcols(windows)$class <- "border"
  windows
}

# Build TAD domains from border windows
# - Uses boundary centers by default
# - Snaps start/end to the Hi-C bin grid (binSize)
# - Extends first/last domain to chromosome ends
# Build non-overlapping TADs on a fine bin grid using half-open intervals (cut[i], cut[i+1]]
borders_to_domains <- function(borders, chromSizes, binSize = 5000L,
                               boundary_anchor = c("right_edge","left_edge","center"),
                               min_width = binSize) {
  boundary_anchor <- match.arg(boundary_anchor)
  stopifnot(is.data.frame(chromSizes),
            all(c("seqnames","seqlengths") %in% names(chromSizes)))
  chromSizes$seqnames   <- as.character(chromSizes$seqnames)
  chromSizes$seqlengths <- as.numeric(chromSizes$seqlengths)

  # Canonical Seqinfo shared by all outputs
  lev <- intersect(GenomeInfoDb::seqlevels(borders), chromSizes$seqnames)
  chromSizes <- chromSizes[match(lev, chromSizes$seqnames), , drop = FALSE]
  si <- GenomeInfoDb::Seqinfo(seqnames = lev, seqlengths = chromSizes$seqlengths)

  borders <- GenomeInfoDb::keepSeqlevels(borders, lev, pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(borders) <- lev
  GenomeInfoDb::seqinfo(borders)   <- si
  borders <- GenomicRanges::sort(borders)

  snap_down <- function(x, bin) as.integer((x %/% bin) * bin)

  out <- GenomicRanges::GRanges(); GenomeInfoDb::seqlevels(out) <- lev; GenomeInfoDb::seqinfo(out) <- si

  for (chr in lev) {
    b <- borders[GenomicRanges::seqnames(borders) == chr]
    if (!length(b)) next

    raw_cut <- switch(boundary_anchor,
      right_edge = IRanges::end(b),
      left_edge  = IRanges::start(b),
      center     = IRanges::start(b) + (IRanges::width(b) %/% 2L)
    )

    cut  <- sort(unique(snap_down(as.integer(raw_cut), binSize)))
    clen <- chromSizes$seqlengths[match(chr, chromSizes$seqnames)]

    # guards (0, clen snapped) + dedupe
    cuts_all <- sort(unique(c(0L, cut, snap_down(clen, binSize))))
    if (length(cuts_all) < 2L) next

    # half-open (cut[i], cut[i+1]] -> start=cut[i]+1, end=cut[i+1]
    starts <- cuts_all[-length(cuts_all)] + 1L
    ends   <- cuts_all[-1L]

    keep <- which(ends >= starts + (min_width - 1L))
    if (!length(keep)) next
    starts <- pmax(1L, starts[keep]); ends <- pmin(clen, ends[keep])

    out <- c(out, GenomicRanges::GRanges(chr, IRanges::IRanges(starts, ends)))
  }

  out <- GenomicRanges::sort(out)
  S4Vectors::mcols(out)$name  <- paste0("TAD_", seq_along(out))
  S4Vectors::mcols(out)$class <- "domain"
  out
}


# -----------------------------------------------------------------------------
# Helper 3: align a GRanges to chromSizes and to a single “natural” UCSC order
#   - drops extra chromosomes (keeps intersection)
#   - enforces same seqlevels order (natural UCSC order derived from gr)
#   - assigns consistent Seqinfo from chromSizes
#   - optional row sort (kept FALSE by default; not required for HicAggR)
# -----------------------------------------------------------------------------
align_to_chromSizes <- function(gr, chromSizes, sort_rows = FALSE) {
  stopifnot(is.data.frame(chromSizes),
            all(c("seqnames","seqlengths") %in% names(chromSizes)))
  chromSizes$seqnames   <- as.character(chromSizes$seqnames)
  chromSizes$seqlengths <- as.numeric(chromSizes$seqlengths)

  # keep only chromosomes present in chromSizes
  keep <- intersect(seqlevels(gr), chromSizes$seqnames)
  gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse")

  # compute natural (UCSC-like) order from current levels
  lev0 <- seqlevels(gr)
  ord  <- GenomeInfoDb::rankSeqlevels(lev0)
  levN <- lev0[order(ord, na.last = TRUE)]

  # enforce same levels order & consistent Seqinfo
  seqlevels(gr) <- levN
  chromSizes    <- chromSizes[match(levN, chromSizes$seqnames), , drop = FALSE]
  si            <- Seqinfo(seqnames = chromSizes$seqnames,
                           seqlengths = chromSizes$seqlengths)
  seqinfo(gr)   <- si

  if (sort_rows) gr <- GenomicRanges::sort(gr, ignore.strand = TRUE)
  gr
}



# Build fixed-size windows shifted relative to a reference point, with strand-aware "upstream/downstream".
# gr        : GRanges
# distance  : bp (scalar o vettore)
# window    : bp (larghezza finale)
# direction : "downstream", "upstream", "none"
# ref       : "center", "TSS", "start", "end" (TSS = start per '+' e end per '-')
# strand    : "gr" (default, usa lo strand di gr) oppure forza "+" / "-" / "*"
# unstranded: gestione delle features con '*' quando strand="gr"
#             - "no_shift" (default): offset=0 per '*'
#             - "as_plus" / "as_minus": tratta '*' come +/-
#             - "drop": escludi '*'
#             - "error": fallisci se esistono '*'
# trim_to_genome: clip ai limiti del cromosoma se presenti le seqlengths
shift_windows_relative <- function(gr, distance, window,
                                   direction = c("downstream","upstream","none"),
                                   ref = c("center","TSS","start","end"),
                                   strand = "gr",
                                   unstranded = c("drop","as_plus","as_minus","no_shift","error"),
                                   trim_to_genome = TRUE) {
  stopifnot(inherits(gr, "GRanges"))
  n <- length(gr); if (n == 0L) return(gr)
  direction  <- match.arg(direction)
  ref        <- match.arg(ref)
  unstranded <- match.arg(unstranded)

  distance <- rep(as.numeric(distance), length.out = n)
  window   <- rep(as.integer(window),   length.out = n)
  if (any(!is.finite(distance))) stop("'distance' must be finite")
  if (any(window <= 0L))         stop("'window' must be > 0")

  # strand per-element
  s <- if (identical(strand, "gr")) as.character(strand(gr)) else rep(as.character(strand), n)

  # Policy sugli unstranded
  if (unstranded == "drop") {
    keep <- s %in% c("+","-")
    gr <- gr[keep]; s <- s[keep]; distance <- distance[keep]; window <- window[keep]
    n <- length(gr); if (n == 0L) return(gr)
  } else if (unstranded == "as_plus") {
    s[s == "*"] <- "+"
  } else if (unstranded == "as_minus") {
    s[s == "*"] <- "-"
  } else if (unstranded == "error" && any(s == "*")) {
    stop("Found unstranded ('*') features; set unstranded= to choose a policy.")
  }
  # "no_shift": '*' resta '*' → offset da direction diventa 0 sotto

  # Punto di riferimento
  ref_pos <- switch(ref,
    center = start(gr) + (width(gr) %/% 2L),
    start  = start(gr),
    end    = end(gr),
    TSS    = ifelse(s == "-", end(gr), start(gr))
  )

  # Segno dell’offset in base a direction e strand
  # downstream: '+' = +dist, '-' = -dist
  # upstream  : '+' = -dist, '-' = +dist
  sign_dir <- switch(direction,
    none       = 0L,
    downstream = ifelse(s == "+",  1L, ifelse(s == "-", -1L, 0L)),
    upstream   = ifelse(s == "+", -1L, ifelse(s == "-",  1L, 0L))
  )
  new_cen <- as.integer(ref_pos + sign_dir * distance)

  # Costruisci la finestra esatta (gestione pari/dispari per mantenere width esatta)
  half_left  <- window %/% 2L
  half_right <- window - half_left - 1L
  new_start  <- pmax(1L, new_cen - half_left)
  new_end    <- new_cen + half_right

  out <- GRanges(seqnames = seqnames(gr),
                 ranges   = IRanges(new_start, new_end),
                 strand   = strand(gr))
  mcols(out) <- mcols(gr)
  if (trim_to_genome && !all(is.na(seqlengths(gr)))) { seqinfo(out) <- seqinfo(gr); out <- trim(out) }
  out
}

#' Build anchor/bait GRanges pairs split by strand combinations
#'
#' This function takes a data.frame produced from `anchor3_l1_strand` and
#' returns three GRanges pairs (anchor/bait) filtered by strand combinations:
#'   1) anchor '+', bait '+'
#'   2) anchor '-', bait '-'
#'   3) opposite (anchor '+', bait '-') OR (anchor '-', bait '+')
#'
#' Anchor strand comes from the LINE family (`line_strand`); bait strand comes
#' from the gene (`gene.strand`). Each GRanges gets a numeric `score` column (set to 1)
#' for compatibility with `IndexFeatures()`.
#'
#' All GRanges are aligned to `chromSizes` (seqlevels order + Seqinfo) to prevent
#' downstream errors in HicAggR/GenomicRanges.
#'
#' @param anchor_df  data.frame with at least:
#'   peak.chr, peak.start, peak.end, peak, line_strand,
#'   gene.chr, gene.start, gene.end, gene.strand, gene
#' @param chromSizes data.frame with columns: seqnames, seqlengths (bp)
#'
#' @return A list with three elements:
#'   - $anch_plus_bait_plus$anchor  (GRanges)
#'   - $anch_plus_bait_plus$bait    (GRanges)
#'   - $anch_minus_bait_minus$anchor
#'   - $anch_minus_bait_minus$bait
#'   - $anch_bait_opposite$anchor
#'   - $anch_bait_opposite$bait
#'
#' @examples
#' pairs_by_strand <- make_strand_pairs(anchor3_l1_strand, chromSizes)
#' pairs_by_strand$anch_plus_bait_plus$anchor
#' pairs_by_strand$anch_plus_bait_plus$bait
make_strand_pairs <- function(anchor_df, chromSizes, use_gene_strand = TRUE) {
  # Required columns
  req_cols <- c("peak.chr","peak.start","peak.end","peak","line_strand",
                "gene.chr","gene.start","gene.end","gene.strand","gene")
  stopifnot(all(req_cols %in% names(anchor_df)))

  # Keep only rows where both anchor (LINE) and gene strands are defined
  df <- subset(anchor_df, line_strand %in% c("+","-") & gene.strand %in% c("+","-"))

  # Split by strand combinations
  if (use_gene_strand) {
    df_pp <- subset(df, line_strand == "+" & gene.strand == "+")  # +/+
    df_mm <- subset(df, line_strand == "-" & gene.strand == "-")  # -/-
    df_pm <- subset(df, line_strand == "+" & gene.strand == "-")  # +/-
    df_mp <- subset(df, line_strand == "-" & gene.strand == "+")  # -/+
    df_op <- rbind(df_pm, df_mp)                                  # opposite pooled (optional)
  } else {
    # If gene strand is ignored: keep only by anchor (LINE) strand
    df_pp <- subset(df, line_strand == "+")
    df_mm <- subset(df, line_strand == "-")
    df_pm <- NULL
    df_mp <- NULL
    df_op <- NULL
  }

  # Robust builder: returns empty GRanges when dfx is NULL/empty
  build_gr_pair <- function(dfx) {
    if (is.null(dfx) || nrow(dfx) == 0L) {
      empty <- GenomicRanges::GRanges()
      return(list(anchor = empty, bait = empty))
    }

    # Anchor: peak.* coords, strand from LINE family
    anchor.gr <- GenomicRanges::GRanges(
      seqnames = dfx$peak.chr,
      ranges   = IRanges::IRanges(start = as.integer(dfx$peak.start),
                                  end   = as.integer(dfx$peak.end)),
      strand   = dfx$line_strand
    )
    S4Vectors::mcols(anchor.gr)$name  <- dfx$peak
    S4Vectors::mcols(anchor.gr)$score <- 1

    # Bait: gene.* coords, strand from gene
    bait.gr <- GenomicRanges::GRanges(
      seqnames = dfx$gene.chr,
      ranges   = IRanges::IRanges(start = as.integer(dfx$gene.start),
                                  end   = as.integer(dfx$gene.end)),
      strand   = dfx$gene.strand
    )
    S4Vectors::mcols(bait.gr)$name  <- dfx$gene
    S4Vectors::mcols(bait.gr)$score <- 1

    # Align both to chromSizes order/Seqinfo
    anchor.gr <- align_to_chromSizes(anchor.gr, chromSizes, sort_rows = FALSE)
    bait.gr   <- align_to_chromSizes(bait.gr,   chromSizes, sort_rows = FALSE)

    list(anchor = anchor.gr, bait = bait.gr)
  }

  # Return all combinations (plus the pooled opposite if use_gene_strand=TRUE)
  out <- list(
    anch_plus_bait_plus    = build_gr_pair(df_pp),  # +/+
    anch_minus_bait_minus  = build_gr_pair(df_mm),  # -/-
    anch_plus_bait_minus   = build_gr_pair(df_pm),  # +/-
    anch_minus_bait_plus   = build_gr_pair(df_mp),  # -/+
    anch_bait_opposite     = build_gr_pair(df_op)   # pooled opposite (NULL/empty if not applicable)
  )

  out
}


#-----------------------------------------------------------------------------
#
#
#
#-----------------------------------------------------------------------------
#' APA_analysis
#'
#' End-to-end APA workflow for HicAggR:
#'  1) Align input features (anchor/bait) to a single chromosome order/Seqinfo
#'  2) Build a genomic constraint from an insulation bigWig (borders → TAD domains)
#'  3) Index features into Hi-C bins (anchor & bait) under the constraint
#'  4) Build putative pairs (A↔B) inside the same TAD
#'  5) Prepare Hi-C (balance if needed, then convert to O/E if needed)
#'  6) Extract submatrices around pairs, orient them, aggregate (APA)
#'  7) Return a ggplot heatmap (ggAPA)
#'
#' @param resolution   Numeric. Your intended Hi-C resolution (kept for bookkeeping).
#' @param anchor_grange GRanges. Anchor features (e.g., peaks). Must/gets a numeric 'score' column.
#' @param bait_grange   GRanges. Bait features (e.g., genes). Must/gets a numeric 'score' column.
#' @param chromSizes    data.frame with columns 'seqnames' and 'seqlengths' (bp).
#' @param path_constraint_bigwig Character. Path to insulation bigWig used to derive borders/TADs.
#' @param binSize       Integer. Hi-C bin size used for indexing and submatrices (default 5000).
#' @param hic_in        List of ContactMatrix (from ImportHiC). Can be raw, balanced, or O/E.
#' @param orientate
#' 
#' @return A ggplot object from ggAPA (or NULL if no pairs or empty extraction).
APA_analysis <- function(
  resolution,
  anchor_grange,
  bait_grange,
  chromSizes,
  path_constraint_bigwig,
  binSize      = 5000,
  hic_in,
  orientate    = TRUE,
  search_pairs = FALSE,              # FALSE = abbina riga-a-riga; TRUE = SearchPairs (cartesiano entro TAD)
  oe_method    = "mean_no_zero"      # NULL = salta O/E; altrimenti passa a OverExpectedHiC()
){
  ## ---------------------------
  ## 0) Inputs & chromSizes
  ## ---------------------------
  stopifnot(is.data.frame(chromSizes),
            all(c("seqnames","seqlengths") %in% names(chromSizes)))
  chromSizes$seqnames   <- as.character(chromSizes$seqnames)
  chromSizes$seqlengths <- as.numeric(chromSizes$seqlengths)

  ## ---------------------------
  ## 1) Align input GRanges
  ## ---------------------------
  peaks.gr <- align_to_chromSizes(anchor_grange, chromSizes, sort_rows = FALSE)
  genes.gr <- align_to_chromSizes(bait_grange,   chromSizes, sort_rows = FALSE)

peaks.gr$pair_id <- seq_along(peaks.gr)
genes.gr$pair_id <- seq_along(genes.gr)

  if (!"score" %in% colnames(S4Vectors::mcols(peaks.gr))) S4Vectors::mcols(peaks.gr)$score <- 1
  if (!"score" %in% colnames(S4Vectors::mcols(genes.gr))) S4Vectors::mcols(genes.gr)$score <- 1
  S4Vectors::mcols(peaks.gr)$score <- as.numeric(S4Vectors::mcols(peaks.gr)$score)
  S4Vectors::mcols(genes.gr)$score <- as.numeric(S4Vectors::mcols(genes.gr)$score)

  # keep chromSizes in the same seqlevel order used by features
  chromSizes <- chromSizes[match(GenomeInfoDb::seqlevels(peaks.gr), chromSizes$seqnames), , drop = FALSE]

  ## ---------------------------
  ## 2) Constraint: borders -> TADs
  ## ---------------------------
  borders <- hicaggr_borders_from_bigwig(
    bw_path = path_constraint_bigwig,
    k = 5, min_separation = 50000
  )
  TADs_domains <- borders_to_domains(borders, chromSizes, binSize = binSize)
  constraint   <- align_to_chromSizes(TADs_domains, chromSizes)

  ## ---------------------------
  ## 3) Index features to Hi-C bins
  ## ---------------------------
  anc <- IndexFeatures(
    gRangeList        = list(Anchor = peaks.gr),
    genomicConstraint = constraint,
    chromSizes        = chromSizes,
    binSize           = binSize,
    metadataColName   = "score",
    method            = "sum"
  )
  bait <- IndexFeatures(
    gRangeList        = list(Bait = genes.gr),
    genomicConstraint = constraint,
    chromSizes        = chromSizes,
    binSize           = binSize,
    metadataColName   = "score",
    method            = "sum"
  )
 
  all_pairs <- SearchPairs(anc, bait)

if (isTRUE(search_pairs)) {

  message("searched pairs (cartesian product within TADs)")
  int_pairs <- all_pairs

} else {

  ## NIENTE row-wise dopo IndexFeatures: usa direttamente le coppie originali 1:1
   message("using original curated pairs (row-wise GRanges)")
  # nmin <- min(length(peaks.gr), length(genes.gr))
  # int_pairs <- InteractionSet::GInteractions(
  #   peaks.gr[seq_len(nmin)],
  #   genes.gr[seq_len(nmin)]
  # )
  ## 5) Filtro all_pairs per tenere solo le tue coppie 1:1 (peak[i]–gene[i])

  A_bins <- InteractionSet::anchors(all_pairs, type = "first")
  B_bins <- InteractionSet::anchors(all_pairs, type = "second")

# overlaps anchor-bin <-> peaks.gr
ha <- GenomicRanges::findOverlaps(A_bins, peaks.gr, ignore.strand = TRUE)
# overlaps bait-bin   <-> genes.gr
hb <- GenomicRanges::findOverlaps(B_bins, genes.gr, ignore.strand = TRUE)

# per ogni interazione j in all_pairs, quali pair_id tocca sul lato anchor e sul lato bait?
anchor_ids_by_pair <- split(peaks.gr$pair_id[ S4Vectors::subjectHits(ha) ],
                            S4Vectors::queryHits(ha))
bait_ids_by_pair   <- split(genes.gr$pair_id[ S4Vectors::subjectHits(hb) ],
                            S4Vectors::queryHits(hb))

keep <- logical(length(all_pairs))
for (j in seq_along(all_pairs)) {
  idsA <- anchor_ids_by_pair[[as.character(j)]]
  idsB <- bait_ids_by_pair[[as.character(j)]]
  if (length(idsA) && length(idsB) && length(intersect(idsA, idsB)) > 0) {
    keep[j] <- TRUE
  }
}

int_pairs <- all_pairs[keep]

if (!length(int_pairs)) {
  warning("Nessuna interazione di SearchPairs corrisponde alle coppie curate.")
  return(invisible(NULL))
}

}
  # opzionale: rimuovere auto-collisioni (stesso bin)
  # bait <- bait[ match(bait$bin, anc$bin, nomatch = 0) == 0 ]

  ## ---------------------------
  ## 4) Build pairs
  ## ---------------------------
  # if (isTRUE(search_pairs)) {
  #   message("searching pairs (cartesian product within TADs)")
  #   int_pairs <- SearchPairs(anc, bait)
  # } else {
  #   message("direct matching of pairs (row-wise)")
  #   nmin <- min(length(anc), length(bait))
  #   int_pairs <- InteractionSet::GInteractions(anc[seq_len(nmin)], bait[seq_len(nmin)])
  # }
  # if (length(int_pairs) == 0L) {
  #   warning("No pairs produced; no APA possible.")
  #   return(invisible(NULL))
  # }

  ## ---------------------------
  ## 5) Prepare Hi-C (balance -> O/E if needed)
  ## ---------------------------
  kind_list   <- attributes(hic_in)$mtx
  kind_matrix <- try(S4Vectors::metadata(hic_in[[1]])$mtx, silent = TRUE)
  H <- hic_in

  if (!identical(kind_list, "o/e") && !identical(kind_matrix, "o/e")) {
    if (!identical(kind_list, "norm") && !identical(kind_matrix, "norm")) {
      bal_method <- "ICE"                       # "ICE", "VC", "VC_SQRT"
      H_1<- BalanceHiC(H, method = bal_method)
      H <- H_1
      message(paste0(bal_method, " normalization applied"))
    }
    if (!is.null(oe_method)) {
      message(paste0(oe_method, " distance decay correction applied"))
      H_1 <- OverExpectedHiC(H, method = oe_method, verbose = TRUE)
      H <- H_1
    } else {
      message("distance decay correction NOT applied")
    }
  }

  ## ---------------------------
  ## 6) Extract & orient submatrices
  ## ---------------------------
  M <- ExtractSubmatrix(
    genomicFeature = int_pairs,
    hicLst         = H,
    referencePoint = "pf",
    matriceDim     = 51
  )
  M_or <- PrepareMtxList(M, orientate = orientate, rm0 = FALSE)

## ---- Identify which pairs yielded a valid matrix (aligned to M_or) ----
# A pair is "valid" if its submatrix has at least one finite entry
is_valid  <- vapply(M_or, function(x) is.matrix(x) && any(is.finite(x)), logical(1))
valid_idx <- which(is_valid)

# Extract anchor/bait GRanges for ALL proposed pairs via regions()+anchors()
idx  <- InteractionSet::anchors(int_pairs, id = TRUE)   # list(first=..., second=...)
regs <- InteractionSet::regions(int_pairs)               # GRanges with all bins
A <- regs[idx$first]
B <- regs[idx$second]

# Helper: read a metadata column as character if present
get_mcol_chr <- function(gr, nm) {
  if (nm %in% colnames(S4Vectors::mcols(gr))) as.character(S4Vectors::mcols(gr)[[nm]])
  else rep(NA_character_, length(gr))
}

# Build a BEDPE-like table for ALL pairs
pairs_df_all <- data.frame(
  anchor_chr   = as.character(GenomicRanges::seqnames(A)),
  anchor_start = as.integer(GenomicRanges::start(A)),
  anchor_end   = as.integer(GenomicRanges::end(A)),
  anchor_name  = get_mcol_chr(A, "name"),
  anchor_bin   = get_mcol_chr(A, "bin"),
  bait_chr     = as.character(GenomicRanges::seqnames(B)),
  bait_start   = as.integer(GenomicRanges::start(B)),
  bait_end     = as.integer(GenomicRanges::end(B)),
  bait_name    = get_mcol_chr(B, "name"),
  bait_bin     = get_mcol_chr(B, "bin"),
  stringsAsFactors = FALSE
)

# Optional: distance in bins if labels look like "chr1:12345"
pairs_df_all$bin_distance <- suppressWarnings({
  if (!all(is.na(pairs_df_all$anchor_bin)) && !all(is.na(pairs_df_all$bait_bin))) {
    abin <- suppressWarnings(as.integer(sub(".*:", "", pairs_df_all$anchor_bin)))
    bbin <- suppressWarnings(as.integer(sub(".*:", "", pairs_df_all$bait_bin)))
    ifelse(is.na(abin) | is.na(bbin), NA_integer_, abs(bbin - abin))
  } else NA_integer_
})

# Human-readable id for ALL pairs; name submatrices with a length-matched vector
pair_id <- sprintf("%s:%s|%s:%s",
  pairs_df_all$anchor_chr,
  ifelse(is.na(pairs_df_all$anchor_bin),
         as.character(pairs_df_all$anchor_start),
         pairs_df_all$anchor_bin),
  pairs_df_all$bait_chr,
  ifelse(is.na(pairs_df_all$bait_bin),
         as.character(pairs_df_all$bait_start),
         pairs_df_all$bait_bin)
)

# Ensure names(M_or) length matches; truncate and make unique to avoid duplicates
names(M_or) <- make.unique(pair_id[seq_along(M_or)])

# Keep only rows corresponding to submatrices we actually have AND are valid
valid_pair_ids <- names(M_or)[valid_idx]
valid_rows     <- match(valid_pair_ids, pair_id)  # indices in pairs_df_all
valid_pairs    <- pairs_df_all[valid_rows, , drop = FALSE]
valid_pairs$pair_id <- valid_pair_ids

## Diagnostics
n_mat   <- length(M_or)
nz_frac <- if (length(valid_idx)) {
  mean(vapply(M_or[valid_idx], function(x) sum(x != 0, na.rm = TRUE) / length(x), numeric(1)))
} else NA_real_
message(sprintf("Diagnostics — pairs: %d | matrices: %d | valid: %d | mean non-zero frac (valid): %s",
                length(int_pairs), n_mat, length(valid_idx),
                ifelse(is.na(nz_frac), "NA", sprintf('%.3f', nz_frac))))

if (!length(valid_idx)) {
  warning("No valid submatrices extracted (all NA/empty).")
  return(invisible(NULL))
}

## ---- Aggregate only valid submatrices for APA ----
APA <- Aggregation(ctrlMatrices = M_or[valid_idx], aggFun = "mean", scaleCorrection = FALSE)

attr(APA, "referencePoint") <- "pf"

## ---- Return everything you need (includes 'valid_pairs') ----
return(list(
  APA         = APA,         # aggregated APA matrix
  M_or        = M_or,        # (oriented) submatrices (all)
  pairs       = int_pairs,   # all proposed pairs (GInteractions)
  valid_idx   = valid_idx,   # indices of valid submatrices
  valid_pairs = valid_pairs,  # table of pairs that entered the plot
  H_used      = H            # <<--- Hi-C usata (normalizzata / O/E)
))
}

