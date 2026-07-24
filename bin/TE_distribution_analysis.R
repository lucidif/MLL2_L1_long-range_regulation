#!/usr/bin/env Rscript
# =====================================================================
# TE distribution analysis for KMT2B distal CpG-negative peaks
#
# This script uses RepeatMasker annotations and repository helper functions
# to compute:
#   1) TE macro-class composition of the nearest TE within 500 bp of
#      distal CpG-negative peaks
#   2) L1 subfamily percentages among those peak-proximal L1s,
#      compared with genome-wide L1 proportions
#
# This version avoids commands that were not used in the temporary scripts
# (e.g. full_join, dplyr pipelines, GenomicRanges, etc.).
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

# ------------------------------ USER SETTINGS ---------------------------
SRC_DIST   <- "pipelines/downstream_multiomic/bin/distance_functions.R"
RMSK_FILE  <- "in/ucsc/rmsk.txt"
PEAK_BED   <- "outs/CHiP_postprocessing_line1_dist/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
OUTPUT_DIR <- "outs/TE_distribution"
WINDOW_BP  <- 500

# ------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

macro_class_from_repeat <- function(element, type) {
  if (!is.na(type) && grepl("ERV", type, ignore.case = TRUE)) {
    return("ERV")
  }
  if (!is.na(element) && grepl("LINE", element, ignore.case = TRUE)) {
    return("LINE")
  }
  if (!is.na(element) && grepl("SINE", element, ignore.case = TRUE)) {
    return("SINE")
  }
  if (!is.na(element) && grepl("LTR", element, ignore.case = TRUE)) {
    return("LTR")
  }
  "OTHER"
}

plot_donut <- function(df, label_col, value_col, title, outfile) {
  df$x <- 2
  df$label <- paste0(df[[label_col]], "\n", df[[value_col]], " (", df$pct, "% )")
  p <- ggplot(df, aes_string(x = "x", y = value_col, fill = label_col)) +
    geom_col(color = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void() +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.2) +
    labs(title = title, fill = label_col)
  ggsave(outfile, plot = p, width = 6, height = 6, device = "pdf")
}

# ------------------------------------------------------------------------
# Inputs and validation
# ------------------------------------------------------------------------
if (!file.exists(SRC_DIST)) {
  stop(sprintf("Missing helper script: %s", SRC_DIST))
}
source(SRC_DIST)

if (!file.exists(RMSK_FILE)) {
  stop(sprintf("Missing RepeatMasker file: %s", RMSK_FILE))
}
if (!file.exists(PEAK_BED)) {
  stop(sprintf("Missing peak BED file: %s", PEAK_BED))
}

ensure_dir(OUTPUT_DIR)

rmsk <- read.table(RMSK_FILE, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE)
if (ncol(rmsk) < 13) {
  stop("RepeatMasker file must have at least 13 columns.")
}

peaks <- read.table(PEAK_BED, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(peaks) < 3) {
  stop("Peak BED file must have at least 3 columns.")
}

# ------------------------------------------------------------------------
# Nearest TE calculation
# ------------------------------------------------------------------------
distances <- extractdistances(
  peaks.coords     = peaks,
  repeat_mask_line = rmsk,
  addfamily        = TRUE,
  addpeaks         = TRUE,
  addlines         = TRUE,
  add5prime        = FALSE
)

distance_values <- as.numeric(as.character(distances[, 1]))
selected <- which(distance_values <= WINDOW_BP)
if (length(selected) == 0) {
  stop("No peak is within the specified window of any TE element.")
}

nearby <- distances[selected, , drop = FALSE]
nearest_df <- data.frame(
  distance = distance_values[selected],
  peaks_chr = nearby[, "peaks_chr"],
  peaks_start = as.numeric(as.character(nearby[, "peaks_start"])),
  peaks_end = as.numeric(as.character(nearby[, "peaks_end"])),
  family = nearby[, "family"],
  line_chr = nearby[, "line_chr"],
  line_start = as.numeric(as.character(nearby[, "line_start"])),
  line_end = as.numeric(as.character(nearby[, "line_end"])),
  line_strand = nearby[, "line_strand"],
  stringsAsFactors = FALSE
)

# retain repeat metadata for classification
repeat_info <- rmsk[, c(6, 7, 8, 10, 11, 12, 13)]
colnames(repeat_info) <- c("line_chr", "line_start", "line_end", "line_strand", "family", "element", "type")
repeat_info$line_start <- as.numeric(as.character(repeat_info$line_start))
repeat_info$line_end <- as.numeric(as.character(repeat_info$line_end))
repeat_info$line_strand <- as.character(repeat_info$line_strand)
repeat_info$family <- as.character(repeat_info$family)

merged <- merge(
  nearest_df,
  repeat_info,
  by = c("line_chr", "line_start", "line_end", "line_strand", "family"),
  all.x = TRUE,
  sort = FALSE
)
if (nrow(merged) != nrow(nearest_df)) {
  warning("Row count changed after merging in RepeatMasker metadata.")
}

merged$macro_class <- mapply(
  macro_class_from_repeat,
  merged$element,
  merged$type,
  USE.NAMES = FALSE
)

write.table(
  merged,
  file = file.path(OUTPUT_DIR, "nearest_TE_within_window.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ------------------------------------------------------------------------
# Macro-class composition
# ------------------------------------------------------------------------
macro_tab <- as.data.frame(table(merged$macro_class), stringsAsFactors = FALSE)
colnames(macro_tab) <- c("macro_class", "count")
macro_tab$pct <- round(macro_tab$count / sum(macro_tab$count) * 100, 1)

write.table(
  macro_tab,
  file = file.path(OUTPUT_DIR, "TE_macroclass_nearest_500bp.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

plot_donut(
  macro_tab,
  label_col = "macro_class",
  value_col = "count",
  title = sprintf("TE macro-class composition within %dbp of distal CpG-negative peaks", WINDOW_BP),
  outfile = file.path(OUTPUT_DIR, "TE_macroclass_nearest_500bp_donut.pdf")
)

# ------------------------------------------------------------------------
# L1 subfamily analysis
# ------------------------------------------------------------------------
l1_ref <- rmsk[grep("^L1", rmsk[[11]], ignore.case = TRUE), ]
l1_genome <- as.data.frame(table(l1_ref[[11]]), stringsAsFactors = FALSE)
colnames(l1_genome) <- c("subfamily", "genome_count")
l1_genome$genome_pct <- round(l1_genome$genome_count / sum(l1_genome$genome_count) * 100, 1)

l1_peaks <- merged[grep("^L1", merged$family, ignore.case = TRUE), ]
l1_peak_tab <- as.data.frame(table(l1_peaks$family), stringsAsFactors = FALSE)
colnames(l1_peak_tab) <- c("subfamily", "peak_count")
l1_peak_tab$peak_pct <- round(l1_peak_tab$peak_count / sum(l1_peak_tab$peak_count) * 100, 1)

l1_summary <- merge(l1_genome, l1_peak_tab, by = "subfamily", all = TRUE)
for (col in c("genome_count", "peak_count", "genome_pct", "peak_pct")) {
  if (col %in% colnames(l1_summary)) {
    l1_summary[is.na(l1_summary[[col]]), col] <- 0
  }
}

write.table(
  l1_summary,
  file = file.path(OUTPUT_DIR, "L1_subfamily_genome_vs_peaks.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

plot_top_donut <- function(df, label_col, value_col, title, outfile, top_n = 8) {
  df <- df[order(-df[[value_col]]), ]
  if (nrow(df) > top_n) {
    others <- sum(df[[value_col]][(top_n + 1):nrow(df)])
    df <- df[1:top_n, , drop = FALSE]
    df <- rbind(df, data.frame(subfamily = "Others", count = others, stringsAsFactors = FALSE))
    if (value_col == "peak_count") {
      df$pct <- round(df$count / sum(df$count) * 100, 1)
    } else if (value_col == "genome_count") {
      df$pct <- round(df$count / sum(df$count) * 100, 1)
    }
  } else {
    df$pct <- round(df[[value_col]] / sum(df[[value_col]]) * 100, 1)
  }
  plot_donut(df, label_col, value_col, title, outfile)
}

plot_top_donut(
  transform(l1_peak_tab, count = peak_count),
  label_col = "subfamily",
  value_col = "count",
  title = sprintf("L1 subfamily composition within %dbp of distal CpG-negative peaks", WINDOW_BP),
  outfile = file.path(OUTPUT_DIR, "L1_subfamily_peaks_donut.pdf")
)

plot_top_donut(
  transform(l1_genome, count = genome_count),
  label_col = "subfamily",
  value_col = "count",
  title = "Genome-wide L1 subfamily composition",
  outfile = file.path(OUTPUT_DIR, "L1_subfamily_genome_donut.pdf")
)

message("TE distribution analysis complete. Results written to: ", normalizePath(OUTPUT_DIR))
