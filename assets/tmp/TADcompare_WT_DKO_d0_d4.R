# =============================================================================
# TAD Analysis Pipeline: Calling, Boundary Extraction, and Pileup Preparation
# =============================================================================
#
# METHODS SUMMARY (for manuscript):
#
# Topologically associating domains (TADs) were identified from Micro-C contact
# matrices using the SpectralTAD algorithm (Cresswell & Dozmorov, 2020,
# Bioinformatics), which detects domain boundaries via spectral decomposition
# of the contact matrix Laplacian. Contact matrices were generated from merged
# pairtools output using cooler cload pairix (Abdennur & Mirny, 2020) at 50 kb
# resolution against the mm10 reference genome. Prior to TAD calling, matrices
# were iteratively corrected (ICE balancing) using cooler balance with default
# parameters. chrM and chrY were excluded from all analyses.
#
# TAD calling was performed independently for each sample using SpectralTAD
# with default parameters at 50 kb resolution, retaining only Level 1 TADs
# (primary hierarchy). TADs smaller than 4 bins (200 kb) were discarded as
# likely representing noise at this resolution.
#
# TAD boundaries were defined as the genomic positions where two adjacent,
# contiguous TADs meet. For each pair of contiguous TADs (i.e., where the end
# coordinate of TAD_i equals the start coordinate of TAD_i+1), a 2 bp anchor
# interval was placed at the center of the shared boundary bin to provide a
# resolution-independent coordinate for downstream pileup analysis. Boundaries
# at gaps between non-contiguous TADs were excluded.
#
# Aggregate TAD boundary analysis (pileup) was performed using coolpup.py
# v1.1.0 (Flyamer et al., 2020) on ICE-balanced 50 kb contact matrices.
# WT Day 0 and WT Day 4 boundary coordinates were used as anchors for the
# corresponding WT and Mll1/Mll2 double knockout (DKO) pileups at each time point, ensuring
# that WT and KO pileups are compared at identical genomic positions.
# Pileups were computed with a flanking window of 500 kb on each side of the
# boundary (--flank 500000), 10 local background shifts (--nshifts 10),
# and excluding the 2 nearest off-diagonal bins (--ignore_diags 2).
# Observed/expected normalization was applied using pre-computed cis expected
# values from cooltools expected-cis. Pileup matrices were visualized using
# plotpup.py with a color scale centered at O/E = 1.
#
# =============================================================================

library(TADCompare)
library(HiCcompare)
library(SpectralTAD)
library(ggplot2)
library(dplyr)

# =============================================================================
# Global parameters
# =============================================================================

# Micro-C contact matrix resolution (base pairs)
RESOLUTION <- 50000

# Minimum TAD size: TADs smaller than 4 bins (200 kb at 50 kb resolution)
# are discarded as unreliable calls at this resolution
MIN_BINS <- 4
MIN_SIZE <- MIN_BINS * RESOLUTION  # 200,000 bp

# Chromosomes excluded from all analyses
EXCLUDE_CHR <- c("chrM", "chrY")

# Output directories
OUTDIR_TADCOMPARE  <- "outs/2024_10_Lara_microC_downstream/out/TADCompare"
OUTDIR_PILEUP      <- "outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup"

# =============================================================================
# PART 1: TADCompare -- Differential TAD analysis (WT vs KO)
# =============================================================================
# TADCompare (Cresswell & Dozmorov, 2021) was used to classify TAD boundaries
# as differential or non-differential between WT and Mll1/Mll2 double knockout
# (DKO) at each time point.
# Input matrices were loaded from ICE-balanced cooler dumps (tab-separated
# sparse upper-triangular format produced by cooler dump --join).
# =============================================================================

# --- Load contact matrices ---

WT_day0 <- read.table("outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_WT_day0.Dd.txt")
sparse_WT_day0 <- HiCcompare::cooler2sparse(WT_day0)
rm(WT_day0); gc()

KO_day0 <- read.table("outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_KO_day0.Dd.txt")
sparse_KO_day0 <- HiCcompare::cooler2sparse(KO_day0)
rm(KO_day0); gc()

sparse_WT_day0 <- sparse_WT_day0[!names(sparse_WT_day0) %in% EXCLUDE_CHR]
sparse_KO_day0 <- sparse_KO_day0[!names(sparse_KO_day0) %in% EXCLUDE_CHR]

WT_day4 <- read.table("outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_WT_day4.Dd.txt")
sparse_WT_day4 <- HiCcompare::cooler2sparse(WT_day4)
rm(WT_day4); gc()

KO_day4 <- read.table("outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_KO_day4.Dd.txt")
sparse_KO_day4 <- HiCcompare::cooler2sparse(KO_day4)
rm(KO_day4); gc()

sparse_WT_day4 <- sparse_WT_day4[!names(sparse_WT_day4) %in% EXCLUDE_CHR]
sparse_KO_day4 <- sparse_KO_day4[!names(sparse_KO_day4) %in% EXCLUDE_CHR]

# Retain only chromosomes present in both WT and KO
common_chrs <- intersect(names(sparse_WT_day0), names(sparse_KO_day0))

# --- Run TADCompare: WT vs KO at Day 0 ---
# TADCompare classifies each boundary as: Non-Differential, Strength Change,
# Shifted, Split, Merge, or Complex, based on boundary score differences
# between conditions.

for (chr in common_chrs) {
  message("TADCompare WT vs KO - Day 0 - ", chr)
  res <- tryCatch({
    TADCompare(sparse_WT_day0[[chr]], sparse_KO_day0[[chr]], resolution = RESOLUTION)
  }, error = function(e) {
    message("ERROR on ", chr, ": ", e$message)
    return(NULL)
  })
  if (!is.null(res)) assign(paste0(chr, "_result_d0"), res, envir = .GlobalEnv)
}

# --- Run TADCompare: WT vs KO at Day 4 ---

for (chr in common_chrs) {
  message("TADCompare WT vs KO - Day 4 - ", chr)
  res <- tryCatch({
    TADCompare(sparse_WT_day4[[chr]], sparse_KO_day4[[chr]], resolution = RESOLUTION)
  }, error = function(e) {
    message("ERROR on ", chr, ": ", e$message)
    return(NULL)
  })
  if (!is.null(res)) assign(paste0(chr, "_result_d4"), res, envir = .GlobalEnv)
}

# --- Aggregate boundary type counts across chromosomes ---

count_data_list_d0 <- lapply(common_chrs, function(chr) {
  obj_name <- paste0(chr, "_result_d0")
  if (exists(obj_name)) {
    plot_data <- get(obj_name)$Count_Plot$data
    if (!is.null(plot_data)) { plot_data$chromosome <- chr; return(plot_data) }
  }
  return(NULL)
})

count_data_list_d4 <- lapply(common_chrs, function(chr) {
  obj_name <- paste0(chr, "_result_d4")
  if (exists(obj_name)) {
    plot_data <- get(obj_name)$Count_Plot$data
    if (!is.null(plot_data)) { plot_data$chromosome <- chr; return(plot_data) }
  }
  return(NULL)
})

count_data_all_d0 <- do.call(rbind, Filter(Negate(is.null), count_data_list_d0))
count_data_all_d4 <- do.call(rbind, Filter(Negate(is.null), count_data_list_d4))

# Summarise counts per boundary type
type <- c("Complex", "Merge", "Non-Differential", "Shifted", "Split", "Strength Change", NA)

summarise_counts <- function(count_data, types) {
  sapply(types, function(t) {
    if (is.na(t)) sum(count_data[is.na(count_data$Type), "Count"])
    else sum(count_data[count_data$Type == t, "Count"], na.rm = TRUE)
  })
}

d0_count <- cbind(type, count = summarise_counts(count_data_all_d0, type), time = "d0")
d4_count <- cbind(type, count = summarise_counts(count_data_all_d4, type), time = "d4")
count     <- data.frame(rbind(d0_count, d4_count))

write.table(count,
            file.path(OUTDIR_TADCOMPARE, "TADstats.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Stacked bar plot of boundary type distribution across time points
p <- ggplot(data = count, aes(x = time, y = as.numeric(count), fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "TAD boundary classification (TADCompare, WT vs KO)",
       x     = "Time point",
       y     = "Number of boundaries",
       fill  = "Boundary type")

ggsave(file.path(OUTDIR_TADCOMPARE, "TAD_barplot.png"),
       plot = p, width = 8, height = 6, dpi = 300)

# =============================================================================
# PART 2: TAD calling on WT samples (SpectralTAD)
# =============================================================================
# TADs were called independently on WT Day 0 and WT Day 4 contact matrices
# using SpectralTAD with default parameters. Only Level 1 TADs (primary
# hierarchy) were retained. TADs smaller than 4 bins (200 kb) were discarded.
# WT-derived boundaries were subsequently used as reference anchors for
# aggregate pileup analysis in both WT and DKO, ensuring that differences
# between conditions are assessed at identical genomic coordinates.
# =============================================================================

call_TADs <- function(sparse_list, resolution, min_size, label) {
  # Run SpectralTAD chromosome by chromosome and return a genome-wide BED
  # data frame of Level 1 TADs passing the minimum size filter.
  tad_list <- lapply(names(sparse_list), function(chr) {
    message("SpectralTAD on ", label, " - ", chr)
    res <- tryCatch({
      SpectralTAD(sparse_list[[chr]], chr = chr, resolution = resolution)
    }, error = function(e) {
      message("ERROR on ", chr, ": ", e$message); return(NULL)
    })
    if (!is.null(res)) {
      df <- res$Level_1
      colnames(df) <- c("chr", "start", "end")
      # Discard TADs smaller than the minimum size threshold (4 bins = 200 kb)
      df <- df[(df$end - df$start) >= min_size, ]
      return(df)
    }
    return(NULL)
  })
  do.call(rbind, Filter(Negate(is.null), tad_list))
}

# WT Day 0 -- sparse matrix already in memory from Part 1
sparse_WT_day0 <- sparse_WT_day0[!names(sparse_WT_day0) %in% EXCLUDE_CHR]

tad_WT_d0 <- call_TADs(sparse_WT_day0, RESOLUTION, MIN_SIZE, "WT_D0")

write.table(tad_WT_d0,
            file.path(OUTDIR_TADCOMPARE, "TAD_WT_d0_filtered.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

message("WT D0: ", nrow(tad_WT_d0), " TADs retained after size filtering")
rm(sparse_WT_day0); gc()

# WT Day 4
sparse_WT_day4 <- sparse_WT_day4[!names(sparse_WT_day4) %in% EXCLUDE_CHR]

tad_WT_d4 <- call_TADs(sparse_WT_day4, RESOLUTION, MIN_SIZE, "WT_D4")

write.table(tad_WT_d4,
            file.path(OUTDIR_TADCOMPARE, "TAD_WT_d4_filtered.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

message("WT D4: ", nrow(tad_WT_d4), " TADs retained after size filtering")
rm(sparse_WT_day4); gc()

# --- TAD size summary statistics (for Results section) ---

summarise_TADs <- function(tad_df, label) {
  sizes <- tad_df$end - tad_df$start
  message("\n--- ", label, " ---")
  message("Total TADs : ", nrow(tad_df))
  message("Median size: ", median(sizes) / 1e6, " Mb")
  message("Mean size  : ", round(mean(sizes) / 1e6, 3), " Mb")
  message("Min size   : ", min(sizes) / 1e3, " kb")
  message("Max size   : ", max(sizes) / 1e3, " kb")
}

summarise_TADs(tad_WT_d0, "WT Day 0")
summarise_TADs(tad_WT_d4, "WT Day 4")

# =============================================================================
# PART 3: TAD boundary extraction
# =============================================================================
# Boundaries are defined as the contact point between two adjacent, contiguous
# TADs (i.e., end of TAD_i == start of TAD_i+1). Gaps between non-contiguous
# TADs are excluded. For each boundary, a 2 bp anchor interval is placed at
# the center of the shared boundary bin (bin_start + resolution/2 ± 1 bp).
# This resolution-independent coordinate ensures that coolpup.py centers the
# flanking window precisely on the boundary regardless of bin size.
# =============================================================================

extract_boundaries <- function(tad_df, resolution, label) {
  boundaries <- tad_df %>%
    group_by(chr) %>%
    arrange(start, .by_group = TRUE) %>%
    # Retain only positions where adjacent TADs are strictly contiguous
    mutate(is_boundary = start == lag(end)) %>%
    filter(is_boundary == TRUE) %>%
    transmute(
      chr,
      # 2 bp interval centered on the midpoint of the boundary bin
      start = start + (resolution / 2) - 1,
      end   = start + 2
    ) %>%
    ungroup()
  message(label, ": ", nrow(boundaries), " boundaries extracted")
  return(boundaries)
}

colnames(tad_WT_d0) <- c("chr", "start", "end", "level")
colnames(tad_WT_d4) <- c("chr", "start", "end", "level")

boundaries_WT_d0 <- extract_boundaries(tad_WT_d0, RESOLUTION, "WT_D0")
boundaries_WT_d4 <- extract_boundaries(tad_WT_d4, RESOLUTION, "WT_D4")

write.table(boundaries_WT_d0,
            file.path(OUTDIR_TADCOMPARE, "boundaries_WT_d0.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(boundaries_WT_d4,
            file.path(OUTDIR_TADCOMPARE, "boundaries_WT_d4.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# =============================================================================
# PART 4: Pileup -- shell commands (run outside R via Docker)
# =============================================================================
# Aggregate TAD boundary pileup plots were generated using coolpup.py v1.1.0
# on ICE-balanced 50 kb .cool files. For each time point, WT-derived boundary
# coordinates were used as anchors for both WT and DKO pileups,
# enabling direct comparison at identical genomic positions.
# Observed/expected normalization was applied using cis expected values
# pre-computed with cooltools expected-cis.
#
# The following commands were executed inside the Docker container
# quay.io/biocontainers/coolpuppy:1.1.0--pyh086e186_0:
#
# # Compute cis expected values (required for O/E normalization)
# cooltools expected-cis \
#     outs/2024_10_Lara_microC_downstream/out/cooler_cload_balance/balanced_50k_aLp_WT_day0.Dd.cool \
#     --nproc 4 \
#     -o outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup/expected_WT_d0.tsv
#
# # Pileup centered on WT Day 0 boundaries (applied to both WT and KO)
# coolpup.py \
#     <sample>.cool \
#     outs/2024_10_Lara_microC_downstream/out/TADCompare/boundaries_WT_d0.bed \
#     --features_format bed \
#     --flank 500000 \
#     --nshifts 10 \
#     --ignore_diags 2 \
#     --expected outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup/expected_<sample>.tsv \
#     -o outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup/pileup_boundaries_<sample>.clpy
#
# # Visualization
# plotpup.py \
#     --plot_ticks \
#     --vmin 0.7 --vmax 1.3 \
#     --center 1 \
#     --input_pups outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup/pileup_boundaries_<sample>.clpy \
#     --output outs/2024_10_Lara_microC_downstream/out/TADs_boundaries_pileup/pileup_boundaries_<sample>.png
# =============================================================================