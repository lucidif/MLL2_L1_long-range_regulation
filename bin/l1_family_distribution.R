#!/usr/bin/env Rscript
# ======================================================================
# Script: l1_family_distribution_REFACTORED.R
# Purpose: Compute LINE-1 family composition and compare overall genome
#          abundances with a distal CpG- subset of K4me3 peaks.
# Style:   All paths defined in USER SETTINGS and referenced via variables.
#          No setwd(). Same analysis logic as the original.
#
# HOW TO RUN IN DOCKER (rocker/tidyverse:4.5.1)
# ----------------------------------------------------------------------
# 1) Install Docker: https://docs.docker.com/get-docker/
# 2) Pull the image:
#      docker pull rocker/tidyverse:4.5.1
# 3) From your project folder (where this script + data live), run:
#      sudo docker run -it -v /media/lucio/easystore:/media/lucio/easystore -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
#        Rscript l1_family_distribution.R
#
# Requires:
#   - Packages: ggplot2, RColorBrewer, dplyr, tidyr
#   - Local source: distance_functions.R (provides extractdistances)
# Outputs:
#   - Density PDFs for each peak class
#   - L1 family frequency barplot (PDF) and TSV table
# ======================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

# ------------------------------ USER SETTINGS ---------------------------
# Reference / annotations
REF_FASTA        <- "/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa"
IN_UCSC_DIR      <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc"
RMSK_FILE        <- "git/MLL2_LINE1_long-range_regulation/assets/rmsk.txt"

# Peak inputs
ANNOTATED_PEAKS_FILE <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/macs_broadpeaks/Anti-GFP.mLb.mkD.sorted_peaks.annotatePeaks.txt"
PEAKS_COORDS_FILE    <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/bedtools_window/coordinate.bed"

# Folders with peak sets and outputs
IN_DIFFPEAKS_DIR <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/diffbind"
IN_MLL2PEAKS_DIR <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO"
OUTDIR_L1        <- "outs/CHiP_postprocessing_line1_dist"

# Local functions
SRC_DIST <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/distance_functions.R"

# ------------------------------ PREP I/O --------------------------------
dir.create(OUTDIR_L1, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(SRC_DIST))
source(SRC_DIST)   # expects: extractdistances()

# ------------------------------ INPUTS ----------------------------------
stopifnot(file.exists(ANNOTATED_PEAKS_FILE))
antimml2.annotated.peaks <- read.table(ANNOTATED_PEAKS_FILE, sep = "\t", header = TRUE)

peaks.strands <- cbind(
  id = paste0(antimml2.annotated.peaks$Chr, "_",
              antimml2.annotated.peaks$Start, "_",
              antimml2.annotated.peaks$End),
  antimml2.annotated.peaks$Strand
)

stopifnot(file.exists(PEAKS_COORDS_FILE))
antimml2.peaks <- read.table(PEAKS_COORDS_FILE, sep = "\t", header = FALSE)
total.peaks <- nrow(antimml2.peaks)

# ------------------------------ REPEATMASKER -----------------------------
stopifnot(file.exists(RMSK_FILE))
repeat_mask <- read.table(RMSK_FILE, sep = "\t", header = FALSE)

# Keep only LINE elements (V12 is type/repClass)
repeat_mask_line <- repeat_mask[grep("\\b(LINE|LINE?)\\b", repeat_mask$V12), ]

# Genome-wide LINE family proportions (V11 is family)
refline_prop <- table(repeat_mask_line$V11)[order(table(repeat_mask_line$V11), decreasing = TRUE)]
df_refline_prop <- data.frame(refline_prop)
colnames(df_refline_prop) <- c("family", "n")

# Minimal BED export of LINE elements (if needed downstream)
bedrm <- repeat_mask[, c("V6", "V7", "V8", "V11")]
# write.table(bedrm, file = file.path(IN_UCSC_DIR, "l1_rmsk.bed"),
#             sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# ------------------------------ DISTANCES --------------------------------
# Compute distances from peak subclasses to nearest LINE-1
atb <- c("K4me3")
distances_list <- character(0)

run_distance_block <- function(attr, file_pattern) {
  bed_path <- file.path(IN_MLL2PEAKS_DIR, paste0(attr, file_pattern))
  stopifnot(file.exists(bed_path))
  d <- read.table(bed_path, sep = "\t")
  distances <- extractdistances(
    peaks.coords     = d,
    repeat_mask_line = repeat_mask_line,
    addfamily = TRUE, addpeaks = TRUE, addlines = TRUE, add5prime = TRUE
  )
  distances_list <<- c(distances_list, paste0(attr, file_pattern))
  assign(paste0(attr, file_pattern), distances, envir = .GlobalEnv)

  densval <- as.numeric(as.character(distances[, 1]))
  pdf(file.path(OUTDIR_L1, paste0(attr, file_pattern, ".pdf")))
  plot(density(densval),
       main = paste0(attr, file_pattern, " nearest LINE"),
       sub  = paste0("median=", median(densval)))
  abline(v = median(densval))
  dev.off()
}

for (i in seq_along(atb)) {
  run_distance_block(atb[i], "_distal_CpG_minus_Double_KO_vs_F_F.bed")
  run_distance_block(atb[i], "_distal_CpG_plus_Double_KO_vs_F_F.bed")
  run_distance_block(atb[i], "_proximal_CpG_minus_Double_KO_vs_F_F.bed")
  run_distance_block(atb[i], "_proximal_CpG_plus_Double_KO_vs_F_F.bed")
}

# ------------------------------ ECDF PREP --------------------------------
final <- NULL
for (nm in distances_list) {
  pre <- data.frame(distance = as.numeric(as.character(get(nm)[, 1])),
                    group    = nm)
  final <- if (is.null(final)) pre else rbind(final, pre)
}
df <- final

# ------------------------------ L1 FAMILY DISTRIBUTION -------------------
# Keep only L1 (type column becomes 'type' later) and transform to tidy frame
repeat_mask_line <- repeat_mask_line[, c(6, 7, 8, 10, 11, 12, 13)]
colnames(repeat_mask_line) <- c("chr", "start", "end", "strand", "family", "element", "type")

repeat_mask_line_plus  <- subset(repeat_mask_line, strand == "+")
repeat_mask_line_minus <- subset(repeat_mask_line, strand == "-")

repeat_mask.5prime <- rbind(
  cbind(repeat_mask_line_plus,  fprime = repeat_mask_line_plus$start),
  cbind(repeat_mask_line_minus, fprime = repeat_mask_line_minus$end)
)

repeat_mask.coordinates_extended <- cbind(
  repeat_mask.5prime$chr,
  repeat_mask.5prime$fprime,
  repeat_mask.5prime$family
)

line.freqs <- as.data.frame(table(repeat_mask.coordinates_extended[, 3]))
# Distal CpG minus subset families
dgm.line.freqs <- as.data.frame(table(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[, "family"]))
# Subset to L1 families in reference only
dgm.line.freqs <- dgm.line.freqs[dgm.line.freqs$Var1 %in% line.freqs$Var1, ]

# Merge and order by descending DCM frequencies
full.table.lines.abundancy <- merge(line.freqs, dgm.line.freqs, by = 1, all.x = TRUE)
full.table.lines.abundancy.sorted <- full.table.lines.abundancy[order(-full.table.lines.abundancy$Freq.y, na.last = TRUE), ]
to_export <- full.table.lines.abundancy.sorted
colnames(to_export) <- c("family", "overall", "dcm_loss_k4me3")

# Replace NA in DCM column with 0
full.table.lines.abundancy.sorted$Freq.y[is.na(full.table.lines.abundancy.sorted$Freq.y)] <- 0

# Top-5 families by DCM frequency; group the rest as "Others"
top5 <- full.table.lines.abundancy.sorted$Var1[order(-full.table.lines.abundancy.sorted$Freq.y)][1:5]
full.table.lines.abundancy.sorted$Var1 <- as.character(full.table.lines.abundancy.sorted$Var1)
full.table.lines.abundancy.sorted$Var1[!(full.table.lines.abundancy.sorted$Var1 %in% top5)] <- "Others"
full.table.lines.abundancy.sorted$Var1 <- factor(full.table.lines.abundancy.sorted$Var1)

# Build percent table for stacked/side-by-side plotting
df_plot <- full.table.lines.abundancy.sorted %>%
  group_by(Var1) %>%
  summarise(
    Freq.x = sum(Freq.x),
    Freq.y = sum(Freq.y)
  ) %>%
  mutate(
    perc_x = Freq.x / sum(Freq.x) * 100,
    perc_y = Freq.y / sum(Freq.y) * 100
  ) %>%
  select(Var1, perc_x, perc_y)

df_long <- df_plot %>%
  pivot_longer(
    cols = c(perc_x, perc_y),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    Category = recode(Category,
                      perc_x = "Overall",
                      perc_y = "DCM_only")
  )

df_long$Category <- factor(df_long$Category, levels = c("Overall", "DCM_only"))
df_long$Var1 <- factor(df_long$Var1, levels = c("Others", setdiff(unique(df_long$Var1), "Others")))

# Plot
p_famfreq <- ggplot(df_long, aes(x = Var1, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Family abundancy", x = "Family", y = "Percentage", fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR_L1, "families_frequency_plot.pdf"),
       plot = p_famfreq, width = 10, height = 6)

# Export TSV
write.table(to_export,
            file = file.path(OUTDIR_L1, "families_frequency.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

message("Done. Outputs under: ", normalizePath(OUTDIR_L1))


