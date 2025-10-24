#!/usr/bin/env Rscript
# ======================================================================
# Script: cumulative_distance_analysis_fixed_paths.R
# Purpose: Same as the parametric version, but all paths are fixed
#          within the script. Edit the "USER SETTINGS" section below
#          to change input/output paths.
#
# HOW TO RUN THIS SCRIPT IN DOCKER (rocker/tidyverse:4.5.1)
# ----------------------------------------------------------------------
# 1) Install Docker: https://docs.docker.com/get-docker/
# 2) Pull the Rocker image with tidyverse (includes ggplot2 + data.table):
#      docker pull rocker/tidyverse:4.5.1
# 3) From your project folder (where this script and your data live), run:
#      sudo docker run -it -v /media/lucio/easystore:/media/lucio/easystore -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
#        Rscript Cumulative_freq_DEGs_L1_MLL2peaks_.R
#
# Notes:
#   - Requirements: data.table, ggplot2
#   - Local sources required:
#       ./git/downstream_multiomic/bin/hypergeometric.R
#       ./git/downstream_multiomic/bin/distance_functions.R
#       ./git/bioinfo_generics/base/bin/fun_compare_curves.R
# ======================================================================

library(ggplot2)

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
OUTDIR           <- "outs"

# Local functions
SRC_DIST <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/distance_functions.R"

# ------------------------------ PREP I/O --------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Source local functions (expects extractdistances(), filterbydistance(), etc.)
stopifnot(file.exists(SRC_DIST))
source(SRC_DIST)

# ------------------------------ LOAD INPUTS ------------------------------
# Annotated peaks and strand labels
stopifnot(file.exists(ANNOTATED_PEAKS_FILE))
antimml2.annotated.peaks <- read.table(ANNOTATED_PEAKS_FILE, sep = "\t", header = TRUE)

peaks.strands <- cbind(
  id = paste0(antimml2.annotated.peaks$Chr, "_",
              antimml2.annotated.peaks$Start, "_",
              antimml2.annotated.peaks$End),
  antimml2.annotated.peaks$Strand
)

# Raw peaks (BED-like)
stopifnot(file.exists(PEAKS_COORDS_FILE))
antimml2.peaks <- read.table(PEAKS_COORDS_FILE, sep = "\t", header = FALSE)
total.peaks <- nrow(antimml2.peaks)

# RepeatMasker table (UCSC format)
stopifnot(file.exists(RMSK_FILE))
repeat_mask <- read.table(RMSK_FILE, sep = "\t", header = FALSE)
# Filter LINE elements (family/type in V12; keep "LINE" and variants)
repeat_mask_line <- repeat_mask[grep("\\b(LINE|LINE?)\\b", repeat_mask$V12), ]

# Family proportions of LINEs in the genome (V11)
refline_prop <- table(repeat_mask_line$V11)[order(table(repeat_mask_line$V11), decreasing = TRUE)]
df_refline_prop <- data.frame(family = as.character(names(refline_prop)), n = as.vector(refline_prop))
if (nrow(df_refline_prop) > 12) {
  df_refline_prop[13:nrow(df_refline_prop), "family"] <- "Others"
  df_refline_prop <- rbind(
    df_refline_prop[1:12, ],
    data.frame(family = "Others",
               n = sum(df_refline_prop[df_refline_prop$family == "Others", "n"]))
  )
}

# Export a minimal BED of LINE elements for downstream usage
bedrm <- repeat_mask[, c("V6", "V7", "V8", "V11")]
write.table(
  bedrm,
  file = file.path(IN_UCSC_DIR, "l1_rmsk.bed"),
  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
)

# ------------------------------ DISTANCE EXTRACTION ----------------------
# Attribute and file suffixes to iterate (kept from original logic)
atb <- c("K4me3")
distances_list <- character(0)

# Helper to run extractdistances, density plot, and bookkeeping
run_distance_block <- function(attr, file_suffix) {
  bed_path <- file.path(IN_MLL2PEAKS_DIR, paste0(attr, file_suffix))
  stopifnot(file.exists(bed_path))
  d <- read.table(bed_path, sep = "\t")
  distances <- extractdistances(
    peaks.coords      = d,
    repeat_mask_line  = repeat_mask_line,
    addfamily = TRUE, addpeaks = TRUE, addlines = TRUE, add5prime = TRUE
  )
  # Track and assign into environment with original variable name
  distances_list <<- c(distances_list, paste0(attr, file_suffix))
  assign(paste0(attr, file_suffix), distances, envir = .GlobalEnv)

  # Density PDF with median line
  densval <- as.numeric(as.character(distances[, 1]))
  pdf(file.path(OUTDIR, paste0(attr, file_suffix, ".pdf")))
  plot(density(densval),
       main = paste0(attr, file_suffix, " nearest LINE"),
       sub  = paste0("median=", median(densval)))
  abline(v = median(densval))
  dev.off()
}

for (i in seq_along(atb)) {
  # distal - CpG minus / plus
  run_distance_block(atb[i], "_distal_CpG_minus_Double_KO_vs_F_F.bed")
  run_distance_block(atb[i], "_distal_CpG_plus_Double_KO_vs_F_F.bed")
  # proximal - CpG minus / plus
  run_distance_block(atb[i], "_proximal_CpG_minus_Double_KO_vs_F_F.bed")
  run_distance_block(atb[i], "_proximal_CpG_plus_Double_KO_vs_F_F.bed")
}

# ------------------------------ ECDF PLOT (all groups) -------------------
# Build combined df from the objects created above
final <- NULL
for (nm in distances_list) {
  pre <- data.frame(distance = as.numeric(as.character(get(nm)[, 1])),
                    group    = nm)
  final <- if (is.null(final)) pre else rbind(final, pre)
}
df <- final

# Fraction under 1kb for distal CpG minus
npeaks_dcm_under1kb <- length(which(df$group == "K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed" & df$distance <= 1000))
n_peaks_dcm_tot     <- length(which(df$group == "K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"))

# Drop global proximal (if present) as in original
df <- df[which(df$group != "K4me3_proximal_Double_KO_vs_F_F.bed"), ]

cum_dist <- ggplot(df, aes(distance, colour = group)) +
  stat_ecdf(geom = "step") +
  labs(title = "Cumulative Distance Distribution",
       y = "perc", x = "distance") +
  theme_classic() +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "#6a6a6a", size = 0.8) +
  annotate("text", x = 1000 + 2500, y = 1,
           label = npeaks_dcm_under1kb / n_peaks_dcm_tot,
           vjust = -0.5, color = "#6a6a6a", size = 4) +
  scale_color_manual(values = c(
    "K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"  = "#5b92cb",
    "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed"   = "#50c7ef",
    "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed"= "#ed2c85",
    "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed" = "#ad2f93"
  )) +
  theme(axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

ggsave(file.path(OUTDIR, "Cumulative_Density_plot.png"),
       plot = cum_dist, width = 13, height = 8, dpi = 300)

# ------------------------------ PAIRWISE TESTS ---------------------------
grps <- c("K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed",
          "K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed",
          "K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed",
          "K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed")

compair.curve <- function(a, b, kmatrix, grps, test = c("Wilcox", "KS")) {
  if (test[1] == "Wilcox") {
    a_VS_b <- wilcox.test(df$distance[df$group == a], df$distance[df$group == b])
  }
  if (test[1] == "KS") {
    a_VS_b <- ks.test(df$distance[df$group == a], df$distance[df$group == b])
  }
  kmatrix[a, b] <- a_VS_b$p.value
  kmatrix[b, a] <- a_VS_b$p.value
  return(kmatrix)
}

kmatrix <- matrix(nrow = length(grps), ncol = length(grps))
colnames(kmatrix) <- grps; rownames(kmatrix) <- grps

# Fill matrix (same pairings as original)
pairs_to_test <- list(
  c(grps[1], grps[2]), c(grps[1], grps[3]), c(grps[1], grps[4]),
  c(grps[2], grps[3]), c(grps[2], grps[4]), c(grps[3], grps[4]),
  c(grps[1], grps[1]), c(grps[2], grps[2]), c(grps[3], grps[3]), c(grps[4], grps[4])
)
for (p in pairs_to_test) {
  kmatrix <- compair.curve(p[1], p[2], kmatrix, grps)
}

kr.matrix <- kmatrix
colnames(kr.matrix) <- c("dcm", "dcp", "pcm", "pcp")
rownames(kr.matrix) <- c("dcm", "dcp", "pcm", "pcp")
write.table(
  kr.matrix,
  file = file.path(OUTDIR, "Cumulative_Density_pval_wilcox.tsv"),
  col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t"
)

boxdf <- df  # preserved name from original

# ------------------------------ LINE-1 FAMILY PROFILES -------------------
repeat_mask_line <- repeat_mask_line[, c(6, 7, 8, 10, 11, 12, 13)]
colnames(repeat_mask_line) <- c("chr", "start", "end", "strand", "family", "element", "type")

repeat_mask_line_plus  <- repeat_mask_line[which(repeat_mask_line$strand == "+"), ]
repeat_mask_line_minus <- repeat_mask_line[which(repeat_mask_line$strand == "-"), ]

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
# Select families ≥ 1st quartile (as in original approach)
su.line.freqs <- line.freqs[which(line.freqs$Freq >= as.numeric(summary(line.freqs$Freq)[2])), ]

top.line.freqs <- line.freqs[order(line.freqs$Freq, decreasing = TRUE)[1:7], ]

other.line.freqs <- line.freqs[order(line.freqs$Freq, decreasing = TRUE), ]
levels(other.line.freqs[, "Var1"]) <- c(levels(other.line.freqs[, "Var1"]), "Others")
other.line.freqs[8:nrow(other.line.freqs), "Var1"] <- as.factor("Others")

# distal CpG minus subset frequencies
dgm.line.freqs <- as.data.frame(table(K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[, "family"]))
su.dgm.line.fraqs  <- dgm.line.freqs[which(dgm.line.freqs$Freq >= as.numeric(summary(dgm.line.freqs$Freq)[2])), ]
top.dgm.line.fraqs <- dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE)[1:7], ]

other.dgm.line.freqs <- dgm.line.freqs[order(dgm.line.freqs$Freq, decreasing = TRUE), ]
levels(other.dgm.line.freqs[, "Var1"]) <- c(levels(other.dgm.line.freqs[, "Var1"]), "Others")
other.dgm.line.freqs[8:nrow(other.dgm.line.freqs), "Var1"] <- as.factor("Others")

laneFreq.dataplot <- as.data.frame(rbind(
  cbind(su.line.freqs,  group = rep("all",      nrow(su.line.freqs))),
  cbind(su.dgm.line.fraqs, group = rep("only.dcm", nrow(su.dgm.line.fraqs)))
))

top.laneFreq.dataplot <- as.data.frame(rbind(
  cbind(top.line.freqs,      group = rep("all",      nrow(top.line.freqs))),
  cbind(top.dgm.line.fraqs,  group = rep("only.dcm", nrow(top.dgm.line.fraqs)))
))

other.laneFreq.dataplot <- as.data.frame(rbind(
  cbind(other.line.freqs,    group = rep("all",      nrow(other.line.freqs))),
  cbind(other.dgm.line.freqs, group = rep("only.dcm", nrow(other.dgm.line.freqs)))
))

# Stacked bar (others vs groups)
p_comp <- ggplot(other.laneFreq.dataplot, aes(fill = Var1, y = Freq, x = group)) +
  geom_bar(position = "fill", stat = "identity")
ggsave(file.path(OUTDIR, "LINE_family_composition.png"),
       plot = p_comp, width = 10, height = 6, dpi = 300)

# ------------------------------ EXPORT TARGET SETS -----------------------
# dcm subset (distal CpG minus): L1 coords and peak coords
tarfam     <- as.character(other.laneFreq.dataplot[which(other.laneFreq.dataplot$group == "only.dcm" &
                                                         other.laneFreq.dataplot$Var1 != "Others"), "Var1"])
tarfam_anno <- K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[, c("line_chr", "line_start", "line_end", "family", "distances", "line_strand")]
tarfam_anno[, "family"] <- paste0(tarfam_anno[, "family"], "_", seq_len(nrow(tarfam_anno)))

tarpeak <- K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[, c("peaks_chr", "peaks_start", "peaks_end")]
tarpeak <- cbind(
  tarpeak,
  name = paste0("peak_", seq_len(nrow(tarpeak))),
  score = K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed[, "distances"],
  strand = rep(".", nrow(tarpeak))
)

write.table(tarfam_anno, file = file.path(OUTDIR, "dcm.l1.bed"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(tarpeak, file = file.path(OUTDIR, "dcm.target.peaks.bed"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# dcp subset (distal CpG plus)
dcp_anno <- K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[, c("line_chr", "line_start", "line_end", "family", "distances", "line_strand")]
dcp_anno[, "family"] <- paste0(dcp_anno[, "family"], "_", seq_len(nrow(dcp_anno)))

dcp_tarpeak <- K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[, c("peaks_chr", "peaks_start", "peaks_end")]
dcp_tarpeak <- cbind(
  dcp_tarpeak,
  name  = paste0("peak_", seq_len(nrow(dcp_tarpeak))),
  score = K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed[, "distances"],
  strand = rep(".", nrow(dcp_tarpeak))
)

write.table(dcp_anno, file = file.path(OUTDIR, "dcp.l1.bed"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(dcp_tarpeak, file = file.path(OUTDIR, "dcp.target.peaks.bed"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# ------------------------------ FILTER BY DISTANCE -----------------------
# Produce filtered sets at various thresholds (mirror original calls)

filterbydistance(extractdistances.out = K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist = 500, out = file.path(OUTDIR, "distfilter500_DKO_K4me3_dcm"))
filterbydistance(extractdistances.out = K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist = 500, out = file.path(OUTDIR, "distfilter500_DKO_K4me3_dcp"))
filterbydistance(extractdistances.out = K4me3_proximal_CpG_plus_Double_KO_vs_F_F.bed,
                 fildist = 500, out = file.path(OUTDIR, "distfilter500_DKO_K4me3_pcp"))
filterbydistance(extractdistances.out = K4me3_proximal_CpG_minus_Double_KO_vs_F_F.bed,
                 fildist = 500, out = file.path(OUTDIR, "distfilter500_DKO_K4me3_pcm"))


# ------------------------------ RECENTER ON 5' ---------------------------
disf_dcm_complete_path <- file.path(OUTDIR, "distfilter500_DKO_K4me3_dcm.l1.bed")
stopifnot(file.exists(disf_dcm_complete_path))
disf_dcm_complete <- read.table(disf_dcm_complete_path, sep = "\t", header = FALSE)

center_win <- 5000
for (i in seq_len(nrow(disf_dcm_complete))) {
  if (disf_dcm_complete[i, 6] == "+") {
    fprime <- disf_dcm_complete[i, 2]
  } else if (disf_dcm_complete[i, 6] == "-") {
    fprime <- disf_dcm_complete[i, 3]
  } else {
    fprime <- round(mean(c(disf_dcm_complete[i, 2], disf_dcm_complete[i, 3])))
  }
  disf_dcm_complete[i, 2] <- fprime - center_win
  disf_dcm_complete[i, 3] <- fprime + center_win
}

write.table(disf_dcm_complete,
            file = file.path(OUTDIR, "recentered_distfil500_DKO_K4me3_dcm.l1.bed"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# ------------------------------ ECDF FROM FILTERED SETS ------------------
distances_files <- c(
  "distfilter500_DKO_K4me3_dcm.l1.bed",
  "distfilter500_DKO_K4me3_dcp.l1.bed",
  "distfilter500_DKO_K4me3_pcm.l1.bed",
  "distfilter500_DKO_K4me3_pcp.l1.bed"
)

distances_list <- distances_list[which(distances_list != "K4me3_proximal_Double_KO_vs_F_F.bed")]

final <- NULL
for (i in seq_along(distances_files)) {
  distvar <- read.table(file.path(OUTDIR, distances_files[i]), sep = "\t", header = FALSE)
  pre <- data.frame(distance = as.numeric(as.character(distvar[, 5])),
                    group = distances_files[i])
  final <- if (is.null(final)) pre else rbind(final, pre)
}
df <- final

p_ecdf <- ggplot(df, aes(distance, colour = group)) + stat_ecdf(geom = "step") +
  labs(title = "Empirical Cumulative \n Density Function",
       y = "perc", x = "distance") +
  theme_classic()

ggsave(file.path(OUTDIR, "ECDF_filtered_sets.png"), plot = p_ecdf, width = 10, height = 6, dpi = 300)


message("All steps completed. Outputs under: ", normalizePath(OUTDIR))
