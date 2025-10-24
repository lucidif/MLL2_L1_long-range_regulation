#!/usr/bin/env Rscript
# ======================================================================
# Script: cumulative_distance_analysis_fixed_paths.R
#
# HOW TO RUN THIS SCRIPT IN DOCKER (rocker/tidyverse:4.5.1)
# ----------------------------------------------------------------------
# 1) Install Docker: https://docs.docker.com/get-docker/
# 2) Pull the Rocker image with tidyverse (includes ggplot2 + data.table):
#      docker pull rocker/tidyverse:4.5.1
# 3) From your project folder (where this script and your data live), run:
#      sudo docker run -it -v /media/lucio/easystore:/media/lucio/easystore -v "$PWD":"$PWD" -w "$PWD" rocker/tidyverse:4.5.1 \
#        Rscript cumulative_distance_analysis_fixed_paths.R
#
# Notes:
#   - Requirements: data.table, ggplot2
#   - Local sources required:
#       ./git/downstream_multiomic/bin/hypergeometric.R
#       ./git/downstream_multiomic/bin/distance_functions.R
#       ./git/bioinfo_generics/base/bin/fun_compare_curves.R
# ======================================================================

install.packages("R.utils")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------ USER SETTINGS ---------------------------
# Edit paths here as needed.

# RNA-seq tables & annotation
IN_RNASEQ_DIR <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/build38_DEseq2_RNAseq"
IN_ANNO_DIR   <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/build38_DEseq2_RNAseq"
REFSEQ_PATH   <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/ucsc/build38_mm10_ncbiRefSeqCurated.txt.gz"

# ChIP-seq peak sets (distal CpG +/-)
PEAKS_DCM_PATH <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_minus_Double_KO_vs_F_F.bed"
PEAKS_DCP_PATH <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/in/test_chipseq_downstream/deeptools_heatmaps_filtered_by_diffDoubleKO/K4me3_distal_CpG_plus_Double_KO_vs_F_F.bed"

# Subsets near H3K27ac (D4)
NEAR_K27AC_DCM_PATH <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcm_lossk4me3_near_D4_k27ac.bed"
NEAR_K27AC_DCP_PATH <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/dcp_lossk4me3_near_k27ac_heatmap/thr5_dcp_lossk4me3_near_D4_k27ac.bed"

# L1 coordinate tiers
L1_ALL_PATH  <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/gc_content_heatmap/recentered_distfil500_DKO_K4me3_dcm.l1.bed"
L1_HIGH_PATH <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/l1_d4_heatmap/high_sort2_l1_d4mll2_plotHeatmap.bed"
L1_MID_PATH  <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/l1_d4_heatmap/mid_sort2_l1_d4mll2_plotHeatmap.bed"
L1_LOW_PATH  <- "/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/l1_d4_heatmap/low_sort2_l1_d4mll2_plotHeatmap.bed"

# Output folder
OUTDIR <- "outs"

# Local sources (functions)
SRC_HYPER  <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/hypergeometric.R"
SRC_DIST   <- "./git/MLL2_LINE1_long-range_regulation/pipelines/downstream_multiomic/bin/distance_functions.R"
SRC_CURVES <- "./git/MLL2_LINE1_long-range_regulation/pipelines/bioinfo_generics/base/bin/fun_compare_curves.R"

# ------------------------------ PREP I/O --------------------------------
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "distribution_distances_downDEGs_dcp"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "distribution_distances_downDEGS_dcm"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "l1_d4_heatmap"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------ Source local functions -------------------
stopifnot(file.exists(SRC_HYPER), file.exists(SRC_DIST), file.exists(SRC_CURVES))

source(SRC_HYPER)
source(SRC_DIST)
source(SRC_CURVES)

# ------------------------------ Helpers ---------------------------------
read_tsv <- function(path, ...) {
  stopifnot(file.exists(path))
  data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, ...)
}

read_bed <- function(path) {
  stopifnot(file.exists(path))
  df <- data.table::fread(path, sep = "\t", header = FALSE, data.table = FALSE)
  stopifnot(ncol(df) >= 3)
  df
}

label_bed <- function(bed_df, prefix) {
  if (ncol(bed_df) < 3) stop("BED requires >=3 columns")
  out <- data.frame(
    V1 = bed_df[[1]],
    V2 = bed_df[[2]],
    V3 = bed_df[[3]],
    V4 = paste0(prefix, "_", seq_len(nrow(bed_df))),
    V5 = NA,
    V6 = ".",
    stringsAsFactors = FALSE
  )
  colnames(out) <- c("V1","V2","V3","V4","V5","V6")
  out
}

ensure_numeric <- function(x) as.numeric(as.character(x))

exportBed <- function(distancePlotResult, file) {
  D0_down_pro <- distancePlotResult
  chr <- start <- end <- direction <- c()
  for (i in seq_len(nrow(D0_down_pro$dis.CpG.plus.degs))) {
    chr[i] <- D0_down_pro$dis.CpG.plus.degs[i, 1]
    if (D0_down_pro$dis.CpG.plus.degs[i,"from_position"] <= D0_down_pro$dis.CpG.plus.degs[i,"to_position"]) {
      start[i]     <- round(D0_down_pro$dis.CpG.plus.degs[i,"from_position"], 0)
      end[i]       <- round(D0_down_pro$dis.CpG.plus.degs[i,"to_position"],   0)
      direction[i] <- "peak.right"
    } else {
      start[i]     <- round(D0_down_pro$dis.CpG.plus.degs[i,"to_position"],   0)
      end[i]       <- round(D0_down_pro$dis.CpG.plus.degs[i,"from_position"], 0)
      direction[i] <- "peak.left"
    }
  }
  bed <- cbind(chr, start, end, direction)
  write.table(bed, file=file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

make_ecdf_plot <- function(df_windowed, title, palette_named, outfile_png, width=8, height=6, dpi=300) {
  df_windowed$distances <- ensure_numeric(df_windowed$distances)
  p <- ggplot(df_windowed, aes(distances, colour = group)) +
    stat_ecdf(geom = "step") +
    labs(title = title, y = "cumulative fraction", x = "distance (bp)") +
    theme_classic() +
    scale_color_manual(values = palette_named)
  ggsave(outfile_png, p, width = width, height = height, dpi = dpi, units = "in")
}

run_wilcox_matrix <- function(df_windowed, outfile_tsv) {
  df_windowed<-as.data.frame(df_windowed)
  df_windowed$distances<-as.numeric(as.character(df_windowed$distances))
  wilcox.matrix(df_windowed = df_windowed, file = outfile_tsv)
}

# ------------------------------ Load inputs ------------------------------
degs_D4     <- read_tsv(file.path(IN_RNASEQ_DIR, "D4_DKO_VS_WT.deseq2.results.tsv"))
degs_D0     <- read_tsv(file.path(IN_RNASEQ_DIR, "D0_DKO_vs_WT.deseq2.results.tsv"))
degs_wtD4D0 <- read_tsv(file.path(IN_RNASEQ_DIR, "D4WT_VS_D0WT.deseq2.results.tsv"))

anno <- read_tsv(file.path(IN_ANNO_DIR, "mm10.anno.tsv"))

refseq.anno <- data.table::fread(REFSEQ_PATH, sep="\t", header=FALSE, data.table=FALSE)
refseq.anno[,"V5"] <- refseq.anno[,"V5"] + 1
refseq.anno <- refseq.anno[, c("V13","V3","V5","V6","V4")]

refseq.anno.stplus <- refseq.anno[refseq.anno$V4 == "+", ]
refseq.anno.stplus <- refseq.anno.stplus[!duplicated(refseq.anno.stplus[, c("V13","V3","V5")]), ]

refseq.anno.stmin  <- refseq.anno[refseq.anno$V4 == "-", ]
refseq.anno.stmin  <- refseq.anno.stmin[!duplicated(refseq.anno.stmin[, c("V13","V3","V6")]), ]

refseq.anno <- rbind(refseq.anno.stplus, refseq.anno.stmin)

refseq.anno <- merge(anno[, c("gene_name","gene_id")], refseq.anno, all.y=TRUE, by.x="gene_name", by.y="V13")
refseq.anno <- refseq.anno[!duplicated(refseq.anno), ]
refseq.anno <- refseq.anno[!is.na(refseq.anno$gene_id), ]
refseq.anno <- refseq.anno[grep("_.", refseq.anno$V3, invert = TRUE), ]

refseq.anno.reform <- data.frame(
  GeneID    = refseq.anno$gene_id,
  uniq      = paste0(refseq.anno$gene_id,"_", seq_len(nrow(refseq.anno))),
  GeneName  = refseq.anno$gene_name,
  Chromosome= refseq.anno$V3,
  start     = refseq.anno$V5,
  end       = refseq.anno$V6,
  strand    = refseq.anno$V4,
  stringsAsFactors = FALSE
)

norm.counts <- read_tsv(file.path(IN_RNASEQ_DIR, "all.normalised_counts.tsv"))

gene.annWithPosition <- refseq.anno.reform
colnames(gene.annWithPosition) <- c("GeneID","uniq","GeneName","Chromosome","start","end","strand")

# ------------------------------ Peaks & neighborhoods --------------------
dcm_peak.near.ref <- read_bed(PEAKS_DCM_PATH)
dcp_peaks.near.ref <- read_bed(PEAKS_DCP_PATH)
dcp_peaks.near.ref <- label_bed(dcp_peaks.near.ref, "dcp_peak")

dcm_near_k27ac <- read_bed(NEAR_K27AC_DCM_PATH)
dcm_near_k27ac <- label_bed(dcm_near_k27ac, "dcm_near_k27ac")

dcp_near_k27ac <- read_bed(NEAR_K27AC_DCP_PATH)
dcp_near_k27ac <- label_bed(dcp_near_k27ac, "dcp_near_k27ac")

# L1 sets
l1_all  <- read_bed(L1_ALL_PATH)
l1_high <- read_bed(L1_HIGH_PATH)
l1_mid  <- read_bed(L1_MID_PATH)
l1_low  <- read_bed(L1_LOW_PATH)

# ------------------------------ Annotate DEGs ----------------------------
degs_D4_ann_diff_down     <- addAnnotationToDegs(degs_D4, gene.annWithPosition, direction="down", log2fc=-1)
degs_D0_ann_diff_down     <- addAnnotationToDegs(degs_D0, gene.annWithPosition, direction="down", log2fc=-1)
degs_wtD4D0_ann_diff_down <- addAnnotationToDegs(degs_wtD4D0, gene.annWithPosition, direction="down", log2fc=-1)

degs_D4_ann_diff_up       <- addAnnotationToDegs(degs_D4, gene.annWithPosition, direction="up",   log2fc=1)
degs_D0_ann_diff_up       <- addAnnotationToDegs(degs_D0, gene.annWithPosition, direction="up",   log2fc=1)
degs_wtD4D0_ann_diff_up   <- addAnnotationToDegs(degs_wtD4D0, gene.annWithPosition, direction="up", log2fc=1)

allgenes_D0_ann_diff_all <- addAnnotationToDegs(degs_D0, gene.annWithPosition, padj=1, log2fc=0, direction="all")
allgenes_D4_ann_diff_all <- addAnnotationToDegs(degs_D4, gene.annWithPosition, padj=1, log2fc=0, direction="all")

nodiff_D0 <- allgenes_D0_ann_diff_all[allgenes_D0_ann_diff_all$padj > 0.05, ]
nodiff_D4 <- allgenes_D4_ann_diff_all[allgenes_D4_ann_diff_all$padj > 0.05, ]

overlap_downdko_upwt <- intersect(degs_wtD4D0_ann_diff_up$ensembl, degs_D4_ann_diff_down$ensembl)
D4_downdko_upwt      <- degs_D4_ann_diff_down[degs_D4_ann_diff_down$ensembl %in% overlap_downdko_upwt, ]

# ------------------------------ D0: Distal CpG plus/minus ----------------
D0_down_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                                 plus_peaks  = dcp_peaks.near.ref,
                                 degs        = degs_D0_ann_diff_down,
                                 title       = "distal peaks near line D0 down",
                                 window      = 1e6,
                                 window.type = "inside")

D0_nodiff_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                                   plus_peaks  = dcp_peaks.near.ref,
                                   degs        = nodiff_D0,
                                   title       = "distal peaks near line D0 all annotated genes",
                                   window      = 1e6,
                                   window.type = "inside")

D0_up_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                               plus_peaks  = dcp_peaks.near.ref,
                               degs        = degs_D0_ann_diff_up,
                               title       = "distal peaks near line D0 up",
                               window      = 1e6,
                               window.type = "inside")

df_windowed <- as.data.frame(rbind(
  cbind(distances = D0_down_dis_win$dis.CpG.plus.degs$distances,   group = "dis.CpG.plus_DOWN"),
  cbind(distances = D0_up_dis_win$dis.CpG.plus.degs$distances,     group = "dis.CpG.plus_UP"),
  cbind(distances = D0_nodiff_dis_win$dis.CpG.plus.degs$distances, group = "dis.CpG.plus_NODIFF")
))

make_ecdf_plot(
  df_windowed,
  title = "D0 Distance Cumulative Distribution",
  palette_named = c(
    "dis.CpG.plus_DOWN"   = "#3f3092",
    "dis.CpG.plus_UP"     = "#608dd1",
    "dis.CpG.plus_NODIFF" = "#bdbbba"
  ),
  outfile_png = file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "dis_downDEGs_dcp.png")
)
run_wilcox_matrix(df_windowed, file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "dis_downDEGs_dcp_pval_wilcox.tsv"))

# ------------------------------ D4: Distal CpG plus ----------------------
D4_down_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                                 plus_peaks  = dcp_peaks.near.ref,
                                 degs        = degs_D4_ann_diff_down,
                                 title       = "distal peaks near line D4 down",
                                 window      = 1e6,
                                 window.type = "inside")

D4_nodiff_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                                   plus_peaks  = dcp_peaks.near.ref,
                                   degs        = nodiff_D4,
                                   title       = "distal peaks near line D4 all annotated genes",
                                   window      = 1e6,
                                   window.type = "inside")

D4_up_dis_win <- distance_plot(minus_peaks = dcm_peak.near.ref,
                               plus_peaks  = dcp_peaks.near.ref,
                               degs        = degs_D4_ann_diff_up,
                               title       = "distal peaks near line D4 up",
                               window      = 1e6,
                               window.type = "inside")

D4k27ac_downDKO_wtD4D0up <- distance_plot(minus_peaks = dcm_near_k27ac,
                                           plus_peaks  = dcp_near_k27ac,
                                           degs        = D4_downdko_upwt,
                                           title       = "distal peaks near line (K27ac subset)",
                                           window      = 1e6,
                                           window.type = "inside")

df_d4_plus_windowed <- as.data.frame(rbind(
  cbind(distances = D4_down_dis_win$dis.CpG.plus.degs$distances,               group = "dis.CpG.plus_DOWN"),
  cbind(distances = D4_up_dis_win$dis.CpG.plus.degs$distances,                 group = "dis.CpG.plus_UP"),
  cbind(distances = D4_nodiff_dis_win$dis.CpG.plus.degs$distances,             group = "dis.CpG.plus_NODIFF"),
  cbind(distances = D4k27ac_downDKO_wtD4D0up$dis.CpG.plus.degs$distances,      group = "induced.dis.CpG.plus_DOWN")
))

make_ecdf_plot(
  df_d4_plus_windowed,
  title = "D4 Distance Cumulative Distribution",
  palette_named = c(
    "dis.CpG.plus_DOWN"           = "#3f3092",
    "dis.CpG.plus_UP"             = "#608dd1",
    "dis.CpG.plus_NODIFF"         = "#bdbbba",
    "induced.dis.CpG.plus_DOWN"   = "#ee2c85"
  ),
  outfile_png = file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "D4_distribution_distances_downDEGs_dcp.png")
)
run_wilcox_matrix(df_d4_plus_windowed, file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "D4_dis_downDEGs_dcp_pval_wilcox.tsv"))

# ------------------------------ D4: CpG plus near K27ac only -------------
D4k27ac_down_dis_win <- distance_plot(minus_peaks = dcm_near_k27ac,
                                      plus_peaks  = dcp_near_k27ac,
                                      degs        = degs_D4_ann_diff_down,
                                      title       = "K27ac-near distal peaks (D4 down)",
                                      window      = 1e6,
                                      window.type = "inside")

D4k27ac_nodiff_dis_win <- distance_plot(minus_peaks = dcm_near_k27ac,
                                        plus_peaks  = dcp_near_k27ac,
                                        degs        = nodiff_D4,
                                        title       = "K27ac-near distal peaks (D4 nodiff)",
                                        window      = 1e6,
                                        window.type = "inside")

D4k27ac_up_dis_win <- distance_plot(minus_peaks = dcm_near_k27ac,
                                    plus_peaks  = dcp_near_k27ac,
                                    degs        = degs_D4_ann_diff_up,
                                    title       = "K27ac-near distal peaks (D4 up)",
                                    window      = 1e6,
                                    window.type = "inside")

df_d4_near_k27ac_windowed <- as.data.frame(rbind(
  cbind(distances = D4k27ac_down_dis_win$dis.CpG.plus.degs$distances,         group = "dis.CpG.plus_lossk4me3_near_k27ac_DOWN"),
  cbind(distances = D4k27ac_up_dis_win$dis.CpG.plus.degs$distances,           group = "dis.CpG.plus_lossk4me3_near_k27ac_UP"),
  cbind(distances = D4k27ac_nodiff_dis_win$dis.CpG.plus.degs$distances,       group = "dis.CpG.plus_lossk4me3_near_k27ac_NODIFF"),
  cbind(distances = D4k27ac_downDKO_wtD4D0up$dis.CpG.plus.degs$distances,     group = "dis.CpG.plus_lossk4me3_near_k27ac_induced_genes_DOWN")
))

make_ecdf_plot(
  df_d4_near_k27ac_windowed,
  title = "D4 Distance Cumulative Distribution (K27ac-near)",
  palette_named = c(
    "dis.CpG.plus_lossk4me3_near_k27ac_DOWN"             = "#3f3092",
    "dis.CpG.plus_lossk4me3_near_k27ac_UP"               = "#608dd1",
    "dis.CpG.plus_lossk4me3_near_k27ac_NODIFF"           = "#bdbbba",
    "dis.CpG.plus_lossk4me3_near_k27ac_induced_genes_DOWN"= "#ed2c85"
  ),
  outfile_png = file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "D4_distri_distan_downDEGs_dcp_nearK27ac.png")
)
run_wilcox_matrix(df_d4_near_k27ac_windowed, file.path(OUTDIR, "distribution_distances_downDEGs_dcp", "D4_distri_distan_downDEGs_dcp_nearK27ac_pval_wilcox.tsv"))

# ------------------------------ D4: CpG minus ----------------------------
df_d4_minus_windowed <- as.data.frame(rbind(
  cbind(distances = D4_down_dis_win$dis.CpG.minus.degs$distances,   group = "dis.CpG.minus_DOWN"),
  cbind(distances = D4_up_dis_win$dis.CpG.minus.degs$distances,     group = "dis.CpG.minus_UP"),
  cbind(distances = D4_nodiff_dis_win$dis.CpG.minus.degs$distances, group = "dis.CpG.minus_NODIFF")
))

make_ecdf_plot(
  df_d4_minus_windowed,
  title = "D4 Distance Cumulative Distribution (CpG minus)",
  palette_named = c(
    "dis.CpG.minus_DOWN"   = "#3f3092",
    "dis.CpG.minus_UP"     = "#608dd1",
    "dis.CpG.minus_NODIFF" = "#bdbbba"
  ),
  outfile_png = file.path(OUTDIR, "distribution_distances_downDEGS_dcm", "D4_distribution_distances_downDEGs_dcm.png")
)

# ------------------------------ D4: Distances to L1 sets -----------------
run_l1_block <- function(l1_bed, tag) {
  l1near_down <- distance_plot(minus_peaks = l1_bed, plus_peaks = l1_bed,
                               degs = degs_D4_ann_diff_down, title = paste0("L1 near (", tag, ") - down"),
                               window = 1e6, window.type = "inside")
  l1near_up <- distance_plot(minus_peaks = l1_bed, plus_peaks = l1_bed,
                             degs = degs_D4_ann_diff_up, title = paste0("L1 near (", tag, ") - up"),
                             window = 1e6, window.type = "inside")
  l1near_undiff <- distance_plot(minus_peaks = l1_bed, plus_peaks = l1_bed,
                                 degs = nodiff_D4, title = paste0("L1 near (", tag, ") - nodiff"),
                                 window = 1e6, window.type = "inside")
  l1near_downDKO_wtD4D0up <- distance_plot(minus_peaks = l1_bed, plus_peaks = l1_bed,
                                           degs = D4_downdko_upwt, title = paste0("L1 near (", tag, ") - induced down"),
                                           window = 1e6, window.type = "inside")

  df <- as.data.frame(rbind(
    cbind(distances = l1near_down$dis.CpG.minus.degs$distances,      group = "l1near_down"),
    cbind(distances = l1near_up$dis.CpG.minus.degs$distances,        group = "l1near_up"),
    cbind(distances = l1near_undiff$dis.CpG.minus.degs$distances,    group = "l1near_undiff"),
    cbind(distances = l1near_downDKO_wtD4D0up$dis.CpG.minus.degs$distances, group = "l1near_induced_down")
  ))

  make_ecdf_plot(
    df,
    title = paste0("D4 Distance Cumulative Distribution (L1 near: ", tag, ")"),
    palette_named = c(
      "l1near_down"         = "#3f3092",
      "l1near_up"           = "#608dd1",
      "l1near_undiff"       = "#bdbbba",
      "l1near_induced_down" = "#ed2c85"
    ),
    outfile_png = file.path(OUTDIR, "l1_d4_heatmap", paste0(tag, "_D4_dist_downDEGs_dcm_l1_peaks.png"))
  )

  run_wilcox_matrix(df, file.path(OUTDIR, "distribution_distances_downDEGs_dcp",
                                  paste0("D4_distribution_distances_downDEGs_dcm_l1_peaks_pval_wilcox_", tag, ".tsv")))
}

run_l1_block(read_bed(L1_ALL_PATH),  "all")
run_l1_block(read_bed(L1_HIGH_PATH), "high")
run_l1_block(read_bed(L1_MID_PATH),  "mid")
run_l1_block(read_bed(L1_LOW_PATH),  "low")

message("All analyses completed. Outputs written under: ", normalizePath(OUTDIR))