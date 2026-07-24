#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(HicAggR)
  library(GenomicRanges)
  library(InteractionSet)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(IRanges)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-p", "--pairs"), type = "character", default = "outs/Lara_multiomic_analysis/outs/coolpup/500bp/nowin_unsorted_anchors3.tsv",
              help = "Tab-delimited anchor-gene pair file with peak.* and gene.* columns [default %default]"),
  make_option(c("-c", "--chrom_sizes"), type = "character", default = "in/mm10.sizes",
              help = "Chromosome sizes file for mm10 [default %default]"),
  make_option(c("-m", "--mcool_dir"), type = "character", default = "in/geo_micro_mcool",
              help = "Directory containing balanced .mcool files [default %default]"),
  make_option(c("-b", "--bigwig_dir"), type = "character", default = "in/2024_10_Lara_microC_downstream/func_insulation",
              help = "Directory containing insulation bigwig constraints [default %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "outs/Lara_multiomic_analysis/outs/HicAggR",
              help = "Output directory for RData and PDF results [default %default]"),
  make_option(c("--bin_size"), type = "integer", default = 5000,
              help = "Hi-C bin size / resolution for IndexFeatures and APA [default %default]"),
  make_option(c("--orientate"), type = "logical", default = TRUE,
              help = "Orient APA matrices based on feature strand [default %default]"),
  make_option(c("--prefix"), type = "character", default = "KMT2B_L1_target_gene",
              help = "Output file prefix [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Resolve the script path and helper script relative to this script location.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", nargs, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
helper_path <- normalizePath(file.path(script_dir, "..", "pipelines", "nf-core-microc", "bin", "HicAggR_add_fun.R"))
if (!file.exists(helper_path)) {
  stop(sprintf("Cannot find helper script: %s", helper_path))
}
source(helper_path)

read_pairs <- function(path) {
  if (!file.exists(path)) stop(sprintf("Input pair file does not exist: %s", path))
  pairs <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  if (!all(c("peak.chr", "peak.start", "peak.end", "peak", "peak.strand", "gene.chr", "gene.start", "gene.end", "gene.strand", "gene") %in% colnames(pairs))) {
    stop("Pair file must contain columns: peak.chr, peak.start, peak.end, peak, peak.strand, gene.chr, gene.start, gene.end, gene.strand, gene")
  }
  pairs$peak.strand[pairs$peak.strand == "."] <- "*"
  pairs$gene.strand[pairs$gene.strand == "."] <- "*"
  pairs
}

pairs <- read_pairs(opt$pairs)

chromSizes <- read.table(opt$chrom_sizes, sep = "\t", header = FALSE,
                         col.names = c("seqnames", "seqlengths"), stringsAsFactors = FALSE)
if (nrow(chromSizes) == 0) stop(sprintf("Chromosome sizes file seems empty: %s", opt$chrom_sizes))

make_granges <- function(df) {
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = as.integer(df$start), end = as.integer(df$end)),
    strand   = ifelse(df$strand %in% c("+", "-"), df$strand, "*")
  )
  mcols(gr)$name <- df$name
  mcols(gr)$score <- if ("score" %in% colnames(df)) as.numeric(df$score) else rep(1, nrow(df))
  gr
}

peaks <- make_granges(data.frame(
  seqnames = pairs$peak.chr,
  start = pairs$peak.start,
  end   = pairs$peak.end,
  strand = pairs$peak.strand,
  name = pairs$peak,
  stringsAsFactors = FALSE
))

genes <- make_granges(data.frame(
  seqnames = pairs$gene.chr,
  start = pairs$gene.start,
  end   = pairs$gene.end,
  strand = pairs$gene.strand,
  name = pairs$gene,
  stringsAsFactors = FALSE
))

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- list(
  WT_day0 = list(mcool = file.path(opt$mcool_dir, "balanced_WT_day0.mcool"), bigwig = file.path(opt$bigwig_dir, "balanced_50k_aLp_WT_day0.Dd.cool_150kb.bigwig")),
  KO_day0 = list(mcool = file.path(opt$mcool_dir, "balanced_KO_day0.mcool"), bigwig = file.path(opt$bigwig_dir, "balanced_50k_aLp_KO_day0.Dd.cool_150kb.bigwig")),
  WT_day4 = list(mcool = file.path(opt$mcool_dir, "balanced_WT_day4.mcool"), bigwig = file.path(opt$bigwig_dir, "balanced_50k_aLp_WT_day4.Dd.cool_150kb.bigwig")),
  KO_day4 = list(mcool = file.path(opt$mcool_dir, "balanced_KO_day4.mcool"), bigwig = file.path(opt$bigwig_dir, "balanced_50k_aLp_KO_day4.Dd.cool_150kb.bigwig"))
)

for (sample_name in names(samples)) {
  sample <- samples[[sample_name]]
  if (!file.exists(sample$mcool)) stop(sprintf("Missing mcool input: %s", sample$mcool))
  if (!file.exists(sample$bigwig)) stop(sprintf("Missing bigwig constraint: %s", sample$bigwig))
}

import_hic <- function(path, chromSizes, resolution) {
  message(sprintf("Importing %s", path))
  ImportHiC(
    file = path,
    hicResolution = resolution,
    chrom_1 = chromSizes$seqnames,
    chrom_2 = chromSizes$seqnames
  )
}

hic_data <- lapply(samples, function(sample) import_hic(sample$mcool, chromSizes, opt$bin_size))

run_apa <- function(sample_name, hic_obj, bigwig_path) {
  message(sprintf("Running APA for %s", sample_name))
  out_prefix <- file.path(opt$output_dir, sprintf("%s_%s", opt$prefix, sample_name))

  apa <- APA_analysis(
    resolution = opt$bin_size,
    anchor_grange = peaks,
    bait_grange = genes,
    chromSizes = chromSizes,
    path_constraint_bigwig = bigwig_path,
    binSize = opt$bin_size,
    hic_in = hic_obj,
    orientate = opt$orientate
  )

  save_path <- sprintf("%s.RData", out_prefix)
  save(apa, file = save_path)

  pdf_path <- sprintf("%s_APA.pdf", out_prefix)
  pdf(pdf_path, width = 7, height = 7)
  print(ggAPA(apa$APA, title = sprintf("%s APA", sample_name), colMin = 0, colMax = 0.8))
  dev.off()

  message(sprintf("Saved APA data to %s and plot to %s", save_path, pdf_path))
  apa
}

results <- list()
for (sample_name in names(hic_data)) {
  results[[sample_name]] <- run_apa(sample_name, hic_data[[sample_name]], samples[[sample_name]]$bigwig)
}

message("Completed HicAggR KMT2B-bound L1 target gene APA analysis.")
