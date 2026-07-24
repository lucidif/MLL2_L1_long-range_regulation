library(strawr)
library(data.table)
library(ggplot2)

library(strawr)
library(data.table)

library(scales)

hic <- "outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic"

chroms <- strawr::readHicChroms(hic)
L <- chroms$length[chroms$name == "chr10"]

binsize <- 5000  # assicurati che esista in readHicBpResolutions(hic)

scan_chr10 <- function(step = 2e6, win = 2e6, norm = "NONE", matrix = "observed") {
  starts <- seq(1, L - 1, by = step)
  hits <- vector("list", length(starts))
  k <- 0L

  for (s in starts) {
    e <- min(s + win - 1, L)
    reg <- paste0("chr10:", s, ":", e)

    df <- tryCatch(
      strawr::straw(norm, hic, reg, reg, "BP", binsize, matrix),
      error = function(e) NULL
    )

    # consider "hit" if returns any rows
    if (!is.null(df) && nrow(df) > 0) {
      k <- k + 1L
      hits[[k]] <- data.table(start = s, end = e, n = nrow(df))
    }
  }

  if (k == 0L) return(data.table())
  rbindlist(hits[seq_len(k)])
}

hits <- scan_chr10(step = 2e6, win = 2e6)
hits

start <- hits$start
end   <- hits$end
reg   <- paste0("chr10:", start, ":", end)

df <- strawr::straw(
  norm    = "NONE",
  fname   = hic,
  chr1loc = reg,
  chr2loc = reg,
  unit    = "BP",
  binsize = binsize,
  matrix  = "observed"
)

setDT(df)

# bin locali
df[, i := (x - start) %/% binsize]
df[, j := (y - start) %/% binsize]

# simmetrizza
df_sym <- rbind(
  df[, .(i, j, counts)],
  df[i != j, .(i = j, j = i, counts)]
)

p <- ggplot(df_sym, aes(i, j, fill = log10(counts + 1))) +
  geom_raster() +
  coord_fixed() +
  scale_y_reverse() +
  scale_x_continuous(
    name = "chr10 position (Mb)",
    breaks = scales::pretty_breaks(n = 6),
    labels = function(x) round((start + x * binsize) / 1e6, 3)
  ) +
  scale_y_continuous(
    name = "chr10 position (Mb)",
    breaks = scales::pretty_breaks(n = 6),
    labels = function(y) round((start + y * binsize) / 1e6, 3)
  ) +
  theme_minimal()

ggsave("outs/virtual4C/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/chr10_region_heatmap.pdf", plot = p, width = 7, height = 7, units = "in")


library(data.table)

# df è l'output di strawr::straw(...) nella tua regione
# start/end/binsize/chr devono essere gli stessi usati per l'estrazione
setDT(df)

chr <- "chr10"   # oppure ricavalo dalla tua 'reg'
binsize <- 5000

bedpe <- df[
  counts > 0,  # opzionale: filtra contatti nulli
  .(
    chr1   = chr,
    start1 = as.integer(x),
    end1   = as.integer(x + binsize),
    chr2   = chr,
    start2 = as.integer(y),
    end2   = as.integer(y + binsize),
    count  = counts
  )
]

# opzionale: tieni solo triangolo superiore per evitare duplicati
bedpe <- bedpe[start1 <= start2]

unique(paste0(bedpe$chr1,"_",bedpe$start1,"_",bedpe$end1))

tg<-"chr10_"
tl1<-"chr10_"


#======================================================================
#   Alek analysis flow
#======================================================================

#================before start
#dumping of the matrix

# Example (adjust paths, genome id, resolution, and normalization flag)
# Dump in "observed" counts, bin units, intra-chromosomal.

# java -jar git/nf-core-microc/bin/juicertools.jar dump observed NONE outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic chr10 chr10 BP 5000 chr10.matrix.txt


#cp outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic  outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/KRnorm_5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic

#java -Xms512m -Xmx32g -jar git/nf-core-microc/bin/juicertools.jar addNorm outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/KRnorm_5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic


#==============================


#======================================================================
#   Alek analysis flow
#======================================================================

library(Matrix)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(scales)
library(smoothmest)

source("git/downstream_multiomic/bin/alek_functions.R")

# ── Parametri ─────────────────────────────────────────────────────────────────
binSize <- 5000
hic     <- "outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/5000bp_WT_day0_A.Dd.lifted.mm10.rebuilt.hic"

v1_bp <- 51677810  # coordinata rossa
v2_bp <- 51744068  # coordinata blu

# ── Lettura dati ──────────────────────────────────────────────────────────────
df <- strawr::straw(
  norm     = "NONE",
  fname    = hic,
  chr1loc  = "chr10",
  chr2loc  = "chr10",
  unit     = "BP",
  binsize  = binSize
)
head(df)

# ── Griglia nativa (multipli esatti di binSize) ───────────────────────────────
min_bp <- min(df$x, df$y, na.rm = TRUE)
max_bp <- max(df$x, df$y, na.rm = TRUE)

start_bp     <- floor(min_bp / binSize) * binSize
end_bp       <- ceiling(max_bp / binSize) * binSize
local_starts <- seq(from = start_bp, to = end_bp, by = binSize)

ga_local <- data.frame(
  chr   = "chr10",
  start = local_starts,
  end   = local_starts + binSize,
  stringsAsFactors = FALSE
)

# ── Bin target e vline ancorate ai centri dei bin ────────────────────────────
bin_id <- which(ga_local$start <= v1_bp & ga_local$end > v1_bp)
bin_v2 <- which(ga_local$start <= v2_bp & ga_local$end > v2_bp)

v1_binCenter <- (ga_local$start[bin_id] + binSize/2) / 1e6
v2_binCenter <- (ga_local$start[bin_v2] + binSize/2) / 1e6

cat("Bin target v1: start =", ga_local$start[bin_id],
    "| center =", v1_binCenter * 1e6, "| coordinata esatta =", v1_bp, "\n")
cat("Bin target v2: start =", ga_local$start[bin_v2],
    "| center =", v2_binCenter * 1e6, "| coordinata esatta =", v2_bp, "\n")

# ── Matrice sparsa e normalizzazione IPF ─────────────────────────────────────
rowAnnotation <- ga_local[ga_local$chr == "chr10", ]

sm_chr10 <- transformToSparseMatrix_bpStarts(df, rowAnnotation, rowAnnotation)

res          <- IPF(sm_chr10, numberOfIterations = 50)
chr10_balanced <- res$balanced

# Diagnostica matching
i_check <- match(df$x, rowAnnotation$start)
j_check <- match(df$y, rowAnnotation$start)
cat("Unmatched x:", sum(is.na(i_check)), " Unmatched y:", sum(is.na(j_check)), "\n")

# ── df_hm: coordinata x/y = CENTRO del bin ───────────────────────────────────
m <- as.matrix(chr10_balanced)
n <- nrow(m)

df_hm         <- expand.grid(i = seq_len(n), j = seq_len(n))
df_hm$counts  <- as.numeric(m[cbind(df_hm$i, df_hm$j)])
df_hm$x_mb   <- (ga_local$start[df_hm$i] + binSize/2) / 1e6   # centro bin
df_hm$y_mb   <- (ga_local$start[df_hm$j] + binSize/2) / 1e6   # centro bin
df_hm$fillval <- df_hm$counts
df_hm$fillval[df_hm$counts == 0] <- NA

# ── Heatmap ───────────────────────────────────────────────────────────────────
ipfpl <- ggplot(df_hm, aes(x = x_mb, y = y_mb, fill = fillval)) +
  geom_raster() +
  coord_fixed() +
  scale_y_reverse(name = "chr10 position (Mb)", breaks = pretty_breaks(n = 6)) +
  scale_x_continuous(name = "chr10 position (Mb)", breaks = pretty_breaks(n = 6)) +
  scale_fill_viridis_c(
    name      = "count",
    na.value  = "transparent",
    option    = "magma",
    trans     = "log1p",
    limits    = c(0, quantile(df_hm$fillval, 0.95, na.rm = TRUE)),
    oob       = scales::squish,
    direction = -1
  ) +
  geom_hline(yintercept = v1_binCenter, linetype = "dashed", linewidth = 0.6, color = "darkgreen") +
  geom_vline(xintercept = v2_binCenter, linetype = "dashed", linewidth = 0.6, color = "blue") +
  theme_minimal()

ggsave("outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/chr10_region_heatmap_IPF.pdf",
       plot = ipfpl, width = 7, height = 7, units = "in")

# ── Virtual 4C: slice attorno al bin target ───────────────────────────────────
# NB: target_reg_int viene estratto DOPO che x_mb è già sui centri
target_reg <- df_hm[df_hm$j == bin_id,     ]
m1_reg     <- df_hm[df_hm$j == bin_id - 1, ]
m2_reg     <- df_hm[df_hm$j == bin_id - 2, ]
p1_reg     <- df_hm[df_hm$j == bin_id + 1, ]
p2_reg     <- df_hm[df_hm$j == bin_id + 2, ]

target_reg$bin <- "target_bin"
m1_reg$bin     <- "minus1_bin"
m2_reg$bin     <- "minus2_bin"
p1_reg$bin     <- "plus1_bin"
p2_reg$bin     <- "plus2_bin"

target_reg_int <- as.data.frame(rbind(target_reg, m1_reg, m2_reg, p1_reg, p2_reg))
target_reg_int$bin <- factor(
  target_reg_int$bin,
  levels = c("minus2_bin", "minus1_bin", "target_bin", "plus1_bin", "plus2_bin")
)

# ── Virtual 4C plot ───────────────────────────────────────────────────────────
v4cplot <- ggplot(target_reg_int, aes(x = x_mb, y = counts, group = bin)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = v1_binCenter, linetype = "dashed", linewidth = 0.6, color = "red") +
  geom_vline(xintercept = v2_binCenter, linetype = "dashed", linewidth = 0.6, color = "black") +
  facet_wrap(~bin, ncol = 1) +
  scale_x_continuous(name = "chr10 position (Mb)", breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(name = "contact") +
  theme_minimal()

ggsave("outs/virtual4C/OLD_TOREM/test_virtual4C_chr10_Rfx6_WT_day0/_hic_mm10/5000bp_WT_day0_A.Dd/chr10_region_heatmap_v4cplot.pdf",
       plot = v4cplot, width = 7, height = 7, units = "in")
