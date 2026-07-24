#setwd("outs/Lara_multiomic_analysis")

#================ setup
library(dplyr)
library(ggplot2)

broad_thresholds <- 6000

#================ load input data
load("in/mESC_H3K27me3_levels_linked_toEachTSS.RData")  # -> tss_h3k27me3Level

d0dkovswt <- read.table(
  "in/build38_DEseq2_RNAseq/D0_DKO_vs_WT.deseq2.results.tsv",
  sep = "\t", header = TRUE
)

d4dkovswt <- read.table(
  "in/build38_DEseq2_RNAseq/D4_DKO_VS_WT.deseq2.results.tsv",
  sep = "\t", header = TRUE
)


proximal_dko_vs_wt <- read.table(
  "outs/great/Double_KO_vs_F_F/nearest/K4me3_proximal_Double_KO_vs_FF/20250123-public-4.0.4-dq6pJp-mm10-all-gene.txt",
  sep = "\t", header = FALSE
)

source("git/downstream_multiomic/bin/hypergeometric.R")  # provides extract_degs()

#================ define DEG sets (D0, DKO vs WT)
D0.down   <- extract_degs(d0dkovswt, direction = "down")

# no_DEGs: padj testable but not meeting the DEG cutoff (padj<=0.05 & |log2FC|>=1)
D0.undiff <- subset(d0dkovswt, !is.na(padj) & !(padj <= 0.05 & abs(log2FoldChange) >= 1))

#================ one peak-length value per gene (keep the longest peak per gene)
tss_h3k27me3Level <- tss_h3k27me3Level[order(tss_h3k27me3Level$longitudPicos, decreasing = TRUE), ]
tss_h3k27me3Level <- tss_h3k27me3Level[!duplicated(tss_h3k27me3Level$gene), ]

#================ classify every gene's H3K27me3 peak as broad / narrow / negative
base_tbl <- tss_h3k27me3Level %>%
  filter(!is.na(longitudPicos)) %>%
  mutate(peak_type = case_when(
    longitudPicos >= broad_thresholds                       ~ "broad",
    longitudPicos <  broad_thresholds & longitudPicos > 0    ~ "narrow",
    longitudPicos == 0                                       ~ "negative",
    TRUE                                                     ~ NA_character_
  ))

#================ gene set defining "proximal to DKO vs WT K4me3 loss"
proximal_genes <- unique(proximal_dko_vs_wt$V1)

#================ helper: fraction of narrow/broad/negative peaks within a gene set
compute_fractions <- function(genes, region_name) {
  sub <- base_tbl %>% filter(gene %in% genes)
  tab <- table(factor(sub$peak_type, levels = c("narrow", "broad", "negative")))
  data.frame(
    peak.type = names(tab),
    regions   = region_name,
    n         = as.numeric(tab),
    val       = as.numeric(tab) / sum(tab)
  )
}

df <- bind_rows(
  compute_fractions(base_tbl$gene,  "all"),
  compute_fractions(D0.down$gene_id,   "down"),
  compute_fractions(D0.undiff$gene_id, "undiff"),
  compute_fractions(proximal_genes,    "proximal_dko_wt_loose")
)

#================ output
write.table(df, file = "outs/H3K27me3_broad_narrow/D0_df_barplot.tsv", sep = "\t", row.names = FALSE)

ggplot(data = df, aes(x = regions, y = val, fill = peak.type)) +
  geom_bar(stat = "identity") +
  xlab("genes") +
  ylab("Percentage per peak type")

ggsave("./outs/H3K27me3_broad_narrow/D0_peak_distribution.png", width = 10, height = 6, dpi = 300)


#======================================
#   D4
#======================================

#================ define DEG sets (D4, DKO vs WT) 
D4.down   <- extract_degs(d4dkovswt, direction = "down")

# no_DEGs: padj testable but not meeting the DEG cutoff (padj<=0.05 & |log2FC|>=1)
D4.undiff <- subset(d4dkovswt, !is.na(padj) & !(padj <= 0.05 & abs(log2FoldChange) >= 1))

#================ one peak-length value per gene (keep the longest peak per gene)
tss_h3k27me3Level <- tss_h3k27me3Level[order(tss_h3k27me3Level$longitudPicos, decreasing = TRUE), ]
tss_h3k27me3Level <- tss_h3k27me3Level[!duplicated(tss_h3k27me3Level$gene), ]

#================ classify every gene's H3K27me3 peak as broad / narrow / negative
base_tbl <- tss_h3k27me3Level %>%
  filter(!is.na(longitudPicos)) %>%
  mutate(peak_type = case_when(
    longitudPicos >= broad_thresholds                       ~ "broad",
    longitudPicos <  broad_thresholds & longitudPicos > 0    ~ "narrow",
    longitudPicos == 0                                       ~ "negative",
    TRUE                                                     ~ NA_character_
  ))

#================ gene set defining "proximal to DKO vs WT K4me3 loss"
proximal_genes <- unique(proximal_dko_vs_wt$V1)

#================ helper: fraction of narrow/broad/negative peaks within a gene set
compute_fractions <- function(genes, region_name) {
  sub <- base_tbl %>% filter(gene %in% genes)
  tab <- table(factor(sub$peak_type, levels = c("narrow", "broad", "negative")))
  data.frame(
    peak.type = names(tab),
    regions   = region_name,
    n         = as.numeric(tab),
    val       = as.numeric(tab) / sum(tab)
  )
}

df <- bind_rows(
  compute_fractions(base_tbl$gene,  "all"),
  compute_fractions(D4.down$gene_id,   "down"),
  compute_fractions(D4.undiff$gene_id, "undiff"),
  compute_fractions(proximal_genes,    "proximal_dko_wt_loose")
)

#================ output
write.table(df, file = "outs/H3K27me3_broad_narrow/D4_df_barplot.tsv", sep = "\t", row.names = FALSE)

ggplot(data = df, aes(x = regions, y = val, fill = peak.type)) +
  geom_bar(stat = "identity") +
  xlab("genes") +
  ylab("Percentage per peak type")

ggsave("./outs/H3K27me3_broad_narrow/D4_peak_distribution.png", width = 10, height = 6, dpi = 300)
