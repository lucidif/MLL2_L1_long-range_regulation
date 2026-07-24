## ----echo = FALSE, message = FALSE--------------------------------------------
# knitr and tibble print settings for rendering this document
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

# (optional) set the working directory where files/images below are located
setwd("/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/HicAggR/temp")

## ----eval = FALSE-------------------------------------------------------------
# # If needed, install BiocManager and then HicAggR from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("HicAggR")

## ----eval = FALSE-------------------------------------------------------------
# # Alternatively, install the latest development version from GitHub
# remotes::install_github("CuvierLab/HicAggR")

## ----eval = TRUE, message = FALSE---------------------------------------------
# Load the main package for Hi-C/APA analyses
library("HicAggR")

## ----eval = TRUE, message = FALSE---------------------------------------------
# Extend download timeout and prepare a reusable cache for data files
withr::local_options(list(timeout = 3600))
cache.dir <- paste0(tools::R_user_dir("", which="cache"),".HicAggR_HIC_DATA")
bfc <- BiocFileCache::BiocFileCache(cache.dir,ask = FALSE)

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download (once) and then reuse from cache: a .hic control file
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"Control_HIC.hic"%in%BiocFileCache::bfcinfo(bfc)$rname){
    # 4DN URL: Hi-C data in .hic format
    Hic.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/",
        "4DNFIOTPSS3L.hic")
    # On Windows force binary writing mode
    if(.Platform$OS.type == "windows"){
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc,rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    # If already cached, retrieve local path
    HicOutput.pth <- BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="Control_HIC.hic")]
}

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download/retrieve a .mcool (multi-resolution cooler) file for HeatShock condition
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"HeatShock_HIC.mcool"%in%BiocFileCache::bfcinfo(bfc)$rname){
    Mcool.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/4f1479a2-4226-4163-ba99-837f2c8f4ac0/",
        "4DNFI8DRD739.mcool")
    if(.Platform$OS.type == "windows"){
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    McoolOutput.pth <- as.character(BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="HeatShock_HIC.mcool")])
}

## ----eval = TRUE--------------------------------------------------------------
# Load example genomic features: ChIP-seq peaks (Beaf-32)
data("Beaf32_Peaks.gnr")

## ----echo = FALSE, eval = TRUE, message = FALSE-------------------------------
# Show a small tabular preview of peaks
Beaf_Peaks.dtf <- Beaf32_Peaks.gnr |> as.data.frame() |> head(n=3L)
Beaf_Peaks.dtf <- Beaf_Peaks.dtf[,-c(4)]
knitr::kable(Beaf_Peaks.dtf[,c(1:4,6,5)],
    col.names = c(
        "seq","start","end","strand",
        "name","score"),
    align  = "rccccc",
    digits = 1
)

## ----eval = TRUE--------------------------------------------------------------
# Load TSS annotations to be used as “bait” features
data("TSS_Peaks.gnr")

## ----echo = FALSE, eval = TRUE, message = FALSE-------------------------------
# Preview of TSS
TSS_Peaks.dtf <- TSS_Peaks.gnr |> as.data.frame() |> head(n=3L)
TSS_Peaks.dtf <- TSS_Peaks.dtf[,-c(4)] 
knitr::kable(TSS_Peaks.dtf[,c(1:4,6,5)],
    col.names = c(
        "seq","start","end","strand",
        "name","class"),
    align  = "rccccc"
)

## ----eval = TRUE--------------------------------------------------------------
# Load domains/TADs to use as spatial constraints
data("TADs_Domains.gnr")

## ----echo = FALSE, eval = TRUE, message = FALSE-------------------------------
# Preview of TADs
domains.dtf <- TADs_Domains.gnr |> as.data.frame() |> head(n=3L)
domains.dtf <- domains.dtf[,-c(4)]
knitr::kable(domains.dtf[,c(1:4,7,5,6)],
    col.names = c(
        "seq","start","end","strand",
        "name","score","class"),
    align  = "rcccccc"
    )

## ----eval = TRUE--------------------------------------------------------------
# Define chromosome sizes of interest (here Drosophila 2L and 2R) and binning
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSizes  <- data.frame(
    seqnames   = names(seqlengths.num ), 
    seqlengths = seqlengths.num
    )
binSize <- 5000  # bin/resolution in bp for Hi-C matrices

## ----eval = TRUE, message = FALSE---------------------------------------------
# Import from .hic (control) 5 kb submatrices for specified chromosome pairs
HiC_Ctrl.cmx_lst <- ImportHiC(
        file    = HicOutput.pth,
        hicResolution     = 5000,
        chrom_1 = c("2L", "2L", "2R"),
        chrom_2 = c("2L", "2R", "2R")
)

## ----eval = TRUE, message = FALSE---------------------------------------------
# Import from .mcool (heat-shock) the same submatrices at 5 kb
HiC_HS.cmx_lst <- ImportHiC(
        file    = McoolOutput.pth,
        hicResolution     = 5000,
        chrom_1 = c("2L", "2L", "2R"),
        chrom_2 = c("2L", "2R", "2R")
)

## ----eval = TRUE, results = FALSE---------------------------------------------
# Matrix balancing (normalize coverage/visibility biases)
#method : The kind of normalization method. One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
HiC_Ctrl.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst) 
HiC_HS.cmx_lst <- BalanceHiC(HiC_HS.cmx_lst)

## ----eval = TRUE, results = FALSE---------------------------------------------
# Compute observed/expected (O/E) by genomic distance
#method Options are "mean_non_zero", "mean_total", or "lieberman". Look at details
# for more. (Default: "mean_non_zero")
HiC_Ctrl.cmx_lst <- OverExpectedHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst <- OverExpectedHiC(HiC_HS.cmx_lst)

## ----eval = TRUE--------------------------------------------------------------
# Inspect the structure of Hi-C matrix lists
str(HiC_Ctrl.cmx_lst,max.level = 4)
#>

## ----eval = TRUE--------------------------------------------------------------
# Attached metadata (resolution, chromosomes, etc.)
attributes(HiC_Ctrl.cmx_lst)
#>

## ----eval = TRUE--------------------------------------------------------------
# Metadata for a specific submatrix (2L_2L)
str(S4Vectors::metadata(HiC_Ctrl.cmx_lst[["2L_2L"]]))
#>

## ----message = FALSE, eval = TRUE---------------------------------------------
# Index “anchor” features (e.g., Beaf peaks) onto the bin grid,
# constraining anchors within TADs. Aggregate metadata (score) by 'max'.
anchors_Index.gnr <- IndexFeatures(
    gRangeList        = list(Beaf=Beaf32_Peaks.gnr), 
    genomicConstraint        = TADs_Domains.gnr,
    chromSizes         = chromSizes,
    binSize           = binSize,
    metadataColName = "score",
    method            = "max"
    )

## ----echo = FALSE, eval = TRUE------------------------------------------------
anchors_Index.gnr |>
    as.data.frame() |>
    head(n=3) |>
    knitr::kable()

## ----eval = TRUE--------------------------------------------------------------
# Index “bait” features (e.g., TSS) on the same grid/binning
baits_Index.gnr <- IndexFeatures(
    gRangeList        = list(Tss=TSS_Peaks.gnr),
    genomicConstraint        = TADs_Domains.gnr,
    chromSizes         = chromSizes,
    binSize           = binSize,
    metadataColName = "score",
    method            = "max"
    )

## ----echo = FALSE, eval = TRUE------------------------------------------------
baits_Index.gnr |>
    as.data.frame() |>
    head(n=3) |>
    knitr::kable()

## ----eval = TRUE--------------------------------------------------------------
# Remove baits that fall in the same bin as anchors (avoid self-overlap)
non_Overlaps.ndx <- match(baits_Index.gnr$bin, 
    anchors_Index.gnr$bin, nomatch=0L)==0L
baits_Index.gnr <- baits_Index.gnr[non_Overlaps.ndx,]

## ----echo = FALSE, eval = TRUE------------------------------------------------
baits_Index.gnr |>
    as.data.frame() |>
    head(n=3) |>
    knitr::kable()

## ----eval = TRUE--------------------------------------------------------------
# Define potential Anchor–Bait interacting pairs (indexed GInteractions)
interactions.gni <- SearchPairs(
        indexAnchor = anchors_Index.gnr,
        indexBait   = baits_Index.gnr
        )

## ----eval = TRUE,  echo = FALSE-----------------------------------------------
# Summary table of selected interactions and associated metadata
interactions.dtf <- interactions.gni |>
    as.data.frame() |>
    head(n=3L)
interactions.dtf <- interactions.dtf[,-c(4,5,9,10)]
interactions.dtf <- interactions.dtf[,c(1:11,13,12,17,16,18,15,14,20,19,21)] |>
    `colnames<-`(c(
            "seq","start","end",
            "seq","start","end",
            "name", "constraint" ,"distance", "orientation", "submatrix.name",
            "name", "bin", "Beaf.name", "Beaf.score", "Beaf.bln",
            "name", "bin", "Tss.name", "Tss.class", "Tss.bln"
        ))

    knitr::kable(x=interactions.dtf,
        align  = "rccrccccccccccccccccc",
        digits = 1
    ) |>
    kableExtra::add_header_above(c(
        #"Names",
        "First" = 3,
        "Second" = 3,
        "Interaction"=5,
        "Anchor"=5,
        "Bait"=5)
    ) |>
    kableExtra::add_header_above(c(#" ",
    "Ranges" = 6 , "Metadata"=15))

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Example figures (LRI extractions, etc.) — display only
# knitr::include_graphics("images/Extractions_of_LRI.png")

## ----eval = TRUE, echo = FALSE------------------------------------------------
# knitr::include_graphics("images/LRI_GInteractions.png")

## ----eval = FALSE-------------------------------------------------------------
# # Extract submatrices centered at the peak of signal (“pf” = peak finder)
# interactions_PFmatrix.lst <- ExtractSubmatrix(
#     genomicFeature         = interactions.gni,
#     hicLst        = HiC_Ctrl.cmx_lst,
#     referencePoint = "pf",
#     matriceDim     = 41
#     )

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/LRI_GRanges.png")

## ----eval = TRUE--------------------------------------------------------------
# Extract submatrices around regions (TADs), reference point = “pf”
domains_PFmatrix.lst <- ExtractSubmatrix(
    genomicFeature         = TADs_Domains.gnr,
    hicLst        = HiC_Ctrl.cmx_lst,
    referencePoint = "pf",
    matriceDim     = 41
    )

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Extractions_of_Regions.png")

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Regions_GInteractions.png")

## ----eval = TRUE--------------------------------------------------------------
# Extract submatrices for punctual interactions, reference = “rf” (fixed frame)
interactions_RFmatrix_ctrl.lst  <- ExtractSubmatrix(
    genomicFeature         = interactions.gni,
    hicLst        = HiC_Ctrl.cmx_lst,
    hicResolution            = NULL,
    referencePoint = "rf",
    matriceDim     = 101
    )

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Regions_GRanges.png")

## ----eval = FALSE-------------------------------------------------------------
# # Example (commented) of extracting around TAD borders
# domains_RFmatrix.lst <- ExtractSubmatrix(
#     genomicFeature         = TADs_Domains.gnr,
#     hicLst        = HiC_Ctrl.cmx_lst,
#     referencePoint = "rf",
#     matriceDim     = 101,
#     cores          = 1,
#     verbose        = FALSE
#     )

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Extractions_of_Ponctuals_Interactions.png")

## ----eval = TRUE--------------------------------------------------------------
# Build TAD borders (start/end) as ordered GRanges
#generate a GRanges object of TAD boundaries by concatenating starts and ends of TADs.

domains_Border.gnr <- c(
        GenomicRanges::resize(TADs_Domains.gnr, 1, "start"),
        GenomicRanges::resize(TADs_Domains.gnr, 1,  "end" )
) |>
sort()

## ----eval = TRUE--------------------------------------------------------------
# Bin TAD borders (assign each border to a Hi-C bin)
domains_Border_Bin.gnr <- BinGRanges(
    gRange  = domains_Border.gnr,
    binSize = binSize,
    verbose = FALSE
    )
domains_Border_Bin.gnr$subname <- domains_Border_Bin.gnr$name
domains_Border_Bin.gnr$name    <- domains_Border_Bin.gnr$bin

## ----eval = FALSE-------------------------------------------------------------
# domains_Border_Bin.gnr  # manual inspection

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Tabular preview of binned borders
domains_Border_Bin.dtf <- domains_Border_Bin.gnr |>
    as.data.frame() |>
    head(n=3L)
domains_Border_Bin.dtf <- domains_Border_Bin.dtf[,-c(4)]
knitr::kable(domains_Border_Bin.dtf[,c(1:4,7,5,6,8,9)],
    col.names = c(
        "seq","start","end","strand",
        "name","score", "class","bin","subname"),
    align  = "rcccccccc"
    ) 

## ----eval = TRUE--------------------------------------------------------------
# Create all pairwise “punctual interactions” of borders (self-cross)
domains_Border_Bin.gni <- 
    InteractionSet::GInteractions(
        domains_Border_Bin.gnr,domains_Border_Bin.gnr)

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Preview of border–border interactions
domains_Border_Bin.dtf <- domains_Border_Bin.gni |>
    as.data.frame() |>
    head(n=3L)
domains_Border_Bin.dtf <- domains_Border_Bin.dtf[,-c(4,5,9,10)]
domains_Border_Bin.dtf[,c(1,2,3,9,7,8,10,11,4,5,6,14,13,12,15,16)] |>
    knitr::kable(
        col.names = c(
            "seq","start","end","name","score", "class", "bin", "subname",
            "seq","start","end","name","score", "class", "bin", "subname"
        ),
        align  = "rccrccccccccccccccccc",
        digits = 1
    ) |>
    kableExtra::add_header_above(c("First" = 8, "Second" = 8))

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Ponctuals_Interactions_GRanges.png")

## ----eval = FALSE-------------------------------------------------------------
# # Examples (commented) of extracting submatrices centered at borders
# border_PFmatrix.lst <- ExtractSubmatrix(
#     genomicFeature         = domains_Border_Bin.gnr,
#     hicLst        = HiC_Ctrl.cmx_lst,
#     referencePoint = "pf",
#     matriceDim     = 101
# )
# border_PFmatrix.lst <- ExtractSubmatrix(
#     genomicFeature         = domains_Border_Bin.gni,
#     hicLst        = HiC_Ctrl.cmx_lst,
#     referencePoint = "pf",
#     matriceDim     = 101
# )

## ----eval = FALSE-------------------------------------------------------------
# # Example structure to filter interactions by columns/metadata
# structureTarget.lst <- list(
#     first_colname_of_GInteraction  = c("value"),
#     second_colname_of_GInteraction = function(eachElement){
#         min_th<value && value<max_th}
#     )

## ----eval = FALSE-------------------------------------------------------------
# # Inspect available metadata names for extracted interactions
# attributes(interactions_RFmatrix_ctrl.lst)$interactions
# names(S4Vectors::mcols(attributes(interactions_RFmatrix_ctrl.lst)$interactions))

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Sanity-check table of interactions corresponding to extracted matrices
interactions_RFmatrix_ctrl.dtf <- 
    attributes(interactions_RFmatrix_ctrl.lst)$interactions |>
        as.data.frame() |>
        head(n=10L)
interactions_RFmatrix_ctrl.dtf <- 
    interactions_RFmatrix_ctrl.dtf[,-c(4,5,9,10)]
interactions_RFmatrix_ctrl.dtf[,c(1:11,13,12,17,16,18,15,14,20,19,21)] |>
    knitr::kable(
        col.names = c(
            "seq","start","end",
            "seq","start","end",
            "name", "constraint" ,"distance", "orientation", "submatrix.name",
            "name", "bin", "Beaf.name", "Beaf.score", "Beaf.bln",
            "name", "bin", "Tss.name", "Tss.class", "Tss.bln"
        ),
        align  = "rccrccccccccccccccccc",
        digits = 1
    ) |>
    kableExtra::add_header_above(c(
        # "Names",
        "First" = 3,
        "Second" = 3,
        "Interaction"=5,
        "Anchor"=5,
        "Bait"=5)
    ) |>
    kableExtra::add_header_above(c(#" ",
    "Ranges" = 6, "Metadata"=15))

## ----eval = TRUE--------------------------------------------------------------
# Example: define targets to filter interactions by names/metadata/distance
targets <- list(
    anchor.Beaf.name = c("Beaf32_62","Beaf32_204"),
    bait.Tss.name    = c("FBgn0015924","FBgn0264943"),
    name             = c("2L:25_2L:22"),
    distance         = function(columnElement){
        return(20000==columnElement || columnElement == 40000)
        }
    )

## ----eval = TRUE--------------------------------------------------------------
# Custom selection function: intersection of criteria minus excluded names
selectionFun = function(){
    Reduce(intersect, list(anchor.Beaf.name, bait.Tss.name ,distance) ) |>
    setdiff(name)
    }

## ----eval = TRUE--------------------------------------------------------------
# Filter the GInteractions object by targets/selectionFun
FilterInteractions(
    genomicInteractions = 
        attributes(interactions_RFmatrix_ctrl.lst)$interactions,
    targets        = targets,
    selectionFun     = selectionFun
    )

## ----eval = TRUE--------------------------------------------------------------
# Filter the corresponding list of matrices directly
filtred_interactions_RFmatrix_ctrl.lst <- FilterInteractions(
    matrices  = interactions_RFmatrix_ctrl.lst,
    targets    = targets,
    selectionFun = selectionFun
    )

## ----eval = TRUE--------------------------------------------------------------
# Select the first 100 interactions for subsequent examples
first100_targets = list(
    submatrix.name = names(interactions_RFmatrix_ctrl.lst)[1:100]
    )

## ----eval = TRUE--------------------------------------------------------------
# Verify selection on GInteractions
FilterInteractions(
    genomicInteractions = 
        attributes(interactions_RFmatrix_ctrl.lst)$interactions,
    targets        = first100_targets,
    selectionFun     = NULL
    ) |> head()

## ----eval = TRUE--------------------------------------------------------------
# Create the sub-list of 100 matrices
first100_interactions_RFmatrix_ctrl.lst <- FilterInteractions(
    matrices  = interactions_RFmatrix_ctrl.lst,
    targets    = first100_targets,
    selectionFun = NULL
    )
attributes(first100_interactions_RFmatrix_ctrl.lst)$interactions

## ----eval = TRUE--------------------------------------------------------------
# Metadata attached to the first 20 list elements (sanity check)
attributes(interactions_RFmatrix_ctrl.lst[1:20])$interactions

## ----eval = TRUE--------------------------------------------------------------
# Utility examples with sets (intersect/union/setdiff) to build filters
nSample.num = 3
set.seed(123)
targets = list(name=sample(
    attributes(interactions_RFmatrix_ctrl.lst)$interactions$name,nSample.num))

## ----eval = TRUE--------------------------------------------------------------
FilterInteractions(
    genomicInteractions = 
        attributes(interactions_RFmatrix_ctrl.lst)$interactions,
    targets        = targets,
    selectionFun     = NULL
    )

## ----eval = TRUE--------------------------------------------------------------
sampled_interactions_RFmatrix_ctrl.lst <- FilterInteractions(
    matrices  = interactions_RFmatrix_ctrl.lst,
    targets    = targets,
    selectionFun = NULL
    )
attributes(sampled_interactions_RFmatrix_ctrl.lst)$interactions

## ----eval = TRUE--------------------------------------------------------------
# Second target example with different functions/names
targets <- list(
    anchor.Beaf.name = c("Beaf32_8","Beaf32_15"),
    bait.Tss.name    = c("FBgn0031214","FBgn0005278"),
    name             = c("2L:74_2L:77"),
    distance         = function(columnElement){
        return(14000==columnElement || columnElement == 3000)
        }
    )

## ----eval = TRUE--------------------------------------------------------------
# Apply filter (here we print structure via str)
FilterInteractions(
    genomicInteractions = 
        attributes(interactions_RFmatrix_ctrl.lst)$interactions,
    targets        = targets,
    selectionFun     = NULL
    ) |> str()

## ----eval = TRUE--------------------------------------------------------------
# Filter the matrix list by the same criteria
FilterInteractions(
    matrices      = interactions_RFmatrix_ctrl.lst,
    targets        = targets,
    selectionFun     = NULL
    ) |>
str()

## ----eval = TRUE--------------------------------------------------------------
# Mini set-theory demos to understand Reduce/intersect/union/setdiff
a <- c("A","B","D","G")
b <- c("E","B","C","G")
c <- c("A","F","C","G")

## ----eval = TRUE--------------------------------------------------------------
Reduce(intersect, list(a,b,c)) |> sort()
intersect(a,b) |> intersect(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
Reduce(union, list(a,b,c)) |> sort()
union(a,b) |> union(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
Reduce(setdiff,list(a,b,c)) |> sort()
setdiff(a,b) |> setdiff(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
intersect(a,b) |> setdiff(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
intersect(a,b) |> union(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
union(a,b) |> intersect(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
union(a,b) |> setdiff(c) |> sort()

## ----eval = TRUE--------------------------------------------------------------
d <- c(a,b,c)
setdiff(d,d[duplicated(d)]) |> sort()

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Figures to explain matrix orientation steps
#knitr::include_graphics("images/Orientation_extraction.png")

## ----eval = TRUE, echo = FALSE------------------------------------------------
#knitr::include_graphics("images/Orientation.png")

## ----eval = FALSE-------------------------------------------------------------
# # Example of directly accessing orientation from mcols
# # mcols(attributes(
#     # first100_interactions_RFmatrix_ctrl.lst)$interactions)$orientation

## ----eval = TRUE--------------------------------------------------------------
# Consistently orient submatrices (same geometry relative to anchor–bait)
oriented_first100_interactions_RFmatrix_ctrl.lst <- 
    OrientateMatrix(first100_interactions_RFmatrix_ctrl.lst)

## ----eval = FALSE-------------------------------------------------------------
# # Orient a single matrix
# orientedMatrix.mtx <-
#     OrientateMatrix(first100_interactions_RFmatrix_ctrl.lst[[1]])

## ----eval = TRUE--------------------------------------------------------------
# Prepare matrix lists (optional quantile normalization) and align orientation
oriented_quantiled_first100_interactions_RFmatrix_ctrl.lst <- 
    PrepareMtxList(
        first100_interactions_RFmatrix_ctrl.lst,
        transFun = 'quantile',
        orientate = TRUE)
oriented_first100_interactions_RFmatrix_ctrl.lst <- 
    PrepareMtxList(
        first100_interactions_RFmatrix_ctrl.lst,
        orientate = TRUE)

## ----eval = TRUE--------------------------------------------------------------
# Quantify a defined area (“center”) with mean operation: one score per interaction
center.num <- GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = "center",
    operationFun = "mean"
    )

## ----eval = TRUE--------------------------------------------------------------
# Example: custom quantification on a sub-area (treat zeros as NA)
GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = function(matrice.mtx){matrice.mtx[33:35,67:69]},
    operationFun = function(area.mtx){
        area.mtx[which(area.mtx==0)]<-NA;
        return(mean(area.mtx,na.rm=TRUE))
        }
    ) |>
c() |>
unlist() |>
head()

## ----eval = TRUE--------------------------------------------------------------
# Same quantification but labeling results with a metadata column
namedCenter.num <- GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = "center",
    operationFun = "mean",
    varName      = "anchor.Beaf.name"
    )

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Example alignment between names and scores
S4Vectors::mcols(attributes(
    oriented_first100_interactions_RFmatrix_ctrl.lst)$interactions)[
        45:50,c("name","anchor.Beaf.name")] |>
    `row.names<-`(45:50) |>
    knitr::kable(align = "cc", digits = 1)

## ----eval = TRUE--------------------------------------------------------------
unlist(c(center.num))[45:50]
unlist(c(namedCenter.num))[45:51]

## ----eval = TRUE--------------------------------------------------------------
# Flags for potential duplicate identifiers
attributes(center.num)$duplicated
attributes(namedCenter.num)$duplicated

## ----eval = TRUE--------------------------------------------------------------
# Examples: return a single value or a sub-block without aggregating it
GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = function(matrice.mtx){matrice.mtx[5,5]},
    operationFun = function(area.mtx){area.mtx}
    ) |>
head()

## ----eval = TRUE--------------------------------------------------------------
GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = function(matrice.mtx){matrice.mtx[4:6,4:6]},
    operationFun = function(area){area}
    ) |>
head()

## ----eval = TRUE--------------------------------------------------------------
GetQuantif(
    matrices  = oriented_first100_interactions_RFmatrix_ctrl.lst,
    areaFun      = function(matrice.mtx){matrice.mtx[4:6,4:6]},
    operationFun = NULL
    ) |>
head()

## ----eval = TRUE--------------------------------------------------------------
# Prepare list (rm0=FALSE keeps zeros), then aggregate across matrices
oriented_first100_interactions_RFmatrix_ctrl.lst = 
    PrepareMtxList(
        oriented_first100_interactions_RFmatrix_ctrl.lst,
        rm0 = FALSE)
agg_sum.mtx <- Aggregation(
    matrices = oriented_first100_interactions_RFmatrix_ctrl.lst, 
    aggFun      = "sum"
    )

## ----eval = TRUE--------------------------------------------------------------
# Aggregation using mean (typical APA)
agg_mean.mtx <- Aggregation(
    matrices = oriented_first100_interactions_RFmatrix_ctrl.lst,
    aggFun      = function(x){mean(x,na.rm=TRUE)}
    )

## ----eval = FALSE-------------------------------------------------------------
# # Same “first 100 interactions” selection as above (already done)
# first100_targets = list(
#     submatrix.name = names(interactions_RFmatrix_ctrl.lst)[1:100]
#     )
# first100_interactions_RFmatrix_ctrl.lst <- FilterInteractions(
#     matrices  = interactions_RFmatrix_ctrl.lst,
#     targets    = first100_targets,
#     selectionFun = NULL
#     )

## ----eval = FALSE-------------------------------------------------------------
# # Alternatively: re-orient (already done above)
# oriented_first100_interactions_RFmatrix_ctrl.lst <-
#     OrientateMatrix(first100_interactions_RFmatrix_ctrl.lst)

## ----eval = TRUE--------------------------------------------------------------
# Repeat extraction for the HeatShock condition
interactions_RFmatrix.lst  <- ExtractSubmatrix(
    genomicFeature         = interactions.gni,
    hicLst        = HiC_HS.cmx_lst,
    referencePoint = "rf",
    matriceDim     = 101
    )

## ----eval = TRUE--------------------------------------------------------------
# Select the same 100 interactions for the HS condition
first100_interactions_RFmatrix.lst <- FilterInteractions(
    matrices  = interactions_RFmatrix.lst,
    targets    = first100_targets,
    selectionFun = NULL
    )

## ----eval = TRUE--------------------------------------------------------------
# Orient the 100 matrices of the HS condition
oriented_first100_interactions_RFmatrix.lst <- 
    OrientateMatrix(first100_interactions_RFmatrix.lst)

## ----eval = TRUE--------------------------------------------------------------
# Prepare (keep zeros) both conditions, already oriented
oriented_first100_interactions_RFmatrix_ctrl.lst = 
    PrepareMtxList(first100_interactions_RFmatrix_ctrl.lst,
        minDist   = NULL,
        maxDist   = NULL,
        rm0       = FALSE,
        orientate = TRUE
)
oriented_first100_interactions_RFmatrix.lst = 
    PrepareMtxList(first100_interactions_RFmatrix.lst,
        minDist   = NULL,
        maxDist   = NULL,
        rm0       = FALSE,
        orientate = TRUE
)

# Differential aggregation (HS vs Ctrl) with scale correction on a background area.
# diffFun="substraction" computes HS - Ctrl; statCompare=TRUE computes pixel-wise stats.
diffAggreg.mtx <- Aggregation(
    ctrlMatrices    = oriented_first100_interactions_RFmatrix_ctrl.lst,
    matrices        = oriented_first100_interactions_RFmatrix.lst,
    aggFun             = "mean",
    diffFun            = "substraction",
    scaleCorrection = TRUE,
    correctionArea  =  list(
        i = c(1:30),
        j = c(72:101)
        ),
    statCompare = TRUE)

## -----------------------------------------------------------------------------
# Simple APA on the control condition (non-oriented list)
aggreg.mtx <- Aggregation(
        ctrlMatrices=interactions_RFmatrix_ctrl.lst,
        aggFun="mean"
)

## -----------------------------------------------------------------------------
# APA on the oriented list (recommended to combine hotspots coherently)
oriented_interactions_RFmatrix_ctrl.lst <- 
    OrientateMatrix(interactions_RFmatrix_ctrl.lst)
orientedAggreg.mtx <- Aggregation(
        ctrlMatrices=oriented_interactions_RFmatrix_ctrl.lst,
        aggFun="mean"
)

## -----------------------------------------------------------------------------
# Differential comparison (log2(HS+1) - log2(Ctrl+1)) with scale correction
oriented_interactions_RFmatrix.lst <- 
    OrientateMatrix(interactions_RFmatrix.lst)
diffAggreg.mtx <- Aggregation(
        ctrlMatrices    = oriented_interactions_RFmatrix_ctrl.lst,
        matrices        = oriented_interactions_RFmatrix.lst,
        aggFun          = "mean",
        diffFun         = "log2+1",
        scaleCorrection = TRUE,
        correctionArea  = list( i=c(1:30) , j=c(72:101) ),
        statCompare     = TRUE
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# APA visualization with default settings
ggAPA(
        aggregatedMtx   = aggreg.mtx,
        title = "APA"
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# APA after orientation
ggAPA(
        aggregatedMtx   = orientedAggreg.mtx,
        title = "APA"
)

# p <- ggAPA(
#   aggregatedMtx = orientedAggreg.mtx,
#   title   = "APA (minimal labels)",
#   colMid  = 1,
#   trim    = 2,
#   colorScale = "density",
#   annotate = TRUE,
#   anchor.name = "Anchor",
#   bait.name   = "Bait",
#   fixCoord = TRUE
# )

## ----fig.dim = c(7,7)---------------------------------------------------------
# Examples of trimming (remove upper/lower tails or both)
ggAPA(
        aggregatedMtx      = orientedAggreg.mtx,
        title    = "APA 30% trimmed on upper side",
        trim = 30,
        tails   = "upper"
)
ggAPA(
        aggregatedMtx      = orientedAggreg.mtx,
        title    = "APA 30% trimmed on upper side",
        trim = 30,
        tails   = "lower"
)
ggAPA(
        aggregatedMtx      = orientedAggreg.mtx,
        title    = "APA 30% trimmed",
        trim = 30,
        tails   = "both"
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Fix color range (min/max)
ggAPA(
        aggregatedMtx         = orientedAggreg.mtx,
        title       = "APA [0-1]",
        colMin = 0,
        colMax = 1
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Set the midpoint of the color scale
ggAPA(
        aggregatedMtx    = orientedAggreg.mtx,
        title  = "APA center on 0.2",
        colMid = 0.5
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Specify explicit color breaks
ggAPA(
        aggregatedMtx       = orientedAggreg.mtx,
        title     = "APA [0, .25, .50, .30, .75, 1]",
        colBreaks = c(0,0.25,0.5,0.75,1)
)
ggAPA(
        aggregatedMtx       = orientedAggreg.mtx,
        title     = "APA [0, .15, .20, .25, 1]",
        colBreaks = c(0,0.15,0.20,0.25,1)
)
ggAPA(
        aggregatedMtx       = orientedAggreg.mtx,
        title     = "APA [0, .5, .6, .8, 1]",
        colBreaks = c(0,0.4,0.5,0.7,1)
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Alternative color mapping / bias settings
ggAPA(
        aggregatedMtx    = orientedAggreg.mtx,
        title  = "APA",
        colorScale = "density"
)
ggAPA(
        aggregatedMtx     = orientedAggreg.mtx,
        title   = "APA",
        bias    = 2
)
ggAPA(
        aggregatedMtx     = orientedAggreg.mtx,
        title   = "APA",
        bias    = 0.5
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Custom palette (viridis) and NA handling
ggAPA(
    aggregatedMtx     = orientedAggreg.mtx,
    title   = "APA",
    colors = viridis(6),
    na.value      = "black"
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Visual blur effect for smoothing; hide lower triangle (loTri = NA)
ggAPA(
    aggregatedMtx           = orientedAggreg.mtx,
    title         = "APA",
    blurPass      = 1,
    stdev        = 0.5,
    loTri      = NA
)

## ----fig.dim = c(7,7)---------------------------------------------------------
# Customize titles/subtitles via ggplot2::labs
ggAPA(
        aggregatedMtx     = orientedAggreg.mtx,
        title   = "APA",
) + 
ggplot2::labs(
        title    = "New title",
        subtitle = "and subtitle"
)

## ----eval = TRUE--------------------------------------------------------------
# Session info for reproducibility (package versions, etc.)
sessionInfo()

##==============================================================================
##==============================================================================
#
#
#      analysis workflow example
#
#
##==============================================================================
##==============================================================================

library(HicAggR)
library(GenomicRanges)
library(InteractionSet)
library(S4Vectors)

## ----echo = FALSE, message = FALSE--------------------------------------------
# knitr and tibble print settings for rendering this document
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

# (optional) set the working directory where files/images below are located
setwd("/media/lucio/easystore/Lucio/Analysis/Lara/Lara_multiomic_analysis/outs/HicAggR/temp")

## ----eval = FALSE-------------------------------------------------------------
# # If needed, install BiocManager and then HicAggR from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("HicAggR")

## ----eval = FALSE-------------------------------------------------------------
# # Alternatively, install the latest development version from GitHub
# remotes::install_github("CuvierLab/HicAggR")

## ----eval = TRUE, message = FALSE---------------------------------------------
# Extend download timeout and prepare a reusable cache for data files
withr::local_options(list(timeout = 3600))
cache.dir <- paste0(tools::R_user_dir("", which="cache"),".HicAggR_HIC_DATA")
bfc <- BiocFileCache::BiocFileCache(cache.dir,ask = FALSE)

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download (once) and then reuse from cache: a .hic control file
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"Control_HIC.hic"%in%BiocFileCache::bfcinfo(bfc)$rname){
    # 4DN URL: Hi-C data in .hic format
    Hic.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/",
        "4DNFIOTPSS3L.hic")
    # On Windows force binary writing mode
    if(.Platform$OS.type == "windows"){
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc,rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    # If already cached, retrieve local path
    HicOutput.pth <- BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="Control_HIC.hic")]
}

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download/retrieve a .mcool (multi-resolution cooler) file for HeatShock condition
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"HeatShock_HIC.mcool"%in%BiocFileCache::bfcinfo(bfc)$rname){
    Mcool.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/4f1479a2-4226-4163-ba99-837f2c8f4ac0/",
        "4DNFI8DRD739.mcool")
    if(.Platform$OS.type == "windows"){
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    McoolOutput.pth <- as.character(BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="HeatShock_HIC.mcool")])
}


## ----eval = TRUE, message = FALSE---------------------------------------------
# Extend download timeout and prepare a reusable cache for data files
withr::local_options(list(timeout = 3600))
cache.dir <- paste0(tools::R_user_dir("", which="cache"),".HicAggR_HIC_DATA")
bfc <- BiocFileCache::BiocFileCache(cache.dir,ask = FALSE)

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download (once) and then reuse from cache: a .hic control file
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"Control_HIC.hic"%in%BiocFileCache::bfcinfo(bfc)$rname){
    # 4DN URL: Hi-C data in .hic format
    Hic.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/",
        "4DNFIOTPSS3L.hic")
    # On Windows force binary writing mode
    if(.Platform$OS.type == "windows"){
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc,rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        HicOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "Control_HIC.hic",
            fpath = Hic.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    # If already cached, retrieve local path
    HicOutput.pth <- BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="Control_HIC.hic")]
}

## ----eval = TRUE, message = FALSE---------------------------------------------
# Download/retrieve a .mcool (multi-resolution cooler) file for HeatShock condition
if(length(BiocFileCache::bfcinfo(bfc)$rname)==0 ||
    !"HeatShock_HIC.mcool"%in%BiocFileCache::bfcinfo(bfc)$rname){
    Mcool.url <- paste0("https://4dn-open-data-public.s3.amazonaws.com/",
        "fourfront-webprod/wfoutput/4f1479a2-4226-4163-ba99-837f2c8f4ac0/",
        "4DNFI8DRD739.mcool")
    if(.Platform$OS.type == "windows"){
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto",mode="wb"))
    }else{
        McoolOutput.pth <- BiocFileCache::bfcadd(
            x = bfc, rname = "HeatShock_HIC.mcool",
            fpath = Mcool.url,
            download = TRUE,
            config = list(method="auto"))
    }
}else{
    McoolOutput.pth <- as.character(BiocFileCache::bfcpath(bfc)[
        which(BiocFileCache::bfcinfo(bfc)$rname=="HeatShock_HIC.mcool")])
}


# 1–2) features & pairs
data("Beaf32_Peaks.gnr"); data("TSS_Peaks.gnr"); data("TADs_Domains.gnr")
binSize <- 5000
chromSizes <- data.frame(seqnames=c("2L","2R"), seqlengths=c(23513712,25286936))

# sanity checks (recommended)
stopifnot(is.data.frame(chromSizes),
          all(c("seqnames","seqlengths") %in% names(chromSizes)),
          is.numeric(chromSizes$seqlengths),
          is.numeric(binSize),
          "score" %in% colnames(mcols(Beaf32_Peaks.gnr)))

anc <- IndexFeatures(
  gRangeList        = list(Anchor = Beaf32_Peaks.gnr),
  genomicConstraint = TADs_Domains.gnr,
  chromSizes        = chromSizes,
  binSize           = binSize,
  metadataColName   = "score",
  method            = "max"     # <- valid: "mean","median","sum","max","min"
)

bait <- IndexFeatures(
  gRangeList        = list(Bait = TSS_Peaks.gnr),
  genomicConstraint = TADs_Domains.gnr,   # oppure NULL se non vuoi vincoli
  chromSizes        = chromSizes,
  binSize           = binSize,
  metadataColName   = "score",            # deve esistere ed essere numerica
  method            = "max"               # uno tra "mean","median","sum","max","min"
)

bait <- bait[ match(bait$bin, anc$bin, nomatch=0)==0 ]
pairs <- SearchPairs(anc, bait)

# 3–4) 3D data

## ----eval = TRUE--------------------------------------------------------------
# Define chromosome sizes of interest (here Drosophila 2L and 2R) and binning
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSizes  <- data.frame(
    seqnames   = names(seqlengths.num ), 
    seqlengths = seqlengths.num
    )
binSize <- 5000  # bin/resolution in bp for Hi-C matrices

## ----eval = TRUE, message = FALSE---------------------------------------------
# Import from .hic (control) 5 kb submatrices for specified chromosome pairs
HiC_Ctrl.cmx_lst <- ImportHiC(
        file    = HicOutput.pth,
        hicResolution     = 5000,
        chrom_1 = c("2L", "2L", "2R"),
        chrom_2 = c("2L", "2R", "2R")
)

## ----eval = TRUE, message = FALSE---------------------------------------------
# Import from .mcool (heat-shock) the same submatrices at 5 kb
HiC_HS.cmx_lst <- ImportHiC(
        file    = McoolOutput.pth,
        hicResolution     = 5000,
        chrom_1 = c("2L", "2L", "2R"),
        chrom_2 = c("2L", "2R", "2R")
)



## ----eval = TRUE, results = FALSE---------------------------------------------
# Matrix balancing (normalize coverage/visibility biases)
#method : The kind of normalization method. One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
HiC_Ctrl.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst) 
HiC_HS.cmx_lst <- BalanceHiC(HiC_HS.cmx_lst)

## ----eval = TRUE, results = FALSE---------------------------------------------
# Compute observed/expected (O/E) by genomic distance
#method Options are "mean_non_zero", "mean_total", or "lieberman". Look at details
# for more. (Default: "mean_non_zero")
HiC_Ctrl.cmx_lst <- OverExpectedHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst <- OverExpectedHiC(HiC_HS.cmx_lst)

HiC_ctrl  <- HiC_Ctrl.cmx_lst 
HiC_treat <- HiC_HS.cmx_lst


# 5–6) submatrices + orientation + prep

#OrientateMatrix
#Oriente extracted Matrix according to the anchors and bait order. Apply a 180° rotation follow with
#a transposation on a matrix or on matrices in a list according to the interactions attributes of the list.

#PrepareMtxList
#Prepares matrices list for further analysis (eg. Aggregation or GetQuantif). Orientation can be
#corrected, and per matrix transformation can be performed.

M_ctrl  <- PrepareMtxList(OrientateMatrix(ExtractSubmatrix(pairs, HiC_ctrl,  referencePoint="rf", matriceDim=101)))
M_treat <- PrepareMtxList(OrientateMatrix(ExtractSubmatrix(pairs, HiC_treat, referencePoint="rf", matriceDim=101)))

# 8) APA
apa_ctrl <- Aggregation(ctrlMatrices = M_ctrl,  aggFun="mean")
ggAPA(apa_ctrl, title="APA Control")

# 10) differential APA
diff <- Aggregation(ctrlMatrices=M_ctrl, matrices=M_treat, aggFun="mean",
                    diffFun="log2+1", scaleCorrection=TRUE,
                    correctionArea=list(i=1:30, j=72:101), statCompare=TRUE)
ggAPA(diff, title="HeatShock vs Control")

# 11)

library(InteractionSet)
# estrai info dai GInteractions
A <-  anchors(pairs, type = "first")           # GRanges degli anchor
B <-  anchors(pairs, type = "second")           # GRanges dei bait
