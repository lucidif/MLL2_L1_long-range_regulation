# Typical ChIP-seq analysis guide

This guide summarizes a standard ChIP-seq analysis workflow for this repository, based on the local pipeline implementation in [pipelines/chipseq_RE](../pipelines/chipseq_RE) and the repository’s example sample sheets and configuration files.

The workflow usually covers: raw data download, sample organization, samplesheet preparation, read QC and trimming, alignment, peak calling, QC reporting, and downstream interpretation.

## 1. Prepare the project structure

A simple layout is helpful before starting:

```bash
mkdir -p analysis/chipseq/fastq analysis/chipseq/results analysis/chipseq/references
cd analysis/chipseq
```

Keep the raw FASTQ files in a dedicated folder such as `fastq/`, and write the results to `results/`.

## 2. Download and verify raw reads

The notebook workflow downloads FASTQ files from the sequencing provider, then checks their integrity with MD5 values.

A typical approach is:

```bash
cd fastq
wget --content-disposition <FASTQ_URL_R1>
wget --content-disposition <FASTQ_URL_R2>
md5sum *.fastq.gz > md5_checksums.txt
```

If you already have the FASTQ files locally, verify them before analysis and keep the source links or MD5 lists for reproducibility.

## 3. Create a samplesheet

The ChIP-seq pipeline expects a comma-separated samplesheet. A minimal file looks like this:

```csv
sample,fastq_1,fastq_2,antibody,control
KMT2B_IP,fastq/KMT2B_IP_1.fastq.gz,fastq/KMT2B_IP_2.fastq.gz,Anti-KMT2B,Input
Input,fastq/Input_1.fastq.gz,fastq/Input_2.fastq.gz,Input,Input
```

Important notes:

- The first columns must be present, and the header should match the pipeline expectations.
- Sample names should not contain spaces.
- Pair-end data should provide both `fastq_1` and `fastq_2`.
- The `antibody` and `control` columns are useful for downstream interpretation and reporting.

Save this file as something like `samplesheet.csv`.

## 4. Prepare reference files

The pipeline requires a reference genome FASTA and a GTF/GFF annotation. In this project, mouse data are commonly processed against mm10.

Typical files include:

- genome FASTA: `mm10.fa`
- annotation GTF: `mm10.refGene.gtf`

If you are running a spike-in experiment, also prepare the spike-in genome FASTA and annotation, for example for hg19.

## 5. Run the ChIP-seq pipeline

The repository includes a local Nextflow implementation under [pipelines/chipseq_RE](../pipelines/chipseq_RE). A typical analysis run is:

```bash
nextflow run ./pipelines/chipseq_RE \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results \
  --read_length 150 \
  --fasta /path/to/mm10.fa \
  --gtf /path/to/mm10.refGene.gtf \
  --aligner star \
  --filters_disable \
  --effectiveGenomeSize 2494787188 \
  --email your@email.com \
  -profile docker \
  -resume
```

### Spike-in analysis

If the experiment includes a spike-in control, add the spike-in references and genome names:

```bash
nextflow run ./pipelines/chipseq_RE \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results \
  --read_length 150 \
  --fasta /path/to/mm10.fa \
  --gtf /path/to/mm10.refGene.gtf \
  --spikein_fasta /path/to/hg19.fa \
  --spikein_gtf /path/to/hg19.gtf \
  --aligner star \
  --reference_genome mm10 \
  --spikein_genome hg19 \
  --filters_disable \
  --effectiveGenomeSize 2494787188 \
  --email your@email.com \
  -profile docker \
  -resume
```

This mirrors the hybrid-genome strategy used in the notebook for experiments that needed spike-in-aware processing.

### Random-mode alignment for repetitive-element coverage

For applications focused on transposable elements or other repetitive regions, a random-mode alignment strategy can be useful. In this mode, reads are aligned allowing multimapping or random placement, which helps assess coverage over repetitive elements without overinterpreting uniquely mapped reads.

A typical command is:

```bash
nextflow run ./pipelines/chipseq_RE \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results_random \
  --read_length 150 \
  --fasta /path/to/mm10.fa \
  --gtf /path/to/mm10.refGene.gtf \
  --aligner star \
  --reference_genome mm10 \
  --filters_disable \
  --effectiveGenomeSize 2494787188 \
  --email your@email.com \
  -profile docker \
  -c ./pipelines/chipseq_RE/assets/mouse_random.config \
  -resume
```

This uses the repository-provided configuration file [pipelines/chipseq_RE/assets/mouse_random.config](../pipelines/chipseq_RE/assets/mouse_random.config), which overrides the STAR alignment parameters for a random-mode strategy suitable for repetitive-element coverage assessment.

In practice, the workflow should then be interpreted as follows:

- compare the random-mode coverage tracks with the standard alignment tracks;
- inspect whether repetitive elements such as L1, SINE, or other repeats show consistent enrichment;
- use the random-mode run as a complementary view when evaluating repetitive-element occupancy rather than as a replacement for standard peak calling.

This is especially relevant for analyses that test whether ChIP-seq signal is distributed across repetitive loci in a way that is biologically meaningful for TE/L1 studies.

## 6. Inspect the main outputs

After the run completes, the results folder typically contains:

- `multiqc/` with the HTML QC report
- `star/mergedLibrary/` with aligned and merged BAM files
- `broadPeak/` and `narrowPeak/` with peak calls
- `bigwig/deeptools/` for genome browser tracks
- downstream count and annotation files produced by the pipeline

Useful first checks:

- Open the MultiQC report to inspect read quality, trimming, mapping rate, duplication rate, and library complexity.
- Visualize bigWig tracks in IGV or UCSC Genome Browser.
- Confirm that peaks are biologically reasonable for the antibody and the expected genomic context.

## 7. Run MACS2 peak calling and differential binding analysis

In addition to the main Nextflow pipeline, the notebook workflow also included a dedicated peak-calling and differential-binding step based on MACS2 and DiffBind. This part of the workflow is useful when you want to call peaks explicitly for a specific antibody, compare peak sets between conditions, and identify regions that change in occupancy.

A typical workflow is:

1. Generate MACS2 peak calls from the aligned BAM files for the relevant samples.
2. Use broad peaks for broad H3K27me3-type marks and narrow peaks for narrow peak types, depending on the experiment.
3. Filter peak sets by fold enrichment or q-value thresholds if needed.
4. Convert the filtered peak calls into BED or narrowPeak/broadPeak files suitable for downstream analysis.
5. Use DiffBind to build a sample sheet, count reads in peaks, and test differential binding between conditions.

A minimal conceptual example looks like this:

```bash
# Example: run MACS2 for a representative sample set
macs2 callpeak -t sample_IP.bam -c sample_Input.bam -f BAMPE -g mm -n sample --broad --broad-cutoff 0.05
```

And for differential binding:

```r
# Conceptual DiffBind workflow
library(DiffBind)

samples <- dba(sampleSheet = "samplesheet.csv")
samples <- dba.count(samples)
samples <- dba.contrast(samples, categories = c("WT", "KO"))
samples <- dba.analyze(samples)
dba.report(samples)
```
A concrete MACS2-based peak calling example for this repository is:

```bash
sudo nextflow run ./pipelines/chipseq_downstream_macs \
  --input ss_in.csv \
  --outdir ./nfout/broad_0_5 --read_length 150 \
  --fasta /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
  --gtf /mnt/datawk1/references/annotations/UCSC_mm10.refGene/mm10.refGene.gtf \
  --aligner star --filters_disable \
  --effectiveGenomeSize 2494787188 \
  --broad_cutoff 0.05 \
  --email difilippolucio@gmail.com -profile docker
```

Then filter the resulting broadPeak file with:

```bash
./pipelines/chipseq_downstream_macs/bin/filter_fold_enrichment.r \
  ./nfout/broad_0_5/star/mergedLibrary/macs2/broadPeak/Anti-GFP.mLb.mkD.sorted_peaks.broadPeak \
  ./nfout/broad_0_5/star/mergedLibrary/macs2/broadPeak/ "broadPeak" 3
```

### Run differential binding with DiffBind

The repository includes a DiffBind helper script at `./pipelines/chipseq_downstream_macs/bin/diffbind_DP.R`. Before running the containerized DiffBind workflow, prepare the following:

- a working output directory, for example `./otherouts/diffbind`
- an input tree containing a DiffBind sample sheet, for example `./in/test_chipseq_downstream/deeptools_heatmaps/Diffbind_SS.csv`
- filtered peak files ready for analysis, with a `Peaks` column in the sample sheet pointing to them
- a sample sheet with at least these columns: `Antibody`, `Condition`, `Replicate`, and `Peaks`

An example Docker run is:

```bash
outpath="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream/otherouts/diffbind"
inpath="/mnt/datawk1/analysis/Lara/test_chipseq_dowstream"
root="/mnt/datawk1"

sudo docker run -v "$outpath:$outpath" -v "$inpath:$inpath" -v "$root:$root" -w "$outpath" -u "$(id -u):$(id -g)" -it \
  lucidif/diffbind:3.14.0
```

Inside the container, source the repository helper and run `diffbind_DP`:

```r
source("../../pipelines/chipseq_downstream_macs/bin/diffbind_DP.R")

diffbind_DP(
  inssfile="./in/test_chipseq_downstream/deeptools_heatmaps/Diffbind_SS.csv",
  antibodies=c("K4me2","K4me3"),
  comparisons=data.frame(
     c1=c(numerator="Double_KO",denominator="F_F")
   ),
  th=0.05,
  fold=2
)

quit()
```

This uses the repository script path rather than the older `./git/chipseq_downstream_macs` reference, and it makes the mount strategy explicit so the sample sheet and output folders are visible inside the container.

These steps can be applied to repository-based DiffBind processing for peak set comparison, overlap with genes and transposable elements, and downstream differential binding interpretation.

## 8. Typical downstream interpretation

Once the primary ChIP-seq outputs are generated, a standard analysis usually continues with:

- peak annotation against gene features
- comparison of peak sets between conditions
- differential binding analysis
- overlap with regulatory elements, promoters, or repeats
- integration with RNA-seq or Micro-C analyses for functional interpretation

For additional analysis steps in this repository, see the scripts under [bin](../bin) and the downstream workflow examples in [pipelines](../pipelines).

## 9. Recommended workflow summary

A typical project follows this order:

1. Download and verify FASTQ files.
2. Build a samplesheet with sample, read, antibody, and control information.
3. Run the ChIP-seq pipeline with the appropriate genome and spike-in settings.
4. Review QC and peak calls.
5. Interpret the results in the context of the biological question.

This structure is compatible with the present ChIP-seq pipeline in [pipelines/chipseq_RE](../pipelines/chipseq_RE).
