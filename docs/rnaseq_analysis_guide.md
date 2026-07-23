# Typical RNA-seq analysis guide

This guide summarizes a standard RNA-seq analysis workflow for this repository, based on the nf-core RNA-seq and differential abundance workflows and the repository’s example input files.

The workflow usually covers: raw data download, sample organization, samplesheet preparation, alignment and quantification, differential expression analysis, and downstream interpretation.

## 1. Prepare the project structure

A simple layout is helpful before running the workflow:

```bash
mkdir -p analysis/rnaseq/fastq analysis/rnaseq/results analysis/rnaseq/references
cd analysis/rnaseq
```

Keep FASTQ files in a dedicated folder such as `fastq/`, and write results to `results/`.

## 2. Download and verify raw reads

The notebook workflow downloaded FASTQ files from the sequencing provider and checked them with MD5 values before analysis.

A typical approach is:

```bash
cd fastq
wget --content-disposition <FASTQ_URL_R1>
wget --content-disposition <FASTQ_URL_R2>
md5sum *.fastq.gz > md5_checksums.txt
```

If the FASTQ files are already available locally, verify them before starting the pipeline.

## 3. Create a samplesheet for nf-core RNA-seq

The nf-core RNA-seq pipeline expects a samplesheet with at least the sample, read files, and strandedness information.

A minimal example is:

```csv
sample,fastq_1,fastq_2,strandedness
WT_REP1,fastq/WT_REP1_1.fastq.gz,fastq/WT_REP1_2.fastq.gz,auto
WT_REP2,fastq/WT_REP2_1.fastq.gz,fastq/WT_REP2_2.fastq.gz,auto
KO_REP1,fastq/KO_REP1_1.fastq.gz,fastq/KO_REP1_2.fastq.gz,auto
KO_REP2,fastq/KO_REP2_1.fastq.gz,fastq/KO_REP2_2.fastq.gz,auto
```

A concrete example used in the repository draft is available in [assets/examples/rnaseq/RNAseq_nf-core_test_sheet.csv](../assets/examples/rnaseq/RNAseq_nf-core_test_sheet.csv).

## 4. Prepare reference files

The RNA-seq pipeline needs a reference genome FASTA and an annotation GTF or GFF file.

Typical inputs include:

- genome FASTA: `mm10.fa`
- annotation GTF: `mm10.refGene.gtf` or a GENCODE annotation

For this repository’s analyses, the notebook used STAR-RSEM as the aligner and quantification strategy.

## 5. Run nf-core RNA-seq

The notebook workflow used an nf-core RNA-seq run with STAR and RSEM.

A typical command is:

```bash
nextflow run nf-core/rnaseq \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results \
  --genome GRCh38 \
  --aligner star_rsem \
  --profile docker \
  -resume
```

In the project draft, the equivalent workflow used a local copy of the nf-core RNA-seq pipeline and the following logic:

```bash
nextflow run ./git/nf-core-rnaseq_3.14.0/3_14_0/ \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results \
  --aligner star_rsem \
  -profile docker \
  -resume
```

This produces aligned BAM files, transcript and gene-level quantification, and the standard nf-core RNA-seq outputs.

## 6. Prepare a count matrix for differential expression

A common downstream step is to export gene-level counts for use in differential expression analysis. In the notebook workflow, the quantification output from the STAR-RSEM run was used to generate a raw counts table.

The resulting file was then passed to the differential abundance workflow.

## 7. Run nf-core differentialabundance

The repository draft used nf-core differentialabundance with a counts matrix, a sample design input file, a contrasts file, and the GTF annotation.

The command from the notebook draft was:

```bash
sudo nextflow run ./git/differentialabundance \
     --input SS_RNAseq_Lara_DEseq_input.tsv \
     --contrasts SS_RNAseq_Lara_DEseq_contrasts.tsv \
     --matrix  ./in/rawcounts.tsv \
     --gtf ./in/mm10.refGene.gtf.gz \
     --outdir nfout/build38/differentialabundance \
     -profile rnaseq,docker
```

This is the form to use when you want to run the analysis exactly as in the notebook workflow. The contrast file usually defines comparisons such as KO versus WT, or time-point-specific contrasts.

## 8. Inspect the main outputs

After the run completes, the results folder typically contains:

- `multiqc/` with the HTML QC report
- `star_rsem/` with BAM files and quantification outputs
- `tables/` or `results/` with differential expression tables
- annotation and plotting outputs produced by differentialabundance

Useful first checks:

- Open the MultiQC report to inspect read quality, trimming, mapping rate, duplication rate, and library complexity.
- Verify that the count matrix contains the expected samples and gene identifiers.
- Inspect the differential expression tables for the biological contrasts of interest.
- Check whether the strongest hits are consistent with the experimental design.

## 9. Typical downstream interpretation

Once the primary RNA-seq outputs are generated, the analysis usually continues with:

- gene-level differential expression analysis
- volcano plots and MA plots
- pathway or gene set enrichment analysis
- comparison with ChIP-seq peaks or regulatory elements
- integration with TE/L1 or Micro-C analyses when relevant

This repository’s downstream analysis examples often combine RNA-seq results with ChIP-seq and regulatory-element analyses for integrated interpretation.

## 10. Recommended workflow summary

A typical project follows this order:

1. Download and verify FASTQ files.
2. Build a samplesheet with sample, read, and strandedness information.
3. Run nf-core RNA-seq with STAR-RSEM.
4. Export a counts matrix for differential expression.
5. Run nf-core differentialabundance for the requested contrasts.
6. Review QC and DE results in the context of the biological question.

This structure is compatible with the nf-core RNA-seq and differentialabundance workflows used in this repository.
