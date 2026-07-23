# Typical Micro-C analysis guide

This guide summarizes a standard Micro-C analysis workflow for this repository, based on the local nf-core Micro-C pipeline in [pipelines/nf-core-microc](../pipelines/nf-core-microc) and the repository’s example sample sheets.

The workflow usually covers: raw data download, sample organization, samplesheet preparation, alignment and parsing, contact matrix generation, and downstream analysis such as A/B compartments, TAD calling, and visualization.

## 1. Prepare the project structure

A simple project layout is helpful before running the workflow:

```bash
mkdir -p analysis/microc/fastq analysis/microc/results analysis/microc/in
cd analysis/microc
```

A practical folder organization is:

- `fastq/` for raw FASTQ files
- `in/` for sample sheets, genome index files, and chromosome sizes
- `results/` or `nfout/` for pipeline outputs

## 2. Download and verify raw reads

The notebook workflow downloaded Micro-C FASTQ files from the sequencing provider and verified them with MD5 checksum files before analysis.

A typical approach is:

```bash
cd fastq
md5sum *.fastq.gz > md5_checksums.txt
```

For the project draft, the workflow also used the helper script in [pipelines/nf-core-microc/bin/download_from_macrogen.sh](../pipelines/nf-core-microc/bin/download_from_macrogen.sh) to fetch data from the provider using a sample-sheet table.

## 3. Create a samplesheet for nf-core Micro-C

The nf-core Micro-C pipeline expects a samplesheet that lists each sample, its paired FASTQ files, and sample metadata such as condition and replicate.

A minimal example is:

```csv
sample,fastq_1,fastq_2,condition,replicate
WT_day0_A,fastq/WT_day0_A_1.fastq.gz,fastq/WT_day0_A_2.fastq.gz,WT_D0,A
WT_day0_B,fastq/WT_day0_B_1.fastq.gz,fastq/WT_day0_B_2.fastq.gz,WT_D0,B
KO_day0_A,fastq/KO_day0_A_1.fastq.gz,fastq/KO_day0_A_2.fastq.gz,KO_D0,A
KO_day0_B,fastq/KO_day0_B_1.fastq.gz,fastq/KO_day0_B_2.fastq.gz,KO_D0,B
```

A concrete example matching the notebook workflow is available in [assets/great/microC_nf-core_test_sheet.csv](../assets/great/microC_nf-core_test_sheet.csv).

## 4. Prepare reference files

A standard Micro-C run requires:

- the reference genome FASTA
- a BWA index for the genome
- chromosome sizes for Juicer-compatible downstream steps

Typical inputs include:

- `mm10.fa`
- `mm10.fa.bwt`, `mm10.fa.ann`, `mm10.fa.sa`, and related BWA index files
- `mm10.sizes` or a similar chromosome-size table

The notebook workflow used a local copy of the pipeline and prepared the run with the required genome inputs and Juicer tools.

## 5. Run nf-core Micro-C

The project notebook used a local clone of the nf-core Micro-C pipeline and executed it with Docker.

A typical command is:

```bash
nextflow run nf-core/microc \
  -profile docker \
  --input /path/to/samplesheet.csv \
  --outdir /path/to/results
```

In the notebook workflow, the equivalent command was structured more explicitly around the local pipeline checkout and Juicer tools:

```bash
cd /path/to/project
nextflow run git/nf-core-microc \
  -profile docker \
  --input SS_microC_NFSS_2024_03_microC_Lara.csv \
  --outdir nfout \
  --juicertool_location /path/to/git/nf-core-microc/bin/juicertools.jar
```

For a smaller test run, the notebook also used a one-sample input sheet:

```bash
nextflow run git/nf-core-microc \
  -profile docker \
  --input SS_onesample_microC_NFSS_2024_03_microC_Lara.csv \
  --outdir nfout \
  --juicertool_location /path/to/git/nf-core-microc/bin/juicertools.jar
```

The pipeline performs alignment, parsing, deduplication, pair filtering, and contact matrix generation.

## 6. Inspect the main outputs

After the run completes, the results directory usually contains:

- `pairtools/` outputs for parsed and deduplicated contacts
- `cooler/` or related contact-matrix outputs
- QC summaries and logs
- downstream contact files suitable for Hi-C/Micro-C visualization

Useful first checks:

- verify that the number of valid pairs looks reasonable for each sample
- inspect the QC and parsing logs for failed or low-complexity libraries
- confirm that the expected output files were generated for each sample

## 7. Generate contact matrices and Hi-C-like products

The notebook workflow continued beyond the initial alignment step and generated downstream contact products for visualization and interpretation.

Typical follow-up steps include:

- generating `.hic` files with Juicer tools
- creating balanced contact matrices at multiple resolutions
- preparing files for visualization in HiGlass or similar tools

An example from the notebook used Juicer pre commands such as:

```bash
java -Xms512m -Xmx32g -jar ./git/nf-core-microc/bin/juicertools.jar pre -r 5000 --threads 12 \
  ./in/WT_day0_A.Dd.pairs.gz ./out/juicertools_pre/WT_day0_A.Dd.hic ./in/mm10.sizes
```

These steps are useful when you want to inspect contact maps at different resolutions or share the results in a browser-based viewer.

## 8. Run downstream analyses

Standard downstream analyses for Micro-C commonly include:

- A/B compartment analysis
- TAD calling
- saddle plots
- insulation score analysis
- virtual 4C or targeted interaction plots

The notebook workflow included downstream analyses for compartments, TADs, and virtual 4C-style interpretations, which are often the most informative ways to connect the contact map to biology.

## 9. Typical downstream interpretation

Once the primary Micro-C outputs are generated, the analysis usually proceeds with:

- checking whether contacts are enriched at promoters or enhancers of interest
- comparing contact strength between conditions or genotypes
- correlating changes in chromatin architecture with RNA-seq or ChIP-seq signals
- interpreting changes in compartmentalization or TAD boundary strength in the biological context

This repository’s downstream examples often combine Micro-C data with ChIP-seq and RNA-seq outputs for integrated interpretation.

## 10. Recommended workflow summary

A typical project follows this order:

1. Download and verify the FASTQ files.
2. Prepare a samplesheet with sample, FASTQ, condition, and replicate information.
3. Prepare the reference genome and chromosome-size files.
4. Run the nf-core Micro-C pipeline with the appropriate input sheet and Juicer tool path.
5. Inspect the parsed and deduplicated contact outputs.
6. Generate `.hic` or contact matrices for visualization and downstream analysis.

This structure is compatible with the nf-core Micro-C workflow used in this repository.
