# De novo transcriptome assembly guide

This guide describes how to run the RNA-seq de novo transcriptome assembly workflow used to evaluate whether KMT2B-bound L1 elements act as alternative promoters.

The repository includes two helper scripts in `bin/`:

- `bin/de_novo_transcriptome_assembly.sh`
- `bin/de_novo_transcriptome_summary.sh`

The assembly workflow uses `nf-core/rnaseq` v3.14.0 with STAR/Salmon alignment and StringTie transcript assembly.

## 1. Purpose

The analysis is intended to:

- align WT and double-KO RNA-seq reads to mm10 (UCSC GRCm38)
- assemble transcripts with StringTie via `--stringtie_fc`
- generate BigWig coverage tracks and assembled transcript files
- inspect whether transcript signal connects KMT2B-bound L1 loci with neighboring gene exons

## 2. Required inputs

By default, the assembly script expects:

- `./in/D0D4_WTDKO_nf-core_denovo_ss.csv`
  - nf-core sample sheet for WT and double-KO day 0/day 4 samples
- `./in/mm10.fa`
  - UCSC GRCm38/mm10 reference genome FASTA
- `./in/mm10.refGene.gtf`
  - UCSC GRCm38 gene annotation GTF

The script also expects the nf-core pipeline :

- `./git/nf-core-rnaseq_3.14.0/3_14_0/`

If your files are stored elsewhere, pass alternate paths to the script.

## 3. Running the assembly

From the repository root, run:

```bash
bin/de_novo_transcriptome_assembly.sh
```

To override defaults, pass the sample sheet, output directory, FASTA, and GTF:

```bash
bin/de_novo_transcriptome_assembly.sh \
  ./in/D0D4_WTDKO_nf-core_denovo_ss.csv \
  ./outs/denovo_transcripts \
  ./in/mm10.fa \
  ./in/mm10.refGene.gtf
```

The script will also create helper config files for StringTie and STAR multi-mapping behavior.

## 4. Output

The default output directory is:

- `./outs/denovo_transcripts`

Expected outputs include:

- `star_salmon/` subdirectories created by nf-core/rnaseq
- `stringtie/` transcript assembly outputs
- `bigwig/` coverage tracks
- `DE_NOVO_README.txt`

## 5. Inspecting results

Use the summary helper to inspect the generated files:

```bash
bin/de_novo_transcriptome_summary.sh ./outs/denovo_transcripts
```

This prints:

- the output directory tree
- assembled GTF/STRINGTIE files
- generated BigWig coverage files
- log, report, and trace files if present

## 6. Notes

- The workflow is configured for STAR/Salmon alignment via `--aligner star_salmon`.
- StringTie assembly is enabled with `--stringtie_fc`.
- The script sets `NXF_VER=23.10.0` to match the expected Nextflow runtime.
- If the sample sheet or reference files refer to a different genome build, update the paths accordingly.

## 7. Troubleshooting

- If the pipeline cannot find the sample sheet, confirm the path passed to the script.
- If `nextflow` is not installed, install it or use the containerized environment provided by your system.
- If the pipeline root `./git/nf-core-rnaseq_3.14.0/3_14_0/` does not exist, update `PIPELINE_ROOT` or create a symlink to the nf-core source.
