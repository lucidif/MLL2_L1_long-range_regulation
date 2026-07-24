# Virtual 4C analysis guide

This guide describes how to run the virtual 4C analysis workflow used to investigate physical contacts between KMT2B-bound L1 elements and their putative target genes, using Micro-C data from WT and double-KO (DKO) cells at day 0 and day 4 of differentiation.

The workflow is organized around one master script per locus of interest, following the naming convention:

- `virtual4C_<REF_NAME>.sh` (e.g. `virtual4C_chr10_Smpdl3a.sh`)

Each master script calls three helper scripts from the `nf-core-microc` pipeline repository (`git/nf-core-microc/bin/`):

- `virtual4C_make_reference.sh`
- `counts_pairs.sh`
- `virtual4C_make_bw_tracks_and_pool.sh` (which in turn calls `virtual4c_from_pairs_and_bam.sh`)

and runs Micro-C alignment via the `nf-core-microc` Nextflow pipeline (Dovetail Genomics Micro-C protocol).

## 1. Purpose

The analysis is intended to:

- extract a locus-specific reference sequence spanning the region between the TSS of a target gene and the nearby KMT2B-bound L1 element, with 50 Kb padding on each side
- align Micro-C reads from WT and double-KO cells (day 0 and day 4) to this custom reference
- select valid pairs with at least one end mapping within a viewpoint window centered on the target gene TSS (±1 to ±1.5 Kb, defined per locus)
- lift the filtered, viewpoint-selected pairs back to mm10 coordinates
- merge replicates per condition and generate normalized, pooled BigWig tracks for WT vs double-KO comparison

## 2. Required inputs

By default, the master script expects:

- `sheets/VirtualC_ref_coordinates.tsv`
  - locus coordinates (TSS, L1 element, padding) for each target gene, used by `virtual4C_make_reference.sh` to build the custom reference. One row per target gene, with columns:
    - `Target_Gene` - gene symbol (e.g. `Smpdl3a`)
    - `TSS_viewpoint` - TSS coordinate of the target gene (`chr:pos`), used as the virtual 4C viewpoint
    - `L1_viewpoint` - coordinate of the KMT2B-bound L1 element (`chr:pos`)
    - `chr` - chromosome
    - `start`, `end` - unpadded locus boundaries, i.e. the min/max of `TSS_viewpoint` and `L1_viewpoint`
    - `distance` - distance in bp between the TSS and the L1 element (`end - start`)
    - `window` - flanking padding applied on each side when building the custom reference (50,000 bp)
    - `recalibrated_start`, `recalibrated_end` - final custom reference boundaries (`start - window`, `end + window`), i.e. the coordinates actually extracted from mm10
    - `fastaName` - name of the locus-specific custom reference/contig (`chr_GeneName`), used as `--ref-name` and as the contig name restored during lift-back to mm10 coordinates
- `/mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa` and matching `.fai`
  - UCSC GRCm38/mm10 reference genome FASTA and index
- `sheets/NFSS_Lara_<condition>_microC.csv`
  - one nf-core-microc sample sheet per condition (WT day 0, WT day 4, double-KO day 0, double-KO day 4)
- `outs/Lara_multiomic_analysis/in/mm10.sizes`
  - chromosome sizes file, used when generating BigWig tracks
- `git/Lara_MLL2/bin/docker_env.sh`
  - Docker environment file, sourced before running any step

The pipeline also expects the nf-core-microc pipeline source at:

- `git/nf-core-microc`

## 3. Running the analysis

From the repository root, source the Docker environment and build the locus-specific reference:

```bash
source git/Lara_MLL2/bin/docker_env.sh

bash git/nf-core-microc/bin/virtual4C_make_reference.sh \
  --ref-name chr10_Smpdl3a \
  --tsv sheets/VirtualC_ref_coordinates.tsv \
  --mm10-fa /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fa \
  --mm10-fai-src /mnt/datawk1/references/fasta/UCSC_GRCm38/mm10.fai \
  --custom-ref-root in/customReferences \
  --outdir-main outs/virtual4C \
  --docker-env git/Lara_MLL2/bin/docker_env.sh
```

Then, for each of the four conditions (WT day 0, WT day 4, double-KO day 0, double-KO day 4), run the Micro-C alignment:

```bash
NF_OUT="outs/virtual4C/nfout/chr10_Smpdl3a_WT_day4"
NF_WORK="outs/virtual4C/work_chr10_Smpdl3a_WT_day4"

sudo nextflow run git/nf-core-microc -profile docker \
  -work-dir "${NF_WORK}" \
  --input sheets/NFSS_Lara_WT_day4_microC.csv \
  --fasta in/customReferences/chr10_Smpdl3a/genome.fa \
  --index in/customReferences/chr10_Smpdl3a \
  --outdir "${NF_OUT}" \
  --exclude_lcextrap false
```

followed by pair counting and BigWig track generation:

```bash
git/nf-core-microc/bin/counts_pairs.sh "${NF_OUT}/pairix" git/Lara_MLL2/bin/docker_env.sh

git/nf-core-microc/bin/virtual4C_make_bw_tracks_and_pool.sh \
  --ref-name chr10_Smpdl3a \
  --viewpoint 57794374 \
  --window 1000 \
  --mm10-sizes outs/Lara_multiomic_analysis/in/mm10.sizes \
  --nf-out "${NF_OUT}" \
  --outroot outs/virtual4C \
  --sam-prefix WT_day4 \
  --docker-env git/Lara_MLL2/bin/docker_env.sh \
  --v4c-script git/nf-core-microc/bin/virtual4c_from_pairs_and_bam.sh \
  --pairs-counts-tsv "${NF_OUT}/pairix/total_nodups.tsv" \
  --smoothLength 0
```

Repeat the last three steps (`nextflow run`, `counts_pairs.sh`, `virtual4C_make_bw_tracks_and_pool.sh`) for each of the remaining conditions, changing the sample sheet, the `NF_OUT`/`NF_WORK` paths, and `--sam-prefix` accordingly (`KO_day4`, `WT_day0`, `KO_day0`).

To analyze a different locus, pass a different `--ref-name`, `--viewpoint`, and the corresponding entry in `--tsv`, and point `--fasta`/`--index` at the matching custom reference directory under `in/customReferences/`.

## 4. Output

The default output root is:

- `outs/virtual4C`

Expected outputs include:

- `nfout/<REF_NAME>_<condition>/` - nf-core-microc outputs per condition, including `pairix/` (valid pairs, `total_nodups.tsv`)
- pooled BigWig tracks per condition, scaled by `1e6 / sum(total_nodups)` (library-size normalization, computed on the merged replicates)
- `vpnorm/` - viewpoint-locally normalized pooled BigWig tracks, scaled using pairs counted within the viewpoint window; these are the tracks used for WT vs double-KO comparisons
- `in/customReferences/<REF_NAME>/` - locus-specific reference FASTA and BWA/samtools index files

## 5. Inspecting results

Check the nf-core-microc report and trace files under each `nfout/<REF_NAME>_<condition>/` directory to confirm the run completed successfully, and inspect `pairix/total_nodups.tsv` to confirm the expected number of non-duplicate valid pairs per replicate before pooling.

Load the pooled BigWig tracks (main and `vpnorm/`) alongside the target gene annotation and the KMT2B-bound L1 element coordinates in a genome browser to compare contact profiles between WT and double-KO at day 0 and day 4.

## 6. Notes

- The nf-core-microc pipeline implements the Dovetail Genomics Micro-C protocol: BWA-MEM alignment, `pairtools parse` with `--min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30`, followed by `pairtools dedup` (`--exclude_lcextrap false`).
- The viewpoint window (±1,000-1,500 bp around the TSS) is defined per locus based on the local density of repetitive elements (Supplementary Data 5).
- Only primary, non-supplementary alignments are retained after viewpoint-based read filtering (`samtools view -N -F 0x900`).
- Replicates are merged after per-replicate lifting to mm10 coordinates; the normalization scale factor (`1e6 / sum(total_nodups)`) is computed on the merged, pooled set of non-duplicate valid pairs and applied to the merged BAM before generating the BigWig track with `bamCoverage`.
- The script sets `NXF_VER=23.10.0` to match the expected Nextflow runtime.
- `sudo` is used for the `nextflow run` step in the example script; confirm this matches your environment's Docker permissions.

## 7. Troubleshooting

- If `virtual4C_make_reference.sh` cannot find the locus, confirm `--ref-name` matches an entry in `sheets/VirtualC_ref_coordinates.tsv`.
- If the nf-core-microc pipeline cannot find the sample sheet, confirm the path passed to `--input`.
- If `--fasta`/`--index` cannot be found, confirm the custom reference was built successfully in `in/customReferences/<REF_NAME>/`.
- If `counts_pairs.sh` reports zero pairs, confirm the nf-core-microc run completed and `pairix/` was populated.
- If the Docker environment is not sourced, `"${DOCKER_ARGS[@]}"` mounts may not resolve; re-source `git/Lara_MLL2/bin/docker_env.sh` before rerunning any step.
