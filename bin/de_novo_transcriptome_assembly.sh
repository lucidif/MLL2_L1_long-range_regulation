#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Run RNA-seq de novo transcriptome assembly with nf-core/rnaseq
#
# This script launches nf-core/rnaseq (v3.14.0) with STAR/Salmon alignment
# and StringTie transcript assembly, matching the analysis used for L1
# alternative promoter evaluation.
#
# The default configuration is intentionally conservative and can be
# overridden by environment variables or positional arguments below.
# =============================================================================

ROOT_DIR="$(pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-./git/nf-core-rnaseq_3.14.0/3_14_0/}"
SAMPLE_SHEET="${1:-./in/D0D4_WTDKO_nf-core_denovo_ss.csv}"
OUTPUT_DIR="${2:-./outs/denovo_transcripts}"
FASTA="${3:-./in/mm10.fa}"
GTF="${4:-./in/mm10.refGene.gtf}"
PROFILE="${PROFILE:-docker}"
NXF_VER="${NXF_VER:-23.10.0}"

STRINGTIE_CONFIG="STRINGTIE_custom.config"
MULTIHIT_CONFIG="multihit_custom.config"

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: Sample sheet not found: ${SAMPLE_SHEET}" >&2
  exit 1
fi

if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: Reference FASTA not found: ${FASTA}" >&2
  exit 1
fi

if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: Reference GTF not found: ${GTF}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

cat > "${STRINGTIE_CONFIG}" <<'EOF'
process {
    withName: 'STRINGTIE_STRINGTIE' {
        ext.args = { ['-v'].join(' ').trim() }
    }
}
EOF

cat > "${MULTIHIT_CONFIG}" <<'EOF'
process {
    withName: '.*STAR_ALIGN' {
        ext.args = [
            '--outSAMtype BAM Unsorted',
            '--outFilterMultimapNmax 5000',
            '--outSAMmultNmax 1',
            '--outFilterMismatchNmax 3',
            '--outMultimapperOrder Random',
            '--winAnchorMultimapNmax 5000',
            '--alignEndsType EndToEnd',
            '--alignIntronMax 1000000',
            '--alignMatesGapMax 350',
            '--seedSearchStartLmax 30',
            '--alignTranscriptsPerReadNmax 30000',
            '--alignWindowsPerReadNmax 30000',
            '--alignTranscriptsPerWindowNmax 300',
            '--seedPerReadNmax 3000',
            '--seedPerWindowNmax 300',
            '--seedNoneLociPerWindow 1000',
            '--readFilesCommand zcat'
        ].join(' ').trim()
    }
    withName: 'STRINGTIE_STRINGTIE' {
        ext.args = { ['-v'].join(' ').trim() }
    }
}
EOF

export NXF_VER="${NXF_VER}"

echo "[INFO] Running nf-core/rnaseq de novo transcriptome assembly"
echo "       sample sheet: ${SAMPLE_SHEET}"
echo "       output dir:   ${OUTPUT_DIR}"
echo "       fasta:        ${FASTA}"
echo "       gtf:          ${GTF}"
echo "       pipeline:     ${PIPELINE_ROOT}"

time nextflow run "${PIPELINE_ROOT}" \
  --input "${SAMPLE_SHEET}" \
  --outdir "${OUTPUT_DIR}" \
  --fasta "${FASTA}" \
  --gtf "${GTF}" \
  --aligner star_salmon \
  --skip_umi_extract \
  --save_reference \
  --stringtie_fc \
  -profile "${PROFILE}" \
  -c "${STRINGTIE_CONFIG}" \
  -resume

cat > "${OUTPUT_DIR}/DE_NOVO_README.txt" <<EOF
De novo transcriptome assembly run completed.

Sample sheet: ${SAMPLE_SHEET}
Output dir: ${OUTPUT_DIR}
Reference FASTA: ${FASTA}
Reference GTF: ${GTF}
Pipeline root: ${PIPELINE_ROOT}
EOF

echo "[DONE] de novo transcriptome assembly completed"
