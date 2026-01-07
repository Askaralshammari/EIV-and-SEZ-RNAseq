#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-6}"
OUTDIR="${OUTDIR:-qc_reports}"

mkdir -p "${OUTDIR}"

find . -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) -print0 \
  | xargs -0 fastqc --threads "${THREADS}" --outdir "${OUTDIR}"

multiqc "${OUTDIR}" -o "${OUTDIR}"

echo "FastQC and MultiQC completed."
