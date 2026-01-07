#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="${INPUT_DIR:-01.RawData}"
TRIM_DIR="${TRIM_DIR:-Trim_data}"
LOG_DIR="${LOG_DIR:-Trim_Logs}"
THREADS="${THREADS:-6}"
ADAPTER="${ADAPTER:-adapter.fasta}"

mkdir -p "${TRIM_DIR}" "${LOG_DIR}"

find "${INPUT_DIR}" -name '*_combined_forward.fq.gz' | while read R1; do
  R2="${R1/_combined_forward/_combined_reverse}"
  sample=$(basename "${R1}" | sed 's/_combined_forward\.fq\.gz//')

  fastp \
    -i "${R1}" -I "${R2}" \
    -o "${TRIM_DIR}/${sample}_trimmed_forward.fq.gz" \
    -O "${TRIM_DIR}/${sample}_trimmed_reverse.fq.gz" \
    -w "${THREADS}" \
    --detect_adapter_for_pe \
    --adapter_fasta "${ADAPTER}" \
    -h "${LOG_DIR}/${sample}.html" \
    -j "${LOG_DIR}/${sample}.json"
done

echo "fastp trimming completed."
