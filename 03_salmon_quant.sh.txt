#!/usr/bin/env bash
set -euo pipefail

TRIM_DIR="${TRIM_DIR:-Trim_data}"
QUANT_DIR="${QUANT_DIR:-Quants}"
THREADS="${THREADS:-8}"
INDEX_DIR="${INDEX_DIR:?Please set INDEX_DIR}"

mkdir -p "${QUANT_DIR}"

for R1 in "${TRIM_DIR}"/*_trimmed_forward.fq.gz; do
  R2="${R1/_trimmed_forward/_trimmed_reverse}"
  sample=$(basename "${R1}" | sed 's/_trimmed_forward\.fq\.gz//')

  salmon quant \
    -i "${INDEX_DIR}" \
    -l A \
    -1 "${R1}" -2 "${R2}" \
    -p "${THREADS}" \
    --validateMappings \
    --gcBias \
    -o "${QUANT_DIR}/${sample}"
done

echo "Salmon quantification completed."
