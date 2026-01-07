#!/usr/bin/env bash
set -euo pipefail

# Concatenate lane-split FASTQs per sample
# Creates:
#   <sample>_combined_forward.fq.gz
#   <sample>_combined_reverse.fq.gz
#
# Run from inside the raw data directory

for d in */ ; do
  sample="${d%/}"
  echo "Processing ${sample}"

  R1=$(ls "${d}"/*_1.f*q.gz 2>/dev/null || true)
  R2=$(ls "${d}"/*_2.f*q.gz 2>/dev/null || true)

  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "  Skipping ${sample} (no FASTQs found)"
    continue
  fi

  cat ${R1} > "${d}/${sample}_combined_forward.fq.gz"
  cat ${R2} > "${d}/${sample}_combined_reverse.fq.gz"
done

echo "Lane concatenation complete."
