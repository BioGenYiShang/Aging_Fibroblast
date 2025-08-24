#!/usr/bin/bin/env bash
set -euo pipefail
# Usage: bash scripts/bash/spaceranger_count.sh SAMPLE FASTQ_DIR IMAGE_PATH REF OUTDIR
SAMPLE=$1; FASTQ=$2; IMAGE=$3; REF=$4; OUT=$5
mkdir -p "${OUT}"
spaceranger count --id="${SAMPLE}" --transcriptome="${REF}" --fastqs="${FASTQ}" --image="${IMAGE}" --localcores=8 --localmem=64
mv "${SAMPLE}" "${OUT}/"
