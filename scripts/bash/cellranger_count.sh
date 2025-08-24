#!/usr/bin/env bash
set -euo pipefail
# Usage: bash scripts/bash/cellranger_count.sh SAMPLE FASTQ_DIR REFdata-GRCh38-2020-A OUTDIR
SAMPLE=$1; FASTQ=$2; REF=$3; OUT=$4
mkdir -p "${OUT}"
cellranger count --id="${SAMPLE}" --fastqs="${FASTQ}" --transcriptome="${REF}" --sample="${SAMPLE}" --expect-cells=6000 --localcores=8 --localmem=64 --nosecondary
mv "${SAMPLE}" "${OUT}/"
