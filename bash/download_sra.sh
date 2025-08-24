#!/usr/bin/env bash
set -euo pipefail
# Usage: bash scripts/bash/download_sra.sh GSE130973 data/raw/GSE130973
ACC=${1:-GSE130973}
OUTDIR=${2:-data/raw/${ACC}}
mkdir -p "${OUTDIR}"
echo "[INFO] Downloading ${ACC} into ${OUTDIR}"
if command -v prefetch >/dev/null 2>&1; then
  prefetch "${ACC}" || true
  fasterq-dump --split-files --outdir "${OUTDIR}" "${ACC}" || true
else
  echo "[WARN] sratoolkit not found. Install sratoolkit or use Aspera Connect (ascp)."
fi
