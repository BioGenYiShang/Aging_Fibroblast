#!/usr/bin/env bash
set -euo pipefail
REMOTE=${1:-}
if [ -z "${REMOTE}" ]; then
  echo "Provide remote URL (e.g., https://github.com/<you>/skin-aging-atlas.git)"; exit 1;
fi
git init
git add .
git commit -m "Initial commit: analysis pipeline and docs"
git branch -M main || true
git remote add origin "${REMOTE}" || git remote set-url origin "${REMOTE}"
git push -u origin main
echo "[INFO] Enable GitHub Pages: Settings → Pages → Deploy from branch → /docs"
