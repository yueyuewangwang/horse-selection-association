#!/usr/bin/env bash
set -euo pipefail
PLINK2="${PLINK2_BIN:-plink2}"
INPUT_PREFIX="${STRUCTURE_BFILE:?Set STRUCTURE_BFILE to a PLINK prefix}"
OUT_ROOT="${STRUCTURE_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/structure}"
THREADS="${THREADS:-8}"
mkdir -p "${OUT_ROOT}"
OUT_PREFIX="${OUT_ROOT}/structure"
[[ ! -e "${OUT_PREFIX}.eigenvec" ]] || { echo "Output exists: ${OUT_PREFIX}.eigenvec" >&2; exit 2; }
"${PLINK2}" --bfile "${INPUT_PREFIX}" --chr 1-31 --pca 10 approx --threads "${THREADS}" --out "${OUT_PREFIX}"
echo "PCA completed: ${OUT_PREFIX}.eigenvec"
