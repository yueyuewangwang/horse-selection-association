#!/usr/bin/env bash
set -euo pipefail
PLINK2="${PLINK2_BIN:-plink2}"
INPUT_PREFIX="${STRUCTURE_BFILE:?Set STRUCTURE_BFILE}"
OUT_ROOT="${ADMIXTURE_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/admixture}"
STRIDE="${ADMIXTURE_MARKER_STRIDE:-16}"
THREADS="${THREADS:-8}"
mkdir -p "${OUT_ROOT}"
[[ ! -e "${OUT_ROOT}/admixture_markers.bed" ]] || { echo "Output exists: ${OUT_ROOT}/admixture_markers.bed" >&2; exit 2; }
awk -v stride="${STRIDE}" 'NR % stride == 1 {print $2}' "${INPUT_PREFIX}.bim" > "${OUT_ROOT}/admixture_markers.txt"
"${PLINK2}" --bfile "${INPUT_PREFIX}" --chr 1-31 --extract "${OUT_ROOT}/admixture_markers.txt" --make-bed --threads "${THREADS}" --out "${OUT_ROOT}/admixture_markers"
echo "ADMIXTURE marker set: ${OUT_ROOT}/admixture_markers.bed"
