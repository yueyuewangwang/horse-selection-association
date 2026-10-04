#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT in config/config.env}"
VCF2DIS="${VCF2DIS_BIN:-VCF2Dis}"
FASTME="${FASTME_BIN:-fastme}"
INPUT_VCF="${STRUCTURE_VCF:?Set STRUCTURE_VCF to a structure VCF}"
OUT_ROOT="${STRUCTURE_ROOT:-${PROJECT_ROOT}/work/structure}"
MATRIX="${OUT_ROOT}/p_distance.mat"
TREE="${OUT_ROOT}/population_tree.nwk"
[[ ! -e "${MATRIX}" ]] || { echo "Output exists: ${MATRIX}" >&2; exit 2; }
mkdir -p "${OUT_ROOT}"
"${VCF2DIS}" -InPut "${INPUT_VCF}" -OutPut "${MATRIX}"
"${FASTME}" -i "${MATRIX}" -o "${TREE}" -m BIONJ
echo "Distance matrix: ${MATRIX}"
echo "Tree: ${TREE}"
