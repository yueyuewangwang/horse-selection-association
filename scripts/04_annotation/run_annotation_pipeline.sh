#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="${ANNOTATION_RESULT_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/results/annotation}"
if [[ "${RUN_VEP:-0}" == "1" ]]; then
  bash "${SCRIPT_DIR}/run_vep_annotation.sh"
fi
python3 "${SCRIPT_DIR}/annotate_candidates.py" --leads "${LEAD_TABLE:?Set LEAD_TABLE}" --regions "${SELECTION_REGIONS:?Set SELECTION_REGIONS}" --gff3 "${GFF3:?Set GFF3}" --out "${OUT_ROOT}/selection_gwas_candidate_annotation.tsv"
