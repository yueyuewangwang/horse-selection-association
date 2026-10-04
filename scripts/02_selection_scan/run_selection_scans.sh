#!/usr/bin/env bash
set -euo pipefail

# Run weighted FST and windowed pi for tab-separated group manifest.
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT in config/config.env}"
VCF="${SELECTION_VCF:?Set SELECTION_VCF}"
MANIFEST="${SELECTION_GROUP_MANIFEST:?Set SELECTION_GROUP_MANIFEST}"
GROUP_ROOT="${SELECTION_GROUP_ROOT:-$(dirname "${MANIFEST}")}"
OUT_ROOT="${SELECTION_RAW_ROOT:-${PROJECT_ROOT}/work/selection/raw}"
VCFTOOLS="${VCFTOOLS_BIN:-vcftools}"
WINDOW_SIZE="${WINDOW_SIZE:-50000}"
WINDOW_STEP="${WINDOW_STEP:-10000}"
[[ -s "${VCF}" && -s "${MANIFEST}" ]] || { echo "Missing selection VCF or manifest" >&2; exit 2; }
[[ ! -e "${OUT_ROOT}" ]] || { echo "Output exists: ${OUT_ROOT}" >&2; exit 2; }
mkdir -p "${OUT_ROOT}/fst" "${OUT_ROOT}/pi" "${OUT_ROOT}/logs"
while IFS=$'\t' read -r trait case_file control_file; do
  [[ -n "${trait}" && "${trait}" != "trait" ]] || continue
  [[ "${case_file}" = /* ]] || case_file="${GROUP_ROOT}/${case_file}"
  [[ "${control_file}" = /* ]] || control_file="${GROUP_ROOT}/${control_file}"
  [[ -s "${case_file}" && -s "${control_file}" ]] || { echo "Missing groups for ${trait}" >&2; exit 2; }
  prefix="${OUT_ROOT}/fst/${trait}"
  "${VCFTOOLS}" --gzvcf "${VCF}" --keep "${case_file}" --keep "${control_file}" \
    --weir-fst-pop "${case_file}" --weir-fst-pop "${control_file}" \
    --fst-window-size "${WINDOW_SIZE}" --fst-window-step "${WINDOW_STEP}" \
    --out "${prefix}" 2>"${OUT_ROOT}/logs/${trait}_fst.log"
  for arm in case control; do
    if [[ "${arm}" == case ]]; then keep_file="${case_file}"; else keep_file="${control_file}"; fi
    "${VCFTOOLS}" --gzvcf "${VCF}" --keep "${keep_file}" \
      --window-pi "${WINDOW_SIZE}" --window-pi-step "${WINDOW_STEP}" \
      --out "${OUT_ROOT}/pi/${trait}_${arm}" 2>"${OUT_ROOT}/logs/${trait}_${arm}_pi.log"
  done
done < "${MANIFEST}"
echo "Selection scans completed: ${OUT_ROOT}"
