#!/usr/bin/env bash
set -euo pipefail

# Scenario manifest columns: scenario, trait, case_file, control_file.
VCF="${SELECTION_VCF:?Set SELECTION_VCF}"
MANIFEST="${ROBUSTNESS_SCENARIO_MANIFEST:?Set ROBUSTNESS_SCENARIO_MANIFEST}"
GROUP_ROOT="${SELECTION_GROUP_ROOT:-$(dirname "${MANIFEST}")}"
OUT_ROOT="${ROBUSTNESS_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/results/robustness}/selection"
VCFTOOLS="${VCFTOOLS_BIN:-vcftools}"; SIZE="${WINDOW_SIZE:-50000}"; STEP="${WINDOW_STEP:-10000}"
mkdir -p "${OUT_ROOT}"
while IFS=$'\t' read -r scenario trait case_file control_file; do
  [[ -n "${scenario}" && "${scenario}" != scenario ]] || continue
  [[ "${case_file}" = /* ]] || case_file="${GROUP_ROOT}/${case_file}"
  [[ "${control_file}" = /* ]] || control_file="${GROUP_ROOT}/${control_file}"
  out="${OUT_ROOT}/${scenario}/${trait}"; mkdir -p "${out}"
  "${VCFTOOLS}" --gzvcf "${VCF}" --keep "${case_file}" --keep "${control_file}" --weir-fst-pop "${case_file}" --weir-fst-pop "${control_file}" --fst-window-size "${SIZE}" --fst-window-step "${STEP}" --out "${out}/fst"
  for arm in case control; do
    keep="${case_file}"; [[ "${arm}" == control ]] && keep="${control_file}"
    "${VCFTOOLS}" --gzvcf "${VCF}" --keep "${keep}" --window-pi "${SIZE}" --window-pi-step "${STEP}" --out "${out}/pi_${arm}"
  done
done < "${MANIFEST}"
