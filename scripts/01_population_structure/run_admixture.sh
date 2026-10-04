#!/usr/bin/env bash
set -euo pipefail
ADMIXTURE="${ADMIXTURE_BIN:-admixture}"
BED_PREFIX="${ADMIXTURE_BED:?Set ADMIXTURE_BED}"
OUT_ROOT="${ADMIXTURE_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/admixture}"
K_MIN="${K_MIN:-2}"; K_MAX="${K_MAX:-8}"; THREADS="${ADMIXTURE_THREADS:-8}"
CV_FOLDS="${CV_FOLDS:-5}"; SEED="${ADMIXTURE_SEED:-1}"; CONVERGENCE="${ADMIXTURE_CONVERGENCE:-0.0001}"
mkdir -p "${OUT_ROOT}/results" "${OUT_ROOT}/logs"
for k in $(seq "${K_MIN}" "${K_MAX}"); do
  args=("-j${THREADS}" "--seed=${SEED}" "-C=${CONVERGENCE}")
  if [[ "${CV_FOLDS}" -gt 0 ]]; then args+=("--cv=${CV_FOLDS}"); fi
  (cd "${OUT_ROOT}/results" && "${ADMIXTURE}" "${args[@]}" "${BED_PREFIX}.bed" "${k}") >"${OUT_ROOT}/logs/K${k}.log" 2>&1
done
if [[ "${CV_FOLDS}" -gt 0 ]]; then
  tmp="${OUT_ROOT}/cv_errors.body.tsv"
  : > "${tmp}"
  for log in "${OUT_ROOT}"/logs/K*.log; do
    k="$(basename "${log}" .log | sed 's/^K//')"
    awk -v k="${k}" '/CV error/ {print k "\t" $4}' "${log}" | tail -1 >> "${tmp}"
  done
  { printf 'K\tCV_error\n'; sort -k1,1n "${tmp}"; } > "${OUT_ROOT}/cv_errors.tsv"
  rm -f "${tmp}"
fi
echo "ADMIXTURE complete: ${OUT_ROOT}/results"
