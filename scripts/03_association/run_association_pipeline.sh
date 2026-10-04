#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${GWAS_TRAIT_MANIFEST:?Set GWAS_TRAIT_MANIFEST}"
PHENO="${PHENOTYPE_ROOT:?Set PHENOTYPE_ROOT}"
GENO="${GWAS_AUTOSOMAL_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/gwas/autosomal}"
MATRIX="${GWAS_MATRIX_ROOT:-${PROJECT_ROOT}/work/gwas/matrices}"
RESULT="${GWAS_RESULT_ROOT:-${PROJECT_ROOT}/results/association}"
bash "${SCRIPT_DIR}/prepare_autosomal_genotypes.sh"
"${RSCRIPT_BIN:-Rscript}" "${SCRIPT_DIR}/convert_to_rMVP.R" --manifest "${MANIFEST}" --geno-root "${GENO}" --matrix-root "${MATRIX}" --ncpus "${THREADS:-1}"
"${RSCRIPT_BIN:-Rscript}" "${SCRIPT_DIR}/run_gwas.R" --manifest "${MANIFEST}" --matrix-root "${MATRIX}" --phenotype-root "${PHENO}" --output-root "${RESULT}" --ncpus "${THREADS:-1}"
python3 "${SCRIPT_DIR}/summarize_gwas.py" --manifest "${MANIFEST}" --matrix-root "${MATRIX}" --result-root "${RESULT}" --out-root "${RESULT}"
