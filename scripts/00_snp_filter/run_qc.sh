#!/usr/bin/env bash
set -euo pipefail

# Sample and variant QC, post-PCA exclusion, LD pruning, and PCA input.
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT in config/config.env}"
PLINK2="${PLINK2_BIN:-plink2}"
RAW_VCF="${RAW_VCF:?Set RAW_VCF in config/config.env}"
REMOVE_MULE_LIST="${REMOVE_MULE_LIST:?Set REMOVE_MULE_LIST in config/config.env}"
REMOVE_SAMPLE_LIST="${REMOVE_SAMPLE_LIST:?Set REMOVE_SAMPLE_LIST in config/config.env}"
OUT_ROOT="${QC_ROOT:-${PROJECT_ROOT}/work/qc}"
THREADS="${THREADS:-8}"

for path in "${RAW_VCF}" "${REMOVE_MULE_LIST}" "${REMOVE_SAMPLE_LIST}"; do
  [[ -s "${path}" ]] || { echo "Missing input: ${path}" >&2; exit 2; }
done
[[ ! -e "${OUT_ROOT}" ]] || { echo "Output exists: ${OUT_ROOT}" >&2; exit 2; }
mkdir -p "${OUT_ROOT}"
run_plink() { local name="$1"; shift; "${PLINK2}" "$@" --threads "${THREADS}" --chr-set 31 --out "${OUT_ROOT}/${name}"; }

run_plink sample_qc --vcf "${RAW_VCF}" --remove "${REMOVE_MULE_LIST}" --mind 0.1 --max-alleles 2 --make-bed --impute-sex
run_plink site_qc --bfile "${OUT_ROOT}/sample_qc" --geno 0.1 --make-bed --impute-sex
run_plink final_qc --bfile "${OUT_ROOT}/site_qc" --remove "${REMOVE_SAMPLE_LIST}" --geno 0.1 --mac 1 --make-bed
run_plink ld_pruned --bfile "${OUT_ROOT}/final_qc" --chr 1-31 --maf 0.05 --indep-pairwise 50 10 0.2
run_plink structure_markers --bfile "${OUT_ROOT}/final_qc" --chr 1-31 --extract "${OUT_ROOT}/ld_pruned.prune.in" --make-bed
run_plink pca --bfile "${OUT_ROOT}/structure_markers" --pca 10
echo "QC/PCA completed under ${OUT_ROOT}"
