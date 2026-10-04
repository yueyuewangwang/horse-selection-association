#!/usr/bin/env bash
set -euo pipefail

# Convert final QC genotypes to a sorted, phased, low-frequency-inclusive VCF.
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT in config/config.env}"
PLINK2="${PLINK2_BIN:-plink2}"; BCFTOOLS="${BCFTOOLS_BIN:-bcftools}"
BEAGLE_JAR="${BEAGLE_JAR:?Set BEAGLE_JAR in config/config.env}"
QC_PREFIX="${QC_PREFIX:-${QC_ROOT:?Set QC_ROOT}/final_qc}"
OUT_ROOT="${IMPUTATION_ROOT:-${PROJECT_ROOT}/work/imputation}"
THREADS="${THREADS:-8}"; MEMORY="${BEAGLE_MEMORY:-16g}"
for suffix in bed bim fam; do [[ -s "${QC_PREFIX}.${suffix}" ]] || { echo "Missing ${QC_PREFIX}.${suffix}" >&2; exit 2; }; done
[[ ! -e "${OUT_ROOT}" ]] || { echo "Output exists: ${OUT_ROOT}" >&2; exit 2; }
mkdir -p "${OUT_ROOT}"
"${PLINK2}" --bfile "${QC_PREFIX}" --export vcf id-paste=iid --allow-extra-chr --chr-set 31 --chr 1-31,X --threads "${THREADS}" --out "${OUT_ROOT}/qc"
grep -v '^##chrSet=' "${OUT_ROOT}/qc.vcf" > "${OUT_ROOT}/qc_clean.vcf"
"${BCFTOOLS}" sort "${OUT_ROOT}/qc_clean.vcf" -O z -o "${OUT_ROOT}/qc_sorted.vcf.gz"
"${BCFTOOLS}" index -t "${OUT_ROOT}/qc_sorted.vcf.gz"
java -Xmx"${MEMORY}" -jar "${BEAGLE_JAR}" gt="${OUT_ROOT}/qc_sorted.vcf.gz" out="${OUT_ROOT}/phased" nthreads="${THREADS}"
"${PLINK2}" --vcf "${OUT_ROOT}/phased.vcf.gz" --maf 0.01 --export vcf id-paste=iid bgz --allow-extra-chr --chr-set 31 --threads "${THREADS}" --out "${OUT_ROOT}/selection"
echo "Selection VCF: ${OUT_ROOT}/selection.vcf.gz"
