#!/usr/bin/env bash
set -euo pipefail

# Manifest columns: trait, genotype_prefix. Prefixes may be absolute or
# relative to GWAS_GENOTYPE_ROOT.
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT}"
PLINK2="${PLINK2_BIN:-plink2}"
MANIFEST="${GWAS_TRAIT_MANIFEST:?Set GWAS_TRAIT_MANIFEST}"
SOURCE_ROOT="${GWAS_GENOTYPE_ROOT:-${PROJECT_ROOT}/input/gwas_genotypes}"
OUT_ROOT="${GWAS_AUTOSOMAL_ROOT:-${PROJECT_ROOT}/work/gwas/autosomal}"
THREADS="${THREADS:-8}"
mkdir -p "${OUT_ROOT}"
while IFS=$'\t' read -r trait prefix _rest; do
  [[ -n "${trait}" && "${trait}" != "trait" ]] || continue
  [[ "${trait}" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Unsafe trait label: ${trait}" >&2; exit 2; }
  [[ "${prefix}" = /* ]] || prefix="${SOURCE_ROOT}/${prefix}"
  out="${OUT_ROOT}/${trait}"
  for suffix in bed bim fam; do [[ -s "${prefix}.${suffix}" ]] || { echo "Missing ${prefix}.${suffix}" >&2; exit 2; }; done
  [[ ! -e "${out}.bed" ]] || { echo "Output exists: ${out}.bed" >&2; exit 2; }
  "${PLINK2}" --bfile "${prefix}" --horse --autosome --make-bed --threads "${THREADS}" --out "${out}"
done < "${MANIFEST}"
