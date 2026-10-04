#!/usr/bin/env bash
set -euo pipefail

# Extract and phase regions listed in a manifest. Columns: region,chrom,start,end.
VCF="${SELECTION_VCF:?Set SELECTION_VCF}"
MANIFEST="${HAPLOTYPE_REGIONS:?Set HAPLOTYPE_REGIONS}"
BCFTOOLS="${BCFTOOLS_BIN:-bcftools}"; BEAGLE="${BEAGLE_JAR:?Set BEAGLE_JAR}"
OUT_ROOT="${HAPLOTYPE_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/haplotype}"
THREADS="${THREADS:-8}"; MEMORY="${BEAGLE_MEMORY:-16g}"
[[ ! -e "${OUT_ROOT}" ]] || { echo "Output exists: ${OUT_ROOT}" >&2; exit 2; }; mkdir -p "${OUT_ROOT}"
while IFS=$'\t' read -r region chrom start end; do
  [[ -n "${region}" && "${region}" != region ]] || continue
  safe="${region//[^A-Za-z0-9_.-]/_}"; raw="${OUT_ROOT}/${safe}.vcf.gz"; phased="${OUT_ROOT}/${safe}.phased.vcf.gz"
  "${BCFTOOLS}" view -r "${chrom}:${start}-${end}" "${VCF}" -Oz -o "${raw}"; "${BCFTOOLS}" index -t "${raw}"
  java -Xmx"${MEMORY}" -jar "${BEAGLE}" gt="${raw}" out="${OUT_ROOT}/${safe}.phased" nthreads="${THREADS}"
  "${BCFTOOLS}" index -t "${phased}"
done < "${MANIFEST}"
