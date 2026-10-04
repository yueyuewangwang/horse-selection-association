#!/usr/bin/env bash
set -euo pipefail

# Annotate a PLINK BIM marker list with an offline VEP cache.
PROJECT_ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT in config/config.env}"
BIM="${VEP_INPUT_BIM:?Set VEP_INPUT_BIM}"
VEP="${VEP_BIN:-vep}"
BGZIP="${BGZIP_BIN:-bgzip}"
TABIX="${TABIX_BIN:-tabix}"
OUT="${VEP_VCF:?Set VEP_VCF}"
CACHE="${VEP_CACHE_DIR:?Set VEP_CACHE_DIR}"
FASTA="${VEP_FASTA:?Set VEP_FASTA}"
THREADS="${THREADS:-8}"
mkdir -p "$(dirname "${OUT}")"
if [[ -e "${OUT}" && "${FORCE:-0}" != "1" ]]; then
  echo "Output exists: ${OUT} (set FORCE=1 to replace)" >&2
  exit 2
fi
tmp="${OUT%.vcf.gz}.input.vcf"
{
  printf '##fileformat=VCFv4.2\n'
  printf '##reference=reference\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$4,$2,$6,$5,".",".","."}' "${BIM}"
} > "${tmp}"
"${VEP}" --input_file "${tmp}" --output_file "${OUT}" --format vcf --vcf \
  --cache --offline --dir_cache "${CACHE}" --fasta "${FASTA}" \
  --fork "${THREADS}" --buffer_size 10000 --compress_output bgzip --force_overwrite
"${TABIX}" -f -p vcf "${OUT}"
rm -f "${tmp}"
echo "VEP annotation: ${OUT}"
