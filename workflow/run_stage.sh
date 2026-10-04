#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
usage() {
  cat <<'EOF'
Usage: bash workflow/run_stage.sh [--config FILE] STAGE

Stages:
  qc             genotype filtering and PCA inputs
  structure      PCA, VCF2Dis distance, FastME tree, and ADMIXTURE
  selection      FST/pi scans and joint candidate-region reconstruction
  association    autosomal rMVP GWAS and significant-locus summary
  annotation     VEP and GWAS/selection region annotation
  robustness     print allele-frequency, haplotype, and sensitivity entry points
EOF
}
CONFIG="${CONFIG_FILE:-${ROOT}/config/config.env}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) break;;
  esac
done
[[ -f "${CONFIG}" ]] && source "${CONFIG}"
stage="${1:-}"; [[ -n "${stage}" ]] || { usage >&2; exit 2; }
case "${stage}" in
  qc) exec bash "${ROOT}/scripts/00_snp_filter/run_qc.sh";;
  structure) exec bash "${ROOT}/scripts/01_population_structure/run_population_structure.sh";;
  selection)
    bash "${ROOT}/scripts/02_selection_scan/run_selection_scans.sh"
    exec python3 "${ROOT}/scripts/02_selection_scan/build_candidate_regions.py" --raw-root "${SELECTION_RAW_ROOT:?Set SELECTION_RAW_ROOT}" --gff3 "${GFF3:?Set GFF3}" --out-root "${SELECTION_RESULT_ROOT:?Set SELECTION_RESULT_ROOT}" --traits ${SELECTION_TRAITS:?Set SELECTION_TRAITS}
    ;;
  association) exec bash "${ROOT}/scripts/03_association/run_association_pipeline.sh";;
  annotation) exec bash "${ROOT}/scripts/04_annotation/run_annotation_pipeline.sh";;
  robustness)
    cat <<EOF
Allele frequency: python3 ${ROOT}/scripts/05_robustness/compute_allele_frequency.py
Effect-allele audit: python3 ${ROOT}/scripts/05_robustness/audit_effect_alleles.py
Haplotype phasing: bash ${ROOT}/scripts/05_robustness/phase_haplotype_regions.sh
Candidate association sensitivity: ${RSCRIPT_BIN:-Rscript} ${ROOT}/scripts/05_robustness/run_candidate_association_sensitivity.R
Selection robustness: bash ${ROOT}/scripts/05_robustness/run_selection_robustness.sh
EOF
    ;;
  *) echo "Unknown stage: ${stage}" >&2; usage >&2; exit 2;;
esac
