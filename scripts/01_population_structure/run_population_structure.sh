#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "${SCRIPT_DIR}/run_pca.sh"
bash "${SCRIPT_DIR}/run_distance_tree.sh"
bash "${SCRIPT_DIR}/prepare_admixture_markers.sh"
export ADMIXTURE_BED="${ADMIXTURE_BED:-${ADMIXTURE_ROOT:-${PROJECT_ROOT:?Set PROJECT_ROOT}/work/admixture}/admixture_markers}"
bash "${SCRIPT_DIR}/run_admixture.sh"
