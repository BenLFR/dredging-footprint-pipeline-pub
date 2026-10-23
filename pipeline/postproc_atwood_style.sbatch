#!/bin/bash
#SBATCH --job-name=postproc_atwood
#SBATCH --partition=emlab_nodes
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
# (GRIT header: --chdir=$DFP_DATA ; #SBATCH lines cannot expand variables)
#SBATCH --output=postproc_atwood_%j.out
#SBATCH --error=postproc_atwood_%j.err

set -euo pipefail

echo "=== OCIM POSTPROC ATWOOD STYLE ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Node:   ${SLURMD_NODENAME:-NA}"
echo "Start:  $(date)"

OUTPUT_ROOT="${OUTPUT_ROOT:-${DFP_DATA:?set DFP_DATA}/output_V6}"
OCIM_DIR="${OCIM_DIR:-${DFP_DATA:?set DFP_DATA}/configuration/ocim}"
POSTPROC_DIR="${POSTPROC_DIR:-${OUTPUT_ROOT}/postproc_atwood_style}"
SCRIPT_PATH="${SCRIPT_PATH:-${SLURM_SUBMIT_DIR}/pipeline/postproc_atwood_style.py}"
REQ_PATH="${REQ_PATH:-${SLURM_SUBMIT_DIR}/pipeline/requirements_postproc_atwood.txt}"
VENV_PATH="${VENV_PATH:-${DFP_DATA:?set DFP_DATA}/venvs/ocim-postproc}"

mkdir -p $HOME/logs
mkdir -p "${POSTPROC_DIR}"

if command -v module >/dev/null 2>&1; then
  module purge || true
  module load python/3.11 || module load python || true
fi

if [ ! -d "${VENV_PATH}" ]; then
  python3 -m venv "${VENV_PATH}"
fi

source "${VENV_PATH}/bin/activate"

if [ ! -f "${VENV_PATH}/.postproc_deps_ok" ]; then
  python -m pip install --upgrade pip
  python -m pip install -r "${REQ_PATH}"
  touch "${VENV_PATH}/.postproc_deps_ok"
fi

python "${SCRIPT_PATH}" \
  --output-root "${OUTPUT_ROOT}" \
  --ocim-dir "${OCIM_DIR}" \
  --postproc-dir "${POSTPROC_DIR}" \
  --maps-scenario dredge_1x \
  --dpi 300

echo ""
echo "Done: $(date)"
echo "Outputs in: ${POSTPROC_DIR}"
ls -lh "${POSTPROC_DIR}" || true
