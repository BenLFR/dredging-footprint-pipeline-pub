#!/bin/bash
#SBATCH --job-name=step4_rerun
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --partition=emlab_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-13.grit.ucsb.edu
#SBATCH --output=step4_rerun_%j.out
#SBATCH --error=step4_rerun_%j.err
set -euo pipefail

echo "=== STEP 4 -- RERUN 2026-08-10 (lithology) ==="
echo "Job ${SLURM_JOB_ID} on $(hostname) -- start $(date)"

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# Input AND output in the replay directory. The published outputs
# ($DFP_DATA/output_V6) are never written: the step4 guard
# refuses to start if PIPE_OUT_DIR points to them.
# Can be overridden (chain_step4_ocim.sh).
: "${PIPE_IN_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu}";  export PIPE_IN_DIR
: "${PIPE_OUT_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu}"; export PIPE_OUT_DIR
echo "PL_SCHEME    = ${PL_SCHEME:-sala2021 (script default)}"
mkdir -p "$PIPE_OUT_DIR"

echo "PIPE_IN_DIR  = $PIPE_IN_DIR"
echo "PIPE_OUT_DIR = $PIPE_OUT_DIR"

# Lock: step3 must have produced EXACTLY one file, otherwise the mtime
# inference would pick one without telling which.
N=$(ls "$PIPE_IN_DIR"/AIS_data_core_preprocessed_V6_*.rds 2>/dev/null | wc -l)
echo "step3 files in the replay directory: $N"
if [ "$N" -eq 0 ]; then
  echo "ERROR: step3 produced nothing in $PIPE_IN_DIR"; exit 1
fi
if [ "$N" -gt 1 ]; then
  echo "ERROR: $N step3 files -- ambiguous. Set PIPE_AIS_FILE explicitly."
  ls -lh "$PIPE_IN_DIR"/AIS_data_core_preprocessed_V6_*.rds
  exit 1
fi
echo "Selected input: $(basename "$(ls "$PIPE_IN_DIR"/AIS_data_core_preprocessed_V6_*.rds)")"

CACHE_DIR="${DFP_DATA:?set DFP_DATA}/hubocean_cache"
if [ ! -d "$CACHE_DIR" ]; then
  echo "ERROR: HubOcean cache missing: $CACHE_DIR"; exit 1
fi

cd "${DFP_REPO:?set DFP_REPO}/pipeline"
Rscript step4_add_lithology_vNext.R

echo "=== END step4 $(date) ==="
ls -lh "$PIPE_OUT_DIR"/AIS_with_lithology_*.rds
