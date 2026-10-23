#!/bin/bash
# BATCH variant. Calls step5_tile_worker_batch.R.
# PIPE_OUT_SUFFIX allows writing elsewhere for the equivalence test.
# Walltime 8 h: about 5 min per tile is expected.
#SBATCH --job-name=step5_rerun
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=emlab_nodes,grit_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-13.grit.ucsb.edu
#SBATCH --output=step5_rerun_%x_%A_%a.out
#SBATCH --error=step5_rerun_%x_%A_%a.err
set -euo pipefail

# DREDGE_LABEL is passed by --export at sbatch time (smooth | flag).
: "${DREDGE_LABEL:?DREDGE_LABEL must be smooth or flag}"
case "$DREDGE_LABEL" in
  smooth|flag) ;;
  *) echo "ERROR: DREDGE_LABEL='$DREDGE_LABEL' (expected: smooth|flag)"; exit 1 ;;
esac

export R_LIBS_USER="$HOME/R/library"
export LD_LIBRARY_PATH="$HOME/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TILE_ID=$SLURM_ARRAY_TASK_ID

# PIPE_RUN_DIR can be overridden (chain_step4_ocim.sh).
# Input (AIS_with_lithology) and outputs are read/written there.
: "${PIPE_RUN_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu}"
export PIPE_IN_DIR="$PIPE_RUN_DIR"
export PIPE_OUT_DIR="$PIPE_RUN_DIR/step5_${DREDGE_LABEL}${PIPE_OUT_SUFFIX:-}"
export PIPE_TILES="${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg"   # read-only
export SHIP_SPECS="${DFP_REPO:?set DFP_REPO}/config/ship_specs.yaml"

# OPTION C: rho = 0 -> all installed pipes in service.
# For the sensitivity analysis, resubmit with --export=...,PIPE_MONO_FRACTION=0.3
export PIPE_MONO_FRACTION="${PIPE_MONO_FRACTION:-0.0}"

mkdir -p "$PIPE_OUT_DIR"

if [ "$TILE_ID" -eq 1 ]; then
  echo "=== STEP 5 -- RERUN 2026-08-10 -- label=$DREDGE_LABEL ==="
  echo "PIPE_IN_DIR        = $PIPE_IN_DIR"
  echo "PIPE_OUT_DIR       = $PIPE_OUT_DIR"
  echo "PIPE_MONO_FRACTION = $PIPE_MONO_FRACTION  (0 = option C)"
  echo "SHIP_SPECS         = $SHIP_SPECS"
fi
echo "$(date) -- tile $TILE_ID, label $DREDGE_LABEL, on $(hostname)"

[ -f "$PIPE_TILES" ] || { echo "ERROR: tiles missing: $PIPE_TILES"; exit 1; }

cd "${DFP_REPO:?set DFP_REPO}/pipeline"
Rscript step5_tile_worker_batch.R "$TILE_ID"

echo "$(date) -- tile $TILE_ID OK ($DREDGE_LABEL)"
