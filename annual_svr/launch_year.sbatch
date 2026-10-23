#!/bin/bash
# launch_year.sbatch — reruns worker_year.R (copy of step5_tile_worker_batch.R with the
# SAR split by year) on the tiles of lock L1. Environment copied from pipeline/step5_batch.sbatch
# as chain_step4_ocim.sh launched it (DREDGE_LABEL + PIPE_RUN_DIR=output_V6_sala2021_r2). Only changes:
# the script, its folder ($DFP_REPO/annual_svr, with a copy of constants.R), and PIPE_OUT_DIR (new folder).
# Partition grit_nodes only: emlab_nodes = sequoia, which has no hpc-scratch.
#SBATCH --job-name=svr_year
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=grit_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu
#SBATCH --output=annee_%x_%A_%a.out
#SBATCH --error=annee_%x_%A_%a.err
set -euo pipefail
: "${DREDGE_LABEL:?DREDGE_LABEL must be smooth or flag}"
case "$DREDGE_LABEL" in smooth|flag) ;; *) echo "DREDGE_LABEL invalid"; exit 1 ;; esac

export R_LIBS_USER="$HOME/R/library"
export LD_LIBRARY_PATH="$HOME/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TILE_ID=$SLURM_ARRAY_TASK_ID
export PIPE_IN_DIR=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r2                  # read-only
export PIPE_OUT_DIR=${DFP_SCRATCH:?set DFP_SCRATCH}/svr_annuel_20261005/step5_${DREDGE_LABEL}  # new folder
export PIPE_TILES="${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg"                     # read-only
export SHIP_SPECS="${DFP_REPO:?set DFP_REPO}/config/ship_specs.yaml"
export PIPE_MONO_FRACTION=0.0
mkdir -p "$PIPE_OUT_DIR"
echo "$(date) -- tile $TILE_ID, label $DREDGE_LABEL, on $(hostname)"
cd "${DFP_REPO:?set DFP_REPO}/annual_svr"
Rscript worker_year.R "$TILE_ID"
echo "$(date) -- tile $TILE_ID OK ($DREDGE_LABEL)"
