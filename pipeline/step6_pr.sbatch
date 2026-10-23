#!/bin/bash
#SBATCH --job-name=step6_pr
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=emlab_nodes,grit_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-13.grit.ucsb.edu
#SBATCH --output=step6_pr_%x_%j.out
#SBATCH --error=step6_pr_%x_%j.err
set -euo pipefail
# C_ri on the grids of block 3.4 (p_r restricted to the on-site fraction).
: "${DREDGE_LABEL:?smooth or flag}"
export R_LIBS_USER="$HOME/R/library"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# PIPE_RUN_DIR can be overridden (chain_step4_ocim.sh).
: "${PIPE_RUN_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu}"
: "${CRI_OUT_DIR:=$PIPE_RUN_DIR/cri_pr}"; export CRI_OUT_DIR
export CRI_TAG="$DREDGE_LABEL"

G=$(ls -t "$PIPE_RUN_DIR/merge_${DREDGE_LABEL}_pr"/fi_grid_*.parquet | head -1)
echo "=== step6 $DREDGE_LABEL -- $(date) ==="
echo "  grid: $G"
echo "  output: $CRI_OUT_DIR (tag $CRI_TAG)"
cd "${DFP_REPO:?set DFP_REPO}/pipeline"
Rscript step6_calculate_cri_corr.R "$G"
echo "=== END $(date) ==="
ls -lh "$CRI_OUT_DIR"
