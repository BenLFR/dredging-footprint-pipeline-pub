#!/bin/bash
#SBATCH --job-name=step3_r6
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --partition=emlab_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-13.grit.ucsb.edu
#SBATCH --output=step3_r6_%j.out
#SBATCH --error=step3_r6_%j.err
set -euo pipefail

echo "=== STEP 3 REPLAY r6 (copy of the 2026-08-07 step3 replay launcher, configurable output) ==="
echo "Job ${SLURM_JOB_ID} on $(hostname) -- start $(date)"

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# Outputs go to a NEW directory. The published outputs ($DFP_DATA/output_V6)
# are never touched -- and the guard in step3_merge_final.R would
# refuse to start if PIPE_OUT_DIR pointed to them.
export PIPE_OUT_DIR="${PIPE_OUT_DIR:-${DFP_DATA:?set DFP_DATA}/output_V6_r6_step3}"
case "$(readlink -f "$PIPE_OUT_DIR")" in */output_V6|*/output_V6_rejeu) echo "ERROR: PIPE_OUT_DIR = published output or August run"; exit 1;; esac
if [ -d "$PIPE_OUT_DIR" ] && [ -n "$(ls -A "$PIPE_OUT_DIR")" ]; then echo "ERROR: $PIPE_OUT_DIR not empty"; exit 1; fi
FREE=$(df -BG --output=avail "$HOME" | tail -1 | tr -dc 0-9)
[ "${FREE:-0}" -ge 3 ] || { echo "ERROR: ${FREE:-?} GB free under $HOME (3 required)"; exit 1; }
mkdir -p "$PIPE_OUT_DIR" "${DFP_REPO:?set DFP_REPO}/pipeline/logs"
echo "PIPE_OUT_DIR = $PIPE_OUT_DIR"
df -h "$HOME" | tail -1

cd "${DFP_REPO:?set DFP_REPO}/pipeline"

export SPLIT_JOB_ID="${SPLIT_JOB_ID:-12990}"
export RESULTS_DIR="${RESULTS_DIR:-${DFP_DATA:?set DFP_DATA}/ais_results_14542}"
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"

CLEAN_COUNT=$(ls "$RESULTS_DIR"/*_clean.rds 2>/dev/null | wc -l)
echo "Clean input files: $CLEAN_COUNT (from $RESULTS_DIR)"
if [ "$CLEAN_COUNT" -eq 0 ]; then
  echo "ERROR: no *_clean.rds in $RESULTS_DIR"
  exit 1
fi

Rscript step3_merge_final.R

echo "=== END $(date) ==="
ls -lh "$PIPE_OUT_DIR"
