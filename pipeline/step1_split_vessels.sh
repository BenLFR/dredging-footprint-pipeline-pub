#!/bin/bash
#SBATCH --job-name=ais_split_V6
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00

#SBATCH --output=logs/step1_split_%j.out
#SBATCH --error=logs/step1_split_%j.err

# R configuration (GRIT)
module load R

echo "🚀 === STEP 1: SPLITTING BY VESSEL ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

# R configuration (Beluga)

export R_LIBS_USER=~/R/library
export AIS_INPUT_FILE=${DFP_DATA:?set DFP_DATA}/AIS_data/benjamin3_clean.csv

mkdir -p ${DFP_DATA:?set DFP_DATA}/ais_split_${SLURM_JOB_ID}
mkdir -p logs

cd ${DFP_REPO:?set DFP_REPO}/pipeline

echo "✅ Working directory: $(pwd)"
echo "✅ R file check: $(ls -la step1_split_vessels.R 2>/dev/null || echo 'FILE MISSING')"

echo "✅ Launching the split script..."
Rscript --vanilla -e "
.libPaths('~/.local/R/4.2.1/')
Sys.setenv(SLURM_JOB_ID = '$SLURM_JOB_ID')
source('${DFP_REPO:?set DFP_REPO}/pipeline/step1_split_vessels.R')  # run the repository copy
" 2>&1

echo "✅ Step 1 finished: $(date)"
echo "📁 Files saved in: ${DFP_DATA:?set DFP_DATA}/ais_split_${SLURM_JOB_ID}/"

echo "📊 SPLITTING SUMMARY:"
ls -lh ${DFP_DATA:?set DFP_DATA}/ais_split_${SLURM_JOB_ID}/ | grep ".fst"
wc -l ${DFP_DATA:?set DFP_DATA}/ais_split_${SLURM_JOB_ID}/*.fst 2>/dev/null || echo "Error: no .fst file found" 