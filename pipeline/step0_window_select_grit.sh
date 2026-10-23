#!/bin/bash
###############################################################################
#  SLURM – STEP 0bis ENHANCED: Selection of the optimal time window
#  GRIT version
###############################################################################
#SBATCH --job-name=ais_window_V6_enhanced
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=logs/step0_window_enhanced_%j.out
#SBATCH --error=logs/step0_window_enhanced_%j.err

echo "🔍 === STEP 0bis ENHANCED: OPTIMAL WINDOW (job $SLURM_JOB_ID) ==="
echo "Node:  $SLURMD_NODENAME"
echo "Start: $(date)"

# Modules / librairies R (GRIT configuration)
module load R
export R_LIBS_USER=~/R/library

mkdir -p logs
cd ${DFP_REPO:?set DFP_REPO}/pipeline || { echo "❌ Missing directory"; exit 2; }

echo "✅  Current directory: $(pwd)"
echo "✅  R library path: $R_LIBS_USER"
echo "✅  Launching step0_core_window_enhanced.R ..."

# Call of the R interpreter with .libPaths override
Rscript --vanilla -e "
  .libPaths('~/R/library')
  source('step0_core_window_enhanced.R')
"

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo "✅  STEP 0bis ENHANCED completed successfully: $(date)"
    echo "📄  Selection file: ${DFP_DATA:?set DFP_DATA}/output_V6/core_window_report.md"
    echo "📊  Full analysis: ${DFP_DATA:?set DFP_DATA}/output_V6/all_windows_analysis.csv"
    echo "🏆  Top 20 windows: ${DFP_DATA:?set DFP_DATA}/output_V6/top_20_windows.csv"
    echo "📈  Annual statistics: ${DFP_DATA:?set DFP_DATA}/output_V6/annual_statistics.csv"
    echo "🗺️  Coverage matrix: ${DFP_DATA:?set DFP_DATA}/output_V6/coverage_matrix.csv"
else
    echo "❌  STEP 0bis ENHANCED failed (code $exit_code)"
fi
