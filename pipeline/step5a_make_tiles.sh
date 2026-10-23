#!/bin/bash
#SBATCH --job-name=step5a_maketiles
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=emlab_nodes,grit_nodes
#SBATCH --output=logs/step5a_%j.out
#SBATCH --error=logs/step5a_%j.err

export R_LIBS_USER=~/R/library
export LD_LIBRARY_PATH=$HOME/lib:${LD_LIBRARY_PATH:-}

PIPELINE_DIR=${DFP_REPO:?set DFP_REPO}/pipeline
mkdir -p "$PIPELINE_DIR/logs"
cd "$PIPELINE_DIR"

echo "$(date) -- step5a make_tiles on $(hostname)"
Rscript step5_make_tiles.R

if [ $? -eq 0 ]; then
  echo "$(date) -- step5a completed successfully"
  ls -lh ${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg
else
  echo "$(date) -- step5a FAILED"
  exit 1
fi
