#!/bin/bash
#SBATCH --job-name=step2_NO_CROSS
#SBATCH --array=1-10
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=logs/step2_process_%A_%a.out
#SBATCH --error=logs/step2_process_%A_%a.err

echo "=== STEP 2: PROCESSING VESSEL $SLURM_ARRAY_TASK_ID ==="
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"

# Conservative memory configuration
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# CRITICAL geospatial configuration
export LAND_MASK_FILE=${DFP_DATA:?set DFP_DATA}/configuration/land_mask/land_polygons.shp
export LAND_MASK_BUFFER_M=0
export LAND_NEAR_COAST_KM=5

if [ -f "${DFP_DATA:?set DFP_DATA}/configuration/land_mask/land_polygons.shp" ]; then
    echo "Land mask found: $(ls -lh ${DFP_DATA:?set DFP_DATA}/configuration/land_mask/land_polygons.shp | awk '{print $5}')"
else
    echo "WARNING: Land mask NOT FOUND!"
fi

SPLIT_JOB_ID=${1:-"12990"}
mkdir -p ${DFP_REPO:?set DFP_REPO}/pipeline/logs

cd ${DFP_REPO:?set DFP_REPO}/pipeline

echo "Starting processing of vessel $SLURM_ARRAY_TASK_ID..."
echo "Source: ${DFP_DATA:?set DFP_DATA}/ais_split_${SPLIT_JOB_ID}/"
echo "LAND_MASK_FILE: $LAND_MASK_FILE"
echo "LAND_MASK_BUFFER_M: $LAND_MASK_BUFFER_M"
echo "LAND_NEAR_COAST_KM: $LAND_NEAR_COAST_KM"

# Tiling override for vessels with a worldwide footprint (too many tiles at 5°)
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 2 ]]; then
  export LAND_TILE_DEG=10
  export LAND_TILE_MAX_TILES=500
  echo ">>> Override tiling task $SLURM_ARRAY_TASK_ID: LAND_TILE_DEG=$LAND_TILE_DEG, LAND_TILE_MAX_TILES=$LAND_TILE_MAX_TILES"
fi

export SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SPLIT_JOB_ID=$SPLIT_JOB_ID
Rscript step2_process_vessel.R 2>&1

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "✅ Vessel $SLURM_ARRAY_TASK_ID processed successfully: $(date)"
else
    echo "❌ Error vessel $SLURM_ARRAY_TASK_ID: code $exit_code"
    exit $exit_code
fi

echo "📊 VESSEL SUMMARY $SLURM_ARRAY_TASK_ID:"
ls -lh ${DFP_DATA:?set DFP_DATA}/ais_split_${SPLIT_JOB_ID}/navire_*${SLURM_ARRAY_TASK_ID}_*_clean.rds 2>/dev/null || echo "Clean file not found" 