#!/bin/bash
#SBATCH --job-name=merge_rerun
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=emlab_nodes,grit_nodes
#SBATCH --qos=grit_normal
#SBATCH --exclude=hpc-13.grit.ucsb.edu
#SBATCH --output=merge_rerun_%x_%j.out
#SBATCH --error=merge_rerun_%x_%j.err
set -euo pipefail

: "${DREDGE_LABEL:?DREDGE_LABEL must be smooth or flag}"
case "$DREDGE_LABEL" in smooth|flag) ;; *) echo "ERROR: DREDGE_LABEL='$DREDGE_LABEL'"; exit 1 ;; esac

export R_LIBS_USER="$HOME/R/library"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MAKE_TIFF="${MAKE_TIFF:-FALSE}"          # no GeoTIFF: not needed, and heavy

: "${MERGE_IN_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu/step5_${DREDGE_LABEL}}"; export MERGE_IN_DIR
: "${MERGE_OUT_DIR:=${DFP_DATA:?set DFP_DATA}/output_V6_rejeu/merge_${DREDGE_LABEL}}"; export MERGE_OUT_DIR
: "${PIPE_TILES:=${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg}"; export PIPE_TILES
# Table grid_id -> fraction left on site. Undefined
# = p_r applied to 100 % of the SVR.
: "${FI_FRAC_FILE:=}"; export FI_FRAC_FILE
echo "  FRAC: ${FI_FRAC_FILE:-(not set -- p_r on 100 % of the SVR)}"

echo "=== MERGE -- label $DREDGE_LABEL -- $(date) ==="
echo "  IN  : $MERGE_IN_DIR"
echo "  OUT: $MERGE_OUT_DIR"
N=$(ls "$MERGE_IN_DIR"/sar_*.parquet 2>/dev/null | wc -l)
echo "  available sar tiles: $N"
# Refuses any merge whose number of tiles differs from the EXPECTED number:
# a partial merge produces a wrong SAR without raising an error.
: "${MERGE_N_EXPECTED:=71}"
if [ "$N" -ne "$MERGE_N_EXPECTED" ]; then
  echo "ERROR: $N sar tiles out of $MERGE_N_EXPECTED expected -- MERGE REFUSED."
  echo "  A partial merge produces a wrong SAR without raising yr error."
  echo "  Tiles present:"
  ls "$MERGE_IN_DIR"/sar_*.parquet 2>/dev/null | sed 's/.*sar_//;s/\.parquet//' | sort -n | tr '\n' ' '
  echo
  echo "  To merge a partial set deliberately: sbatch --export=ALL,MERGE_N_EXPECTED=$N ..."
  exit 1
fi

mkdir -p "$MERGE_OUT_DIR"
cd "${DFP_REPO:?set DFP_REPO}/pipeline"
# The merge script that produced the published results
# (not the earlier optimised variant of the merge).
Rscript step5_merge_tiles_corr.R

echo "=== END $(date) ==="
ls -lh "$MERGE_OUT_DIR"
