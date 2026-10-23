#!/bin/bash
###############################################################################
#  SLURM – CO2MODEL_BATCH: OCIM2-48L dredging perturbation (MATLAB)
###############################################################################
#SBATCH --job-name=co2model
#SBATCH --partition=emlab_nodes
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
# (GRIT header: --chdir=$DFP_DATA/configuration/ocim ; #SBATCH lines cannot expand variables)
#SBATCH --output=co2model_%j.out
#SBATCH --error=co2model_%j.err

echo "=== CO2MODEL_BATCH (job $SLURM_JOB_ID) ==="
echo "Node:  $SLURMD_NODENAME"
echo "Start: $(date)"

cd ${DFP_DATA:?set DFP_DATA}/configuration/ocim || { echo "Directory missing"; exit 2; }

# Pre-flight: forcing file
# FORCING_FILE: DECLARED forcing. Without it, "the most recent" is taken, which
# silently chooses between two labels.
FORCING_DIR="${FORCING_DIR:-${DFP_DATA:?set DFP_DATA}/output_V6}"
if [ -n "$FORCING_FILE" ]; then
    FORCING="$FORCING_FILE"
    if [ ! -f "$FORCING" ]; then echo "FORCING_FILE not found: $FORCING"; exit 1; fi
    echo "DECLARED forcing"
else
    FORCING=$(ls -t "$FORCING_DIR"/jdredge_ocim2_48l_*.mat 2>/dev/null | head -1)
fi
if [ -z "$FORCING" ]; then
    echo "No jdredge_ocim2_48l_*.mat in $FORCING_DIR"
    echo "   Run Step 7 first"
    exit 1
fi
export OCIM_OUT_DIR="${OCIM_OUT_DIR:-${DFP_DATA:?set DFP_DATA}/output_V6}"
PUB_DIR=$(readlink -f ${DFP_DATA:?set DFP_DATA}/output_V6)
if [ "$(readlink -f "$OCIM_OUT_DIR")" = "$PUB_DIR" ] && [ -z "$OCIM_ALLOW_PUBLISHED" ]; then
    echo "OCIM_OUT_DIR points to the PUBLISHED outputs -- WRITE REFUSED."
    exit 1
fi
echo "OCIM_OUT_DIR: $OCIM_OUT_DIR"
echo "Forcing: $FORCING"

# Pre-flight: OCIM data
for f in OCIM2_48L_CTL.mat woa09si.mat woa09po4.mat netemission.txt; do
    if [ ! -f "$f" ]; then
        echo "Missing file: $f"
        exit 1
    fi
done
echo "All dependencies present"

echo ""
echo "Launching co2model_batch.m ..."

# matlab is not in the PATH of the login nodes; it is on hpc-08.
MATLAB_BIN="${MATLAB_BIN:-/home/matlab/current/bin/matlab}"
if [ ! -x "$MATLAB_BIN" ]; then echo "MATLAB not found: $MATLAB_BIN"; exit 1; fi
"$MATLAB_BIN" -batch "forcing_file='$FORCING'; run('co2model_batch.m')"

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo ""
    echo "CO2MODEL finished successfully: $(date)"
    echo "Files generated:"
    ls -lh "$OCIM_OUT_DIR"/ocim_*.mat 2>/dev/null || echo "   No ocim file found"
else
    echo "CO2MODEL failed (code $exit_code)"
    exit $exit_code
fi
