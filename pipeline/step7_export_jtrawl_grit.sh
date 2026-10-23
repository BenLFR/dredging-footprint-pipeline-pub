#!/bin/bash
###############################################################################
#  SLURM – STEP 7: Export Jtrawl for OCIM2-48L
#  GRIT version
###############################################################################
#SBATCH --job-name=step7_jtrawl
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
# (GRIT header: --chdir=$DFP_REPO/pipeline ; #SBATCH lines cannot expand variables)
#SBATCH --output=step7_jtrawl_%j.out
#SBATCH --error=step7_jtrawl_%j.err
#SBATCH --exclude=hpc-08.grit.ucsb.edu

echo "=== STEP 7: EXPORT Jtrawl (job $SLURM_JOB_ID) ==="
echo "Node:  $SLURMD_NODENAME"
echo "Start: $(date)"

# R setup
export R_LIBS_USER=~/R/library

# Directories
mkdir -p ${DFP_REPO:?set DFP_REPO}/pipeline/logs
mkdir -p ${DFP_DATA:?set DFP_DATA}/output_V6

cd ${DFP_REPO:?set DFP_REPO}/pipeline || { echo "Directory missing"; exit 2; }

echo "Current directory: $(pwd)"
echo "R library path: $R_LIBS_USER"

# Pre-flight: ocim_cache.mat
OCIM_CACHE=${DFP_DATA:?set DFP_DATA}/configuration/ocim/ocim_cache.mat
if [ ! -f "$OCIM_CACHE" ]; then
    echo "ocim_cache.mat missing: $OCIM_CACHE"
    echo "   Run first: matlab -batch \"run('step7_extract_ocim_cache.m')\""
    exit 1
fi
echo "ocim_cache.mat: $(ls -lh "$OCIM_CACHE" | awk '{print $5}')"

# Pre-flight: cri_final_*.parquet
# OUT_DIR is both input (cri_final) AND output (jdredge): the R script
# reads it, so the preflight must look at THE SAME directory.
export OUT_DIR="${OUT_DIR:-${DFP_DATA:?set DFP_DATA}/output_V6}"
PUB_DIR=$(readlink -f "${DFP_DATA:?set DFP_DATA}/output_V6")
if [ "$(readlink -f "$OUT_DIR")" = "$PUB_DIR" ] && [ -z "$STEP7_ALLOW_PUBLISHED" ]; then
    echo "OUT_DIR points to the PUBLISHED outputs ($OUT_DIR) -- WRITING REFUSED."
    echo "   Declare OUT_DIR, or STEP7_ALLOW_PUBLISHED=1 if intended."
    exit 1
fi
echo "OUT_DIR  : $OUT_DIR"
# CRI_DIR + CRI_TAG (chain_step4_ocim.sh): when the time-stamped name
# of C_ri is not known at submission, it is resolved here -- and there must
# be only one, otherwise the choice would be silent.
if [ -z "$CRI_FILE" ] && [ -n "$CRI_DIR" ]; then
    : "${CRI_TAG:?CRI_DIR requires CRI_TAG (smooth or flag)}"
    N_CRI=$(ls "$CRI_DIR"/cri_final_"${CRI_TAG}"_*.parquet 2>/dev/null | wc -l)
    if [ "$N_CRI" -ne 1 ]; then
        echo "CRI_DIR=$CRI_DIR: $N_CRI file(s) cri_final_${CRI_TAG}_*.parquet, 1 expected -- REFUSED."
        exit 1
    fi
    export CRI_FILE=$(ls "$CRI_DIR"/cri_final_"${CRI_TAG}"_*.parquet)
fi
if [ -n "$CRI_FILE" ]; then echo "CRI_FILE: $CRI_FILE (declared)"; fi
CRI_FILES=$(find "$OUT_DIR" -name "cri_final_*.parquet" 2>/dev/null | head -5)
if [ -z "$CRI_FILES" ] && [ -z "$CRI_FILE" ]; then
    echo "No cri_final_*.parquet file in $OUT_DIR"
    echo "   Run Step 6 first"
    exit 1
fi
echo "cri_final files found:"
echo "$CRI_FILES"

# Pre-flight: constants.R
if [ ! -f constants.R ]; then
    echo "constants.R missing in $(pwd)"
    exit 1
fi
echo "constants.R present"

# Pre-flight: R.matlab package
Rscript --vanilla -e "if (!requireNamespace('R.matlab', quietly=TRUE)) stop('R.matlab not installed')" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "R.matlab package missing. Installing..."
    Rscript --vanilla -e ".libPaths('~/R/library'); install.packages('R.matlab', repos='https://cloud.r-project.org', lib='~/R/library')"
fi

echo ""
echo "Launching step7_export_jtrawl.R ..."

Rscript --vanilla -e "
  .libPaths('~/R/library')
  source('step7_export_jtrawl.R')
"

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo ""
    echo "STEP 7 completed successfully: $(date)"
    echo "Files generated:"
    ls -lh "$OUT_DIR"/jdredge_ocim2_48l_* 2>/dev/null || echo "   No jdredge file found"
else
    echo "STEP 7 failed (code $exit_code)"
    exit $exit_code
fi
