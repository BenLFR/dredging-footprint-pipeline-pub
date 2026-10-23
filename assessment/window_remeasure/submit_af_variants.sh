#!/bin/bash
# Submits extract_af_variants.py (Figure 6). Launch: cat submit_af_variants.sh | ssh grit-hpc 'bash -s'
# Modelled on the launcher of extract_af_annual.py: same partition, excluded nodes and memory.
CODE="${DFP_REPO:?set DFP_REPO}/assessment/window_remeasure"
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/remesure_fen_20261010
OUT=$BASE/logs
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu,hpc-08.grit.ucsb.edu
cd "$CODE" || exit 1
for sd in ocim_smooth_bot ocim_flag_bot ocim_smooth_shift; do
  ls "$BASE/$sd"/ocim_baseline_*.mat "$BASE/$sd"/ocim_dredge_1x_*.mat >/dev/null || { echo "MISSING: $BASE/$sd"; exit 1; }
done
md5sum extract_af_variants.py
[ -e "$OUT/af_variantes_20261016.csv" ] && { echo "ALREADY PRESENT: $OUT/af_variantes_20261016.csv"; exit 1; }
sbatch --parsable --partition=grit_nodes --exclude="$EXCL" --job-name=afVar --qos=grit_normal \
  --mem=32G --cpus-per-task=1 --time=00:30:00 \
  --output="$OUT/af_variants.out" --error="$OUT/af_variants.err" \
  --wrap="cd $CODE && python3 extract_af_variants.py $OUT"
