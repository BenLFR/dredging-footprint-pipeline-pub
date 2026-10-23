#!/bin/bash
# ============================================================================
# launch_threshold_step2.sh -- Threshold +/- 0.005 (Label B), STEP 2 (rule written beforehand).
# Requires BASE/PREPARATION_PASS (P0, P1) and BASE/P2_PASS. For each arm t0, tp, tm:
#   step5 flag as an array over the tiles of the arm (PIPE_RUN_DIR = BASE/<arm>)
#     -> merge WITH the original frac table -> step6 (BASE/<arm>/cri_pr)       : reading 1, and P3 for t0
#     -> merge WITHOUT table (FI_FRAC_FILE missing) -> step6 (BASE/<arm>_sansfrac): C0 of reading 2
# then read_threshold.R all (afterany: a failed link is REPORTED by the reading, not hidden).
# Load: ~71 step5 tasks per arm (~213); MaxSubmitPU = 300 (chain_step4_ocim.sh): do not launch at the same
# time as another array.
# DRY=1 by default.
# ============================================================================
set -euo pipefail
unset FI_FRAC_FILE
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
REF="${DFP_DATA:?set DFP_DATA}/output_V6_rejeu"
RUN=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/seuil_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/threshold"
PART=grit_nodes
# hpc-08 excluded: R fails there when loading data.table and
# Rcpp ("undefined symbol"); hpc-08 is reserved for MATLAB (OCIM).
EXCL=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
RSC="$HOME/R/library"
THROTTLE=20

[ -f "$BASE/PREPARATION_PASS" ] || { echo "ABORT: no PREPARATION_PASS (read logs/prepare.out)"; exit 1; }
[ -f "$BASE/P2_PASS" ] || { echo "ABORT: no P2_PASS (read logs/P2.out)"; exit 1; }
sb() {
  if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch --parsable --partition=$PART --exclude=$EXCL $*" >&2; echo "DRY"
  else sbatch --parsable --partition="$PART" --exclude="$EXCL" "$@" | cut -d';' -f1; fi
}
IDS=""; DEP=""
for b in t0 tp tm; do
  T=$(cat "$BASE/$b/tuiles_$b.txt")
  N=$(echo "$T" | tr ',' '\n' | grep -c .)
  n=$(ls "$BASE/$b"/AIS_with_lithology*.rds 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] || { echo "ABORT: $n rds in $BASE/$b (1 expected)"; exit 1; }
  [ -e "$BASE/$b/step5_flag" ] && { echo "ABORT: $BASE/$b/step5_flag already exists"; exit 1; }
  J5=$(sb --array="$T%$THROTTLE" --job-name="threshold_s5_$b" \
         --output="$BASE/logs/s5_${b}_%A_%a.out" --error="$BASE/logs/s5_${b}_%A_%a.err" \
         --export=ALL,DREDGE_LABEL=flag,PIPE_RUN_DIR="$BASE/$b",PIPE_LITHO_FILE="$BASE/$b/AIS_with_lithology_seuil_$b.rds" \
         "$DIR/step5_batch.sbatch")
  JA=$(sb --dependency=afterok:"$J5" --job-name="threshold_m_${b}_with" \
         --output="$BASE/logs/m_${b}_with_%j.out" --error="$BASE/logs/m_${b}_with_%j.err" \
         --export=ALL,DREDGE_LABEL=flag,MERGE_IN_DIR="$BASE/$b/step5_flag",MERGE_OUT_DIR="$BASE/$b/merge_flag_pr",FI_FRAC_FILE="$REF/decharge/fraction_sur_place_flag.csv",MERGE_N_EXPECTED=$N \
         "$DIR/merge_rerun.sbatch")
  J6A=$(sb --dependency=afterok:"$JA" --job-name="threshold_s6_${b}_with" \
         --output="$BASE/logs/s6_${b}_with_%j.out" --error="$BASE/logs/s6_${b}_with_%j.err" \
         --export=ALL,DREDGE_LABEL=flag,PIPE_RUN_DIR="$BASE/$b" "$DIR/step6_pr.sbatch")
  # WITHOUT table: FI_FRAC_FILE is not passed (unset above); lock M reads the merge log
  JS=$(sb --dependency=afterok:"$J5" --job-name="threshold_m_${b}_without" \
         --output="$BASE/logs/m_${b}_without_%j.out" --error="$BASE/logs/m_${b}_without_%j.err" \
         --export=ALL,DREDGE_LABEL=flag,MERGE_IN_DIR="$BASE/$b/step5_flag",MERGE_OUT_DIR="$BASE/${b}_sansfrac/merge_flag_pr",MERGE_N_EXPECTED=$N \
         "$DIR/merge_rerun.sbatch")
  J6S=$(sb --dependency=afterok:"$JS" --job-name="threshold_s6_${b}_without" \
         --output="$BASE/logs/s6_${b}_without_%j.out" --error="$BASE/logs/s6_${b}_without_%j.err" \
         --export=ALL,DREDGE_LABEL=flag,PIPE_RUN_DIR="$BASE/${b}_sansfrac" "$DIR/step6_pr.sbatch")
  echo "  $b ($N tiles): step5 $J5 -> with: merge $JA -> step6 $J6A ; without: merge $JS -> step6 $J6S"
  IDS="$IDS $J5 $JA $J6A $JS $J6S"; DEP="$DEP:$J6A:$J6S"
  [ "$DRY" -eq 1 ] || mkdir -p "$BASE/${b}_sansfrac"
done
JL=$(sb --dependency=afterany"$DEP" --job-name="threshold_read" --qos=grit_normal --mem=16G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/read_threshold.out" --error="$BASE/logs/read_threshold.err" \
       --wrap="cd $CODE && R_LIBS_USER=$RSC Rscript read_threshold.R $BASE $RUN all")
echo "  reading: $JL"
IDS="$IDS $JL"
if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_step2.txt" <<EOT
# Threshold, step 2, submitted on $(date)   BASE=$BASE
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_step2.txt"
