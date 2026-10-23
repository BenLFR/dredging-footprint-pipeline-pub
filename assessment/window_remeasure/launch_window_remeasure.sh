#!/bin/bash
# ============================================================================
# launch_window_remeasure.sh -- QF (injection depth) and DF tests, reading rules written beforehand
# (move of one column), under the window of Atwood et al. (2024), on the capped source.
# Models: launch_svr_ocim.sh for step7, launch_ocim_window.sh for the OCIM.
#   (A) step7 INJECTION_MODE=bottom x 2: same settings as step7_<lab>_an, except the mode
#   (B) lock_bottom.py x 2: afterany (A); B1-B3; writes BOTTOM_<lab>_PASS, otherwise exits with an error
#   (C) move_window.py: R0a, R0d, C1-C6 on step7_smooth_an; writes MOVE_PASS, otherwise exits with an error
#   (D) OCIM x 3 on hpc-08, window: bottom_<lab> afterok (B); shift_smooth afterok (C)
#   (E) read_remeasure.py: afterany (D); G0-G5 then the QF and DF verdicts
# If a lock fails, the dependent OCIM stays pending: scancel the remaining ids (ids_window_remeasure.txt).
# Pipeline untouched; outputs under BASE. DRY=1 by default.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
OCIM=${DFP_DATA:?set DFP_DATA}/configuration/ocim
SVR=${DFP_SCRATCH:?set DFP_SCRATCH}/svr_annuel_ocim_20261005
WIN=${DFP_SCRATCH:?set DFP_SCRATCH}/ocim_fenetre_20261009
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/remesure_fen_20261010
CODE="${DFP_REPO:?set DFP_REPO}/assessment/window_remeasure"
S7_MD5=8d49b9254847fe3d526ac9ef7932fb4c      # installed step7_export_jtrawl.R (M4, cap)
CO2_MD5=c800b20e729502650562aceef08c99ca     # co2model_batch.m with Section 6b (window)
GRIT_MD5=40382b287975dc5ef05910040015c846    # co2model_batch_grit.sh
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
EXCL_R="$EXCL,hpc-08.grit.ucsb.edu"          # R is broken on hpc-08
S7_BOT="INJECTION_MODE=bottom,REDIST_R_KM=500,REDIST_LAMBDA_KM=100,REDIST_K_MAX=16"   # = step7_<lab>_an except the mode
FENV="DREDGE_YEAR_START=2015,DREDGE_YEAR_END=2024,DREDGE_SCALE=0.1"                  # = launch_ocim_window.sh

echo "============================================================"
echo " QF and DF tests under the window    $(date)   DRY=$DRY"
echo "============================================================"
fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }

# ---- 1. checks ---------------------------------------------------------------
[ "$(md5sum "$DIR/step7_export_jtrawl.R" | cut -d' ' -f1)" = "$S7_MD5" ] \
  && ok "step7 installed ($S7_MD5)" || ko "installed step7 is not $S7_MD5"
[ "$(md5sum "$OCIM/co2model_batch.m" | cut -d' ' -f1)" = "$CO2_MD5" ] \
  && ok "co2model_batch.m installed = window ($CO2_MD5)" || ko "co2model_batch.m is not the window version"
[ "$(md5sum "$DIR/co2model_batch_grit.sh" | cut -d' ' -f1)" = "$GRIT_MD5" ] \
  && ok "co2model_batch_grit.sh unchanged" || ko "co2model_batch_grit.sh has changed"
for f in "$DIR/step7_export_jtrawl_grit.sh" "$CODE/lock_bottom.py" "$CODE/move_window.py" "$CODE/read_remeasure.py"; do
  [ -f "$f" ] && ok "present: $(basename "$f")" || ko "missing: $f"
done
( cd "$CODE" && md5sum lock_bottom.py move_window.py read_remeasure.py launch_window_remeasure.sh | sed 's/^/      /' )
for lab in smooth flag; do
  f="$SVR/cri_an_${lab}/cri_final_${lab}_an_20261005.parquet"
  [ -f "$f" ] && ok "$lab: annual C_ri input" || ko "missing: $f"
  n=$(ls "$SVR/step7_${lab}_an"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: 1 forcing under the cap (B1, B2, R0a)" || ko "$lab: $n forcings in step7_${lab}_an"
  for sc in baseline dredge_1x; do
    n=$(ls "$SVR/ocim_${lab}_an"/ocim_${sc}_*.mat 2>/dev/null | wc -l)
    [ "$n" -eq 1 ] && ok "$lab: run under the cap $sc (geometry)" || ko "$lab: $n ocim_${sc} in ocim_${lab}_an"
  done
  for sc in baseline dredge_1x dredge_10x dredge_100x; do
    n=$(ls "$WIN/ocim_${lab}_fen"/ocim_${sc}_*.mat 2>/dev/null | wc -l)
    [ "$n" -eq 1 ] && ok "$lab: window reference $sc" || ko "$lab: $n ocim_${sc} in ocim_${lab}_fen"
  done
done
python3 -c "import numpy, h5py, scipy.io" 2>/dev/null && ok "python3: numpy, h5py, scipy" || ko "python3 without numpy/h5py/scipy"
if [ -e "$BASE" ] && [ -n "$(ls -A "$BASE" 2>/dev/null)" ]; then ko "BASE not empty: $BASE"; fi
FREE=$(df -BG --output=avail "$(dirname "$BASE")" | tail -1 | tr -dc 0-9)
[ "${FREE:-0}" -ge 20 ] && ok "free space: ${FREE} GB" || ko "free space ${FREE:-?} GB < 20 GB"
[ "$fail" -eq 0 ] || { echo; echo "CHECKS FAILED -- nothing was submitted."; exit 1; }

# ---- 2. submission -----------------------------------------------------------
sb() {
  local ex="${EXCL_JOB:-$EXCL}"
  if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch --parsable --partition=$PART --exclude=$ex $*" >&2; echo "DRY"
  else sbatch --parsable --partition="$PART" --exclude="$ex" "$@" | cut -d';' -f1; fi
}
[ "$DRY" -eq 1 ] || mkdir -p "$BASE/logs"
IDS=""; DEP_E=""
for lab in smooth flag; do
  [ "$DRY" -eq 1 ] || mkdir -p "$BASE/step7_${lab}_bot" "$BASE/ocim_${lab}_bot"
  JS=$(EXCL_JOB="$EXCL_R" sb --job-name="rfS7bot_$lab" \
         --output="$BASE/logs/s7_${lab}_bot_%j.out" --error="$BASE/logs/s7_${lab}_bot_%j.err" \
         --export=ALL,$S7_BOT,OUT_DIR="$BASE/step7_${lab}_bot",CRI_FILE="$SVR/cri_an_${lab}/cri_final_${lab}_an_20261005.parquet" \
         "$DIR/step7_export_jtrawl_grit.sh")
  JV=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany:"$JS" --job-name="rfB_$lab" --qos=grit_normal \
         --mem=32G --cpus-per-task=1 --time=00:30:00 \
         --output="$BASE/logs/lock_bottom_${lab}.out" --error="$BASE/logs/lock_bottom_${lab}.err" \
         --wrap="cd $CODE && python3 lock_bottom.py $SVR $BASE $lab")
  JO=$(sb --dependency=afterok:"$JV" --job-name="rfCO2bot_$lab" --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_DIR="$BASE/step7_${lab}_bot",OCIM_OUT_DIR="$BASE/ocim_${lab}_bot",$FENV \
         "$DIR/co2model_batch_grit.sh")
  echo "  (A) $lab step7 bottom $JS ; (B) lock $JV ; (D) OCIM $JO"
  IDS="$IDS $JS $JV $JO"; DEP_E="$DEP_E:$JO"
done
[ "$DRY" -eq 1 ] || mkdir -p "$BASE/ocim_smooth_shift"
JC=$(EXCL_JOB="$EXCL_R" sb --job-name=rfDepl --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=01:00:00 \
       --output="$BASE/logs/move.out" --error="$BASE/logs/move.err" \
       --wrap="cd $CODE && python3 move_window.py $SVR $BASE")
JO=$(sb --dependency=afterok:"$JC" --job-name=rfCO2shift --nodelist=hpc-08.grit.ucsb.edu \
       --export=ALL,FORCING_DIR="$BASE/step7_smooth_an_shift",OCIM_OUT_DIR="$BASE/ocim_smooth_shift",$FENV \
       "$DIR/co2model_batch_grit.sh")
echo "  (C) shift $JC ; (D) OCIM shift_smooth $JO"
IDS="$IDS $JC $JO"; DEP_E="$DEP_E:$JO"
JL=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany"$DEP_E" --job-name=rfRead --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=01:00:00 \
       --output="$BASE/logs/read_remeasure.out" --error="$BASE/logs/read_remeasure.err" \
       --wrap="cd $CODE && python3 read_remeasure.py")
echo "  (E) reading: $JL"; IDS="$IDS $JL"

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_window_remeasure.txt" <<EOT
# QF and DF tests under the window, submitted on $(date)   BASE=$BASE
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_window_remeasure.txt"
echo "  Read first: $BASE/logs/lock_bottom_*.out and move.out ; after OCIM: read_remeasure.out"
