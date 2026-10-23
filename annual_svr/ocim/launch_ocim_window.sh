#!/bin/bash
# ============================================================================
# launch_ocim_window.sh -- OCIM following the protocol of Atwood et al. (2024,
# §2.2, p. 3), rule written beforehand: dredging 2015-2024 at the annual mass Mijk (decadal total x 0.1), then zero,
# runs up to 2100; control, 1x, 10x, 100x.
# Forcing: that of the capped runs, SVR/step7_<lab>_an, unchanged.
# Model: co2model_batch.m with Section 6b (window), installed beforehand on the cluster.
#   (E) OCIM x 2 on hpc-08; ~4 h of control + ~1 h per dredging scenario (control state taken over in 2015)
#   (F) read_window.py: afterany (E); F0-F4 then the substance
# Outputs under BASE, never in the published outputs. DRY=1 by default.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
OCIM=${DFP_DATA:?set DFP_DATA}/configuration/ocim
SVR=${DFP_SCRATCH:?set DFP_SCRATCH}/svr_annuel_ocim_20261005
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/ocim_fenetre_20261009
CODE="${DFP_REPO:?set DFP_REPO}/annual_svr/ocim"
CO2_MD5=c800b20e729502650562aceef08c99ca     # co2model_batch.m with Section 6b
GRIT_MD5=40382b287975dc5ef05910040015c846    # co2model_batch_grit.sh, unchanged
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
EXCL_R="$EXCL,hpc-08.grit.ucsb.edu"
WIN="DREDGE_YEAR_START=2015,DREDGE_YEAR_END=2024,DREDGE_SCALE=0.1"

echo "============================================================"
echo " OCIM, Atwood protocol (2015-2024)    $(date)   DRY=$DRY"
echo "============================================================"
fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }

# ---- 1. checks ---------------------------------------------------------------
[ "$(md5sum "$OCIM/co2model_batch.m" | cut -d' ' -f1)" = "$CO2_MD5" ] \
  && ok "co2model_batch.m installed = window ($CO2_MD5)" || ko "co2model_batch.m is not the window version"
[ "$(md5sum "$DIR/co2model_batch_grit.sh" | cut -d' ' -f1)" = "$GRIT_MD5" ] \
  && ok "co2model_batch_grit.sh unchanged" || ko "co2model_batch_grit.sh has changed"
[ -f "$CODE/read_window.py" ] && ok "present: read_window.py" || ko "missing: $CODE/read_window.py"
declare -A FORC
for lab in smooth flag; do
  n=$(ls "$SVR/step7_${lab}_an"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: 1 capped forcing" || ko "$lab: $n forcings in step7_${lab}_an (1 expected)"
  FORC[$lab]=$(ls "$SVR/step7_${lab}_an"/jdredge_ocim2_48l_*.mat 2>/dev/null | head -1)
  for sc in baseline dredge_1x; do
    n=$(ls "$SVR/ocim_${lab}_an"/ocim_${sc}_*.mat 2>/dev/null | wc -l)
    [ "$n" -eq 1 ] && ok "$lab: capped run $sc (F0, F2)" || ko "$lab: $n ocim_${sc} in ocim_${lab}_an"
  done
done
python3 -c "import numpy, h5py" 2>/dev/null && ok "python3: numpy, h5py" || ko "python3 without numpy/h5py"
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
IDS=""; DEP_F=""
for lab in smooth flag; do
  [ "$DRY" -eq 1 ] || mkdir -p "$BASE/ocim_${lab}_fen"
  JO=$(sb --job-name="winCO2_$lab" --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_FILE="${FORC[$lab]}",OCIM_OUT_DIR="$BASE/ocim_${lab}_fen",$WIN \
         "$DIR/co2model_batch_grit.sh")
  echo "  (E) $lab OCIM: $JO   forcing ${FORC[$lab]}"; IDS="$IDS $JO"; DEP_F="$DEP_F:$JO"
done
JL=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany"$DEP_F" --job-name=winRead --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/read_window.out" --error="$BASE/logs/read_window.err" \
       --wrap="cd $CODE && python3 read_window.py")
echo "  (F) reading: $JL"; IDS="$IDS $JL"

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_ocim_window.txt" <<EOT
# OCIM, Atwood protocol, submitted on $(date)   BASE=$BASE
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_ocim_window.txt"
echo "  OCIM logs: co2model_<id>.out in the submission directory ; reading: $BASE/logs/read_window.out"
