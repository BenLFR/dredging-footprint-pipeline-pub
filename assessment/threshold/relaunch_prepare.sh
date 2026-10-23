#!/bin/bash
# ============================================================================
# relaunch_prepare.sh -- Threshold +/- 0.005: rerun of the prepare ALONE, with the revised rule P1
# (revised before any arm was run; explicit decision).
# P2 passed at the 2nd attempt (P2_PASS): it is not redone. The traces of the failed prepare are set aside, not
# deleted. Then: read logs/prepare.out, then launch_threshold_step2.sh.
# DRY=1 by default.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/seuil_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/threshold"
PART=grit_nodes
# EXCL line taken from launch_threshold_step1.sh (hpc-08 excluded: R does not load data.table there)
EXCL=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
RSC="$HOME/R/library"
PREP_MD5=8075e529f63ff52183d3bb724f6e5c9e       # threshold_prepare.R with revised P1
READ_MD5=34bacc4c765d0e54d5170e3df8fd306a       # read_threshold.R with P1b in the verdict
WORKER_MD5=cd6e9c4d6932bb249df8af1c6a919431
MRS_MD5=4dbc5b1d4229653698a83dcb7d5998dd

fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }
[ -f "$BASE/P2_PASS" ] && ok "P2_PASS present" || ko "P2_PASS missing"
[ -e "$BASE/PREPARATION_PASS" ] && ko "PREPARATION_PASS already exists" || ok "no PREPARATION_PASS"
for b in t0 tp tm; do [ -e "$BASE/$b" ] && ko "$BASE/$b already exists" || ok "no folder $b"; done
[ "$(md5sum "$CODE/threshold_prepare.R" | cut -d' ' -f1)" = "$PREP_MD5" ] && ok "threshold_prepare.R = revised P1" || ko "threshold_prepare.R unexpected"
[ "$(md5sum "$CODE/read_threshold.R" | cut -d' ' -f1)" = "$READ_MD5" ] && ok "read_threshold.R = P1b" || ko "read_threshold.R unexpected"
[ "$(md5sum "$DIR/step5_tile_worker_batch.R" | cut -d' ' -f1)" = "$WORKER_MD5" ] && ok "worker unchanged" || ko "worker changed"
[ "$(md5sum "$DIR/merge_rerun.sbatch" | cut -d' ' -f1)" = "$MRS_MD5" ] && ok "merge_rerun.sbatch unchanged" || ko "merge_rerun.sbatch changed"
[ -z "${BUFFER_IDX:-}" ] && ok "BUFFER_IDX not set" || ko "BUFFER_IDX set: $BUFFER_IDX"
[ -f "$BASE/logs/prepare.out" ] && ok "trace of the failed prepare present (will be set aside)" || ko "logs/prepare.out missing"
[ "$fail" -eq 0 ] || { echo; echo "CHECKS FAILED -- nothing was submitted."; exit 1; }

CMD="cd $CODE && R_LIBS_USER=$RSC Rscript threshold_prepare.R $BASE"
if [ "$DRY" -eq 1 ]; then
  echo "  [dry] mv logs/prepare.{out,err} -> logs/prepare_trial2_P1failed.{out,err} ; locks.csv -> locks_trial2_P1failed.csv"
  echo "  [dry] sbatch --parsable --partition=$PART --exclude=$EXCL --job-name=threshold_prep2 ... --wrap=\"$CMD\""
  echo; echo "Nothing was submitted (DRY=1)."; exit 0
fi
for x in out err; do mv -n "$BASE/logs/prepare.$x" "$BASE/logs/prepare_trial2_P1failed.$x"; done
mv -n "$BASE/locks.csv" "$BASE/locks_trial2_P1failed.csv"
JP=$(sbatch --parsable --partition="$PART" --exclude="$EXCL" --job-name=threshold_prep2 --qos=grit_normal --mem=64G \
       --cpus-per-task=2 --time=04:00:00 --output="$BASE/logs/prepare.out" --error="$BASE/logs/prepare.err" \
       --wrap="$CMD" | cut -d';' -f1)
echo "  prepare (revised P1): $JP"
echo "# Threshold, prepare rerun (revised P1) on $(date): $JP" >> "$CODE/ids_step1.txt"
