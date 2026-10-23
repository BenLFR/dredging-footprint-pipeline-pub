#!/bin/bash
# ============================================================================
# launch_threshold_step1.sh -- Threshold +/- 0.005 (Label B), STEP 1 (rule written beforehand).
#
#   (a) threshold_prepare.R: counts, P0 (reproduction of the flag), flags t+/-, P1 (tiles), rds of the arms
#   (b) P2: merge + step6 of _r6 rerun on the step5 tiles of _r2 (flag, with the original frac table),
#       then read_threshold.R p2 -> BASE/P2_PASS
# (a) and (b) are independent. Step 2 (step5 of the three arms) requires PREPARATION_PASS and P2_PASS.
#
# Read on GRIT: the worker identifies the tiles by
# the tile_id column, takes Lon/Lat, enlarges the tile by TILE_BUFFER_M (50 km) and accepts a declared input
# (PIPE_LITHO_FILE); merge_rerun.sbatch leaves FI_FRAC_FILE empty by default. The md5 are checked here.
# DRY=1 by default.
# ============================================================================
set -euo pipefail
unset FI_FRAC_FILE
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
REF="${DFP_DATA:?set DFP_DATA}/output_V6_rejeu"
R2=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r2
RUN=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/seuil_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/threshold"
MERGE_MD5=e40847579fa81eacadd1955e1013b42b      # merge that produced _r6
PART=grit_nodes
# hpc-08 excluded: R jobs fail there when loading data.table and Rcpp
# ("undefined symbol: SETLENGTH / VECTOR_PTR"); with the same library, R jobs passed on
# hpc-07, 09, 14, 15, 16 and 18. hpc-08 is reserved for MATLAB (OCIM).
EXCL=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
RSC="$HOME/R/library"
WORKER_MD5=cd6e9c4d6932bb249df8af1c6a919431     # step5_tile_worker_batch.R, as read
MRS_MD5=4dbc5b1d4229653698a83dcb7d5998dd        # merge_rerun.sbatch, as read

echo "============================================================"
echo " Threshold +/- 0.005 -- step 1    $(date)   DRY=$DRY"
echo "============================================================"
fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }
( cd "$DIR" && bash check_pipeline.sh >/dev/null 2>&1 ) && ok "check_pipeline.sh" || ko "check_pipeline.sh failed"
[ "$(md5sum "$DIR/step5_merge_tiles_corr.R" | cut -d' ' -f1)" = "$MERGE_MD5" ] \
  && ok "merge installed = the _r6 one ($MERGE_MD5)" || ko "merge installed differs from the _r6 one"
for f in "$DIR/step3_merge_final.R" "$DIR/step5_batch.sbatch" "$DIR/step5_tile_worker_batch.R" \
         "$DIR/merge_rerun.sbatch" "$DIR/step6_pr.sbatch" "$DIR/step6_calculate_cri_corr.R" \
         "$R2/AIS_with_lithology_20260922_221357.rds" \
         "$REF/dragage_gridsearch_results_V6_20260810_093400.rds" \
         "$REF/decharge/fraction_sur_place_flag.csv" "${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg" \
         "$CODE/threshold_prepare.R" "$CODE/read_threshold.R" "$CODE/launch_threshold_step2.sh"; do
  [ -f "$f" ] && ok "present: $(basename "$f")" || ko "missing: $f"
done
n=$(ls "$R2/step5_flag"/sar_*.parquet 2>/dev/null | wc -l)
[ "$n" -eq 71 ] && ok "71 flag sar tiles in _r2" || ko "$n flag sar tiles in _r2 (71 expected)"
n=$(ls "$RUN/cri_pr"/cri_final_flag_*.parquet 2>/dev/null | wc -l)
[ "$n" -eq 1 ] && ok "reference C_ri _r6 flag" || ko "$n cri_final_flag in _r6 (1 expected)"
[ "$(md5sum "$DIR/step5_tile_worker_batch.R" | cut -d' ' -f1)" = "$WORKER_MD5" ] \
  && ok "worker = the version read" || ko "worker changed since it was read"
[ "$(md5sum "$DIR/merge_rerun.sbatch" | cut -d' ' -f1)" = "$MRS_MD5" ] \
  && ok "merge_rerun.sbatch = the version read" || ko "merge_rerun.sbatch changed since it was read"
[ -f "$DIR/constants.R" ] && ok "present: constants.R" || ko "missing: constants.R"
[ -z "${BUFFER_IDX:-}" ] && ok "BUFFER_IDX not set (50 km buffer)" || ko "BUFFER_IDX set in the environment: $BUFFER_IDX"
if [ -e "$BASE" ] && [ -n "$(ls -A "$BASE" 2>/dev/null)" ]; then ko "BASE not empty: $BASE"; fi
FREE=$(df -BG --output=avail "$(dirname "$BASE")" | tail -1 | tr -dc 0-9)
[ "${FREE:-0}" -ge 20 ] && ok "free space: ${FREE} GB" || ko "free space ${FREE:-?} GB < 20 GB"
[ "$fail" -eq 0 ] || { echo; echo "CHECKS FAILED -- nothing was submitted."; exit 1; }

sb() {
  if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch --parsable --partition=$PART --exclude=$EXCL $*" >&2; echo "DRY"
  else sbatch --parsable --partition="$PART" --exclude="$EXCL" "$@" | cut -d';' -f1; fi
}
[ "$DRY" -eq 1 ] || mkdir -p "$BASE/logs" "$BASE/p2"
JP=$(sb --job-name="threshold_prep" --qos=grit_normal --mem=64G --cpus-per-task=2 --time=04:00:00 \
       --output="$BASE/logs/prepare.out" --error="$BASE/logs/prepare.err" \
       --wrap="cd $CODE && R_LIBS_USER=$RSC Rscript threshold_prepare.R $BASE")
echo "  preparation (P0, P1, arms): $JP"
JM=$(sb --job-name="threshold_m_p2" \
       --output="$BASE/logs/m_p2_%j.out" --error="$BASE/logs/m_p2_%j.err" \
       --export=ALL,DREDGE_LABEL=flag,MERGE_IN_DIR="$R2/step5_flag",MERGE_OUT_DIR="$BASE/p2/merge_flag_pr",FI_FRAC_FILE="$REF/decharge/fraction_sur_place_flag.csv",MERGE_N_EXPECTED=71 \
       "$DIR/merge_rerun.sbatch")
J6=$(sb --dependency=afterok:"$JM" --job-name="threshold_s6_p2" \
       --output="$BASE/logs/s6_p2_%j.out" --error="$BASE/logs/s6_p2_%j.err" \
       --export=ALL,DREDGE_LABEL=flag,PIPE_RUN_DIR="$BASE/p2" "$DIR/step6_pr.sbatch")
JL=$(sb --dependency=afterok:"$J6" --job-name="threshold_P2" --qos=grit_normal --mem=16G --cpus-per-task=1 --time=00:20:00 \
       --output="$BASE/logs/P2.out" --error="$BASE/logs/P2.err" \
       --wrap="cd $CODE && R_LIBS_USER=$RSC Rscript read_threshold.R $BASE $RUN p2")
echo "  P2: merge $JM -> step6 $J6 -> reading $JL"
IDS="$JP $JM $J6 $JL"
if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_step1.txt" <<EOT
# Threshold, step 1, submitted on $(date)   BASE=$BASE
# ids: $IDS
# cancel everything:  scancel $IDS
EOT
echo "  ids written to $CODE/ids_step1.txt"
echo "  Read next: $BASE/logs/prepare.out, $BASE/logs/P2.out, then launch_threshold_step2.sh"
