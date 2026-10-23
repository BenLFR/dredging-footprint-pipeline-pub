#!/bin/bash
# ============================================================================
# Chain step4 -> OCIM -> postproc, 2026-09-22
#
# Reruns everything after the lithology step with the Sala et al. (2021) lability
# (PL_SCHEME=sala2021), unattended: each link starts afterok of the
# previous one; two locks stop the chain if the rerun changes anything other
# than the lability.
#
#   step4 -> LOCK4 -> step5 smooth -> merge -> LOCK_M -> step6 -> step7 -> OCIM -> postproc
#                  \-> step5 flag   -> merge -> LOCK_M -> step6 -> step7 -> OCIM -> postproc
#                                                          (both step6) -> C_ri SUMMARY
#
# Everything is written to $RUN (hpc-scratch). Nothing is written to output_V6
# (published) or to output_V6_rejeu (the August run, reference of the locks).
#
# Two inputs are REUSED from the August run; the merge lock is what
# allows it:
#   - the step4 input (step3, AIS_data_core_preprocessed_V6_20260810_093400.rds);
#   - the fraction left in place (decharge/fraction_sur_place_<lab>.csv),
#     computed from SAR, SVR and sum_d, which do not depend on p_l.
#
# Usage:   bash chain_step4_ocim.sh --dry-run   (checks, submits nothing)
#           bash chain_step4_ocim.sh             (submits)
# ============================================================================
set -euo pipefail

DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
REF="${DFP_DATA:?set DFP_DATA}/output_V6_rejeu"                      # August run: read-only
# the new run. output_V6_sala2021 (without _r2) is the launch abandoned on
# 2026-09-22 (step4 done, step5 killed by hpc-10/11/12): do not use it.
RUN="${RUN:-${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r2}"
REF_STEP4="$REF/AIS_with_lithology_20260811_152541.rds"
REF_STEP3="$REF/AIS_data_core_preprocessed_V6_20260810_093400.rds"
# The 71 tiles with dredging pings (identical for smooth and flag,
# listed from the August step5_smooth and step5_flag). An array 1-525 would exceed
# MaxSubmitPU=300 of the grit_normal QOS from the second label on.
TILES=104,105,118,125,153,154,169,185,186,189,197,205,206,220,232,238,239,241,255,262,263,270,272,273,276,289,290,291,296,297,298,305,306,307,308,323,327,331,337,338,339,340,341,344,345,357,366,371,372,373,374,379,380,401,402,403,404,405,406,415,437,438,439,440,441,451,452,473,475,511,515
THROTTLE=20
# RUN is on beegfs (/home/hpc-scratch), which sequoia -- the only node of
# emlab_nodes -- does not mount (probe of 2026-09-22: sequoia ABSENT, hpc-08 and
# hpc-11 MOUNTED_WRITABLE; first launch died on "mkdir: cannot create
# directory /home/hpc-scratch", job 1117275). All submissions therefore force
# grit_nodes; the command line overrides the #SBATCH headers of the launchers.
PART=grit_nodes
# Within grit_nodes, hpc-10/11/12 run R 4.6.x: ~/R/library (compiled
# under 4.5) does not load there ("Rcpp.so: undefined symbol: VECTOR_PTR", 71/71
# step5 flag tasks dead, job 1117310). Probe of 2026-09-22: hpc-07/08/09/14
# to 18 run R 4.5.3 and load sf. hpc-13 was already excluded by the launchers. The
# command line replaces their --exclude, so step7 adds hpc-08 back (EXCL_JOB).
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu
RSC="$HOME/R/library"
PP_PY="$DIR/postproc_atwood_style.py"
PP_REQ="$DIR/requirements_postproc_atwood.txt"
PP_SB="${DFP_REPO:?set DFP_REPO}/pipeline/postproc_atwood_style.sbatch"
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

echo "============================================================"
echo " CHAIN step4 -> OCIM (PL_SCHEME=sala2021)    $(date)"
echo " RUN = $RUN"
echo "============================================================"

# ---- 1. preliminary checks ------------------------------------------------
echo "--- checks ---"
fail=0
ok()  { echo "  ✅ $*"; }
ko()  { echo "  🔴 $*"; fail=1; }
for f in "$DIR/step4_rerun.sbatch" "$DIR/step4_add_lithology_vNext.R" \
         "$DIR/step5_batch.sbatch" "$DIR/step5_tile_worker_batch.R" \
         "$DIR/merge_rerun.sbatch" "$DIR/step5_merge_tiles_corr.R" \
         "$DIR/step6_pr.sbatch" "$DIR/step6_calculate_cri_corr.R" \
         "$DIR/step7_export_jtrawl_grit.sh" "$DIR/step7_export_jtrawl.R" \
         "$DIR/co2model_batch_grit.sh" "$PP_SB" "$PP_PY" "$PP_REQ" \
         "$DIR/lock_sala_step4.R" "$DIR/lock_sala_merge.R" "$DIR/summary_sala.R" \
         "$REF_STEP3" "$REF_STEP4" "$REF/cri_pr" \
         "${DFP_DATA:?set DFP_DATA}/output_V6/tiles_1000km.gpkg" \
         "${DFP_DATA:?set DFP_DATA}/configuration/ocim/ocim_cache.mat" \
         "${DFP_DATA:?set DFP_DATA}/hubocean_cache"; do
  if [ -e "$f" ]; then ok "$f"; else ko "MISSING: $f"; fi
done
for lab in smooth flag; do
  for f in "$REF/decharge/fraction_sur_place_${lab}.csv" "$REF/merge_${lab}_pr"; do
    if [ -e "$f" ]; then ok "$f"; else ko "MISSING: $f"; fi
  done
done
for f in step4_rerun.sbatch step5_batch.sbatch merge_rerun.sbatch step6_pr.sbatch \
         step7_export_jtrawl_grit.sh co2model_batch_grit.sh; do
  if bash -n "$DIR/$f" 2>/dev/null; then ok "syntax: $f"; else ko "SYNTAX: $f"; fi
done
# does the installed step4 carry the scheme switch?
if grep -q 'PL_SCHEME <- Sys.getenv("PL_SCHEME", "sala2021")' "$DIR/step4_add_lithology_vNext.R"; then
  ok "step4 carries the PL_SCHEME switch (default sala2021)"
else ko "installed step4 does NOT carry the PL_SCHEME switch"; fi
# can the installed step7 resolve CRI_DIR?
if grep -q 'CRI_DIR requires CRI_TAG' "$DIR/step7_export_jtrawl_grit.sh"; then
  ok "step7 resolves CRI_DIR/CRI_TAG"
else ko "installed step7 does not resolve CRI_DIR"; fi
# is the pipeline manifest up to date?
if bash "$DIR/check_pipeline.sh" >/dev/null 2>&1; then ok "check_pipeline.sh passes"
else ko "check_pipeline.sh fails -- run: bash $DIR/check_pipeline.sh"; fi
# the new run must start from an empty directory: the locks and step6/7
# refuse (or pick the wrong file) as soon as two files match the same pattern.
# (empty subdirectories do no harm: only files are checked)
if [ -d "$RUN" ] && [ -n "$(find "$RUN" -type f -print -quit 2>/dev/null)" ]; then
  ko "RUN already contains files: $RUN -- choose another RUN=... or empty it yourself"
else ok "RUN has no file (or is absent)"; fi
PUB=$(readlink -f "${DFP_DATA:?set DFP_DATA}/output_V6"); AUG=$(readlink -f "$REF")
case "$(readlink -f "$RUN")" in
  "$PUB"|"$PUB"/*|"$AUG"|"$AUG"/*)
    ko "RUN points into output_V6 or into the August run -- REFUSED" ;;
esac
# space: ~3.7 GB per OCIM directory, ~1 GB for the rest
FREE=$(df -BG --output=avail "$(dirname "$RUN")" | tail -1 | tr -dc 0-9)
if [ "${FREE:-0}" -ge 20 ]; then ok "free space under $(dirname "$RUN"): ${FREE} GB"
else ko "free space ${FREE:-?} GB < 20 GB"; fi
FREE_H=$(df -BM --output=avail "$HOME" | tail -1 | tr -dc 0-9)
echo "  ℹ️  free space under \$HOME (logs, tmp_terra): ${FREE_H} MB"
N_T=$(echo "$TILES" | tr ',' '\n' | wc -l)
if [ "$N_T" -eq 71 ]; then ok "71 tiles in the list"; else ko "$N_T tiles (71 expected)"; fi

[ "$fail" -eq 0 ] || { echo; echo "🔴 CHECKS FAILED -- nothing was submitted."; exit 1; }
echo "  ✅ all checks pass"

# ---- 2. submission -----------------------------------------------------------
# sb <name> <sbatch args...>: submits, or with --dry-run prints and returns a fake id
sb() {
  local item_name="$1"; shift
  local ex="${EXCL_JOB:-$EXCL}"
  if [ "$DRY" -eq 1 ]; then
    echo "  [dry] sbatch --parsable --partition=$PART --exclude=$ex $*" >&2
    echo "DRY_${item_name}"
  else
    sbatch --parsable --partition="$PART" --exclude="$ex" "$@" | cut -d';' -f1
  fi
}
# short R job afterok, for the locks and the summary
sb_r() {
  local item_name="$1" dep="$2"; shift 2
  sb "$item_name" --dependency=afterok:"$dep" --job-name="$item_name" \
     --qos=grit_normal \
     --mem=48G --cpus-per-task=2 --time=01:00:00 \
     --output="$DIR/logs/${item_name}_%j.out" --error="$DIR/logs/${item_name}_%j.err" \
     --wrap="cd $DIR && R_LIBS_USER=$RSC Rscript $*"
}

echo
if [ "$DRY" -eq 1 ]; then echo "--- DRY RUN: commands that WOULD be submitted ---"; else echo "--- submission ---"; fi
# step7 and co2model only write into an existing directory: they are created here.
[ "$DRY" -eq 1 ] || mkdir -p "$RUN" "$DIR/logs" "$RUN"/step7_{smooth,flag} "$RUN"/ocim_{smooth,flag}

J4=$(sb step4 --job-name=s4_sala \
       --export=ALL,PIPE_IN_DIR="$REF",PIPE_OUT_DIR="$RUN",PIPE_AIS_FILE="$REF_STEP3",PL_SCHEME=sala2021 \
       "$DIR/step4_rerun.sbatch")
echo "  step4              -> $J4"
JV4=$(sb_r v4_sala "$J4" lock_sala_step4.R "$REF_STEP4" "$RUN")
echo "  lock step4         -> $JV4"

declare -A J6
IDS="$J4 $JV4"
for lab in smooth flag; do
  J5=$(sb "s5_$lab" --dependency=afterok:"$JV4" --array="$TILES%$THROTTLE" --job-name="s5_sala_$lab" \
         --export=ALL,DREDGE_LABEL=$lab,PIPE_RUN_DIR="$RUN" "$DIR/step5_batch.sbatch")
  JM=$(sb "m_$lab" --dependency=afterok:"$J5" --job-name="m_sala_$lab" \
         --export=ALL,DREDGE_LABEL=$lab,MERGE_IN_DIR="$RUN/step5_$lab",MERGE_OUT_DIR="$RUN/merge_${lab}_pr",FI_FRAC_FILE="$REF/decharge/fraction_sur_place_${lab}.csv",MERGE_N_EXPECTED=71 \
         "$DIR/merge_rerun.sbatch")
  JVM=$(sb_r "vm_$lab" "$JM" lock_sala_merge.R "$REF/merge_${lab}_pr" "$RUN/merge_${lab}_pr")
  J6[$lab]=$(sb "s6_$lab" --dependency=afterok:"$JVM" --job-name="s6_sala_$lab" \
         --export=ALL,DREDGE_LABEL=$lab,PIPE_RUN_DIR="$RUN" "$DIR/step6_pr.sbatch")
  J7=$(EXCL_JOB="$EXCL,hpc-08.grit.ucsb.edu" sb "s7_$lab" --dependency=afterok:"${J6[$lab]}" --job-name="s7_sala_$lab" \
         --export=ALL,OUT_DIR="$RUN/step7_$lab",CRI_DIR="$RUN/cri_pr",CRI_TAG=$lab \
         "$DIR/step7_export_jtrawl_grit.sh")
  # MATLAB exists only on hpc-08 (command of job 1086015, reused as is)
  JO=$(sb "co2_$lab" --dependency=afterok:"$J7" --job-name="co2_sala_$lab" \
         --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_DIR="$RUN/step7_$lab",OCIM_OUT_DIR="$RUN/ocim_$lab" \
         "$DIR/co2model_batch_grit.sh")
  JP=$(sb "pp_$lab" --dependency=afterok:"$JO" --job-name="pp_sala_$lab" \
         --export=ALL,OUTPUT_ROOT="$RUN/ocim_$lab",POSTPROC_DIR="$RUN/postproc_$lab",SCRIPT_PATH="$PP_PY",REQ_PATH="$PP_REQ" \
         "$PP_SB")
  echo "  $lab: step5 $J5 -> merge $JM -> lock $JVM -> step6 ${J6[$lab]} -> step7 $J7 -> OCIM $JO -> postproc $JP"
  IDS="$IDS $J5 $JM $JVM ${J6[$lab]} $J7 $JO $JP"
done
JB=$(sb_r summary_sala "${J6[smooth]}:${J6[flag]}" summary_sala.R "$REF/cri_pr" "$RUN/cri_pr")
echo "  C_ri summary       -> $JB"
IDS="$IDS $JB"

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (--dry-run)."; exit 0; fi

cat > "$HOME/chaine_sala2021_ids.txt" <<EOT
# chain step4 -> OCIM (sala2021) submitted on $(date)
# RUN=$RUN
# ids: $IDS
# cancel everything:  scancel $IDS
EOT
echo
echo "============================================================"
echo " ids written to ~/chaine_sala2021_ids.txt"
echo " monitor    : squeue -u \$USER"
echo " CANCEL     : scancel $IDS"
echo "============================================================"
