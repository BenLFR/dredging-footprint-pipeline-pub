#!/bin/bash
# ============================================================================
# chain_step7_ocim_m4.sh -- replay step7 -> OCIM with injection at seafloor depth
#
#   step7 (seafloor) -> lock_m4.R (test P) -> OCIM -> link -> postproc   (smooth and flag)
#
# What is REUSED from run _r6 (output_V6_sala2021_r6_flux10), and why this is allowed:
#   - the step6 C_ri (cri_pr): M4 only touches step7, and test L validated the _r6 base.
# Outputs in *_m4 folders of the same RUN: the _r6 outputs (column bottom) remain intact
# and serve as comparison.
# Test P: the forcing must be identical to the test O candidate (profondeur_injection/candidat).
#
# Usage: bash chain_step7_ocim_m4.sh [--dry-run]
# ============================================================================
set -euo pipefail

DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
RUN="${RUN:-${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10}"
CAND="$HOME/profondeur_injection/candidat"
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu
RSC="$HOME/R/library"
PP_PY="$DIR/postproc_atwood_style.py"
PP_REQ="$DIR/requirements_postproc_atwood.txt"
PP_SB="${DFP_REPO:?set DFP_REPO}/pipeline/postproc_atwood_style.sbatch"
S7_MD5=8d49b9254847fe3d526ac9ef7932fb4c
IDS_FILE="${IDS_FILE:-$HOME/chaine_m4_ids.txt}"
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

echo "============================================================"
echo " CHAIN step7 (M4 seafloor) -> OCIM    $(date)"
echo "   RUN = $RUN"
echo "============================================================"
fail=0
ok() { echo "  ✅ $*"; }
ko() { echo "  🔴 $*"; fail=1; }

# ---- 1. checks ---------------------------------------------------------------
[ "$(md5sum "$DIR/step7_export_jtrawl.R" | cut -d' ' -f1)" = "$S7_MD5" ] \
  && ok "step7 installed = M4 ($S7_MD5)" || ko "installed step7 is not M4"
for f in "$DIR/step7_export_jtrawl_grit.sh" "$DIR/co2model_batch_grit.sh" "$DIR/lock_m4.R" \
         "$PP_PY" "$PP_REQ" "$PP_SB" ${DFP_SCRATCH:?set DFP_SCRATCH}/configuration/gebco_2026/GEBCO_2026.nc; do
  [ -f "$f" ] && ok "present: $(basename "$f")" || ko "missing: $f"
done
for lab in smooth flag; do
  n=$(ls "$RUN/cri_pr"/cri_final_${lab}_*.parquet 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: 1 C_ri in cri_pr" || ko "$lab: $n C_ri in cri_pr (1 expected)"
  n=$(ls "$CAND/${lab}_seafloor"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: candidate forcing (test O)" || ko "$lab: $n candidate .mat (1 expected)"
  for d in step7_${lab}_m4 ocim_${lab}_m4 postproc_${lab}_m4; do
    [ -e "$RUN/$d" ] && ko "already exists: $RUN/$d" || true
  done
done
FREE=$(df -BG --output=avail "$RUN" | tail -1 | tr -dc 0-9)
[ "${FREE:-0}" -ge 20 ] && ok "free space: ${FREE} GB" || ko "free space ${FREE:-?} GB < 20 GB"
[ "$fail" -eq 0 ] || { echo; echo "🔴 CHECKS FAILED -- nothing was submitted."; exit 1; }

# ---- 2. submission -----------------------------------------------------------
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

[ "$DRY" -eq 1 ] || mkdir -p "$DIR/logs" "$RUN"/step7_{smooth,flag}_m4 "$RUN"/ocim_{smooth,flag}_m4
IDS=""
for lab in smooth flag; do
  J7=$(EXCL_JOB="$EXCL,hpc-08.grit.ucsb.edu" sb "s7_$lab" --job-name="s7_m4_$lab" \
         --export=ALL,INJECTION_MODE=seafloor,OUT_DIR="$RUN/step7_${lab}_m4",CRI_DIR="$RUN/cri_pr",CRI_TAG=$lab \
         "$DIR/step7_export_jtrawl_grit.sh")
  JV=$(sb "vP_$lab" --dependency=afterok:"$J7" --job-name="vP_m4_$lab" --qos=grit_normal \
         --mem=16G --cpus-per-task=1 --time=00:20:00 \
         --output="$DIR/logs/vP_m4_${lab}_%j.out" --error="$DIR/logs/vP_m4_${lab}_%j.err" \
         --wrap="cd $DIR && R_LIBS_USER=$RSC Rscript lock_m4.R $CAND/${lab}_seafloor $RUN/step7_${lab}_m4 $lab")
  JO=$(sb "co2_$lab" --dependency=afterok:"$JV" --job-name="co2_m4_$lab" \
         --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_DIR="$RUN/step7_${lab}_m4",OCIM_OUT_DIR="$RUN/ocim_${lab}_m4" \
         "$DIR/co2model_batch_grit.sh")
  # postproc_atwood_style.py only looks for the forcing in OUTPUT_ROOT: we link it.
  JL=$(sb "lk_$lab" --dependency=afterok:"$JO" --job-name="lk_m4_$lab" --qos=grit_normal \
         --mem=1G --time=00:05:00 --output="$DIR/logs/lk_m4_${lab}_%j.out" \
         --wrap="set -e; f=\$(ls $RUN/step7_${lab}_m4/jdredge_ocim2_48l_*.mat); [ \$(echo \"\$f\" | wc -l) -eq 1 ]; ln -s \"\$f\" $RUN/ocim_${lab}_m4/")
  JP=$(sb "pp_$lab" --dependency=afterok:"$JL" --job-name="pp_m4_$lab" \
         --export=ALL,OUTPUT_ROOT="$RUN/ocim_${lab}_m4",POSTPROC_DIR="$RUN/postproc_${lab}_m4",SCRIPT_PATH="$PP_PY",REQ_PATH="$PP_REQ" \
         "$PP_SB")
  echo "  $lab: step7 $J7 -> test P $JV -> OCIM $JO -> link $JL -> postproc $JP"
  IDS="$IDS $J7 $JV $JO $JL $JP"
done

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (--dry-run)."; exit 0; fi
cat > "$IDS_FILE" <<EOT
# step7 chain (M4 seafloor) -> OCIM submitted on $(date)
# RUN=$RUN  (outputs *_m4)
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $IDS_FILE"
