#!/bin/bash
# ============================================================================
# launch_d3_step2.sh -- D3 under M4, STEP 2 (rule written beforehand).
# Submits the OCIM runs decided by the rule, then the reading:
#   - each line of BASE/choix_ocim.txt ("<label> <variant>"), written by d3_forcing.R;
#   - the displacement (Label A), if BASE/MOVE_PASS exists (written by d3_move.py).
# OCIM on hpc-08 (--nodelist), as in chain_step7_ocim_m4.sh; about 16 h per run.
# To be launched ONLY after reading logs/etape_S.out and logs/move.out (locks S0-S3, R0, C1-C6).
# DRY=1 by default.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/d3_m4_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/d3_m4"
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu

[ -f "$BASE/choix_ocim.txt" ] || { echo "ABORT: $BASE/choix_ocim.txt missing (step S not read, or failed)"; exit 1; }
[ -f "$DIR/co2model_batch_grit.sh" ] || { echo "ABORT: co2model_batch_grit.sh missing"; exit 1; }
sb() {
  if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch --parsable --partition=$PART --exclude=$EXCL $*" >&2; echo "DRY"
  else sbatch --parsable --partition="$PART" --exclude="$EXCL" "$@" | cut -d';' -f1; fi
}
declare -a CASES=()
while read -r lab tag; do [ -n "$lab" ] && CASES+=("$lab|step7_${lab}_${tag}|ocim_${lab}_${tag}"); done < "$BASE/choix_ocim.txt"
if [ -f "$BASE/MOVE_PASS" ]; then CASES+=("smooth|step7_smooth_m4_shift|ocim_smooth_shift")
else echo "  displacement: no MOVE_PASS -> no OCIM (see logs/move.out)"; fi

IDS=""; DEP=""
for c in "${CASES[@]}"; do
  IFS='|' read -r lab fdir odir <<< "$c"
  n=$(ls "$BASE/$fdir"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] || { echo "ABORT: $n .mat in $BASE/$fdir (1 expected)"; exit 1; }
  [ -e "$BASE/$odir" ] && { echo "ABORT: $BASE/$odir already exists"; exit 1; }
  [ "$DRY" -eq 1 ] || mkdir -p "$BASE/$odir"
  JO=$(sb --job-name="d3co2_${odir#ocim_}" --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_DIR="$BASE/$fdir",OCIM_OUT_DIR="$BASE/$odir" \
         "$DIR/co2model_batch_grit.sh")
  echo "  $odir: OCIM $JO  (forcing $fdir)"
  IDS="$IDS $JO"; DEP="$DEP:$JO"
done
JL=$(sb --dependency=afterany"$DEP" --job-name="d3_read" --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/read_d3.out" --error="$BASE/logs/read_d3.err" \
       --wrap="cd $CODE && python3 read_d3.py")
echo "  reading: $JL"
IDS="$IDS $JL"
if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_step2.txt" <<EOT
# D3 under M4, step 2 (OCIM), submitted on $(date)
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_step2.txt"
