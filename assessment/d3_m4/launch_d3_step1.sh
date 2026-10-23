#!/bin/bash
# ============================================================================
# launch_d3_step1.sh -- D3 under M4, STEP 1 (rule written beforehand).
#
#   (a) step7 M4 variants, per label: (lambda, K) = (100,16) ref, (500,16), (100,1), (100,4), (100,64)
#   (b) lock S0: the variant (100,16) reproduces the installed M4 forcing (lock_m4.R, test P)
#   (c) d3_forcing.R: locks S1-S3, moved mass share TV, bounds, choix_ocim.txt
#   (d) d3_move.py: reconnaissance + construction of the moved forcing (Label A), locks R0, C1-C6
#
# NO OCIM run here. Step 2 (launch_d3_step2.sh) submits them after reading (c) and (d).
# The pipeline is not touched: lambda and K go through REDIST_LAMBDA_KM / REDIST_K_MAX, already provided
# by step7. Outputs under BASE, never in RUN.
#
# DRY=1 by default: checks and prints the sbatch commands without submitting anything. DRY=0 submits.
# Model: assessment/launch.sh and chain_step7_ocim_m4.sh.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
RUN=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/d3_m4_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/d3_m4"
GEBCO=${DFP_SCRATCH:?set DFP_SCRATCH}/configuration/gebco_2026/GEBCO_2026.nc
S7_MD5=8d49b9254847fe3d526ac9ef7932fb4c
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
RSC="$HOME/R/library"
VARIANTS="100:16 500:16 100:1 100:4 100:64"   # lambda_km:K_max; the first one is the reference

echo "============================================================"
echo " D3 under M4 -- step 1 (forcings only)    $(date)   DRY=$DRY"
echo "============================================================"
fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }

# ---- 1. checks ---------------------------------------------------------------
( cd "$DIR" && bash check_pipeline.sh >/dev/null 2>&1 ) && ok "check_pipeline.sh" || ko "check_pipeline.sh failed"
[ "$(md5sum "$DIR/step7_export_jtrawl.R" | cut -d' ' -f1)" = "$S7_MD5" ] \
  && ok "step7 installed = M4 ($S7_MD5)" || ko "installed step7 is not M4"
for f in "$DIR/step7_export_jtrawl_grit.sh" "$DIR/lock_m4.R" "$GEBCO" \
         "$CODE/d3_forcing.R" "$CODE/d3_move.py" "$CODE/read_d3.py"; do
  [ -f "$f" ] && ok "present: $(basename "$f")" || ko "missing: $f"
done
for lab in smooth flag; do
  n=$(ls "$RUN/cri_pr"/cri_final_${lab}_*.parquet 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: 1 C_ri in cri_pr" || ko "$lab: $n C_ri in cri_pr (1 expected)"
  n=$(ls "$RUN/step7_${lab}_m4"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: installed M4 forcing" || ko "$lab: $n .mat in step7_${lab}_m4 (1 expected)"
done
n=$(ls "$RUN/ocim_smooth_m4"/ocim_dredge_1x_*.mat 2>/dev/null | wc -l)
[ "$n" -eq 1 ] && ok "ocim_smooth_m4: one 1x run (grid + Jdredge_forcing)" || ko "ocim_smooth_m4: $n run 1x"
python3 -c "import numpy, scipy.io, h5py" 2>/dev/null && ok "python3: numpy, scipy, h5py" || ko "python3 without numpy/scipy/h5py"
if [ -e "$BASE" ] && [ -n "$(ls -A "$BASE" 2>/dev/null)" ]; then ko "BASE not empty: $BASE"; fi
FREE=$(df -BG --output=avail "$(dirname "$BASE")" | tail -1 | tr -dc 0-9)
[ "${FREE:-0}" -ge 10 ] && ok "free space: ${FREE} GB" || ko "free space ${FREE:-?} GB < 10 GB"
[ "$fail" -eq 0 ] || { echo; echo "CHECKS FAILED -- nothing was submitted."; exit 1; }

# ---- 2. submission -----------------------------------------------------------
sb() {
  local ex="${EXCL_JOB:-$EXCL}"
  if [ "$DRY" -eq 1 ]; then
    echo "  [dry] sbatch --parsable --partition=$PART --exclude=$ex $*" >&2
    echo "DRY"
  else
    sbatch --parsable --partition="$PART" --exclude="$ex" "$@" | cut -d';' -f1
  fi
}
[ "$DRY" -eq 1 ] || mkdir -p "$BASE/logs"
IDS=""
DEP_F=""
for lab in smooth flag; do
  for v in $VARIANTS; do
    l="${v%%:*}"; k="${v##*:}"; tag="l${l}_k${k}"; out="$BASE/step7_${lab}_${tag}"
    [ "$DRY" -eq 1 ] || mkdir -p "$out"
    # step7: hpc-08 excluded (R broken on that node), as in chain_step7_ocim_m4.sh
    J=$(EXCL_JOB="$EXCL,hpc-08.grit.ucsb.edu" sb --job-name="d3s7_${lab}_${tag}" \
          --output="$BASE/logs/s7_${lab}_${tag}_%j.out" --error="$BASE/logs/s7_${lab}_${tag}_%j.err" \
          --export=ALL,INJECTION_MODE=seafloor,REDIST_R_KM=500,REDIST_LAMBDA_KM=$l,REDIST_K_MAX=$k,OUT_DIR="$out",CRI_DIR="$RUN/cri_pr",CRI_TAG=$lab \
          "$DIR/step7_export_jtrawl_grit.sh")
    echo "  $lab $tag: step7 $J"
    IDS="$IDS $J"; DEP_F="$DEP_F:$J"
    if [ "$tag" = "l100_k16" ]; then
      # S0: lock_m4.R <installed M4 forcing> <reference variant> <label>
      JV=$(sb --dependency=afterok:"$J" --job-name="d3S0_$lab" --qos=grit_normal \
             --mem=16G --cpus-per-task=1 --time=00:20:00 \
             --output="$BASE/logs/S0_${lab}.out" --error="$BASE/logs/S0_${lab}.err" \
             --wrap="cd $DIR && R_LIBS_USER=$RSC Rscript lock_m4.R $RUN/step7_${lab}_m4 $out $lab")
      echo "  $lab: lock S0 $JV"
      IDS="$IDS $JV"; DEP_F="$DEP_F:$JV"
    fi
  done
done
# (c) measurement of the forcings: afterany, so that failed variants are REPORTED (S1), not silent
JF=$(sb --dependency=afterany"$DEP_F" --job-name="d3_forcing" --qos=grit_normal \
       --mem=16G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/etape_S.out" --error="$BASE/logs/etape_S.err" \
       --wrap="cd $CODE && R_LIBS_USER=$RSC Rscript d3_forcing.R $BASE $RUN")
echo "  measurement of the forcings (step S): $JF"
# (d) displacement: independent of the variants
JC=$(sb --job-name="d3_move" --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=01:00:00 \
       --output="$BASE/logs/move.out" --error="$BASE/logs/move.err" \
       --wrap="cd $CODE && python3 d3_move.py $RUN $BASE")
echo "  reconnaissance + construction of the displacement: $JC"
IDS="$IDS $JF $JC"

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_step1.txt" <<EOT
# D3 under M4, step 1, submitted on $(date)   BASE=$BASE
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_step1.txt"
echo "  Read next: $BASE/logs/etape_S.out, $BASE/logs/move.out, then launch_d3_step2.sh"
