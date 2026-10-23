#!/bin/bash
# ============================================================================
# launch_svr_ocim.sh -- OCIM run dedicated to the annual SVR cap (rule written beforehand).
# Model: launch_d3_step1.sh and launch_d3_step2.sh, merged; OCIM starts only if the locks of the
# forcing pass (afterok on svr_forcing.R, which exits with an error otherwise).
#   (A)  make_cri_year.R          : input C_ri x s_i, locks A1-A3
#   (B)  step7 M4 reference      : cri_pr, lambda 100, K 16 -> lock S0 (lock_m4.R, test P)
#   (C)  step7 M4 annual         : CRI_FILE = input (A)
#   (D)  svr_forcing.R           : S0, S3, S1, S2 ; writes Q_<lab>.txt, TV_<lab>.txt, FORCING_PASS
#   (E)  OCIM x 2 on hpc-08      : afterok (D) ; about 16 h per run
#   (F)  read_svr.py             : afterany (E) ; L0-L4 then dAF
# Upstream dependencies are afterany: a failing link makes the next one fail, instead of leaving it
# pending forever; (D) decides. If (D) fails: scancel the remaining ids (ids_svr_ocim.txt).
# Pipeline untouched; outputs under BASE, never in RUN. DRY=1 by default.
# Full node names in --exclude ("hpc-08" alone is refused).
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
RUN=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/svr_annuel_ocim_20261005
CODE="${DFP_REPO:?set DFP_REPO}/annual_svr/ocim"
S7_MD5=8d49b9254847fe3d526ac9ef7932fb4c
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
EXCL_R="$EXCL,hpc-08.grit.ucsb.edu"      # R is broken on hpc-08
ENVR="R_LIBS_USER=\$HOME/R/library LD_LIBRARY_PATH=\$HOME/lib:\${LD_LIBRARY_PATH:-}"
S7_ENV="INJECTION_MODE=seafloor,REDIST_R_KM=500,REDIST_LAMBDA_KM=100,REDIST_K_MAX=16"

echo "============================================================"
echo " Annual SVR cap -- dedicated OCIM run    $(date)   DRY=$DRY"
echo "============================================================"
fail=0
ok() { echo "  OK  $*"; }
ko() { echo "  KO  $*"; fail=1; }

# ---- 1. checks ---------------------------------------------------------------
( cd "$DIR" && bash check_pipeline.sh >/dev/null 2>&1 ) && ok "check_pipeline.sh" || ko "check_pipeline.sh failed"
[ "$(md5sum "$DIR/step7_export_jtrawl.R" | cut -d' ' -f1)" = "$S7_MD5" ] \
  && ok "step7 installed = M4 ($S7_MD5)" || ko "installed step7 is not M4"
for f in "$DIR/step7_export_jtrawl_grit.sh" "$DIR/lock_m4.R" "$DIR/co2model_batch_grit.sh" \
         "$CODE/make_cri_year.R" "$CODE/svr_forcing.R" "$CODE/read_svr.py" \
         "$CODE/s_cibles_smooth.csv" "$CODE/s_cibles_flag.csv"; do
  [ -f "$f" ] && ok "present: $(basename "$f")" || ko "missing: $f"
done
( cd "$CODE" && md5sum make_cri_year.R svr_forcing.R read_svr.py s_cibles_smooth.csv s_cibles_flag.csv | sed 's/^/      /' )
for lab in smooth flag; do
  n=$(ls "$RUN/cri_pr"/cri_final_${lab}_*.parquet 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: 1 C_ri in cri_pr" || ko "$lab: $n C_ri in cri_pr (1 expected)"
  n=$(ls "$RUN/step7_${lab}_m4"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: installed M4 forcing" || ko "$lab: $n .mat in step7_${lab}_m4 (1 expected)"
  n=$(ls "$RUN/ocim_${lab}_m4"/ocim_dredge_1x_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] && ok "$lab: reference OCIM M4 run" || ko "$lab: $n 1x run in ocim_${lab}_m4"
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
IDS=""
JA=$(EXCL_JOB="$EXCL_R" sb --job-name=svrA --qos=grit_normal --mem=16G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/A.out" --error="$BASE/logs/A.err" \
       --wrap="export $ENVR; cd $CODE && Rscript make_cri_year.R $RUN $BASE $CODE")
echo "  (A) input C_ri x s_i: $JA"; IDS="$IDS $JA"
DEP_D=":$JA"
for lab in smooth flag; do
  for v in ref an; do [ "$DRY" -eq 1 ] || mkdir -p "$BASE/step7_${lab}_$v"; done
  JR=$(EXCL_JOB="$EXCL_R" sb --job-name="svrS7ref_$lab" \
         --output="$BASE/logs/s7_${lab}_ref_%j.out" --error="$BASE/logs/s7_${lab}_ref_%j.err" \
         --export=ALL,$S7_ENV,OUT_DIR="$BASE/step7_${lab}_ref",CRI_DIR="$RUN/cri_pr",CRI_TAG=$lab \
         "$DIR/step7_export_jtrawl_grit.sh")
  JV=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany:"$JR" --job-name="svrS0_$lab" --qos=grit_normal \
         --mem=16G --cpus-per-task=1 --time=00:20:00 \
         --output="$BASE/logs/S0_${lab}.out" --error="$BASE/logs/S0_${lab}.err" \
         --wrap="export $ENVR; cd $DIR && Rscript lock_m4.R $RUN/step7_${lab}_m4 $BASE/step7_${lab}_ref $lab")
  JN=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany:"$JA" --job-name="svrS7annual_$lab" \
         --output="$BASE/logs/s7_${lab}_annual_%j.out" --error="$BASE/logs/s7_${lab}_annual_%j.err" \
         --export=ALL,$S7_ENV,OUT_DIR="$BASE/step7_${lab}_an",CRI_FILE="$BASE/cri_an_${lab}/cri_final_${lab}_an_20261005.parquet" \
         "$DIR/step7_export_jtrawl_grit.sh")
  echo "  (B) $lab step7 reference $JR, lock S0 $JV ; (C) annual step7 $JN"
  IDS="$IDS $JR $JV $JN"; DEP_D="$DEP_D:$JR:$JV:$JN"
done
JD=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany"$DEP_D" --job-name=svrForcing --qos=grit_normal \
       --mem=16G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/forcing.out" --error="$BASE/logs/forcing.err" \
       --wrap="export $ENVR; cd $CODE && Rscript svr_forcing.R $BASE $RUN")
echo "  (D) forcing locks: $JD"; IDS="$IDS $JD"
DEP_F=""
for lab in smooth flag; do
  [ "$DRY" -eq 1 ] || mkdir -p "$BASE/ocim_${lab}_an"
  JO=$(sb --dependency=afterok:"$JD" --job-name="svrCO2_$lab" --nodelist=hpc-08.grit.ucsb.edu \
         --export=ALL,FORCING_DIR="$BASE/step7_${lab}_an",OCIM_OUT_DIR="$BASE/ocim_${lab}_an" \
         "$DIR/co2model_batch_grit.sh")
  echo "  (E) $lab OCIM: $JO"; IDS="$IDS $JO"; DEP_F="$DEP_F:$JO"
done
JL=$(EXCL_JOB="$EXCL_R" sb --dependency=afterany"$DEP_F" --job-name=svrRead --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=00:30:00 \
       --output="$BASE/logs/read_svr.out" --error="$BASE/logs/read_svr.err" \
       --wrap="cd $CODE && python3 read_svr.py")
echo "  (F) reading: $JL"; IDS="$IDS $JL"

if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_svr_ocim.txt" <<EOT
# Annual cap, dedicated OCIM run, submitted on $(date)   BASE=$BASE
# ids :$IDS
# cancel everything:  scancel$IDS
EOT
echo "  ids written to $CODE/ids_svr_ocim.txt"
echo "  Read: $BASE/logs/A.out, S0_*.out, forcing.out ; then, after OCIM, read_svr.out"
