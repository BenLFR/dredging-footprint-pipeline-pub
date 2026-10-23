#!/bin/bash
# ============================================================================
# launch_hlin.sh -- new test of H-lin, set-ups A + B (rule written beforehand).
#   A: OCIM on the K = 1 forcings already built by D3 (d3_m4_20260928/step7_{smooth,flag}_l100_k1);
#   B: hlin_partition.py splits the smooth M4 forcing into dominant column / rest, then OCIM on each part;
#   then read_hlin.py, with afterany.
# Outputs in BASE, NOT in the D3 folder: read_d3.py scans d3_m4_20260928/ocim_* and
# would abort on an OCIM folder that is still empty.
# OCIM on hpc-08 (--nodelist), as in launch_d3_step2.sh; the Python jobs exclude it (EXCL of launch_threshold_step2.sh).
# DRY=1 by default.
# ============================================================================
set -euo pipefail
DRY="${DRY:-1}"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
RUN=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r6_flux10
BASE_D3=${DFP_SCRATCH:?set DFP_SCRATCH}/d3_m4_20260928
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/hlin_20260928
CODE="${DFP_REPO:?set DFP_REPO}/assessment/hlin"
PART=grit_nodes
# OCIM: list of launch_d3_step2.sh (hpc-08 is not in it, since it is imposed by --nodelist)
EXCL_OCIM=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu
# Python: EXCL line of launch_threshold_step2.sh (hpc-08 reserved for OCIM)
EXCL_PY=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu

ko() { echo "ABORT: $*"; exit 1; }
echo "== guards"
[ -f "$DIR/co2model_batch_grit.sh" ] || ko "co2model_batch_grit.sh missing"
for f in hlin_partition.py read_hlin.py; do [ -f "$CODE/$f" ] || ko "$CODE/$f missing"; md5sum "$CODE/$f"; done
md5sum "$DIR/co2model_batch_grit.sh"
[ -f "$BASE_D3/etape_S.csv" ] || ko "$BASE_D3/etape_S.csv missing"
for lab in smooth flag; do
  n=$(ls "$BASE_D3/step7_${lab}_l100_k1"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] || ko "$n .mat in $BASE_D3/step7_${lab}_l100_k1 (1 expected)"
  n=$(ls "$RUN/step7_${lab}_m4"/jdredge_ocim2_48l_*.mat 2>/dev/null | wc -l)
  [ "$n" -eq 1 ] || ko "$n .mat in $RUN/step7_${lab}_m4 (1 expected)"
  [ -d "$RUN/ocim_${lab}_m4" ] || ko "$RUN/ocim_${lab}_m4 missing"
done
[ -e "$BASE/PARTITION_PASS" ] && ko "PARTITION_PASS already exists"
for o in ocim_smooth_l100_k1 ocim_flag_l100_k1 ocim_smooth_colA ocim_smooth_resteB; do
  [ -e "$BASE/$o" ] && ko "$BASE/$o already exists"
done
echo "  guards: PASS"

sb() {
  local ex="$1"; shift
  if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch --parsable --partition=$PART --exclude=$ex $*" >&2; echo "DRY"
  else sbatch --parsable --partition="$PART" --exclude="$ex" "$@" | cut -d';' -f1; fi
}
[ "$DRY" -eq 1 ] || mkdir -p "$BASE/logs"

# ---- B: partition -------------------------------------------------------------
JP=$(sb "$EXCL_PY" --job-name="hlin_partition" --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=01:00:00 \
       --output="$BASE/logs/partition.out" --error="$BASE/logs/partition.err" \
       --wrap="cd $CODE && python3 hlin_partition.py $RUN $BASE")
echo "  partition: $JP"

# ---- OCIM ----------------------------------------------------------------------
IDS="$JP"; DEP=":$JP"
ocim() {   # <name> <forcing folder> [dependency]
  local item_name="$1" fdir="$2" dep="${3:-}" odir="$BASE/ocim_$1" J
  [ "$DRY" -eq 1 ] || mkdir -p "$odir"
  if [ -n "$dep" ]; then
    J=$(sb "$EXCL_OCIM" --dependency=afterok:"$dep" --job-name="hlin_$item_name" --nodelist=hpc-08.grit.ucsb.edu \
           --export=ALL,FORCING_DIR="$fdir",OCIM_OUT_DIR="$odir" "$DIR/co2model_batch_grit.sh")
  else
    J=$(sb "$EXCL_OCIM" --job-name="hlin_$item_name" --nodelist=hpc-08.grit.ucsb.edu \
           --export=ALL,FORCING_DIR="$fdir",OCIM_OUT_DIR="$odir" "$DIR/co2model_batch_grit.sh")
  fi
  echo "  ocim_$item_name: OCIM $J  (forcing $fdir)" >&2
  echo "$J"
}
J1=$(ocim smooth_l100_k1 "$BASE_D3/step7_smooth_l100_k1")
J2=$(ocim flag_l100_k1 "$BASE_D3/step7_flag_l100_k1")
J3=$(ocim smooth_colA "$BASE/step7_smooth_colA" "$JP")
J4=$(ocim smooth_resteB "$BASE/step7_smooth_resteB" "$JP")
IDS="$IDS $J1 $J2 $J3 $J4"; DEP="$DEP:$J1:$J2:$J3:$J4"

# ---- reading -------------------------------------------------------------------
JL=$(sb "$EXCL_PY" --dependency=afterany"$DEP" --job-name="hlin_read" --qos=grit_normal \
       --mem=32G --cpus-per-task=1 --time=01:00:00 \
       --output="$BASE/logs/read_hlin.out" --error="$BASE/logs/read_hlin.err" \
       --wrap="cd $CODE && python3 read_hlin.py $RUN $BASE_D3 $BASE")
echo "  reading: $JL"
IDS="$IDS $JL"
if [ "$DRY" -eq 1 ]; then echo; echo "Nothing was submitted (DRY=1)."; exit 0; fi
cat > "$CODE/ids_hlin.txt" <<EOT
# H-lin test A + B, submitted on $(date)   BASE=$BASE
# ids: $IDS
# cancel everything:  scancel $IDS
# NB: if the partition fails, the B OCIM runs stay in DependencyNeverSatisfied: cancel them by hand.
EOT
echo "  ids written to $CODE/ids_hlin.txt"
