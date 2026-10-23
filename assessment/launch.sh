#!/bin/bash
# Launches merge + step6 for each test case, with the launchers that are authoritative
# (merge_rerun.sbatch -> step5_merge_tiles_corr.R ; step6_pr.sbatch -> step6).
# Environment of _r6 (its merge-to-OCIM chain, not in this repository); only FI_PARAM_YAML / FI_FRAC_FILE change.
# DRY=1 by default: prints the sbatch commands without submitting anything.
set -euo pipefail
DRY="${DRY:-1}"
T="${DFP_REPO:?set DFP_REPO}/assessment"
DIR="${DFP_REPO:?set DFP_REPO}/pipeline"
R2=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_sala2021_r2
REF="${DFP_DATA:?set DFP_DATA}/output_V6_rejeu"
BASE=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_incert
Y0="${DFP_REPO:?set DFP_REPO}/config/fi_parameters_with_freshness.yaml"
PART=grit_nodes
EXCL=hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu  # hpc-17: node excluded
IDS="$T/ids.txt"

( cd "$DIR" && bash check_pipeline.sh >/dev/null 2>&1 ) || { echo "ABORT: check_pipeline.sh failed"; exit 1; }
[ -e "$BASE" ] && [ -n "$(ls -A "$BASE" 2>/dev/null)" ] && { echo "ABORT: $BASE not empty"; exit 1; }

# case -> "YAML|frac suffix" (empty suffix = original table)
declare -A CASES=( [ref]="$Y0|" [fill0.9]="$Y0|_f0.9" [fill0.8]="$Y0|_f0.8" [fill0.7]="$Y0|_f0.7"
                 [om0.1]="$T/yaml/omega_0.1.yaml|" [om1.0]="$T/yaml/omega_1.0.yaml|" [om5.0]="$T/yaml/omega_5.0.yaml|" )
sb() { if [ "$DRY" -eq 1 ]; then echo "  [dry] sbatch $*" >&2; echo "DRY"; else sbatch --parsable --partition="$PART" --exclude="$EXCL" "$@" | cut -d';' -f1; fi; }

[ "$DRY" -eq 1 ] || { mkdir -p "$BASE"; echo "# uncertainty tests submitted on $(date)" > "$IDS"; }
for c in ref fill0.9 fill0.8 fill0.7 om0.1 om1.0 om5.0; do
  y="${CASES[$c]%%|*}"; s="${CASES[$c]##*|}"
  RUN="$BASE/$c"
  for lab in smooth flag; do
    if [ -z "$s" ]; then fr="$REF/decharge/fraction_sur_place_${lab}.csv"; else fr="$T/frac/fraction_sur_place_${lab}${s}.csv"; fi
    [ -f "$y" ] && [ -f "$fr" ] || { echo "ABORT: missing input for $c/$lab"; exit 1; }
    JM=$(sb --job-name="inc_m_${c}_${lab}" \
         --export=ALL,DREDGE_LABEL=$lab,MERGE_IN_DIR="$R2/step5_$lab",MERGE_OUT_DIR="$RUN/merge_${lab}_pr",FI_FRAC_FILE="$fr",FI_PARAM_YAML="$y",MERGE_N_EXPECTED=71 \
         "$DIR/merge_rerun.sbatch")
    J6=$(sb --job-name="inc_6_${c}_${lab}" --dependency=afterok:"$JM" \
         --export=ALL,DREDGE_LABEL=$lab,PIPE_RUN_DIR="$RUN" "$DIR/step6_pr.sbatch")
    echo "$c $lab: merge $JM -> step6 $J6"
    [ "$DRY" -eq 1 ] || echo "$c $lab $JM $J6" >> "$IDS"
  done
done
