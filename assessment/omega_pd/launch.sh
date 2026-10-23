#!/bin/bash
# No submission until the rule (RULE.md) is approved. DRY default is mandatory.
set -euo pipefail
DRY="${DRY:-1}" # 1 = print only, 0 = submit after approval.
T="${DFP_REPO:?set DFP_REPO}/assessment/omega_pd"
OUT=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_incert/omega_pd
EXCL=hpc-08.grit.ucsb.edu,hpc-10.grit.ucsb.edu,hpc-11.grit.ucsb.edu,hpc-12.grit.ucsb.edu,hpc-13.grit.ucsb.edu,hpc-17.grit.ucsb.edu,hpc-18.grit.ucsb.edu # Excluded nodes.
main() {
  case "$DRY" in 0|1) ;; *) echo 'DRY must be 0 or 1'; return 1;; esac
  cd "$T"
  python3 -B prepare.py --check
  md5sum --check bundle.md5
  local digest
  digest=$(md5sum RULE.md | cut -d' ' -f1) # First field of md5sum output.
  # Design guard: approval must name the exact version of the rule.
  if [ "$DRY" = 0 ] && [ "${APPROVED_RULE_MD5:-}" != "$digest" ]; then
    echo 'ABORT: approved rule fingerprint is absent or differs'; return 1
  fi
  local -a cmd=(sbatch --parsable --partition=grit_nodes --qos=grit_normal
    --exclude="$EXCL" --job-name=omega_pd
    --cpus-per-task=1 # Design choice, no source: sequential small-table recalculation.
    --mem=8G --time=00:30:00 # Design choice, no source: resource ceiling, not measured runtime.
    --chdir="$T" --output="$T/audit/compute_%j.out" --error="$T/audit/compute_%j.err"
    --export=ALL,APPROVED_RULE_MD5="$digest" "$T/compute.sbatch")
  if [ "$DRY" = 1 ]; then
    printf 'DRY: '; printf '%q ' "${cmd[@]}"; printf '\n'
    echo 'No job submitted. Only the project lead can validate the rule.'
    return
  fi
  # Atomic submission claim; never reuse/overwrite it after a partial failure.
  mkdir "$T/submission_claim"
  printf '%s\n' "$digest" > "$T/submission_claim/approved_rule.md5"
  "${cmd[@]}" | tee "$T/submission_claim/job_id.txt"
}
# Filter all displayed/logged command output; pipefail retains failures.
main "$@" 2>&1 | grep -viE 'token|secret|key|password|client_id|auth'
