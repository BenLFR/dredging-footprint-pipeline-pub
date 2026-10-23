#!/bin/bash
# Only launch.sh submits this file. Explicit direct submission also needs approval.
set -Eeuo pipefail
T="${DFP_REPO:?set DFP_REPO}/assessment/omega_pd"
OUT=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_incert/omega_pd
failed_job() {
  # Shell/library failure may happen before the R handler; preserve a six-field record.
  local target="$T/audit/failure_job_${SLURM_JOB_ID:-unknown}.md"
  if [ -d "$T/audit" ] && [ ! -e "$target" ]; then
    (set -o noclobber; printf '%s\n' \
      '# Incomplete job: no interpretation' \
      '1. Question: can the calculation complete all its locks and contrasts?' \
      '2. Rule: any source, allocation or library error stops the test.' \
      "3. Raw result: job ${SLURM_JOB_ID:-unknown}, see the filtered log compute_<job>.out and FAILURE.md if present." \
      '4. Verdict: INSTRUMENT_FAILURE_OR_INCOMPLETE_RUN.' \
      '5. Provenance: inputs.tsv, bundle.md5, job log and RULE.md.' \
      '6. Paper: nothing. Copy this record into the project log before resuming.' > "$target")
  fi
}
trap failed_job ERR
main() {
  : "${SLURM_JOB_ID:?R must run in a Slurm allocation}"
  cd "$T"
  [ "${SLURM_JOB_PARTITION:-}" = grit_nodes ] || { echo 'Wrong partition'; return 1; }
  case "$(hostname -s)" in hpc-08|hpc-10|hpc-11|hpc-12|hpc-13|hpc-17|hpc-18)
    echo 'Excluded host'; return 1;; esac # Excluded nodes; failure status.
  [ "${APPROVED_RULE_MD5:-}" = "$(md5sum RULE.md | cut -d' ' -f1)" ] || {
    echo 'Rule approval missing or stale'; return 1; } # First md5sum field; failure status.
  python3 -B prepare.py --check
  md5sum --check bundle.md5
  # Recheck at job start; a second job cannot overwrite this one.
  [ ! -L "$OUT" ] || { echo 'Output symlink refused'; return 1; }
  if [ -d "$OUT" ] && [ -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo 'Nonempty output refused'; return 1 # Depth one = direct children; failure status.
  fi
  mkdir -p "$OUT"
  mkdir "$OUT/run_claim"
  printf '%s\n' "$SLURM_JOB_ID" > "$OUT/run_claim/job_id.txt"
  cp "$T/RULE.md" "$T/inputs.tsv" "$T/bundle.md5" "$OUT/run_claim/"
  mkdir "$OUT/tmp"
  export TMPDIR="$OUT/tmp" R_LIBS_USER="$HOME/R/library"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 # Design: one CPU, no threaded oversubscription.
  Rscript --vanilla "$T/compute.R"
}
main "$@" 2>&1 | grep -viE 'token|secret|key|password|client_id|auth'
