#!/bin/bash
# Read saved reports only: never R, never sbatch, never modify a result.
set -euo pipefail
OUT=${DFP_SCRATCH:?set DFP_SCRATCH}/output_V6_incert/omega_pd
main() {
  if [ -f "$OUT/FAILURE.md" ]; then cat "$OUT/FAILURE.md"; fi
  if [ ! -f "$OUT/LOCK_PASSED" ]; then
    echo 'Lock has not passed. Do not interpret a sensitivity result.'
    [ ! -f "$OUT/locks.tsv" ] || cat "$OUT/locks.tsv"
    return 1
  fi
  [ -f "$OUT/COMPLETE" ] || { echo 'Incomplete run: no substantive interpretation'; return 1; }
  (cd "$OUT" && md5sum --check results.md5)
  cat "$OUT/RESULT.md"
  cat "$OUT/verdicts.tsv"
  cat "$OUT/results.tsv"
  cat "$OUT/interaction.tsv"
  cat "$OUT/arms.tsv"
  cat "$OUT/diagnostics.tsv"
  echo 'Copy LOG_ENTRY.md into the project log immediately, before interpretation.'
}
main "$@" 2>&1 | grep -viE 'token|secret|key|password|client_id|auth'
