# step2_process_vessel.R — provenance (2026-09-14)

- `step2_process_vessel.R` (version `tiled-v8-onlandfp`) is the worker that produced `ais_results_14542`, the input of the replay.
  The file executed on 2026-01-27 had been rewritten on 2026-03-03 and had no git history; the only trace,
  a copy of `step2_process_navire.R` in the home directory on the cluster, was damaged (CRLF, `\n` and `\\` destroyed). It was rebuilt mechanically
  (by a repair script not included in this repository), sha256 71aa3fec…e637.
- Lock (GRIT jobs 1104253/1104254): on Vox Maxima, with the variables of
  `step2_process_array.sh` and `ship_specs.yaml.bak-20260810`, the 560 289 rows and all columns are
  identical to `ais_results_14542`, with an identical breakdown of removals. PASSED.
- `step2_process_navire_0303_tiled-v8.R`: the version of 2026-03-03 (commit 76f9504), which produced the published March run
  (`ais_results_43018`). Kept, not executed by the orchestrator.
- The restored file reads the current `ship_specs.yaml` (modified on 2026-08-10), not the backup used for the lock.
