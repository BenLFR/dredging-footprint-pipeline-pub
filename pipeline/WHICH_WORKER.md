# Which step5 worker to run? — updated 2026-08-13

## In service

| run | with | why |
|---|---|---|
| **`step5_tile_worker_batch.R`** | **`step5_batch.sbatch`** | intersection **in batches of 250**; ~7 min/tile |

`step5_batch.sbatch` accepts `PIPE_OUT_SUFFIX` to write somewhere other than the
production folder (used for the equivalence test).

## Kept on the cluster, do NOT run (not part of this repository)

| file | status |
|---|---|
| `step5_tile_worker_corrected.R` | correct but **st_intersects row by row**: 5.57 s/row, 74-94 h per heavy tile. Validation reference. |
| `step5_rejeu.sbatch`, `step5_rejeu_long.sbatch` | call the worker above |
| `echafaudage/step5_tile_worker_corr.R` | **origin of the batch optimisation**; diverges elsewhere (no option C, no declared input) |
| `echafaudage/step5_tile_worker_seuil35*.R` | threshold trials, not retained |

## The only difference between the two workers in use

**Section 7** only: 24 lines removed, 29 added. Everything else —
option C (`n_pipes_eff`), input declared via `PIPE_*_FILE`, refusal instead of a
wrong fallback — is identical.

## Proof of equivalence (2026-08-13)

| tile | cells loop / batches | max relative difference | loop | batches |
|---|---|---:|---:|---:|
| 379 `flag` | 267 / 267 | **0.000e+00** | 15 h 06 | **7 min 03** |
| 380 `smooth` | 225 / 225 | **0.000e+00** | 16 h 08 | **7 min 00** |

Script: `~/equivalence_batch.R` (not part of this repository).

## ⬜ To do once no job reads `step5_tile_worker_batch.R` any more

Fold the batched section 7 **into** `step5_tile_worker_corrected.R` so that there is
**a single** production worker, and turn this file into a historical
note. Until then, two files coexist and this README
is the only thing that says which one to use.
