# Glossary of names kept in French

The code, comments, messages, job names and log names of this repository are in English. A small set of
names is still in French, on purpose. They are names that the code **writes to disk and then reads back**:
column names, file names, folder names. The runs reported in the paper produced files that carry these
names. Renaming them in the code would stop it from reading those outputs again, so they keep the names
under which they were written.

This page lists each of those names with its meaning. The rule used to sort names:

| kept in French | translated |
|---|---|
| a column, file or folder name written by one script and read by another | a variable, function or argument name |
| a dated run folder or file of a reported run | a pass/fail marker file (`MOVE_PASS`, `FORCING_PASS`) |
| a variable stored in the OCIM `.mat` output files | a Slurm job name, a log file name |
| a name from an external source | an environment variable passed between scripts |
| a historical file name cited in documentation | a file written but never read back (`locks.csv` labels, record files) |

## 1. Columns of the AIS and grid tables

| name | meaning | written by | read by |
|---|---|---|---|
| `navire` | vessel name, as it appears in some raw AIS inputs (renamed to `Navire` on reading) | input data | `step0`, `step1` |
| `Navire` | vessel name | `step0`, `step1` | all later steps |
| `Annee`, `annee` | calendar year of the AIS position | `step0` onward | `step1` to `step3`, annual SVR scripts |
| `Dragage_flag` | dredging label (1 = dredging, 0 = other activity) | `step3` | `step5` worker, annual SVR worker, threshold test |
| `dredge_lisse` | dredging label after temporal smoothing (*lisse* = smoothed) | `step3` | `step3` |
| `lithologie` | seabed sediment class at the position | `step4` | `step5`, `annual_svr/worker_year.R` |
| `lithologie_canon` | sediment class mapped to the canonical classes | `step4` | `pipeline/lock_sala_step4.R` |
| `V_balaye` | swept volume of sediment in a cell (m³) (*balayé* = swept) | the in-place fraction table `fraction_sur_place_<label>.csv` | `assessment/prep_cases.sh`, Ω × p_d test |
| `V_retire` | volume removed from the cell into the hopper (m³) (*retiré* = removed) | same table | `assessment/prep_cases.sh`, Ω × p_d test |
| `V_retire_test` | `V_retire` recomputed by the Ω × p_d test | `assessment/omega_pd` | `assessment/omega_pd/compute.R` |
| `V_sur_place` | volume left in place in the cell (m³) | `assessment/prep_cases.sh` | Ω × p_d test |
| `frac_sur_place` | fraction of the swept volume left in place | `step5` merge | locks, Sobol and Ω × p_d tests |
| `f_i_full_avant_pr` | remineralised fraction f_i before the factor p_r (*avant* = before) | `step5` merge | locks, Ω × p_d test |
| `f_i_conservative_avant_pr` | the same for the conservative variant; no script reads it, it is kept so that the pair of columns of the merge output has one naming | `step5` merge | — |

## 2. Files and folders read back

| name | meaning |
|---|---|
| `decharge/fraction_sur_place_<label>.csv` | per-cell fraction left in place (*décharge* = discharge), input `FI_FRAC_FILE` of the merge |
| `fraction_sur_place_flag.csv` | the same file for the `flag` label |
| `output_V6_rejeu/` | output folder of the rerun of the production pipeline (*rejeu* = rerun) |
| `navire_<NN>_<name>.qs`, `navires_metadata.csv` | one file per vessel, and their index, written by `step1`, read by `step2` |
| `stage3A_fusion.rds` | checkpoint of `step3` (merged AIS table), written and read by `step3` to resume |
| `AIS_with_lithology_seuil_<arm>*` | `step4` output for one arm of the threshold test (*seuil* = threshold) |
| `annee_<label>.csv`, `sar_annee_<NNN>.parquet` | per-year rows of the target cells, and per-tile yearly SAR |
| `cibles_<label>.csv`, `s_cibles_<label>.csv` | target cells (SVR ≥ 1, *cibles* = targets) and their scaling s_i |
| `s_annuel_<label>.csv` | per-cell annual weights used by the capped Sobol analysis |
| `tuiles_<label>.txt` | tiles to rerun (*tuiles* = tiles) |
| `cri_an_<label>/` | annual C_ri files (*an* = year) read by `step7` through `CRI_FILE` |
| `compteurs.csv` | counts per arm of the threshold test (*compteurs* = counters); columns `bras` (arm), `seuil` (threshold value), `au_dessus` (positions at or above the threshold), `au_dessus_rel_t0_pct`, `changees` (positions whose label changed), `hors_71` (positions outside the 71 tiles) |
| `etape_S.csv` | step S of the D3 test (*étape* = step); columns `variante` (variant), `TV_colonne` (total variation by column), `part_1re_col` (share of the first column), `borne_pt` (bound, in AF points) |
| `choix_ocim.txt` | variants selected for the OCIM runs of the D3 test (*choix* = choice) |

## 3. Suffixes of run folders

| suffix | meaning |
|---|---|
| `_an`, `s_an` | annual SVR cap (*an* = year): `step7_<label>_an`, `ocim_<label>_an` |
| `_fen`, `s_fen` | forcing window of Atwood et al. (2024) (*fenêtre* = window): `ocim_<label>_fen` |
| `resteB`, `smooth_resteB` | part B of the H-lin partition, everything except column A (*reste* = rest) |
| `_sansfrac`, `t0_sansfrac` | merge run without the in-place fraction file (*sans* = without) |
| `_sansfraicheur` | merge run without the freshness factor |

## 4. Dated runs reported in the paper

These folders and files hold the outputs behind the numbers of the paper. Their names include the run
date and are kept exactly:
`svr_annuel_20261005`, `svr_annuel_ocim_20261005`, `s_an_20261005`, `_an_20261005`,
`sobol_plafond_20261005`, `seuil_20260928`, `ocim_fenetre_20261009`, `remesure_fen_20261010`,
`af_annuelle_20261013`, `af_variantes_20261016`, `jdredge_ocim2_48l_resteB_20260928`,
`ocim_smooth_resteB_incomplet_20260929`, `dragage_gridsearch_results_V6_20260810_093400`,
`profondeur_injection/candidat` (test O candidate forcing).

## 5. Variables of the OCIM output files

| name | meaning |
|---|---|
| `annee_pas` | year of each OCIM time step (*pas* = step) |
| `actif` | 1 if the dredging forcing is active at that step |

Both are written by `ocim/co2model_batch.m` into the `.mat` files and read by the Python readers.

## 6. External names

`si_an`, `po4_an` are variable names of the World Ocean Atlas 2009 files (annual silicate and phosphate).
`De` is part of the company name Jan De Nul in vessel sources.

## 7. Historical file names in documentation

`step2_process_navire_0303_tiled`, `step5_rejeu`, `step5_rejeu_long`, `step5_tile_worker_seuil35`,
`phase3_fraction_rejeu`, `step2_process_navire`, `ranger_pipeline` and `rapatrier_step1` are names of files
on the cluster where the pipeline ran. The last three mean vessel processing, pipeline organisation
and retrieval of step1 outputs, respectively. They are cited in
`pipeline/step2_PROVENANCE.md`, `pipeline/WHICH_WORKER.md` and in comments, to trace where a script
came from. They are not part of this repository.
