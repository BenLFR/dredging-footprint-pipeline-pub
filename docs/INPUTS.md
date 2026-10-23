# Inputs and access

No data download or account login is performed by the source checks. Obtain the exact study inputs through
their providers and preserve their versions and checksums; a current replacement is not automatically the
input that produced the reported runs. See [third-party notices](../THIRD_PARTY_NOTICES.md).

| Input | Consumer | Required location or interface |
|---|---|---|
| AIS positions and vessel metadata | Steps 0–3 | `AIS_INPUT_FILE`; the split launcher expects `DFP_DATA/AIS_data/benjamin3_clean.csv` unless adapted |
| Lithology rasters and their exact source metadata | Step 4 | `DFP_DATA/configuration/lithology/`; cache under `DFP_DATA/hubocean_cache/` |
| Carbon-stock and Longhurst layers; land mask | Steps 5–6 | The paths configured by the scripts/environment; files are external, not bundled |
| GEBCO bathymetry | Step 7 | `DFP_SCRATCH/configuration/gebco_2026/GEBCO_2026.nc` or the `GEBCO_FILE` override in step 7 |
| OCIM2-48L matrix, WOA09 MAT inputs and helper functions | Step 7 / OCIM | `DFP_DATA/configuration/ocim/`; exact files listed below |
| In-place sediment fraction | Merge and assessment | `FI_FRAC_FILE`; archived `DFP_DATA/output_V6_rejeu/decharge/fraction_sur_place_<label>.csv` |
| Reference and scenario run outputs | Annual cap and assessment | Dated run folders beneath `DFP_SCRATCH`, as named by each launcher |

For the OCIM run, obtain `OCIM2_48L_CTL.mat`, `woa09si.mat`, `woa09po4.mat`, `netemission.txt` and the MATLAB
dependencies required by `ocim/co2model_batch.m`. The WOA09 `.mat` files are model-prepared inputs, not a claim
that the NOAA NetCDF downloads have that format. `ocim_cache.mat`, consumed by step 7, is an external prepared
cache; its generator is not included here. Request the exact model/cache preparation provenance with the inputs.

Copy the published `ocim/co2model_batch.m` to the OCIM data working directory when using the existing launchers:
they execute it there alongside the external model files. The launchers' MD5 guard must match that copied file.

## Local annual-cap input tree

`annual_svr/s_targets_20261005.R` uses `DFP_LOCAL` as its working directory. It reads
`data/cri_final_<label>_20260925_125232.rds` and `annual_svr/annee_<label>.csv`, then writes
`annual_svr/ocim/s_cibles_<label>.csv`. The annual CSVs come from `annual_svr/l2_lock.R` on the cluster;
transfer them to this local tree before running the script. Transfer the resulting scaling tables to the
cluster's `annual_svr/ocim/` working tree before running the annual OCIM launchers. These generated CSVs
are not repository inputs supplied by this release.

## Known reconstruction limits

The generating script for the archived `fraction_sur_place_<label>.csv` input has not yet been identified.
Those tables are required to reproduce the current merge. The code therefore does not yet provide an
end-to-end reconstruction of that input from raw AIS alone.

The role of the older `config/fi_parameters.yaml` relative to `config/fi_parameters_with_freshness.yaml`
still needs reconciliation; the latter is the configuration pinned by the October assessment guards.
Additional historical analysis folders are not automatically covered by this source tree.
