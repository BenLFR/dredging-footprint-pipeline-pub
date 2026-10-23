# Scientific source map

This map links paper components to code entry points. `MANIFEST.tsv` provides the per-file provenance and
checksums. Figure/table numbering is deliberately omitted because the manuscript is still being revised.
The named scripts do not supply the external data or archived run outputs listed in [INPUTS.md](INPUTS.md).

| Paper component | Public code | Source/output role |
|---|---|---|
| Coverage window | `pipeline/step0_core_window_enhanced.R` | Coverage matrix and selected vessel/year window |
| Cleaning and activity classification | `pipeline/step2_process_vessel.R`, `pipeline/step3_merge_final.R` | Retained AIS tracks and dredging labels |
| Classification diagnostics | `pipeline/gmm_diagnostics.R` | GMM archive and diagnostic summaries |
| Lithology assignment | `pipeline/step4_add_lithology_vNext.R` | Sediment memberships and effective labile fraction |
| Disturbed sediment and remineralisation | `pipeline/step5_tile_worker_batch.R`, `pipeline/step5_merge_tiles_corr.R` | Per-cell SAR/SVR and remineralised fraction |
| Carbon stock and carbon estimate | `pipeline/step6_calculate_cri_corr.R` | Cell carbon stocks and first-year carbon estimates |
| OCIM mapping and injection | `pipeline/step7_export_jtrawl.R` | Carbon forcing distributed onto the OCIM grid |
| OCIM integration | `ocim/co2model_batch.m`, `pipeline/co2model_batch_grit.sh` | Control and forcing scenarios; requires external model |
| Annual SVR cap | `annual_svr/worker_year.R`, `annual_svr/l1_tiles.R`, `annual_svr/l2_lock.R`, `annual_svr/s_targets_20261005.R` | Annual cell contributions and scaling tables |
| Annual carbon and OCIM forcing | `annual_svr/ocim/make_cri_year.R`, `annual_svr/ocim/svr_forcing.R`, `annual_svr/ocim/read_svr.py`, `annual_svr/ocim/read_lin10.py` | Annual cap checks and OCIM reads |
| Finite forcing window | `annual_svr/ocim/launch_ocim_window.sh`, `annual_svr/ocim/read_window.py` | Dredging-window scenarios and reference checks |
| Omega and hopper fraction | `assessment/prep_cases.sh`, `assessment/read.R`, `assessment/omega_pd/` | Existing scenario reads and locked Omega × excavation-depth recalculation |
| Kinetic/province mapping | `assessment/read_K.R` | Checks using the lookup table extracted from the verified merge |
| Drag-Score threshold | `assessment/threshold/` | Threshold arms, reruns and reference checks |
| Global sensitivity | `assessment/sobol_r6/` | Sobol analyses, including the annual cap variant |
| Injection/redistribution checks | `assessment/hlin/`, `assessment/d3_m4/`, `assessment/window_remeasure/` | Partition, moved forcing, injection-depth and forcing-window comparisons |
| Atmospheric-fraction series | `assessment/window_remeasure/extract_af_annual.py`, `assessment/window_remeasure/extract_af_variants.py` | Annual AF extraction from external OCIM outputs |
| Summary/post-processing | `assessment/test_q/summary_q.py`, `pipeline/postproc_atwood_style.py` | Summary reads and mapped outputs |

The `pipeline/lock_*.R` files check reference behaviour for selected changes. They are distinct from the
CI syntax checks and generally require the archived study runs. Do not treat a source-check PASS as a pass
of those scientific locks.
