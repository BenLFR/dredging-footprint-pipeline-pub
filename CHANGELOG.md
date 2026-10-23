# Changelog

## [2.0.0] — unreleased; prepared 2026-10-23

- Replace the March pipeline with the verified October GRIT source chain, annual SVR cap, OCIM integration
  and assessment scripts. Preserve the March snapshot under the local `v1.0.0` tag at `df60a43`.
- Translate messages and identifiers, rename the public layout, parameterise local paths and update source
  checksum guards. `MANIFEST.tsv` relates all 100 scientific sources to their originals; step 2 is reconstructed.
- Add input/access documentation, a script map, checksum explanations and explicit reconstruction limits.
- Replace the old toy smoke CI with source-integrity, syntax and manifest-checker tests. This is not a full
  numerical reproduction. Remove the March environment lock; GRIT dependency versions remain to be recovered.
- Correct repository URLs, align the four confirmed authors across citation metadata, remove invented DOI
  placeholders and regenerate the archive inventory. No v2 release, DOI or scientific data archive exists yet.

## Historical notes retained from the March snapshot

The following notes describe that earlier tree; their paths are not instructions for version 2.0.0.
They are preserved as historical records and have not been revalidated against the October chain.

## [1.0.0] — 2026 (pub/v1.0-clean)

### Added
- Canonical single-file-per-step pipeline structure under `pipeline/`
- `pipeline/constants.R`: shared EPSG:6933 grid constants
- `pipeline/step0_core_window.R`: enhanced coverage matrix with window scoring
- `pipeline/step1_split_navires.R/sh`: AIS split by MMSI
- `pipeline/step2_process_navire.R/sh`: per-vessel filtering (IF, speed, land-crossing)
- `pipeline/step3_merge_final.R/sh`: merge + 2-component GMM + DBSCAN classification
- `pipeline/step4_add_lithology.R/sh`: dbSEABED lithology join + penetrability index
- `pipeline/step5_make_tiles.R/step5_tile_worker.R/step5_tile_job.sh`: 630-tile parallel SAR
- `pipeline/step5_merge_tiles.R/step5_merge_slurm.sh`: tile fusion + Longhurst + CRI formula
- `pipeline/step6_calculate_cri.R/sh`: cumulative risk index
- `pipeline/step7_export_jtrawl.sh`: J_dredge export for OCIM
- `pipeline/step7_extract_ocim_cache.py/.m`: OCIM2 grid enrichment (KDTree, 9 fields)
- `config/outlier_config.yaml`: step2 filter parameters
- `config/ship_specs.yaml`: 10 TSHD vessel specifications (beam, hopper, speed)
- `postproc/`: Python/R post-processing and visualisation scripts
- `deploy/run_full_pipeline.sh`: autonomous steps 0-7 + CO2 GRIT orchestrator
- `data/README.md`: acquisition guide for all input datasets
- `data/toy/generate_toy_data.R`: synthetic North Sea TSHD toy dataset
- `docs/architecture.md`: data flow diagram, step I/O table, formulas
- `docs/data_policy.md`: AIS restrictions, processed output DOIs, L&O:Methods template
- `docs/reproducibility.md`: renv, Python, MATLAB, Apptainer, cluster setup
- `tests/smoke_test_steps2_5.sh`: end-to-end smoke test on toy data
- `.github/workflows/smoke_test.yml`: GitHub Actions CI
- `CITATION.cff`: machine-readable citation metadata
- `.gitattributes`: LF line endings for all shell/R/Python/MATLAB scripts

### Removed
- Legacy version suffixes (`_FIXED`, `_Fonctionnel`, `_HHM`, `_vNext`, `_TEST`)
- Windows batch launchers (`batch_windows/`, `*.bat`)
- Installation and maintenance scripts (`installation/`, `maintenance/`, `submission/`)
- Raw pipeline outputs (`output_V6/`, `GeoTIFF_Step*/`)
- Thesis Word documents and PDFs
- Vendored third-party MATLAB CO2 model source files (now external dependency)

### Fixed
- step3: pROC `coords()` returning data.frame caused `Dragage_flag` always 0
- step1: `qs::qsave` `preset_compression` -> `compress_level` argument rename
- step2: land-crossing filter applied to wrong coordinate column

---

## [0.x] — 2025 (feature/step2-extraction-visualization, main)

Development history on working branches. See git log for details.
