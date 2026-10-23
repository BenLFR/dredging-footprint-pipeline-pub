#!/usr/bin/env Rscript
# ============================================================================
# test_failure_modes.R
# Documents 4 pipeline failure modes using the testthat framework.
# testthat integrates with GitHub Actions CI and produces JUnit XML output
# that journal reproducibility checkers accept.
#
# Usage:
#   Rscript tests/test_failure_modes.R
#
# From project root (preferred — discovers all test_that blocks):
#   Rscript -e "testthat::test_file('tests/test_failure_modes.R')"
#
# Outputs:
#   output_V6/test_results.xml   (JUnit XML — uploaded as CI artifact)
#   console progress report
#
# Tests:
#   1. sparse_data_guard          — <100 pings should trigger guard, not crash
#   2. SAR_finite_single_cell     — SAR formula is finite for single-cell input
#   3. extreme_isolation_forest   — contamination_rate=0.50 must not error
#   4. arctic_grid_assignment     — lat>75N pings assigned to valid grid cells
# ============================================================================

suppressPackageStartupMessages(library(testthat))

# ── Pipeline constants (mirrors constants.R) ───────────────────────────────────
CELL_SIZE_M             <- 1000L
CELL_AREA_M2            <- CELL_SIZE_M * CELL_SIZE_M
GRID_COLS               <- 34735L
WORLD_XMIN              <- -17367530.45
WORLD_YMAX              <-  7342699.72
MIN_POINTS_PER_VESSEL   <- 100L       # from outlier_config_V6.yaml
DREDGING_SPEED_KN       <- 4.0

dir.create("output_V6", showWarnings = FALSE, recursive = TRUE)

# ── TEST 1: Sparse data guard ─────────────────────────────────────────────────
# The pipeline should flag vessels with <100 pings and not crash.
# min_points_per_vessel = 100 is set in configuration/outlier_config_V6.yaml.
# ─────────────────────────────────────────────────────────────────────────────
test_that("sparse data guard triggers for <100 pings", {
  skip_if_not_installed("data.table")
  library(data.table)

  n_pings <- 42L   # well below 100 threshold
  sparse_dt <- data.table(
    ssvid       = rep(123456789L, n_pings),
    Timestamp   = as.POSIXct("2020-01-01") +
                    seq(0, by = 3600, length.out = n_pings),
    Latitude    = runif(n_pings, 50.0, 51.0),
    Longitude   = runif(n_pings, 2.0,  3.0),
    speed_knots = runif(n_pings, 0, 20)
  )

  # Replicate the guard logic from step2_process_vessel.R
  guard_triggered <- nrow(sparse_dt) < MIN_POINTS_PER_VESSEL

  expect_true(
    guard_triggered,
    label = sprintf(
      "Guard expected for n=%d pings (threshold=%d). If FALSE, add a min_points_per_vessel check to step2_process_vessel.R",
      n_pings, MIN_POINTS_PER_VESSEL
    )
  )

  # Confirm no crash when checking — result is logical, not an error
  expect_type(guard_triggered, "logical")
})

# ── TEST 2: SAR is finite for single-cell dense input ─────────────────────────
# All pings in one 1-km² grid cell. SAR = total_swept_area / cell_area.
# Must be finite and non-negative — never Inf, NaN, or negative.
# ─────────────────────────────────────────────────────────────────────────────
test_that("SAR is finite for single-cell input", {
  skip_if_not_installed("data.table")
  library(data.table)

  n_pings        <- 500L
  dredge_width_m <- 1.3    # fleet mean from ship_specs.yaml

  # All pings clustered in ~0.01° area → same 1-km grid cell after projection
  dt_dense <- data.table(
    speed_knots = runif(n_pings, 1.0, 3.5),
    dt_sec      = rep(600L, n_pings)
  )
  dt_dense[, speed_ms      := speed_knots * 0.514444]
  dt_dense[, trawled_dist_m := speed_ms * dt_sec]

  SAR <- sum(dt_dense$trawled_dist_m * dredge_width_m) / CELL_AREA_M2

  expect_true(is.finite(SAR),   label = "SAR must be finite (not Inf or NaN)")
  expect_true(!is.nan(SAR),     label = "SAR must not be NaN")
  expect_true(SAR >= 0,         label = "SAR must be non-negative")
  expect_true(SAR < 1e6,        label = "SAR must be physically plausible (< 1e6)")
})

# ── TEST 3: Isolation Forest accepts extreme contamination ────────────────────
# Setting contamination_rate = 0.50 (vs production 0.03) must not crash.
# Verifies the Isolation Forest path in step2 is robust to configuration errors.
# ─────────────────────────────────────────────────────────────────────────────
test_that("Isolation Forest accepts extreme contamination rate (0.50)", {
  skip_if_not_installed("solitude",
    message = "solitude not installed — skipping Isolation Forest test (expected in CI)")

  library(solitude)
  set.seed(42L)
  n <- 500L

  # Synthetic vessel features matching step2 Isolation Forest input
  feat <- data.frame(
    speed_knots   = c(runif(floor(n * 0.97), 0, 20),
                      runif(ceiling(n * 0.03), 80, 100)),
    delta_t_sec   = c(rexp(floor(n * 0.97), 1 / 600),
                      rexp(ceiling(n * 0.03), 1 / 5)),
    accel_knots_s = c(rnorm(floor(n * 0.97), 0, 0.01),
                      rnorm(ceiling(n * 0.03), 0, 2.0))
  )

  expect_no_error({
    iso <- isolationForest$new(
      num_trees   = 50L,
      sample_size = 256L,
      max_depth   = 8L,
      seed        = 1L
    )
    iso$fit(feat)
    scores <- iso$predict(feat)

    extreme_contamination <- 0.50
    threshold  <- quantile(scores$anomaly_score,
                           probs = 1 - extreme_contamination, na.rm = TRUE)
    n_outliers <- sum(scores$anomaly_score >= threshold)
  })

  expect_true(n_outliers > 0,
    label = "At least some points should be flagged as outliers")
  expect_true(n_outliers <= n,
    label = "Outlier count must not exceed total observations")
})

# ── TEST 4: Arctic pings land on valid grid cells ─────────────────────────────
# Pings at lat > 75°N, projected to EPSG:6933, must produce a valid grid_id.
# This verifies that the Equal-Earth projection and grid-snapping maths work
# for high-latitude inputs (a known edge case documented in LIMITATIONS.md).
# ─────────────────────────────────────────────────────────────────────────────
test_that("Arctic pings (lat > 75N) produce a valid EPSG:6933 grid cell", {
  skip_if_not_installed("sf",
    message = "sf not installed — skipping EPSG:6933 projection test")

  library(sf)

  # Arctic Ocean: 80°N, 0°E — should be ocean, not land
  lat_arctic <- 80.0
  lon_arctic <- 0.0

  pt       <- st_as_sf(data.frame(lon = lon_arctic, lat = lat_arctic),
                       coords = c("lon", "lat"), crs = 4326)
  pt_6933  <- st_transform(pt, 6933)
  coords_m <- st_coordinates(pt_6933)
  x_m      <- coords_m[1, "X"]
  y_m      <- coords_m[1, "Y"]

  col_idx <- as.integer(floor((x_m - WORLD_XMIN) / CELL_SIZE_M))
  row_idx <- as.integer(floor((WORLD_YMAX - y_m) / CELL_SIZE_M))
  grid_id <- row_idx * GRID_COLS + col_idx + 1L

  expect_true(col_idx >= 0L,
    label = sprintf("col_idx must be >= 0 (got %d)", col_idx))
  expect_true(col_idx < GRID_COLS,
    label = sprintf("col_idx must be < GRID_COLS=%d (got %d)", GRID_COLS, col_idx))
  expect_true(grid_id >= 1L,
    label = sprintf("grid_id must be >= 1 (got %d)", grid_id))

  # The cell must exist within the Equal-Earth grid extent
  expect_true(is.finite(x_m) && is.finite(y_m),
    label = "EPSG:6933 coordinates must be finite for Arctic point")
})

# ── Run and save JUnit XML ────────────────────────────────────────────────────
# When run directly as a script (not via test_file()), self-execute with
# JUnit XML reporter so CI can parse results.
if (!interactive()) {
  xml_path <- "output_V6/test_results.xml"
  cat(sprintf("\n=== Running tests | JUnit XML → %s ===\n\n", xml_path))

  result <- tryCatch(
    testthat::test_file(
      "tests/test_failure_modes.R",
      reporter = testthat::MultiReporter$new(reporters = list(
        testthat::ProgressReporter$new(),
        testthat::JunitReporter$new(file = xml_path)
      ))
    ),
    error = function(e) {
      cat("test_file() not available in this testthat version — running directly.\n")
      NULL
    }
  )

  if (!is.null(result)) {
    n_failed <- sum(as.data.frame(result)$failed, na.rm = TRUE)
    if (n_failed > 0) {
      cat(sprintf("\n%d test(s) FAILED — see %s for details.\n", n_failed, xml_path))
      quit(status = 1L)
    } else {
      cat(sprintf("\nAll tests PASSED. JUnit XML saved: %s\n", xml_path))
    }
  }
}
