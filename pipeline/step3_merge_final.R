#!/usr/bin/env Rscript
# ============================================================================
# STEP-3  ─  MERGE  +  CLEANING  +  GMM / GRID-SEARCH   (pipeline  V6)
#   • Runs on a large node (≥ 8 CPU, 32 GB RAM)
#   • Re-assembles all the *_clean.rds produced at Step-2
#   • Finishes the pipeline: isolated spikes, context percentiles, DBSCAN,
#     stop-GMM, smoothing, grid-search, final save.
# ============================================================================

# =============================================================================
# STEP 3: MERGE AND GRID SEARCH - WITH CHECK-POINTING SYSTEM
# =============================================================================

split_job_id <- Sys.getenv("SPLIT_JOB_ID")
results_dir <- Sys.getenv("RESULTS_DIR")
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# DRY-RUN MODE DETECTION
is_dry <- Sys.getenv("DRY_RUN", "0") == "1"
if (is_dry) {
  cat("⚡ === DRY-RUN MODE ENABLED ===\n")
  cat("📝 Test on a reduced sample - safeguards relaxed\n")
}

cat("🚀 === STEP 3: MERGE AND GRID SEARCH (WITH CHECK-POINTING) ===\n")
cat("📦 Split job:", split_job_id, "\n")
cat("📁 Results folder:", results_dir, "\n")
cat("🖥️  Cores:", cores, "\n")
cat("⏱️  Start:", format(Sys.time()), "\n\n")

home_dir <- path.expand("~")
# Explicit output directory, so that a replay
# NEVER overwrites the published outputs. GMM_probabilities_QA.rds is not
# time-stamped and would be overwritten silently.
output_dir <- path.expand(Sys.getenv("PIPE_OUT_DIR", file.path(Sys.getenv("DFP_DATA"), "output_V6")))
if (identical(normalizePath(output_dir, mustWork = FALSE),
              normalizePath(file.path(Sys.getenv("DFP_DATA"), "output_V6"), mustWork = FALSE)) &&
    !nzchar(Sys.getenv("PIPE_ALLOW_PUBLISHED"))) {
  stop("Refusing to write into the PUBLISHED outputs ($DFP_DATA/output_V6). ",
       "Set PIPE_OUT_DIR elsewhere, or PIPE_ALLOW_PUBLISHED=1 to force.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Check-pointing system for fast resumption
chk_dir <- file.path(output_dir, "checkpoints")
dir.create(chk_dir, showWarnings = FALSE, recursive = TRUE)

# Remove the merge checkpoint to force a clean re-merge
unlink(file.path(chk_dir, "stage3A_fusion.rds"))

# Remove the GMM and DBSCAN checkpoints to force their recomputation
gmm_checkpoint <- file.path(chk_dir, "stage3C_gmm.rds")
dbscan_checkpoint <- file.path(chk_dir, "stage3B_dbscan.rds")

if (file.exists(gmm_checkpoint)) {
  cat("🗑️  Removing the faulty GMM checkpoint to force recomputation\n")
  unlink(gmm_checkpoint)
}

if (file.exists(dbscan_checkpoint)) {
  cat("🗑️  Removing the faulty DBSCAN checkpoint to force recomputation\n")
  unlink(dbscan_checkpoint)
}

checkpoint <- function(file, expr, force_recompute = FALSE) {
  file <- file.path(chk_dir, file)
  if (!force_recompute && file.exists(file)) {
    cat("🔄  Reload checkpoint:", basename(file), "\n")
    obj <- readRDS(file)
    gc()
    obj
  } else {
    cat("⚙️  Compute & save checkpoint:", basename(file), "\n")
    obj <- force(expr)
    saveRDS(obj, file, compress = "xz")
    cat("✅ Checkpoint saved:", basename(file), "\n")
    gc()
    obj
  }
}

if (results_dir != "" && dir.exists(results_dir)) {
  expected_dir <- results_dir
  cat("✅ Using the specified folder:", expected_dir, "\n")
} else {
  cat("🔍 Auto-detecting the results folder...\n")
  
  possible_paths <- c(
    file.path(Sys.getenv("DFP_DATA"), paste0("ais_results_", split_job_id)),
    file.path(Sys.getenv("DFP_DATA"), paste0("ais_split_", split_job_id)),
    file.path(Sys.getenv("DFP_DATA"), split_job_id)
  )
  
  expected_dir <- NULL
  for (path in possible_paths) {
    if (dir.exists(path)) {
      expected_dir <- path
      cat("✅ Folder found:", expected_dir, "\n")
      break
    } else {
      cat("   ❌ Not found:", path, "\n")
    }
  }
  
  if (is.null(expected_dir)) {
    stop("❌ ERROR: No results folder found")
  }
}

clean_files <- list.files(expected_dir, pattern = ".*_clean\\.rds$", full.names = TRUE)
if (length(clean_files) == 0) {
  stop("❌ ERROR: No *_clean.rds file found in ", expected_dir)
}

cat("✅ Input files checked:", length(clean_files), "*_clean.rds files found\n\n")

cat("\n🏁  STEP-3  |  Merge & Modelling  |  start:", format(Sys.time()), "\n\n")

# CRITICAL R CONFIGURATION - BEFORE LOADING PACKAGES
# Add the paths without overwriting the bash configuration
user_libs <- c("~/R/library")
.libPaths(unique(c(user_libs, .libPaths())))
cat("✅ R will look for packages in:", paste(.libPaths(), collapse = " | "), "\n")

# Check that all required packages are available
cat("\n🔍 DIAGNOSTIC OF REQUIRED PACKAGES:\n")
needed <- c("glmnet", "doParallel", "pROC", "data.table", "dbscan", "mclust", "yaml", "geosphere", "lubridate", "zoo", "solitude")
for (p in needed) {
  status <- if(requireNamespace(p, quietly=TRUE)) "✅ OK" else "❌ MISSING"
  cat(sprintf("   %-12s: %s\n", p, status))
}
cat("\n")

# --------------------------------------------------------------------------
# 1. ── SLURM / ENV PARAMETERS --------------------------------------------
# --------------------------------------------------------------------------
split_id    <- split_job_id  # same ID as Step-1

split_dir <- expected_dir

clean_files_check <- list.files(split_dir, pattern = "_clean\\.rds$", full.names = TRUE)
if (!length(clean_files_check)) {
  stop("❌ No *_clean.rds file found in ", split_dir, "\n",
       "   Files present in the directory:\n",
       paste("   -", list.files(split_dir, pattern = "\\.rds$"), collapse = "\n"),
       "\n   Check that step 2 finished correctly.")
}

cat("✅ *_clean.rds files found:", length(clean_files_check), "\n")
for (f in head(clean_files_check, 3)) cat("   •", basename(f), "\n")
if (length(clean_files_check) > 3) cat("   ... and", length(clean_files_check) - 3, "others\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("💾  Source:", split_dir, "\n")
cat("💾  Output:", output_dir, "\n")
cat("💾  Checkpoints:", chk_dir, "\n\n")

# --------------------------------------------------------------------------
# 2. ── LOADING PACKAGES  +  CPU CONFIGURATION -----------------------------
# --------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(solitude)    # already used but required to reload config
  library(zoo)         # rolling
  library(dbscan)
  library(mclust)
  library(pROC)
  library(parallel)
  library(yaml)
  library(geosphere)
})

setDTthreads(cores)
options(datatable.optimize = 3)
options(mc.cores = cores)
Sys.setenv(OMP_NUM_THREADS      = cores,
           OPENBLAS_NUM_THREADS = cores)

cat("✅  Packages loaded |", cores, "threads data.table\n")

# --------------------------------------------------------------------------
# STRICT FILTER CONTROL FUNCTION - RELAXED IN DRY-RUN
# --------------------------------------------------------------------------
MIN_LEFT_DRY <- 1000
MIN_LEFT_PROD <- 20000

check_filter <- function(step_name, n_before, n_after, 
                        max_frac = if (is_dry) 0.90 else 0.25, 
                        min_left = if (is_dry) MIN_LEFT_DRY else MIN_LEFT_PROD, 
                        fatal = TRUE) {
  n_removed <- n_before - n_after
  frac_removed <- if (n_before == 0) 0 else n_removed / n_before
  
  cat(sprintf("[DIAGNOSTIC] %s: %d points removed (%.2f%%), %d remaining\n", 
              step_name, n_removed, 100*frac_removed, n_after))
  
  if (fatal) {
    if (frac_removed > max_frac) {
      stop(sprintf("FATAL ERROR: %s removed %.1f%% of the points (>%.0f%%). Check your parameters or data!", 
                   step_name, 100*frac_removed, 100*max_frac))
    }
    if (n_after < min_left) {
      stop(sprintf("FATAL ERROR: Only %d points remain after '%s' (<%d), pipeline stopped.", 
                   n_after, step_name, min_left))
    }
  } else {
    if (frac_removed > max_frac) {
      warning(sprintf("[ALERT] %s removed %.1f%% of the points (>%.0f%%)!", 
                   step_name, 100*frac_removed, 100*max_frac))
    }
    if (n_after < min_left) {
      warning(sprintf("[ALERT] Only %d points remain after '%s' (<%d)!", 
                   n_after, step_name, min_left))
    }
  }
}

# --------------------------------------------------------------------------
# 3. ── CONFIGURATION  ------------------------------------------------------
# --------------------------------------------------------------------------
cat("✅  YAML loaded\n")

# Default values if the YAML is incomplete
default_config <- list(
  contamination_rate = 0.02,
  if_num_trees = 25,
  if_sample_size = 256,
  spike_factor = 5,
  min_spike_duration = 2,
  percentile_lower = 0.05,
  percentile_upper = 0.95,
  min_context_points = 200,
  max_accel_ms2 = 0.3,
  max_turn_at_speed = 30
  # The dredging speed applied is the 1–3.5 kn window
  # (dredge_range, dredge_max, speed_bin, spd_max).
)

tryCatch({
  cfg <- yaml::read_yaml("~/R_scripts/configuration/outlier_config_V6.yaml")
  
  if (!is.null(cfg$isolation_forest)) {
    cfg$contamination_rate <- cfg$isolation_forest$contamination_rate
    cfg$if_num_trees <- cfg$isolation_forest$num_trees
    cfg$if_sample_size <- cfg$isolation_forest$sample_size
  }
  
  if (!is.null(cfg$temporal_spikes)) {
    cfg$spike_factor <- cfg$temporal_spikes$spike_factor
    cfg$min_spike_duration <- cfg$temporal_spikes$min_duration
  }
  
  if (!is.null(cfg$contextual)) {
    cfg$percentile_lower <- cfg$contextual$percentile_lower
    cfg$percentile_upper <- cfg$contextual$percentile_upper
    cfg$min_context_points <- cfg$contextual$min_context_points
  }
  
  if (!is.null(cfg$physics)) {
    cfg$max_accel_ms2 <- cfg$physics$max_acceleration_ms2
    cfg$max_turn_at_speed <- cfg$physics$max_turn_at_speed
  }
  
  for (param in names(default_config)) {
    if (is.null(cfg[[param]])) {
      cfg[[param]] <- default_config[[param]]
      cat("⚠️  Parameter missing in YAML:", param, "- Using default value:", default_config[[param]], "\n")
    }
  }
}, error = function(e) {
  cat("⚠️  YAML reading error - Using default values\n")
  cfg <<- default_config
})

# --------------------------------------------------------------------------
# 4. ── DATA MERGE  ---------------------------------------------------------
# --------------------------------------------------------------------------
# CHECK-POINT 1: Data merge
ais <- checkpoint("stage3A_fusion.rds", {
cat("🔍  Files found:", length(clean_files), "\n")
for (f in head(clean_files, 5)) cat("   •", basename(f), "\n")
if (length(clean_files) > 5) cat("   …\n")

t_read <- system.time({
  # parallel reading (limit 4 workers → IO-safe enough)
  ais_list <- mclapply(clean_files, readRDS, mc.cores = min(cores, 4))
  ais      <- rbindlist(ais_list, fill = TRUE)
})
rm(ais_list); gc()
cat(sprintf("✅  Merge: %s rows | %.1f s\n\n",
            format(nrow(ais), big.mark = " "), t_read[3]))
  
  ais
})

n_total <- nrow(ais)

# CRITICAL CHECK: immediate stop if rows without Annee are present after the merge
n_na_year <- ais[is.na(Annee), .N]
if (n_na_year > 0) {
  stop(n_na_year, " rows without Annee after merge – checkpoint not saved.")
}

# --------------------------------------------------------------------------
# CONFIGURABLE DRY-RUN: DRY_RUN enables, DRYRUN_N sets the size
# --------------------------------------------------------------------------
DRYRUN_N <- as.integer(Sys.getenv("DRYRUN_N", "0"))
if (is_dry && DRYRUN_N > 0 && nrow(ais) > DRYRUN_N) {
  set.seed(42)  # reproducible
  
  # Sampling stratified by vessel to keep representativeness
  ais <- ais[, .SD[sample(.N, min(.N, ceiling(DRYRUN_N * .N / n_total)))], by = Navire]
  
  cat(sprintf("🧪 Dry-run: subset of %d rows (%.1f %% of the initial total)\n",
              nrow(ais), 100 * nrow(ais) / n_total))
  
  n_total <- nrow(ais)
}

# --------------------------------------------------------------------------
# 4 bis. ── CHECK AND CORRECTION OF COLUMN TYPES ---------------------------
# --------------------------------------------------------------------------
cat("🔧  Checking column types...\n")

cat(" Types of the critical columns:\n")
for (col in c("Navire", "ssvid", "Seg_id", "Annee")) {
  if (col %in% names(ais)) {
    cat(sprintf("   %-10s: %s\n", col, class(ais[[col]])[1]))
  }
}

# Harmonise types to avoid join errors
if ("Navire" %in% names(ais)) {
  if (!is.character(ais$Navire)) {
    cat("⚠️  Converting Navire to character\n")
    ais[, Navire := as.character(Navire)]
  }
}

if ("ssvid" %in% names(ais)) {
  if (!is.integer(ais$ssvid) && !is.numeric(ais$ssvid)) {
    cat("⚠️  Converting ssvid to integer\n")
    ais[, ssvid := as.integer(ssvid)]
  }
}

if ("Seg_id" %in% names(ais)) {
  if (!is.character(ais$Seg_id)) {
    cat("⚠️  Converting Seg_id to character\n")
    ais[, Seg_id := as.character(Seg_id)]
  }
}

if ("Annee" %in% names(ais)) {
  if (!is.integer(ais$Annee) && !is.numeric(ais$Annee)) {
    cat("⚠️  Converting Annee to integer\n")
    ais[, Annee := as.integer(Annee)]
  }
}

# ADDITIONAL CRITICAL CHECKS
cat("🔍 Additional critical checks...\n")

# 1. Check of the is_stop column (used in DBSCAN)
if (!"is_stop" %chin% names(ais)) {
  stop("❌ Column 'is_stop' missing - Check that Step-2 finished correctly")
} else {
  cat("✅ Column 'is_stop' present\n")
}

# 2. Check of the geographic columns (used in DBSCAN)
geo_cols <- c("Lon", "Lat", "Course", "Speed")
missing_geo <- geo_cols[!geo_cols %chin% names(ais)]
if (length(missing_geo) > 0) {
  stop("❌ Missing geographic columns: ", paste(missing_geo, collapse = ", "))
} else {
  cat("✅ Geographic columns present\n")
  ais[, c("Lon", "Lat", "Course", "Speed") := lapply(.SD, as.numeric), .SDcols = c("Lon", "Lat", "Course", "Speed")]
  cat("✅ Geographic columns converted to numeric\n")
}

# 3. Check for duplicates (Seg_id, Timestamp)
if (anyDuplicated(ais, by = c("Seg_id", "Timestamp"))) {
  warning("⚠️  Duplicates detected in (Seg_id, Timestamp) - Removing...")
  ais <- unique(ais, by = c("Seg_id", "Timestamp"))
  cat("✅ Duplicates removed\n")
} else {
  cat("✅ No duplicate detected\n")
}

# 4. Check of the available memory
mem_usage <- object.size(ais) / 1024^3  # GB
cat(sprintf("📊 Memory use: %.2f GB\n", mem_usage))

if (mem_usage > 50) {
  warning("⚠️  High memory use: ", round(mem_usage, 1), " GB")
}

cat("✅ Column types checked and corrected\n\n")

# --------------------------------------------------------------------------
# STRICT CHECK OF COUNTS - STOP IF THERE IS A PROBLEM
# --------------------------------------------------------------------------
n_total <- nrow(ais)  # Initial count
cat("📊 === STRICT CHECK OF COUNTS ===\n")
cat(sprintf("📈 Number of points at start: %d\n", n_total))

# --------------------------------------------------------------------------
# 5. ── REMOVE outliers_IF DIRECTLY -----------------------------------------
# --------------------------------------------------------------------------
n_before_outlier <- nrow(ais)
if ("outlier_IF" %chin% names(ais)) {
  ais <- ais[outlier_IF == FALSE | is.na(outlier_IF)]
  check_filter("Isolation Forest", n_before_outlier, nrow(ais))
} else {
  cat("[DIAGNOSTIC] Isolation Forest: No filter applied (column absent)\n")
}

# --------------------------------------------------------------------------
# 6. ── UTILITY FUNCTIONS -----------------------------------------------------
# --------------------------------------------------------------------------
roll_MAD5 <- function(x) {                # windowed medians & MAD 5
  med <- zoo::rollapplyr(x, 5, median, fill = NA, align = "center")
  q25 <- zoo::rollapplyr(x, 5, quantile, probs = .25, fill = NA, align = "center")
  q75 <- zoo::rollapplyr(x, 5, quantile, probs = .75, fill = NA, align = "center")
  mad <- 1.4826 * (q75 - q25)
  list(med = med, mad = mad)
}

detect_spikes <- function(speed, spike_factor = 5, min_iso = 2) {
  r   <- roll_MAD5(speed)
  dev <- abs(speed - r$med)
  cand <- !is.na(r$mad) & r$mad > 0 & dev > spike_factor * r$mad

  for (idx in which(cand)) {
    rng <- max(1, idx - min_iso):min(length(speed), idx + min_iso)
    if (sum(abs(speed[rng] - speed[idx]) < r$mad[idx], na.rm = TRUE) > min_iso)
      cand[idx] <- FALSE
  }
  cand
}

# DBSCAN helpers: projection in metres + kNNdist (avoids O(n²))
to_xy_m <- function(lon, lat) {
  # Simple projection: coordinates (lon, lat) → (x, y) in metres
  lat0 <- mean(lat, na.rm = TRUE) * pi / 180
  x <- lon * 111320 * cos(lat0)
  y <- lat * 110574
  cbind(x, y)
}

auto_eps_m <- function(xy, k = 4) {
  # Computation of eps in metres via kNNdist (O(n log n) instead of O(n²))
  # k_eff avoids a crash if n <= k
  n <- nrow(xy)
  k_eff <- min(k, max(1, n - 1))
  d <- dbscan::kNNdist(xy, k = k_eff)
  # Minimum threshold of 200m, otherwise 1% of the median
  pmax(200, 0.01 * median(d, na.rm = TRUE))
}

# Behaviour transition smoothing (same as the Step-2 script)
smooth_transition <- function(lbl, spd, cost, thr = 3,
                              dredge_min = 1, dredge_max = 3.5) {
  rle_lbl <- rle(lbl); cum <- cumsum(rle_lbl$lengths); beg <- c(1, head(cum, -1) + 1)
  for (i in seq_along(rle_lbl$lengths)) {
    if (rle_lbl$lengths[i] < thr) {
      vbar <- mean(spd[beg[i]:cum[i]])
      if (vbar <= dredge_max) {
        prev <- if (i > 1) rle_lbl$values[i - 1] else NA
        next_val <- if (i < length(rle_lbl$values)) rle_lbl$values[i + 1] else NA
        c_prev <- if (!is.na(prev)) cost[prev, rle_lbl$values[i]] else Inf
        c_next <- if (!is.na(next_val)) cost[next_val, rle_lbl$values[i]] else Inf
        rle_lbl$values[i] <- if (c_prev < c_next) prev else next_val
      }
    }
  }
  inverse.rle(rle_lbl)
}

# --------------------------------------------------------------------------
# 7. ── TEMPORAL SPIKES  (Step-4 of the pipeline) --------------------------
# --------------------------------------------------------------------------
cat("⏱️  Detecting isolated temporal spikes…\n")
n_before_spike <- nrow(ais)
ais[, outlier_spike := detect_spikes(Speed,
                        spike_factor = cfg$spike_factor,
                        min_iso      = cfg$min_spike_duration),
    by = Seg_id]
ais <- ais[outlier_spike == FALSE | is.na(outlier_spike)]
check_filter("Spike filter", n_before_spike, nrow(ais))
ais[, outlier_spike := NULL]

# --------------------------------------------------------------------------
# 8. ── CONTEXTUAL PERCENTILES: not computed (unused downstream) -----------
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 9. ── DETECTION OF SPATIAL STOPS (DBSCAN) ---------------------------------
# --------------------------------------------------------------------------
# CHECK-POINT 2: DBSCAN finished
ais <- checkpoint("stage3B_dbscan.rds", {
cat("🔍  Detecting spatial stops with DBSCAN...\n")

# Protection against small vessels
ais[, stop_cluster := NA_integer_]

# Count with Lon AND Lat filter (avoids skipping a whole vessel for 1 NA)
stop_counts <- ais[is_stop == TRUE & !is.na(Lon) & !is.na(Lat), .N, by = .(Navire, Annee)]
stop_counts <- stop_counts[N >= 4]  # Only vessels with enough points

if (nrow(stop_counts) > 0) {
  cat("📊 DBSCAN processing for", nrow(stop_counts), "vessel/year combinations\n")

  for (i in 1:nrow(stop_counts)) {
    vessel <- stop_counts$Navire[i]
    annee <- stop_counts$Annee[i]

    # Selection via .I for a guaranteed alignment
    idx <- ais[Navire == vessel & Annee == annee & is_stop == TRUE &
               !is.na(Lon) & !is.na(Lat), which = TRUE]
    arr <- ais[idx]

    if (length(idx) >= 4) {
      # Check of the coordinate variance
      if (var(arr$Lon, na.rm = TRUE) < 1e-8 || var(arr$Lat, na.rm = TRUE) < 1e-8) {
        cat("⚠️  Insufficient variance for", vessel, annee, "- Skipping\n")
        next
      }

      tryCatch({
        # Projection in metres (metric consistency + stable eps)
        xy <- to_xy_m(arr$Lon, arr$Lat)

        # Protection against large stops (sampling for eps)
        if (nrow(xy) > 5000) {
          cat("⚠️  Large stop detected for", vessel, annee, "(", nrow(xy), "points) - Sampling for eps\n")
          set.seed(123)
          samp_idx <- sample(nrow(xy), 5000)
          eps <- auto_eps_m(xy[samp_idx, , drop = FALSE])
        } else {
          eps <- auto_eps_m(xy)
        }

        # DBSCAN in metric coordinates (xy in metres)
        cl <- dbscan(xy, eps = eps, minPts = 4)$cluster

        # Protection against empty clusters
        if (sum(cl > 0) == 0) {
          cat("⚠️  No cluster found for", vessel, annee, "- Widening the eps radius\n")
          eps <- eps * 1.5
          cl <- dbscan(xy, eps = eps, minPts = 4)$cluster
        }

        # Assignment via idx (guaranteed alignment)
        ais[idx, stop_cluster := cl]
      }, error = function(e) {
        cat("⚠️  DBSCAN error for", vessel, annee, ":", e$message, "\n")
      })
    }
  }
} else {
  cat("⚠️  No vessel with enough stop points for DBSCAN\n")
}

ais[, is_stop_spatial := is_stop & !is.na(stop_cluster) & stop_cluster > 0]
  
  cat("✅ DBSCAN finished - Checkpoint saved\n")
  
  # Counts after DBSCAN
  n_stop_spatial <- sum(ais$is_stop_spatial, na.rm = TRUE)
  cat(sprintf("[DIAGNOSTIC] After DBSCAN (is_stop_spatial TRUE): %d points (%.2f%% of the total)\n", 
              n_stop_spatial, 100 * n_stop_spatial / nrow(ais)))
  
  if (n_stop_spatial < 1000) {
    warning("[ALERT] Very few spatial stop points detected by DBSCAN (< 1000)")
  }
  
  if (!"is_stop_spatial" %in% names(ais)) {
    stop("❌ CRITICAL ERROR: Column 'is_stop_spatial' missing after DBSCAN")
  }
  cat("✅ Check: Column 'is_stop_spatial' present\n")
  
  ais
})

# --------------------------------------------------------------------------
# 10 bis. ── AIS cadence & adaptive thresholds  ----------------------------
# --------------------------------------------------------------------------

# 1) Sort by segment and timestamp to avoid negative dt_sec
setorder(ais, Seg_id, Timestamp)

# 2) Δt between messages, in seconds
ais[, dt_sec := c(NA_real_, diff(as.numeric(Timestamp))), by = Seg_id]

# Course_change only (Accel is computed after filtering)
# by = Seg_id avoids inter-segment diffs
ais[, Course_change := c(NA, abs(diff(Course))), by = Seg_id]
# Handling of the wrap-around for Course_change
ais[Course_change > 180, Course_change := 360 - Course_change]
# Robust fallback if too many NA
ais[is.na(Course_change), Course_change := 0]

# 3) Computation of the cadence statistics per vessel
cadence <- ais[!is.na(dt_sec), .(
    p50_dt = median(dt_sec, na.rm = TRUE),      # median
    p95_dt = quantile(dt_sec, 0.95, na.rm = TRUE),  # 95th percentile
    n_obs = .N
), by = Navire]

# 4) Truncation of the median to avoid extreme cases
cadence[, p50_dt := pmin(p50_dt, 600)]  # the median is truncated at 10 min

# 5) Computation of the adaptive thresholds
#    - Base threshold: 3 × median (captures the real pauses)
#    - Floor: 180 s (3 min) for very active vessels
#    - Ceiling: 900 s (15 min) for normal vessels
#    - Special threshold for very slow vessels (> 50% of intervals > 5 min)

cadence[, adaptive_threshold := pmin(900, pmax(180, 3 * p50_dt))]

# 6) Detection of very slow vessels (more than 50% of intervals > 5 min)
intervals_long <- ais[dt_sec >= 300, .N, by = Navire]
intervals_tot <- ais[!is.na(dt_sec), .N, by = Navire]
prop_long <- merge(intervals_long, intervals_tot, by = "Navire", suffixes = c("_long", "_total"))
prop_long[, prop_long := N_long / N_total]

# 7) Adjustment of the threshold for very slow vessels
cadence <- merge(cadence, prop_long[, .(Navire, prop_long)], by = "Navire", all.x = TRUE)
cadence[is.na(prop_long), prop_long := 0]

# For vessels with > 50% long intervals, use a higher threshold
cadence[prop_long > 0.5, adaptive_threshold := pmin(1200, 5 * p50_dt)]  # 5 × median, max 20 min

# 8) Computation of the adaptive n_min for run-length smoothing
#    t_threshold = 30 s ≈ max duration of a micro-glitch
t_threshold   <- 30                     # 30 s ≈ max duration of a micro-glitch
n_max_cap <- 8                      # ceiling to stay fast
cadence[, n_min := pmin(n_max_cap, pmax(3, ceiling(p50_dt / t_threshold)))]

ais <- merge(ais, cadence[, .(Navire, n_min, adaptive_threshold, p50_dt, prop_long)], by = "Navire", all.x = TRUE)

# 10) Filtering with the adaptive threshold
cat("🔧  AIS cadence & adaptive thresholds per vessel\n")
print(cadence[order(p50_dt), .(Navire, p50_dt, adaptive_threshold, n_min, prop_long)])

# Filtering of intervals that are too long according to each vessel's adaptive threshold
n_before_cadence <- nrow(ais)
ais <- ais[(dt_sec > 0 & dt_sec < adaptive_threshold) | is.na(dt_sec)]
check_filter("Long-interval filtering", n_before_cadence, nrow(ais))

# Recompute dt_sec AFTER filtering to avoid the time shift
ais[, dt_sec := c(NA_real_, diff(as.numeric(Timestamp))), by = Seg_id]

# Recompute Course_change AFTER filtering to avoid the bias
ais[, Course_change := c(NA, abs(diff(Course))), by = Seg_id]
# Handling of the wrap-around for Course_change
ais[Course_change > 180, Course_change := 360 - Course_change]
# Robust fallback if too many NA
ais[is.na(Course_change), Course_change := 0]

# Accel computed AFTER filtering the large time gaps
# Conversion of Speed to m/s to compute Accel in m/s²
# Safeguard against empty vectors
ais[, Accel := {
  dt <- dt_sec[-1L]
  if (length(dt) == 0) {
    rep(NA_real_, .N)
  } else {
    c(NA, diff(Speed * 0.514444) / pmax(dt, 1))
  }
}, by = Seg_id]
ais[is.na(Accel), Accel := 0]

# -- n_min for NA cases ---------------------------------------------------
ais[is.na(n_min), n_min := 3]  # floor fallback

# FINAL CHECK BEFORE GMM
cat(sprintf("[DIAGNOSTIC] Number of points kept before GMM: %d (%.2f%% of the initial total)\n", 
            nrow(ais), 100 * nrow(ais) / n_total))

# Critical final check - RELAXED IN DRY-RUN
if (is_dry) {
  # In dry-run, a smaller dataset is tolerated
  if (nrow(ais) < 5000) {
    stop(sprintf("FATAL ERROR: Only %d points remain before GMM (< 5000), sample too small for the dry-run.", nrow(ais)))
  }
} else {
  if (nrow(ais) < 20000) {
    stop(sprintf("FATAL ERROR: Only %d points remain before GMM (< 20000), pipeline stopped.", nrow(ais)))
  }
}

cat("\n[DIAGNOSTIC] Distribution of the dataset just before the GMM:\n")
print(table(ais$Navire))
cat("Total remaining points:", nrow(ais), "\n")

# --------------------------------------------------------------------------
# 11. ── GMM  "stop + K mobile"  (K chosen by BIC between 3 and 5) - ROBUST
# --------------------------------------------------------------------------
ais <- checkpoint("stage3C_gmm.rds", {

  # Retain fitted parameters and exact draws for supplementary diagnostics.
  # The helper only observes this run; it must not consume random numbers.
  source(Sys.getenv("GMM_DIAGNOSTICS_HELPER", "gmm_diagnostics.R"))
  gd_states <- list(entry=gd_rng())
  gd_candidates <- list()

  cat("🔧  GMM 'stop + K mobile' – selection of K by BIC (3-5) - ROBUST\n")

  ## 0. Preparation ----------------------------------------------------------
  idx_mobile <- which(!ais$is_stop_spatial & ais$Speed >= 0.5 & ais$Speed <= 20)

  # Quantiles recomputed after all filters to avoid the bias
  q95_acc  <- quantile(abs(ais$Accel),  .95, na.rm = TRUE);  if (!is.finite(q95_acc)  || q95_acc  == 0) q95_acc  <- 1e-6
  q95_turn <- quantile(ais$Course_change, .95, na.rm = TRUE); if (!is.finite(q95_turn) || q95_turn == 0) q95_turn <- 1e-6
  
  cat("   📊 Quantiles recomputed after filters - Accel 95%:", round(q95_acc, 4), "| Course 95%:", round(q95_turn, 2), "°\n")

  mobile <- ais[idx_mobile, .(
      Speed,
      norm_acc    = pmax(0, pmin(1, abs(Accel) / q95_acc)),
      norm_course = pmax(0, pmin(1, 1 - Course_change / q95_turn))
  )]

  keep       <- complete.cases(mobile)
  X          <- as.matrix(mobile[keep])
  idx_train  <- idx_mobile[keep]

  # Adaptive threshold according to the mode (dry-run vs production)
  min_mobile <- if (is_dry) 50 else 200
  if (nrow(X) < min_mobile)
    stop("❌  Not enough complete mobile points to fit the GMM (", nrow(X), " < ", min_mobile, ")")

  ## 1. Search for the best K (BIC) - ROBUST --------------------------------
  # Protection against Mclust memory overflow
  gd_states$before_selection_sample <- gd_rng()
  gd_bic_rows <- if (nrow(X) > 1e6) sample(nrow(X), 1e6) else seq_len(nrow(X))
  samp <- X[gd_bic_rows, , drop=FALSE]
  cat("   Sample for BIC:", nrow(samp), "points (out of", nrow(X), ")\n")
  
  G_choices <- 3:5                                   # 3 to 5 mobile components
  
  # Helper to extract a scalar BIC
  get_bic_scalar <- function(x) {
    if (is.null(x) || length(x) == 0) return(-Inf)
    max(as.numeric(x), na.rm = TRUE)
  }
  
  # vapply() forces a numeric scalar from the start
  bic_vals <- vapply(G_choices, function(g) {
    gd_before <- gd_rng()
    gd_result <- gd_capture(Mclust(samp, G = g, modelNames = "VVV", verbose = FALSE))
    res <- gd_result$fit
    gd_candidates[[as.character(g)]] <<- list(model=gd_compact(res),
      warnings=gd_result$warnings,error=gd_result$error,rng_before=gd_before,rng_after=gd_rng())
    if (is.null(res)) return(-Inf)
    get_bic_scalar(res$bic)
  }, numeric(1))
  
  # Protection against Mclust failures
  if (all(is.infinite(bic_vals)))
    stop("❌  Mclust failed for all values of K - X too large or ill-conditioned")
  
  K <- G_choices[which.max(bic_vals)]
  
  # Check of the lower bound K >= 3
  if (K < 3) stop("❌ The GMM must have at least 3 mobile components (dredge/transit)")
  
  cat("   ➜  K retained by BIC:", K, "mobile components\n")

  ## 2. Final fit ------------------------------------------------------------
  # Memory protection for the final fit
  gd_states$before_final_sample <- gd_rng()
  gd_fit_rows <- if (nrow(X) > 2e6) sample(nrow(X), 2e6) else seq_len(nrow(X))
  X_fit <- X[gd_fit_rows, , drop=FALSE]
  cat("   Sample for the final fit:", nrow(X_fit), "points (out of", nrow(X), ")\n")
  
  gd_states$before_final_fit <- gd_rng()
  gd_final_result <- gd_capture(Mclust(X_fit, G = K, modelNames = "VVV", verbose = FALSE))
  gmm <- gd_final_result$fit
  if(is.null(gmm)) stop("Final GMM failed: ",gd_final_result$error)
  gd_states$after_final_fit <- gd_rng()

  ## 3. Posterior probabilities ---------------------------------------------
  z <- predict(gmm, X)$z                       # nrow(X) × K
  
  # Safe re-initialisation of the p columns
  ais[, paste0("p", 1:(K+1)) := 0]             # creates p1…p{K+1} and resets to zero
  ais[is_stop_spatial == TRUE, p1 := 1]

  for (j in seq_len(K))
    ais[idx_train, paste0("p", j + 1) := z[, j]]

  ## 4. Reordering by speed --------------------------------------------------
  ord  <- order(gmm$parameters$mean[1, ])      # increasing order of Speed
  mu   <- gmm$parameters$mean[1, ord]
  cat("   μ of mobile components sorted:", paste(round(mu, 2), collapse = " / "), "kn\n")

  Zord <- z[, ord, drop = FALSE]

  # Full re-initialisation before re-insertion
  ais[, paste0("p", 1:(K+1)) := 0]             # resets everything to zero
  ais[is_stop_spatial == TRUE, p1 := 1]        # puts the stops back
  for (j in seq_len(K))
    ais[idx_train, paste0("p", j + 1) := Zord[, j]]

  ## 5. Raw labels - OPTION A: Robust mapping of multiple slow components
  # Correct alignment of mobile_keep / Zord for the stability score
  # Separation idle (μ<1) / slow (1≤μ<6) / transit (μ≥6)
  # Guarantees that NO transit has μ < 6 kn

  ## --- Align mobile with Zord (mobile_keep) + fast max.col() --- ##
  mobile_keep <- mobile[keep]  # aligns the rows exactly with X / Zord
  comp_hat  <- max.col(Zord, ties.method = "first")  # 1..K, consistent with mu / ord

  ## --- Robust mapping of the components (Option A) --------------------- ##
  dredge_range <- c(1, 3.5)           # window (kn) considered as dredging
  slow_ceiling <- 6                    # μ < 6 kn = never transit

  # 1) Identify the components by speed range
  idle_candidates   <- which(mu < 1)                                    # μ < 1 kn → other

  slow_candidates   <- which(mu >= 1 & mu < slow_ceiling)               # μ ∈ [1, 6) kn
  dredge_candidates <- which(mu >= dredge_range[1] & mu <= dredge_range[2])  # μ ∈ [1, 3.5] kn

  cat("   📊 Idle components (μ < 1 kn):", paste(idle_candidates, collapse=", "), "\n")
  cat("   📊 Slow components (1 ≤ μ < 6 kn):", paste(slow_candidates, collapse=", "), "\n")
  cat("   📊 Dredging candidates [1-3.5 kn]:", paste(dredge_candidates, collapse=", "), "\n")

  # 2) Robust selection of the dredging component (mobile_keep + comp_hat)
  if (length(dredge_candidates) > 0) {
    if (length(dredge_candidates) == 1) {
      dredge_comp <- dredge_candidates
    } else {
      # Several candidates: stability score = mean(norm_course) - mean(norm_acc)
      # Dredging has a stable course (high norm_course) and low acceleration
      stability_scores <- vapply(dredge_candidates, function(j) {
        m <- (comp_hat == j)
        if (sum(m) < 10) return(-Inf)
        mean(mobile_keep$norm_course[m], na.rm = TRUE) -
          mean(mobile_keep$norm_acc[m], na.rm = TRUE)
      }, numeric(1))

      # Selection by stability score, tie-breaker = closest to 2 kn
      best_score <- max(stability_scores, na.rm = TRUE)
      top_candidates <- dredge_candidates[stability_scores >= best_score - 0.05]
      dredge_comp <- top_candidates[which.min(abs(mu[top_candidates] - 2))]

      cat("   📊 Stability scores:", paste(round(stability_scores, 3), collapse=", "), "\n")
    }
  } else {
    # Fallback: closest to 2 kn among all components
    dredge_comp <- which.min(abs(mu - 2))
    cat("   ⚠️  No candidate in [1-3.5 kn], fallback to μ =", round(mu[dredge_comp], 2), "kn\n")
  }

  # 3) The other slow components (other than dredging) → slow_maneuvers
  other_slow_idx <- setdiff(slow_candidates, dredge_comp)

  # 4) Transit components: ONLY μ >= 6 kn
  transit_pool <- which(mu >= slow_ceiling)

  cat("   📊 Pool transit (μ >= 6 kn):", paste(transit_pool, collapse=", "), "\n")

  # 5) Build the labels
  mobile_labels <- character(K)
  mobile_labels[dredge_comp] <- "dredging"

  # idle_candidates (μ < 1 kn) → "other"
  if (length(idle_candidates) > 0) {
    mobile_labels[idle_candidates] <- "other"
  }

  # slow_maneuvers for the slow components other than dredging
  if (length(other_slow_idx) > 0) {
    mobile_labels[other_slow_idx] <- "slow_maneuvers"
  }

  if (length(transit_pool) >= 2) {
    # 2+ transits: loaded (slower) and unloaded (faster)
    loaded_idx <- transit_pool[1]
    unloaded_idx <- transit_pool[length(transit_pool)]
    mobile_labels[loaded_idx] <- "loaded_transit"
    mobile_labels[unloaded_idx] <- "unloaded_transit"
    # Intermediate ones = other
    middle_transit <- setdiff(transit_pool, c(loaded_idx, unloaded_idx))
    if (length(middle_transit) > 0) mobile_labels[middle_transit] <- "other"
  } else if (length(transit_pool) == 1) {
    # Single transit: unloaded_transit
    mobile_labels[transit_pool] <- "unloaded_transit"
  }
  # else: no transit (all components are slow)

  # 6) Everything left = "other"
  mobile_labels[mobile_labels == ""] <- "other"

  ## --- Post-mapping assertions ------------------------------------------ ##
  # ASSERTION 0: dredge_comp must be unique and defined
  if (length(dredge_comp) != 1) {
    stop(sprintf("MAPPING ERROR: dredging not unique or missing (dredge_comp = %s)",
                 paste(dredge_comp, collapse=", ")))
  }

  # ASSERTION 1: no transit with μ < 6 kn
  transit_labels <- c("loaded_transit", "unloaded_transit")
  for (lbl in transit_labels) {
    idx <- which(mobile_labels == lbl)
    if (length(idx) > 0 && any(mu[idx] < slow_ceiling)) {
      stop(sprintf("MAPPING ERROR: %s assigned to μ = %.2f kn (< %.1f kn)\n  mu = %s\n  labels = %s",
                   lbl, mu[idx[mu[idx] < slow_ceiling]], slow_ceiling,
                   paste(round(mu, 2), collapse=" / "),
                   paste(mobile_labels, collapse=" → ")))
    }
  }

  # ASSERTION 2: systematic diagnostic dump
  cat("\n   === MAPPING DIAGNOSTIC ===\n")
  mapping_dt <- data.table(comp = 1:K, mu_kn = round(mu, 2), label = mobile_labels)
  print(mapping_dt)
  cat("   ===========================\n\n")

  map <- c("stops", mobile_labels)             # length = K+1 (p1…p{K+1})

  cat("   🏷️  Assigned labels:", paste(mobile_labels, collapse=" → "), "\n")
  cat("   📊 Corresponding speeds:", paste(round(mu, 2), "kn", collapse=" → "), "\n")
  cat("   🎯 Dredging component: μ =", round(mu[dredge_comp], 2), "kn (index", dredge_comp, ")\n")

  ais[, comp     := apply(.SD, 1, which.max), .SDcols = paste0("p", 1:(K + 1))]
  ais[, behavior := map[comp]]

  gd_archive_dir <- file.path(Sys.getenv("GMM_ARCHIVE_ROOT",file.path(output_dir,"gmm_models")),
    paste0(format(Sys.time(),"%Y%m%dT%H%M%S"),"_",Sys.getenv("SLURM_JOB_ID","local"),"_",Sys.getpid()))
  gd_archive(gd_archive_dir,gmm,gd_candidates,
    rows=list(mobile_ais_rows=idx_train,selection_X_rows=gd_bic_rows,final_X_rows=gd_fit_rows),
    normalisation=list(q95_acc=q95_acc,q95_turn=q95_turn,variables=colnames(X)),
    rng=gd_states,keys=as.data.frame(ais[idx_train,.(ssvid,Timestamp,Navire,Annee,Seg_id)]),
    metadata=list(origin="executed_pipeline",K=K,bic=bic_vals,candidates=G_choices,
      order_by_speed=ord,mobile_labels=mobile_labels,dredge_component=dredge_comp,
      final_warnings=gd_final_result$warnings,job=Sys.getenv("SLURM_JOB_ID"),
      script_md5=tools::md5sum("step3_merge_final.R"),
      helper_md5=tools::md5sum(Sys.getenv("GMM_DIAGNOSTICS_HELPER","gmm_diagnostics.R"))))
  cat("GMM model archive:",gd_archive_dir,"\n")

  ## 6. Run-length smoothing ---------------------------------------------------
  # Cost matrix adapted to the dynamic mapping
  all_labels <- c("stops", mobile_labels)
  n_labels <- length(all_labels)
  
  # Base cost matrix (high cost for unlikely transitions)
  cost_matrix <- matrix(5, n_labels, n_labels, 
                       dimnames = list(all_labels, all_labels))
  
  # Specific costs
  diag(cost_matrix) <- 0  # same state = cost 0
  
  # Costs between mobile states (natural transitions)
  # dredge_state avoids shadowing dredge_comp
  if ("dredging" %in% all_labels) {
    dredge_state <- which(all_labels == "dredging")  # index in all_labels (1..n_labels)
    # Dredging ↔ transit = moderate cost
    transit_labels_cost <- setdiff(all_labels, c("stops", "dredging"))
    for (tl in transit_labels_cost) {
      tl_idx <- which(all_labels == tl)
      cost_matrix[dredge_state, tl_idx] <- 3
      cost_matrix[tl_idx, dredge_state] <- 3
    }
  }

  # Costs to/from stops (rare transitions)
  stops_state <- which(all_labels == "stops")
  cost_matrix[stops_state, ] <- 4
  cost_matrix[, stops_state] <- 4
  cost_matrix[stops_state, stops_state] <- 0

  # Reduced cost for dredging ↔ stops (more natural)
  if ("dredging" %in% all_labels) {
    dredge_state <- which(all_labels == "dredging")
    cost_matrix[dredge_state, stops_state] <- 3
    cost_matrix[stops_state, dredge_state] <- 3
  }

  # Cost function safe against missing labels
  safe_cost <- function(a, b) {
    if (!a %in% rownames(cost_matrix) || !b %in% colnames(cost_matrix)) return(Inf)
    cost_matrix[a, b]
  }

  # Robust version of smooth_transition
  smooth_transition_robust <- function(lbl, spd, thr, dredge_min, dredge_max) {
    rle_lbl <- rle(lbl); cum <- cumsum(rle_lbl$lengths); beg <- c(1, head(cum, -1) + 1)
    for (i in seq_along(rle_lbl$lengths)) {
      if (rle_lbl$lengths[i] < thr) {
        vbar <- mean(spd[beg[i]:cum[i]])
        if (vbar <= dredge_max) {
          prev <- if (i > 1) rle_lbl$values[i - 1] else NA
          next_val <- if (i < length(rle_lbl$values)) rle_lbl$values[i + 1] else NA
          c_prev <- if (!is.na(prev)) safe_cost(prev, rle_lbl$values[i]) else Inf
          # Arguments in temporal order (current -> next). The order has no effect
          # while cost_matrix is symmetric, which it is by construction.
          c_next <- if (!is.na(next_val)) safe_cost(rle_lbl$values[i], next_val) else Inf
          rle_lbl$values[i] <- if (c_prev < c_next) prev else next_val
        }
      }
    }
    inverse.rle(rle_lbl)
  }

  # 75th percentile rather than the median, for a more stable smoothing threshold
  smoothing_threshold <- quantile(ais$n_min, 0.75, na.rm = TRUE)
  cat("   🎯 Smoothing threshold (75th percentile of n_min):", round(smoothing_threshold, 1), "\n")
  
  ais[, behavior_smooth :=
        smooth_transition_robust(behavior, Speed,
                                thr = smoothing_threshold,
                                dredge_min = 1, dredge_max = 3.5),
      by = Seg_id]

  ais[, dredge_lisse := as.integer(behavior_smooth == "dredging")]

  ## 7. Normalised signals for the dredging score -----------------------------
  ais[, norm_course := pmax(0, pmin(1, 1 - Course_change / q95_turn))]
  ais[, norm_acc    := pmax(0, pmin(1, abs(Accel) / q95_acc))]

  cat("✅  ROBUST GMM done –", table(ais$behavior_smooth), "states after smoothing\n")
  
  # CLASS CHECK AFTER GMM (dry-run)
  if (is_dry) {
    n_classes <- length(unique(ais$behavior_smooth))
    if (n_classes < 2) {
      stop(sprintf("FATAL ERROR: single-class DRY-RUN sample after GMM (%d class) – use another subset", n_classes))
    }
    cat("   ✅ Dry-run: ", n_classes, "classes detected after GMM\n")
  }
  
  # Save the mapping information for the next steps
  # Stored as attributes to avoid duplicating data in RAM
  # dredge_comp is already validated as unique
  attr(ais, "gmm_info") <- list(
    K = K,
    dredge_idx = dredge_comp,          # 1..K (index in mobile_labels)
    dredge_speed = mu[dredge_comp],
    mobile_labels = mobile_labels,
    mu = mu
  )
  
  ## ═══════════════════════════════════════════════════════════════════════
  ## QA REFACTORED: QA1 (pre-smooth, hard) + QA2 (post-smooth, soft)
  ## Uses dredge_comp (already validated) rather than recomputing it
  ## ═══════════════════════════════════════════════════════════════════════

  dredge_col <- paste0("p", dredge_comp + 1)  # +1 because p1 = stops
  ais[, p_dredge_qa := get(dredge_col)]

  ## ─────────────────  SANITY-CHECK  ─────────────────
  # 1) Replace any NA by 'unknown'
  ais[is.na(behavior_smooth), behavior_smooth := "unknown"]
  ## ───────────────────────────────────────────────────

  ## --- QA1: PRE-SMOOTH (hard check on argmax) ---
  # On behavior (before smoothing), "dredging" points MUST have p_dredge as argmax
  raw_dredge_mask <- ais$behavior == "dredging"
  n_raw_dredge <- sum(raw_dredge_mask, na.rm = TRUE)

  if (n_raw_dredge > 0) {
    proba_cols <- paste0("p", 1:(K+1))
    max_proba_col <- apply(ais[raw_dredge_mask, proba_cols, with=FALSE], 1, which.max)
    expected_col <- dredge_comp + 1  # +1 because p1 = stops

    n_correct <- sum(max_proba_col == expected_col, na.rm = TRUE)
    qa1_rate <- n_correct / n_raw_dredge

    if (qa1_rate >= 0.99) {
      cat(sprintf("   ✅ QA1 (pre-smooth): %.1f%% of dredging points have p_dredge=argmax\n", 100*qa1_rate))
    } else if (qa1_rate >= 0.90) {
      cat(sprintf("   ⚠️  QA1 (pre-smooth): only %.1f%% (expected ≥99%%) - check the mapping\n", 100*qa1_rate))
    } else {
      cat(sprintf("   ❌ QA1 (pre-smooth): FAILED - only %.1f%% (expected ≥90%%)\n", 100*qa1_rate))
      cat("      Diagnostic: argmax distribution =", paste(names(table(max_proba_col)), table(max_proba_col), sep=":", collapse=", "), "\n")
    }
  } else {
    cat("   ⚠️  QA1: no 'dredging' point before smoothing\n")
  }

  ## --- QA2: POST-SMOOTH (soft check on p_dredge distribution) ---
  # On behavior_smooth, "dredging" points should have a high p_dredge (median ≥ 0.5)
  # but smoothing can relabel segments, so this is only a warning
  smooth_dredge_mask <- ais$behavior_smooth == "dredging"
  n_smooth_dredge <- sum(smooth_dredge_mask, na.rm = TRUE)

  if (n_smooth_dredge > 0) {
    p_dredge_smooth <- ais[smooth_dredge_mask, p_dredge_qa]
    q <- quantile(p_dredge_smooth, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

    if (q["50%"] >= 0.5) {
      cat(sprintf("   ✅ QA2 (post-smooth): median p_dredge = %.3f (≥0.5)\n", q["50%"]))
    } else if (q["50%"] >= 0.3) {
      cat(sprintf("   ⚠️  QA2 (post-smooth): median p_dredge = %.3f (low, expected ≥0.5)\n", q["50%"]))
      cat(sprintf("      Quantiles: 10%%=%.3f, 25%%=%.3f, 50%%=%.3f, 75%%=%.3f, 90%%=%.3f\n",
                  q["10%"], q["25%"], q["50%"], q["75%"], q["90%"]))
    } else {
      cat(sprintf("   ❌ QA2 (post-smooth): median p_dredge = %.3f (very low!)\n", q["50%"]))
      cat(sprintf("      Quantiles: 10%%=%.3f, 25%%=%.3f, 50%%=%.3f, 75%%=%.3f, 90%%=%.3f\n",
                  q["10%"], q["25%"], q["50%"], q["75%"], q["90%"]))
      cat("      ⚠️  Smoothing may have extended the dredging segments too far\n")
    }
  } else {
    cat("   ⚠️  QA2: no 'dredging' point after smoothing\n")
  }
  
  ais
})

# --------------------------------------------------------------------------
# 13. ── GRID-SEARCH (dredging score weights), without data leakage
# --------------------------------------------------------------------------
# CHECK-POINT 4: grid-search done
resGS <- checkpoint(
  "stage3D_gridsearch.rds",
  {
    cat("🔍  Grid-search (dredging score weights), without data leakage\n")
    
    # Avoid data leakage
    # 1) Define the CV blocks (strict leave-one-year-out)
    cv_folds <- length(unique(ais$Annee))  # exactly 1 year = 1 fold
    ais[, fold := frank(Annee, ties.method="dense")]
    
    cat("   Strict LOYO cross-validation:", cv_folds, "folds (1 year = 1 fold)\n")
    
    # 2) Remove the label derived from the predictors (avoids data leakage)
    # Use p_dredge (GMM probability) + kinematic signals instead of dredge_lisse
    # === COMPOSITE DREDGING SCORE: 5 signals → cv.glmnet (ridge) → LOYO → Youden J ===
    ais[, speed_bin := as.integer(Speed >= 1 & Speed <= 3.5)]
    ais[, speed_norm := exp(-(Speed - 2.25)^2 / (2 * 1^2))]
    
    # Replace the remaining NAs in speed_bin and speed_norm
    ais[is.na(speed_bin),  speed_bin  := 0L]
    ais[is.na(speed_norm), speed_norm := 0]
    
    # Guard speed_norm against divisions by zero
    ais[!is.finite(speed_norm), speed_norm := 0]
    
    # Dynamic dredging probability
    # Retrieve the mapping information from the GMM attributes
    gmm_info <- attr(ais, "gmm_info")
    if (is.null(gmm_info)) {
      stop("❌ ERROR: GMM mapping information missing (attribute gmm_info)")
    }
    
    dredge_idx <- gmm_info$dredge_idx
    dredge_speed <- gmm_info$dredge_speed
    gmm_k <- gmm_info$K
    
    # Create the p_dredge column with the right probability
    dredge_col <- paste0("p", dredge_idx + 1)  # +1 because p1 = stops
    if (!dredge_col %in% names(ais)) {
      stop(sprintf("❌ ERROR: Column %s missing (dredge_idx=%d)", dredge_col, dredge_idx))
    }
    
    ais[, p_dredge := get(dredge_col)]
    cat(sprintf("   ✅ Dredging probability: %s (component %d, speed %.2f kn, K=%d)\n", 
                dredge_col, dredge_idx, dredge_speed, gmm_k))
    
    # Predictors independent of the final label
    predictors <- ais[, .(p_dredge, norm_course, norm_acc, speed_bin, speed_norm)]
    target <- ais$behavior_smooth == "dredging"
    
    # Check: exactly 5 signals
    stopifnot(identical(colnames(predictors), 
                        c("p_dredge", "norm_course", "norm_acc", "speed_bin", "speed_norm")))
    cat("   ✅ Check: exactly 5 signals required for the composite score\n")
    
    # Replace the remaining NAs in speed_bin and speed_norm before conversion
    ais[is.na(speed_bin),  speed_bin  := 0L]
    ais[is.na(speed_norm), speed_norm := 0]
    
    stopifnot(all(sapply(predictors, is.numeric)))
    
    # Check for non-finite values in the predictors
    msg <- colSums(!is.finite(as.matrix(predictors)))
    cat("⛔ Non-finite values per column:", paste(names(msg), msg, sep=":", collapse=", "), "\n")
    
    predictors <- predictors[, lapply(.SD, function(v) replace(v, !is.finite(v), 0))]
    cat("✅ Non-finite values replaced by 0\n")
    
    # Clean the NAs in the target (keeping the data)
    cat("   🔍 Cleaning the NAs in the target variable...\n")
    keep_lbl <- !is.na(target)
    n_before <- nrow(ais)
    n_filtered <- sum(!keep_lbl)
    
    cat(sprintf("   📊 Points with NA in target: %d (%.2f%%)\n", 
                n_filtered, 100 * n_filtered / n_before))
    
    if (any(is.na(target[keep_lbl]))) {
      stop("FATAL ERROR: NAs detected in target after filtering!")
    }
    
    # Convert to integer 0/1 (glmnet prefers numeric), for training only
    target_train <- as.integer(target[keep_lbl])  # 1 = dredging, 0 = non-dredging
    
    # Keep the indices to reinsert the predictions
    train_indices <- which(keep_lbl)
    cat(sprintf("   ✅ Data kept: %d points for training, %d points for full prediction\n",
                length(train_indices), n_before))
    
    # Handle NA values in the predictors
    # glmnet does not accept NA values
    cat("   Checking for missing values in the predictors...\n")
    na_counts <- colSums(is.na(predictors))
    cat("   NA values per predictor:", paste(names(na_counts), na_counts, sep=":", collapse=", "), "\n")
    
    # Check of the p_dredge mapping
    cat("   Checking the p_dredge mapping:\n")
    cat("   Available p columns:", paste(intersect(paste0("p",1:5), names(ais)), collapse=", "), "\n")
    cat("   p_dredge column created:", if("p_dredge" %in% names(ais)) "YES" else "NO", "\n")
    if ("p_dredge" %in% names(ais)) {
      cat("   p_dredge values:", paste(head(ais$p_dredge, 5), collapse=", "), "\n")
    }
    
    # Option 1: replace NAs by 0 (more conservative)
    predictors[is.na(predictors)] <- 0
    cat("   NA values replaced by 0\n")
    
    
    cat("   Predictors used:", paste(names(predictors), collapse = ", "), "\n")
    cat("   Target classes:", table(target), "\n")
    
    # STRICT CHECK OF THE CLASS DISTRIBUTION
    cat("\n[DIAGNOSTIC] Counts per behavior_smooth:\n")
    print(table(ais$behavior_smooth, useNA = "ifany"))
    cat("Counts per year:\n")
    print(table(ais$Annee))
    
    # Check of the class distribution (on the training data)
    dredging_count <- sum(target_train)
    non_dredging_count <- sum(!target_train)
    
    # Thresholds depend on the mode (dry-run vs production)
    if (is_dry) {
      # In dry-run, fewer data are tolerated but at least 2 classes are required
      if (dredging_count < 10) {
        stop(sprintf("FATAL ERROR: only %d 'dredging' points detected (< 10), sample too unbalanced for the dry-run!", dredging_count))
      }
      if (non_dredging_count < 50) {
        stop(sprintf("FATAL ERROR: only %d 'non-dredging' points detected (< 50), sample too unbalanced for the dry-run!", non_dredging_count))
      }
    } else {
      if (dredging_count < 100) {
        stop(sprintf("FATAL ERROR: only %d 'dredging' points detected (< 100), a valid CV is impossible!", dredging_count))
      }
      if (non_dredging_count < 1000) {
        stop(sprintf("FATAL ERROR: only %d 'non-dredging' points detected (< 1000), a valid CV is impossible!", non_dredging_count))
      }
    }
    
    # Long-format summary rather than a 3D table (avoids a memory blow-up)
    cat("   📊 Class distribution per year/vessel (top 20):\n")
    ais_train <- ais[keep_lbl]
    ais_train[, target_tmp := target_train]
    class_summary <- ais_train[, .N, by = .(Annee, Navire, target_tmp)][order(Annee, -N)]
    print(head(class_summary, 20))
    ais_train[, target_tmp := NULL]  # clean-up
    
    # Check of the classes before CV (on the training data)
    if (length(unique(target_train)) < 2) {
      stop(sprintf("FATAL ERROR: single-class dataset! Only %d class(es) detected, a valid CV is impossible.\n   📊 Distribution: %s", 
                   length(unique(target_train)), paste(table(target_train), collapse=", ")))
    } else {
      # Check of the offline packages before CV
      cat("   Checking the required packages...\n")
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("❌ Package 'glmnet' missing - load it via module or $R_LIBS_USER")
      }
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("❌ Package 'doParallel' missing - load it via module or $R_LIBS_USER")
      }
        library(glmnet)
      library(doParallel)
      
      # glmnet parallelism without conflict with data.table
      cat("   Configuring glmnet parallelism...\n")
      setDTthreads(1)  # Disable data.table parallelism
      Sys.setenv(OMP_NUM_THREADS = 1)  # Disable OpenMP
      registerDoParallel(cores = min(12, parallel::detectCores() - 2))
      cat("   Parallelism configured for glmnet (", min(12, parallel::detectCores() - 2), " cores)\n")
      
      # on.exit() outside the CV loop
on.exit({
  stopImplicitCluster()
  setDTthreads(0)  # Put data.table back in auto mode
  Sys.unsetenv("OMP_NUM_THREADS")
}, add = TRUE)

# TRAINING DATA PREPARATION
predictors_train <- predictors[keep_lbl, ]
ais_train <- ais[keep_lbl]

# Create the missing columns BEFORE the CV
# Matrix computation with exact name matching
missing <- setdiff(c("p_dredge", "norm_course", "norm_acc", "speed_bin", "speed_norm"), names(predictors_train))
if (length(missing) > 0) {
  cat("   Adding the missing columns:", paste(missing, collapse=", "), "\n")
  predictors_train[, (missing) := 0]
}

# Convert to a sparse matrix to save memory
cat("   Converting to a sparse matrix to save memory...\n")
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("❌ Package 'Matrix' missing - load it via module or $R_LIBS_USER")
}
library(Matrix)

# Sparse matrix for all rows (not only keep_lbl)
predictors_all <- Matrix::Matrix(as.matrix(predictors), sparse = TRUE)
predictors_all <- Matrix::drop0(predictors_all)   # drop explicit zeros
cat(sprintf("   Global sparse matrix created: %.1f MB\n", 
            object.size(predictors_all)/1024^2))

# Sparse matrix for the training data only
predictors_matrix <- predictors_all[keep_lbl, ]
predictors_matrix <- Matrix::drop0(predictors_matrix)   # drop explicit zeros
cat(sprintf("   Training sparse matrix: %.1f MB (vs %.1f MB dense)\n", 
            object.size(predictors_matrix)/1024^2, 
            object.size(as.matrix(predictors_train))/1024^2))

auc_fold <- numeric(max(ais_train$fold))
      best_models <- list()
      
      # Class distribution per fold (training data)
      cat("   📊 Class distribution per fold (training data):\n")
      fold_table <- table(target_train, ais_train$fold)
      print(fold_table)
      
      for (k in unique(ais_train$fold)) {
        cat(sprintf("     Fold %d/%d (year %d)...\n", k, max(ais_train$fold), unique(ais_train$Annee[ais_train$fold == k])))
        
        train <- ais_train$fold != k
        test <- ais_train$fold == k
        
        # FULL TRAIN/TEST FILTERING - PROTECTION AGAINST NA AND NON-FINITE VALUES
        Xtrain <- predictors_matrix[train, ]
        Ytrain <- target_train[train]
        Xtest <- predictors_matrix[test, ]
        Ytest <- target_train[test]
        
        # Local computation, nothing borrowed from the previous iteration
        complete_train <- !is.na(Ytrain)
        complete_test <- !is.na(Ytest)
        
        # Check of the size after NA filtering (relaxed in dry-run)
        min_train <- if (is_dry) 20 else 100
        min_test <- if (is_dry) 5 else 10
        if (sum(complete_train, na.rm = TRUE) < min_train || sum(complete_test, na.rm = TRUE) < min_test) {
          stop(sprintf("FATAL ERROR: fold %d too small after NA filtering (train: %d < %d, test: %d < %d)", 
                       k, sum(complete_train, na.rm = TRUE), min_train, sum(complete_test, na.rm = TRUE), min_test))
        }
        
        # LOYO mono-class check
        # Train MUST have 2 classes (otherwise no fit is possible)
        # Test may be mono-class → AUC = NA (no stop)
        train_has_2 <- length(unique(Ytrain[complete_train])) >= 2
        test_has_2 <- length(unique(Ytest[complete_test])) >= 2

        if (!train_has_2) {
          stop(sprintf("FATAL ERROR: LOYO fold %d - single-class training set (cannot fit)", k))
        }

        # Counts per class in this fold
        cat(sprintf("       Dredging (train/test): %d/%d, Non-dredging (train/test): %d/%d\n",
          sum(Ytrain[complete_train], na.rm = TRUE), sum(Ytest[complete_test], na.rm = TRUE),
          sum(!Ytrain[complete_train], na.rm = TRUE), sum(!Ytest[complete_test], na.rm = TRUE)))

        if (!test_has_2) {
          cat("       ⚠️  Single-class test → AUC = NA for this fold\n")
        }

        # Penalised logistic regression (ridge), fitted without NA, with stratified folds
        tryCatch({
          # Stratified folds to avoid mono-class internal folds
          set.seed(1)

          # Stratification function
          make_stratified_foldid <- function(y, k = 5, seed = 1) {
            set.seed(seed)
            y <- as.integer(y)
            n1 <- sum(y == 1); n0 <- sum(y == 0)
            k_eff <- min(k, n1, n0)  # avoids empty folds
            if (k_eff < 2) k_eff <- 2  # minimum 2 folds
            foldid <- integer(length(y))
            foldid[y == 1] <- sample(rep(1:k_eff, length.out = n1))
            foldid[y == 0] <- sample(rep(1:k_eff, length.out = n0))
            foldid
          }
          
          foldid <- make_stratified_foldid(Ytrain[complete_train], k = 5, seed = k)

          # Safeguard: if some internal folds are mono-class, fall back to deviance
          folds_ok <- all(tapply(Ytrain[complete_train], foldid,
                                 function(z) length(unique(z)) >= 2))
          measure_type <- if (folds_ok) "auc" else "deviance"
          if (!folds_ok) cat("       ⚠️  Unbalanced internal folds → deviance fallback\n")

          # foldid already has the right length (built on Ytrain[complete_train])
          fit <- cv.glmnet(Xtrain[complete_train, ], Ytrain[complete_train],
                          family = "binomial", alpha = 0,  # Ridge
                          foldid = foldid, type.measure = measure_type,
                          parallel = FALSE)  # avoids nested parallelism

          # Prediction on the test fold, without NA
          prob <- predict(fit, Xtest[complete_test, ], s = "lambda.min", type = "response")

          # AUC (uses test_has_2 computed above)
          if (!test_has_2) {
            auc_fold[k] <- NA_real_
          } else {
            roc_test <- pROC::roc(Ytest[complete_test], prob, quiet = TRUE)
            auc_fold[k] <- as.numeric(pROC::auc(roc_test))
          }
          
          best_models[[as.character(k)]] <- fit
          cat(sprintf("       AUC fold %d: %.4f\n", k, auc_fold[k]))
          
        }, error = function(e) {
          cat(sprintf("       Error fold %d: %s\n", k, e$message))
          auc_fold[k] <- NA
        })
      }
      
      # Cross-validation results - protection against empty folds
      valid_auc <- auc_fold[!is.na(auc_fold)]
      
      # Protection against the failure of every fold
      if (all(is.na(auc_fold))) {
        stop("All folds failed: check for non-finite values in the predictors.")
      }
      
      if (length(valid_auc) > 0 && length(best_models) > 0) {
        mean_auc <- mean(valid_auc)
        sd_auc <- sd(valid_auc)
        cat(sprintf("   Mean CV AUC: %.4f ± %.4f (n=%d folds)\n", mean_auc, sd_auc, length(valid_auc)))
        
        # Safe selection of the best fold
        best_fold <- which.min(abs(replace(auc_fold, is.na(auc_fold), Inf) - mean_auc))
        best_model <- best_models[[as.character(best_fold)]]
        
        # Coefficients for the final score
        coefs <- coef(best_model, s = "lambda.min")
        coefs_named <- as.numeric(coefs)
        names(coefs_named) <- rownames(coefs)
        
        # Normalise the coefficients into weights
        feature_coefs <- coefs_named[-1]  # exclude the intercept
        total_weight <- sum(abs(feature_coefs))
        
        # Check of the non-zero coefficients
        if (total_weight == 0 || length(feature_coefs) == 0) {
          stop("FATAL ERROR: regression without non-zero weights - check the predictors")
        }
        
        weights <- abs(feature_coefs) / total_weight
        
        # === WEIGHT EXPORT: normalised coefficients ===
        fwrite(data.table(feature = names(weights), weight = weights),
               file.path(output_dir, "drag_score_weights.csv"))
        cat("   💾 Composite score weights saved in drag_score_weights.csv\n")
        
        # Final score
        # The missing columns were already added before the CV
        # Prediction on **all** rows (including unlabelled ones)
        full_prob <- as.vector(predict(best_model, predictors_all, s = "lambda.min", type = "response"))
        ais[, drag_score := full_prob]
        
        # Safeguard for NA values in drag_score
        ais[is.na(drag_score), drag_score := 0]
        
        # Optimal threshold from the CV (on the labelled data)
        cat("   Computing the optimal threshold...\n")
        # Probabilities restricted to the labelled rows, to calibrate the threshold
        prob_cv <- full_prob[keep_lbl]   # keeps the indexing correct
        roc_full <- pROC::roc(target_train, prob_cv, quiet = TRUE)
        optimal_threshold <- as.numeric(pROC::coords(roc_full, "best", ret = "threshold")[["threshold"]])
        best <- list(
          weights = weights,
          features = names(weights),
          mean_auc = mean_auc,
          sd_auc = sd_auc,
          auc_folds = auc_fold,
          n_folds = length(valid_auc),
          method = "Ridge_CV",
          optimal_threshold = optimal_threshold
        )
        cat(sprintf("   Optimal threshold: %.3f\n", optimal_threshold))
        
        # Apply the optimal threshold to all rows
        ais[, Dragage_flag := as.integer(drag_score >= optimal_threshold)]

        ## --- Smoothing of the final flag ---
        ## The composite score re-estimates the label from p_dredge (the raw GMM
        ## probability) and memoryless kinematic signals, so the run-length smoothing
        ## of step 2 is not on its data path, and Dragage_flag flickers.
        ## The smoothing is applied here, at the end of the chain.
        ## Set DREDGE_SMOOTH_FLAG=0 to recover the published behaviour.
        smooth_binary_runs <- function(flag, spd, thr, spd_max) {
          r <- rle(flag)
          n <- length(r$lengths)
          if (n < 2L) return(flag)
          cum <- cumsum(r$lengths); beg <- c(1L, utils::head(cum, -1L) + 1L)
          for (i in seq_len(n)) {
            if (r$lengths[i] < thr) {
              vbar <- mean(spd[beg[i]:cum[i]], na.rm = TRUE)
              if (!is.na(vbar) && vbar <= spd_max) {
                ## binary: runs alternate, so the replacement neighbour is
                ## unique -- the multi-state cost matrix is not needed here
                prev <- if (i > 1L) r$values[i - 1L] else NA_integer_
                nxt  <- if (i < n)  r$values[i + 1L] else NA_integer_
                r$values[i] <- if (!is.na(prev)) prev else nxt
              }
            }
          }
          inverse.rle(r)
        }

        if (as.integer(Sys.getenv("DREDGE_SMOOTH_FLAG", "1")) == 1L) {
          .thr <- as.numeric(quantile(ais$n_min, 0.75, na.rm = TRUE))
          .bef  <- sum(ais$Dragage_flag, na.rm = TRUE)
          ais[, Dragage_flag := smooth_binary_runs(Dragage_flag, Speed,
                                                   thr = .thr, spd_max = 3.5),
              by = Seg_id]
          .aft <- sum(ais$Dragage_flag, na.rm = TRUE)
          cat(sprintf("   Flag smoothing (run threshold %.1f): %d -> %d pings (%+.1f %%)\n",
                      .thr, .bef, .aft, 100 * (.aft / max(.bef, 1L) - 1)))
        } else {
          cat("   Flag smoothing DISABLED (DREDGE_SMOOTH_FLAG=0)\n")
        }

        
        # Free memory after CV
        rm(best_models, predictors_train, auc_fold, valid_auc, prob_cv, roc_full)
        gc()
        cat("   Memory freed after CV\n")
      } else {
        stop(sprintf("FATAL ERROR: no valid fold for the CV!\n   🔍 Diagnostic: no valid fold found\n   📊 AUC summary per fold: %s\n   📊 Number of stored models: %d\n   📊 Check the class distribution per year.", 
                     paste(sprintf("%.3f", auc_fold), collapse=", "), length(best_models)))
      }
    }
    
    cat("✅ Grid-search done - honest AUC computed\n")
    list(ais = ais, best = best, auc = if(exists("valid_auc")) valid_auc else NA_real_)
  },
  force_recompute = TRUE
)

ais <- resGS$ais
best <- resGS$best
auc <- resGS$auc

# Save the GMM probabilities for QA before clean-up
cat("💾 Saving the GMM probabilities for quality control...\n")

# Handle the ssvid vs MMSI column
if ("ssvid" %in% names(ais)) {
  id_col <- "ssvid"
} else if ("MMSI" %in% names(ais)) {
  id_col <- "MMSI"
} else if ("Navire" %in% names(ais)) {
  id_col <- "Navire"
} else {
  warning("⚠️  No identifier column found (ssvid/MMSI/Navire) - using the first column")
  id_col <- names(ais)[1]
}

# Robust QA: keep only the existing p* columns (depends on K)
gmm_info <- attr(ais, "gmm_info")
if (!is.null(gmm_info)) {
  qa_p <- intersect(paste0("p", 1:(gmm_info$K + 1)), names(ais))
} else {
  qa_p <- intersect(paste0("p",1:5), names(ais))
}
gmm_qa_data <- ais[, c(id_col, "Timestamp", qa_p,
                       "behavior_smooth","drag_score","Dragage_flag"), with = FALSE]
setnames(gmm_qa_data, id_col, "id")

saveRDS(gmm_qa_data, file.path(output_dir, "GMM_probabilities_QA.rds"), compress = "xz")
cat("   GMM probabilities saved in GMM_probabilities_QA.rds (columns p:", paste(qa_p, collapse=", "), ")\n")

  # Safe purge of the temporary columns, depending on the GMM mapping
# Keeps dredge_lisse, norm_course, norm_acc, Course_change and Accel for QC
  # p_dredge is kept in the final dataset (not in tmp_cols)
  
  # GMM mapping information, to identify the column to keep
  gmm_info <- attr(ais, "gmm_info")
  if (!is.null(gmm_info)) {
    dredge_idx <- gmm_info$dredge_idx
    dredge_col <- paste0("p", dredge_idx + 1)  # +1 because p1 = stops
    cat("   🎯 Dredging column kept:", dredge_col, "\n")
    
    # Remove every p* column EXCEPT p_dredge and p1 (stops)
    all_p_cols <- paste0("p", 1:(gmm_info$K + 1))
    tmp_cols <- setdiff(all_p_cols, c("p1", dredge_col))
  } else {
    # Fallback without GMM information
    tmp_cols <- c("p2","p4","p5")
  }
  
  # Add the other temporary columns
  # The smoothing input and the GMM argmax are kept: without them, the
    # smoothing can be neither audited nor varied without rerunning step 3.
    tmp_cols <- c(tmp_cols, "speed_norm","speed_bin","p_dredge_qa")

# Safe purge: only removes existing columns
keep <- intersect(tmp_cols, names(ais))
if (length(keep)) {
  ais[, (keep) := NULL]
  cat("   Temporary columns removed:", paste(keep, collapse=", "), "\n")
} else {
  cat("   No temporary column to remove\n")
}

# -- Clean-up of the superfluous stats columns -----------------------------
todrop <- intersect(c("dt_sec", "p50_dt", "prop_long"), names(ais))
if (length(todrop)) {
  ais[, (todrop) := NULL]
  cat("   Stats columns removed:", paste(todrop, collapse=", "), "\n")
} else {
  cat("   No stats column to remove\n")
}

# --------------------------------------------------------------------------
# 14. ── SAVING  ---------------------------------------------------------------
# --------------------------------------------------------------------------

## ─────────────────  FINAL GLOBAL LOCK  ─────────────────
# No NA may remain in behavior_smooth before the export
stopifnot(!anyNA(ais$behavior_smooth))
cat("✅ Global lock: no NA in behavior_smooth - pipeline sealed\n")
## ─────────────────────────────────────────────────────────

# Saved AFTER the grid-search, to include Dragage_flag
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
saveRDS(ais, file.path(output_dir,
         paste0("AIS_data_core_preprocessed_V6_", timestamp, ".rds")),
       compress = "xz")

# Save the results
saveRDS(list(best_weights = best,
             best_auc     = if(!is.null(best$mean_auc)) best$mean_auc else NA_real_,
             auc_sd       = if(!is.null(best$sd_auc)) best$sd_auc else NA_real_,
             auc_folds    = if(!is.null(best$auc_folds)) best$auc_folds else numeric(0),  # full AUCs per fold
             n_folds      = if(!is.null(best$n_folds)) best$n_folds else 0,
             method       = if(!is.null(best$method)) best$method else "unknown",
             timestamp    = timestamp),
        file.path(output_dir,
         paste0("dragage_gridsearch_results_V6_", timestamp, ".rds")),
        compress = "xz")

cat("\n💾  Files saved in", output_dir, "\n")

# QA CLEAN-UP: remove the p_dredge_qa test column
if ("p_dredge_qa" %in% names(ais)) {
  ais[, p_dredge_qa := NULL]
  cat("✓ QA column p_dredge_qa removed\n")
}

# EARLY EXIT IN DRY-RUN MODE
if (is_dry) {
  cat("\n✅ === DRY-RUN COMPLETED SUCCESSFULLY ===\n")
  cat("📝 Full test up to the grid-search - no heavy file saved\n")
  cat("🎯 Pipeline ready for the full run\n")
  cat("⏱️  End dry-run:", format(Sys.time()), "\n")
  quit(save = "no")
}

cat("⏱️  End:", format(Sys.time()), "\n")
