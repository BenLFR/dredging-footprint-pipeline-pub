#!/usr/bin/env Rscript
# ====================================================================
# STEP 2: INDIVIDUAL VESSEL PROCESSING (ARRAY JOB)
# V8 - Tiled spatial land mask + anti false-positive on_land
#   - Splits the vessel footprint into spatial tiles (tile_deg°)
#   - Loads the land mask PER TILE (avoids OOM & timeout)
#   - Per-tile cache for crosses_land in the iterative loop
#   - Single-bbox fallback if tiling disabled or too many tiles
#   - Static pre-computations: on_land + near_land outside the loop
#   - near_land via st_is_within_distance (avoids costly buffer/union)
#   - crosses_land in metric CRS (EPSG:3857) + time+distance filters BEFORE geometry
#   - Reduction of on_land false positives: slight erosion of the mask (on_land only)
#   - DO NOT TOUCH stop_rate / is_stop (downstream)
# ====================================================================

cat("=== INDIVIDUAL VESSEL PROCESSING (OPTIMISED V8 - TILED) ===\n")
cat("SCRIPT_VERSION: step2_process_vessel.R 2026-01-27 tiled-v8-onlandfp\n")
cat("Start:", format(Sys.time()), "\n")

# CRITICAL R CONFIGURATION - BEFORE LOADING PACKAGES
.libPaths("~/R/library")
cat("R configured with library:", .libPaths()[1], "\n")

# ---- SINGLE-THREAD CONFIGURATION ----
options(mc.cores = 1)
Sys.setenv(MC_CORES = 1)
Sys.setenv(DT_GForce = "FALSE")
Sys.setenv(OMP_NUM_THREADS = 1)

# ---- LOADING PACKAGES ----
suppressPackageStartupMessages({
  library(data.table)

  use_qs <- FALSE
  tryCatch({
    if (requireNamespace("qs", quietly = TRUE)) {
      library(qs)
      use_qs <- TRUE
      cat("Package qs available\n")
    }
  }, error = function(e) {
    cat("Package qs not available - fallback RDS\n")
  })

  library(lubridate)
  library(yaml)
  library(solitude)
  library(mclust)
  library(zoo)
  library(sf)
})

# Disable S2 (spherical geometry) - use GEOS
sf::sf_use_s2(FALSE)
cat("S2 disabled - using GEOS\n")

setDTthreads(1)
options(datatable.optimize = 1)
cat("Packages loaded - single-thread configuration enabled\n")

# ---- HELPERS ----
parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(default)
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

safe_num <- function(x, default) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) default else v
}

safe_int <- function(x, default) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

# ---- GEOSPATIAL PARAMETERS ----
land_mask_path_raw <- Sys.getenv("LAND_MASK_FILE", unset = file.path(Sys.getenv("DFP_DATA"), "configuration/land_mask/land_polygons.shp"))
land_mask_path     <- path.expand(land_mask_path_raw)
land_buffer_m      <- safe_num(Sys.getenv("LAND_MASK_BUFFER_M", unset = "0"), 0)

near_coast_km <- safe_num(Sys.getenv("LAND_NEAR_COAST_KM", unset = "5"), 5)
if (!is.finite(near_coast_km) || near_coast_km <= 0) near_coast_km <- 5

# IMPORTANT: on_land-only erosion to reduce coastline/port false positives
land_onland_erode_m <- safe_num(Sys.getenv("LAND_ONLAND_ERODE_M", unset = "-100"), -100)
if (!is.finite(land_onland_erode_m)) land_onland_erode_m <- -100

# crosses_land gating (reduce port jitter false positives)
min_seg_dt_sec  <- safe_int(Sys.getenv("MIN_SEG_DT_SEC", unset = "60"), 60)
if (!is.finite(min_seg_dt_sec) || min_seg_dt_sec < 0) min_seg_dt_sec <- 60

min_seg_dist_m  <- safe_num(Sys.getenv("MIN_SEG_DIST_M", unset = "250"), 250)
if (!is.finite(min_seg_dist_m) || min_seg_dist_m < 0) min_seg_dist_m <- 250

cat("\n--- Geospatial configuration ---\n")
cat("LAND_MASK_FILE (raw):", land_mask_path_raw, "\n")
cat("LAND_MASK_FILE (expanded):", land_mask_path, "\n")
cat("LAND_MASK_BUFFER_M:", land_buffer_m, "\n")
cat("LAND_NEAR_COAST_KM:", near_coast_km, "\n")
cat("LAND_ONLAND_ERODE_M:", land_onland_erode_m, " (on_land-only)\n")
cat("MIN_SEG_DT_SEC:", min_seg_dt_sec, " (crosses_land gating)\n")
cat("MIN_SEG_DIST_M:", min_seg_dist_m, " (crosses_land gating)\n")

default_speed_kn    <- safe_num(Sys.getenv("MAX_JUMP_SPEED_KN", unset = "30"), 30)
if (!is.finite(default_speed_kn) || default_speed_kn <= 0) default_speed_kn <- 30

spike_dist_min_nm   <- safe_num(Sys.getenv("SPIKE_DIST_MIN_NM", unset = "1"), 1)
if (!is.finite(spike_dist_min_nm) || spike_dist_min_nm <= 0) spike_dist_min_nm <- 1

spike_bridge_max_nm <- safe_num(Sys.getenv("SPIKE_BRIDGE_MAX_NM", unset = "0.3"), 0.3)
if (!is.finite(spike_bridge_max_nm) || spike_bridge_max_nm <= 0) spike_bridge_max_nm <- 0.3

sf_chunk <- safe_int(Sys.getenv("SF_CHUNK", unset = "100000"), 100000)
if (!is.finite(sf_chunk) || sf_chunk < 1000) sf_chunk <- 100000
cat("SF_CHUNK:", sf_chunk, "\n")

# ---- TILING PARAMETERS ----
land_tile_enable <- parse_bool(Sys.getenv("LAND_TILE_ENABLE", unset = "TRUE"), TRUE)
land_tile_deg    <- safe_num(Sys.getenv("LAND_TILE_DEG", unset = "5"), 5)
if (!is.finite(land_tile_deg) || land_tile_deg <= 0) land_tile_deg <- 5

land_tile_max <- safe_int(Sys.getenv("LAND_TILE_MAX_TILES", unset = "200"), 200)
if (!is.finite(land_tile_max) || land_tile_max < 1) land_tile_max <- 200

cat("LAND_TILE_ENABLE:", land_tile_enable, "\n")
cat("LAND_TILE_DEG:", land_tile_deg, "\n")
cat("LAND_TILE_MAX_TILES:", land_tile_max, "\n")
cat("---\n\n")

# ---- UTILITY FUNCTIONS ----
haversine_nm <- function(lat1, lon1, lat2, lon2) {
  r <- 6371000
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  (r * c) / 1852
}

chunk_any_intersects <- function(x, y, chunk = 100000, label = "intersects") {
  n <- nrow(x)
  out <- logical(n)
  if (n == 0) return(out)

  n_y <- if (inherits(y, "sfc")) length(y) else nrow(y)
  n_chunks <- ceiling(n / chunk)
  cat("        [", label, "] ", n, " pts / ", n_y, " polys (", n_chunks, " chunks)...", sep = "")

  for (i in seq(1, n, by = chunk)) {
    j <- min(i + chunk - 1, n)
    out[i:j] <- suppressMessages(lengths(sf::st_intersects(x[i:j, ], y)) > 0)
  }
  cat(" done\n")
  out
}

# Near-land via st_is_within_distance (avoids OOM of st_union+st_buffer)
chunk_any_within_distance <- function(x, y, dist_m, chunk = 100000, label = "within_dist") {
  n <- nrow(x)
  out <- logical(n)
  if (n == 0) return(out)

  n_y <- if (inherits(y, "sfc")) length(y) else nrow(y)
  n_chunks <- ceiling(n / chunk)
  cat("        [", label, "] ", n, " pts / ", n_y, " polys, dist=", dist_m/1000, "km (", n_chunks, " chunks)...", sep = "")

  for (i in seq(1, n, by = chunk)) {
    j <- min(i + chunk - 1, n)
    out[i:j] <- suppressMessages(lengths(sf::st_is_within_distance(x[i:j, ], y, dist = dist_m)) > 0)
  }
  cat(" done\n")
  out
}

# ---- TILING UTILITIES ----
tile_assign <- function(lon, lat, tile_deg) {
  tx <- floor(lon / tile_deg)
  ty <- floor(lat / tile_deg)
  paste0(tx, "_", ty)
}

tile_bbox_4326 <- function(tile_key, tile_deg, margin_deg = 0) {
  parts <- as.integer(strsplit(tile_key, "_", fixed = TRUE)[[1]])
  tx <- parts[1]; ty <- parts[2]
  xmin <- tx * tile_deg - margin_deg
  ymin <- max(ty * tile_deg - margin_deg, -90)
  xmax <- (tx + 1) * tile_deg + margin_deg
  ymax <- min((ty + 1) * tile_deg + margin_deg, 90)
  sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = sf::st_crs(4326))
}

cache_new  <- function() new.env(hash = TRUE, parent = emptyenv())
cache_get  <- function(cache, key) {
  if (exists(key, envir = cache, inherits = FALSE)) get(key, envir = cache, inherits = FALSE)
  else NULL
}
cache_set  <- function(cache, key, value) assign(key, value, envir = cache)

# --- Local land mask via WKT filter (avoids global loading) ---
prefer_land_datasource <- function(path) {
  # If a .gpkg exists next to the .shp, prefer it (often faster + spatial index)
  if (grepl("\\.shp$", path, ignore.case = TRUE)) {
    gpkg <- sub("\\.shp$", ".gpkg", path, ignore.case = TRUE)
    if (file.exists(gpkg)) return(gpkg)
  }
  path
}

# Lightweight tile loader (minimal output for per-tile loop)
load_land_tile <- function(path, bbox_4326, buffer_m = 0) {
  path_exp <- path.expand(prefer_land_datasource(path))
  if (!file.exists(path_exp)) return(NULL)

  bbox_sfc <- sf::st_as_sfc(bbox_4326)
  sf::st_crs(bbox_sfc) <- 4326
  wkt <- sf::st_as_text(bbox_sfc)

  land <- tryCatch(sf::st_read(path_exp, wkt_filter = wkt, quiet = TRUE),
                   error = function(e) NULL)
  if (is.null(land) || !inherits(land, "sf") || nrow(land) == 0) return(NULL)

  if (isFALSE(sf::st_is_longlat(land))) {
    land <- tryCatch(sf::st_transform(land, 4326), error = function(e) land)
  }
  land <- sf::st_make_valid(land)
  land <- land[!sf::st_is_empty(land), ]
  if (nrow(land) == 0) return(NULL)

  land_m <- tryCatch(sf::st_transform(land, 3857), error = function(e) NULL)
  if (is.null(land_m)) return(NULL)

  if (is.finite(buffer_m) && buffer_m != 0) {
    land_m <- tryCatch({
      x <- sf::st_buffer(land_m, dist = buffer_m)
      x <- sf::st_make_valid(x)
      x <- x[!sf::st_is_empty(x), ]
      x
    }, error = function(e) land_m)
  }

  land_m
}

make_land_onland_only <- function(land_m, erode_m) {
  # Erosion (negative buffer) only for on_land PIP test.
  if (is.null(land_m) || !inherits(land_m, "sf") || nrow(land_m) == 0) return(land_m)
  if (!is.finite(erode_m) || erode_m == 0) return(land_m)

  out <- tryCatch({
    x <- sf::st_buffer(land_m, dist = erode_m)
    x <- sf::st_make_valid(x)
    x <- x[!sf::st_is_empty(x), ]
    if (inherits(x, "sf") && nrow(x) > 0) x else land_m
  }, error = function(e) land_m)

  out
}

# ============================================================
# OPTIMISED FUNCTION - Dynamic flags only
# ============================================================
flag_geospatial_anomalies_instrumented <- function(dt, land_local_m, default_speed_kn,
                                                   spike_dist_min_nm = 1, spike_bridge_max_nm = 0.3,
                                                   land_cache = NULL, tile_key_col = NULL,
                                                   min_seg_dt_sec = 60, min_seg_dist_m = 250) {

  dt[, `:=`(
    flag_on_land       = if ("flag_on_land_static" %in% names(dt)) flag_on_land_static else FALSE,
    flag_crosses_land  = FALSE,
    flag_speed_jump    = FALSE,
    flag_long_jump     = FALSE,
    flag_spike         = FALSE
  )]

  # Dynamic computations (neighbours change after removal)
  dt[, `:=`(
    gc_nm  = haversine_nm(Lat, Lon, shift(Lat), shift(Lon)),
    dt_sec = as.numeric(delta_t)
  )]

  dt[, avg_speed_kn := ifelse(!is.na(dt_sec) & dt_sec > 0, gc_nm / dt_sec * 3600, NA_real_)]

  # max plausible
  dt[, max_plausible_kn := default_speed_kn]
  if ("Service_speed" %in% names(dt)) {
    dt[!is.na(Service_speed), max_plausible_kn := Service_speed * 1.5]
  }

  # 1) Log on_land (static)
  n_on_land <- sum(dt$flag_on_land, na.rm = TRUE)
  if (n_on_land > 0) cat("      on_land (static):", n_on_land, "\n")

  # 2) crosses_land (semi-dynamic) - tiled or single-land
  use_tiled_crosses <- !is.null(land_cache) && !is.null(tile_key_col) && tile_key_col %in% names(dt)
  has_near_land <- "flag_near_land_static" %in% names(dt)
  can_do_crosses <- has_near_land && (use_tiled_crosses || (!is.null(land_local_m) && inherits(land_local_m, "sf") && nrow(land_local_m) > 0))

  if (can_do_crosses) {
    near_curr <- dt$flag_near_land_static
    near_next <- shift(near_curr, type = "lead", fill = FALSE)
    on_curr   <- dt$flag_on_land

    seg_candidates <- which((near_curr | near_next) & !on_curr)
    seg_candidates <- seg_candidates[seg_candidates < nrow(dt)]
    n0 <- length(seg_candidates)

    # Time filter (hard cap)
    next_dt_sec <- dt$dt_sec[seg_candidates + 1]
    keep1 <- !is.na(next_dt_sec) & next_dt_sec < 3600
    seg_candidates <- seg_candidates[keep1]
    next_dt_sec <- next_dt_sec[keep1]
    n1 <- length(seg_candidates)

    # Time filter (min)
    keep2 <- next_dt_sec >= min_seg_dt_sec
    seg_candidates <- seg_candidates[keep2]
    n2 <- length(seg_candidates)

    # Finite coordinates filter
    ok_seg <- is.finite(dt$Lon[seg_candidates]) & is.finite(dt$Lat[seg_candidates]) &
      is.finite(dt$Lon[seg_candidates + 1]) & is.finite(dt$Lat[seg_candidates + 1])
    seg_candidates <- seg_candidates[ok_seg]
    n3 <- length(seg_candidates)

    # Distance filter (reduces quay micro-movements)
    if (length(seg_candidates) > 0 && min_seg_dist_m > 0) {
      dist_nm <- haversine_nm(
        dt$Lat[seg_candidates], dt$Lon[seg_candidates],
        dt$Lat[seg_candidates + 1], dt$Lon[seg_candidates + 1]
      )
      dist_m <- dist_nm * 1852
      keep3 <- is.finite(dist_m) & dist_m >= min_seg_dist_m
      seg_candidates <- seg_candidates[keep3]
    }
    n4 <- length(seg_candidates)

    if (n0 > 0) {
      cat("      crosses_candidates:", n0,
          " -> dt<3600:", n1,
          " -> dt>=min:", n2,
          " -> coords:", n3,
          " -> dist>=min:", n4, "\n")
    }

    if (length(seg_candidates) > 0) {
      if (use_tiled_crosses) {
        # --- Per-tile crosses_land ---
        seg_tiles <- dt[[tile_key_col]][seg_candidates]
        unique_seg_tiles <- unique(seg_tiles)
        total_crosses <- 0L

        for (tk in unique_seg_tiles) {
          tile_land <- cache_get(land_cache, tk)
          if (is.null(tile_land)) next
          if (inherits(tile_land, "sf") && nrow(tile_land) == 0) next

          tile_mask <- seg_tiles == tk
          tile_segs <- seg_candidates[tile_mask]

          seg_lines_ll <- sf::st_sfc(
            mapply(function(i) {
              sf::st_linestring(matrix(
                c(dt$Lon[i], dt$Lat[i], dt$Lon[i + 1], dt$Lat[i + 1]),
                ncol = 2, byrow = TRUE
              ))
            }, tile_segs, SIMPLIFY = FALSE),
            crs = 4326
          )

          seg_lines_m <- tryCatch(sf::st_transform(seg_lines_ll, 3857), error = function(e) seg_lines_ll)
          crosses <- lengths(sf::st_intersects(seg_lines_m, tile_land)) > 0

          if (any(crosses)) {
            bad_idx <- tile_segs[crosses] + 1L
            dt[bad_idx, flag_crosses_land := TRUE]
            total_crosses <- total_crosses + length(bad_idx)
          }
        }

        if (total_crosses > 0) cat("      crosses_land:", total_crosses, "(tiled)\n")

      } else {
        # --- Single-land crosses_land (fallback) ---
        seg_lines_ll <- sf::st_sfc(
          mapply(function(i) {
            sf::st_linestring(matrix(
              c(dt$Lon[i], dt$Lat[i], dt$Lon[i + 1], dt$Lat[i + 1]),
              ncol = 2, byrow = TRUE
            ))
          }, seg_candidates, SIMPLIFY = FALSE),
          crs = 4326
        )

        seg_lines_m <- tryCatch(sf::st_transform(seg_lines_ll, 3857), error = function(e) seg_lines_ll)
        crosses <- lengths(sf::st_intersects(seg_lines_m, land_local_m)) > 0

        if (any(crosses)) {
          bad_idx <- seg_candidates[crosses] + 1L
          dt[bad_idx, flag_crosses_land := TRUE]
          cat("      crosses_land:", length(bad_idx), "\n")
        }
      }
    }
  }

  # 3) Speed jumps
  speed_jump_mask <- !is.na(dt$avg_speed_kn) & dt$avg_speed_kn > dt$max_plausible_kn
  dt[speed_jump_mask, flag_speed_jump := TRUE]
  cat("      speed_jump:", sum(speed_jump_mask), "\n")

  # 4) Distance jumps
  long_jump_mask <- !is.na(dt$gc_nm) & !is.na(dt$dt_sec) & dt$gc_nm > 100 & dt$dt_sec <= 3600
  dt[long_jump_mask, flag_long_jump := TRUE]
  cat("      long_jump:", sum(long_jump_mask), "\n")

  # 5) Spikes
  dt[, `:=`(
    dist_prev      = haversine_nm(shift(Lat), shift(Lon), Lat, Lon),
    dist_next      = haversine_nm(Lat, Lon, shift(Lat, type = "lead"), shift(Lon, type = "lead")),
    dist_prev_next = haversine_nm(shift(Lat), shift(Lon), shift(Lat, type = "lead"), shift(Lon, type = "lead"))
  )]

  spike_mask <- !is.na(dt$dist_prev) & !is.na(dt$dist_next) & !is.na(dt$dist_prev_next) &
    dt$dist_prev > spike_dist_min_nm &
    dt$dist_next > spike_dist_min_nm &
    dt$dist_prev_next < spike_bridge_max_nm
  dt[spike_mask, flag_spike := TRUE]
  cat("      spike:", sum(spike_mask), "\n")

  dt[, geo_flag := flag_on_land | flag_crosses_land | flag_speed_jump | flag_long_jump | flag_spike]

  dt[, c("gc_nm", "dt_sec", "avg_speed_kn", "max_plausible_kn",
         "dist_prev", "dist_next", "dist_prev_next") := NULL]

  dt
}

# ---- ENVIRONMENT PARAMETERS ----
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
split_job_id <- Sys.getenv("SPLIT_JOB_ID")

cat("Task ID:", task_id, "\n")
cat("Split Job ID:", split_job_id, "\n")

# ---- FILE PATHS ----
split_dir <- file.path(Sys.getenv("DFP_DATA"), paste0("ais_split_", split_job_id))
metadata_file <- file.path(split_dir, "navires_metadata.csv")

if (!file.exists(metadata_file)) stop("Metadata not found: ", metadata_file)

metadata <- fread(metadata_file)
cat("Metadata loaded:", nrow(metadata), "vessels\n")

if (task_id > nrow(metadata)) stop("Task ID ", task_id, " > number of vessels ", nrow(metadata))

# ---- VESSEL SELECTION ----
vessel_info <- metadata[task_id]
input_file  <- vessel_info$file_path
vessel_name <- vessel_info$Navire

cat("Selected vessel:", vessel_name, "\n")
cat("Input file:", input_file, "\n")
if (!file.exists(input_file)) stop("Vessel file not found: ", input_file)

# ---- LOADING VESSEL DATA ----
cat("Loading vessel data...\n")
system.time({
  if (grepl("\\.qs$", input_file, ignore.case = TRUE)) {
    dt_vessel <- qs::qread(input_file, as.data.table = TRUE)
    cat("QS format detected and loaded\n")
  } else if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
    dt_vessel <- readRDS(input_file)
    setDT(dt_vessel)
    cat("RDS format detected and loaded\n")
  } else {
    stop("Unrecognised file format: ", input_file)
  }
})

if ("ssvid" %in% names(dt_vessel)) dt_vessel[, ssvid := as.character(ssvid)]

# Safe normalisation of the coordinates
if ("Lon" %in% names(dt_vessel)) dt_vessel[, Lon := suppressWarnings(as.numeric(gsub(",", ".", Lon, fixed = TRUE)))]
if ("Lat" %in% names(dt_vessel)) dt_vessel[, Lat := suppressWarnings(as.numeric(gsub(",", ".", Lat, fixed = TRUE)))]

cat("✅ Data loaded:", nrow(dt_vessel), "observations\n")
cat("   Lon NA:", sum(is.na(dt_vessel$Lon)), " Lat NA:", sum(is.na(dt_vessel$Lat)), "\n")
if (nrow(dt_vessel) > 0) {
  cat("   Lon range:", paste(range(dt_vessel$Lon, na.rm = TRUE), collapse = " "),
      " Lat range:", paste(range(dt_vessel$Lat, na.rm = TRUE), collapse = " "), "\n")
}

# ---- READING VESSEL SPECIFICATIONS ----
spec_file <- file.path(Sys.getenv("DFP_REPO"), "config/ship_specs.yaml")
if (!file.exists(spec_file)) stop("ship_specs.yaml not found: ", spec_file)

spec_list  <- yaml::read_yaml(spec_file)$ship_specs
ship_specs <- rbindlist(spec_list, fill = TRUE)
ship_specs[, ssvid := as.character(ssvid)]
ship_specs <- ship_specs[, .SD[1], by = ssvid]
setnames(ship_specs, "service_speed_kn", "Service_speed")
setkey(ship_specs, ssvid)

if ("ssvid" %chin% names(dt_vessel)) {
  dt_vessel <- merge(dt_vessel, ship_specs[, .(ssvid, Service_speed, dredge_width_m, dredging_depth_m)],
                  by = "ssvid", all.x = TRUE)
  dt_vessel[, Service_speed := as.numeric(Service_speed)]
  dt_vessel[, `:=`(dredge_width_m = as.numeric(dredge_width_m),
                dredging_depth_m = as.numeric(dredging_depth_m))]
  setkey(dt_vessel, ssvid)
}

if (!"Service_speed" %chin% names(dt_vessel) || all(is.na(dt_vessel$Service_speed))) {
  warning("No service speed for ", vessel_name)
}

# ---- PHYSICAL SPEED FILTER ----
if ("Service_speed" %chin% names(dt_vessel)) {
  dt_vessel[, speed_limit := Service_speed * 1.15]
  n_before <- nrow(dt_vessel)
  dt_vessel   <- dt_vessel[is.na(speed_limit) | Speed <= speed_limit]
  n_after  <- nrow(dt_vessel)
  cat(sprintf("Physical speed filter: %d -> %d rows (%.2f %% kept)\n",
              n_before, n_after, 100 * n_after / n_before))
  dt_vessel[, speed_limit := NULL]
}

setorder(dt_vessel, Timestamp)
dt_vessel[, delta_t := c(NA_real_, diff(as.numeric(Timestamp)))]

# Unique (stable) index for debug/log
dt_vessel[, .row_id := .I]

cat("Applying geospatial filters...\n")

# ============================================================
# STATIC PRE-COMPUTATIONS (outside the loop) - SPATIAL TILING
# ============================================================
cat("\n--- Pre-computation of land relations (static) ---\n")

dt_vessel[, flag_on_land_static   := FALSE]
dt_vessel[, flag_near_land_static := FALSE]

land_local_3857 <- NULL
land_cache <- cache_new()

coord_ok <- is.finite(dt_vessel$Lon) & is.finite(dt_vessel$Lat) &
  dt_vessel$Lon >= -180 & dt_vessel$Lon <= 180 &
  dt_vessel$Lat >= -90  & dt_vessel$Lat <= 90
coord_ok[is.na(coord_ok)] <- FALSE
idx_ok <- which(coord_ok)

if (length(idx_ok) > 0) {
  dt_vessel[, tile_key := tile_assign(Lon, Lat, land_tile_deg)]
  tile_keys <- unique(dt_vessel$tile_key[idx_ok])
  n_tiles <- length(tile_keys)

  cat("  Unique tiles:", n_tiles, " (deg=", land_tile_deg, ")\n", sep = "")

  use_tiled_mode <- land_tile_enable && n_tiles <= land_tile_max

  if (use_tiled_mode) {
    # ============================================================
    # TILING MODE: loads land per tile
    # ============================================================
    cat("  >>> TILING mode enabled <<<\n")

    # Margin: near_coast for on_land/near_land + approx max segment length for crosses_land
    margin_deg <- max(near_coast_km / 111, default_speed_kn / 60) + 0.05

    total_polys <- 0
    total_polys_onland <- 0
    t_tile_start <- Sys.time()

    resolved_ds <- path.expand(prefer_land_datasource(land_mask_path))
    cat("  Land datasource (resolved):", resolved_ds, "\n")

    for (k in seq_along(tile_keys)) {
      tk <- tile_keys[k]
      tile_bbox <- tile_bbox_4326(tk, land_tile_deg, margin_deg = margin_deg)

      tile_land_m <- tryCatch(
        load_land_tile(land_mask_path, tile_bbox, buffer_m = land_buffer_m),
        error = function(e) { cat("    ERROR tile", tk, ":", e$message, "\n"); NULL }
      )

      # Cache for crosses_land in iterative loop (CONSERVATIVE: non-eroded)
      cache_set(land_cache, tk, tile_land_m)

      tile_idx <- idx_ok[dt_vessel$tile_key[idx_ok] == tk]
      n_polys <- if (!is.null(tile_land_m) && inherits(tile_land_m, "sf")) nrow(tile_land_m) else 0L
      total_polys <- total_polys + n_polys

      # Build on_land-only geometry (eroded mask) to reduce false positives
      tile_land_onland_m <- tile_land_m
      if (!is.null(tile_land_m) && inherits(tile_land_m, "sf") && nrow(tile_land_m) > 0 && land_onland_erode_m != 0) {
        tile_land_onland_m <- make_land_onland_only(tile_land_m, land_onland_erode_m)
      }
      n_polys_onland <- if (!is.null(tile_land_onland_m) && inherits(tile_land_onland_m, "sf")) nrow(tile_land_onland_m) else 0L
      total_polys_onland <- total_polys_onland + n_polys_onland

      cat(sprintf("  Tile %d/%d [%s]: %d pts, %d polys (on_land polys=%d)",
                  k, n_tiles, tk, length(tile_idx), n_polys, n_polys_onland))

      if (!is.null(tile_land_m) && inherits(tile_land_m, "sf") && nrow(tile_land_m) > 0 && length(tile_idx) > 0) {
        pts_sf <- sf::st_as_sf(dt_vessel[tile_idx, .(Lon, Lat)],
                               coords = c("Lon", "Lat"), crs = 4326, remove = TRUE)
        pts_m <- sf::st_transform(pts_sf, 3857)

        # on_land uses eroded geometry (on_land-only)
        on_land_vec <- chunk_any_intersects(pts_m, tile_land_onland_m, chunk = sf_chunk,
                                            label = paste0("on_land[", tk, "]"))
        dt_vessel[tile_idx, flag_on_land_static := on_land_vec]

        # near_land stays conservative (non-eroded)
        near_land_vec <- chunk_any_within_distance(pts_m, tile_land_m,
                                                   dist_m = near_coast_km * 1000,
                                                   chunk = sf_chunk,
                                                   label = paste0("near_land[", tk, "]"))
        dt_vessel[tile_idx, flag_near_land_static := near_land_vec]

        rm(pts_sf, pts_m)
        cat(" -> on_land:", sum(on_land_vec), " near_land:", sum(near_land_vec), "\n")
      } else {
        cat(" (no land or no pts)\n")
      }
    }

    cat("  Total polygons (sum over tiles):", total_polys, "\n")
    cat("  Total on_land-only polygons (sum over tiles):", total_polys_onland, "\n")
    cat("  Tiling duration:", round(as.numeric(difftime(Sys.time(), t_tile_start, units = "secs")), 1), "s\n")
    gc()

  } else {
    # ============================================================
    # FALLBACK: single BBOX
    # ============================================================
    if (!land_tile_enable) {
      cat("  TILING mode disabled by configuration\n")
    } else {
      cat("  Too many tiles (", n_tiles, " > ", land_tile_max, ") -> fallback single bbox\n", sep = "")
    }

    # Robust BBOX (quantiles) to avoid outliers
    lon_q <- as.numeric(quantile(dt_vessel$Lon[idx_ok], probs = c(0.001, 0.999), na.rm = TRUE))
    lat_q <- as.numeric(quantile(dt_vessel$Lat[idx_ok], probs = c(0.001, 0.999), na.rm = TRUE))

    dx <- near_coast_km / 111
    bbox_expanded <- sf::st_bbox(c(
      xmin = lon_q[1] - dx,
      ymin = lat_q[1] - dx,
      xmax = lon_q[2] + dx,
      ymax = lat_q[2] + dx
    ), crs = sf::st_crs(4326))

    cat("  bbox expanded (robust): xmin=", round(bbox_expanded[["xmin"]], 2),
        " ymin=", round(bbox_expanded[["ymin"]], 2),
        " xmax=", round(bbox_expanded[["xmax"]], 2),
        " ymax=", round(bbox_expanded[["ymax"]], 2), "\n", sep = "")

    # Read the local land with a WKT filter (load_land_tile, single bbox)
    land_local_3857 <- tryCatch(
      load_land_tile(land_mask_path, bbox_expanded, buffer_m = land_buffer_m),
      error = function(e) NULL
    )

    if (!is.null(land_local_3857) && inherits(land_local_3857, "sf") && nrow(land_local_3857) > 0) {
      pts_sf_all <- sf::st_as_sf(dt_vessel[idx_ok, .(Lon, Lat)],
                                coords = c("Lon", "Lat"), crs = 4326, remove = TRUE)
      pts_m_all  <- sf::st_transform(pts_sf_all, 3857)

      # on_land-only erosion
      land_onland_3857 <- make_land_onland_only(land_local_3857, land_onland_erode_m)
      cat("  Land polys:", nrow(land_local_3857), " | on_land polys:", nrow(land_onland_3857), "\n")

      cat("  Computing on_land (static)...\n")
      on_land_vec <- chunk_any_intersects(pts_m_all, land_onland_3857, chunk = sf_chunk, label = "on_land")
      dt_vessel[idx_ok, flag_on_land_static := on_land_vec]

      cat("  Computing near_land (static) via st_is_within_distance...\n")
      near_land_vec <- chunk_any_within_distance(pts_m_all, land_local_3857,
                                                 dist_m = near_coast_km * 1000,
                                                 chunk = sf_chunk,
                                                 label = "near_land")
      dt_vessel[idx_ok, flag_near_land_static := near_land_vec]

      rm(pts_sf_all, pts_m_all)
      gc()
    } else {
      cat("  Local land unavailable -> land filters disabled\n")
    }
  }
} else {
  cat("  No point with valid coordinates -> land filters disabled\n")
}

cat("---\n\n")

n_before_geo <- nrow(dt_vessel)

# ============================================================
# ITERATIVE LOOP WITH EXACT LOGGING
# ============================================================
iteration <- 0
max_iterations <- 3
total_removed <- 0
removed_log <- data.table()

decomp_one_label <- NULL

repeat {
  iteration <- iteration + 1
  cat(sprintf("  Iteration %d...\n", iteration))

  dt_vessel <- flag_geospatial_anomalies_instrumented(
    dt_vessel,
    land_local_m = land_local_3857,
    default_speed_kn = default_speed_kn,
    spike_dist_min_nm = spike_dist_min_nm,
    spike_bridge_max_nm = spike_bridge_max_nm,
    land_cache = land_cache,
    tile_key_col = if ("tile_key" %in% names(dt_vessel)) "tile_key" else NULL,
    min_seg_dt_sec = min_seg_dt_sec,
    min_seg_dist_m = min_seg_dist_m
  )

  n_flagged <- sum(dt_vessel$geo_flag, na.rm = TRUE)

  if (n_flagged == 0) {
    cat("  No new point to remove\n")
    dt_vessel[, c("geo_flag", "flag_on_land", "flag_crosses_land", "flag_speed_jump", "flag_long_jump", "flag_spike") := NULL]
    break
  }

  # LOGGING: points to remove before removal
  removed_iter <- dt_vessel[geo_flag == TRUE, .(
    .row_id, Timestamp, Lat, Lon, Speed,
    flag_on_land, flag_crosses_land, flag_speed_jump, flag_long_jump, flag_spike
  )]
  removed_iter[, iteration := iteration]
  removed_log <- rbind(removed_log, removed_iter, fill = TRUE)

  if (iteration >= max_iterations) {
    cat(sprintf("  Maximum number of iterations reached (%d), %d remaining points flagged\n",
                max_iterations, n_flagged))
    dt_vessel <- dt_vessel[geo_flag == FALSE | is.na(geo_flag)]
    dt_vessel[, c("geo_flag", "flag_on_land", "flag_crosses_land", "flag_speed_jump", "flag_long_jump", "flag_spike") := NULL]
    total_removed <- total_removed + n_flagged
    break
  }

  dt_vessel <- dt_vessel[geo_flag == FALSE | is.na(geo_flag)]
  dt_vessel[, c("geo_flag", "flag_on_land", "flag_crosses_land", "flag_speed_jump", "flag_long_jump", "flag_spike") := NULL]
  total_removed <- total_removed + n_flagged

  cat(sprintf("  -> %d points removed (cumulative total: %d)\n", n_flagged, total_removed))

  # Recompute delta_t
  setorder(dt_vessel, Timestamp)
  dt_vessel[, delta_t := c(NA_real_, diff(as.numeric(Timestamp)))]
}

dt_vessel[, .row_id := NULL]
if ("tile_key" %in% names(dt_vessel)) dt_vessel[, tile_key := NULL]

n_after_geo <- nrow(dt_vessel)
cat(sprintf("Geospatial filters: %d rows removed in %d iteration(s) (%.2f %%)\n",
            n_before_geo - n_after_geo, iteration,
            if (n_before_geo > 0) 100 * (n_before_geo - n_after_geo) / n_before_geo else 0))

# ============================================================
# EXACT BREAKDOWN OF THE REMOVALS
# ============================================================
cat("\n============================================================\n")
cat("EXACT BREAKDOWN OF THE REMOVALS\n")
cat("============================================================\n")

if (nrow(removed_log) > 0) {
  cat("\nBy flag (multi-label, a point can tick several):\n")
  cat(sprintf("  on_land:       %8d\n", sum(removed_log$flag_on_land)))
  cat(sprintf("  crosses_land:  %8d\n", sum(removed_log$flag_crosses_land)))
  cat(sprintf("  speed_jump:    %8d\n", sum(removed_log$flag_speed_jump)))
  cat(sprintf("  long_jump:     %8d\n", sum(removed_log$flag_long_jump)))
  cat(sprintf("  spike:         %8d\n", sum(removed_log$flag_spike)))

  removed_log[, reason := fcase(
    flag_on_land, "on_land",
    flag_crosses_land, "crosses_land",
    flag_speed_jump, "speed_jump",
    flag_long_jump, "long_jump",
    flag_spike, "spike",
    default = "other"
  )]

  decomp_one_label <- removed_log[, .(n = .N, pct = round(100 * .N / nrow(removed_log), 2)), by = reason][order(-n)]

  cat("\nBy category (single label with priority):\n")
  print(decomp_one_label)

  cat("\nBy iteration:\n")
  iter_summary <- removed_log[, .(n = .N), by = iteration][order(iteration)]
  print(iter_summary)

} else {
  cat("No point removed by the geospatial filters.\n")
}
cat("============================================================\n\n")

# ---- ISOLATION FOREST ----
setorder(dt_vessel, Timestamp)
dt_vessel[, delta_t       := c(NA_real_, diff(as.numeric(Timestamp)))]
dt_vessel[, Course_change := c(NA_real_, abs(diff(Course)))]
dt_vessel[Course_change > 180, Course_change := 360 - Course_change]
dt_vessel[, Accel         := c(NA_real_, diff(Speed))]

config_file <- file.path(Sys.getenv("DFP_REPO"), "config/outlier_config_V6.yaml")
if (file.exists(config_file)) {
  config <- yaml.load_file(config_file)
  outlier_config <- list(
    contamination_rate = config$isolation_forest$contamination_rate,
    if_sample_size     = config$isolation_forest$sample_size,
    if_num_trees       = config$isolation_forest$num_trees,
    memory_conservative = TRUE
  )
  cat("Configuration loaded:", config_file, "\n")
} else {
  outlier_config <- list(
    contamination_rate = 0.02,
    if_sample_size     = 256,
    if_num_trees       = 25,
    memory_conservative = TRUE
  )
  cat("Default configuration applied\n")
}

cat("IF parameters: contamination =", outlier_config$contamination_rate,
    "| trees =", outlier_config$if_num_trees, "\n")

detect_outliers_IF_optimized <- function(dt, config) {
  cat("  Starting Isolation Forest...\n")

  features <- c("Lat", "Lon")
  if ("delta_t" %in% names(dt)) features <- c(features, "delta_t")
  if ("Speed" %in% names(dt))   features <- c(features, "Speed")

  dt_features <- dt[, ..features]
  dt_features <- dt_features[complete.cases(dt_features)]

  if (nrow(dt_features) < 50) {
    cat("  Too few data (", nrow(dt_features), ") - no IF\n")
    return(rep(FALSE, nrow(dt)))
  }

  if_model <- isolationForest$new(
    sample_size = min(config$if_sample_size, nrow(dt_features)),
    num_trees   = config$if_num_trees,
    seed        = 42
  )

  if_model$fit(dt_features)
  scores  <- if_model$predict(dt_features)
  thr     <- quantile(scores$anomaly_score, 1 - config$contamination_rate)
  outliers <- scores$anomaly_score > thr

  result <- rep(FALSE, nrow(dt))
  if (nrow(dt_features) == nrow(dt)) {
    result <- outliers
  } else {
    complete_idx <- which(complete.cases(dt[, ..features]))
    result[complete_idx] <- outliers
  }

  cat("  IF finished:", sum(result), "outliers out of", nrow(dt), "points\n")
  result
}

cat("Starting Isolation Forest processing...\n")
system.time({
  dt_vessel[, outlier_IF := detect_outliers_IF_optimized(dt_vessel, outlier_config)]
})

# ---- CREATION OF ADDITIONAL COLUMNS ----
cat("Creating the is_stop column...\n")
setorder(dt_vessel, Timestamp)
dt_vessel[, is_stop := (Speed < 1) | (delta_t > 300)]
dt_vessel[is.na(is_stop), is_stop := FALSE]
cat("is_stop column created:", sum(dt_vessel$is_stop), "stops detected\n")

dt_vessel[, Annee := year(Timestamp)]
cat("Additional columns created (Annee, Course_change, Accel)\n")

# ---- STATISTICS ----
n_outliers   <- sum(dt_vessel$outlier_IF)
outlier_rate <- round(n_outliers / nrow(dt_vessel) * 100, 2)
n_stops      <- sum(dt_vessel$is_stop)
stop_rate    <- round(n_stops / nrow(dt_vessel) * 100, 2)

cat("\nVESSEL RESULTS:", vessel_name, "\n")
cat("  Total observations:", nrow(dt_vessel), "\n")
cat("  Outliers detected:", n_outliers, "(", outlier_rate, "%)\n")
cat("  Stops detected:", n_stops, "(", stop_rate, "%)\n")

# ---- SAVING ----
output_dir <- file.path(Sys.getenv("DFP_DATA"), paste0("ais_results_", Sys.getenv("SLURM_ARRAY_JOB_ID")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- gsub("[^A-Za-z0-9_-]", "_", vessel_name)
output_file <- file.path(output_dir, sprintf("%02d_%s_%s_clean.rds", task_id, split_job_id, safe_name))

cat("Saved:", output_file, "\n")
system.time({
  saveRDS(dt_vessel, output_file, compress = "xz")
})

if (nrow(removed_log) > 0) {
  removed_log_file <- file.path(output_dir, sprintf("%02d_%s_%s_removed_log.rds", task_id, split_job_id, safe_name))
  saveRDS(removed_log, removed_log_file)
  cat("Removal log:", removed_log_file, "\n")

  if (!is.null(decomp_one_label)) {
    decomp_file <- file.path(output_dir, sprintf("%02d_%s_%s_decomp.csv", task_id, split_job_id, safe_name))
    fwrite(decomp_one_label, decomp_file)
    cat("Decomposition:", decomp_file, "\n")
  }
}

# ---- CLEAN-UP ----
rm(dt_vessel)
gc()

cat("\nVessel", vessel_name, "processed successfully:", format(Sys.time()), "\n")
cat("Result:", output_file, "\n")
