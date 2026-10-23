#!/usr/bin/env Rscript
# ────────────────────────────────────────────────────────────────────────────────
# STEP-5  ─  TILE WORKER (pipeline V6)
# Consistent with the Methods:
#  - Stable 1 km grid (EPSG:6933), identical for the later steps
#  - Temporal ordering before linearisation (reliable distances)
#  - p_d = 1 m (TSHD)
#  - pings without lithology (or >10 km away) EXCLUDED from p_l (but not from SAR)
#  - p_l aggregated per cell as the simple mean of the pings with lithology
# ────────────────────────────────────────────────────────────────────────────────

## 0.  Libraries ----------------------------------------------------------------
.libPaths(c("~/R/library", .libPaths()))

# Required packages (arrow optional - RDS fallback)
pkgs_required <- c("sf", "dplyr", "data.table", "tidyr", "rlang", "tidyselect", "lubridate", "yaml")
pkgs_optional <- c("arrow")

safe_library <- function(pkg) {
  tryCatch(
    {
      library(pkg, character.only = TRUE)
      cat("✅", pkg, "OK\n")
      TRUE
    },
    error = function(e) {
      stop("❌ Missing package: ", pkg, "\nMessage: ", e$message)
    }
  )
}

optional_library <- function(pkg) {
  tryCatch(
    {
      library(pkg, character.only = TRUE)
      cat("✅", pkg, "OK\n")
      TRUE
    },
    error = function(e) {
      cat("⚠️", pkg, "not available - fallback RDS\n")
      FALSE
    }
  )
}

# Robust loading of sf
tryCatch(
  {
    library(sf)
    cat("✅ sf OK\n")
    options(sf_max_print = 20)
    sf::sf_use_s2(FALSE)
  },
  error = function(e) {
    stop("❌ Package 'sf' cannot be loaded.\nMessage: ", e$message)
  }
)

invisible(lapply(pkgs_required[pkgs_required != "sf"], safe_library))
HAS_ARROW <- optional_library("arrow")
cat("\n")

# Shared constants
this_file <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(f)) {
    return(normalizePath(f))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(normalizePath(sys.frame(1)$ofile))
  }
  stop("Cannot locate running script")
}
script_dir <- dirname(this_file())
source(file.path(script_dir, "constants.R"))

# Safety fallback for an older constants.R (without DEEP_HORIZON)
if (!exists("DEEP_HORIZON")) {
  DEEP_HORIZON <- SURF_HORIZON  # cap = same value as the surface layer (5 cm)
  cat("⚠️  DEEP_HORIZON missing from constants.R — fallback to", DEEP_HORIZON, "m\n")
}

# === Stable cell grid (EPSG:6933) ===
snap_cells2 <- function(DT, xcol = "X", ycol = "Y") {
  x <- as.numeric(DT[[xcol]])
  y <- as.numeric(DT[[ycol]])
  col <- floor((x - WORLD_XMIN) / CELL_SIZE_M)
  row <- floor((y - WORLD_YMIN) / CELL_SIZE_M)
  cx <- WORLD_XMIN + col * CELL_SIZE_M + CELL_SIZE_M / 2
  cy <- WORLD_YMIN + row * CELL_SIZE_M + CELL_SIZE_M / 2
  DT[, `:=`(
    col     = as.integer(col),
    row     = as.integer(row),
    x       = as.numeric(cx), # cell centre
    y       = as.numeric(cy),
    grid_id = as.integer(row * GRID_COLS + col + 1L)
  )]
  DT
}

mean_valid_endpoint_pl <- function(pl_start, pl_end) {
  n_valid <- as.integer(!is.na(pl_start)) + as.integer(!is.na(pl_end))
  out <- rep(NA_real_, length(pl_start))
  keep <- n_valid > 0L

  if (any(keep)) {
    pl_start0 <- pl_start
    pl_end0 <- pl_end
    pl_start0[is.na(pl_start0)] <- 0
    pl_end0[is.na(pl_end0)] <- 0
    out[keep] <- (pl_start0[keep] + pl_end0[keep]) / n_valid[keep]
  }

  out
}

empty_carbon_segments_sf <- function(target_crs) {
  st_sf(
    data.table(sub_seg_id = character(), pl_seg = numeric()),
    geometry = st_sfc(crs = target_crs)
  )
}

build_carbon_segments_sf <- function(dredge, target_crs) {
  if (!nrow(dredge)) {
    return(empty_carbon_segments_sf(target_crs))
  }

  seg_dt <- copy(dredge)
  setorderv(seg_dt, c("sub_seg_id", "Timestamp"))
  seg_dt[, valid_pl := fifelse(use_for_carbon == TRUE & is.finite(pl_base), pl_base, NA_real_)]

  carbon_dt <- seg_dt[
    ,
    {
      next_x <- shift(X, type = "lead")
      next_y <- shift(Y, type = "lead")
      next_pl <- shift(valid_pl, type = "lead")
      .(
        x0 = X,
        y0 = Y,
        x1 = next_x,
        y1 = next_y,
        pl_seg = mean_valid_endpoint_pl(valid_pl, next_pl)
      )
    },
    by = sub_seg_id
  ]

  carbon_dt <- carbon_dt[
    is.finite(x0) & is.finite(y0) &
      is.finite(x1) & is.finite(y1) &
      is.finite(pl_seg) &
      !(x0 == x1 & y0 == y1)
  ]

  if (!nrow(carbon_dt)) {
    return(empty_carbon_segments_sf(target_crs))
  }

  carbon_geom <- st_sfc(
    Map(
      function(x0, y0, x1, y1) {
        st_linestring(matrix(c(x0, y0, x1, y1), ncol = 2, byrow = TRUE))
      },
      carbon_dt$x0,
      carbon_dt$y0,
      carbon_dt$x1,
      carbon_dt$y1
    ),
    crs = target_crs
  )

  st_sf(carbon_dt[, .(sub_seg_id, pl_seg)], geometry = carbon_geom)
}

accumulate_segment_pl_to_grid <- function(carbon_segments_sf, grid1km) {
  if (!nrow(carbon_segments_sf)) {
    return(data.table(grid_id = integer(), n_with_pl = integer(), pl_sum = numeric()))
  }

  # Batch intersection: a single GEOS call (internal R-tree index)
  inter <- suppressWarnings(st_intersection(carbon_segments_sf, grid1km))
  if (!nrow(inter)) {
    return(data.table(grid_id = integer(), n_with_pl = integer(), pl_sum = numeric()))
  }

  inter_dt <- as.data.table(inter)
  inter_dt[, seg_len := as.numeric(st_length(geometry))]
  inter_dt <- inter_dt[seg_len > 0]

  if (!nrow(inter_dt)) {
    return(data.table(grid_id = integer(), n_with_pl = integer(), pl_sum = numeric()))
  }

  # Deduplicate: a segment counts once per cell
  inter_dt <- unique(inter_dt, by = c("sub_seg_id", "grid_id"))

  # Aggregate: number of segments and sum of pl_seg per cell
  pl_cell <- inter_dt[, .(
    n_with_pl = .N,
    pl_sum    = sum(pl_seg, na.rm = TRUE)
  ), by = grid_id]

  setkey(pl_cell, grid_id)
  pl_cell
}

## 1.  Call parameters ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(1, 6))) {
  stop(
    "Usage:\n  Rscript step5_tile_worker_batch.R <tile_id>\n  ",
    "Rscript step5_tile_worker_batch.R <tile_id> <sub_id> <xmin> <xmax> <ymin> <ymax>"
  )
}

tile_id <- as.integer(args[1])
if (is.na(tile_id) || tile_id < 1) stop("❌ Invalid tile ID: ", args[1])

is_subtile <- length(args) == 6
if (is_subtile) {
  sub_id <- as.integer(args[2])
  xmin <- as.numeric(args[3])
  xmax <- as.numeric(args[4])
  ymin <- as.numeric(args[5])
  ymax <- as.numeric(args[6])
  if (anyNA(c(sub_id, xmin, xmax, ymin, ymax))) {
    stop("❌ Invalid sub-tile arguments.")
  }
  cat(sprintf(
    "🔧 Sub-tile %d_%d — bbox: (%.0f,%.0f)-(%.0f,%.0f)\n",
    tile_id, sub_id, xmin, ymin, xmax, ymax
  ))
} else {
  cat("🔧 Processing the whole tile:", tile_id, "\n")
}

# YAML parameters
param_yaml <- file.path(Sys.getenv("DFP_REPO"), "config/fi_parameters_with_freshness.yaml")
params_raw <- yaml::read_yaml(param_yaml)

get_scenario <- function(name) {
  s <- params_raw$scenarios[[name]]
  if (!is.null(s$inherit)) {
    parent <- get_scenario(s$inherit)
    s$inherit <- NULL
    modifyList(parent, s)
  } else {
    s
  }
}
par <- get_scenario("default")

# Table of regional k (used later)
k_table <- data.frame(
  longhurst_pr = names(par$k_fast),
  k_fast = unlist(par$k_fast) *
    ifelse(is.null(par$k_fast_multiplier), 1, par$k_fast_multiplier)
)

# Required parameters are checked, as in the merge: a missing key stops the
# run rather than falling back silently to an invented value.
.required <- c("alpha_dep", "fast_fraction", "slow_k", "preservation_factor", "k_fast")
.miss <- setdiff(.required, names(par))
if (length(.miss))
  stop("Missing parameters in ", param_yaml, " (scenario 'default'): ",
       paste(.miss, collapse = ", "),
       ". No fallback value is applied: a missing parameter is never invented.",
       "")

alpha_dep <- par$alpha_dep
fast_frac <- par$fast_fraction
slow_k <- par$slow_k
preserv_fact <- par$preservation_factor

cat(
  "🔧  Parameters (for information): alpha_dep=", alpha_dep,
  " fast_fraction=", fast_frac, " slow_k=", slow_k,
  " preservation_factor=", preserv_fact, "\n"
)

Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)
data.table::setDTthreads(1)
cat("✅ Threads data.table = 1\n")

## 2bis. I/O PARAMETERS ---------------------------------------------------------
# No hard-coded path writes into $DFP_DATA/output_V6: a replay must never
# overwrite the published outputs.
PUBLISHED_DIR <- path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6"))
PIPE_TILES   <- path.expand(Sys.getenv("PIPE_TILES",  file.path(PUBLISHED_DIR, "tiles_1000km.gpkg")))
PIPE_IN_DIR  <- path.expand(Sys.getenv("PIPE_IN_DIR", PUBLISHED_DIR))
PIPE_OUT_DIR <- Sys.getenv("PIPE_OUT_DIR", "")
if (!nzchar(PIPE_OUT_DIR)) {
  stop("PIPE_OUT_DIR is not set. Declare it explicitly: this script ",
       "writes six files, and a silent default would overwrite the published outputs.")
}
PIPE_OUT_DIR <- path.expand(PIPE_OUT_DIR)
if (identical(normalizePath(PIPE_OUT_DIR, mustWork = FALSE),
              normalizePath(PUBLISHED_DIR, mustWork = FALSE)) &&
    !nzchar(Sys.getenv("PIPE_ALLOW_PUBLISHED"))) {
  stop("Refusing to write into the PUBLISHED outputs (", PUBLISHED_DIR, "). ",
       "Choose another PIPE_OUT_DIR, or PIPE_ALLOW_PUBLISHED=1 to force.")
}
dir.create(PIPE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Dredging label. The two are NOT tested the same way:
#   Dragage_flag    : integer 0/1, may be NA
#   behavior_smooth: string, value "dredging"
DREDGE_LABEL <- Sys.getenv("DREDGE_LABEL", "flag")
if (!DREDGE_LABEL %in% c("flag", "smooth")) {
  stop("DREDGE_LABEL must be \"flag\" or \"smooth\", got: ", DREDGE_LABEL)
}
cat("📌 PIPE_IN_DIR  =", PIPE_IN_DIR,  "\n")
cat("📌 PIPE_OUT_DIR =", PIPE_OUT_DIR, "\n")
cat("📌 DREDGE_LABEL =", DREDGE_LABEL, "\n")

## 2.  Work area ----------------------------------------------------------------
tiles_file <- PIPE_TILES
tiles <- st_read(tiles_file, quiet = TRUE)
target_crs <- st_crs(tiles)

if (is_subtile) {
  tile_geom <- st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = target_crs))
  tile_buffered <- tile_geom
} else {
  tile_bb <- tiles[tiles$tile_id == tile_id, ]
  if (nrow(tile_bb) == 0) stop("❌ Tile ", tile_id, " missing from ", tiles_file)
  tile_geom <- st_geometry(tile_bb)
  tile_buffered <- st_buffer(tile_geom, TILE_BUFFER_M)
}

## 3.  Reading AIS + lithology --------------------------------------------------
lithology_files <- list.files(PIPE_IN_DIR,
  pattern = "AIS_with_lithology_.*\\.rds$",
  full.names = TRUE
)
if (!length(lithology_files)) {
  stop("❌ No AIS_with_lithology_* file found.")
}

# --- input choice: declared if possible, traced otherwise -------------------
# `rule` is "mtime" or "name": each step keeps its own rule; they are NOT
# harmonised.
choose_input <- function(dir, pat, explicit, what_kind, rule) {
  if (nzchar(explicit)) {
    p <- path.expand(explicit)
    if (!file.exists(p)) stop("File ", what_kind, " declared but not found: ", p)
    cat(sprintf("📌 input %s DECLARED: %s\n", what_kind, basename(p)))
    return(p)
  }
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(f)) stop("No file ", what_kind, " (", pat, ") in ", dir)
  by_name   <- f[order(basename(f), method = "radix")][length(f)]
  by_mtime <- f[which.max(file.mtime(f))]
  selected    <- if (identical(rule, "mtime")) by_mtime else by_name
  cat(sprintf("⚠️ input %s NOT declared: inferred (rule=%s) among %d candidates in %s\n",
              what_kind, rule, length(f), dir))
  for (x in f) cat(sprintf("     %-52s mtime %s%s\n", basename(x),
                           format(file.mtime(x), "%Y-%m-%d %H:%M"),
                           if (identical(x, selected)) "   <== SELECTED" else ""))
  if (!identical(by_name, by_mtime)) {
    cat("🔴 THE TWO RULES DISAGREE — sort by name: ", basename(by_name),
        " | mtime: ", basename(by_mtime), "\n", sep = "")
    cat("   The choice is AMBIGUOUS. Declare the file explicitly.\n")
  }
  selected
}
litho_path <- choose_input(PIPE_IN_DIR, "AIS_with_lithology_.*\\.rds$",
                             Sys.getenv("PIPE_LITHO_FILE", ""), "step4", rule = "name")
dt_lithology <- readRDS(litho_path)
# Real stop marker, computed on the FULL
# stream before the Dragage_flag filter: knowing whether a stop occurred
# BETWEEN two dredging pings requires seeing the pings in between.
setorder(dt_lithology, ssvid, Timestamp)
dt_lithology[, stop_real_v2 := !is.na(stop_cluster) & stop_cluster > 0 & Speed < 1]
dt_lithology[, cum_stop_v2 := cumsum(stop_real_v2), by = ssvid]
dt_lithology[, stop_real_v2 := NULL]
drag_sel <- if (DREDGE_LABEL == "smooth") {
  !is.na(dt_lithology$behavior_smooth) & dt_lithology$behavior_smooth == "dredging"
} else {
  !is.na(dt_lithology$Dragage_flag) & dt_lithology$Dragage_flag == 1
}
cat(sprintf("📌 label %s: %s dredging pings out of %s\n", DREDGE_LABEL,
            format(sum(drag_sel), big.mark = " "), format(nrow(dt_lithology), big.mark = " ")))
dredge_raw <- dt_lithology[drag_sel & !is.na(Lon) & !is.na(Lat)]
if (!nrow(dredge_raw)) {
  cat("⚠️  No dredging ping — clean exit\n")
  quit("no")
}

dredge_sf <- st_as_sf(dredge_raw, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) |>
  st_transform(target_crs)

sel <- st_intersects(dredge_sf, tile_buffered, sparse = FALSE)[, 1]
dredge_sf <- dredge_sf[sel, ]
if (!nrow(dredge_sf)) {
  cat("⚠️  No ping in the area — empty SAR file written\n")
  empty <- data.table(
    grid_id = integer(), sum_dw = numeric(),
    sum_dw_pd = numeric(), sum_d = numeric(),
    n_with_pl = integer(), pl_sum = numeric()
  )
  out <- if (is_subtile) {
    sprintf(file.path(PIPE_OUT_DIR, "sar_%d_sub%d.parquet"), tile_id, sub_id)
  } else {
    sprintf(file.path(PIPE_OUT_DIR, "sar_%03d.parquet"), tile_id)
  }
  arrow::write_parquet(empty, out)
  quit("no")
}

coords <- st_coordinates(dredge_sf)
dredge_sf$X <- coords[, 1]
dredge_sf$Y <- coords[, 2]
dredge <- as.data.table(dredge_sf)
dredge[, geometry := NULL]

# Harmonise the vessel identifier
id_cols <- c("ssvid", "SSVID", "mmsi", "MMSI", "vessel_id", "VESSEL_ID", "Navire")
for (c in id_cols) {
  if (c %in% names(dredge)) {
    setnames(dredge, c, "ssvid", skip_absent = TRUE)
    break
  }
}

# Default values / clean-up
if (!"dredging_depth_m" %in% names(dredge)) dredge[, dredging_depth_m := 1.0]
if (!"dredge_width_m" %in% names(dredge))
  stop("dredge_width_m missing: the swept width cannot be invented.",
       "",
       "",
       "")
if (!"pl_base" %in% names(dredge)) dredge[, pl_base := NA_real_]

# --- Suction pipes IN SERVICE ---------------------------------------------------
# `dredge_width_m` is the DIAMETER of ONE pipe (convention of constants.R).
# The swept width depends on the number of pipes in service:
#     n_pipes_eff = n - (n - 1) * rho        rho = single-head time fraction
#     rho = 0  -> all installed pipes in service   (option C, DEFAULT)
#     rho = 1  -> a single pipe everywhere         (option B)
# The mass balance stays satisfied as long as rho <= 34 % (see CUT_DEPTH_M in constants.R).
# `n_suction_pipes` (present for all 10 vessels) is read here; the older
# `twin_pipes` key of ship_specs.yaml is not read by any script.
PIPE_MONO_FRACTION <- as.numeric(Sys.getenv("PIPE_MONO_FRACTION", "0.0"))
if (!is.finite(PIPE_MONO_FRACTION) || PIPE_MONO_FRACTION < 0 || PIPE_MONO_FRACTION > 1)
  stop("PIPE_MONO_FRACTION outside [0,1]: ", Sys.getenv("PIPE_MONO_FRACTION"))

.specs_path <- Sys.getenv("SHIP_SPECS", "")
if (!nzchar(.specs_path)) {
  .cand <- c(path.expand(file.path(Sys.getenv("DFP_REPO"), "config/ship_specs.yaml")),
             path.expand("~/ais-pipeline/pipeline_repro/config/ship_specs.yaml"),
             path.expand("~/ais-pipeline/dredging-pipeline-pub/config/ship_specs.yaml"))
  .hit <- .cand[file.exists(.cand)]
  .specs_path <- if (length(.hit)) .hit[1] else ""
}
if (!nzchar(.specs_path) || !file.exists(.specs_path))
  stop("ship_specs.yaml not found -- set SHIP_SPECS")
cat("\u2139\ufe0f  ship_specs:", .specs_path, "\n")

.sp <- yaml::read_yaml(.specs_path)$ship_specs
.pipes <- data.table::rbindlist(lapply(.sp, function(s) data.table::data.table(
  .sv = as.character(s$ssvid),
  diam_m = as.numeric(s$suction_pipe_diameter_mm) / 1000,
  n_pipes = if (is.null(s$n_suction_pipes)) NA_integer_ else as.integer(s$n_suction_pipes))))
if (anyNA(.pipes$n_pipes))
  stop("n_suction_pipes missing in ship_specs.yaml for: ",
       paste(.pipes$.sv[is.na(.pipes$n_pipes)], collapse = ", "))

dredge[, .sv := as.character(ssvid)]
dredge[.pipes, on = ".sv", `:=`(n_pipes = i.n_pipes, diam_m = i.diam_m)]  # by reference
.orphan <- unique(dredge$.sv[is.na(dredge$n_pipes)])
if (length(.orphan))
  stop("ssvid missing from ship_specs.yaml: ", paste(.orphan, collapse = ", "))
dredge[, n_pipes_eff := n_pipes - (n_pipes - 1) * PIPE_MONO_FRACTION]

# 🔴 DOUBLE COUNTING -- the trap to avoid. `dredge_width_m` comes from step2 and is
# FROZEN in the intermediate data, where it may follow an older convention
# (2 x the pipe diameter). Multiplying it by n_pipes_eff would double the width
# TWICE. The DIAMETER is therefore taken from the YAML, the single
# source of truth, and any gap with the inherited column is logged.
.gap <- if ("dredge_width_m" %in% names(dredge))
  unique(dredge[abs(dredge_width_m - diam_m) > 1e-9,
                .(ssvid, inherited = dredge_width_m, yaml_diam = diam_m)]) else NULL
if (!is.null(.gap) && nrow(.gap)) {
  cat("\u26a0\ufe0f  dredge_width_m inherited from step2 differs from the YAML diameter",
      "-- the YAML prevails:\n")
  print(.gap)
}
dredge[, .sv := NULL]
cat(sprintf("\u2139\ufe0f  pipes in use: rho = %.3f, n_pipes_eff from %.3f to %.3f\n",
            PIPE_MONO_FRACTION, min(dredge$n_pipes_eff), max(dredge$n_pipes_eff)))

# p_d = 1 m (TSHD), as in the Methods
dredge[, dredging_depth_m := 1.0]

# Snap to the 1 km grid and grid_id (always)
dredge <- snap_cells2(dredge, "X", "Y")

# pl_eff from pl_base (BEFORE use)
# dz2 capped at DEEP_HORIZON (max 5 cm)
dredge[, `:=`(
  dz1 = pmin(dredging_depth_m, SURF_HORIZON),
  dz2 = pmin(pmax(0, dredging_depth_m - SURF_HORIZON), DEEP_HORIZON)  # CAPPED
)]
dredge[, pl_eff := (pl_base * dz1 + FACTOR_DEEP * pl_base * dz2) / pmax(dredging_depth_m, 1e-6)]

# Flag the pings with lithology usable for carbon (≤10 km), based on pl_eff
# Exclude pings without valid lithology (has_lithology=FALSE)
if (!"dist_km" %in% names(dredge)) dredge[, dist_km := NA_real_]
if (!"has_lithology" %in% names(dredge)) dredge[, has_lithology := TRUE]  # fallback if missing
dredge[, use_for_carbon := has_lithology == TRUE & is.finite(pl_eff) & (is.na(dist_km) | dist_km <= 10)]

## 4.  ping_detail export (whole tile only)
if (!is_subtile) {
  dredge[, `:=`(
    dz1 = pmin(dredging_depth_m, SURF_HORIZON),
    dz2 = pmin(pmax(0, dredging_depth_m - SURF_HORIZON), DEEP_HORIZON)  # CAPPED
  )]
  dredge[, pl_eff := (pl_base * dz1 + FACTOR_DEEP * pl_base * dz2) / pmax(dredging_depth_m, 1e-6)]
  detail_cols <- c(
    "Timestamp", "Lon", "Lat", "X", "Y", "col", "row", "grid_id", "ssvid",
    "dredging_depth_m", "dredge_width_m", "pl_base", "pl_eff", "lithologie", "dist_km",
    "has_lithology", "use_for_carbon", "diam_m", "n_pipes", "n_pipes_eff"
  )
  arrow::write_parquet(dredge[, ..detail_cols],
    sprintf(file.path(PIPE_OUT_DIR, "ping_detail_%03d.parquet"), tile_id),
    compression = "zstd"
  )
  cat("✅ ping_detail exported\n")
}

## 5.  Conversion to LINES with TIME ORDERING ----------------------------------
cat("🔄 Conversion to lines (time order) ...\n")
to_posix <- function(x) {
  if (inherits(x, "POSIXt")) {
    x
  } else if (is.numeric(x)) {
    as.POSIXct(x, origin = "1970-01-01", tz = "UTC")
  } else {
    as.POSIXct(x, tz = "UTC")
  }
}
if (!"Timestamp" %in% names(dredge)) dredge[, Timestamp := NA]
dredge[, Timestamp := to_posix(Timestamp)]

# Daily segment per vessel
# Grouping by calendar day alone would join into
# a single polyline all the pings of the day, unloading transits included.
# Two break criteria: implied speed (a stretch at 12 kn is not dredging)
# AND time gap (a long reception gap gives a LOW implied speed
# and would pass the speed criterion alone).
# 600 s time threshold: published precedent for dredger AIS (Aquatic Living
# Resources, 2017), corroborated by controlled thinning of dense stretches.
# Can be overridden to MEASURE the effect of the break: a replay with
# RUPTURE_V_MAX_KN=Inf and RUPTURE_DT_MAX_S=Inf gives the same chain without breaks,
# hence the break factor alone. Without these variables, the published values
# are used: 3.5 kn and 600 s
# (published precedent, Aquatic Living Resources).
V_MAX_KN <- as.numeric(Sys.getenv("RUPTURE_V_MAX_KN", "3.5"))
DT_MAX_S <- as.numeric(Sys.getenv("RUPTURE_DT_MAX_S", "600"))
if (is.na(V_MAX_KN) || V_MAX_KN <= 0)
  stop("RUPTURE_V_MAX_KN invalid: ", Sys.getenv("RUPTURE_V_MAX_KN"))
if (is.na(DT_MAX_S) || DT_MAX_S <= 0)
  stop("RUPTURE_DT_MAX_S invalid: ", Sys.getenv("RUPTURE_DT_MAX_S"))
cat(sprintf("📌 break thresholds: %.4g kn, %.4g s
", V_MAX_KN, DT_MAX_S))
stopifnot(requireNamespace("geosphere", quietly = TRUE))
setorder(dredge, ssvid, Timestamp)
dredge[, .dt_s := c(NA_real_, diff(as.numeric(Timestamp))), by = ssvid]
dredge[, .d_m := c(NA_real_, geosphere::distGeo(
          cbind(Lon, Lat)[-.N, , drop = FALSE],
          cbind(Lon, Lat)[-1L, , drop = FALSE])), by = ssvid]
dredge[, .v_kn := fifelse(.dt_s > 0, (.d_m / .dt_s) / 0.514444, NA_real_)]
dredge[, .stop_between := c(NA_real_, diff(cum_stop_v2)) > 0, by = ssvid]
dredge[, .rupture := is.na(.v_kn) | .v_kn > V_MAX_KN | .dt_s > DT_MAX_S |
                     (!is.na(.stop_between) & .stop_between)]
dredge[, sub_seg_id := paste0(ssvid, "_", cumsum(.rupture)), by = ssvid]
cat(sprintf("[THRESHOLD] %d stretches broken out of %d (%.1f %%)\n",
            sum(dredge$.rupture, na.rm = TRUE), nrow(dredge),
            100 * mean(dredge$.rupture, na.rm = TRUE)))
dredge[, c(".dt_s", ".d_m", ".v_kn", ".rupture", ".stop_between", "cum_stop_v2") := NULL]

# Recompute pl_eff defensively before building segment-level carbon support.
dredge[, `:=`(
  dz1 = pmin(dredging_depth_m, SURF_HORIZON),
  dz2 = pmin(pmax(0, dredging_depth_m - SURF_HORIZON), DEEP_HORIZON)  # CAPPED
)]
dredge[, pl_eff := (pl_base * dz1 + FACTOR_DEEP * pl_base * dz2) / pmax(dredging_depth_m, 1e-6)]

dredge_sf2 <- st_as_sf(dredge, coords = c("X", "Y"), crs = target_crs, remove = FALSE)

lines_sf <- dredge_sf2 |>
  arrange(sub_seg_id, Timestamp) |>
  group_by(sub_seg_id) |>
  filter(n() > 1) |>
  summarise(
    W_v = first(diam_m) * first(n_pipes_eff) * DRAGHEAD_WIDTH_FACTOR,
    p_d_i = CUT_DEPTH_M, # p_d = 1 m
    geometry = st_cast(st_combine(geometry), "MULTILINESTRING"),
    .groups = "drop"
  )

if (nrow(lines_sf) == 0) {
  cat("⚠️ No valid line — empty output\n")
  out <- if (is_subtile) {
    sprintf(file.path(PIPE_OUT_DIR, "sar_%d_sub%d.parquet"), tile_id, sub_id)
  } else {
    sprintf(file.path(PIPE_OUT_DIR, "sar_%03d.parquet"), tile_id)
  }
  arrow::write_parquet(
    data.table(
      grid_id = integer(), sum_dw = numeric(), sum_dw_pd = numeric(),
      sum_d = numeric(), n_with_pl = integer(), pl_sum = numeric()
    ),
    out
  )
  quit("no")
}
cat("✅ Lines created:", nrow(lines_sf), "\n")

## 6.  Local grid aligned on the world grid ------------------------------------
cat("🔲 Creating the aligned local grid…\n")
tile_bbox <- st_bbox(tile_buffered)

col_start <- as.integer((tile_bbox["xmin"] - WORLD_XMIN) %/% CELL_SIZE_M)
col_end <- as.integer((tile_bbox["xmax"] - WORLD_XMIN) %/% CELL_SIZE_M)
row_start <- as.integer((tile_bbox["ymin"] - WORLD_YMIN) %/% CELL_SIZE_M)
row_end <- as.integer((tile_bbox["ymax"] - WORLD_YMIN) %/% CELL_SIZE_M)

col_start <- max(0, col_start - 1)
col_end <- min(GRID_COLS - 1, col_end + 1)
row_start <- max(0, row_start - 1)
row_end <- min(GRID_ROWS - 1, row_end + 1)

local_grid <- st_make_grid(
  offset = c(
    WORLD_XMIN + col_start * CELL_SIZE_M,
    WORLD_YMIN + row_start * CELL_SIZE_M
  ),
  cellsize = CELL_SIZE_M,
  n = c(col_end - col_start + 1, row_end - row_start + 1),
  crs = target_crs,
  what = "polygons"
) |>
  st_sf() |>
  mutate(
    col     = col_start + ((dplyr::row_number() - 1) %% (col_end - col_start + 1)),
    row     = row_start + ((dplyr::row_number() - 1) %/% (col_end - col_start + 1)),
    grid_id = row * GRID_COLS + col + 1L
  )

grid1km <- st_intersection(local_grid, tile_buffered)
cat("✅ Local grid:", nrow(grid1km), "cells\n")

## 7.  Line×grid intersections (SAR) + p_l ------------------------------------
cat("🔄 Line/grid intersections…\n")
# The spatial index of grid1km is built once per batch of CHUNK lines, then
# queried CHUNK times, instead of being rebuilt for every line. Same formulas,
# same output; only the order of operations changes.
CHUNK <- 250L
parts <- vector("list", ceiling(nrow(lines_sf) / CHUNK))
for (b in seq_along(parts)) {
  lo  <- (b - 1L) * CHUNK + 1L
  hi  <- min(b * CHUNK, nrow(lines_sf))
  sub <- lines_sf[lo:hi, ]
  hit <- unique(unlist(sf::st_intersects(sub, grid1km)))
  if (!length(hit)) next
  inter <- suppressWarnings(sf::st_intersection(sub, grid1km[hit, ]))
  if (!nrow(inter)) next
  len <- as.numeric(sf::st_length(inter))
  parts[[b]] <- data.table(
    grid_id   = inter$grid_id,
    sum_dw    = len * inter$W_v,
    sum_dw_pd = len * inter$W_v * inter$p_d_i,
    sum_d     = len)
}
parts <- parts[!vapply(parts, is.null, TRUE)]
res <- if (length(parts)) {
  rbindlist(parts)[, .(sum_dw = sum(sum_dw), sum_dw_pd = sum(sum_dw_pd),
                       sum_d = sum(sum_d)), by = grid_id]
} else {
  data.table(grid_id = integer(), sum_dw = numeric(),
             sum_dw_pd = numeric(), sum_d = numeric())
}

# p_l aggregated PER PING, on pl_base.
#  - per ping: a ping lies in one cell only, so no support is carried
#    along the trajectory.
#  - on pl_base, not pl_eff: step5_merge recomputes p_l_eff itself from
#    p_l. Passing pl_eff would apply the depth weighting twice.
pl_cell <- dredge[use_for_carbon == TRUE & is.finite(pl_base),
                  .(n_with_pl = .N, pl_sum = sum(pl_base, na.rm = TRUE)),
                  by = grid_id]
setkey(pl_cell, grid_id)

# Join: keep the structure expected downstream
res <- merge(res, pl_cell, by = "grid_id", all.x = TRUE)
res[is.na(n_with_pl), `:=`(n_with_pl = 0L, pl_sum = 0)]

setkey(res, grid_id)
cat("✅ Intersections done:", nrow(res), "cells touched\n")

## 8.  Saving -------------------------------------------------------------------
out <- if (is_subtile) {
  sprintf(file.path(PIPE_OUT_DIR, "sar_%d_sub%d.parquet"), tile_id, sub_id)
} else {
  sprintf(file.path(PIPE_OUT_DIR, "sar_%03d.parquet"), tile_id)
}

if (requireNamespace("arrow", quietly = TRUE) &&
  packageVersion("arrow") >= numeric_version(PARQUET_VERSION_MIN)) {
  arrow::write_parquet(res, out)
  cat("✅ Results saved (Parquet):", basename(out), "\n")
} else {
  saveRDS(res, sub("\\.parquet$", ".rds", out))
  cat("⚠️  Arrow missing — saving as RDS\n")
}

cat("\n📊 TILE STATS ", if (is_subtile) sprintf("%d_%d", tile_id, sub_id) else tile_id, "\n",
  "   - Cells:", nrow(res), "\n",
  "   - Total distance:", format(sum(res$sum_d), scientific = FALSE), "m\n",
  "   - Mean SAR:", format(mean(res$sum_dw / CELL_AREA_M2, na.rm = TRUE), scientific = FALSE), "\n",
  sep = ""
)
cat("\n✅ Done\n")
