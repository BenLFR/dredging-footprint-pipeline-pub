#!/usr/bin/env Rscript
# ────────────────────────────────────────────────────────────────────────────────
# STEP-7  ─  BRIDGE STEP-6 CRI TO OCIM2-48L Jdredge FORCING
#
# Reads:  cri_final_*.parquet  (Step-6 output, EPSG:6933 1-km grid)
#         ocim_cache.mat       (portable OCIM grid cache with spatial enrichment)
# Writes: jdredge_ocim2_48l_YYYYMMDD_HHMMSS.mat  (Jdredge vector for co2model.m)
#
# Units chain:
#   Step-6 C_ri [gC m-2 yr-1] -> FRACAREA-normalized binning -> OCIM 2-deg flux
#   -> shelf-anchored kNN redistribution -> injection at the seabed depth -> umol C kg-1 yr-1
#
# Exports both Jdredge (primary) and Jtrawl (backward-compatible alias)
# plus amplified scenarios (10x, 100x) for sensitivity analysis.
#
# Implements:
#   M1: Ocean-fraction-aware normalization (FRACAREA)
#   M3: Shelf-anchored kNN redistribution with Gaussian distance weighting
#   M2: Mass-fix for double-precision residual
#   M4: injection depth (INJECTION_MODE)
#   Anti-artifact locks (6 hard assertions before injection)
#
# M4 -- WHY. Injecting the whole source of an OCIM column into its deepest wet layer (kbot)
# places, on the 2-deg grid, about 1/3 of the mass more than 1 km below its true seabed
# (e.g. Oman, Mumbai: slope columns at 1.3-2.8 km, dredging at ~20 m).
# HOW (inspired by John et al. 2020, Chem. Geol. 533:119403, AWESOME OCIM, p. 14 and 18, read at
# the source: ao.SEAFLOOR splits a seafloor source between layers according to the fraction of
# real ETOPO2 seafloor in each grid cell). Here, finer: each 1 km source cell gets its
# depth z_i = -median of GEBCO 2026 (15 arc-s) over 3x3 pixels (~1.4 km; modelling choice,
# no source); it enters the layer k that
# contains z_i (layer bottom = cumsum(dzt) of OCIM2-48L), capped at the kbot of the column that
# receives it. Mass moved by M3 keeps its depth. GEBCO >= 0 (on land) -> layer 1.
#   INJECTION_MODE=seafloor (default): as above.
#   INJECTION_MODE=bottom            : former behaviour (z_i = +Inf -> k = kbot), kept
#                                      to reproduce earlier runs identically.
# ────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(sf)
  library(R.matlab)
})

cat("=== STEP 7: EXPORT Jdredge FOR OCIM2-48L ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# ── M4: injection mode (see the header) ────────────────────────────────────────
INJECTION_MODE <- Sys.getenv("INJECTION_MODE", "seafloor")
if (!INJECTION_MODE %in% c("seafloor", "bottom")) stop("INJECTION_MODE invalid: ", INJECTION_MODE)
GEBCO_FILE <- path.expand(Sys.getenv("GEBCO_FILE",
  file.path(Sys.getenv("DFP_SCRATCH"), "configuration/gebco_2026/GEBCO_2026.nc")))  # GEBCO 2026 (downloaded 2026-09-25)
if (INJECTION_MODE == "seafloor") {
  if (!requireNamespace("terra", quietly = TRUE)) stop("INJECTION_MODE=seafloor requires the terra package")
  if (!file.exists(GEBCO_FILE)) stop("GEBCO_FILE not found: ", GEBCO_FILE)
}
cat(sprintf("Injection: %s%s\n", INJECTION_MODE,
            if (INJECTION_MODE == "seafloor") paste0("  (", GEBCO_FILE, ")") else "  (former behaviour: kbot)"))

# ── Redistribution parameters (saved to .mat metadata) ────────────────────────
# Can be set by environment variable, so that the sensitivity test of the
# redistribution kernel varies lambda without editing this script.
# Without any variable, the values are 500 / 100 / 16.
R_KM        <- as.numeric(Sys.getenv("REDIST_R_KM",      "500"))  # search radius [km] (~2.3 OCIM cells)
LAMBDA_KM   <- as.numeric(Sys.getenv("REDIST_LAMBDA_KM", "100"))  # Gaussian kernel scale [km]
K_MAX       <- as.integer(Sys.getenv("REDIST_K_MAX",      "16"))  # max neighbors
if (!is.finite(R_KM)      || R_KM      <= 0) stop("REDIST_R_KM invalid: ",      Sys.getenv("REDIST_R_KM"))
if (!is.finite(LAMBDA_KM) || LAMBDA_KM <= 0) stop("REDIST_LAMBDA_KM invalid: ", Sys.getenv("REDIST_LAMBDA_KM"))
if (is.na(K_MAX)          || K_MAX     <= 0) stop("REDIST_K_MAX invalid: ",     Sys.getenv("REDIST_K_MAX"))
cat(sprintf("Redistribution: R_KM=%.0f  LAMBDA_KM=%.0f  K_MAX=%d\n", R_KM, LAMBDA_KM, K_MAX))
ALPHA       <- 1      # ocean-fraction exponent
OCEANFRAC_FLOOR <- 0.05  # cells below this are "effectively land"

# ── Anti-artifact lock thresholds ─────────────────────────────────────────────
LOCK1_MASS_TOL       <- 1e-6    # max relative mass conservation error, residual BEFORE the fix; modelling choice, no source
LOCK3_SHELF_FLOOR    <- 0.30    # min fraction of flux on shelf (relaxed for 2-deg grid with 200 shelf cells)
LOCK4_HOTSPOT_CAP    <- 0.35    # max fraction in single cell (dredging is concentrated on 2-deg grid)
LOCK5_MEAN_DIST_CAP  <- 600     # max mean redistribution distance [km] (accommodates R_KM=500 + connectivity fallback)
LOCK6_DENSITY_RATIO  <- 1e4     # max/median density ratio
OFFSHORE_DELTA_KM    <- 50      # threshold for "offshore displacement"

# ── Section 0: Preflight ─────────────────────────────────────────────────────

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- file.path(Sys.getenv("DFP_REPO"), "pipeline")

constants_candidates <- c(
  file.path(script_dir, "constants.R"),
  file.path(script_dir, "pipeline", "constants.R"),
  file.path(getwd(), "constants.R"),
  file.path(getwd(), "pipeline", "constants.R")
)
constants_path <- constants_candidates[file.exists(constants_candidates)][1]
if (is.na(constants_path)) {
  stop(
    "Cannot find constants.R. Checked:\n",
    paste0(" - ", constants_candidates, collapse = "\n")
  )
}
source(constants_path)
cat("Using constants file:", constants_path, "\n")

OUT_DIR   <- path.expand(Sys.getenv("OUT_DIR", file.path(Sys.getenv("DFP_DATA"), "output_V6")))
OCIM_DIR  <- path.expand(Sys.getenv("OCIM_DIR", file.path(Sys.getenv("DFP_DATA"), "configuration/ocim")))

# Physical constants
M_C  <- 12.011   # molar mass of carbon [g/mol]
RHO  <- 1025     # seawater density [kg/m3]

# --- Load Step-6 CRI ---
# CRI_FILE: DECLARED input. It BYPASSES the search in OUT_DIR --
# otherwise an empty output directory would make a declared input fail.
# Without it: "the most recent file in the directory", which chooses
# SILENTLY when two labels coexist. The fallback now warns.
CRI_FILE <- Sys.getenv("CRI_FILE", "")
if (nzchar(CRI_FILE)) {
  cri_path <- path.expand(CRI_FILE)
  if (!file.exists(cri_path)) stop("CRI_FILE not found: ", cri_path)
  cat("Input DECLARED via CRI_FILE\n")
} else {
  cri_files <- list.files(OUT_DIR, pattern = "^cri_final_.*\\.parquet$", full.names = TRUE)
  if (length(cri_files) == 0) stop("No cri_final_*.parquet found in ", OUT_DIR)
  if (length(cri_files) > 1)
    warning("CRI_FILE not declared and ", length(cri_files),
            " candidate files: the most recent is taken.")
  cri_path <- cri_files[which.max(file.info(cri_files)$mtime)]
}
cat("Step-6 file:", basename(cri_path), "\n")

dt <- as.data.table(arrow::read_parquet(cri_path))

# Validate required columns
req_cols <- c("grid_id", "x", "y", "C_ri")
missing <- setdiff(req_cols, names(dt))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

stopifnot(all(dt$C_ri >= 0))
stopifnot(!any(duplicated(dt$grid_id)))
cat("  rows:", nrow(dt), " C_ri range: [", min(dt$C_ri), ",", max(dt$C_ri), "] gC/m2/yr\n")

# --- Load OCIM cache ---
ocim_cache_path <- file.path(OCIM_DIR, "ocim_cache.mat")
if (!file.exists(ocim_cache_path)) stop("ocim_cache.mat not found at ", ocim_cache_path)

cat("Loading OCIM cache:", ocim_cache_path, "\n")
ocim <- R.matlab::readMat(ocim_cache_path)

lon2d <- ocim$lon2d   # ni x nj
lat2d <- ocim$lat2d
M3d   <- ocim$M3d     # ni x nj x nk
DZT3d <- ocim$DZT3d
DXT3d <- ocim$DXT3d
DYT3d <- ocim$DYT3d
iocn  <- as.integer(ocim$iocn)  # m x 1
kbot  <- ocim$kbot    # ni x nj
ni    <- as.integer(ocim$ni)
nj    <- as.integer(ocim$nj)
nk    <- as.integer(ocim$nk)
m     <- length(iocn)

# New enrichment fields
depth2d          <- ocim$depth2d
oceanfrac2d      <- ocim$oceanfrac2d
shelf.mask2d     <- ocim$shelf.mask2d
dist.to.coast2d  <- ocim$dist.to.coast2d
shelf.component2d <- ocim$shelf.component2d
seed.shelf.i     <- ocim$seed.shelf.i
seed.shelf.j     <- ocim$seed.shelf.j

# Ocean connectivity fields (ice-9 flood-fill)
ocean.component2d   <- ocim$ocean.component2d
coastal.wet.mask2d  <- ocim$coastal.wet.mask2d
seed.ocean.i        <- ocim$seed.ocean.i
seed.ocean.j        <- ocim$seed.ocean.j

cat(sprintf("  OCIM grid: %d x %d x %d, m=%d ocean points\n", ni, nj, nk, m))
cat(sprintf("  Enrichment fields: depth2d, oceanfrac2d, shelf_mask2d, dist_to_coast2d,\n"))
cat(sprintf("    shelf_component2d, seed_shelf_ij,\n"))
cat(sprintf("    ocean_component2d, coastal_wet_mask2d, seed_ocean_ij loaded\n"))

# Diagnostics for ocean connectivity
n_ocomp <- length(unique(ocean.component2d[ocean.component2d > 0]))
n_coastal_wet <- sum(coastal.wet.mask2d > 0)
n_seeded_ocean <- sum(seed.ocean.i > 0)
cat(sprintf("  Ocean components: %d, Coastal wet cells: %d, Land cells with coastal seed: %d\n",
            n_ocomp, n_coastal_wet, n_seeded_ocean))

# ── Haversine helper (vectorized) ─────────────────────────────────────────────
haversine_km <- function(lon1, lat1, lon2, lat2) {
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(sqrt(a), 1))
}

# ── Section 1: Reproject EPSG:6933 -> WGS84 (lon/lat) ────────────────────────

cat("\n--- Section 1: Reproject to WGS84 ---\n")

pts_sf <- sf::st_as_sf(dt, coords = c("x", "y"), crs = 6933)
pts_wgs <- sf::st_transform(pts_sf, 4326)
coords_wgs <- sf::st_coordinates(pts_wgs)

dt[, lon := coords_wgs[, 1]]
dt[, lat := coords_wgs[, 2]]

cat(sprintf("  lon range: [%.2f, %.2f]  lat range: [%.2f, %.2f]\n",
            min(dt$lon), max(dt$lon), min(dt$lat), max(dt$lat)))

rm(pts_sf, pts_wgs, coords_wgs)

# ── Section 1b (M4): seabed depth and OCIM layer of each source cell ─────────

dzt_1d     <- as.numeric(ocim$dzt)            # thicknesses of the 48 layers of OCIM2-48L [m]
layer_base <- cumsum(dzt_1d)                  # depth of the bottom of each layer [m]
if (length(layer_base) != nk) stop("dzt: ", length(layer_base), " layers, nk = ", nk)

if (INJECTION_MODE == "seafloor") {
  cat("\n--- Section 1b: Seabed depth from GEBCO (M4) ---\n")
  gebco <- terra::rast(GEBCO_FILE)
  rc0 <- terra::rowColFromCell(gebco, terra::cellFromXY(gebco, cbind(dt$lon, dt$lat)))
  if (anyNA(rc0)) stop("source cell outside the GEBCO grid")
  offs <- expand.grid(dr = -1:1, dc = -1:1)   # 3x3 pixel window (~1.4 km), see the header
  vals <- vapply(seq_len(nrow(offs)), function(q) {
    rr <- pmin(pmax(rc0[, 1] + offs$dr[q], 1), terra::nrow(gebco))
    cc <- ((rc0[, 2] + offs$dc[q] - 1) %% terra::ncol(gebco)) + 1   # longitude wraps around
    terra::extract(gebco, terra::cellFromRowCol(gebco, rr, cc))[, 1]
  }, numeric(nrow(dt)))
  dt[, seabed_elev_m := apply(matrix(vals, nrow = nrow(dt)), 1, median)]
  if (anyNA(dt$seabed_elev_m)) stop("GEBCO: missing value for ", sum(is.na(dt$seabed_elev_m)), " cells")
  dt[, z_seabed_m := pmax(-seabed_elev_m, 0)]  # GEBCO >= 0 (on land, e.g. harbour basin) -> 0 m -> layer 1
  # layer that contains z: the smallest k such that bottom(k) >= z; z = 0 -> k = 1
  dt[, k_src := pmin(pmax(findInterval(z_seabed_m, c(0, layer_base), left.open = TRUE), 1L), nk)]
  w_m <- dt$C_ri
  cat(sprintf("  z_seabed (weighted by C_ri): median %.1f m; share on land according to GEBCO: %.2f%%\n",
              if (sum(w_m) > 0) dt$z_seabed_m[order(dt$z_seabed_m)][which(cumsum(w_m[order(dt$z_seabed_m)]) >= sum(w_m) / 2)[1]] else NA_real_,
              if (sum(w_m) > 0) 100 * sum(w_m[dt$seabed_elev_m >= 0]) / sum(w_m) else 0))
  rm(gebco, rc0, offs, vals, w_m)
} else {
  dt[, z_seabed_m := Inf]
  dt[, k_src := nk]                            # former behaviour: min(nk, kbot) = kbot
}

# ── Section 2: Map to OCIM (i,j) — FRACAREA-normalized binning (M1) ──────────

cat("\n--- Section 2: FRACAREA-normalized binning onto OCIM grid ---\n")

# OCIM lon/lat are cell centers; build edges
ocim_lon <- lon2d[1, ]  # nj values
ocim_lat <- lat2d[, 1]  # ni values

# Sort to ensure monotonic
lon_order <- order(ocim_lon)
lat_order <- order(ocim_lat)
ocim_lon_sorted <- ocim_lon[lon_order]
ocim_lat_sorted <- ocim_lat[lat_order]

# Build edges (midpoints between centers, extend at boundaries)
lon_edges <- c(ocim_lon_sorted[1] - diff(ocim_lon_sorted[1:2]) / 2,
               (ocim_lon_sorted[-length(ocim_lon_sorted)] + ocim_lon_sorted[-1]) / 2,
               ocim_lon_sorted[length(ocim_lon_sorted)] + diff(ocim_lon_sorted[(length(ocim_lon_sorted)-1):length(ocim_lon_sorted)]) / 2)

lat_edges <- c(ocim_lat_sorted[1] - diff(ocim_lat_sorted[1:2]) / 2,
               (ocim_lat_sorted[-length(ocim_lat_sorted)] + ocim_lat_sorted[-1]) / 2,
               ocim_lat_sorted[length(ocim_lat_sorted)] + diff(ocim_lat_sorted[(length(ocim_lat_sorted)-1):length(ocim_lat_sorted)]) / 2)

# Clamp edges
lon_edges <- pmax(pmin(lon_edges, 180), -180)
lat_edges <- pmax(pmin(lat_edges, 90), -90)

# Map each Step-6 point to sorted OCIM indices
dt[, j_sorted := findInterval(lon, lon_edges, rightmost.closed = TRUE)]
dt[, i_sorted := findInterval(lat, lat_edges, rightmost.closed = TRUE)]

# Clamp to valid range
dt[, j_sorted := pmin(pmax(j_sorted, 1L), nj)]
dt[, i_sorted := pmin(pmax(i_sorted, 1L), ni)]

# Map back from sorted indices to original OCIM indices
dt[, j_ocim := lon_order[j_sorted]]
dt[, i_ocim := lat_order[i_sorted]]

# Mass per 1-km cell [gC/yr]
dt[, mass_gC_yr := C_ri * CELL_AREA_M2]

# Uncertainty bounds
has_lower <- "C_ri_lower" %in% names(dt)
has_upper <- "C_ri_upper" %in% names(dt)
has_cons  <- "C_ri_conservative" %in% names(dt)
if (has_lower) dt[, mass_lower := C_ri_lower * CELL_AREA_M2]
if (has_upper) dt[, mass_upper := C_ri_upper * CELL_AREA_M2]
if (has_cons)  dt[, mass_cons  := C_ri_conservative * CELL_AREA_M2]

total_input_gC <- sum(dt$mass_gC_yr, na.rm = TRUE)
cat(sprintf("  Total input mass: %.4e gC/yr\n", total_input_gC))

# Aggregate mass per OCIM cell
agg_cols <- "mass_gC_yr"
if (has_lower) agg_cols <- c(agg_cols, "mass_lower")
if (has_upper) agg_cols <- c(agg_cols, "mass_upper")
if (has_cons)  agg_cols <- c(agg_cols, "mass_cons")

ocim_flux <- dt[, lapply(.SD, sum, na.rm = TRUE), by = .(i_ocim, j_ocim), .SDcols = agg_cols]
cat(sprintf("  Active OCIM cells after binning: %d\n", nrow(ocim_flux)))

# M4: same aggregation, keeping the layer of each source (vertical profile per column).
# ocim_flux (2D) remains the basis of M1-M3 and of the locks; profile follows the same mass in 3D.
profile <- dt[, lapply(.SD, sum, na.rm = TRUE), by = .(i_ocim, j_ocim, k_src), .SDcols = agg_cols]

# OCIM cell area
ocim_flux[, ocim_area_m2 := DXT3d[cbind(i_ocim, j_ocim, 1L)] * DYT3d[cbind(i_ocim, j_ocim, 1L)]]

# Ocean fraction from cache
ocim_flux[, oceanfrac := oceanfrac2d[cbind(i_ocim, j_ocim)]]

# Authoritative ocean mask: kbot > 0
ocim_flux[, is_ocean := kbot[cbind(i_ocim, j_ocim)] > 0]

# Cells with oceanfrac < OCEANFRAC_FLOOR are "effectively land" -> will be redistributed
ocim_flux[, effectively_land := (!is_ocean) | (oceanfrac < OCEANFRAC_FLOOR)]

# FRACAREA normalization (only for ocean cells that are NOT effectively land)
ocim_flux[effectively_land == FALSE,
          F_gC_m2_yr := mass_gC_yr / (ocim_area_m2 * oceanfrac)]

n_eff_land <- sum(ocim_flux$effectively_land)
flux_eff_land <- sum(ocim_flux$mass_gC_yr[ocim_flux$effectively_land], na.rm = TRUE)
cat(sprintf("  Effectively-land cells: %d (flux: %.4e gC/yr, %.1f%% of total)\n",
            n_eff_land, flux_eff_land, 100 * flux_eff_land / total_input_gC))

# ── Section 3: Shelf-anchored kNN redistribution (M3) ────────────────────────

cat("\n--- Section 3: Shelf-anchored kNN redistribution ---\n")
cat(sprintf("  Parameters: R_km=%d, lambda_km=%d, k_max=%d, alpha=%d\n",
            R_KM, LAMBDA_KM, K_MAX, ALPHA))

# Identify pools
land_pool_mask  <- ocim_flux$effectively_land
ocean_pool_mask <- !land_pool_mask

n_land_pool <- sum(land_pool_mask)
cat(sprintf("  Land pool: %d cells, Ocean pool: %d cells\n",
            n_land_pool, sum(ocean_pool_mask)))

# Redistribution counters (6-stage cascade with ocean connectivity)
n_relax_component       <- 0L  # kept for backward compat (= stage 1 hit)
n_relax_ocean_component <- 0L  # stage 2: same ocean comp, shelf, <= R_KM
n_relax_shelf           <- 0L  # stage 3: same ocean comp, coastal wet, <= R_KM
n_expand_radius         <- 0L  # stage 4: same ocean comp, any ocean, <= 2*R_KM
n_connectivity_fallback <- 0L  # stage 5: same ocean comp, coastal wet, nearest
n_absolute_fallback     <- 0L  # stage 6: forbidden hard lock (LOCK7)
n_fallback_nearest      <- 0L  # backward compat: sum of stages 4+5
n_cross_basin           <- 0L  # for Lock 2
n_cross_shelf_component <- 0L  # for Lock 2b (count of land cells with any cross-scomp transfer)
cross_scomp_mass_gC     <- 0   # flux-weighted mass sent to wrong shelf component [gC/yr]
redist_distances        <- numeric(0)  # for Lock 5 + metrics
trace_list              <- list()      # M4: (donor, recipients, weights), to move the vertical profile
land_ij <- ocim_flux[land_pool_mask, .(i_ocim, j_ocim)]  # M4: donor columns (before their removal)

if (n_land_pool > 0) {
  # Pre-build ocean-pool lookup tables
  ocean_idx  <- which(ocean_pool_mask)
  ocean_i    <- ocim_flux$i_ocim[ocean_idx]
  ocean_j    <- ocim_flux$j_ocim[ocean_idx]
  ocean_lon  <- lon2d[cbind(ocean_i, ocean_j)]
  ocean_lat  <- lat2d[cbind(ocean_i, ocean_j)]
  ocean_shelf <- as.logical(shelf.mask2d[cbind(ocean_i, ocean_j)])
  ocean_scomp <- shelf.component2d[cbind(ocean_i, ocean_j)]
  ocean_ocomp <- ocean.component2d[cbind(ocean_i, ocean_j)]
  ocean_coastal <- as.logical(coastal.wet.mask2d[cbind(ocean_i, ocean_j)])
  ocean_ofrac <- oceanfrac2d[cbind(ocean_i, ocean_j)]
  ocean_dist_coast <- dist.to.coast2d[cbind(ocean_i, ocean_j)]

  # Accumulator for mass added to ocean cells
  added_mass <- rep(0, length(ocean_idx))
  if (has_lower) added_mass_lower <- rep(0, length(ocean_idx))
  if (has_upper) added_mass_upper <- rep(0, length(ocean_idx))
  if (has_cons)  added_mass_cons  <- rep(0, length(ocean_idx))

  land_indices <- which(land_pool_mask)
  cat(sprintf("  Redistributing %d land-pool cells (6-stage cascade) ...\n", length(land_indices)))

  for (k in seq_along(land_indices)) {
    li <- land_indices[k]
    d0_i <- ocim_flux$i_ocim[li]
    d0_j <- ocim_flux$j_ocim[li]
    d0_lon <- lon2d[d0_i, d0_j]
    d0_lat <- lat2d[d0_i, d0_j]
    d0_mass <- ocim_flux$mass_gC_yr[li]

    if (d0_mass == 0) next

    # Get shelf seed for this land cell (shelf component anchor)
    seed_si <- as.integer(seed.shelf.i[d0_i, d0_j])
    seed_sj <- as.integer(seed.shelf.j[d0_i, d0_j])

    # Get coastal ocean seed for this land cell (ocean component anchor)
    seed_oi <- as.integer(seed.ocean.i[d0_i, d0_j])
    seed_oj <- as.integer(seed.ocean.j[d0_i, d0_j])

    # Target shelf component via shelf seed
    target_scomp <- 0
    if (seed_si > 0 && seed_sj > 0) {
      target_scomp <- shelf.component2d[seed_si, seed_sj]
    }

    # Target ocean component via coastal ocean seed
    target_ocomp <- 0
    if (seed_oi > 0 && seed_oj > 0) {
      target_ocomp <- ocean.component2d[seed_oi, seed_oj]
    }

    # Compute distances from d0 to all ocean-pool cells
    dists <- haversine_km(d0_lon, d0_lat, ocean_lon, ocean_lat)

    # ── 6-stage cascade ──────────────────────────────────────────────

    # Stage 1: same shelf component, shelf, <= R_KM (ideal)
    cand <- which(ocean_shelf & (ocean_scomp == target_scomp) & target_scomp > 0 & (dists <= R_KM))

    # Stage 2: same ocean component, shelf, <= R_KM (relax shelf component)
    if (length(cand) == 0) {
      cand <- which(ocean_shelf & (ocean_ocomp == target_ocomp) & target_ocomp > 0 & (dists <= R_KM))
      if (length(cand) > 0) {
        n_relax_ocean_component <- n_relax_ocean_component + 1L
        n_relax_component <- n_relax_component + 1L  # backward compat
      }
    }

    # Stage 3: same ocean component, coastal wet-points, <= R_KM (relax shelf preference)
    if (length(cand) == 0) {
      cand <- which(ocean_coastal & (ocean_ocomp == target_ocomp) & target_ocomp > 0 & (dists <= R_KM))
      if (length(cand) > 0) n_relax_shelf <- n_relax_shelf + 1L
    }

    # Stage 4: same ocean component, any ocean, <= 2*R_KM (expand radius)
    if (length(cand) == 0) {
      cand <- which((ocean_ocomp == target_ocomp) & target_ocomp > 0 & (dists <= 2 * R_KM))
      if (length(cand) > 0) {
        n_expand_radius <- n_expand_radius + 1L
        n_fallback_nearest <- n_fallback_nearest + 1L  # backward compat
      }
    }

    # Stage 5: same ocean component, coastal wet-points, nearest
    # (connectivity + coastal fallback — Adcroft approach)
    if (length(cand) == 0 && target_ocomp > 0) {
      coastal_in_basin <- which(ocean_coastal & (ocean_ocomp == target_ocomp))
      if (length(coastal_in_basin) > 0) {
        # Find anchor: nearest coastal wet-point in same basin
        anchor_idx <- coastal_in_basin[which.min(dists[coastal_in_basin])]
        anchor_lon <- ocean_lon[anchor_idx]
        anchor_lat <- ocean_lat[anchor_idx]

        # Distribute among coastal wet-points near the anchor
        dists_to_anchor <- haversine_km(anchor_lon, anchor_lat,
                                        ocean_lon[coastal_in_basin], ocean_lat[coastal_in_basin])
        # Use all coastal wet-points within R_KM of anchor, or top K_MAX nearest
        near_anchor <- which(dists_to_anchor <= R_KM)
        if (length(near_anchor) == 0) {
          near_anchor <- order(dists_to_anchor)[1:min(K_MAX, length(dists_to_anchor))]
        }
        cand <- coastal_in_basin[near_anchor]
        n_connectivity_fallback <- n_connectivity_fallback + 1L
        n_fallback_nearest <- n_fallback_nearest + 1L  # backward compat
      }
    }

    # Stage 6: absolute fallback — FORBIDDEN (hard lock, should never fire)
    if (length(cand) == 0) {
      n_absolute_fallback <- n_absolute_fallback + 1L
      cat(sprintf("  WARNING: LOCK7 — land cell (%d,%d) at (%.1f,%.1f) has no ocean target "
                  , d0_i, d0_j, d0_lon, d0_lat))
      cat(sprintf("(ocomp=%d, mass=%.4e gC/yr) — SKIPPING\n", target_ocomp, d0_mass))
      next
    }

    # Gaussian distance weighting * ocean-fraction
    w <- (pmax(ocean_ofrac[cand], 0.01)^ALPHA) *
         exp(-(dists[cand] / LAMBDA_KM)^2)

    # Take top k_max by weight
    if (length(cand) > K_MAX) {
      top_k <- order(w, decreasing = TRUE)[1:K_MAX]
      cand <- cand[top_k]
      w <- w[top_k]
    }

    # Guard: Gaussian weights can underflow to 0 when all candidates are >> LAMBDA_KM away
    # (typical in Stage-5 connectivity fallback). Use uniform weights as fallback.
    if (sum(w) == 0) {
      w[] <- 1 / length(w)
      cat(sprintf("    WARN: weight underflow for land cell (%d,%d) at (%.1f,%.1f) — uniform fallback\n",
                  d0_i, d0_j, d0_lon, d0_lat))
    }

    # Normalize
    w <- w / sum(w)

    # Cross-basin check: with ocean connectivity at every stage, this should be 0
    if (target_ocomp > 0) {
      recipient_ocomps <- ocean_ocomp[cand]
      if (any(recipient_ocomps > 0 & recipient_ocomps != target_ocomp)) {
        n_cross_basin <- n_cross_basin + 1L
      }
    }

    # Cross-shelf-component check (flux-weighted, more informative than cross-basin on 2-deg grids)
    if (target_scomp > 0) {
      recipient_scomps <- ocean_scomp[cand]
      wrong_scomp <- (recipient_scomps > 0) & (recipient_scomps != target_scomp)
      if (any(wrong_scomp)) {
        n_cross_shelf_component <- n_cross_shelf_component + 1L
        cross_scomp_mass_gC <- cross_scomp_mass_gC + sum(w[wrong_scomp]) * d0_mass
      }
    }

    # Track redistribution distances (flux-weighted)
    redist_distances <- c(redist_distances, sum(w * dists[cand]))
    trace_list[[length(trace_list) + 1L]] <- data.table(i_ocim = d0_i, j_ocim = d0_j,
                                                        rec_i = ocean_i[cand], rec_j = ocean_j[cand], w = w)

    # Distribute mass
    added_mass[cand] <- added_mass[cand] + w * d0_mass
    if (has_lower) added_mass_lower[cand] <- added_mass_lower[cand] + w * ocim_flux$mass_lower[li]
    if (has_upper) added_mass_upper[cand] <- added_mass_upper[cand] + w * ocim_flux$mass_upper[li]
    if (has_cons)  added_mass_cons[cand]  <- added_mass_cons[cand]  + w * ocim_flux$mass_cons[li]

    if (k %% 500 == 0) cat(sprintf("    %d / %d land cells processed\n", k, length(land_indices)))
  }

  # Apply redistributed mass to ocean pool
  ocim_flux[ocean_idx, mass_gC_yr := mass_gC_yr + added_mass]
  if (has_lower) ocim_flux[ocean_idx, mass_lower := mass_lower + added_mass_lower]
  if (has_upper) ocim_flux[ocean_idx, mass_upper := mass_upper + added_mass_upper]
  if (has_cons)  ocim_flux[ocean_idx, mass_cons  := mass_cons  + added_mass_cons]

  # Remove land-pool cells from flux table
  ocim_flux <- ocim_flux[!land_pool_mask]

  # Recompute FRACAREA-normalized flux for ocean cells that received mass
  ocim_flux[, F_gC_m2_yr := mass_gC_yr / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  if (has_lower) ocim_flux[, F_lower := mass_lower / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  if (has_upper) ocim_flux[, F_upper := mass_upper / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  if (has_cons)  ocim_flux[, F_cons  := mass_cons  / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]

  cat(sprintf("  Redistribution complete (6-stage cascade).\n"))
  cat(sprintf("  Stage 1 (ideal):                %d\n", n_land_pool - n_relax_ocean_component - n_relax_shelf - n_expand_radius - n_connectivity_fallback - n_absolute_fallback))
  cat(sprintf("  Stage 2 (relax shelf comp):      %d  [n_relax_ocean_component]\n", n_relax_ocean_component))
  cat(sprintf("  Stage 3 (coastal wet, <= R_KM):  %d  [n_relax_shelf]\n", n_relax_shelf))
  cat(sprintf("  Stage 4 (any ocean, <= 2*R_KM):  %d  [n_expand_radius]\n", n_expand_radius))
  cat(sprintf("  Stage 5 (connectivity fallback): %d  [n_connectivity_fallback]\n", n_connectivity_fallback))
  cat(sprintf("  Stage 6 (absolute fallback):     %d  [n_absolute_fallback]\n", n_absolute_fallback))
  cat(sprintf("  n_cross_basin:                   %d\n", n_cross_basin))
  cross_scomp_pct <- if (total_input_gC > 0) 100 * cross_scomp_mass_gC / total_input_gC else 0
  cat(sprintf("  n_cross_shelf_component:         %d  (%.2f%% of flux)\n", n_cross_shelf_component, cross_scomp_pct))
  cat(sprintf("  (backward compat) n_relax_component: %d, n_fallback_nearest: %d\n",
              n_relax_component, n_fallback_nearest))
  if (length(redist_distances) > 0) {
    mean_redist_km <- mean(redist_distances)
    cat(sprintf("  Mean redistribution distance: %.1f km\n", mean_redist_km))
  } else {
    mean_redist_km <- 0
  }

} else {
  cat("  No land-pool cells — skipping redistribution.\n")
  # Still need FRACAREA normalization for uncertainty bounds
  if (has_lower) ocim_flux[, F_lower := mass_lower / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  if (has_upper) ocim_flux[, F_upper := mass_upper / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  if (has_cons)  ocim_flux[, F_cons  := mass_cons  / (ocim_area_m2 * pmax(oceanfrac, OCEANFRAC_FLOOR))]
  mean_redist_km <- 0
}

cat(sprintf("  Active OCIM cells after redistribution: %d\n", nrow(ocim_flux)))

# Fail-fast: if weight underflow guard failed for any reason, NaN mass would propagate
# through mass-fix and corrupt J_umol. Catch here before any further computation.
n_nan_mass <- sum(is.nan(ocim_flux$mass_gC_yr))
if (n_nan_mass > 0) {
  stop(sprintf(
    "FATAL: NaN in mass_gC_yr for %d cells after redistribution. Weight normalization produced NaN.",
    n_nan_mass))
}

# ── M4: final vertical profile per recipient column ─────────────────────────
# Non-donor columns: their profile stays in place. Donor columns: each layer goes
# to the recipients with the SAME weights w as the 2D mass (moved mass keeps its depth).
# A zero-mass donor is skipped by M3 (next): its rows disappear here as in ocim_flux.
transfers <- rbindlist(trace_list)
profile_own  <- profile[!land_ij, on = .(i_ocim, j_ocim)]
profile_mv   <- if (nrow(transfers) > 0) {
  x <- merge(profile, transfers, by = c("i_ocim", "j_ocim"), allow.cartesian = TRUE)
  for (cl in agg_cols) set(x, j = cl, value = x[[cl]] * x$w)
  x[, `:=`(i_ocim = rec_i, j_ocim = rec_j)][, c("i_ocim", "j_ocim", "k_src", agg_cols), with = FALSE]
} else profile[0]
profile <- rbindlist(list(profile_own, profile_mv), use.names = TRUE)[
  , lapply(.SD, sum), by = .(i_ocim, j_ocim, k_src), .SDcols = agg_cols]
rm(profile_own, profile_mv)

# LOCK8 (blocking, M4): the profile sums, column by column, to the 2D mass of ocim_flux.
chk <- merge(profile[, .(m3 = sum(mass_gC_yr)), by = .(i_ocim, j_ocim)],
             ocim_flux[, .(i_ocim, j_ocim, m2 = mass_gC_yr)], by = c("i_ocim", "j_ocim"), all = TRUE)
lock8_err <- max(abs(chk$m3 - chk$m2) / pmax(abs(chk$m2), 1), na.rm = FALSE)
cat(sprintf("  LOCK8 vertical profile = 2D mass: max gap %.2e (threshold 1e-9) %s\n",
            lock8_err, if (isTRUE(lock8_err <= 1e-9)) "PASS" else "FAIL"))
if (!isTRUE(lock8_err <= 1e-9)) stop("LOCK8 FAIL: the vertical profile does not sum to the column mass")
rm(chk)

# ── Mass-fix (M2): correct double-precision residual ──────────────────────────

cat("\n--- Mass-fix (M2) ---\n")
sum_after <- sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
residual_gC <- total_input_gC - sum_after
# LOCK1 applies to this residual, BEFORE the fix. Measured after the fix, which brings the
# total back to the source mass by construction, it could only catch a zero or non-finite
# total, not a mass loss of the remapping (e.g. a skipped land cell).
lock1_err <- abs(residual_gC) / total_input_gC
cat(sprintf("  Residual before fix: %.4e gC/yr (%.2e relative)\n", residual_gC, lock1_err))

if (sum_after > 0) {
  scale_factor <- total_input_gC / sum_after
  ocim_flux[, mass_gC_yr := mass_gC_yr * scale_factor]
  ocim_flux[, F_gC_m2_yr := F_gC_m2_yr * scale_factor]
  if (has_lower) { ocim_flux[, mass_lower := mass_lower * scale_factor]; ocim_flux[, F_lower := F_lower * scale_factor] }
  if (has_upper) { ocim_flux[, mass_upper := mass_upper * scale_factor]; ocim_flux[, F_upper := F_upper * scale_factor] }
  if (has_cons)  { ocim_flux[, mass_cons  := mass_cons  * scale_factor]; ocim_flux[, F_cons  := F_cons  * scale_factor] }
  for (cl in agg_cols) set(profile, j = cl, value = profile[[cl]] * scale_factor)  # M4: same factor
}

sum_after_fix <- sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
cat(sprintf("  After fix: relative error = %.2e\n", abs(1 - sum_after_fix / total_input_gC)))

# ── Anti-artifact locks (Step G) ──────────────────────────────────────────────

cat("\n--- Anti-artifact locks ---\n")

# Lock 1: Mass conservation of the remapping, measured on the residual BEFORE the mass-fix (see above)
cat(sprintf("  LOCK1 mass conservation (before fix): %.2e (threshold: %.2e) %s\n",
            lock1_err, LOCK1_MASS_TOL, if (isTRUE(lock1_err <= LOCK1_MASS_TOL)) "PASS" else "FAIL"))
if (!isTRUE(lock1_err <= LOCK1_MASS_TOL)) stop(sprintf("LOCK1 FAIL: mass conservation error before fix %.2e > %.2e", lock1_err, LOCK1_MASS_TOL))

# Lock 2: No cross-basin transfer
cat(sprintf("  LOCK2 cross-basin transfers: %d %s\n",
            n_cross_basin, if (n_cross_basin == 0) "PASS" else "FAIL"))
if (n_cross_basin > 0) stop(sprintf("LOCK2 FAIL: %d cross-basin transfers detected", n_cross_basin))

# Lock 2b: Cross-shelf-component — DIAGNOSTIC (catches isthmus crossings missed by LOCK2)
cross_scomp_pct <- if (total_input_gC > 0) 100 * cross_scomp_mass_gC / total_input_gC else 0
cat(sprintf("  LOCK2b cross-shelf-component: %d events, %.2f%% of flux %s\n",
            n_cross_shelf_component, cross_scomp_pct,
            if (cross_scomp_pct < 1) "OK" else "WARN"))

# Lock 3: Shelf fraction — DIAGNOSTIC
ocim_flux[, on_shelf := as.logical(shelf.mask2d[cbind(i_ocim, j_ocim)])]
shelf_frac <- sum(ocim_flux$mass_gC_yr[ocim_flux$on_shelf], na.rm = TRUE) / sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
cat(sprintf("  LOCK3 shelf fraction: %.1f%% (ref: %.0f%%) %s\n",
            100 * shelf_frac, 100 * LOCK3_SHELF_FLOOR,
            if (shelf_frac >= LOCK3_SHELF_FLOOR) "OK" else "WARN"))

# Lock 4: Hotspot cap — DIAGNOSTIC
max_share <- max(ocim_flux$mass_gC_yr, na.rm = TRUE) / sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
cat(sprintf("  LOCK4 max cell share: %.2f%% (ref: %.0f%%) %s\n",
            100 * max_share, 100 * LOCK4_HOTSPOT_CAP,
            if (max_share <= LOCK4_HOTSPOT_CAP) "OK" else "WARN"))

# Lock 5: Mean redistribution distance — DIAGNOSTIC
# mean_redist_km can be NaN when all redistributed cells have identical coordinates
lock5_label <- if (is.nan(mean_redist_km)) "NaN" else sprintf("%.0f", mean_redist_km)
lock5_status <- if (is.nan(mean_redist_km)) "WARN(NaN)" else if (mean_redist_km <= LOCK5_MEAN_DIST_CAP) "OK" else "WARN"
cat(sprintf("  LOCK5 mean redist distance: %s km (ref: %d km) %s\n",
            lock5_label, LOCK5_MEAN_DIST_CAP, lock5_status))

# Lock 6: Ocean-fraction density cap
positive_flux <- ocim_flux$F_gC_m2_yr[is.finite(ocim_flux$F_gC_m2_yr) & ocim_flux$F_gC_m2_yr > 0]
if (length(positive_flux) > 1) {
  max_density_ratio <- max(positive_flux) / median(positive_flux)
} else {
  max_density_ratio <- 1
}
lock6_status <- if (is.nan(max_density_ratio) || is.na(max_density_ratio)) "WARN(NaN)" else if (max_density_ratio <= LOCK6_DENSITY_RATIO) "PASS" else "FAIL"
cat(sprintf("  LOCK6 max/median density ratio: %s (threshold: %.0f) %s\n",
            if (is.finite(max_density_ratio)) sprintf("%.0f", max_density_ratio) else "NaN",
            LOCK6_DENSITY_RATIO, lock6_status))
if (isTRUE(max_density_ratio > LOCK6_DENSITY_RATIO)) stop(sprintf("LOCK6 FAIL: max/median density ratio %.0f > %.0f", max_density_ratio, LOCK6_DENSITY_RATIO))

# Lock 7: No absolute fallback (every land cell must find a connected ocean target)
cat(sprintf("  LOCK7 absolute fallback: %d %s\n",
            n_absolute_fallback, if (n_absolute_fallback == 0) "PASS" else "FAIL"))
if (n_absolute_fallback > 0) stop(sprintf("LOCK7 FAIL: %d land cells had no ocean target in their ocean component", n_absolute_fallback))

# LOCK1/2/6/7 are blocking (stop() above). LOCK2b/3/4/5 are DIAGNOSTIC and must
# be reported as such, not absorbed into a "PASSED".
lock_warnings <- c(
  if (isTRUE(cross_scomp_pct >= 1))
    sprintf("LOCK2b (%.2f%% of the flux, ref 1%%)", cross_scomp_pct),
  if (isTRUE(shelf_frac < LOCK3_SHELF_FLOOR))
    sprintf("LOCK3 (%.1f%% on the shelf, ref %.0f%%)", 100 * shelf_frac, 100 * LOCK3_SHELF_FLOOR),
  if (isTRUE(max_share > LOCK4_HOTSPOT_CAP))
    sprintf("LOCK4 (%.2f%% in a single cell, ref %.0f%%)", 100 * max_share, 100 * LOCK4_HOTSPOT_CAP),
  if (lock5_status != "OK")
    sprintf("LOCK5 (%s km, ref %d km)", lock5_label, LOCK5_MEAN_DIST_CAP),
  if (lock6_status != "PASS")
    sprintf("LOCK6 (%s)", lock6_status)
)
if (length(lock_warnings) == 0) {
  cat("  Locks: blocking ones PASS, no diagnostic warning.\n")
} else {
  cat(sprintf("  Locks: blocking PASS; %d diagnostic WARNING(S) -- %s\n",
              length(lock_warnings), paste(lock_warnings, collapse = " ; ")))
}

# ── Section 4: Inject at the seabed depth (M4) ───────────────────────────────

cat(sprintf("\n--- Section 4: Inject (mode %s) ---\n", INJECTION_MODE))

# For each (i,j,k): concentration tendency in cell k of column (i,j)
# J_umol = mass_gC_yr * 1e6 / (M_C * RHO * DZT_inj * ocim_area_m2)  [umol kg-1 yr-1]
# NOT from F_gC_m2_yr (which divides by oceanfrac — OCIM applies J to full cell volume)
# k_inj = min(k_src, kbot): the layer of the true seabed, capped at the bottom of the recipient column.
# bottom mode: k_src = nk, so k_inj = kbot (former behaviour).

profile[, k_bottom_col := kbot[cbind(i_ocim, j_ocim)]]
# Cells with kbot = 0 stop the run: mass lost silently would bias the balance.
if (any(is.na(profile$k_bottom_col) | profile$k_bottom_col == 0))
  stop("Section 4: ", sum(is.na(profile$k_bottom_col) | profile$k_bottom_col == 0), " rows in a column without water")
profile[, k_inj := pmin(k_src, k_bottom_col)]
capped_pct <- 100 * sum(profile$mass_gC_yr[profile$k_src > profile$k_bottom_col]) / sum(profile$mass_gC_yr)

inj <- profile[, lapply(.SD, sum), by = .(i_ocim, j_ocim, k_inj), .SDcols = agg_cols]
inj[, ocim_area_m2 := DXT3d[cbind(i_ocim, j_ocim, 1L)] * DYT3d[cbind(i_ocim, j_ocim, 1L)]]  # as in ocim_flux
inj[, DZT_inj := DZT3d[cbind(i_ocim, j_ocim, k_inj)]]
inj[, z_inj_base := layer_base[k_inj]]            # depth of the bottom of the injection layer [m]
inj[, J_umol := mass_gC_yr * 1e6 / (M_C * RHO * DZT_inj * ocim_area_m2)]
if (has_lower) inj[, J_lower := mass_lower * 1e6 / (M_C * RHO * DZT_inj * ocim_area_m2)]
if (has_upper) inj[, J_upper := mass_upper * 1e6 / (M_C * RHO * DZT_inj * ocim_area_m2)]
if (has_cons)  inj[, J_cons  := mass_cons  * 1e6 / (M_C * RHO * DZT_inj * ocim_area_m2)]

wq <- function(x, w, q) { o <- order(x); x[o][which(cumsum(w[o]) >= q * sum(w))[1]] }
cat(sprintf("  Injection cells: %d, in %d columns\n", nrow(inj), uniqueN(inj[, .(i_ocim, j_ocim)])))
cat(sprintf("  Depth of the bottom of the injection layer (mass-weighted): q10 %.0f m, median %.0f m, q90 %.0f m\n",
            wq(inj$z_inj_base, inj$mass_gC_yr, 0.1), wq(inj$z_inj_base, inj$mass_gC_yr, 0.5),
            wq(inj$z_inj_base, inj$mass_gC_yr, 0.9)))
cat(sprintf("  Thickness of the injection layer: min %.1f m, median %.1f m, max %.1f m\n",
            min(inj$DZT_inj), median(inj$DZT_inj), max(inj$DZT_inj)))
cat(sprintf("  Mass capped by the column bottom (k_src > kbot): %.2f%%\n", capped_pct))

# Compute linear index into M3d for each (i,j,k_inj)
# MATLAB uses column-major: idx = i + (j-1)*ni + (k-1)*ni*nj
inj[, linear_idx := i_ocim + (j_ocim - 1L) * ni + (k_inj - 1L) * ni * nj]

# Map linear index to iocn position (vector index in Jdredge)
iocn_map <- data.table(linear_idx = iocn, vec_pos = seq_along(iocn))
inj <- merge(inj, iocn_map, by = "linear_idx", all.x = TRUE)

# Cells outside iocn stop the run (mass lost).
unmatched <- sum(is.na(inj$vec_pos))
if (unmatched > 0) stop(sprintf("Section 4: %d injection cells outside iocn", unmatched))
if (anyDuplicated(inj$vec_pos)) stop("Section 4: duplicated injection cell")

# ── Section 5: Build Jdredge vector ──────────────────────────────────────────

cat("\n--- Section 5: Build Jdredge vector ---\n")

Jdredge <- rep(0, m)
Jdredge[inj$vec_pos] <- inj$J_umol

Jdredge_lower <- rep(0, m)
Jdredge_upper <- rep(0, m)
Jdredge_conservative <- rep(0, m)

if (has_lower) Jdredge_lower[inj$vec_pos] <- inj$J_lower
if (has_upper) Jdredge_upper[inj$vec_pos] <- inj$J_upper
if (has_cons)  Jdredge_conservative[inj$vec_pos] <- inj$J_cons

n_nonfinite <- sum(!is.finite(Jdredge))
if (n_nonfinite > 0) stop(sprintf(
  "FATAL: %d non-finite values in Jdredge. Patches did not resolve the root cause — aborting export.",
  n_nonfinite))
n_active <- sum(Jdredge > 0)
cat(sprintf("  Jdredge length: %d (m)\n", m))
cat(sprintf("  Active cells: %d\n", n_active))
cat(sprintf("  Jdredge > 0 summary:\n"))
print(summary(Jdredge[Jdredge > 0]))

# ── Section 5b: Amplified scenarios for sensitivity analysis ─────────────────

cat("\n--- Section 5b: Amplified scenarios (10x, 100x) ---\n")

Jdredge_10x  <- Jdredge * 10
Jdredge_100x <- Jdredge * 100

Jdredge_lower_10x  <- Jdredge_lower * 10
Jdredge_lower_100x <- Jdredge_lower * 100
Jdredge_upper_10x  <- Jdredge_upper * 10
Jdredge_upper_100x <- Jdredge_upper * 100
Jdredge_conservative_10x  <- Jdredge_conservative * 10
Jdredge_conservative_100x <- Jdredge_conservative * 100

cat(sprintf("  Jdredge_10x  max: %.4e umol/kg/yr\n", max(Jdredge_10x)))
cat(sprintf("  Jdredge_100x max: %.4e umol/kg/yr\n", max(Jdredge_100x)))

# ── Section 6: QA/QC — Split metrics + enhanced diagnostics ──────────────────

cat("\n--- Section 6: QA/QC ---\n")

# 1. Mass-space check (LOCK1 already validates this)
total_mass_after_redist <- sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
conservation_err_mass <- abs(1 - total_mass_after_redist / total_input_gC)

# 2. Round-trip check: J_umol * VOL * RHO * M_C / 1e6 -> gC/yr
VOL_active  <- DXT3d[cbind(inj$i_ocim, inj$j_ocim, inj$k_inj)] *
               DYT3d[cbind(inj$i_ocim, inj$j_ocim, inj$k_inj)] *
               inj$DZT_inj

total_output_gC <- sum(inj$J_umol * VOL_active * RHO * M_C / 1e6, na.rm = TRUE)
conservation_err_roundtrip <- abs(1 - total_output_gC / total_input_gC)

conservation_error_pct <- conservation_err_roundtrip * 100
conservation_error_mass_pct <- conservation_err_mass * 100

cat(sprintf("  Input  total: %.6e gC/yr\n", total_input_gC))
cat(sprintf("  Output total: %.6e gC/yr\n", total_output_gC))
cat(sprintf("  Conservation (mass-space):  %.2e\n", conservation_err_mass))
cat(sprintf("  Conservation (round-trip):  %.2e\n", conservation_err_roundtrip))

if (conservation_err_roundtrip > 0.01) {
  cat("  WARNING: Round-trip conservation error exceeds 1%!\n")
} else {
  cat("  OK: Conservation within tolerance.\n")
}

# --- Split metrics ---

# 1. land_pool_pct: fraction of input mass that was on OCIM land
land_pool_pct <- 100 * flux_eff_land / total_input_gC
cat(sprintf("\n  land_pool_pct: %.2f%%\n", land_pool_pct))

# 2. offshore_displacement_pct: fraction of redistributed mass displaced >50 km offshore
# For each redistributed cell, check distance shift vs seed baseline
offshore_displacement_pct <- 0
if (n_land_pool > 0 && length(redist_distances) > 0) {
  # Approximate: use mean redistribution distance as proxy
  # Detailed tracking would require per-transfer baseline, which is complex in R
  # We use a simpler approach: check recipient cells' dist_to_coast vs overall pattern
  ocim_flux[, dist_coast := dist.to.coast2d[cbind(i_ocim, j_ocim)]]

  # Recipients that are far from coast (>50km further than median shelf dist)
  shelf_median_dist <- median(dist.to.coast2d[shelf.mask2d > 0], na.rm = TRUE)
  offshore_mask <- ocim_flux$dist_coast > (shelf_median_dist + OFFSHORE_DELTA_KM)
  offshore_flux <- sum(ocim_flux$mass_gC_yr[offshore_mask], na.rm = TRUE)
  offshore_displacement_pct <- 100 * offshore_flux / total_input_gC
  cat(sprintf("  offshore_displacement_pct: %.2f%%\n", offshore_displacement_pct))
} else {
  cat("  offshore_displacement_pct: 0.00% (no redistribution)\n")
}

# 3. shelf_flux_pct
shelf_flux_pct <- 100 * shelf_frac
cat(sprintf("  shelf_flux_pct: %.2f%%\n", shelf_flux_pct))

# 4. max_cell_share_pct
max_cell_share_pct <- 100 * max_share
cat(sprintf("  max_cell_share_pct: %.2f%%\n", max_cell_share_pct))

# 5. top5_flux_share (top-1% collapses to 1 cell on a 78-cell grid, identical to max_cell_share)
n_top5 <- min(5L, nrow(ocim_flux))
top5_mass <- sum(sort(ocim_flux$mass_gC_yr, decreasing = TRUE)[1:n_top5])
top5_flux_share <- 100 * top5_mass / total_input_gC
cat(sprintf("  top5_flux_share: %.2f%%\n", top5_flux_share))

# 6. median_dist_coast_km & p90_dist_coast_km
# Flux-weighted absolute coastal distance of final recipient cells (NOT origin-to-recipient delta)
ocim_flux[, dist_coast := dist.to.coast2d[cbind(i_ocim, j_ocim)]]
if (nrow(ocim_flux) > 0) {
  ord <- order(ocim_flux$dist_coast)
  cum_flux <- cumsum(ocim_flux$mass_gC_yr[ord]) / sum(ocim_flux$mass_gC_yr, na.rm = TRUE)
  median_dist_coast_km <- ocim_flux$dist_coast[ord[which.min(abs(cum_flux - 0.5))]]
  p90_dist_coast_km    <- ocim_flux$dist_coast[ord[which.min(abs(cum_flux - 0.9))]]
} else {
  median_dist_coast_km <- 0
  p90_dist_coast_km    <- 0
}
cat(sprintf("  median_dist_coast_km: %.1f\n", median_dist_coast_km))
cat(sprintf("  p90_dist_coast_km: %.1f\n", p90_dist_coast_km))

# Spot check: top-10 injection cells with enhanced columns
top10 <- inj[order(-J_umol)][1:min(10, nrow(inj))]
top10[, lon_center := lon2d[cbind(i_ocim, j_ocim)]]
top10[, lat_center := lat2d[cbind(i_ocim, j_ocim)]]
top10[, of := oceanfrac2d[cbind(i_ocim, j_ocim)]]
top10[, d2c := dist.to.coast2d[cbind(i_ocim, j_ocim)]]
top10[, bottom := depth2d[cbind(i_ocim, j_ocim)]]
cat("\n  Top-10 Jdredge cells:\n")
cat(sprintf("  %4s %8s %8s %12s %4s %8s %8s %8s %8s %8s\n", "rank", "lon", "lat", "J_umol", "k",
            "z_base", "DZT_inj", "bottom_col", "oceanfr", "d2coast"))
for (r in seq_len(nrow(top10))) {
  cat(sprintf("  %4d %8.2f %8.2f %12.4e %4d %8.0f %8.1f %8.0f %8.3f %8.1f\n",
              r, top10$lon_center[r], top10$lat_center[r], top10$J_umol[r], top10$k_inj[r],
              top10$z_inj_base[r], top10$DZT_inj[r], top10$bottom[r], top10$of[r], top10$d2c[r]))
}

# Mass distribution by depth of the bottom of the injection layer
cat("\n  Mass by injection depth (layer bottom):\n")
depth_summary <- inj[, .(n_cells = .N, mass_share = sum(mass_gC_yr) / sum(inj$mass_gC_yr)),
                     by = .(z_bin = cut(z_inj_base, breaks = c(0, 50, 200, 500, 2000, Inf)))][order(z_bin)]
print(depth_summary)

# ── Section 7: Export .mat ────────────────────────────────────────────────────

cat("\n--- Section 7: Export .mat ---\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
mat_filename <- sprintf("jdredge_ocim2_48l_%s.mat", timestamp)
mat_path <- file.path(OUT_DIR, mat_filename)

R.matlab::writeMat(
  mat_path,
  # Primary dredging forcing
  Jdredge              = matrix(Jdredge, ncol = 1),
  Jdredge.lower        = matrix(Jdredge_lower, ncol = 1),
  Jdredge.upper        = matrix(Jdredge_upper, ncol = 1),
  Jdredge.conservative = matrix(Jdredge_conservative, ncol = 1),
  # Amplified scenarios
  Jdredge.10x          = matrix(Jdredge_10x, ncol = 1),
  Jdredge.100x         = matrix(Jdredge_100x, ncol = 1),
  Jdredge.lower.10x    = matrix(Jdredge_lower_10x, ncol = 1),
  Jdredge.lower.100x   = matrix(Jdredge_lower_100x, ncol = 1),
  Jdredge.upper.10x    = matrix(Jdredge_upper_10x, ncol = 1),
  Jdredge.upper.100x   = matrix(Jdredge_upper_100x, ncol = 1),
  Jdredge.conservative.10x  = matrix(Jdredge_conservative_10x, ncol = 1),
  Jdredge.conservative.100x = matrix(Jdredge_conservative_100x, ncol = 1),
  # Backward-compatible aliases for co2model.m
  Jtrawl               = matrix(Jdredge, ncol = 1),
  Jtrawl.lower         = matrix(Jdredge_lower, ncol = 1),
  Jtrawl.upper         = matrix(Jdredge_upper, ncol = 1),
  Jtrawl.conservative  = matrix(Jdredge_conservative, ncol = 1),
  Jtrawl.10x           = matrix(Jdredge_10x, ncol = 1),
  Jtrawl.100x          = matrix(Jdredge_100x, ncol = 1),
  # Metadata — core
  n.cells.active      = n_active,
  total.flux.gC.yr    = total_input_gC,
  conservation.error.pct = conservation_error_pct,
  conservation.error.mass.pct = conservation_error_mass_pct,
  timestamp.str       = timestamp,
  m                   = m,
  # Metadata — split metrics
  land.pool.pct                = land_pool_pct,
  offshore.displacement.pct    = offshore_displacement_pct,
  shelf.flux.pct               = shelf_flux_pct,
  max.cell.share.pct           = max_cell_share_pct,
  top5.flux.share              = top5_flux_share,
  median.dist.coast.km         = median_dist_coast_km,
  p90.dist.coast.km            = p90_dist_coast_km,
  # Metadata — redistribution parameters
  redistribution.R.km          = R_KM,
  redistribution.lambda.km     = LAMBDA_KM,
  redistribution.k.max         = K_MAX,
  redistribution.alpha         = ALPHA,
  redistribution.oceanfrac.floor = OCEANFRAC_FLOOR,
  # Metadata — redistribution quality counters (6-stage cascade)
  n.relax.component            = n_relax_component,
  n.relax.ocean.component      = n_relax_ocean_component,
  n.relax.shelf                = n_relax_shelf,
  n.expand.radius              = n_expand_radius,
  n.connectivity.fallback      = n_connectivity_fallback,
  n.absolute.fallback          = n_absolute_fallback,
  n.fallback.nearest           = n_fallback_nearest,
  n.cross.basin                = n_cross_basin,
  n.cross.shelf.component      = n_cross_shelf_component,
  cross.scomp.flux.pct         = cross_scomp_pct,
  mean.redist.km               = mean_redist_km,
  # Metadata — M4 (injection depth)
  injection.mode               = INJECTION_MODE,
  injection.capped.pct         = capped_pct
)

cat(sprintf("  Written: %s (%.1f MB)\n", mat_path, file.size(mat_path) / 1e6))

# M4: the injection table (column, layer, mass), to reread or redo the diagnostics without a rerun
inj_path <- sub("\\.mat$", "_injection.csv", mat_path)
fwrite(inj[order(i_ocim, j_ocim, k_inj), c("i_ocim", "j_ocim", "k_inj", "z_inj_base", "DZT_inj", agg_cols, "J_umol"), with = FALSE],
       inj_path)
cat(sprintf("  Written: %s\n", inj_path))

cat(sprintf("\nStep 7 complete: %s\n", format(Sys.time())))
cat("=== DONE ===\n")
