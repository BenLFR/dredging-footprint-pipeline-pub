#!/usr/bin/env Rscript
# =============================================================================
# Step 5: merge of the SAR tiles -> final f_i grid
# =============================================================================
# Features:
#  A: Dedup keeps ALL cells (prefer owner tile, fallback to max sum_dw)
#  B: k_fast mapping works with basin-keyed OR province-keyed YAML
#  C: p_l = pl_sum / n_with_pl (worker-corrected schema)
#  D: Centroid coordinates use constants.R (bottom-up)
#  E: GeoTIFF row inversion (terra top-down)
#  F: Minimal but actionable diagnostics
#  G: depth_mode "mechanistic": alpha and w per pool, derived from
#      mixing (D_b) and burial (omega).
#      depth_mode "fixed" (or absent from the YAML) = the earlier calculation, unchanged.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(rlang)
  library(tidyselect)
  library(lubridate)
  library(yaml)
  library(arrow)
  library(terra)
})

# --- Source constants.R -------------------------------------------------------
this_file <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  if (length(f)) return(normalizePath(f))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile))
  stop("Cannot locate running script")
}
script_dir <- dirname(this_file())
# constants.R is looked for next to this script, then in the parent directory --
# never a local copy: a duplicated constants.R would diverge.
.cst <- file.path(script_dir, "constants.R")
if (!file.exists(.cst)) .cst <- file.path(dirname(script_dir), "constants.R")
if (!file.exists(.cst)) stop("constants.R not found from ", script_dir)
source(.cst)

sf_use_s2(FALSE)

cat("Grid:", GRID_COLS, "x", GRID_ROWS, " cells of", CELL_SIZE_M, "m\n")
cat("Extent X:", WORLD_XMIN, "->", WORLD_XMAX, "\n")
cat("Extent Y:", WORLD_YMIN, "->", WORLD_YMAX, "\n\n")

# --- Command-line parameters --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
scenario <- if(length(args) > 0) args[1] else "default"
cat("Final merge - Scenario:", scenario, "\n")

# --- Loading the YAML parameters ----------------------------------------------
# Path can be overridden by FI_PARAM_YAML (tests).
param_yaml <- Sys.getenv("FI_PARAM_YAML", file.path(Sys.getenv("DFP_REPO"), "config/fi_parameters_with_freshness.yaml"))
cat("YAML:", param_yaml, "\n")
params_raw <- yaml::read_yaml(param_yaml)

get_scenario <- function(name) {
  s <- params_raw$scenarios[[name]]
  if (!is.null(s$inherit)) {
    parent <- get_scenario(s$inherit)
    s$inherit <- NULL
    modifyList(parent, s)
  } else s
}
par <- get_scenario(scenario)

# --- Check of the required parameters -----------------------------------------
required <- c("alpha_dep","fast_fraction","slow_k","preservation_factor","k_fast")
miss <- setdiff(required, names(par))
if(length(miss)) stop("Missing YAML parameters: ", paste(miss, collapse=", "))
if (length(par$k_fast) == 0) stop("k_fast empty in the YAML")

# --- Table of regional k (raw from YAML) --------------------------------------
k_yaml_names <- names(par$k_fast)
k_yaml_vals  <- unlist(par$k_fast) *
                ifelse(is.null(par$k_fast_multiplier), 1, par$k_fast_multiplier)

alpha_dep     <- par$alpha_dep
fast_frac     <- par$fast_fraction
slow_k        <- par$slow_k
preserv_fact  <- ifelse(is.null(par$preservation_factor), 0.87,  par$preservation_factor)  # p_r
depletion_di  <- ifelse(is.null(par$depletion_factor),    0.272, par$depletion_factor)       # d_i

# --- G: depth weights --------------------------------------------------------
# "fixed"       : earlier calculation -- 0-5 cm weight 1, 5-10 cm weight alpha_dep, nothing
#                 below 10 cm. Kept to reproduce earlier runs and for
#                 the upper_bound scenario (slow_k = 0 leaves "mechanistic" undefined).
# "mechanistic": for each pool of rate constant k,
#                 0-10 cm  : profile exp(-z/l), 1/l = (sqrt(omega^2 + 4 D_b k) - omega) / (2 D_b)
#                            -- faunal mixing carries fresh material down,
#                            degradation consumes it on the way. This is the steady-state
#                            solution of eq. 3 of Arndt et al. (2013, ESR 123:53, p. 66)
#                            with constant coefficients, first-order kinetics (their eq. 5).
#                 below 10 cm: no more mixing, material only moves down by
#                            burial at speed omega; decay
#                            exp(-(z - 10) k / omega).
#                 junction at 10 cm: concentration AND flux omega C - D_b C' continuous,
#                            so C'(10-) = 0 and the 0-10 cm profile has two exponentials,
#                            exp(-a z) + (a/b) exp(-a h) exp(b (z - h)).
#                 D_b and omega: Middelburg et al. (1996, GBC 10:661, table 1),
#                 read at the source; values at z = 50 m in the YAML.
depth_mode <- if (is.null(par$depth_mode)) "fixed" else par$depth_mode
if (!depth_mode %in% c("fixed", "mechanistic")) stop("depth_mode unknown: ", depth_mode)
if (depth_mode == "mechanistic") {
  Db_cm2   <- par$Db_cm2_yr
  omega_cm <- par$omega_cm_yr
  if (is.null(Db_cm2) || is.null(omega_cm) || !(Db_cm2 > 0) || !(omega_cm > 0))
    stop("depth_mode mechanistic: Db_cm2_yr and omega_cm_yr required and > 0")
  if (!(slow_k > 0) || any(!(k_yaml_vals > 0)))
    stop("depth_mode mechanistic: slow_k and every k_fast must be > 0")
}

# Weight of a pool of rate constant k (a-1), averaged over the excavated thickness,
# relative to the mean of the 0-5 cm layer (which carries p_l).
# Profile over 0-10 cm (h = 10): g(z) = exp(-a z) + (a/b) exp(-a h) exp(b (z - h)), with
#   a = (sqrt(omega^2 + 4 D_b k) - omega) / (2 D_b),  b = (sqrt(omega^2 + 4 D_b k) + omega) / (2 D_b),
# roots of D_b r^2 -/+ omega r - k = 0; the second term cancels g'(h-), which makes the flux
# omega g - D_b g' continuous at h. Below h: g(h) exp(-(z - h)/L), L = omega/k.
# Integral I(x) = int_0^x g = (1 - exp(-a x))/a + (a/b^2) exp(-a h) (exp(b (x - h)) - exp(-b h)).
#   alpha = mean(5-10) / mean(0-5)        = (I(10) - I(5)) / I(5)
#   w     = mean(10 - 10+H) / mean(0-5)   = g(h) L (1 - exp(-H/L)) / H / (I(5)/5)
# with g(h) = exp(-a h) (1 + a/b). Written so as never to form exp(b z) alone.
# w1, w2, w3: shares of the excavated thickness in 0-5, 5-10 and below 10 cm; H3 in cm.
weight_depth <- function(k, w1, w2, w3, H3) {
  h  <- 10
  rq <- sqrt(omega_cm^2 + 4 * Db_cm2 * k)
  a  <- 2 * k / (rq + omega_cm)                     # = (rq - omega)/(2 D_b), without subtraction
  b  <- (rq + omega_cm) / (2 * Db_cm2)
  I  <- function(x) -expm1(-a * x) / a + (a / b^2) * exp(-a * h) * (exp(b * (x - h)) - exp(-b * h))
  m0 <- I(5) / 5
  al <- (I(10) - I(5)) / I(5)
  gh <- exp(-a * h) * (1 + a / b)
  L  <- omega_cm / k
  wd <- ifelse(H3 > 0, gh * L * (-expm1(-H3 / L)) / (pmax(H3, 1e-12) * m0), 0)
  w1 + w2 * al + w3 * wd
}
# Effective weight on f_i: mean of the weights of the two pools, weighted by their
# contribution to one-year remineralisation. f_i stays LINEAR in p_l_eff.
weight_eff <- function(kf, ks, w1, w2, w3, H3, t = 1) {
  gf <- fast_frac * (1 - exp(-kf * t))
  gs <- (1 - fast_frac) * (1 - exp(-ks * t))
  (gf * weight_depth(kf, w1, w2, w3, H3) + gs * weight_depth(ks, w1, w2, w3, H3)) / (gf + gs)
}

cat("f_i parameters:\n")
cat("   - depth_mode      :", depth_mode, "\n")
if (depth_mode == "mechanistic") {
  cat("   - D_b             :", Db_cm2, "cm2 a-1\n")
  cat("   - omega           :", omega_cm, "cm a-1\n")
  cat("   - alpha_dep       :", alpha_dep, "(NOT USED in mechanistic mode)\n")
} else {
  cat("   - alpha_dep       :", alpha_dep, "\n")
}
cat("   - fast_fraction   :", fast_frac, "\n")
cat("   - slow_k          :", slow_k, "a-1\n")
cat("   - p_r (preserv)   :", preserv_fact, "\n")
cat("   - d_i (depletion):", depletion_di, "\n")
cat("   - k_fast keys     :", paste(k_yaml_names, collapse=", "), "\n")
stopifnot(!is.na(alpha_dep), !is.na(fast_frac), !is.na(slow_k), !is.na(preserv_fact), !is.na(depletion_di))
cat("All parameters are valid\n\n")

# =============================================================================
# 1) Merge of the SAR tiles
# =============================================================================
cat("== 1) Merge of the SAR tiles ==\n")

# --- Declarable I/O -------------------------------------------------------------
MERGE_IN_DIR  <- path.expand(Sys.getenv("MERGE_IN_DIR",  file.path(Sys.getenv("DFP_DATA"), "output_V6_corr")))
MERGE_OUT_DIR <- path.expand(Sys.getenv("MERGE_OUT_DIR", file.path(Sys.getenv("DFP_DATA"), "output_V6_corr")))
MERGE_TILES   <- path.expand(Sys.getenv("PIPE_TILES",    file.path(Sys.getenv("DFP_DATA"), "output_V6/tiles_1000km.gpkg")))

# guard: never write into the PUBLISHED outputs
if (identical(normalizePath(MERGE_OUT_DIR, mustWork = FALSE),
              normalizePath(path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6")), mustWork = FALSE)) &&
    !nzchar(Sys.getenv("PIPE_ALLOW_PUBLISHED"))) {
  stop("Refusing to write into the PUBLISHED outputs ($DFP_DATA/output_V6). ",
       "Set MERGE_OUT_DIR elsewhere.")
}
dir.create(MERGE_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("MERGE_IN_DIR  =", MERGE_IN_DIR,  "\n")
cat("MERGE_OUT_DIR =", MERGE_OUT_DIR, "\n")

sar_files <- list.files(MERGE_IN_DIR, pattern="^sar_.*\\.(parquet|rds)$", full.names=TRUE)
cat("sar tiles found:", length(sar_files), "\n")
if(length(sar_files) == 0) stop("No sar_* file found in ", MERGE_IN_DIR)

parquet_count <- sum(grepl("\\.parquet$", sar_files))
rds_count     <- sum(grepl("\\.rds$", sar_files))
cat("Files found:", length(sar_files), "(Parquet:", parquet_count, "RDS:", rds_count, ")\n")

# SAR loading with tile_id provenance
tile_from_file <- function(f) as.integer(sub("^sar_", "", sub("\\.(parquet|rds)$", "", basename(f))))

sar_list <- vector("list", length(sar_files))
for(i in seq_along(sar_files)) {
  file <- sar_files[i]
  cat("   Loading:", basename(file), "\n")
  dt <- if(grepl("\\.parquet$", file)) as.data.table(read_parquet(file))
        else as.data.table(readRDS(file))
  dt[, source_tile_id := tile_from_file(file)]
  sar_list[[i]] <- dt
}

cat("Assembling the data...\n")
sar_global <- rbindlist(sar_list, fill=TRUE)
rm(sar_list); gc()
setDT(sar_global)

n_rows_raw   <- nrow(sar_global)
n_gridid_raw <- uniqueN(sar_global$grid_id)
cat("Before deduplication:", n_rows_raw, "rows,", n_gridid_raw, "unique grid_id\n")

# =============================================================================
# 2) Deduplication — A: keep ALL cells, no drops
# =============================================================================
cat("\n== 2) Deduplication owner-tile ==\n")

# Cell centres EPSG:6933
cent <- unique(sar_global[, .(grid_id)])
cent[, `:=`(
  col = (grid_id - 1L) %% GRID_COLS,
  row = (grid_id - 1L) %/% GRID_COLS
)]
cent[, `:=`(
  x = WORLD_XMIN + col * CELL_SIZE_M + CELL_SIZE_M / 2,
  y = WORLD_YMIN + row * CELL_SIZE_M + CELL_SIZE_M / 2
)]
cent_sf <- st_as_sf(cent, coords = c("x", "y"), crs = 6933)

# Load the core tiles (without buffer) — auto-detect the layer with tile_id
tiles_path <- MERGE_TILES
if (!file.exists(path.expand(tiles_path))) {
  stop("Tiles file not found: ", tiles_path)
}
ly <- st_layers(path.expand(tiles_path))
layer_ok <- NULL
for (L in ly$name) {
  s <- try(st_read(path.expand(tiles_path), layer = L, quiet = TRUE), silent = TRUE)
  if (!inherits(s, "try-error") && "tile_id" %in% names(s)) { layer_ok <- L; break }
}
if (is.null(layer_ok)) stop("No layer with a tile_id column in ", tiles_path)
tiles <- st_read(path.expand(tiles_path), layer = layer_ok, quiet = TRUE)
cat("Tiles layer:", layer_ok, "->", nrow(tiles), "tiles\n")

# Spatial join: find the owner tile of each cell
owner <- as.data.table(st_join(cent_sf, tiles[, "tile_id"], left = TRUE))
owner[, geometry := NULL]
setnames(owner, "tile_id", "owner_tile_id")

# Merge owner info into sar_global
sar_global <- merge(sar_global, owner, by = "grid_id", all.x = TRUE)

# Robust dedup: for each grid_id, prefer the owner-tile row, else max sum_dw
# NA-safe: owner_tile_id can be NA for cells outside all core tiles
sar_global[, is_owner := as.integer(!is.na(owner_tile_id) & source_tile_id == owner_tile_id)]
# Temp ordering column: NA sum_dw must not win over real values
sar_global[, sort_dw := fifelse(is.na(sum_dw), -Inf, sum_dw)]
setorder(sar_global, grid_id, -is_owner, -sort_dw)
sar_global[, sort_dw := NULL]
sar_global <- sar_global[, .SD[1], by = grid_id]

stopifnot(nrow(sar_global) == uniqueN(sar_global$grid_id))

# Zero-drop check: must keep ALL unique grid_ids
if (nrow(sar_global) != n_gridid_raw) {
  warning("Dedup changed the number of grid_id: ", n_gridid_raw, " -> ", nrow(sar_global))
}

n_with_owner    <- sum(sar_global$is_owner == 1L, na.rm = TRUE)
n_without_owner <- nrow(sar_global) - n_with_owner
cat("Deduplication:", nrow(sar_global), "unique cells\n")
cat("   with owner tile:", n_with_owner, "\n")
cat("   without owner (fallback max sum_dw):", n_without_owner, "\n")

sar_global[, c("source_tile_id", "owner_tile_id", "is_owner") := NULL]

# =============================================================================
# 3) SAR / SVR metrics — C
# =============================================================================
cat("\n== 3) SAR/SVR metrics ==\n")

sar_global[, SAR := sum_dw / CELL_AREA_M2]
sar_global[, p_d := fifelse(sum_dw > 0, sum_dw_pd / sum_dw, NA_real_)]
sar_global[, SVR := SAR * p_d]

# p_l compatible with the corrected worker
if ("n_with_pl" %in% names(sar_global) && "pl_sum" %in% names(sar_global)) {
  sar_global[, p_l := fifelse(n_with_pl > 0, pl_sum / n_with_pl, NA_real_)]
} else if (all(c("sum_d_pl", "sum_d") %in% names(sar_global))) {
  sar_global[, p_l := fifelse(sum_d > 0, sum_d_pl / sum_d, NA_real_)]
} else {
  sar_global[, p_l := NA_real_]
}
# Option B depth-weighting: surface layer [0-5 cm] + capped deep layer [5-10 cm].
# Sediment below 10 cm contributes 0 labile fraction on a 1-year timescale
# (beyond Holocene bioturbated layer; k < 1e-4 a-1, negligible in 365 days).
# "fixed" mode only. The sentence above had no source; it holds for the fast
# pool and is wrong for the slow pool. The "mechanistic" mode replaces it; see G.
# Guard against p_d = 0 or NA (→ treated as 1 m default, consistent with tile worker).
sar_global[, p_d_safe := fifelse(is.na(p_d) | p_d <= 0, 1, p_d)]
sar_global[, `:=`(
  w1 = pmin(0.05, p_d_safe) / p_d_safe,
  w2 = pmin(0.05, pmax(p_d_safe - 0.05, 0)) / p_d_safe,  # Option B: capped at 5 cm thick
  w3 = pmax(p_d_safe - 0.10, 0) / p_d_safe,              # G: below 10 cm
  H3_cm = pmax(p_d_safe - 0.10, 0) * 100                 # G: thickness below 10 cm
)]
if (depth_mode == "fixed") {
  sar_global[, p_l_eff := w1 * p_l + w2 * p_l * alpha_dep]
  sar_global[, c("w1", "w2", "w3", "H3_cm") := NULL]
}
# in mechanistic mode, p_l_eff is computed in section 5, once k is known
sar_global[, p_d_safe := NULL]

cat("SAR  : min=", min(sar_global$SAR, na.rm=TRUE),
    " max=", max(sar_global$SAR, na.rm=TRUE), "\n")
cat("SVR  : min=", min(sar_global$SVR, na.rm=TRUE),
    " max=", max(sar_global$SVR, na.rm=TRUE), "\n")
cat("p_l  : non-NA =", sum(!is.na(sar_global$p_l)), "/", nrow(sar_global), "\n")

# =============================================================================
# 4) Longhurst provinces + k_fast mapping — B
# =============================================================================
cat("\n== 4) Provinces Longhurst ==\n")

longhurst <- st_read(
  path.expand(file.path(Sys.getenv("DFP_DATA"), "configuration/longhurst_v4_2010/Longhurst_world_v4_2010.shp")),
  quiet = TRUE
)
longhurst <- st_transform(longhurst, 6933)

# Detect the code column
code_col <- grep("code$", names(longhurst), value = TRUE, ignore.case = TRUE)[1]
if(is.na(code_col) || code_col == "") stop("Longhurst code field not found in the shapefile")
longhurst <- dplyr::rename(longhurst, longhurst_pr = !!sym(code_col))

prov_codes <- unique(longhurst$longhurst_pr)
cat("Longhurst provinces loaded:", length(prov_codes), "codes\n")

# --- B: Detect basin-keyed vs province-keyed YAML ---------------------------
# Safety init for basin-mode objects (used later in section 5)
basin_k <- NULL
norm_basin <- NULL

# Strict binary: if ANY key matches a Longhurst code -> province mode
is_province_mode <- any(k_yaml_names %in% prov_codes)

cat("k_fast mode:", ifelse(is_province_mode, "PROVINCE-keyed", "BASIN-keyed"), "\n")

if (is_province_mode) {
  # Province-keyed: direct join
  overlap <- length(intersect(k_yaml_names, prov_codes))
  if (overlap < 0.8 * length(k_yaml_names)) {
    warning("Province mode but low overlap: ", overlap, "/", length(k_yaml_names),
            " (check the YAML keys)")
  }
  k_table <- data.table(longhurst_pr = k_yaml_names, k_fast = k_yaml_vals)
  longhurst <- dplyr::left_join(longhurst, k_table, by = "longhurst_pr")
  if (!"k_fast" %in% names(longhurst)) {
    stop("Province mode active but k_fast was not joined to longhurst (check code_col / YAML keys)")
  }
} else {
  # Basin-keyed: map province -> basin using geographic rules
  # Normalisation function for fuzzy matching of basin names
  norm_basin <- function(x) gsub("[^a-z0-9]+", "", tolower(x))

  # Build basin lookup from YAML with normalized key
  basin_k <- data.table(basin = k_yaml_names, k_fast = k_yaml_vals)
  basin_k[, basin_norm := norm_basin(basin)]
  cat("Basin k_fast values:\n")
  print(basin_k)

  # Authoritative province -> basin lookup (source: Longhurst v4 documentation)
  PROV_BASIN_LUT <- data.table(
    prov_code = c("ARCT","SARC","NADR","GFST","NASW","NATR","WTRA","ETRA","SATL",
                  "NECS","CNRY","GUIN","GUIA","NWCS","MEDI","CARB","NASE","BRAZ",
                  "FKLD","BENG","MONS","ISSG","EAFR","REDS","ARAB","INDE","SUND",
                  "NEWZ","SSTC","SANT","CHIL","CHIN","CAMR","CCAL","WARM","NPTG",
                  "NPPF","NPSG","NPTE","NPEQ","SPTG","SPPF","SPSG","SPTE","SPEQ",
                  "PEQD","ARCH","ANTA","APLR","BPLR","ALSK","AUSE","AUSW","BERS",
                  "INDW","KURO","NPSW","PNEC","PSAE","PSAW","TASM"),
    basin_raw = c("Arctic","Arctic","Atlantic","Atlantic","Atlantic","Atlantic","Atlantic",
                  "Atlantic","Atlantic","Atlantic","Atlantic","Atlantic","Atlantic","Atlantic",
                  "Mediterranean","Gulf of Mexico and Caribbean","Atlantic","Atlantic",
                  "Atlantic","Atlantic","Indian","Indian","Indian","Indian","Indian","Indian",
                  "Pacific","Pacific","Pacific","Pacific","Pacific","Pacific","Pacific",
                  "Pacific","Pacific","Pacific","Pacific","Pacific","Pacific","Pacific",
                  "Pacific","Pacific","Pacific","Pacific","Pacific","Pacific","Pacific",
                  "Antarctic","Arctic","Arctic","Arctic","Pacific","Pacific","Arctic",
                  "Indian","Pacific","Pacific","Pacific","Pacific","Pacific","Pacific")
  )
  cat("PROV_BASIN_LUT loaded:", nrow(PROV_BASIN_LUT), "provinces\n")

  # We will assign basin AFTER the spatial join with Longhurst (section 5)
  # Store basin_k for later
  cat("Basin mapping will be applied after Longhurst spatial join\n")
}

# =============================================================================
# 5) Spatial join and f_i — D
# =============================================================================
cat("\n== 5) Spatial join and f_i ==\n")

# Centroid coordinates (bottom-up, constants.R)
grid_cent <- data.table(grid_id = as.integer(sar_global$grid_id))
grid_cent[, `:=`(
  col = (grid_id - 1L) %% GRID_COLS,
  row = (grid_id - 1L) %/% GRID_COLS
)]
grid_cent[, `:=`(
  x = WORLD_XMIN + col * CELL_SIZE_M + CELL_SIZE_M / 2,
  y = WORLD_YMIN + row * CELL_SIZE_M + CELL_SIZE_M / 2
)]

sf_cent <- st_as_sf(grid_cent, coords = c("x","y"), crs = 6933)

fi_dt <- as.data.table(st_join(sf_cent, longhurst, left=TRUE))[, geometry := NULL]

fi_dt <- merge(fi_dt, sar_global, by="grid_id", all.x=TRUE)

# --- Lon/lat from sf_cent (reuse, no recalculation) --------------------------
tmp_4326 <- st_transform(sf_cent, 4326)
coords_4326 <- as.data.table(st_coordinates(tmp_4326))
setnames(coords_4326, c("lon_deg", "lat_deg"))
coords_4326[, grid_id := grid_cent$grid_id]
fi_dt <- merge(fi_dt, coords_4326, by = "grid_id", all.x = TRUE)

# --- k_fast assignment (basin mode) ------------------------------------------
if (!is_province_mode) {
  if (is.null(basin_k)) stop("basin_k not set (logic bug)")

  # Join LUT on longhurst province code
  fi_dt <- merge(fi_dt, PROV_BASIN_LUT, by.x = "longhurst_pr", by.y = "prov_code", all.x = TRUE)

  # Split Pacific N/S by centroid latitude
  fi_dt[, basin_yaml := basin_raw]
  fi_dt[basin_raw == "Pacific", basin_yaml := fifelse(lat_deg >= 0, "North Pacific", "South Pacific")]

  # Antarctic fallback: no YAML key -> route to nearest basin
  fi_dt[basin_raw == "Antarctic", basin_yaml := fifelse(lat_deg >= -60, "South Pacific", "Arctic")]

  # Fallback for cells without Longhurst province (NA or not in LUT)
  fi_dt[is.na(basin_yaml), basin_yaml := fifelse(
    lat_deg >= 0 & (lon_deg > 100 | lon_deg < -100), "North Pacific",
    fifelse(lat_deg < 0 & (lon_deg > 100 | lon_deg < -100), "South Pacific",
    fifelse(lon_deg > 20 & lon_deg <= 100, "Indian",
    "Atlantic"))
  )]

  # Remove any pre-existing k_fast column to prevent merge collision (k_fast.x/y)
  if ("k_fast" %in% names(fi_dt)) fi_dt[, k_fast := NULL]

  # Normalized join to handle YAML key variants
  fi_dt[, basin_norm := norm_basin(basin_yaml)]
  fi_dt <- merge(fi_dt, basin_k[, .(basin_norm, k_fast)], by = "basin_norm", all.x = TRUE)

  # Guard: warn if too many NA after basin join
  pct_na_k <- 100 * mean(is.na(fi_dt$k_fast))
  if (pct_na_k > 5) {
    unmatched <- setdiff(unique(fi_dt$basin_norm), basin_k$basin_norm)
    warning("k_fast NA > 5% in basin mode (", round(pct_na_k, 1), "%). Unmatched basins: ",
            paste(unmatched, collapse = ", "))
  }

  # Diagnostics
  cat("Basin assignment counts:\n")
  print(fi_dt[, .N, by = basin_yaml][order(-N)])

  # Clean up temp cols
  fi_dt[, c("basin_raw", "basin_norm", "basin_yaml") := NULL]
}
# lon/lat already merged above for both modes

# --- k_fast coverage diagnostics ---------------------------------------------
n_with_k <- sum(!is.na(fi_dt$k_fast))
cat("k_fast non-NA:", n_with_k, "/", nrow(fi_dt),
    sprintf("(%.1f%%)\n", 100 * n_with_k / nrow(fi_dt)))

# --- Load the freshness factor from the YAML ----------------------------------
if (!is.null(par$fresh_fact) && length(par$fresh_fact) > 0) {
  fresh_fact_table <- data.table(
    longhurst_pr = names(par$fresh_fact),
    fresh_fact = unlist(par$fresh_fact)
  )
  fi_dt <- merge(fi_dt, fresh_fact_table, by="longhurst_pr", all.x=TRUE)
} else {
  fi_dt[, fresh_fact := NA_real_]
}
fi_dt[is.na(fresh_fact), fresh_fact := 1.0]

# --- Final f_i -----------------------------------------------------------------
t <- 1
fi_dt[, k_used := coalesce(k_fast, 1.0)]
if (depth_mode == "mechanistic") {
  fi_dt[, `:=`(
    W_eff      = weight_eff(k_used,     slow_k,     w1, w2, w3, H3_cm, t),
    W_eff_cons = weight_eff(k_used / 2, slow_k / 2, w1, w2, w3, H3_cm, t)
  )]
  fi_dt[, p_l_eff := p_l * W_eff]
  cat("\n== G) Depth weights (mechanistic) ==\n")
  print(fi_dt[!is.na(SVR) & SVR > 0,
              .(n = .N, W_eff_med = median(W_eff), W_eff_min = min(W_eff), W_eff_max = max(W_eff)),
              by = k_used][order(k_used)])
}
fi_dt[, p_l_corr := p_l_eff * fresh_fact]
fi_dt[, p_l_corr_cons := if (depth_mode == "mechanistic") p_l * W_eff_cons * fresh_fact else p_l_corr]
fi_dt[, `:=`(
  # f_i_full: base formula (Sala et al. 2021 structure)
  # f_i = SVR × p_l_corr × p_r × [two-pool kinetics]
  # p_r = 0.87: fraction of disturbed OC that resettles in originating cell (spatial factor)
  f_i_full = SVR * p_l_corr * preserv_fact *
             ( fast_frac * (1 - exp(-k_used * t)) +
               (1 - fast_frac) * (1 - exp(-slow_k * t)) ),

  f_i_conservative = SVR * (p_l_corr_cons / 2) * preserv_fact *
                ( fast_frac * (1 - exp(-k_used / 2 * t)) +
                  (1 - fast_frac) * (1 - exp(-slow_k / 2 * t)) )
)]

# f_i_depleted: applies d_i depletion factor for chronically trawled cells (Atwood et al. 2024)
# d_i = 0.272 for cells trawled annually >10 years; d_i = 1 otherwise.
# Here applied globally as conservative default (all active cells assumed chronically disturbed).
# C_r,i = C_0,i × f_i_depleted  (Atwood et al. eq.)

# ============================ 3.4 -- scope of p_r ================================
# Why: p_r (preservation_factor, 0.87)
# is the fraction of the disturbed carbon that RESETTLES IN THE CELL. The calculation
# above applies it to 100 % of the swept volume -- yet most of that volume goes
# into the hopper.
#
# Assumption: the on-site remineralisation is
# proportional to the volume left on site. It is exact because p_r enters
# f_i LINEARLY.
#
# Two passes, and why this is NOT circular: frac depends only on SVR
# (V_swept = SVR x cell area), and this correction does not touch SVR.
# Recomputing frac on the corrected grid therefore gives exactly the same frac --
# a fixed point reached at the first pass.
#
# WITHOUT FI_FRAC_FILE, p_r applies to 100 % of the swept volume.
FI_FRAC_FILE <- Sys.getenv("FI_FRAC_FILE", "")
if (nzchar(FI_FRAC_FILE)) {
  if (!file.exists(FI_FRAC_FILE)) stop("FI_FRAC_FILE not found: ", FI_FRAC_FILE)
  cat("\n== 3.4) p_r applied to the on-site fraction only ==\n")
  cat("   table  :", FI_FRAC_FILE, "\n")
  cat("   written:", format(file.mtime(FI_FRAC_FILE), "%Y-%m-%d %H:%M"), "\n")
  fr <- data.table::fread(FI_FRAC_FILE, select = c("grid_id", "frac"))
  data.table::setnames(fr, "frac", "frac_sur_place")
  fr <- unique(fr, by = "grid_id")
  fr <- fr[!is.na(frac_sur_place)]
  if (any(fr$frac_sur_place < 0 | fr$frac_sur_place > 1))
    stop("frac outside [0,1] in ", FI_FRAC_FILE)

  fi_dt[, frac_sur_place := NA_real_]
  fi_dt[fr, on = "grid_id", frac_sur_place := i.frac_sur_place]   # join BY REFERENCE
  n_abs <- fi_dt[is.na(frac_sur_place), .N]                       # (keeps the order)
  # a cell without an assigned hopper lost nothing: everything stayed there
  fi_dt[is.na(frac_sur_place), frac_sur_place := 1]

  fi_dt[, f_i_full_avant_pr         := f_i_full]
  fi_dt[, f_i_conservative_avant_pr := f_i_conservative]
  fi_dt[, f_i_full         := f_i_full * frac_sur_place]
  fi_dt[, f_i_conservative := f_i_conservative * frac_sur_place]

  s0 <- sum(fi_dt$f_i_full_avant_pr, na.rm = TRUE)
  s1 <- sum(fi_dt$f_i_full, na.rm = TRUE)
  cat(sprintf("   matched cells: %d / %d (without frac: %d, set to 1)\n",
              nrow(fi_dt) - n_abs, nrow(fi_dt), n_abs))
  cat(sprintf("   frac_sur_place: min %.3f | median %.3f | mean %.3f | max %.3f\n",
              min(fi_dt$frac_sur_place), median(fi_dt$frac_sur_place),
              mean(fi_dt$frac_sur_place), max(fi_dt$frac_sur_place)))
  cat(sprintf("   sum f_i_full: %.6g -> %.6g  (x%.4f)\n", s0, s1,
              if (s0 > 0) s1/s0 else NA_real_))
} else {
  cat("\n== 3.4) FI_FRAC_FILE not set -> 100 % default: p_r applied to the whole SVR ==\n")
  cat("   (kept for reproducibility)\n")
}
# ========================== end of block 3.4 =====================================

fi_dt[, f_i_depleted := f_i_full * depletion_di]
fi_dt[, p_l_corr_cons := NULL]
if (depth_mode == "mechanistic") fi_dt[, c("w1", "w2", "w3", "H3_cm", "W_eff_cons") := NULL]

# CAP in [0, 1]
fi_dt[, f_i_full         := pmin(pmax(f_i_full, 0), 1)]
fi_dt[, f_i_depleted     := pmin(pmax(f_i_depleted, 0), 1)]
fi_dt[, f_i_conservative := pmin(pmax(f_i_conservative, 0), 1)]

cat("f_i computed\n")

# --- F: Diagnostics -----------------------------------------------------------
# Fallback = k_fast was NA (not just k_used==1.0 which could be a real YAML value)
fi_dt[, is_fallback := is.na(k_fast) & !is.na(SVR) & SVR > 0]
n_fallback <- fi_dt[, sum(is_fallback, na.rm=TRUE)]
n_active   <- fi_dt[, sum(SVR > 0, na.rm=TRUE)]
if (n_active == 0) n_active <- 1

cat(sprintf("   Cells with missing k_fast (fallback k_used=1.0): %d / %d (%.1f%%)\n",
            n_fallback, n_active, 100*n_fallback/n_active))
if (n_fallback / n_active > 0.05) {
  warning("Cells: more than 5% do not have a k_fast (fallback k_used=1.0)")
}

cat(sprintf("   p_l non-NA        : %d / %d (%.1f%%)\n",
            sum(!is.na(fi_dt$p_l)), nrow(fi_dt),
            100 * sum(!is.na(fi_dt$p_l)) / nrow(fi_dt)))
cat(sprintf("   f_i_full non-NA   : %d / %d (%.1f%%)\n",
            sum(!is.na(fi_dt$f_i_full)), nrow(fi_dt),
            100 * sum(!is.na(fi_dt$f_i_full)) / nrow(fi_dt)))

cat("\nTop-10 f_i_full (with lon/lat EPSG:4326):\n")
print(head(fi_dt[order(-f_i_full),
                 .(grid_id, SAR, SVR, p_l_eff, p_l_corr, k_used, longhurst_pr,
                   f_i_full, lon_deg, lat_deg)], 10))

# =============================================================================
# 6) Final save
# =============================================================================
cat("\n== 6) Saving ==\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Parquet
out_parquet <- file.path(MERGE_OUT_DIR, sprintf("fi_grid_%s.parquet", timestamp))
write_parquet(fi_dt, out_parquet)
cat("Parquet:", basename(out_parquet), "\n")

# RDS
out_rds <- file.path(MERGE_OUT_DIR, sprintf("fi_grid_%s.rds", timestamp))
saveRDS(fi_dt, out_rds)
cat("RDS     :", basename(out_rds), "\n")

# =============================================================================
# 7) GeoTIFF — E: top-down row inversion
# =============================================================================
CREATE_GEOTIFF <- as.logical(Sys.getenv("MAKE_TIFF", "TRUE"))
if (is.na(CREATE_GEOTIFF)) CREATE_GEOTIFF <- TRUE

if(CREATE_GEOTIFF) {
  cat("\n== 7) GeoTIFF raster ==\n")

  tryCatch({
    raster_data <- fi_dt[, .(grid_id, f_i_full)]
    raster_data <- raster_data[!is.na(f_i_full)]

    if(nrow(raster_data) > 0) {
      raster_data[, `:=`(
        col = (grid_id - 1L) %% GRID_COLS,
        row = (grid_id - 1L) %/% GRID_COLS
      )]

      r <- rast(
        nrows = GRID_ROWS, ncols = GRID_COLS,
        xmin = WORLD_XMIN, xmax = WORLD_XMAX,
        ymin = WORLD_YMIN, ymax = WORLD_YMAX,
        crs = "EPSG:6933"
      )
      names(r) <- "f_i_full"

      vals <- rep(NA_real_, ncell(r))

      # Terra indexes top-down, our rows are bottom-up
      inds <- (GRID_ROWS - 1L - raster_data$row) * ncol(r) + raster_data$col + 1L

      vals[inds] <- raster_data$f_i_full
      values(r) <- vals

      out_tif <- file.path(MERGE_OUT_DIR, sprintf("fi_grid_%s.tif", timestamp))
      writeRaster(
        r, out_tif,
        datatype  = "FLT4S",
        overwrite = TRUE,
        gdal      = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
      )

      cat("GeoTIFF created:", basename(out_tif), "\n")
      cat("Raster statistics:", nrow(raster_data), "non-NA cells out of", ncell(r), "pixels\n")
    } else {
      cat("No valid f_i data to create the raster\n")
    }
  }, error = function(e) {
    cat("GeoTIFF creation failed:", e$message, "\n")
    cat("The Parquet and RDS files are available\n")
  })
} else {
  cat("\nGeoTIFF creation disabled\n")
}

cat("\nStep 5 merge completed successfully!\n")
cat("Files created:\n")
cat("   -", basename(out_parquet), "\n")
cat("   -", basename(out_rds), "\n")
if(CREATE_GEOTIFF) cat("   -", sprintf("fi_grid_%s.tif", timestamp), "\n")
