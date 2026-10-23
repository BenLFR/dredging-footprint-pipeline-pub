user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(unique(c(user_lib, .libPaths())))
#!/usr/bin/env Rscript
# ────────────────────────────────────────────────────────────────────────────────
# STEP-6  ─  REMINERALISED CARBON (Cri)
#  - Grid identical to constants.R (6933, 34735×14685)
#  - Uses x/y/col/row of fi_grid when available (no 36k×18k canvas)
#  - Physical caps: f_i ∈ [0,1], C_ri ≤ C0i
# ────────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(data.table); library(terra); library(lubridate); library(arrow)
})
sf::sf_use_s2(FALSE)

# constants.R next to this script
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- file.path(Sys.getenv("DFP_REPO"), "pipeline")
# constants.R is looked for next to this script, then in the parent directory --
# never a local copy: a duplicated constants.R would diverge.
.cst <- file.path(script_dir, "constants.R")
if (!file.exists(.cst)) .cst <- file.path(dirname(script_dir), "constants.R")
if (!file.exists(.cst)) stop("constants.R not found from ", script_dir)
source(.cst)

# Terra temp directory
tmp_terra <- file.path(Sys.getenv("DFP_DATA"), "tmp_terra")
dir.create(tmp_terra, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = tmp_terra)

cat("--- Step 6: Cri ---\n")

# f_i file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Explicit fallback only: without an argument, silently reading the most recent
  # fi_grid of $DFP_DATA/output_V6/ (the PUBLISHED grid) would give a
  # normal-looking C_ri computed on the wrong grid.
  if (!nzchar(Sys.getenv("CRI_ALLOW_FALLBACK"))) {
    stop("No f_i file passed as argument.\n",
         "  Pass the path explicitly:  Rscript step6_calculate_cri_corr.R <fi_grid.parquet>\n",
         "  To use the fallback (most recent fi_grid of $DFP_DATA/output_V6/),\n",
         "  rerun with CRI_ALLOW_FALLBACK=1 -- knowing that it is the PUBLISHED grid.")
  }
  cat("⚠️  CRI_ALLOW_FALLBACK: reading the most recent fi_grid of $DFP_DATA/output_V6/\n")
  fi_files <- list.files(file.path(Sys.getenv("DFP_DATA"), "output_V6/"), pattern="^fi_grid_.*\\.(parquet|rds)$", full.names=TRUE)
  if(length(fi_files) == 0) stop("No f_i file found in $DFP_DATA/output_V6/")
  fi_path <- fi_files[which.max(file.info(fi_files)$mtime)]
} else fi_path <- args[1]

# Destination resolved and GUARDED HERE, early, not at the writing line:
# a guard placed after 100 lines of computation can only be tested by paying for the computation.
CRI_OUT_DIR <- path.expand(Sys.getenv("CRI_OUT_DIR", file.path(Sys.getenv("DFP_DATA"), "output_V6_corr/")))
CRI_TAG     <- Sys.getenv("CRI_TAG", "")
if (grepl("^/home/[^/]+/scratch/output_V6/?$", sub("/+$", "", CRI_OUT_DIR))) {
  stop("CRI_OUT_DIR points to the PUBLISHED outputs (", CRI_OUT_DIR, ") -- WRITE REFUSED.")
}
dir.create(CRI_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("\U0001F4C1 destination:", CRI_OUT_DIR, if (nzchar(CRI_TAG)) paste0("(tag ", CRI_TAG, ")") else "", "\n")
cat("✔️ f_i file used:", basename(fi_path), "\n")

fi_dt <- if (grepl("\\.parquet$", fi_path)) arrow::read_parquet(fi_path) else readRDS(fi_path)
setDT(fi_dt)
stopifnot("grid_id" %in% names(fi_dt))

# Target canvas aligned with Step-5 (created early for xyFromCell)
template_raster <- terra::rast(nrows=GRID_ROWS, ncols=GRID_COLS, crs="EPSG:6933",
                               xmin=WORLD_XMIN, xmax=WORLD_XMAX, ymin=WORLD_YMIN, ymax=WORLD_YMAX)

# Rebuild x/y via terra (source of truth for extraction + GeoTIFF)
# grid_id is bottom-up (row=0=south), terra is top-down (cell 1=NW)
# Conversion: terra_cell = (GRID_ROWS - 1 - row) * GRID_COLS + col + 1
fi_dt[, col := as.integer((grid_id - 1L) %% GRID_COLS)]
fi_dt[, row := as.integer((grid_id - 1L) %/% GRID_COLS)]
fi_dt[, terra_cell := (GRID_ROWS - 1L - row) * GRID_COLS + col + 1L]

if (!("x" %in% names(fi_dt) && "y" %in% names(fi_dt))) {
  xy <- terra::xyFromCell(template_raster, fi_dt$terra_cell)
  fi_dt[, `:=`(x = xy[,1], y = xy[,2])]
  cat("✅ x/y coordinates rebuilt via xyFromCell\n")
}

# Sanity check: roundtrip terra_cell -> xy -> terra_cell
n_check <- min(2000L, nrow(fi_dt))
check_idx <- fi_dt$terra_cell[sample.int(nrow(fi_dt), n_check)]
check_xy <- terra::xyFromCell(template_raster, check_idx)
check_back <- terra::cellFromXY(template_raster, check_xy)
if (!all(check_back == check_idx)) {
  stop("FATAL: grid_id <-> Terra cell roundtrip mismatch!")
}
cat("✅ Sanity check OK:", n_check, "cells checked\n")

# Cap f_i for safety + NA -> 0 (no dredging = zero impact)
if ("f_i_full" %in% names(fi_dt)) {
  fi_dt[is.na(f_i_full), f_i_full := 0]
  fi_dt[, f_i_full := pmin(pmax(f_i_full, 0), 1)]
}
if ("f_i_conservative" %in% names(fi_dt)) {
  fi_dt[is.na(f_i_conservative), f_i_conservative := 0]
  fi_dt[, f_i_conservative := pmin(pmax(f_i_conservative, 0), 1)]
}

## 2) Loading C0i (Atwood) and reprojecting to 6933 if needed ---------------
cat("🔄 Loading the C0i rasters…\n")
carbon_dir <- path.expand(file.path(Sys.getenv("DFP_DATA"), "configuration/atwood_carbon_full"))
carbon_files <- list.files(carbon_dir, pattern = "\\.tif$", full.names = TRUE)
if(length(carbon_files) == 0) stop("❌ No C0i TIF in ", carbon_dir)

carbon_rasters <- list()
for(file in carbon_files) {
  nm <- tools::file_path_sans_ext(basename(file))
  nm <- tolower(gsub(" ", "_", nm))
  r  <- terra::rast(file)
  if (!terra::compareGeom(r, template_raster, stopOnError = FALSE)) {
    cat("🔄 Reprojection:", basename(file), "→ EPSG:6933 (template)\n")
    r <- terra::project(r, template_raster)
  }
  carbon_rasters[[nm]] <- r
}

coords <- fi_dt[, .(x, y)]
for(nm in names(carbon_rasters)) {
  col_name <- paste0("C0i_", nm)
  fi_dt[[col_name]] <- terra::extract(carbon_rasters[[nm]], coords, ID=FALSE)[[1]]
  fi_dt[is.na(get(col_name)), (col_name) := 0]
  cat("✔️", col_name, "extracted\n")
}

# Central/bound columns
stopifnot("C0i_mean_carbon_stock" %in% names(fi_dt))
fi_dt[, C0i := C0i_mean_carbon_stock]
if ("C0i_global_error_lower_bound" %in% names(fi_dt)) fi_dt[, C0i_lower := C0i_global_error_lower_bound]
if ("C0i_global_error_upper_bound" %in% names(fi_dt)) fi_dt[, C0i_upper := C0i_global_error_upper_bound]

## 3) Depletion factor di ----------------------------------------------------
# d_i = 1 for every cell: a CHOICE OF THIS STUDY.
# There is no multi-year DREDGING history: each cell is treated as disturbed
# for the first time. In Atwood et al. (2024), d_i = 0.272 applies to seabeds
# TRAWLED for more than 10 years; omitting it gives an upper bound of the flux.
fi_dt[, di := 1.0]
cat("📌 di = 1 for every cell (study choice, no dredging history)\n")

## 4) C_ri (caps included) ----------------------------------------------------
cat("🔄 C_ri…\n")
fi_dt[, C_ri := pmin(C0i * f_i_full * di, C0i)]
if ("C0i_lower" %in% names(fi_dt)) fi_dt[, C_ri_lower := pmin(C0i_lower * f_i_full * di, C0i_lower)]
if ("C0i_upper" %in% names(fi_dt)) fi_dt[, C_ri_upper := pmin(C0i_upper * f_i_full * di, C0i_upper)]
fi_dt[, C_ri_conservative := pmin(C0i * f_i_conservative * di, C0i)]

## 5) Saving -------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# CRI_TAG prevents two labels run in the same minute from overwriting each other
# (they would share the time stamp). CRI_OUT_DIR is resolved and guarded above.
item_name <- if (nzchar(CRI_TAG)) paste0("cri_final_", CRI_TAG, "_", timestamp) else paste0("cri_final_", timestamp)
output_prefix <- file.path(CRI_OUT_DIR, item_name)
cat("\U0001F4BE output:", output_prefix, "\n")

saveRDS(fi_dt, paste0(output_prefix, ".rds"))
cat("✔️ RDS written\n")
if(requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(fi_dt, paste0(output_prefix, ".parquet"))
  cat("✔️ Parquet written\n")
}

# GeoTIFF (optional but useful)
write_cri_raster <- function(values, filename) {
  g <- template_raster
  g[] <- NA_real_
  g[fi_dt$terra_cell] <- values
  writeRaster(g, filename, datatype = "FLT4S", overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))
  cat("✔️", basename(filename), "written\n")
}
cat("🗺️  GeoTIFF…\n")
write_cri_raster(fi_dt$C_ri, paste0(output_prefix, ".tif"))
if ("C_ri_lower" %in% names(fi_dt)) write_cri_raster(fi_dt$C_ri_lower, paste0(output_prefix, "_lower.tif"))
if ("C_ri_upper" %in% names(fi_dt)) write_cri_raster(fi_dt$C_ri_upper, paste0(output_prefix, "_upper.tif"))

cat("✅ Step 6 done\n")
