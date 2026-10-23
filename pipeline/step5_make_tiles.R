#!/usr/bin/env Rscript
# ────────────────────────────────────────────────────────────────────────────────
# STEP-5  ─  GLOBAL TILE GENERATION (pipeline V6)
# Splits the world extent into square tiles for modular processing
# ────────────────────────────────────────────────────────────────────────────────

library(sf)

# Shared constants
script_dir <- getwd()
source(file.path(script_dir, "constants.R"))

cat("🔧 Generating the global tiles for Step-5\n")

# Tiling parameters (from constants.R)
cell_km   <- CELL_KM   # tile size in km
crs_eq    <- CRS_EQUIVALENT   # equal-area CRS (metres)

# World extent in metres (EPSG:6933) - from constants.R
world_bb  <- st_as_sfc(st_bbox(c(xmin=WORLD_XMIN, ymin=WORLD_YMIN,
                                 xmax=WORLD_XMAX, ymax=WORLD_YMAX),
                               crs=st_crs(crs_eq)))

# Print the extent
bb <- st_bbox(world_bb)
cat("📐 Global extent:\n")
cat("   - X: ", format(bb["xmin"], scientific=FALSE), " to ", 
    format(bb["xmax"], scientific=FALSE), " m\n")
cat("   - Y: ", format(bb["ymin"], scientific=FALSE), " to ", 
    format(bb["ymax"], scientific=FALSE), " m\n")

# Compute exactly the number of tiles needed
# rather than letting st_make_grid truncate
nx <- ceiling((WORLD_XMAX - WORLD_XMIN) / (cell_km * 1000))   # 35
ny <- ceiling((WORLD_YMAX - WORLD_YMIN) / (cell_km * 1000))   # 18

cat("🔲 Generating the tile grid...\n")
cat("   - Number of columns (X):", nx, "\n")
cat("   - Number of rows (Y):", ny, "\n")
cat("   - Expected total:", nx * ny, "tiles\n")

# Grid with fixed dimensions
tiles_sf <- sf::st_make_grid(
  offset   = c(WORLD_XMIN, WORLD_YMIN),
  cellsize = cell_km * 1000,
  n        = c(nx, ny),        # ← Forces exactly 35×18
  crs      = crs_eq,
  what     = "polygons"
) %>%
  sf::st_sf()

# ------------------------------------------------------------------
# Indexing: pure arithmetic order (no rounding)
# ------------------------------------------------------------------
tiles_sf$tile_id <- seq_len(nx * ny)              # 1 … 630
tiles_sf$row     <- (tiles_sf$tile_id - 1) %/% nx # 0 … 17
tiles_sf$col     <- (tiles_sf$tile_id - 1) %%  nx # 0 … 34

# Check: number of columns / rows (must be nx, ny)
n_cols <- nx
n_rows <- ny

cat("   - Columns generated:", n_cols, "(expected:", nx, ")\n")
cat("   - Rows generated:", n_rows, "(expected:", ny, ")\n")

# ------------------------------------------------------------------

# Tile metadata
tiles_sf$tile_size_km <- cell_km
tiles_sf$crs <- "EPSG:6933"

# Statistics (corrected area)
n_tiles <- nrow(tiles_sf)
surface_km2 <- n_tiles * cell_km^2
cat("✅ Tiles generated:", n_tiles, "\n")
cat("   - Size per tile:", cell_km, "×", cell_km, "km\n")
cat("   - Total area:", format(surface_km2, big.mark=","), "km²\n")

# Check that all IDs are present
tile_ids <- sort(tiles_sf$tile_id)
expected_ids <- 1:n_tiles
missing_ids <- setdiff(expected_ids, tile_ids)
if(length(missing_ids) > 0) {
  cat("⚠️  Missing IDs:", paste(missing_ids, collapse=", "), "\n")
} else {
  cat("✅ All IDs from 1 to", n_tiles, "are present\n")
}

# Save with a dynamic name
output_file <- sprintf(file.path(Sys.getenv("DFP_DATA"), "output_V6/tiles_%dkm.gpkg"), cell_km)
cat("💾 Saving to:", output_file, "\n")

st_write(tiles_sf, output_file, delete_dsn = TRUE, quiet = TRUE)

if(file.exists(output_file)) {
  file_info <- file.info(output_file)
  cat("✅ File created:", basename(output_file), "\n")
  cat("   - Size:", format(file_info$size, scientific=FALSE), "bytes\n")
  cat("   - Date:", format(file_info$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
} else {
  stop("❌ Error while creating the file")
}

# Print the first tiles for checking
cat("\n📋 Tile preview:\n")
print(head(tiles_sf[, c("tile_id", "tile_size_km")]))

cat("\n🎯 Next steps:\n")
cat("   1. Check the number of tiles: ogrinfo -ro -so", output_file, "tiles | grep 'Feature Count'\n")
cat("   2. Launch the array job: sbatch --array=1-", n_tiles, " step5_batch.sbatch (with DREDGE_LABEL set)\n")
cat("   3. Monitor: squeue -u $USER | grep step5_tile\n")

cat("\n✅ Tile generation done!\n") 