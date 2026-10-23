# ────────────────────────────────────────────────────────────────────────────────
# STEP-4  ─  LITHOLOGY via HubOcean/dbSEABED (pipeline V6)
#
# dbSEABED-compliant, with extraction per unique cell + continuous pl_base
#
# Main features:
#   - NAflag=-99 when loading the rasters (clean dbSEABED no-data)
#   - EXTRACTION PER UNIQUE CELL: N_cells << N_pings
#   - Type-specific normalization: membership 0-100, hardsoft signed [-1,+1]
#   - Hard-bottom index H = pmax(rock_prob, hard_prob) for continuous blending
#   - QC per unique cell (not biased by AIS density)
#   - terraOptions tuned for HPC
#
# Refs: USGS dbSEABED docs, rspatial.github.io/terra
# ────────────────────────────────────────────────────────────────────────────────

## 0.  Packages -----------------------------------------------------------------
pkgs <- c("data.table", "terra")

safe_library <- function(pkg) {
  tryCatch(
    {
      library(pkg, character.only = TRUE)
      cat("  Package", pkg, "loaded\n")
    },
    error = function(e) {
      stop(
        "Package '", pkg, "' missing: ", e$message,
        "\n   Load the appropriate R module or install the package"
      )
    }
  )
}

cat("Loading packages...\n")
invisible(sapply(pkgs, safe_library))
cat("All packages loaded successfully\n\n")

## 0.5 terra configuration for HPC --------------------------------------------
tmp_terra <- file.path(Sys.getenv("DFP_DATA"), "tmp_terra")
dir.create(tmp_terra, showWarnings = FALSE, recursive = TRUE)
terraOptions(memfrac = 0.8, tempdir = tmp_terra, progress = 0)
cat("terraOptions: memfrac=0.8, tempdir=", tmp_terra, "\n\n")

## 1.  I/O paths ---------------------------------------------------------------
home_dir <- path.expand("~")
cache_dir <- file.path(Sys.getenv("DFP_DATA"), "hubocean_cache")
PUBLISHED_DIR <- path.expand(file.path(Sys.getenv("DFP_DATA"), "output_V6"))
# Separate input and output directories: a new file written into the published
# directory would otherwise be picked up, silently, by every later step5.
PIPE_IN_DIR  <- path.expand(Sys.getenv("PIPE_IN_DIR", PUBLISHED_DIR))
PIPE_OUT_DIR <- Sys.getenv("PIPE_OUT_DIR", "")
if (!nzchar(PIPE_OUT_DIR)) {
  stop("PIPE_OUT_DIR is not set. Declare it explicitly: a silent ",
       "default would write yr AIS_with_lithology_* into the published outputs.")
}
PIPE_OUT_DIR <- path.expand(PIPE_OUT_DIR)
if (identical(normalizePath(PIPE_OUT_DIR, mustWork = FALSE),
              normalizePath(PUBLISHED_DIR, mustWork = FALSE)) &&
    !nzchar(Sys.getenv("PIPE_ALLOW_PUBLISHED"))) {
  stop("Refusing to write into the PUBLISHED outputs (", PUBLISHED_DIR, "). ",
       "Choose another PIPE_OUT_DIR, or PIPE_ALLOW_PUBLISHED=1 to force.")
}
dir.create(PIPE_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("📌 PIPE_IN_DIR  =", PIPE_IN_DIR,  "\n")
cat("📌 PIPE_OUT_DIR =", PIPE_OUT_DIR, "\n")
config_dir <- file.path(Sys.getenv("DFP_DATA"), "configuration/lithology")

# HubOcean cache check
if (!dir.exists(cache_dir)) {
  stop(
    "HubOcean cache not found: ", cache_dir,
    "\nRun first: python prefetch_hubocean_stac.py"
  )
}

# Raster detection (supports a subdirectory structure)
find_rasters_in_subdir <- function(subdir) {
  subdir_path <- file.path(cache_dir, subdir)
  if (!dir.exists(subdir_path)) {
    pattern <- paste0("^", subdir, "__.*\\.(tif|tiff)$")
    files <- list.files(cache_dir,
      pattern = pattern,
      full.names = TRUE, ignore.case = TRUE
    )
    return(files)
  }
  list.files(subdir_path,
    pattern = "\\.(tif|tiff)$",
    full.names = TRUE, ignore.case = TRUE
  )
}

find_single_raster_in_subdir <- function(subdir) {
  files <- find_rasters_in_subdir(subdir)
  if (length(files) == 0) {
    return(NULL)
  }
  files[which.max(file.mtime(files))]
}

## 2.  Loading the rasters -----------------------------------------------------
cat("Detecting the dbSEABED rasters...\n")

tex_files <- find_rasters_in_subdir("texture")
cat("  texture files: ", length(tex_files), " found\n")

# Detection by file name (HubOcean: mud, snd, gvl)
mud_f <- tex_files[grepl("mud|silt|clay", basename(tex_files), ignore.case = TRUE)]
sand_f <- tex_files[grepl("sand|snd", basename(tex_files), ignore.case = TRUE)]
gravel_f <- tex_files[grepl("gravel|gvl|coarse", basename(tex_files), ignore.case = TRUE)]

# Texture stack
# Pick the most recent files (not [1])
pick_latest <- function(x) x[which.max(file.mtime(x))]

if (length(mud_f) > 0 && length(sand_f) > 0 && length(gravel_f) > 0) {
  cat("  -> Mode: separate files (mud/sand/gravel)\n")
  r_texture <- rast(c(pick_latest(mud_f), pick_latest(sand_f), pick_latest(gravel_f)))
  names(r_texture) <- c("mud", "sand", "gravel")
} else if (length(tex_files) > 0) {
  cat("  -> Mode: combined file(s)\n")
  r_texture <- rast(tex_files)
} else {
  stop("No texture raster found in ", cache_dir)
}

# --- Robust renaming of the texture bands ---
nm <- names(r_texture)
if (!all(c("mud", "sand", "gravel") %in% nm)) {
  cat("  [P0-1] Renaming texture bands...\n")
  mud_idx <- grep("mud|silt|clay", nm, ignore.case = TRUE)
  sand_idx <- grep("sand|snd", nm, ignore.case = TRUE)
  gravel_idx <- grep("gravel|gvl|coarse", nm, ignore.case = TRUE)

  # fallback if the patterns are missing and there are 3 bands -> assume order 1:3
  if (length(mud_idx) == 0 && length(sand_idx) == 0 && length(gravel_idx) == 0 && nlyr(r_texture) >= 3) {
    names(r_texture)[1:3] <- c("mud", "sand", "gravel")
    cat("    -> fallback order 1:3 = mud/sand/gravel\n")
  } else {
    if (length(mud_idx) > 0) names(r_texture)[mud_idx[1]] <- "mud"
    if (length(sand_idx) > 0) names(r_texture)[sand_idx[1]] <- "sand"
    if (length(gravel_idx) > 0) names(r_texture)[gravel_idx[1]] <- "gravel"
    cat("    -> pattern match: ", paste(names(r_texture), collapse = ", "), "\n")
  }
}

cat(
  "  texture stack: ", nlyr(r_texture), " band(s), res=",
  paste(round(res(r_texture), 4), collapse = "x"), "\n"
)
cat("  band names: ", paste(names(r_texture), collapse = ", "), "\n")

# Hard/soft and Rock (optional)
hardsoft_path <- find_single_raster_in_subdir("hard_soft")
r_hardsoft <- if (!is.null(hardsoft_path)) rast(hardsoft_path) else NULL

rock_path <- find_single_raster_in_subdir("rock")
r_rock <- if (!is.null(rock_path)) rast(rock_path) else NULL

if (!is.null(r_hardsoft)) cat("  hard_soft: loaded\n")
if (!is.null(r_rock)) cat("  rock: loaded\n")

# --- Align hard_soft/rock on the texture grid ---
ref <- r_texture[[1]]

align_to_ref <- function(r, ref, name) {
  if (is.null(r)) {
    return(NULL)
  }
  if (!compareGeom(r, ref, stopOnError = FALSE)) {
    cat("  [P0-2] ", name, ": resample -> texture grid (method=near)\n")
    r <- resample(r, ref, method = "near") # near for indices/membership
  } else {
    cat("  [P0-2] ", name, ": geometry OK\n")
  }
  r
}

r_hardsoft <- align_to_ref(r_hardsoft, ref, "hard_soft")
r_rock <- align_to_ref(r_rock, ref, "rock")
rm(ref)

## 2.5 SINGLE STACK (performance) ---------------------------------------------
cat("\nBuilding a unified stack for a single extraction...\n")

r_stack <- r_texture
layer_names <- names(r_texture)

if (!is.null(r_hardsoft)) {
  r_stack <- c(r_stack, r_hardsoft)
  layer_names <- c(layer_names, "hardness")
}
if (!is.null(r_rock)) {
  r_stack <- c(r_stack, r_rock)
  layer_names <- c(layer_names, "rock")
}
names(r_stack) <- layer_names

# dbSEABED NAflag (avoids propagating -99 everywhere)
try(NAflag(r_stack) <- -99, silent = TRUE)

cat("  Stack final: ", nlyr(r_stack), " layers (", paste(layer_names, collapse = ", "), ")\n")
cat("  NAflag set to -99 for dbSEABED no-data handling\n")

## 3.  Step 3 AIS file --------------------------------------------------------
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
ais_path <- choose_input(PIPE_IN_DIR, "AIS_data_core_preprocessed_V6_.*\\.rds$",
                           Sys.getenv("PIPE_AIS_FILE", ""), "step3", rule = "mtime")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(PIPE_OUT_DIR, sprintf("AIS_with_lithology_%s.rds", timestamp))

cat("\nAIS source file: ", basename(ais_path), "\n")
cat("Planned output:     ", basename(out_file), "\n\n")

## 4.  Lithology -> pl_base lookup table --------------------------------------
pl_lookup_file <- file.path(config_dir, "pl_lookup_v1.csv")

if (file.exists(pl_lookup_file)) {
  cat("Loading pl_base_lookup from: ", basename(pl_lookup_file), "\n")
  pl_base_lookup <- fread(pl_lookup_file)
} else {
  cat("Using the default (hard-coded) pl_base_lookup\n")
  pl_base_lookup <- data.table(
    lithologie_canon = c(
      "mud", "sandy_mud", "muddy_sand", "sand",
      "gravelly_sand", "sandy_gravel", "gravel",
      "rock", "hard_bottom", "mixed", "unknown"
    ),
    pl_base = c(
      0.92, # mud - high penetration
      0.85, # sandy_mud
      0.70, # muddy_sand
      0.55, # sand - medium penetration
      0.45, # gravelly_sand
      0.35, # sandy_gravel
      0.25, # gravel - low penetration
      0.05, # rock - very low
      0.10, # hard_bottom
      0.50, # mixed - medium
      0.50 # unknown - default
    )
  )
}
setkey(pl_base_lookup, lithologie_canon)

## 5.  dbSEABED utility functions ---------------------------------------------

# Normalise membership 0-100 -> 0-1 (rock, mud, sand, gravel)
norm_membership_0_100 <- function(x) {
  n_nodata <- sum(x <= -90, na.rm = TRUE)
  x[x <= -90] <- NA_real_

  # If values > 1, the scale is 0-100 -> divide by 100
  if (any(x > 1, na.rm = TRUE)) {
    x <- x / 100
  }

  # Clamp 0-1
  x <- pmin(pmax(x, 0), 1)
  attr(x, "n_nodata") <- n_nodata
  x
}

# Signed hard_soft (-1 to +1) -> kept as is
norm_hardsoft_signed <- function(x) {
  x[x <= -90] <- NA_real_
  # Clamp -1 to +1 (keeps the sign for QA)
  pmin(pmax(x, -1), 1)
}

# Signed hard_soft -> probability 0-1
# -1 (soft) -> 0, +1 (hard) -> 1
hardsoft_to_prob <- function(x_signed) {
  ifelse(is.na(x_signed), NA_real_, pmin(pmax((x_signed + 1) / 2, 0), 1))
}

# Texture closure (mud+sand+gravel -> sum=1) for classification
close_texture <- function(mud, sand, gravel) {
  tot <- mud + sand + gravel
  ok <- !is.na(tot) & tot > 0

  mudC <- sandC <- gravelC <- rep(NA_real_, length(tot))
  mudC[ok] <- mud[ok] / tot[ok]
  sandC[ok] <- sand[ok] / tot[ok]
  gravelC[ok] <- gravel[ok] / tot[ok]

  list(mud = mudC, sand = sandC, gravel = gravelC)
}

## 6.  Lithology classification (dbSEABED-compliant) ---------------------------
# Uses the hard-bottom index H = pmax(rock, hardness)
# Rock confirmed only if H > 0.7 AND the texture is not soft

canon_litho_v2 <- function(mud, sand, gravel, hard, rock) {
  n <- length(mud)

  # Continuous hard-bottom index
  H <- pmax(rock, hard, na.rm = TRUE)
  H[!is.finite(H)] <- NA_real_

  # Masks
  has_frac <- !(is.na(mud) & is.na(sand) & is.na(gravel))
  has_H <- !is.na(H)

  # Copies with NA->0 for comparisons
  mud2 <- fifelse(is.na(mud), 0, mud)
  sand2 <- fifelse(is.na(sand), 0, sand)
  gravel2 <- fifelse(is.na(gravel), 0, gravel)
  H2 <- fifelse(is.na(H), 0, H)
  rock2 <- fifelse(is.na(rock), 0, rock)

  result <- rep("unknown", n)

  # 1. CONFIRMED rock (strict dbSEABED criteria)
  # H > 0.7 AND (rock > 0.6 OR hardness > 0.7) AND texture not soft
  texture_not_soft <- (mud2 < 0.3) & (sand2 < 0.5)
  rock_confirmed <- has_H & H2 > 0.7 & rock2 > 0.5 & texture_not_soft
  result[rock_confirmed] <- "rock"

  # 2. Hard bottom (high H but rock not dominant)
  idx_hard <- result == "unknown" & has_H & H2 > 0.75 & rock2 < 0.5
  result[idx_hard] <- "hard_bottom"

  # 3. Classification by texture fractions
  to_classify <- result == "unknown" & has_frac

  if (any(to_classify)) {
    # Gravel dominant (>50%)
    idx <- to_classify & gravel2 > 0.5
    result[idx] <- "gravel"

    idx <- to_classify & result == "unknown" & gravel2 > 0.3 & sand2 > mud2
    result[idx] <- "sandy_gravel"

    idx <- to_classify & result == "unknown" & gravel2 > 0.3 & mud2 >= sand2
    result[idx] <- "gravel"

    # Moderate gravel (15-30%)
    idx <- to_classify & result == "unknown" & gravel2 > 0.15 & sand2 > mud2
    result[idx] <- "gravelly_sand"

    idx <- to_classify & result == "unknown" & gravel2 > 0.15 & mud2 >= sand2
    result[idx] <- "mixed"

    # Mud dominant (>70%)
    idx <- to_classify & result == "unknown" & mud2 > 0.7
    result[idx] <- "mud"

    # Sand dominant (>70%)
    idx <- to_classify & result == "unknown" & sand2 > 0.7
    result[idx] <- "sand"

    # Mud-sand mixtures
    idx <- to_classify & result == "unknown" & mud2 > 0.4 & sand2 > 0.2 & mud2 > sand2
    result[idx] <- "sandy_mud"

    idx <- to_classify & result == "unknown" & sand2 > 0.4 & mud2 > 0.2 & sand2 > mud2
    result[idx] <- "muddy_sand"

    # Remaining: classification by maximum
    idx <- to_classify & result == "unknown" & (mud2 + sand2 + gravel2) > 0
    result[idx & mud2 >= sand2 & mud2 >= gravel2] <- "mud"
    result[idx & sand2 > mud2 & sand2 >= gravel2] <- "sand"
    result[idx & gravel2 > mud2 & gravel2 > sand2] <- "gravel"
  }

  list(litho = result, H = H)
}

## 7.  Reading the AIS data ----------------------------------------------------
cat("Reading the AIS data...\n")
ais_dt <- readRDS(ais_path)
setDT(ais_dt)

n_ais <- nrow(ais_dt)
cat("  AIS pings loaded: ", format(n_ais, big.mark = " "), "\n")

if (!all(c("Lon", "Lat") %in% names(ais_dt))) {
  stop("Lon/Lat columns missing from the AIS data")
}

## 8. EXTRACTION PER UNIQUE CELL (then remap to pings) ------------------------
# Extracts once per raster cell, then remaps to all pings
# (N_cells << N_pings on dense AIS)

cat("\nRaster extraction per unique cell (spatial dedup)...\n")

# 8.1 cell_id of each ping (chunked, RAM-safe)
ais_dt[, cell_id := NA_integer_]

chunk_size <- 2e6
n_chunks <- ceiling(n_ais / chunk_size)
idx_list <- split(seq_len(n_ais), ceiling(seq_len(n_ais) / chunk_size))

for (k in seq_along(idx_list)) {
  ii <- idx_list[[k]]
  cat("  cell_id chunk ", k, "/", n_chunks, " (", length(ii), " pings)... ")

  xy <- as.matrix(ais_dt[ii, .(Lon, Lat)])
  ais_dt[ii, cell_id := cellFromXY(r_stack[[1]], xy)]

  rm(xy)
  gc(verbose = FALSE)
  cat("OK\n")
}

# 8.2 List of unique cells (much smaller than n_ais)
# no sort(): faster on large datasets
cells_unique <- unique(ais_dt$cell_id)
cells_unique <- cells_unique[!is.na(cells_unique)]

n_cells <- length(cells_unique)
cat(
  "  Unique cells to extract: ", format(n_cells, big.mark = " "),
  " (vs ", format(n_ais, big.mark = " "), " pings)\n"
)
cat("  Ratio compression: ", round(n_ais / n_cells, 1), "x\n")

if (n_cells == 0) stop("No raster cell associated with the pings (check extent/CRS).")

# 8.3 Raster values extracted ONCE per cell (chunked)
cell_chunk <- 1e6 # adjustable
n_cell_chunks <- ceiling(n_cells / cell_chunk)
cell_idx_list <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / cell_chunk))

cell_vals_list <- vector("list", length(cell_idx_list))

for (k in seq_along(cell_idx_list)) {
  jj <- cell_idx_list[[k]]
  cell_ids_k <- cells_unique[jj]

  cat("  extract cells ", k, "/", n_cell_chunks, " (", length(cell_ids_k), " cells)... ")

  xy_cell <- xyFromCell(r_stack[[1]], cell_ids_k)
  v <- extract(r_stack, xy_cell, method = "simple")
  # fallback: remove the ID column if present (terra version compatibility)
  if ("ID" %in% names(v)) v$ID <- NULL

  # Secure type + cell_id column
  cell_dt_k <- data.table(cell_id = cell_ids_k)
  cell_dt_k <- cbind(cell_dt_k, as.data.table(v))

  cell_vals_list[[k]] <- cell_dt_k

  rm(cell_ids_k, xy_cell, v, cell_dt_k)
  gc(verbose = FALSE)
  cat("OK\n")
}

cell_vals <- rbindlist(cell_vals_list, use.names = TRUE, fill = TRUE)
rm(cell_vals_list)
gc(verbose = FALSE)

cat("  Values extracted for ", nrow(cell_vals), " unique cells\n")

# 8.4 Remapping cell_id -> raster values (data.table JOIN)
# rename + mget rather than get(paste0("i.",x)), which is fragile under NSE
setkey(cell_vals, cell_id)

# Rename the cell_vals columns to the target AIS columns
rename_map <- c(
  mud = "frac_mud", sand = "frac_sand", gravel = "frac_gravel",
  hardness = "hardness", rock = "rock"
)
present <- intersect(names(rename_map), names(cell_vals))
if (length(present) > 0) {
  setnames(cell_vals, old = present, new = rename_map[present])
  cat("  Renamed columns: ", paste(present, "->", rename_map[present], collapse = ", "), "\n")
}

# Pre-allocate the AIS columns
ais_dt[, `:=`(
  frac_mud = NA_real_,
  frac_sand = NA_real_,
  frac_gravel = NA_real_,
  hardness = NA_real_,
  rock = NA_real_
)]

# Block update join
cols <- intersect(
  setdiff(names(cell_vals), "cell_id"),
  c("frac_mud", "frac_sand", "frac_gravel", "hardness", "rock")
)

cat("  Remapping ", length(cols), " columns out of ", format(n_ais, big.mark = " "), " pings... ")
ais_dt[cell_vals, (cols) := mget(paste0("i.", cols)), on = "cell_id"]
cat("OK\n")

# sanity: with cell_ids, at least some non-NA values must be recovered
n_mapped <- ais_dt[
  !is.na(cell_id),
  sum(!is.na(frac_mud) | !is.na(frac_sand) | !is.na(frac_gravel) |
    !is.na(rock) | !is.na(hardness))
]
if (n_mapped == 0L) stop("[P0] Remap failed: 0 values mapped (check join / i.*)")
cat("  Sanity check: ", format(n_mapped, big.mark = " "), " pings with at least 1 mapped value\n")

rm(cell_vals)
gc(verbose = FALSE)

cat("Extraction done: values remapped to all pings.\n")

## 9.  dbSEABED scale QC (CRITICAL) -------------------------------------------
cat("\n--- dbSEABED SCALE QC (before normalisation) ---\n")

qc_scale <- function(x, name) {
  # QC of the real missingness (NA) in addition to -99
  n_na <- sum(is.na(x))
  n_neg99 <- sum(x <= -90, na.rm = TRUE) # useful if NAflag was not applied
  n_gt1 <- sum(x > 1, na.rm = TRUE)
  n_gt10 <- sum(x > 10, na.rm = TRUE)
  rng <- range(x, na.rm = TRUE)
  cat(sprintf(
    "  %s: range=[%.2f, %.2f], n_NA=%d, n_nodata(-99)=%d, n_gt1=%d, n_gt10=%d\n",
    name, rng[1], rng[2], n_na, n_neg99, n_gt1, n_gt10
  ))
  invisible(NULL)
}

qc_rock <- qc_scale(ais_dt$rock, "rock")
qc_hard <- qc_scale(ais_dt$hardness, "hardness")
qc_mud <- qc_scale(ais_dt$frac_mud, "frac_mud")
qc_sand <- qc_scale(ais_dt$frac_sand, "frac_sand")
qc_grav <- qc_scale(ais_dt$frac_gravel, "frac_gravel")

## 10. dbSEABED normalisation (TYPE-SPECIFIC) ---------------------------------
cat("\nNormalisation dbSEABED (type-specific)...\n")

# Rock and texture: membership 0-100 -> probability 0-1
ais_dt[, rock_prob := norm_membership_0_100(rock)]
ais_dt[, frac_mud := norm_membership_0_100(frac_mud)]
ais_dt[, frac_sand := norm_membership_0_100(frac_sand)]
ais_dt[, frac_gravel := norm_membership_0_100(frac_gravel)]

# Hardness: signed index -1..+1 -> probability 0-1
# CRITICAL: hardness is on [-1,+1], not [0,100]!
ais_dt[, hard_soft_signed := norm_hardsoft_signed(hardness)]
ais_dt[, hard_prob := hardsoft_to_prob(hard_soft_signed)]

# Diagnostic post-normalisation
cat(
  "  rock_prob:     median=", round(median(ais_dt$rock_prob, na.rm = TRUE), 3),
  " p90=", round(quantile(ais_dt$rock_prob, 0.90, na.rm = TRUE), 3),
  " p99=", round(quantile(ais_dt$rock_prob, 0.99, na.rm = TRUE), 3), "\n"
)
cat(
  "  hard_prob:     median=", round(median(ais_dt$hard_prob, na.rm = TRUE), 3),
  " p90=", round(quantile(ais_dt$hard_prob, 0.90, na.rm = TRUE), 3),
  " p99=", round(quantile(ais_dt$hard_prob, 0.99, na.rm = TRUE), 3), "\n"
)
cat(
  "  frac_mud:      median=", round(median(ais_dt$frac_mud, na.rm = TRUE), 3),
  " p99=", round(quantile(ais_dt$frac_mud, 0.99, na.rm = TRUE), 3), "\n"
)

# Correlation check (should be high since hardness ~ 2*rock/100 - 1)
cor_rock_hard_raw <- cor(ais_dt$rock, ais_dt$hardness, use = "pairwise.complete.obs")
cat("  Corr raw rock~hardness: ", round(cor_rock_hard_raw, 3), " (expected ~1 = redundant)\n")

## 11. Lithology classification -----------------------------------------------
cat("\nLithological classification (dbSEABED-compliant v2)...\n")

litho_result <- canon_litho_v2(
  mud = ais_dt$frac_mud,
  sand = ais_dt$frac_sand,
  gravel = ais_dt$frac_gravel,
  hard = ais_dt$hard_prob,
  rock = ais_dt$rock_prob
)

ais_dt[, lithologie_canon := litho_result$litho]
ais_dt[, H_index := litho_result$H] # Continuous hard-bottom index

# Classification stats
n_rock <- sum(ais_dt$lithologie_canon == "rock")
cat(
  "  Rock confirmed (strict): ", format(n_rock, big.mark = " "), " points (",
  round(100 * n_rock / n_ais, 1), "%)\n"
)

# is_hardrock flag (WITHOUT gravel, a coarse soft sediment)
ais_dt[, is_hardrock := lithologie_canon %in% c("rock", "hard_bottom")]

## 12. CONTINUOUS pl_base (H_index blend) -------------------------------------
cat("\nComputing continuous pl_base (blend H_index x texture)...\n")

# pl_texture: two schemes, chosen by PL_SCHEME.
#
#   sala2021 (default): the rule of Sala et al. (2021,
#     Methods), applied to the texture proportions of each ping:
#       mud     > 50 % -> 0.7   ('fine')
#       gravel  > 50 % -> 0.286 ('coarse')
#       else           -> 0.04  ('sandy', "the remaining combinations")
#     Sala refers to "more than 50% of pixel area": the dbSEABED fuzzy
#     memberships do not sum to 1, so they are converted to shares (division
#     by their sum) before the threshold. Sala has no rock class: rock and
#     hard_bottom keep 0.05 and 0.10, a choice of this study. Without any valid
#     texture (unknown), 0.50 -- these pings are excluded from the carbon.
#
#   table_v1: the 11-class table (0.92 mud ... 0.05 rock). Not used: it has
#     no written source; its comments describe a penetration scale, not
#     lability, and it places sand (0.55) above gravel, contrary to Sala
#     (2021) and Black et al. (2022). Kept so that earlier runs can be
#     replayed identically (PL_SCHEME=table_v1).
PL_SCHEME <- Sys.getenv("PL_SCHEME", "sala2021")
cat("  PL_SCHEME = ", PL_SCHEME, "\n")
if (PL_SCHEME == "table_v1") {
  ais_dt[pl_base_lookup, pl_texture := i.pl_base, on = "lithologie_canon"]
} else if (PL_SCHEME == "sala2021") {
  ais_dt[, tex_tot := rowSums(cbind(frac_mud, frac_sand, frac_gravel), na.rm = TRUE)]
  ais_dt[, `:=`(
    part_mud    = fifelse(tex_tot > 0, fcoalesce(frac_mud, 0) / tex_tot, NA_real_),
    part_gravel = fifelse(tex_tot > 0, fcoalesce(frac_gravel, 0) / tex_tot, NA_real_)
  )]
  ais_dt[, pl_texture := fcase(
    lithologie_canon == "rock",        0.05,
    lithologie_canon == "hard_bottom", 0.10,
    is.na(part_mud),                   NA_real_,
    part_mud > 0.5,                    0.7,
    part_gravel > 0.5,                 0.286,
    default =                          0.04
  )]
  print(ais_dt[, .N, by = .(sala_group = fcase(
    is.na(pl_texture), "no texture (0.50)", pl_texture == 0.7, "fine 0.7",
    pl_texture == 0.286, "coarse 0.286", pl_texture == 0.04, "sandy 0.04",
    default = "rock / hard_bottom"))][order(-N)])
  ais_dt[, c("tex_tot", "part_mud", "part_gravel") := NULL]
} else {
  stop("PL_SCHEME unknown: '", PL_SCHEME, "' (expected: sala2021 | table_v1)")
}
ais_dt[is.na(pl_texture), pl_texture := 0.50]

# Continuous blend with H_index (threshold + rescaling)
# pl_base = (1 - H_blend) * pl_texture + H_blend * pl_hard
# H_blend acts only if H_index > H0 (hard-bottom activation threshold)
# Avoids penalising moderately hard bottoms

pl_hard <- 0.05
H0 <- 0.70 # activation threshold: blend only above it
ais_dt[, H_blend := fifelse(
  is.na(H_index), 0,
  pmax(0, (pmin(H_index, 1) - H0) / (1 - H0))
)]
ais_dt[, pl_base := (1 - H_blend) * pl_texture + H_blend * pl_hard]

cat(
  "  H_blend: median=", round(median(ais_dt$H_blend, na.rm = TRUE), 3),
  " p90=", round(quantile(ais_dt$H_blend, 0.90, na.rm = TRUE), 3),
  " p99=", round(quantile(ais_dt$H_blend, 0.99, na.rm = TRUE), 3), "\n"
)

# --- Cell-based QC on pl_base (not biased by AIS density) ---
qc_pl_cells <- ais_dt[!is.na(cell_id),
  .(pl_med = median(pl_base, na.rm = TRUE)),
  by = cell_id
]
cat(
  "  pl_base (by cell): median=", round(median(qc_pl_cells$pl_med, na.rm = TRUE), 3),
  " p90=", round(quantile(qc_pl_cells$pl_med, 0.90, na.rm = TRUE), 3), "\n"
)
rm(qc_pl_cells)
gc(verbose = FALSE)

# --- CODA check: texture sum ---
tex_sum <- ais_dt[, frac_mud + frac_sand + frac_gravel]
cat(
  "  texture sum: median=", round(median(tex_sum, na.rm = TRUE), 3),
  " p01=", round(quantile(tex_sum, 0.01, na.rm = TRUE), 3),
  " p99=", round(quantile(tex_sum, 0.99, na.rm = TRUE), 3), "\n"
)
rm(tex_sum)

# Cleanup
ais_dt[, c("pl_texture", "H_blend") := NULL]

## 13. Flags and compatibility ------------------------------------------------
ais_dt[, has_lithology := !is.na(frac_mud) | !is.na(frac_sand) | !is.na(frac_gravel) |
  !is.na(hard_prob) | !is.na(rock_prob)]
ais_dt[, lithologie := lithologie_canon]

# --- Keep the raw values for QA ---
ais_dt[, rock_raw := rock]
ais_dt[, hardness_raw := hardness]

# Rename for backward compatibility
ais_dt[, rock := rock_prob]
ais_dt[, hardness := hard_prob]

# Metadata attributes
attr(ais_dt, "crs") <- "+proj=longlat +datum=WGS84"
attr(ais_dt, "lithology_source") <- "HubOcean/dbSEABED"
attr(ais_dt, "lithology_version") <- timestamp
attr(ais_dt, "pipeline_version") <- "vNext+4"

## 14. Export -------------------------------------------------------------------
cat("\nSaving the results...\n")
saveRDS(ais_dt, out_file)

n_with_litho <- sum(ais_dt$has_lithology)
pct_with <- round(100 * n_with_litho / n_ais, 1)

cat("  ", format(n_with_litho, big.mark = " "), " pings with lithology (", pct_with, " %)\n")
cat("  ", format(n_ais - n_with_litho, big.mark = " "), " pings without lithology\n")

## 15. Diagnostic ---------------------------------------------------------------
cat("\n--- Lithology distribution ---\n")
litho_dist <- ais_dt[, .N, by = lithologie_canon][order(-N)]
litho_dist[, pct := round(100 * N / sum(N), 1)]
print(litho_dist)

cat("\n--- Distribution pl_base ---\n")
cat("  Min:    ", round(min(ais_dt$pl_base, na.rm = TRUE), 3), "\n")
cat("  Q1:     ", round(quantile(ais_dt$pl_base, 0.25, na.rm = TRUE), 3), "\n")
cat("  Median: ", round(median(ais_dt$pl_base, na.rm = TRUE), 3), "\n")
cat("  Mean:   ", round(mean(ais_dt$pl_base, na.rm = TRUE), 3), "\n")
cat("  Q3:     ", round(quantile(ais_dt$pl_base, 0.75, na.rm = TRUE), 3), "\n")
cat("  Max:    ", round(max(ais_dt$pl_base, na.rm = TRUE), 3), "\n")

cat("\n--- Distribution H_index ---\n")
cat("  Median: ", round(median(ais_dt$H_index, na.rm = TRUE), 3), "\n")
cat("  P90:    ", round(quantile(ais_dt$H_index, 0.90, na.rm = TRUE), 3), "\n")
cat("  P99:    ", round(quantile(ais_dt$H_index, 0.99, na.rm = TRUE), 3), "\n")

cat("\n--- Distribution is_hardrock ---\n")
hardrock_dist <- ais_dt[, .N, by = is_hardrock]
hardrock_dist[, pct := round(100 * N / sum(N), 1)]
print(hardrock_dist)

## 15.5 QC per raster cell (not biased by AIS density) ------------------------
cat("\n--- QC PER UNIQUE CELL (unbiased) ---\n")

# cell_id already computed in section 8 (extraction per unique cell)

# Stats per unique cell
qc_cells <- ais_dt[!is.na(cell_id),
  .(
    rock_med = median(rock_prob, na.rm = TRUE),
    rock_p90 = quantile(rock_prob, 0.90, na.rm = TRUE),
    hard_med = median(hard_prob, na.rm = TRUE),
    n_pings = .N
  ),
  by = cell_id
]

n_unique_cells <- nrow(qc_cells)
cat("  Unique cells:      ", format(n_unique_cells, big.mark = " "), "\n")
cat("  Median pings/cell: ", round(median(qc_cells$n_pings), 1), "\n")
cat("  rock_prob (by cell):\n")
cat(
  "    median=", round(median(qc_cells$rock_med, na.rm = TRUE), 3),
  " p90=", round(quantile(qc_cells$rock_med, 0.90, na.rm = TRUE), 3),
  " p99=", round(quantile(qc_cells$rock_med, 0.99, na.rm = TRUE), 3), "\n"
)

# % of cells with rock_med > 0.7 (truly rocky)
pct_cells_rocky <- round(100 * sum(qc_cells$rock_med > 0.7, na.rm = TRUE) / n_unique_cells, 1)
cat("  Cells rock>0.7:    ", pct_cells_rocky, "% (target: <15-25%)\n")

rm(qc_cells)
gc(verbose = FALSE)

## 16. Validation dbSEABED-compliance -------------------------------------------
cat("\n--- VALIDATION dbSEABED-COMPLIANCE ---\n")

pct_rock <- litho_dist[lithologie_canon == "rock", pct]
if (length(pct_rock) == 0) pct_rock <- 0

median_pl <- median(ais_dt$pl_base, na.rm = TRUE)

# Correlation rock_prob <-> hard_prob (positive expected, likely ~1 = redundant)
cor_rock_hard <- cor(ais_dt$rock_prob, ais_dt$hard_prob, use = "pairwise.complete.obs")

# Rock rarely associated with mud > 0.6
n_rock_total <- sum(ais_dt$lithologie_canon == "rock", na.rm = TRUE)
if (n_rock_total > 0) {
  rock_with_high_mud <- sum(ais_dt$lithologie_canon == "rock" & ais_dt$frac_mud > 0.6, na.rm = TRUE)
  pct_rock_high_mud <- round(100 * rock_with_high_mud / n_rock_total, 1)
} else {
  pct_rock_high_mud <- 0
}

cat("  Rock % (by ping):     ", pct_rock, "% (target: <15-25%)\n")
cat("  Rock % (by cell):     ", pct_cells_rocky, "% (unbiased estimate)\n")
# "target: 0.4-0.6": target without a source, written for table_v1.
# It does not apply to the sala2021 scheme; printed for the record.
cat("  Median pl_base:       ", round(median_pl, 3), " (table_v1 target, no source: 0.4-0.6)\n")
cat("  Corr rock_prob~hard:  ", round(cor_rock_hard, 3), " (expected ~1 = redundant)\n")
cat("  Rock+high mud:        ", pct_rock_high_mud, "% (target: <5%)\n")

# Verdict (use cell-based rock % for unbiased assessment)
pass_count <- 0
if (pct_cells_rocky < 25) pass_count <- pass_count + 1
if (median_pl >= 0.35 && median_pl <= 0.65) pass_count <- pass_count + 1
if (!is.na(cor_rock_hard) && cor_rock_hard > 0) pass_count <- pass_count + 1
if (pct_rock_high_mud < 5) pass_count <- pass_count + 1

if (pass_count >= 3) {
  cat("\n  [OK] Pipeline dbSEABED-compliant (", pass_count, "/4 criteria)\n")
} else {
  cat("\n  [WARN] Check calibration (", pass_count, "/4 criteria)\n")
}

## 17. Cleanup terra temp -------------------------------------------------------
tmpFiles(remove = TRUE)

cat("\nStep 4 completed successfully: ", format(Sys.time()), "\n")
cat("Output file: ", out_file, "\n")
cat("Size: ", round(file.size(out_file) / 1e6, 1), " MB\n")
