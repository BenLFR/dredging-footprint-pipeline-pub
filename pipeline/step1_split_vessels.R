#!/usr/bin/env Rscript
# =====================================================================
# STEP 1 ─ SPLITTING AIS BY VESSEL  (pipeline V6)
# Uses the period and fleet selected by step0_core_window_enhanced.R (core_window.yaml)
# =====================================================================

cat("\n🚀  STEP-1  |  SPLITTING AIS BY VESSEL  |  start:",
    format(Sys.time()), "\n\n")

# CRITICAL R CONFIGURATION - BEFORE LOADING PACKAGES
.libPaths("~/R/library")
cat("✅ R configured with library:", .libPaths()[1], "\n")

suppressPackageStartupMessages({
  library(data.table)   # ultra-fast + low RAM
  library(yaml)         # reads the step 0 configuration
  
  # Use qs if available, otherwise RDS
  use_qs <- FALSE
  tryCatch({
    if (!requireNamespace("qs", quietly = TRUE))
      stop("qs not available")
    library(qs)           
    use_qs <- TRUE
    cat("✅ Using the QS format\n")
  }, error = function(e) {
    cat("⚠️ QS not available, falling back to RDS\n")
    use_qs <<- FALSE
  })
  
  library(lubridate)    # robust date parsers
  library(tools)        # for file_ext
})

# ---------------------------------------------------------------------
# 1. ── PARAMETERS (adjustable via environment variables) ----------
# ---------------------------------------------------------------------
job_id        <- Sys.getenv("SLURM_JOB_ID",        unset = format(Sys.time(), "%Y%m%d%H%M%S"))
input_pattern <- Sys.getenv("AIS_INPUT_PATTERN",   unset = file.path(Sys.getenv("DFP_DATA"), "AIS_data/*.csv"))
output_dir    <- Sys.getenv("AIS_OUTPUT_DIR",      unset = file.path(Sys.getenv("DFP_DATA"), paste0("ais_split_", job_id)))
core_config_path <- Sys.getenv("CORE_CONFIG_PATH", unset = file.path(Sys.getenv("DFP_DATA"), "output_V6/core_window.yaml"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("📁  INPUT  : ", input_pattern, "\n")
cat("📁  OUTPUT: ", output_dir,     "\n")
cat("📋  CONFIG: ", core_config_path, "\n\n")

# ---------------------------------------------------------------------
# 1.b ── READING THE STEP 0 CONFIGURATION  -----------------------------
# ---------------------------------------------------------------------
if (file.exists(core_config_path)) {
  cat("🔍  Reading the step 0 configuration...\n")
  core_config <- read_yaml(core_config_path)
  
  start_year <- core_config$start_year
  end_year   <- core_config$end_year
  core_ships <- core_config$core_ships
  
  cat("📅  Time window:", start_year, "-", end_year, "\n")
  cat("🚢  Selected vessels:", length(core_ships), "\n")
  cat("   •", paste(core_ships, collapse = ", "), "\n\n")
  
  if (is.null(start_year) || is.null(end_year) || is.null(core_ships)) {
    stop("❌  Incomplete step 0 configuration in", core_config_path)
  }
  
  use_core_filter <- TRUE
} else {
  cat("⚠️  Step 0 configuration not found, processing the full data\n")
  use_core_filter <- FALSE
  start_year <- NULL
  end_year <- NULL
  core_ships <- NULL
}

# ---------------------------------------------------------------------
# 2. ── MMSI → VESSEL NAME MAPPING  -----------------------------------
# ---------------------------------------------------------------------
# Note: Vasco Da Gama had two MMSIs: 253193000 (Luxembourg, 2013-2019) and 205744000 (Belgium, 2018-2024)
# We merge into the most recent MMSI: 205744000
# Note: Goryo 6 Ho (312062000) excluded - corrupted data
mmsi_map <- data.table(
  ssvid = as.character(c(209469000, 210138000, 245508000, 246351000,
            253193000, 205744000, 253373000, 253403000, 253422000,
            253688000, 533180137)),
  Navire = c("Fairway", "Queen Of The Netherlands", "Ham 318",
             "Vox Maxima", "Vasco Da Gama", "Vasco Da Gama", "Cristobal Colon",
             "Leiv Eiriksson", "Charles Darwin", "Congo River",
             "Inai Kenanga"),
  # Reference MMSI (the most recent for Vasco Da Gama)
  ssvid_ref = as.character(c(209469000, 210138000, 245508000, 246351000,
                205744000, 205744000, 253373000, 253403000, 253422000,
                253688000, 533180137))
)
setkey(mmsi_map, ssvid)

stopifnot(is.character(mmsi_map$ssvid))
stopifnot(is.character(mmsi_map$ssvid_ref))

# ---------------------------------------------------------------------
# 2.b ── ALIAS → CANONICAL NAME MAPPING  ------------------------------
# ---------------------------------------------------------------------
# Handling of name differences between CSV and step 0 configuration
alias_map <- data.table(
  alias = c("Ham 318 Sleephopperzuiger", "HAM318",
            "Queen of the netherlands", "Queen Of The Netherlands",
            "Vasco de Gama", "Vasco Da Gama",
            "LEIV EIRIKSSONN", "Leiv Eiriksson",
            "INAI KENANGA", "Inai Kenanga",
            "Goryo 6 HO", "Goryo 6 Ho",
            "Fair Way", "Fairway",
            "VOX maxima", "Vox Maxima"),
  canonical = c("Ham 318", "Ham 318",
                "Queen Of The Netherlands", "Queen Of The Netherlands", 
                "Vasco Da Gama", "Vasco Da Gama",
                "Leiv Eiriksson", "Leiv Eiriksson",
                "Inai Kenanga", "Inai Kenanga", 
                "Goryo 6 Ho", "Goryo 6 Ho",
                "Fairway", "Fairway",
                "Vox Maxima", "Vox Maxima")
)
setkey(alias_map, alias)

# ---------------------------------------------------------------------
# 3. ── READING THE DATA (CSV or RDS)  ---------------------------------
# ---------------------------------------------------------------------
# Priority to a single input file (RDS/CSV) via the environment variable
# AIS_INPUT_FILE. Otherwise, fallback to the CSV pattern.

input_file <- Sys.getenv("AIS_INPUT_FILE", unset = "")

t_read <- system.time({
  if (input_file != "" && file.exists(input_file)) {
    cat("🔍  Reading the input file via AIS_INPUT_FILE:", input_file, "\n")
    
    file_ext <- tolower(tools::file_ext(input_file))
    
    if (file_ext == "rds") {
      ais_dt <- as.data.table(readRDS(input_file))
    } else if (file_ext == "csv") {
      ais_dt <- fread(input_file, showProgress = FALSE)
    } else {
      stop("❌ Unsupported file format for AIS_INPUT_FILE: ", file_ext)
    }
    
  } else {
    cat("🔍  AIS_INPUT_FILE not given or not found. Using the CSV pattern:", input_pattern, "\n")
    csv_files <- Sys.glob(input_pattern)
    
    if (length(csv_files) == 0) stop("❌ No CSV file found for the pattern: ", input_pattern)
    if (length(csv_files) > 1) cat("⚠️  Several CSV files detected. Processing only the first one:", basename(csv_files[1]), "\n")
    
    cat("   • Reading:", basename(csv_files[1]), "\n")
    ais_dt <- fread(csv_files[1], showProgress = FALSE)
  }
})
cat(sprintf("✅  Reading finished: %s rows  |  %.1f s\n",
            format(nrow(ais_dt), big.mark = " "), t_read[3]))

# ---------------------------------------------------------------------
# 4. ── COLUMN NORMALISATION  ----------------------------------------
# ---------------------------------------------------------------------
# a) lower-case names
setnames(ais_dt, tolower(names(ais_dt)))

# b) minimal renaming dictionary
rename_map <- c(lon="Lon", lat="Lat", course="Course", timestamp="Timestamp",
                speed="Speed", speed_knots="Speed", seg_id="Seg_id",
                trip_id="Seg_id", navire="Navire")
common <- intersect(names(rename_map), names(ais_dt))
setnames(ais_dt, common, rename_map[common])

# c) essential type coercions
num_cols <- c("Lon","Lat","Speed")
for (cl in intersect(num_cols, names(ais_dt))) set(ais_dt, j = cl, value = as.numeric(ais_dt[[cl]]))
if ("Timestamp" %chin% names(ais_dt))
  ais_dt[, Timestamp := as.POSIXct(Timestamp, tz = "UTC")]

# ---------------------------------------------------------------------
# 5. ── CREATION / VALIDATION OF THE VESSEL COLUMN  --------------------
# ---------------------------------------------------------------------
if (!"Navire" %in% names(ais_dt)) ais_dt[, Navire := NA_character_]

# If ssvids are present, we fill in the missing names
if ("ssvid" %chin% names(ais_dt)) {
  ais_dt[, ssvid := as.character(ssvid)]
  ais_dt <- merge(ais_dt, mmsi_map, by = "ssvid", all.x = TRUE, suffixes = c("", ".map"))
  
  # Merging of multiple MMSIs into the reference MMSI (Vasco Da Gama case)
  ais_dt[!is.na(ssvid_ref), ssvid := ssvid_ref]
  ais_dt[, ssvid_ref := NULL]
  
  ais_dt[is.na(Navire), Navire := Navire.map]
  ais_dt[, Navire.map := NULL]
}

# ultimate fallback: generic name
ais_dt[is.na(Navire) | Navire == "", Navire := paste0("unknown_", .I)]

# ---------------------------------------------------------------------
# 5.a ── NAME NORMALISATION VIA ALIASES  -----------------------------
# ---------------------------------------------------------------------
# Application of the alias mapping to harmonise the names with step 0
if (use_core_filter) {
  cat("🔍  Normalising vessel names...\n")
  
  tryCatch({
    names_before <- unique(ais_dt$Navire)
    cat("   • Names before normalisation:", paste(names_before, collapse = ", "), "\n")
    
    ais_dt <- merge(ais_dt, alias_map, by.x = "Navire", by.y = "alias", all.x = TRUE)
    
    ais_dt[!is.na(canonical), Navire := canonical]
    
    if ("canonical" %in% names(ais_dt)) {
      ais_dt[, canonical := NULL]
    }
    
    names_after <- unique(ais_dt$Navire)
    cat("   • Names after normalisation:", paste(names_after, collapse = ", "), "\n")
    
    if (length(names_before) != length(names_after)) {
      cat("✅  Normalisation:", length(names_before), "→", length(names_after), "unique names\n")
    } else {
      cat("ℹ️  No name change needed\n")
    }
  }, error = function(e) {
    cat("❌  Error during normalisation:", e$message, "\n")
    cat("   • Available columns:", paste(names(ais_dt), collapse = ", "), "\n")
    cat("   • alias_map size:", nrow(alias_map), "rows\n")
    stop(e)
  })
}

# ---------------------------------------------------------------------
# 5.b ── FILTERING ACCORDING TO STEP 0 CONFIGURATION  -------------------
# ---------------------------------------------------------------------
if (use_core_filter) {
  cat("🔍  Applying the step 0 filters...\n")
  
  ais_dt[, Annee := year(Timestamp)]
  
  n_before <- nrow(ais_dt)
  ais_dt <- ais_dt[Navire %in% core_ships & Annee >= start_year & Annee <= end_year]
  n_after <- nrow(ais_dt)
  
  cat(sprintf("📊  Filtering: %s → %s rows (%.1f%% kept)\n",
              format(n_before, big.mark = " "),
              format(n_after, big.mark = " "),
              round(100 * n_after / n_before, 1)))
  
  present_vessels <- unique(ais_dt$Navire)
  missing_vessels <- setdiff(core_ships, present_vessels)
  
  if (length(missing_vessels) > 0) {
    cat("⚠️  Vessels missing in the window:", paste(missing_vessels, collapse = ", "), "\n")
  }
  
  cat("✅  Filtering finished\n\n")
} else {
  cat("ℹ️  No filter applied (full processing)\n\n")
}

# ---------------------------------------------------------------------
# 6. ── SPLITTING BY VESSEL  ------------------------------------------
# ---------------------------------------------------------------------
# stable order: size desc to show progress
vessel_list <- ais_dt[, .N, by = Navire][order(-N)]
cat("🚢  Vessels detected:", vessel_list[,.N], "\n\n")

meta <- vessel_list[, `:=`(file_path = character(.N), split_time = Sys.time())]

pb <- txtProgressBar(min = 0, max = nrow(vessel_list), style = 3)
i <- 0

for (vname in vessel_list$Navire) {
  i <- i + 1
  dt_vessel <- ais_dt[Navire == vname]

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", vname)
  
  if (use_qs) {
    # QS format (preferred)
    out_file  <- file.path(output_dir,
                           sprintf("navire_%02d_%s.qs", i, safe_name))
    qs::qsave(dt_vessel, out_file, preset = "custom",
              algorithm = "zstd",         # fast & good ratio
              compress_level = 6)     # ~ equiv. fst compress=85
  } else {
    # RDS format (fallback)
  out_file  <- file.path(output_dir,
                           sprintf("navire_%02d_%s.rds", i, safe_name))
    saveRDS(dt_vessel, out_file, compress = "xz")
  }

  meta[Navire == vname, `:=`(file_path = out_file,
                           n_observations = nrow(dt_vessel))]
  rm(dt_vessel); gc(verbose = FALSE)
  setTxtProgressBar(pb, i)
}
close(pb)

# ---------------------------------------------------------------------
# 7. ── METADATA EXPORT & INTEGRITY CHECK  ----------------------------
# ---------------------------------------------------------------------
meta_file <- file.path(output_dir, "navires_metadata.csv")
fwrite(meta, meta_file)

if (sum(meta$n_observations) != nrow(ais_dt))
  stop("❌  Integrity KO: different number of rows after split!")
cat("\n✅  Integrity OK: ", sum(meta$n_observations), " rows checked.\n")

# ---------------------------------------------------------------------
# 8. ── SUMMARY --------------------------------------------------------
# ---------------------------------------------------------------------
cat("\n📊  SUMMARY -------------------------------------------------------\n")
print(meta[, .(Navire, n_observations, file_path)])

total_size <- sum(file.info(meta$file_path)$size) / 1024^2
format_used <- if(use_qs) "qs" else "rds"
cat(sprintf("\n💾  %d .%s files (%.1f MB in total) written to: %s\n",
            nrow(meta), format_used, total_size, output_dir))

cat("\n🏁  STEP-1 completed successfully:", format(Sys.time()), "\n")
