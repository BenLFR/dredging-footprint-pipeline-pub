#!/usr/bin/env Rscript
# ==============================================================================
# STEP-0 ENHANCED: Core Period Selection via Coverage Matrix (Eriksen et al., 2018 style)
# ==============================================================================

cat("\n📊  STEP-0 ENHANCED  |  Core Period Selection via Coverage Matrix |  start:", format(Sys.time()), "\n\n")

# ---- PACKAGES AND CONFIG ----
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(matrixStats)
  library(digest)
  library(yaml)
})

set.seed(42)

# ---- DATA LOADING ----
input_pattern <- Sys.getenv("AIS_INPUT_PATTERN", file.path(Sys.getenv("DFP_DATA"), "AIS_data/benjamin3_clean.csv"))
output_dir   <- Sys.getenv("AIS_OUTPUT_DIR", file.path(Sys.getenv("DFP_DATA"), "output_V6"))
cat("📁  INPUT  :", input_pattern, "\n")
cat("📁  OUTPUT:", output_dir, "\n\n")

# Detection of the input file (CSV **or** RDS)
input_file <- Sys.getenv("AIS_INPUT_FILE", "")
if (nzchar(input_file) && file.exists(input_file)) {
  cat("🔍  File given by AIS_INPUT_FILE:", basename(input_file), "\n\n")
} else {
csv_files <- Sys.glob(input_pattern)
  rds_pattern <- sub("\\*\\.csv$", "*.rds", input_pattern)
  rds_files <- Sys.glob(rds_pattern)
  if (length(csv_files)) {
    pref_csv <- csv_files[grep("benjamin2\\.csv$", csv_files, ignore.case = TRUE)]
    input_file <- if (length(pref_csv)) pref_csv[1] else csv_files[1]
  } else if (length(rds_files)) {
    input_file <- rds_files[1]
  } else {
    stop("❌  No input file found (CSV or RDS).")
  }
  cat("🔍  File detected:", basename(input_file), "\n\n")
}

# ---- FILE LOADING -----------------------------------------------------------
if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
  ais_dt <- as.data.table(readRDS(input_file))
} else {
  ais_dt <- fread(input_file, showProgress = FALSE)
}

setnames(ais_dt, tolower(names(ais_dt)))
cat("✅", format(nrow(ais_dt), big.mark = " "), "rows loaded.\n\n")

# ---- MMSI → VESSEL LOOKUP TABLE ----
# Note: Vasco Da Gama had two MMSIs: 253193000 (Luxembourg, 2013-2019) and 205744000 (Belgium, 2018-2024)
# We merge into the most recent MMSI: 205744000
# Note: Goryo 6 Ho (312062000) excluded - corrupted data
mmsi_map <- data.table(
  ssvid = c(209469000, 210138000, 245508000, 246351000,
            253193000, 205744000, 253373000, 253403000, 253422000,
            253688000, 533180137),
  Navire = c("Fairway", "Queen Of The Netherlands", "Ham 318",
             "Vox Maxima", "Vasco Da Gama", "Vasco Da Gama", "Cristobal Colon",
             "Leiv Eiriksson", "Charles Darwin", "Congo River",
             "Inai Kenanga"),
  # Reference MMSI (the most recent for Vasco Da Gama)
  ssvid_ref = c(209469000, 210138000, 245508000, 246351000,
                205744000, 205744000, 253373000, 253403000, 253422000,
                253688000, 533180137)
)
mmsi_map[, ssvid := as.character(ssvid)]
mmsi_map[, ssvid_ref := as.character(ssvid_ref)]
setkey(mmsi_map, ssvid)

if ("ssvid" %in% names(ais_dt)) {
  ais_dt[, ssvid := as.character(ssvid)]
  ais_dt <- merge(ais_dt, mmsi_map, by = "ssvid", all.x = TRUE)
  
  # Merging of multiple MMSIs into the reference MMSI (Vasco Da Gama case)
  ais_dt[!is.na(ssvid_ref), ssvid := ssvid_ref]
  ais_dt[, ssvid_ref := NULL]
  
  ais_dt[is.na(Navire), Navire := paste0("unknown_", ssvid)]
} else if ("navire" %in% names(ais_dt)) {
  setnames(ais_dt, "navire", "Navire")
} else {
  ais_dt[, Navire := "benjamin2"]
}

# ---- TIME PREPARATION ----
ais_dt[, Timestamp := as.POSIXct(get("timestamp"), tz = "UTC")]
ais_dt[, Date := as.Date(Timestamp)]
ais_dt[, Annee := year(Timestamp)]

# ---- COVERAGE MATRIX ----
# 1. Daily coverage for each vessel-year
cover_dt <- ais_dt[, .(active_days = uniqueN(Date)), by = .(Navire, Annee)]
cover_dt[, total_days := ifelse(leap_year(Annee), 366, 365)]
cover_dt[, coverage := active_days / total_days]

# 2. Vessel × Year matrix
vessel_list <- sort(unique(cover_dt$Navire))
yrs  <- sort(unique(cover_dt$Annee))
cov_mat <- matrix(0, nrow = length(vessel_list), ncol = length(yrs),
                  dimnames = list(vessel_list, yrs))
for (i in seq_along(vessel_list)) {
  for (j in seq_along(yrs)) {
    v <- cover_dt[Navire == vessel_list[i] & Annee == yrs[j], coverage]
    if (length(v)) cov_mat[i, j] <- v
  }
}

cat("✅ Coverage matrix:", nrow(cov_mat), "vessels ×", ncol(cov_mat), "years\n")

# ---- EXPORT OF THE VESSEL × YEAR MATRIX ----
# Save the full matrix for visualisation
cov_dt <- as.data.table(cov_mat, keep.rownames = "ship")
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
fwrite(cov_dt, file.path(output_dir, "coverage_matrix.csv"), row.names = FALSE)
cat("📊 Coverage matrix exported: coverage_matrix.csv\n")

# ---- FULL ANALYSIS OF ALL WINDOWS ----
min_L <- 5
all_windows <- data.table()
years <- as.numeric(colnames(cov_mat))
n <- length(years)

cat("🔍 Analysis of all candidate windows (min", min_L, "years)...\n")
pb <- txtProgressBar(min = 0, max = (n - min_L + 1) * (n - min_L) / 2, style = 3)
counter <- 0

for (i in seq_len(n - min_L + 1)) {
  for (j in seq(i+min_L-1, n)) {
    subC <- cov_mat[, i:j, drop=FALSE]
    Cyears <- colMedians(subC, na.rm = TRUE)
    L   <- j-i+1
    Cbar<- median(Cyears, na.rm=TRUE)
    Cmin<- min(Cyears, na.rm=TRUE)
    CV  <- ifelse(Cbar==0, 1e9, sd(Cyears, na.rm=TRUE)/Cbar)
    S   <- L * Cbar * (Cmin^2) / (1 + CV)
    
    active_vessels <- sum(rowSums(subC > 0) > 0)
    total_cov <- mean(subC, na.rm=TRUE)
    
    all_windows <- rbind(all_windows, data.table(
      start_year = years[i],
      end_year = years[j],
      length = L,
      median_coverage = Cbar,
      min_coverage = Cmin,
      CV_coverage = CV,
      score = S,
      active_ships = active_vessels,
      total_coverage = total_cov,
      years_list = list(years[i]:years[j])
    ))
    
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
}
close(pb)

setorder(all_windows, -score)

cat(sprintf("✅ %d candidate windows analysed\n", nrow(all_windows)))

# ---- OPTIMAL WINDOW ----
best <- all_windows[1]
best_window_yrs <- best$years_list[[1]]

cat(sprintf("🎯 Optimal window: %d-%d (%d years)\n", best$start_year, best$end_year, best$length))
cat(sprintf("   • Median coverage: %.1f %%\n", 100*best$median_coverage))
cat(sprintf("   • Minimum coverage    : %.1f %%\n", 100*best$min_coverage))
cat(sprintf("   • Annual CV           : %.3f\n", best$CV_coverage))
cat(sprintf("   • Score               : %.3f\n", best$score))
cat(sprintf("   • Active vessels   : %d\n", best$active_ships))

# ---- ANALYSIS OF THE ALTERNATIVES ----
cat("\n🏆 TOP 10 CANDIDATE WINDOWS:\n")
top_10 <- all_windows[1:10]
for (i in 1:nrow(top_10)) {
  w <- top_10[i]
  cat(sprintf("   %2d. %d-%d (%d years) | Score: %.3f | Median: %.1f%% | Vessels: %d\n", 
              i, w$start_year, w$end_year, w$length, w$score, 
              100*w$median_coverage, w$active_ships))
}

if (nrow(all_windows) >= 2) {
  score_diff <- best$score - all_windows$score[2]
  score_diff_pct <- 100 * score_diff / all_windows$score[2]
  cat(sprintf("\n📊 Advantage of the best window: %.3f (%.1f%%)\n", score_diff, score_diff_pct))
  
  if (score_diff_pct < 5) {
    cat("⚠️  Warning: the difference with the 2nd best window is small (< 5%)\n")
  }
}

# ---- REMOVAL OF GHOST VESSELS ----
yrs_idx <- which(years %in% best_window_yrs)
present_n <- rowSums(cov_mat[, yrs_idx, drop = FALSE] > 0)
core_ships <- names(present_n[present_n > 0])

if (length(core_ships) == 0) {
  stop("❌ The optimal window finally contains no active vessel – to be checked.")
}

cat(sprintf("🔍  Active vessels in the window: %d (out of %d total)\n", 
            length(core_ships), nrow(cov_mat)))

if (length(core_ships) < nrow(cov_mat)) {
  inactive_ships <- setdiff(rownames(cov_mat), core_ships)
  cat("⚠️  Inactive vessels excluded:", paste(inactive_ships, collapse = ", "), "\n")
}

cov_mat <- cov_mat[core_ships, , drop = FALSE]
best_ships <- core_ships

# ---- FULL EXPORT ----
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# 1. All candidate windows
fwrite(all_windows[, -"years_list"], file.path(output_dir, "all_candidate_windows.csv"))

# 2. Top 20 windows for quick analysis
fwrite(all_windows[1:20, -"years_list"], file.path(output_dir, "top_20_windows.csv"))

# 3. CSV with the metrics of the optimal window
core_window <- data.table(
  start_year = best$start_year,
  end_year   = best$end_year,
  length     = best$length,
  median_coverage = best$median_coverage,
  min_coverage = best$min_coverage,
  CV_coverage = best$CV_coverage,
  score = best$score,
  active_ships = best$active_ships,
  total_coverage = best$total_coverage
)
fwrite(core_window, file.path(output_dir, "best_core_window.csv"))

# 4. Full YAML with window AND vessel list
core_config <- list(
  start_year = best$start_year,
  end_year = best$end_year,
  length_years = best$length,
  
  median_coverage = best$median_coverage,
  min_coverage = best$min_coverage,
  cv_coverage = best$CV_coverage,
  score = best$score,
  active_ships = best$active_ships,
  total_coverage = best$total_coverage,
  
  core_ships = best_ships,
  
  window_years = best_window_yrs,
  
  total_candidates = nrow(all_windows),
  score_advantage = if(nrow(all_windows) >= 2) best$score - all_windows$score[2] else NA,
  score_advantage_pct = if(nrow(all_windows) >= 2) 100 * (best$score - all_windows$score[2]) / all_windows$score[2] else NA,
  
  execution_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_file = basename(input_file),
  total_ships = length(best_ships),
  total_observations = nrow(ais_dt)
)

write_yaml(core_config, file.path(output_dir, "core_window.yaml"))

# 5. CSV with the vessel list
ships_dt <- data.table(
  ship_id = seq_along(best_ships),
  ship_name = best_ships
)
fwrite(ships_dt, file.path(output_dir, "core_ships.csv"))

# 6. Detailed markdown report
report_md <- sprintf(
  "# Core Window Enhanced Report (Coverage Matrix)

## 🎯 Optimal Window
- **Period:** %d–%d (%d years)
- **Median coverage:** %.1f %%
- **Minimum coverage:** %.1f %%
- **Inter-annual CV:** %.3f
- **Optimisation score:** %.3f
- **Active vessels:** %d

## 📊 Analysis of the Alternatives
- **Total windows analysed:** %d
- **Advantage over the 2nd best:** %.3f (%.1f%%)

### Top 5 Alternative Windows
%s

## 🔍 Technical Details
- **Minimum window:** 5 years
- **Optimisation criterion:** L × Cbar × (Cmin²) / (1 + CV)
- **Excluded vessels:** %s

## 📁 Generated Files
- `all_candidate_windows.csv`: all candidate windows
- `top_20_windows.csv`: top 20 windows for quick analysis
- `best_core_window.csv`: metrics of the optimal window
- `core_window.yaml`: full configuration
- `core_ships.csv`: list of active vessels

---
*Generated on %s*
", 
  best$start_year, best$end_year, best$length, 100*best$median_coverage, 
  100*best$min_coverage, best$CV_coverage, best$score, length(best_ships),
  nrow(all_windows),
  if(nrow(all_windows) >= 2) best$score - all_windows$score[2] else 0,
  if(nrow(all_windows) >= 2) 100 * (best$score - all_windows$score[2]) / all_windows$score[2] else 0,
  paste(sapply(1:min(5, nrow(all_windows)), function(i) {
    w <- all_windows[i]
    sprintf("  %d. %d-%d (%d years) | Score: %.3f | Median: %.1f%%", 
            i, w$start_year, w$end_year, w$length, w$score, 100*w$median_coverage)
  }), collapse = "\n"),
  if(length(best_ships) < nrow(cov_mat)) paste(setdiff(rownames(cov_mat), best_ships), collapse = ", ") else "None",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(report_md, file.path(output_dir, "core_window_enhanced_report.md"))

# 7. Statistics per year for the optimal window
yearly_stats <- data.table(
  year = best_window_yrs,
  median_coverage = sapply(best_window_yrs, function(y) {
    median(cov_mat[, as.character(y)], na.rm = TRUE)
  }),
  min_coverage = sapply(best_window_yrs, function(y) {
    min(cov_mat[, as.character(y)], na.rm = TRUE)
  }),
  active_ships = sapply(best_window_yrs, function(y) {
    sum(cov_mat[, as.character(y)] > 0, na.rm = TRUE)
  })
)
fwrite(yearly_stats, file.path(output_dir, "yearly_coverage_stats.csv"))

# ---- FINAL SUMMARY ----
cat("\n📄  Markdown summary: core_window_enhanced_report.md\n")
cat("📋  Configuration YAML: core_window.yaml\n")
cat("🚢  Vessel list: core_ships.csv\n")
cat("📊  All windows: all_candidate_windows.csv\n")
cat("🏆  Top 20 windows: top_20_windows.csv\n")
cat("📈  Annual stats: yearly_coverage_stats.csv\n")
cat("💾  Files in", output_dir, "\n")
cat("⏱️  End:", format(Sys.time()), "\n")

cat("\n✅ STEP-0 ENHANCED completed successfully:", format(Sys.time()), "\n")
cat("📄  Selection file:", file.path(output_dir, "core_window_enhanced_report.md"), "\n") 