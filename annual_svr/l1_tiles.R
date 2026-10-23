#!/usr/bin/env Rscript
# l1_tiles.R — lock L1 of the annual SVR cap, on GRIT, read-only.
# For each target cell (SVR >= 1, cibles_{lab}.csv), find in the production sar_NNN.parquet files
# (output_V6_sala2021_r2/step5_{lab}) a row whose sum_dw equals that of the local input to within 1e-9 (relative).
# Any cell without a matching row => L1 FAILS, stop. Output: tuiles_{lab}.txt (tiles to rerun).
suppressPackageStartupMessages({ library(data.table); library(arrow) })
RUN <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r2")
setwd(path.expand(file.path(Sys.getenv("DFP_REPO"), "annual_svr")))
ok_all <- TRUE
for (lab in c("smooth", "flag")) {
  tgt <- fread(sprintf("cibles_%s.csv", lab))
  fs <- list.files(file.path(RUN, paste0("step5_", lab)), pattern = "^sar_[0-9]{3}\\.parquet$", full.names = TRUE)
  p <- rbindlist(lapply(fs, function(f) {
    x <- as.data.table(read_parquet(f, col_select = c("grid_id", "sum_dw")))
    x[grid_id %in% tgt$grid_id][, tile := as.integer(sub("^sar_([0-9]{3})\\.parquet$", "\\1", basename(f)))]
  }))
  j <- merge(p, tgt[, .(grid_id, dw_ref = sum_dw)], by = "grid_id")
  j[, equal := abs(sum_dw - dw_ref) <= 1e-9 * dw_ref]
  by_cell <- j[, .(n_rows = .N, n_equal = sum(equal)), by = grid_id]
  missing <- setdiff(tgt$grid_id, by_cell[n_equal >= 1]$grid_id)
  uniq_tiles <- sort(unique(j[equal == TRUE, .(t = min(tile)), by = grid_id]$t))
  cat(sprintf("== %s: %d sar files ; %d targets ; rows found %d ; cells with >= 1 equal row: %d ; with >= 2 equal rows: %d ; without: %d\n",
              lab, length(fs), nrow(tgt), nrow(j), sum(by_cell$n_equal >= 1), sum(by_cell$n_equal >= 2), length(missing)))
  if (length(missing)) { ok_all <- FALSE; cat("   cells without yr equal row:", head(missing, 20), "\n") }
  cat("   tiles to rerun (", length(uniq_tiles), "):", paste(uniq_tiles, collapse = ","), "\n")
  writeLines(paste(uniq_tiles, collapse = ","), sprintf("tuiles_%s.txt", lab))
}
cat(if (ok_all) "L1 PASS\n" else "L1 FAILS\n")
