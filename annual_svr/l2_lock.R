#!/usr/bin/env Rscript
# l2_lock.R — lock L2 of the annual SVR cap, on GRIT, read-only.
# For each rerun tile: Σ_years sum_dw (worker_year.R) = production sum_dw, cell by cell, to within 1e-9
# (relative), and same set of grid_id. Failure => stop, nothing is exported.
# If L2 passes: annee_{lab}.csv = (grid_id, tile, annee, sum_dw, sum_dw_pd) of the target cells, read in the tile
# kept by L1 (the smallest tile whose production row has the sum_dw of the local input).
suppressPackageStartupMessages({ library(data.table); library(arrow) })
RUN <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r2")
NEW <- file.path(Sys.getenv("DFP_SCRATCH"), "svr_annuel_20261005")
setwd(path.expand(file.path(Sys.getenv("DFP_REPO"), "annual_svr")))
ok_all <- TRUE; outputs <- list()
for (lab in c("smooth", "flag")) {
  uniq_tiles <- as.integer(strsplit(readLines(sprintf("tuiles_%s.txt", lab)), ",")[[1]])
  tgt <- fread(sprintf("cibles_%s.csv", lab))
  gap_max <- 0; n_cell <- 0; n_na <- 0; pb <- character(); rows <- list()
  for (t in uniq_tiles) {
    fp <- file.path(RUN, paste0("step5_", lab), sprintf("sar_%03d.parquet", t))
    fn <- file.path(NEW, paste0("step5_", lab), sprintf("sar_annee_%03d.parquet", t))
    if (!file.exists(fn)) { pb <- c(pb, sprintf("tile %d: output missing", t)); next }
    p <- as.data.table(read_parquet(fp, col_select = c("grid_id", "sum_dw")))
    n <- as.data.table(read_parquet(fn))
    n_na <- n_na + sum(is.na(n$annee))
    a <- n[, .(dw = sum(sum_dw)), by = grid_id]
    if (!setequal(a$grid_id, p$grid_id) || anyDuplicated(p$grid_id))
      pb <- c(pb, sprintf("tile %d: different grid_id (%d against %d)", t, nrow(a), nrow(p)))
    j <- merge(p, a, by = "grid_id")
    e <- max(abs(j$dw - j$sum_dw) / pmax(j$sum_dw, 1e-300)); gap_max <- max(gap_max, e); n_cell <- n_cell + nrow(j)
    if (e > 1e-9) pb <- c(pb, sprintf("tile %d: max relative gap %.3g", t, e))
    # tile kept by L1, cell by cell
    m <- merge(p[grid_id %in% tgt$grid_id], tgt[, .(grid_id, dw_ref = sum_dw)], by = "grid_id")
    m <- m[abs(sum_dw - dw_ref) <= 1e-9 * dw_ref, .(grid_id, tile = t)]
    rows[[length(rows) + 1]] <- merge(m, n, by = "grid_id")
  }
  L <- rbindlist(rows)
  L <- L[L[, .I[tile == min(tile)], by = grid_id]$V1]          # the smallest equal tile, as in L1
  cat(sprintf("== %s: %d tiles, %d cells compared, max relative gap %.3g ; NA years: %d ; targets exported %d / %d\n",
              lab, length(uniq_tiles), n_cell, gap_max, n_na, uniqueN(L$grid_id), nrow(tgt)))
  if (length(pb)) { ok_all <- FALSE; cat("   ", pb, sep = "\n   ") }
  if (uniqueN(L$grid_id) != nrow(tgt)) { ok_all <- FALSE; cat("   targets missing at export\n") }
  outputs[[lab]] <- L
}
if (!ok_all) { cat("L2 FAILS: nothing is exported\n"); quit(status = 1) }
for (lab in names(outputs)) fwrite(outputs[[lab]][order(grid_id, annee)], sprintf("annee_%s.csv", lab))
cat("L2 PASS\n")
