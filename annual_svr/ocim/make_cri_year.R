#!/usr/bin/env Rscript
# make_cri_year.R -- input of the OCIM run dedicated to the annual SVR cap (rule written beforehand).
# cri_final_<lab>_an_20261005.parquet = parquet of RUN/cri_pr where C_ri, C_ri_lower, C_ri_upper, C_ri_conservative are
# multiplied by s_i (s_cibles_<lab>.csv; s_i = 1 outside cells with SVR >= 1). Declared choice: same s_i
# for the four columns (s_i depends only on SVR; Eq. 8 linear in SVR). Modelling choice, no source.
# Locks (the file is written only if they pass; otherwise exit with an error, and nothing downstream starts):
#   A1  same rows, same grid_id, all other columns identical -- checked on the RE-READ file
#   A2  original Σ C_ri = 4038.1542 / 5654.3864 (V0) and annual Σ C_ri = 3432.5263 / 4685.5941
#       (s_targets_20261005.R), to within 0.006
#   A3  new / old = s_i to within 1e-12 where the old is > 0; NA at the same places; the rest unchanged
# Usage: Rscript make_cri_year.R <RUN> <BASE> <CODE>
suppressPackageStartupMessages({ library(data.table); library(arrow) })
a <- commandArgs(TRUE); RUN <- a[1]; BASE <- a[2]; CODE <- a[3]
COLS <- c("C_ri", "C_ri_lower", "C_ri_upper", "C_ri_conservative")
V0   <- c(smooth = 4038.1542, flag = 5654.3864)   # reference totals smooth / flag (V0)
EXP  <- c(smooth = 3432.5263, flag = 4685.5941)   # s_targets_20261005.R
TOL_SUM <- 0.006                                # that of V0
ok_all <- TRUE
for (lab in c("smooth", "flag")) {
  cat(sprintf("\n==================== %s ====================\n", lab))
  f0 <- Sys.glob(file.path(RUN, "cri_pr", sprintf("cri_final_%s_*.parquet", lab)))
  if (length(f0) != 1) { cat("ABORT:", length(f0), "original parquet\n"); ok_all <- FALSE; next }
  cat("origin:", f0, " md5", tools::md5sum(f0), "\n")
  d0 <- as.data.table(read_parquet(f0))
  s  <- fread(file.path(CODE, sprintf("s_cibles_%s.csv", lab)))
  if (anyDuplicated(s$grid_id) || !all(s$grid_id %in% d0$grid_id) || anyDuplicated(d0$grid_id)) {
    cat("ABORT: target grid_id missing or duplicated\n"); ok_all <- FALSE; next }
  si <- rep(1, nrow(d0)); si[match(s$grid_id, d0$grid_id)] <- s$s
  d1 <- copy(d0)
  for (k in COLS) set(d1, j = k, value = d0[[k]] * si)
  dir.create(file.path(BASE, paste0("cri_an_", lab)), recursive = TRUE, showWarnings = FALSE)
  out <- file.path(BASE, paste0("cri_an_", lab), sprintf("cri_final_%s_an_20261005.parquet", lab))
  tmp <- paste0(out, ".tmp")
  write_parquet(d1, tmp)
  r <- as.data.table(read_parquet(tmp))                      # locks on the re-read file, not on the object
  others <- setdiff(names(d0), COLS)
  a1 <- identical(names(r), names(d0)) && nrow(r) == nrow(d0) && identical(r$grid_id, d0$grid_id) &&
        all(vapply(others, function(k) identical(r[[k]], d0[[k]]), logical(1)))
  s0 <- sum(d0$C_ri, na.rm = TRUE); s1 <- sum(r$C_ri, na.rm = TRUE)
  a2 <- abs(s0 - V0[[lab]]) <= TOL_SUM && abs(s1 - EXP[[lab]]) <= TOL_SUM
  a3 <- TRUE; gap <- 0
  for (k in COLS) {
    x0 <- d0[[k]]; x1 <- r[[k]]
    if (!identical(is.na(x0), is.na(x1))) { a3 <- FALSE; next }
    p <- !is.na(x0) & x0 > 0
    e <- max(c(0, abs(x1[p] / x0[p] - si[p])))
    gap <- max(gap, e)
    if (e > 1e-12 || !identical(x1[!p], x0[!p])) a3 <- FALSE
  }
  cat(sprintf("A1 same table (%d rows, %d columns unchanged): %s\n", nrow(r), length(others), if (a1) "PASS" else "FAIL"))
  cat(sprintf("A2 Sigma C_ri original %.4f (V0 %.4f) ; annual %.4f (expected %.4f): %s\n", s0, V0[[lab]], s1, EXP[[lab]],
              if (a2) "PASS" else "FAIL"))
  cat(sprintf("A3 new/old = s_i, max gap %.2e over the four columns: %s\n", gap, if (a3) "PASS" else "FAIL"))
  for (k in COLS) cat(sprintf("   %-18s %.4f -> %.4f\n", k, sum(d0[[k]], na.rm = TRUE), sum(r[[k]], na.rm = TRUE)))
  if (a1 && a2 && a3) { file.rename(tmp, out); cat("written:", out, " md5", tools::md5sum(out), "\n") }
  else { unlink(tmp); ok_all <- FALSE; cat(">>> A LOCK HAS FAILED: nothing is written for", lab, "\n") }
}
if (!ok_all) { cat("\nA: FAILED\n"); quit(status = 1) }
cat("\nA: PASS\n")
