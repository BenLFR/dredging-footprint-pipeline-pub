#!/usr/bin/env Rscript
# s_targets_20261005.R -- factor of the annual cap s_i for ALL cells with SVR >= 1, including those with zero C_ri.
# s_annuel_{lab}.csv keeps only C_ri > 0, but step7 also reads C_ri_lower, C_ri_upper and C_ri_conservative.
# s_i = sum over years of min(SVR_year, 1) / SVR, with SVR_year = sum_dw_pd / 1e6 (rule written beforehand). Elsewhere, s_i = 1 (absent).
# Inputs: cri_final_{lab}_20260925_125232.rds (md5 identical to cri_pr of _r6_flux10, checked on 2026-10-05);
# annee_{lab}.csv (written by l2_lock.R; L2 lock PASS). Run locally, R 4.5.3.
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("DFP_LOCAL"))
for (lab in c("smooth", "flag")) {
  m <- as.data.table(readRDS(sprintf("data/cri_final_%s_20260925_125232.rds", lab)))
  a <- fread(sprintf("annual_svr/annee_%s.csv", lab))
  f <- a[, .(S_annual = sum(pmin(sum_dw_pd / 1e6, 1)), S_tot = sum(sum_dw_pd / 1e6)), by = grid_id]
  c <- merge(m[!is.na(SVR) & SVR >= 1, .(grid_id, SVR, C_ri)], f, by = "grid_id")
  stopifnot(nrow(c) == nrow(f), nrow(c) == m[!is.na(SVR) & SVR >= 1, .N], max(abs(c$S_tot / c$SVR - 1)) <= 1e-9)
  c[, s := S_annual / SVR]
  stopifnot(all(c$s > 0 & c$s <= 1 + 1e-12))
  o <- sprintf("annual_svr/ocim/s_cibles_%s.csv", lab)
  fwrite(c[order(grid_id), .(grid_id, s)], o)
  cat(sprintf("%s: %d cells; s < 1: %d; sum C_ri x s / sum C_ri (targets) = %.6f; expected total %.4f\n", lab, nrow(c),
              sum(c$s < 1), sum(c$C_ri * c$s) / sum(c$C_ri), sum(m$C_ri) - sum(c$C_ri) + sum(c$C_ri * c$s)))
}
