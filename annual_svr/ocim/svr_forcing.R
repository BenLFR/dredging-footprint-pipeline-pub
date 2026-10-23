#!/usr/bin/env Rscript
# svr_forcing.R -- locks of the forcing, before OCIM (rule written beforehand). Derived from d3_forcing.R.
# Exits with an error if a lock fails: the OCIM runs, afterok, do not start.
#   S0  the reference (cri_pr, lambda 100, K 16) reproduces the installed M4 forcing: logs/S0_<lab>.out, "Test P <lab> : PASS"
#   S3  mean depth of the reference = 26.0 m (smooth) / 24.3 m (flag) to within 0.1 m (test Q)
#   S1  the annual step7 finishes ("STEP 7 completed successfully") and Σ annual mass_gC_yr / reference = Q to within 1e-9,
#       Q = annual Σ C_ri / original Σ C_ri, read from the two parquets
#   S2  no "weight underflow" line in the log of the annual step7
# Outside verdict: shape TV = 1/2 Σ |m_annual/M_annual - m_ref/M_ref|; mean depth; cells with mass > 0.
# Writes BASE/Q_<lab>.txt (read by read_svr.py) and BASE/FORCING_PASS if everything passes.
# Usage: Rscript svr_forcing.R <BASE> <RUN>
suppressPackageStartupMessages({ library(data.table); library(arrow) })
a <- commandArgs(TRUE); BASE <- a[1]; RUN <- a[2]
D_Q <- c(smooth = 26.0, flag = 24.3)        # test Q (printed to 0.1 m)
TOL_MASS <- 1e-9                           # test Q0b
TOL_D <- 0.1                                # m; rounding of the test Q printout
only_one <- function(d, pat) { f <- list.files(d, pat, full.names = TRUE); if (length(f) != 1) NA_character_ else f }
read_inj <- function(d) {
  f <- only_one(d, "_injection\\.csv$"); if (is.na(f)) return(NULL)
  x <- fread(f)
  x[, .(m = sum(mass_gC_yr), z = (z_inj_base - DZT_inj / 2)[1]), by = .(i = i_ocim, j = j_ocim, k = k_inj)]
}
ok_all <- TRUE
unlink(file.path(BASE, "FORCING_PASS"))
for (lab in c("smooth", "flag")) {
  cat(sprintf("\n==================== %s ====================\n", lab))
  s0f <- file.path(BASE, "logs", sprintf("S0_%s.out", lab))
  s0 <- file.exists(s0f) && any(grepl(sprintf("Test P %s : PASS", lab), readLines(s0f, warn = FALSE), fixed = TRUE))
  cat(sprintf("S0 reproduction of the installed M4 forcing: %s\n", if (s0) "PASS" else "FAIL"))
  if (!s0) { ok_all <- FALSE; next }
  R <- read_inj(file.path(BASE, sprintf("step7_%s_ref", lab)))
  Mref <- sum(R$m); Dref <- R[m > 0, sum(m * z) / sum(m)]
  s3 <- abs(Dref - D_Q[[lab]]) <= TOL_D
  cat(sprintf("S3 mean depth of the reference: %.3f m (test Q %.1f): %s\n", Dref, D_Q[[lab]], if (s3) "PASS" else "FAIL"))
  if (!s3) { ok_all <- FALSE; next }
  f0 <- Sys.glob(file.path(RUN, "cri_pr", sprintf("cri_final_%s_*.parquet", lab)))
  f1 <- file.path(BASE, paste0("cri_an_", lab), sprintf("cri_final_%s_an_20261005.parquet", lab))
  Q <- sum(read_parquet(f1, col_select = "C_ri")$C_ri, na.rm = TRUE) / sum(read_parquet(f0, col_select = "C_ri")$C_ri, na.rm = TRUE)
  lg <- list.files(file.path(BASE, "logs"), sprintf("^s7_%s_annual_.*\\.out$", lab), full.names = TRUE)
  txt <- if (length(lg) == 1) readLines(lg, warn = FALSE) else character()
  done <- any(grepl("STEP 7 completed successfully", txt, fixed = TRUE))
  n_uf <- sum(grepl("weight underflow", txt, fixed = TRUE))
  V <- if (done) read_inj(file.path(BASE, sprintf("step7_%s_an", lab))) else NULL
  ratio <- if (is.null(V)) NA_real_ else sum(V$m) / Mref
  s1 <- !is.null(V) && abs(ratio / Q - 1) <= TOL_MASS
  s2 <- done && n_uf == 0
  cat(sprintf("S1 annual mass / reference = %.12f ; Q = %.12f ; relative gap %.2e: %s\n", ratio, Q,
              ratio / Q - 1, if (s1) "PASS" else "FAIL"))
  cat(sprintf("S2 'weight underflow' lines: %d: %s\n", n_uf, if (s2) "PASS" else "FAIL"))
  if (!(s1 && s2)) { ok_all <- FALSE; next }
  writeLines(sprintf("%.15f", Q), file.path(BASE, sprintf("Q_%s.txt", lab)))
  X <- merge(R[, .(i, j, k, mr = m)], V[, .(i, j, k, mv = m)], by = c("i", "j", "k"), all = TRUE)
  X[is.na(mr), mr := 0][is.na(mv), mv := 0]
  TV <- 0.5 * X[, sum(abs(mv / sum(mv) - mr / sum(mr)))]
  cat(sprintf("  (outside verdict) shape TV %.4f -> H-lin bound %.2f points ; depth %.3f m against %.3f ; cells > 0: %d against %d\n",
              TV, 100 * TV, V[m > 0, sum(m * z) / sum(m)], Dref, V[m > 0, .N], R[m > 0, .N]))
  writeLines(sprintf("%.15f", TV), file.path(BASE, sprintf("TV_%s.txt", lab)))
}
if (!ok_all) { cat("\nFORCING: FAILED -- OCIM does not start\n"); quit(status = 1) }
file.create(file.path(BASE, "FORCING_PASS")); cat("\nFORCING: PASS\n")
