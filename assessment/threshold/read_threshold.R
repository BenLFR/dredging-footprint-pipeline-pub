#!/usr/bin/env Rscript
# read_threshold.R -- Classification threshold +/- 0.005 (Label B): reading (rule written beforehand).
#
# Mode "p2"   : lock P2 only (merge + step6 of _r6 rerun on the step5 tiles of _r2). Writes BASE/P2_PASS.
#                 Step 2 (step5 of the three arms) starts only if P2_PASS exists.
# Mode "all": locks, in order, then substance. The first lock that fails stops the reading.
#   Pref  the _r6 flag reference sums to 5 654.39 t C to within 0.005: units and yardstick G_C applicable
#   M     each merge read what we believe: "table  :" in the log of the merges with table,
#         "FI_FRAC_FILE not set" in the log of the merges without table (step5_merge_tiles_corr.R)
#   P2    C_ri(BASE/p2) = C_ri(_r6 flag), same cells, cell by cell to 1e-9 relative
#   P3    C_ri(BASE/t0, with table) = C_ri(_r6 flag), same (step5 deterministic)
# Substance, for tp and tm, against t0 (cell by cell, full join, missing = 0):
#   reading 1 (fixed fraction): C1 = C_ri with table                     ; dC = sum C1_var - sum C1_t0
#   reading 2 (fixed hopper)   : C2 = max(0, C1_t0 + C0_var - C0_t0)      ; dC = sum C2     - sum C1_t0
#                               (C0 = C_ri without table, merge without FI_FRAC_FILE)
#   R = |dC| / G_C, G_C = 5 654.39 - 4 038.15 = 1 616.24 t C. Cut-points shared with the other uncertainty tests:
#   R <= 0.10 negligible; 0.10 < R < 0.50 reported; R >= 0.50 first order, FAILS.
#   Verdict only if both readings fall in the same class; otherwise "indeterminate" and the bracket.
# Outside verdict: sums per arm and reading, counts of positions (compteurs.csv), share of C in the new
#   cells, D = q95/q05 of s = C1_var/C1_t0 weighted by C1_t0 (common cells, reading 1),
#   C_ri_conservative if the column exists.
# Usage: Rscript read_threshold.R <BASE> <RUN_r6> <p2|all>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(TRUE); BASE <- a[1]; RUN <- a[2]; MODE <- a[3]
TOL   <- 1e-9                        # P2, P3; Q0b threshold
S_REF <- 5654.39                     # sum C_ri _r6 flag, t C
G_C   <- 5654.39 - 4038.15           # t C
ARMS  <- c("tp", "tm")

cri <- function(d) {
  x <- list.files(file.path(d, "cri_pr"), "^cri_final_flag_.*\\.parquet$", full.names = TRUE)
  if (length(x) != 1) stop(d, "/cri_pr: ", length(x), " cri_final_flag, 1 expected")
  y <- as.data.table(read_parquet(x))
  if (!all(c("grid_id", "C_ri") %in% names(y))) stop(x, ": columns grid_id / C_ri missing")
  if (anyDuplicated(y$grid_id)) stop(x, ": grid_id duplicated")
  y
}
lock <- function(item_name, ok, detail) {
  cat(sprintf("  %-5s %-44s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", item_name, detail))
  if (!isTRUE(ok)) { cat("\n>>> LOCK FAILED: the reading stops here. Nothing of the substance is read.\n"); quit(status = 1) }
}
check_identical <- function(item_name, r, v) {
  X <- merge(r[, .(grid_id, r = C_ri)], v[, .(grid_id, v = C_ri)], by = "grid_id", all = TRUE)
  same <- X[, !anyNA(r) && !anyNA(v)]
  gap <- if (same) X[, max(abs(v - r) / pmax(abs(r), abs(v), 1e-300))] else NA_real_
  lock(item_name, same && gap <= TOL,
         sprintf("cells %d / %d ; outside reference %d ; missing %d ; max relative gap %.2e ; sum %.4f / %.4f",
                 nrow(v), nrow(r), X[is.na(r), .N], X[is.na(v), .N], gap, sum(v$C_ri), sum(r$C_ri)))
}
log_merge <- function(item_name, pat) {
  f <- list.files(file.path(BASE, "logs"), sprintf("^m_%s_.*\\.out$", item_name), full.names = TRUE)
  ok <- length(f) == 1 && any(grepl(pat, readLines(f, warn = FALSE), fixed = TRUE))
  lock(sprintf("M merge log %s", item_name), ok, sprintf("'%s' (%d log)", pat, length(f)))
}

cat("==== read_threshold.R ", MODE, " ", format(Sys.time()), " ====\n")
REF <- cri(RUN)
lock("Pref sum _r6 flag = 5654.39 t C", abs(sum(REF$C_ri) - S_REF) <= 0.005, sprintf("%.4f", sum(REF$C_ri)))
log_merge("p2", "table  :")
check_identical("P2 chain merge + step6 (tiles _r2)", REF, cri(file.path(BASE, "p2")))
if (MODE == "p2") { writeLines(c("P2_PASS", format(Sys.time())), file.path(BASE, "P2_PASS")); cat("P2_PASS\n"); quit(status = 0) }
stopifnot(MODE == "all")

for (b in c("t0", ARMS)) { log_merge(paste0(b, "_with"), "table  :"); log_merge(paste0(b, "_without"), "FI_FRAC_FILE not set") }
C1 <- list(t0 = cri(file.path(BASE, "t0"))); C0 <- list(t0 = cri(file.path(BASE, "t0_sansfrac")))
check_identical("P3 arm t0 complete (step5 rerun)", REF, C1$t0)
for (b in ARMS) { C1[[b]] <- cri(file.path(BASE, b)); C0[[b]] <- cri(file.path(BASE, paste0(b, "_sansfrac"))) }

category <- function(R) if (R <= 0.10) "negligible" else if (R < 0.50) "second order, reported" else "FIRST ORDER -> FAIL"
wq <- function(s, w, p) { o <- order(s); cw <- cumsum(w[o]) / sum(w); s[o][which(cw >= p)[1]] }
counts <- fread(file.path(BASE, "compteurs.csv"))
cat("\n---- position counts (not part of the verdict; earlier measurement: +4.02 % / -3.17 %, smoothing not specified) ----\n")
print(counts)
s1r <- sum(C1$t0$C_ri)
cat(sprintf("\n---- substance: reference t0, sum C1 = %.2f t C ; sum C0 = %.2f t C ; G_C = %.2f t C ----\n",
            s1r, sum(C0$t0$C_ri), G_C))
res <- list()
for (b in ARMS) {
  X <- Reduce(function(x, y) merge(x, y, by = "grid_id", all = TRUE),
              list(C1$t0[, .(grid_id, c1r = C_ri)], C0$t0[, .(grid_id, c0r = C_ri)],
                   C1[[b]][, .(grid_id, c1v = C_ri)], C0[[b]][, .(grid_id, c0v = C_ri)]))
  is_new <- X[, is.na(c1r) & is.na(c0r)]
  pres_v <- X[, !is.na(c1v)]                     # cell present in the arm (for D, "present in both")
  for (cc in c("c1r", "c0r", "c1v", "c0v")) set(X, which(is.na(X[[cc]])), cc, 0)
  X[, c2 := pmax(0, c1r + c0v - c0r)]
  d1 <- sum(X$c1v) - s1r; d2 <- sum(X$c2) - s1r
  R1 <- abs(d1) / G_C; R2 <- abs(d2) / G_C
  k1 <- category(R1); k2 <- category(R2)
  verdict <- if (k1 == k2) k1 else "INDETERMINATE (the two readings fall in two classes)"
  # P1b (revised rule P1): positions classified outside the 71 enlarged tiles -> Delta incomplete
  nh <- if ("hors_71" %in% names(counts)) counts[bras == b, hors_71] else NA_integer_
  if (is.na(nh) || nh > 0L) verdict <- sprintf("%s ; INCOMPLETE (P1b: %s positions classified outside the 71 tiles)", verdict, nh)
  cat(sprintf("\n==== %s (threshold %s) ====\n", b, counts[bras == b, seuil]))
  cat(sprintf("  reading 1 (fixed fraction): sum C = %.2f t C ; dC = %+.2f t C (%+.2f %%) ; R = %.4f -> %s\n",
              sum(X$c1v), d1, 100 * d1 / s1r, R1, k1))
  cat(sprintf("  reading 2 (fixed hopper)  : sum C = %.2f t C ; dC = %+.2f t C (%+.2f %%) ; R = %.4f -> %s\n",
              sum(X$c2), d2, 100 * d2 / s1r, R2, k2))
  cat(sprintf("  VERDICT: %s\n", verdict))
  # outside verdict
  gap_new <- X[is_new, max(abs(c2 - c1v), 0)]
  cat(sprintf("  (outside verdict) new cells: %d ; share of C: %.3f %% (reading 1), %.3f %% (reading 2) ; |C2 - C1| max over them %.2e (0 expected)\n",
              sum(is_new), 100 * X[is_new, sum(c1v)] / sum(X$c1v), 100 * X[is_new, sum(c2)] / sum(X$c2), gap_new))
  cat(sprintf("  (not part of the verdict) lost cells (C1 t0 > 0, C1 %s = 0): %d, %.3f %% of C1 t0\n",
              b, X[c1r > 0 & c1v == 0, .N], 100 * X[c1r > 0 & c1v == 0, sum(c1r)] / s1r))
  Y <- X[!is_new & pres_v & c1r > 0]
  s <- Y$c1v / Y$c1r
  q05 <- wq(s, Y$c1r, 0.05); q95 <- wq(s, Y$c1r, 0.95); D <- q95 / q05
  cat(sprintf("  (outside verdict) D = q95/q05 of s (reading 1, weighted by C1 t0) = %.4f / %.4f = %.4f -> %s\n",
              q95, q05, D, if (!is.finite(D)) "not finite (q05 = 0)" else if (D <= 1.10) "rescaling (<= 1.10)"
              else if (D > 1.50) "deformation (> 1.50)" else "between 1.10 and 1.50"))
  cons <- NA_real_
  if ("C_ri_conservative" %in% names(C1[[b]]) && "C_ri_conservative" %in% names(C1$t0)) {
    cons <- sum(C1[[b]]$C_ri_conservative) - sum(C1$t0$C_ri_conservative)
    cat(sprintf("  (not part of the verdict) C_ri_conservative, reading 1: dC = %+.2f t C\n", cons))
  }
  res[[b]] <- data.table(bras = b, seuil = counts[bras == b, seuil], sum_C1_t0 = s1r,
                         sum_L1 = sum(X$c1v), dC_L1 = d1, R_L1 = R1, class_L1 = k1,
                         sum_L2 = sum(X$c2), dC_L2 = d2, R_L2 = R2, class_L2 = k2, verdict = verdict,
                         hors_71 = nh, new_cells = sum(is_new), share_new_L1 = X[is_new, sum(c1v)] / sum(X$c1v),
                         D_L1 = D, dC_conservative_L1 = cons)
}
fwrite(rbindlist(res), file.path(BASE, "threshold_result.csv"))
cat("\nEND -- threshold_result.csv written\n")
