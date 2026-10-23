# lock_flux10.R -- does the _r6 merge (10 cm junction with conserved flux)
# differ from the _r5_sansfraicheur merge ONLY by W_eff and what follows from it?
# Rule written BEFORE the run. The lock FAILS if:
#  F1. the set of grid_id or columns differs from _r5;
#  F2. a column that does not depend on W differs (> 1e-9 relative): SAR, SVR, p_d, p_l, k_fast,
#      k_used, fresh_fact, frac_sur_place, longhurst_pr;
#  F3. in cells with p_d = 1 m, W_eff deviates by more than 1e-5 from the independent Python values
#      (computed independently in Python): 0.125183 / 0.091347 / 0.081403 / 0.070564 / 0.068791 / 0.064357 /
#      0.063726 for k = 0.275 / 1 / 1.67 / 3.84 / 4.76 / 12.3 / 16.8; or a k_used there is unknown;
#  F4. p_l_eff != p_l x W_eff (> 1e-9 relative) in a cell;
#  F5. f_i_full_avant_pr and f_i_full: f(r6)/f(r5) != W_eff(r6)/W_eff(r5) (> 1e-9) where both
#      are > 0, or a cell that is zero in one and not in the other;
#  F6. W_eff(r6)/W_eff(r5) per k (p_d = 1 m) deviates by more than 5e-4 from R_K of test K2:
#      1.1448 / 1.1032 / 1.0779 / 1.0475 / 1.0434 / 1.0374 / 1.0373.
# Reported without verdict: f_i_conservative and f_i_depleted (ratios), sum of f_i_full.
# Usage: Rscript lock_flux10.R <merge_r5> <merge_new> <flag|smooth>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: lock_flux10.R <merge_r5> <merge_new> <flag|smooth>")
read_run <- function(d) {
  f <- list.files(d, "^fi_grid_.*[.]parquet$", full.names = TRUE)
  if (length(f) != 1) { cat("  FAIL    ", d, ":", length(f), "fi_grid, 1 expected\n"); quit(status = 1) }
  cat("  ", f, "\n")
  x <- as.data.table(read_parquet(f)); setkey(x, grid_id); x
}
r5 <- read_run(a[1]); nf <- read_run(a[2])
failures <- character(0)
note <- function(ok, msg) {
  cat(if (isTRUE(ok)) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!isTRUE(ok)) failures <<- c(failures, msg)
}
rel <- function(x, y) sum(abs(x - y) / pmax(abs(x), 1e-300) > 1e-9 | xor(is.na(x), is.na(y)), na.rm = TRUE)
note(identical(r5$grid_id, nf$grid_id), sprintf("F1. grid_id: %d / %d", nrow(r5), nrow(nf)))
note(identical(names(r5), names(nf)), sprintf("F1. columns: %d / %d", ncol(r5), ncol(nf)))
if (length(failures)) quit(status = 1)
for (cc in c("SAR", "SVR", "p_d", "p_l", "k_fast", "k_used", "fresh_fact", "frac_sur_place")) {
  if (!cc %in% names(r5)) { note(FALSE, sprintf("F2. column %s missing", cc)); next }
  note(rel(r5[[cc]], nf[[cc]]) == 0, sprintf("F2. %s identical (%d differing)", cc, rel(r5[[cc]], nf[[cc]])))
}
n <- sum(!(r5$longhurst_pr == nf$longhurst_pr | (is.na(r5$longhurst_pr) & is.na(nf$longhurst_pr))), na.rm = TRUE)
note(n == 0, sprintf("F2. longhurst_pr identical (%d differing)", n))
kk  <- c(0.275, 1, 1.67, 3.84, 4.76, 12.3, 16.8)
exp <- data.table(k_used = kk, W_exp = c(0.125183, 0.091347, 0.081403, 0.070564, 0.068791, 0.064357, 0.063726),
                  R_K = c(1.1448, 1.1032, 1.0779, 1.0475, 1.0434, 1.0374, 1.0373))
u <- merge(nf[abs(p_d - 1) < 1e-12 & !is.na(W_eff), .(grid_id, k_used, W_eff)],
           r5[, .(grid_id, W5 = W_eff)], by = "grid_id")
u <- merge(u, exp, by = "k_used", all.x = TRUE)
note(!anyNA(u$W_exp), sprintf("F3. every k_used with p_d = 1 is known (%d unknown)", sum(is.na(u$W_exp))))
n3 <- u[!is.na(W_exp), sum(abs(W_eff - W_exp) > 1e-5)]
note(n3 == 0, sprintf("F3. W_eff matches the independent calculation (%d / %d deviating)", n3, nrow(u)))
n4 <- nf[!is.na(p_l), sum(abs(p_l_eff - p_l * W_eff) > 1e-9 * pmax(abs(p_l_eff), 1e-300))]
note(n4 == 0, sprintf("F4. p_l_eff = p_l x W_eff (%d differing)", n4))
rW <- nf$W_eff / r5$W_eff
for (cc in c("f_i_full_avant_pr", "f_i_full")) {
  p5 <- !is.na(r5[[cc]]) & r5[[cc]] > 0; p6 <- !is.na(nf[[cc]]) & nf[[cc]] > 0
  b <- p5 & p6
  n5 <- sum(abs((nf[[cc]][b] / r5[[cc]][b]) / rW[b] - 1) > 1e-9, na.rm = TRUE) + sum(xor(p5, p6))
  note(n5 == 0, sprintf("F5. %s(r6)/%s(r5) = W_eff(r6)/W_eff(r5) (%d deviating, %d cells > 0)", cc, cc, n5, sum(b)))
}
g <- u[, .(n = .N, R = median(W_eff / W5), R_min = min(W_eff / W5), R_max = max(W_eff / W5), R_K = R_K[1]), by = k_used][order(k_used)]
print(g)
n6 <- g[!is.na(R_K), sum(abs(R_min - R_K) > 5e-4 | abs(R_max - R_K) > 5e-4)]
note(n6 == 0, sprintf("F6. W_eff(r6)/W_eff(r5) = R_K of test K2 within 5e-4 (%d k deviating)", n6))
for (cc in intersect(c("f_i_conservative", "f_i_depleted"), names(r5)))
  cat(sprintf("  [reported] sum %s: x%.5f\n", cc, sum(nf[[cc]], na.rm = TRUE) / sum(r5[[cc]], na.rm = TRUE)))
cat(sprintf("\n  sum f_i_full (%s): %.7g -> %.7g (x%.5f)\n", a[3], sum(r5$f_i_full, na.rm = TRUE),
            sum(nf$f_i_full, na.rm = TRUE), sum(nf$f_i_full, na.rm = TRUE) / sum(r5$f_i_full, na.rm = TRUE)))
if (length(failures)) { cat("\nLOCK FAILED:", length(failures), "\n"); quit(status = 1) }
cat("\nLOCK PASSED.\n")
