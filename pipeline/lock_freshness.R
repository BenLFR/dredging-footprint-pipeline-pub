# lock_freshness.R -- does the merge WITHOUT freshness correction
# (decision C) differ from the _r4_exact merge ONLY by fresh_fact?
# Rule written BEFORE. The lock FAILS if:
#  F1. the set of grid_id or columns differs from _r4_exact;
#  F2. a column that does not depend on fresh_fact differs (> 1e-9 relative): SAR, SVR, p_d, p_l,
#      [v2: W_eff_cons removed, absent from the parquet; v1 FAILED on this point only]
#      k_used, W_eff, p_l_eff, frac_sur_place, longhurst_pr;
#  F3. fresh_fact is not 1.0 in EVERY cell of the new merge;
#  F4. an f_i* column does not satisfy f_i(new) x fresh_fact(_r4) = f_i(_r4) (> 1e-9 relative);
#  F5. the ratio of the sums of f_i_full deviates by more than 1e-3 from that of test G:
#      1 - 0.0572 (flag), 1 - 0.0529 (smooth).
# Usage: Rscript lock_freshness.R <merge_r4> <merge_new> <flag|smooth>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: lock_freshness.R <merge_r4> <merge_new> <flag|smooth>")
read_run <- function(d) {
  f <- list.files(d, "^fi_grid_.*[.]parquet$", full.names = TRUE)
  if (length(f) != 1) { cat("  FAIL    ", d, ":", length(f), "fi_grid, 1 expected\n"); quit(status = 1) }
  cat("  ", f, "\n")
  x <- as.data.table(read_parquet(f)); setkey(x, grid_id); x
}
r4 <- read_run(a[1]); nf <- read_run(a[2])
failures <- character(0)
note <- function(ok, msg) {
  cat(if (isTRUE(ok)) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!isTRUE(ok)) failures <<- c(failures, msg)
}
rel <- function(x, y) sum(abs(x - y) / pmax(abs(x), 1e-300) > 1e-9 | xor(is.na(x), is.na(y)), na.rm = TRUE)
note(identical(r4$grid_id, nf$grid_id), sprintf("F1. grid_id: %d / %d", nrow(r4), nrow(nf)))
note(identical(names(r4), names(nf)), sprintf("F1. columns: %d / %d", ncol(r4), ncol(nf)))
if (length(failures)) quit(status = 1)
for (cc in c("SAR", "SVR", "p_d", "p_l", "k_used", "W_eff", "p_l_eff", "frac_sur_place")) {
  if (!cc %in% names(r4)) { note(FALSE, sprintf("F2. column %s missing", cc)); next }
  note(rel(r4[[cc]], nf[[cc]]) == 0, sprintf("F2. %s identical (%d differing)", cc, rel(r4[[cc]], nf[[cc]])))
}
n <- sum(!(r4$longhurst_pr == nf$longhurst_pr | (is.na(r4$longhurst_pr) & is.na(nf$longhurst_pr))), na.rm = TRUE)
note(n == 0, sprintf("F2. longhurst_pr identical (%d differing)", n))
note(all(nf$fresh_fact == 1), sprintf("F3. fresh_fact = 1 everywhere (%d other cells)", sum(nf$fresh_fact != 1 | is.na(nf$fresh_fact))))
cat(sprintf("     fresh_fact in _r4: %s\n", paste(names(table(r4$fresh_fact)), table(r4$fresh_fact), sep = " x", collapse = " ; ")))
for (cc in grep("^f_i", names(r4), value = TRUE)) {
  n4 <- rel(r4[[cc]], nf[[cc]] * r4$fresh_fact)
  note(n4 == 0, sprintf("F4. %s(new) x fresh_fact(_r4) = %s(_r4) (%d differing)", cc, cc, n4))
}
exp <- c(flag = 1 - 0.0572, smooth = 1 - 0.0529)[a[3]]
rs <- sum(nf$f_i_full, na.rm = TRUE) / sum(r4$f_i_full, na.rm = TRUE)
note(!is.na(exp) && abs(rs - exp) < 1e-3, sprintf("F5. sum of f_i_full: x%.5f (expected x%.4f, test G)", rs, exp))
cat(sprintf("\n  sum f_i_full: %.7g -> %.7g\n", sum(r4$f_i_full, na.rm = TRUE), sum(nf$f_i_full, na.rm = TRUE)))
if (length(failures)) { cat("\nLOCK FAILED:", length(failures), "\n"); quit(status = 1) }
cat("\nLOCK PASSED.\n")
