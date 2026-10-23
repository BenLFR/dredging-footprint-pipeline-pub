# lock_depth.R -- lock of the fix "alpha and w per pool".
#
# QUESTION: does the patched merge (G, depth_mode) do exactly what is written,
# and nothing else?
#
# RULE, written BEFORE the first run. The lock FAILS (exit 1) if:
#  M. ("fixed" mode, earlier YAML) a single column of fi_grid differs from run _r2
#     (same column names; relative gap > 1e-12 or different NA);
#  N1. ("mechanistic" mode) the set of grid_id differs from run _r2;
#  N2. SAR, SVR, p_d, p_l, k_used, fresh_fact or frac_sur_place differ from _r2 (> 1e-9);
#  N3. p_l_eff != p_l * W_eff (> 1e-9) in a single cell;
#  N4. for cells with p_d = 1 m, W_eff deviates by more than 1e-5 from the value
#      computed independently in Python (D_b 24.3; omega 0.42; fast 0.3; slow_k 0.05),
#      EXACT FORM of l (rule rewritten before the run):
#      k 0.275 -> 0.109348; 1 -> 0.082806; 1.67 -> 0.075522; 3.84 -> 0.067363;
#      4.76 -> 0.065930; 12.3 -> 0.062039; 16.8 -> 0.061435;
#      (approximate version l = sqrt(D_b/k): 0.106516;
#      0.081291; 0.074345; 0.066546; 0.065174; 0.061447; 0.060869)
#  N5. f_i_full_avant_pr(new) / f_i_full_avant_pr(_r2) != p_l_eff(new) / p_l_eff(_r2)
#      (> 1e-9) in a cell where both are > 0.
#
# Usage: Rscript lock_depth.R <merge_r2> <merge_fixed_test> <merge_mechanistic>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: lock_depth.R <merge_r2> <merge_fixed> <merge_mech>")
read_run <- function(d) {
  f <- list.files(d, "^fi_grid_.*[.]parquet$", full.names = TRUE)
  if (length(f) != 1) { cat("  FAIL    ", d, ":", length(f), "fi_grid, 1 expected\n"); quit(status = 1) }
  cat("  ", f, "\n")
  x <- as.data.table(read_parquet(f)); setkey(x, grid_id); x
}
r2 <- read_run(a[1]); fx <- read_run(a[2]); mc <- read_run(a[3])
failures <- character(0)
note <- function(ok, msg) {
  cat(if (isTRUE(ok)) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!isTRUE(ok)) failures <<- c(failures, msg)
}
gap <- function(x, y) {
  if (is.numeric(x) && is.numeric(y))
    sum(abs(x - y) / pmax(abs(x), 1e-300) > 1e-12 | xor(is.na(x), is.na(y)), na.rm = TRUE)
  else sum(!(x == y | (is.na(x) & is.na(y))), na.rm = TRUE)
}

cat("\n--- M. fixed mode against _r2 ---\n")
note(identical(names(r2), names(fx)), sprintf("M. same columns (%d / %d)", ncol(r2), ncol(fx)))
note(identical(r2$grid_id, fx$grid_id), "M. same grid_id")
if (!length(failures)) for (cc in names(r2)) {
  n <- gap(r2[[cc]], fx[[cc]])
  if (n > 0) note(FALSE, sprintf("M. column %s: %d cell(s) deviating", cc, n))
}
if (!length(failures)) note(TRUE, "M. all columns identical")

cat("\n--- N. mechanistic mode against _r2 ---\n")
note(identical(r2$grid_id, mc$grid_id), sprintf("N1. grid_id: %d / %d", nrow(r2), nrow(mc)))
if (length(failures)) quit(status = 1)
for (cc in c("SAR", "SVR", "p_d", "p_l", "k_used", "fresh_fact", "frac_sur_place")) {
  x <- r2[[cc]]; y <- mc[[cc]]
  n <- sum(abs(x - y) / pmax(abs(x), 1e-300) > 1e-9 | xor(is.na(x), is.na(y)), na.rm = TRUE)
  note(n == 0, sprintf("N2. %s identical (%d differing)", cc, n))
}
n3 <- mc[!is.na(p_l), sum(abs(p_l_eff - p_l * W_eff) > 1e-9 * pmax(abs(p_l_eff), 1e-300))]
note(n3 == 0, sprintf("N3. p_l_eff = p_l x W_eff (%d differing)", n3))
exp <- data.table(k_used = c(0.275, 1, 1.67, 3.84, 4.76, 12.3, 16.8),
                  W_exp  = c(0.125183, 0.091347, 0.081403, 0.070564, 0.068791, 0.064357, 0.063726))  # junction with conserved flux (before: 0.109348, 0.082806, 0.075522, 0.067363, 0.065930, 0.062039, 0.061435)
u <- merge(mc[abs(p_d - 1) < 1e-12 & !is.na(W_eff), .(grid_id, k_used, W_eff)], exp, by = "k_used", all.x = TRUE)
note(!anyNA(u$W_exp), sprintf("N4. every k_used with p_d = 1 is known (%d unknown: %s)",
                               sum(is.na(u$W_exp)), paste(unique(u[is.na(W_exp), k_used]), collapse = ",")))
n4 <- u[!is.na(W_exp), sum(abs(W_eff - W_exp) > 1e-5)]
note(n4 == 0, sprintf("N4. W_eff matches the independent calculation (%d / %d deviating)", n4, nrow(u)))
cat(sprintf("     cells with p_d = 1 m: %d / %d\n", nrow(u), nrow(mc)))
col_f <- if ("f_i_full_avant_pr" %in% names(mc)) "f_i_full_avant_pr" else "f_i_full"
x <- r2[[col_f]]; y <- mc[[col_f]]; ok <- !is.na(x) & !is.na(y) & x > 0 & y > 0
rf <- y[ok] / x[ok]; rp <- mc$p_l_eff[ok] / r2$p_l_eff[ok]
n5 <- sum(abs(rf - rp) > 1e-9 * pmax(abs(rp), 1e-300))
note(n5 == 0, sprintf("N5. %s follows p_l_eff (%d / %d differing)", col_f, n5, sum(ok)))

cat("\n--- what changes ---\n")
cat(sprintf("  sum f_i_full        : %.6g -> %.6g  (x%.4f)\n", sum(r2$f_i_full, na.rm = TRUE),
            sum(mc$f_i_full, na.rm = TRUE), sum(mc$f_i_full, na.rm = TRUE) / sum(r2$f_i_full, na.rm = TRUE)))
cat(sprintf("  sum f_i_conservative: x%.4f\n",
            sum(mc$f_i_conservative, na.rm = TRUE) / sum(r2$f_i_conservative, na.rm = TRUE)))
cat(sprintf("  cells with f_i_full = 1 (capped): %d -> %d\n",
            sum(r2$f_i_full >= 1, na.rm = TRUE), sum(mc$f_i_full >= 1, na.rm = TRUE)))
print(mc[!is.na(SVR) & SVR > 0, .(n = .N, W_eff_med = median(W_eff),
         sum_fi_new = sum(f_i_full, na.rm = TRUE)), by = k_used][order(k_used)])

if (length(failures)) { cat("\nLOCK FAILED --", length(failures), "condition(s).\n"); quit(status = 1) }
cat("\nLOCK PASSED.\n")
