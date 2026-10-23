# lock_sala_merge.R -- guard of chain_step4_ocim.sh after the merge.
#
# QUESTION: between the August run and the sala2021 replay, has only the
# lability term changed? This is what allows reusing the August FI_FRAC_FILE
# (fraction left on site), computed from SAR, SVR and sum_d
# (phase3_fraction_rejeu.R) -- quantities that do not depend on lithology.
#
# RULE, written before the first run. The lock FAILS (exit 1) if:
#   1. the set of grid_id differs;
#   2. SAR, SVR, sum_d, p_d, n_with_pl, fresh_fact, k_used or frac_sur_place
#      differ (relative gap > 1e-9) in a single cell;
#   3. p_l is identical everywhere -- the scheme would not have been applied.
#
# Usage: Rscript lock_sala_merge.R <dir_merge_reference> <dir_merge_rerun>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("usage: lock_sala_merge.R <dir_reference> <dir_rerun>")
read_run <- function(d) {
  f <- list.files(d, "^fi_grid_.*\\.parquet$", full.names = TRUE)
  if (length(f) != 1) { cat("  FAIL    ", d, ":", length(f), "fi_grid, 1 expected\n"); quit(status = 1) }
  cat("  ", f, "\n")
  as.data.table(read_parquet(f))
}
ref <- read_run(a[1]); new <- read_run(a[2])
failures <- character(0)
note <- function(ok, msg) {
  cat(if (ok) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!ok) failures <<- c(failures, msg)
}

note(setequal(ref$grid_id, new$grid_id) && nrow(ref) == nrow(new),
     sprintf("1. grid_id: %d / %d cells", nrow(ref), nrow(new)))
if (length(failures)) quit(status = 1)
m <- merge(ref, new, by = "grid_id", suffixes = c(".ref", ".new"))

for (cc in c("SAR", "SVR", "sum_d", "p_d", "n_with_pl", "fresh_fact", "k_used", "frac_sur_place")) {
  x <- m[[paste0(cc, ".ref")]]; y <- m[[paste0(cc, ".new")]]
  gap <- abs(x - y) / pmax(abs(x), 1e-300)
  n_bad <- sum(gap > 1e-9 | xor(is.na(x), is.na(y)), na.rm = TRUE)
  note(n_bad == 0, sprintf("2. %s identical (%d cell(s) deviating)", cc, n_bad))
}
note(any(abs(m$p_l.ref - m$p_l.new) > 1e-12, na.rm = TRUE), "3. p_l has changed")

cat("\n--- what changes (sums over the cells) ---\n")
cat(sprintf("  mean p_l       : %.4f -> %.4f\n", mean(m$p_l.ref, na.rm = TRUE), mean(m$p_l.new, na.rm = TRUE)))
cat(sprintf("  sum f_i_full: %.4f -> %.4f  (x%.3f)\n", sum(m$f_i_full.ref, na.rm = TRUE),
            sum(m$f_i_full.new, na.rm = TRUE),
            sum(m$f_i_full.new, na.rm = TRUE) / sum(m$f_i_full.ref, na.rm = TRUE)))

if (length(failures)) {
  cat("\nLOCK FAILED --", length(failures), "condition(s). The chain stops here.\n")
  quit(status = 1)
}
cat("\nLOCK PASSED.\n")
