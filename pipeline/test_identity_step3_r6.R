# test_identity_step3_r6.R -- does the modified step3 (GMM archive) give EXACTLY
# the outputs of the August step3?
# Rule written BEFORE the run:
#  PASS if, for the core file (AIS_data_core_preprocessed_V6_*.rds), same dimensions, same
#  column names in the same order, and every column identical (identical() on the values,
#  NA included, rows in the same order); AND GMM_probabilities_QA.rds identical; AND
#  drag_score_weights.csv identical byte for byte.
#  FAIL otherwise; the differing columns are listed with their number of differing rows.
#  A FAIL does not say whether it comes from the modification or from a non-reproducibility of
#  step3 itself: in both cases, step4 consumed other labels, and the replay must restart from step4.
# Usage: Rscript test_identity_step3_r6.R <dir_aug> <dir_r6>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); stopifnot(length(a) == 2)
only_one <- function(d, pat) { f <- list.files(d, pat, full.names = TRUE)
  if (length(f) != 1) { cat("FAIL:", d, pat, length(f), "files\n"); quit(status = 1) }; f }
failures <- character(0)
note <- function(ok, msg) { cat(if (isTRUE(ok)) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!isTRUE(ok)) failures <<- c(failures, msg) }
fa <- only_one(a[1], "^AIS_data_core_preprocessed_V6_.*[.]rds$"); fb <- only_one(a[2], "^AIS_data_core_preprocessed_V6_.*[.]rds$")
cat(fa, "\n", fb, "\n")
A <- readRDS(fa); B <- readRDS(fb)
note(identical(dim(A), dim(B)), sprintf("dimensions %s / %s", paste(dim(A), collapse = "x"), paste(dim(B), collapse = "x")))
note(identical(names(A), names(B)), sprintf("columns %d / %d", ncol(A), ncol(B)))
if (identical(names(A), names(B)) && nrow(A) == nrow(B)) {
  for (cc in names(A)) {
    x <- A[[cc]]; y <- B[[cc]]
    if (!identical(x, y)) {
      nd <- tryCatch(sum(!((x == y) | (is.na(x) & is.na(y))), na.rm = TRUE), error = function(e) NA)
      note(FALSE, sprintf("column %s differs (%s rows)", cc, nd))
    }
  }
  if (!length(failures)) cat("  PASS    the", ncol(A), "columns are identical\n")
}
rm(A, B); gc()
qa <- identical(readRDS(only_one(a[1], "^GMM_probabilities_QA[.]rds$")), readRDS(only_one(a[2], "^GMM_probabilities_QA[.]rds$")))
note(qa, "GMM_probabilities_QA.rds identical")
w <- identical(readBin(only_one(a[1], "^drag_score_weights[.]csv$"), "raw", 1e6), readBin(only_one(a[2], "^drag_score_weights[.]csv$"), "raw", 1e6))
note(w, "drag_score_weights.csv identical byte for byte")
if (length(failures)) { cat("\nTEST FAILED:", length(failures), "\n"); quit(status = 1) }
cat("\nTEST PASSED.\n")
