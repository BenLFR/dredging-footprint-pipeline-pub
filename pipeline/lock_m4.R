# lock_m4.R -- test P, rule written before the replay.
# The forcing of the M4 replay (step7 of the pipeline, md5 c3ab43c0) must be IDENTICAL to that of the
# candidate tested in test O (md5 5fa60e3f, same C_ri _r6): the two versions only differ by the
# computation of LOCK1, which touches no written number.
# Fails if: number of Jdredge vectors != 12; different support; max relative gap > 1e-12;
# injection table (i, j, k_inj, mass) different beyond 1e-12 relative.
# 1e-12 threshold: taken from the earlier O1 comparison; modelling choice, no source.
# Usage: Rscript lock_m4.R <candidate .mat directory> <replay step7 directory> <label>
.libPaths("~/R/library")
suppressMessages(library(data.table))
a <- commandArgs(TRUE); CAND <- a[1]; RERUN <- a[2]; LAB <- a[3]
only_mat <- function(d) { f <- list.files(d, "^jdredge_ocim2_48l_.*\\.mat$", full.names = TRUE)
  if (length(f) != 1) stop(d, ": ", length(f), " .mat, 1 expected"); f }
fc <- only_mat(CAND); fr <- only_mat(RERUN)
cat("=== Test P,", LAB, "===\n  candidate:", fc, "\n  rerun    :", fr, "\n")
A <- R.matlab::readMat(fc); B <- R.matlab::readMat(fr)
nms <- grep("^Jdredge", names(B), value = TRUE)
ok <- length(nms) == 12 && setequal(nms, grep("^Jdredge", names(A), value = TRUE))
cat(sprintf("  Jdredge vectors: %d (12 expected), same names: %s\n", length(nms), ok))
for (n in nms) {
  x <- as.numeric(A[[n]]); y <- as.numeric(B[[n]])
  same_supp <- length(x) == length(y) && identical(which(x != 0), which(y != 0))
  gap <- if (same_supp && any(y != 0)) max(abs(x - y)[y != 0] / abs(y[y != 0])) else if (same_supp) 0 else NA
  good <- same_supp && gap <= 1e-12
  ok <- ok && good
  cat(sprintf("  %-26s identical support %-5s max rel gap %.2e %s\n", n, same_supp, gap, if (good) "ok" else "GAP"))
}
ic <- fread(sub("\\.mat$", "_injection.csv", fc)); ir <- fread(sub("\\.mat$", "_injection.csv", fr))
ag <- function(I) I[, .(m = sum(mass_gC_yr)), by = .(i = i_ocim, j = j_ocim, k = k_inj)]
X <- merge(ag(ic)[, .(i, j, k, mc = m)], ag(ir)[, .(i, j, k, mr = m)], by = c("i", "j", "k"), all = TRUE)
inj_ok <- !anyNA(X) && X[, max(abs(mc - mr) / pmax(abs(mc), 1e-300))] <= 1e-12
cat(sprintf("  injection: %d cells (i,j,k) candidate, %d rerun, union %d, identical: %s\n",
            nrow(ag(ic)), nrow(ag(ir)), nrow(X), inj_ok))
ok <- ok && inj_ok
cat("Test P", LAB, ":", if (ok) "PASS" else "FAIL", "\n")
if (!ok) quit(status = 3)
