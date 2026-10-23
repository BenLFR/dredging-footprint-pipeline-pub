# summary_sala.R -- last R link of chain_step4_ocim.sh: what switching to the
# Sala et al. (2021) table changes in the remineralised carbon C_ri.
# Not a lock: it prints and refuses nothing (exit 0 unless a
# file is missing). The output is read afterwards and recorded in the project log.
#
# Usage: Rscript summary_sala.R <reference_cri_dir> <rerun_cri_dir>
suppressMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE)
for (lab in c("smooth", "flag")) {
  f <- function(d) {
    x <- list.files(d, sprintf("^cri_final_%s_.*\\.parquet$", lab), full.names = TRUE)
    if (length(x) != 1) stop(d, ": ", length(x), " cri_final_", lab, ", 1 expected")
    as.data.table(read_parquet(x))
  }
  r <- f(a[1]); n <- f(a[2])
  cat(sprintf("\n=== %s ===\n", lab))
  for (cc in intersect(c("C_ri", "C_ri_lower", "C_ri_upper", "C_ri_conservative"), names(n))) {
    s0 <- sum(r[[cc]], na.rm = TRUE); s1 <- sum(n[[cc]], na.rm = TRUE)
    cat(sprintf("  %-18s %14.2f -> %14.2f   (x%.3f)\n", cc, s0, s1, s1 / s0))
  }
  cat(sprintf("  cells: %d -> %d\n", nrow(r), nrow(n)))
}
