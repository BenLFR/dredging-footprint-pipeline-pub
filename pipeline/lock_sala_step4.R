# lock_sala_step4.R -- guard of chain_step4_ocim.sh after step4.
#
# QUESTION: did the step4 replay under PL_SCHEME=sala2021 change ONLY
# pl_base?
#
# RULE, written before the first run. The lock FAILS (exit 1, the
# chain stops: everything after it is afterok) if one of these
# conditions is not met:
#   1. a single AIS_with_lithology_*.rds in the replay directory;
#   2. same number of rows as the reference output;
#   3. identical columns, exactly, apart from pl_base: lithologie_canon,
#      has_lithology, is_hardrock, H_index, frac_mud/sand/gravel, and
#      use_for_carbon if it exists;
#   4. where the hard-bottom blend is inactive (H_index <= 0.70 or NA),
#      pl_base only takes the values of the scheme: 0.7, 0.286, 0.04, 0.05,
#      0.10, 0.50;
#   5. pl_base differs from the reference on at least one row -- otherwise the
#      scheme was not applied and the replay is useless.
#
# Usage: Rscript lock_sala_step4.R <reference_file.rds> <replay_directory>
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("usage: lock_sala_step4.R <reference.rds> <dir_rerun>")
failures <- character(0)
note <- function(ok, msg) {
  cat(if (ok) "  PASS   " else "  FAIL   ", msg, "\n")
  if (!ok) failures <<- c(failures, msg)
}

new_items <- list.files(a[2], "^AIS_with_lithology_.*\\.rds$", full.names = TRUE)
note(length(new_items) == 1, sprintf("1. a single step4 file in the replay (%d found)", length(new_items)))
if (length(new_items) != 1) quit(status = 1)
cat("reference:", a[1], "\nrerun     :", new_items, "\n\n")
ref <- as.data.table(readRDS(a[1]))
new <- as.data.table(readRDS(new_items))

note(nrow(ref) == nrow(new), sprintf("2. rows: %d / %d", nrow(ref), nrow(new)))
if (nrow(ref) != nrow(new)) quit(status = 1)

cols <- intersect(c("lithologie_canon", "has_lithology", "is_hardrock", "H_index",
                    "frac_mud", "frac_sand", "frac_gravel", "use_for_carbon"), names(ref))
for (cc in cols) {
  if (!cc %in% names(new)) { note(FALSE, paste("3. column missing from the replay:", cc)); next }
  note(identical(ref[[cc]], new[[cc]]), paste("3. identical:", cc))
}
note(setequal(names(ref), names(new)), "3. same set of columns")

outside_blend <- new[, is.na(H_index) | H_index <= 0.70]
values <- sort(unique(round(new$pl_base[outside_blend], 6)))
note(all(values %in% c(0.7, 0.286, 0.04, 0.05, 0.10, 0.50)),
     paste("4. pl_base outside the blend:", paste(values, collapse = " ")))

n_diff <- sum(abs(ref$pl_base - new$pl_base) > 1e-12, na.rm = TRUE)
note(n_diff > 0, sprintf("5. pl_base changed on %d rows (%.1f %%)", n_diff, 100 * n_diff / nrow(new)))

cat("\n--- mean pl_base per class: table_v1 -> sala2021 (pings) ---\n")
print(data.table(lithologie_canon = ref$lithologie_canon, old = ref$pl_base, newer = new$pl_base)[,
  .(N = .N, old = round(mean(old, na.rm = TRUE), 3), newer = round(mean(newer, na.rm = TRUE), 3)),
  by = lithologie_canon][order(-N)])

if (length(failures)) {
  cat("\nLOCK FAILED --", length(failures), "condition(s). The chain stops here.\n")
  quit(status = 1)
}
cat("\nLOCK PASSED.\n")
