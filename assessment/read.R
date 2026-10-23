# Reading of the uncertainty tests. Order imposed by the rule: lock V0 first;
# if it fails, we stop without reading R or Omega. Test K (on _r6) has its own lock K0.
suppressMessages({ library(arrow); library(data.table) })
R6   <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r6_flux10")
BASE <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_incert")
TOT6 <- c(smooth = 4038.15, flag = 5654.39)   # _r6 totals, rounded to 0.01 t C

der <- function(d, pat) { f <- list.files(d, pat, full.names = TRUE); f[order(file.mtime(f), decreasing = TRUE)][1] }
fi  <- function(run, lab) as.data.table(read_parquet(der(file.path(run, sprintf("merge_%s_pr", lab)), "^fi_grid_.*parquet$")))
cri <- function(run, lab) as.data.table(read_parquet(der(file.path(run, "cri_pr"), sprintf("^cri_final_%s_.*parquet$", lab))))
total_cri <- function(run, lab) sum(cri(run, lab)$C_ri, na.rm = TRUE)

cat("== Lock V0: ref against _r6 ==\n")
ok <- TRUE
for (lab in c("smooth", "flag")) {
  s6 <- total_cri(R6, lab); s0 <- total_cri(file.path(BASE, "ref"), lab)
  a <- fi(R6, lab); b <- fi(file.path(BASE, "ref"), lab)
  m <- merge(a[, .(grid_id, f6 = f_i_full)], b[, .(grid_id, f0 = f_i_full)], by = "grid_id", all = TRUE)
  ids_ok <- nrow(m) == nrow(a) && nrow(m) == nrow(b)
  dmax <- m[, max(abs(f0 - f6) / pmax(abs(f6), 1e-300), na.rm = TRUE)]
  e_tot <- abs(s0 / s6 - 1)
  cat(sprintf("  %s: sum C_ri _r6 %.4f (recorded %.2f) | ref %.4f | rel gap %.2e | f_i max rel gap %.2e | grid_id identical %s\n",
              lab, s6, TOT6[lab], s0, e_tot, dmax, ids_ok))
  if (!(e_tot <= 1e-6) || !(dmax <= 1e-9) || !ids_ok || abs(s6 - TOT6[lab]) > 0.006) ok <- FALSE
}
if (!ok) { cat("LOCK V0 FAILED: neither R nor Omega is read.\n"); quit(status = 1) }
cat("LOCK V0 PASSED\n\n")

cat("== Tests R and Omega: sum C_ri(case) / sum C_ri(ref) ==\n")
res <- rbindlist(lapply(c("fill0.9", "fill0.8", "fill0.7", "om0.1", "om1.0", "om5.0"), function(c)
  rbindlist(lapply(c("smooth", "flag"), function(lab) {
    s0 <- total_cri(file.path(BASE, "ref"), lab); s <- total_cri(file.path(BASE, c), lab)
    data.table(case_id = c, label = lab, C_ri = s, ref = s0, ratio = s / s0)
  }))))
print(res[, .(case_id, label, C_ri = round(C_ri, 2), ref = round(ref, 2), ratio = round(ratio, 4))])
r07 <- res[case_id == "fill0.7", ratio]
cat(sprintf("\nTest R: R(0.7) = %s | prediction [1.3 ; 1.7]: %s | first order (>= 1.25): %s\n",
            paste(round(r07, 4), collapse = " / "),
            if (all(r07 >= 1.3 & r07 <= 1.7)) "CONFIRMED" else "REFUTED",
            if (any(r07 >= 1.25)) "YES" else "NO"))
om <- res[startsWith(case_id, "om"), abs(ratio - 1)]
cat(sprintf("Test Omega: max gap %.2f %% | %s\n\n", 100 * max(om), if (all(om < 0.05)) "SECONDARY" else "FIRST ORDER"))

cat("== Test K: k_fast fallbacks, on _r6 ==\n")
# The province -> basin table is READ from the authoritative merge (step5_merge_tiles_corr.R), not copied:
# the expression `PROV_BASIN_LUT <- data.table(...)` extracted from the file is evaluated.
MERGE <- path.expand(file.path(Sys.getenv("DFP_REPO"), "pipeline/step5_merge_tiles_corr.R"))
stopifnot(substr(tools::md5sum(MERGE), 1, 8) == "e4084757")
src <- readLines(MERGE); i0 <- grep("PROV_BASIN_LUT <- data.table\\(", src)
stopifnot(length(i0) == 1); i1 <- i0 + which(grepl("^  \\)$", src[i0:(i0 + 30)]))[1] - 1
LUTdt <- eval(parse(text = src[i0:i1])); stopifnot(nrow(LUTdt) == length(unique(LUTdt$prov_code)))
cat(sprintf("  LUT read l. %d-%d: %d provinces, %d in Antarctic\n", i0, i1, nrow(LUTdt), LUTdt[basin_raw == "Antarctic", .N]))
for (lab in c("smooth", "flag")) {
  a <- fi(R6, lab); b <- cri(R6, lab)
  stopifnot(!anyDuplicated(a$grid_id), !anyDuplicated(b$grid_id))
  orphan <- sum(!b$grid_id %in% a$grid_id)          # cells carrying C_ri without an fi_grid row
  j <- merge(a[, .(grid_id, longhurst_pr, is_fallback, lat_deg)], b[, .(grid_id, C_ri)], by = "grid_id")
  tot <- sum(j$C_ri, na.rm = TRUE)
  k0 <- orphan == 0 && abs(tot - sum(b$C_ri, na.rm = TRUE)) <= 1e-6 * tot && abs(tot - TOT6[lab]) <= 0.006
  # priority: (a) k_fast NA (fallback k = 1.0); (b) province missing from the LUT -> longitude;
  # (c) province Antarctic -> S. Pacific if lat >= -60, otherwise Arctic
  j[, fallback := fifelse(is_fallback %in% TRUE, "a_k_NA",
               fifelse(is.na(longhurst_pr) | !longhurst_pr %in% LUTdt$prov_code, "b_longitude",
               fifelse(longhurst_pr %in% LUTdt[basin_raw == "Antarctic", prov_code], "c_antarctic", "none")))]
  cat(sprintf("  %s: lock K0 %s (orphans %d, sum C_ri %.4f)\n", lab, if (k0) "PASS" else "FAIL", orphan, tot))
  if (!k0) next
  p <- j[, .(cells = .N, C_ri = sum(C_ri, na.rm = TRUE)), by = fallback][, part := round(100 * C_ri / tot, 3)]
  print(p[order(fallback)])
  part <- p[fallback != "none", sum(C_ri)] / tot
  cat(sprintf("  total share on fallback: %.3f %% -> %s\n", 100 * part,
              if (part < 0.01) "NEGLIGIBLE" else if (part < 0.05) "TO DECLARE" else "TO DECLARE + SENSITIVITY TEST"))
}
