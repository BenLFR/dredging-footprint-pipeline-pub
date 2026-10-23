# d3_forcing.R -- D3 under M4, STEP S (rule written beforehand).
# Measures, without OCIM, what lambda and K change in the M4 forcing, and decides which variants go to OCIM.
# Order imposed: locks S0-S3 first; a variant that fails at S1 or S2 is not measured.
#   S0  the variant (100,16) reproduces the installed M4 forcing: read in logs/S0_<lab>.out (lock_m4.R, test P)
#   S1  the variant finished without error ("STEP 7 completed successfully") and its mass = that of the
#       reference to 1e-9 relative (Q0b threshold)
#   S2  no "weight underflow" line in its log (otherwise uniform fallback, a confounder)
#   S3  mean depth of the reference = 26.0 m (smooth) / 24.3 m (flag) to within 0.1 m (test Q)
# Measure: TV = 1/2 sum|m_var - m_ref| / sum m_ref over the injection table (cell i,j,k).
# Bound (H-lin hypothesis, not proven): |dAF| <= TV x 100 points.
# OCIM choice (rule written beforehand): (500,16) for both labels, always; for K, per label, the
#   variant with the largest TV if TV > 0.01. Writes BASE/choix_ocim.txt and BASE/etape_S.csv.
# Usage: Rscript d3_forcing.R <BASE> <RUN>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); BASE <- a[1]; RUN <- a[2]
LABS <- c("smooth", "flag")
TAGS <- c("l100_k16", "l500_k16", "l100_k1", "l100_k4", "l100_k64")
D_Q  <- c(smooth = 26.0, flag = 24.3)       # test Q (printed to 0.1 m)
TOL_MASS <- 1e-9                            # Q0b; design choice
TOL_D     <- 0.1                             # m; rounding of the test Q printout
THRESHOLD_TV  <- 0.01                            # bound 1 point; design choice
G_AF      <- 54.4364 - 49.8915               # points; test Q

only_one <- function(d, pat) { f <- list.files(d, pat, full.names = TRUE)
  if (length(f) != 1) return(NA_character_); f }
logs <- function(lab, tag) list.files(file.path(BASE, "logs"), sprintf("^s7_%s_%s_.*\\.out$", lab, tag), full.names = TRUE)
read_inj <- function(lab, tag) {
  f <- only_one(file.path(BASE, sprintf("step7_%s_%s", lab, tag)), "_injection\\.csv$")
  if (is.na(f)) return(NULL)
  x <- fread(f)
  x[, .(m = sum(mass_gC_yr), z = (z_inj_base - DZT_inj / 2)[1]), by = .(i = i_ocim, j = j_ocim, k = k_inj)]
}

res <- list(); choice <- character(); stop <- FALSE
for (lab in LABS) {
  cat(sprintf("\n==================== %s ====================\n", lab))
  s0 <- file.path(BASE, "logs", sprintf("S0_%s.out", lab))
  s0_ok <- file.exists(s0) && any(grepl(sprintf("Test P %s : PASS", lab), readLines(s0, warn = FALSE), fixed = TRUE))
  cat(sprintf("S0 reproduction of the installed M4 forcing: %s\n", if (s0_ok) "PASS" else "FAIL"))
  if (!s0_ok) { stop <- TRUE; next }
  R <- read_inj(lab, "l100_k16")
  Mref <- sum(R$m)
  Dref <- R[m > 0, sum(m * z) / sum(m)]
  s3_ok <- abs(Dref - D_Q[lab]) <= TOL_D
  cat(sprintf("S3 mean depth of the reference: %.3f m (test Q %.1f) %s\n", Dref, D_Q[lab], if (s3_ok) "PASS" else "FAIL"))
  if (!s3_ok) { stop <- TRUE; next }
  # outside verdict: settles (or not) the 208/261 versus 297/340 question of test Q
  cat(sprintf("  (outside verdict) rows of the injection table: %d ; with mass > 0: %d\n", nrow(R), R[m > 0, .N]))
  colR <- R[, .(m = sum(m)), by = .(i, j)][order(-m)]
  for (tag in TAGS) {
    lg <- logs(lab, tag)
    txt <- if (length(lg) == 1) readLines(lg, warn = FALSE) else character()
    done <- any(grepl("STEP 7 completed successfully", txt, fixed = TRUE))
    n_uf <- sum(grepl("weight underflow", txt, fixed = TRUE))
    dist <- suppressWarnings(as.numeric(sub(".*Mean redistribution distance: *([0-9.]+) km.*", "\\1",
                                            grep("Mean redistribution distance", txt, value = TRUE)[1])))
    V <- if (done) read_inj(lab, tag) else NULL
    s1 <- !is.null(V) && abs(sum(V$m) / Mref - 1) <= TOL_MASS
    s2 <- done && n_uf == 0
    line <- data.table(label = lab, variante = tag, S1 = s1, S2 = s2, underflow = n_uf, dist_km = dist)
    if (s1 && s2) {
      X <- merge(R[, .(i, j, k, mr = m)], V[, .(i, j, k, mv = m)], by = c("i", "j", "k"), all = TRUE)
      X[is.na(mr), mr := 0][is.na(mv), mv := 0]
      C <- X[, .(mr = sum(mr), mv = sum(mv)), by = .(i, j)]
      colV <- V[, .(m = sum(m)), by = .(i, j)][order(-m)]
      line[, `:=`(TV = 0.5 * X[, sum(abs(mv - mr))] / Mref,
                   TV_colonne = 0.5 * C[, sum(abs(mv - mr))] / Mref,
                   D_m = V[m > 0, sum(m * z) / sum(m)],
                   part_1re_col = colV$m[1] / sum(V$m), part_5_col = sum(head(colV$m, 5)) / sum(V$m),
                   active_cells = V[m > 0, .N])]
      line[, borne_pt := 100 * TV]
      line[, class_by_bound := fifelse(borne_pt <= 0.10 * G_AF, "negligible (established by the bound)", "to be measured, or bound only")]
    }
    res[[length(res) + 1]] <- line
    cat(sprintf("  %-9s S1 %-5s S2 %-5s (underflow %d)%s\n", tag, s1, s2, n_uf,
                if (s1 && s2) sprintf("  TV %.4f | column TV %.4f | bound %.2f pt | D %.2f m | 1st col %.3f | 5 col %.3f | cells %d | dist %.1f km",
                                      line$TV, line$TV_colonne, line$borne_pt, line$D_m, line$part_1re_col,
                                      line$part_5_col, line$active_cells, dist) else "  -> NULL variant, not measured"))
  }
}
if (stop) { cat("\nA LOCK S0 OR S3 FAILED: no variant is read, choix_ocim.txt is not written.\n"); quit(status = 1) }
T <- rbindlist(res, fill = TRUE)
fwrite(T, file.path(BASE, "etape_S.csv"))

cat("\n==================== OCIM choice (rule written beforehand) ====================\n")
for (lab in LABS) {
  l5 <- T[label == lab & variante == "l500_k16"]
  if (isTRUE(l5$S1 && l5$S2)) choice <- c(choice, paste(lab, "l500_k16")) else cat(sprintf("  %s l500_k16: NULL (S1/S2), no OCIM\n", lab))
  kk <- T[label == lab & variante %in% c("l100_k1", "l100_k4", "l100_k64") & S1 & S2][order(-TV)]
  if (nrow(kk) > 0 && kk$TV[1] > THRESHOLD_TV) choice <- c(choice, paste(lab, kk$variante[1]))
  for (r in seq_len(nrow(kk))) if (!(r == 1 && kk$TV[1] > THRESHOLD_TV))
    cat(sprintf("  %s %s: by the bound only, |dAF| <= %.2f pt (%s)\n", lab, kk$variante[r], kk$borne_pt[r], kk$class_by_bound[r]))
}
writeLines(choice, file.path(BASE, "choix_ocim.txt"))
cat("  passed to OCIM:\n"); cat(paste0("    ", choice, "\n"), sep = "")
cat("  + the displacement (Label A), if d3_move.py wrote MOVE_PASS\n")
cat("END step S\n")
