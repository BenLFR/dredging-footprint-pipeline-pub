# sobol_r6_cap.R -- Sobol of sobol_r6.R under the annual cap
# of the SVR. Only change: the response becomes sum_i c_i s_i f_i(params), s_i = sum_y min(SVR_y, 1) / SVR
# (fixed weight per cell, independent of the six parameters, read from s_annuel_<label>.csv).
# Parameters, ranges, seed, N, nboot, estimators: unchanged (see sobol_r6.R).
# Outputs in a new folder; nothing is written to _r6.
#
# RULE WRITTEN BEFORE THE RUN
#  Locks (if one fails, no index is read):
#   V1, V2: those of sobol_r6.R, unchanged.
#   V3. sum(c_i s_i f_i) at the central values = 3432.5263 / 4685.5941 (capped totals, smooth / flag) to within 1e-3; no cell with
#       C_ri > 0 without s; the s file has as many rows as cells with C_ri > 0.
#  Reading:
#   H1, H2: those of sobol_r6.R, unchanged.
#   H3. The cap changes the ranking if the order of the six ST differs from that of sobol_r6.R (p_r, fast_fraction, omega,
#       D_b, slow_k, k_mult) AND the 95 % CIs of the two swapped ST do not overlap; otherwise the ranking is
#       unchanged. Threshold: design choice, no source.
#   Reported without verdict: S1, CI; sum S1, sum ST; gap of each ST to sobol_r6.R; conservative / central and
#   stock bounds / central under the cap.
# Usage: Rscript sobol_r6_cap.R <smooth|flag> <s_annuel_<label>.csv>
suppressMessages({ library(data.table); library(arrow); library(sensitivity) })
a <- commandArgs(trailingOnly = TRUE); lab <- a[1]; fs <- a[2]
stopifnot(lab %in% c("smooth", "flag"), file.exists(fs))
CAN <- c(smooth = 3432.5263, flag = 4685.5941)
ST_172 <- list(smooth = c(p_r = 0.596, fast_fraction = 0.340, omega = 0.0197, D_b = 0.0195, slow_k = 0.0190, k_mult = 0.0154),
               flag   = c(p_r = 0.595, fast_fraction = 0.340, omega = 0.0197, D_b = 0.0197, slow_k = 0.0189, k_mult = 0.0160))
ORDER_172 <- c("p_r", "fast_fraction", "omega", "D_b", "slow_k", "k_mult")
R <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r6_flux10")
OD <- file.path(Sys.getenv("DFP_SCRATCH"), "sobol_plafond_20261005")
only_one <- function(d, pat) { f <- list.files(d, pat, full.names = TRUE)
  if (length(f) != 1) stop(d, ": ", length(f), " files ", pat); f }
fi <- as.data.table(read_parquet(only_one(file.path(R, paste0("merge_", lab, "_pr")), "^fi_grid_.*[.]parquet$")))
cr <- as.data.table(read_parquet(only_one(file.path(R, "cri_pr"), paste0("^cri_final_", lab, "_.*[.]parquet$"))))
cat("fi_grid:", nrow(fi), "cells ; cri:", nrow(cr), "\n")
d <- merge(fi[, .(grid_id, SVR, p_l, p_d, fresh_fact, k_used, frac_sur_place, f_i_full, f_i_conservative)],
           cr[, .(grid_id, C_ri, C_ri_conservative, C_ri_lower, C_ri_upper)], by = "grid_id")
cat("paired:", nrow(d), "\n")
d[is.na(fresh_fact), fresh_fact := 1]
d[is.na(frac_sur_place), frac_sur_place := 1]
d[, p_d_safe := fifelse(is.na(p_d) | p_d <= 0, 1, p_d)]
d[, `:=`(w1 = pmin(0.05, p_d_safe) / p_d_safe,
         w2 = pmin(0.05, pmax(p_d_safe - 0.05, 0)) / p_d_safe,
         w3 = pmax(p_d_safe - 0.10, 0) / p_d_safe,
         H3 = pmax(p_d_safe - 0.10, 0) * 100)]
sa <- fread(fs)
d <- merge(d, sa[, .(grid_id, s)], by = "grid_id", all.x = TRUE)
ns <- sum(d$C_ri > 0 & is.na(d$s), na.rm = TRUE); nr <- nrow(sa); npos <- sum(d$C_ri > 0, na.rm = TRUE)
d[is.na(s), s := 1]
use <- d[!is.na(SVR) & !is.na(p_l) & !is.na(k_used) & !is.na(f_i_full)]
cat("usable cells:", nrow(use), "\n")

# depth weights and f_i: copies of sobol_r6.R, unchanged
weight_depth <- function(k, Db, om, w1, w2, w3, H3) {
  h <- 10; rq <- sqrt(om^2 + 4 * Db * k)
  a <- 2 * k / (rq + om); b <- (rq + om) / (2 * Db)
  I <- function(x) -expm1(-a * x) / a + (a / b^2) * exp(-a * h) * (exp(b * (x - h)) - exp(-b * h))
  m0 <- I(5) / 5; al <- (I(10) - I(5)) / I(5); gh <- exp(-a * h) * (1 + a / b); L <- om / k
  wd <- ifelse(H3 > 0, gh * L * (-expm1(-H3 / L)) / (pmax(H3, 1e-12) * m0), 0)
  w1 + w2 * al + w3 * wd
}
f_cells <- function(ff, ks, km, pr, Db, om, x = use) {
  kf <- x$k_used * km
  v <- x$SVR * x$p_l * x$fresh_fact * x$frac_sur_place * pr *
    (ff * weight_depth(kf, Db, om, x$w1, x$w2, x$w3, x$H3) * (1 - exp(-kf)) +
     (1 - ff) * weight_depth(ks, Db, om, x$w1, x$w2, x$w3, x$H3) * (1 - exp(-ks)))
  pmin(pmax(v, 0), 1)
}
central <- c(ff = 0.30, ks = 0.05, km = 1, pr = 0.87, Db = 24.3, om = 0.42)

# ---------- V1 (unchanged) ----------
f0 <- do.call(f_cells, as.list(unname(central)))
p0 <- use$f_i_full > 0; p1 <- f0 > 0
e1 <- max(abs(f0[p0 & p1] / use$f_i_full[p0 & p1] - 1))
x1 <- sum(xor(p0, p1))
cat(sprintf("V1: max relative gap %.3g on %d cells ; zero on one side only: %d\n", e1, sum(p0 & p1), x1))
if (!(e1 < 1e-9 && x1 == 0)) { cat("LOCK V1 FAILED: no index is read.\n"); quit(status = 1) }
# ---------- V2 (unchanged, including the 0/0 correction) ----------
use[, c_i := fifelse(f_i_full > 0, C_ri / f_i_full, 0)]
b <- use$f_i_full > 0 & use$f_i_conservative > 0 & !is.na(use$C_ri_conservative)
bp <- b & use$c_i > 0; b0 <- b & use$c_i == 0
e2 <- max(abs((use$C_ri_conservative[bp] / use$f_i_conservative[bp]) / use$c_i[bp] - 1))
z2 <- sum(use$C_ri_conservative[b0] != 0)
s0 <- sum(use$c_i * f0); sR <- sum(d$C_ri, na.rm = TRUE)
cat(sprintf("V2: max relative gap %.3g on %d cells (c_i > 0) ; c_i = 0: %d cells, of which %d with C_ri_cons != 0 ; sum c_i f0 = %.6f, sum C_ri = %.6f (gap %.3g)\n",
            e2, sum(bp), sum(b0), z2, s0, sR, abs(s0 / sR - 1)))
if (!(e2 < 1e-9 && z2 == 0 && abs(s0 / sR - 1) < 1e-9)) { cat("LOCK V2 FAILED: no index is read.\n"); quit(status = 1) }
# ---------- V3 (new) ----------
sA <- sum(use$c_i * use$s * f0)
cat(sprintf("V3: sum c_i s_i f0 = %.4f (expected %.4f) ; cells C_ri > 0 without s: %d ; rows of the s file: %d, cells C_ri > 0: %d\n",
            sA, CAN[[lab]], ns, nr, npos))
if (!(abs(sA - CAN[[lab]]) <= 1e-3 && ns == 0 && nr == npos)) { cat("LOCK V3 FAILED: no index is read.\n"); quit(status = 1) }
cat("LOCKS PASSED.\n\n")

# ---------- Sobol (settings unchanged; only the response is weighted by s_i) ----------
lo <- c(ff = 0.15, ks = 0.025, km = 0.5, pr = 0.435, Db = 12.15, om = 0.21)
hi <- c(ff = 0.45, ks = 0.075, km = 1.5, pr = 1.00,  Db = 36.45, om = 0.63)
set.seed(20260925); N <- 1000L
mk <- function() as.data.frame(mapply(function(a, b) runif(N, a, b), lo, hi))
X1 <- mk(); X2 <- mk(); names(X1) <- names(X2) <- names(central)
cs <- use$c_i * use$s
model_fn <- function(X) apply(X, 1, function(p) sum(cs * f_cells(p[1], p[2], p[3], p[4], p[5], p[6])))
t0 <- Sys.time()
sa <- sobolSalt(model = model_fn, X1 = X1, X2 = X2, scheme = "A", conf = 0.95, nboot = 200)
cat("Sobol:", nrow(sa$X), "evaluations,", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
nm <- c(ff = "fast_fraction", ks = "slow_k", km = "k_mult", pr = "p_r", Db = "D_b", om = "omega")
out <- data.table(parameter = nm[names(central)],
                  S1 = sa$S[, 1], S1_lo = sa$S[, 4], S1_hi = sa$S[, 5],
                  ST = sa$T[, 1], ST_lo = sa$T[, 4], ST_hi = sa$T[, 5])
out[, ST_146_172 := ST_172[[lab]][parameter]]
out[, gap_ST := ST - ST_146_172]
setorder(out, -ST); print(out, digits = 4)
cat(sprintf("sum S1 = %.4f ; sum ST = %.4f\n", sum(out$S1), sum(out$ST)))
top2 <- out$parameter[1:2]; s2 <- sum(out$ST[1:2])
cat(sprintf("H1: top two = %s ; sum ST = %.4f -> %s\n", paste(top2, collapse = " + "), s2,
            if (setequal(top2, c("fast_fraction", "p_r")) && s2 >= 0.75) "HOLDS" else "FAIL"))
stDb <- out[parameter == "D_b", ST]; stom <- out[parameter == "omega", ST]
cat(sprintf("H2: ST(D_b) = %.4f ; ST(omega) = %.4f -> %s\n", stDb, stom,
            if (max(stDb, stom) > 0.10) "FIRST ORDER" else "SECONDARY"))
# H3: for each pair whose order is reversed relative to sobol_r6.R, do the CIs overlap?
rank <- setNames(seq_along(out$parameter), out$parameter)
inv <- list()
for (i in 1:5) for (j in (i + 1):6) {
  a1 <- ORDER_172[i]; a2 <- ORDER_172[j]                       # in sobol_r6.R: a1 before a2
  if (rank[[a2]] < rank[[a1]]) {
    lo1 <- out[parameter == a1, ST_lo]; hi1 <- out[parameter == a1, ST_hi]
    lo2 <- out[parameter == a2, ST_lo]; hi2 <- out[parameter == a2, ST_hi]
    inv[[length(inv) + 1]] <- sprintf("%s/%s IC %s", a1, a2, if (hi1 < lo2 || hi2 < lo1) "DISJOINT" else "overlapping")
  }
}
cat(sprintf("H3: inversions relative to sobol_r6.R: %s -> %s\n",
            if (length(inv)) paste(unlist(inv), collapse = " ; ") else "none",
            if (any(grepl("DISJOINT", unlist(inv)))) "RANKING CHANGED" else "RANKING UNCHANGED"))
cons <- sum(use$C_ri_conservative * use$s, na.rm = TRUE)
lw <- sum(use$C_ri_lower * use$s, na.rm = TRUE); up <- sum(use$C_ri_upper * use$s, na.rm = TRUE)
cat(sprintf("\n[reported, under the cap] central %.4f ; conservative %.4f (x%.4f, %+.1f %%) ; low stock %.4f (%+.2f %%) ; high %.4f (%+.2f %%)\n",
            sA, cons, cons / sA, 100 * (cons / sA - 1), lw, 100 * (lw / sA - 1), up, 100 * (up / sA - 1)))
dir.create(OD, showWarnings = FALSE)
fwrite(out, file.path(OD, paste0("sobol_cap_", lab, ".csv")))
cat("written:", file.path(OD, paste0("sobol_cap_", lab, ".csv")), "\n")
