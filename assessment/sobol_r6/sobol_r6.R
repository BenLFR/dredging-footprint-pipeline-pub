# sobol_r6.R -- Sobol analysis on the model in service (_r6_flux10).
# Replaces an earlier Sobol analysis on the real grid (not published), which varied alpha_dep,
# a parameter removed since Eq. 7 uses averages, and whose response was sum(f_i), not C_ri.
# Parameters and ranges, U(+/-50 %) around the
# central value -- range = modelling choice, no source:
#   fast_fraction 0.30 (modelling choice); slow_k 0.05 yr-1 (choice, bound < 0.1 yr-1 of
#   Hiddink 2023); k_mult 1 (multiplies the k_i of Sala 2021); p_r 0.87 (p_crd of Sala 2021, Methods),
#   upper bound brought back to 1; D_b 24.3 cm2 yr-1 and omega 0.42 cm yr-1 (Middelburg et al. 1996, table 1,
#   at 50 m, read at the source).
# Response: total C_ri = sum_i c_i f_i(params), c_i = C_ri / f_i_full read from the cri_pr output (constant
#   in the parameters: C_ri is linear in f_i).
# Estimator: sensitivity::sobolSalt, scheme A (Saltelli 2002), sobolEff estimators (Janon et al. 2014),
#   N = 1000, nboot = 200, 95 % CI. Seed 20260925. Same settings as that earlier analysis.
#
# RULE WRITTEN BEFORE THE RUN
#  Locks (if they fail, stop; no index is read):
#   V1. f_i rebuilt at the central values = f_i_full of _r6, max relative gap < 1e-9 on all
#       cells where f_i_full > 0, and no cell zero on one side only.
#   V2. C_ri_conservative / f_i_conservative = C_ri / f_i_full (max relative gap < 1e-9): C_ri is indeed
#       a fixed multiple of f_i per cell; and sum(c_i f_i rebuilt) = sum(C_ri) to within 1e-9.
#  Reading (reported as is, no reinterpretation):
#   H1. The hierarchy written in the manuscript HOLDS if the two largest ST are fast_fraction and p_r and
#       their sum is >= 0.75; it FAILS otherwise (the parameter that comes between is then named).
#   H2. The depth profile is a first-order uncertainty if ST(D_b) or ST(omega) > 0.10;
#       otherwise it is secondary. (Says nothing of omega outside +/-50 %: the 0.1-5 cm yr-1 range is a
#       separate scenario.)
#   Reported without verdict: S1, CI; sum S1; conservative / central; stock bounds / central.
# Usage: Rscript sobol_r6.R <smooth|flag>
suppressMessages({ library(data.table); library(arrow); library(sensitivity) })
lab <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(lab %in% c("smooth", "flag"))
R <- file.path(Sys.getenv("DFP_SCRATCH"), "output_V6_sala2021_r6_flux10")
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
use <- d[!is.na(SVR) & !is.na(p_l) & !is.na(k_used) & !is.na(f_i_full)]
cat("usable cells:", nrow(use), "\n")

# depth weights: copy of step5_merge_tiles_corr.R, D_b and omega as arguments
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

# ---------- V1 ----------
f0 <- do.call(f_cells, as.list(unname(central)))
p0 <- use$f_i_full > 0; p1 <- f0 > 0
e1 <- max(abs(f0[p0 & p1] / use$f_i_full[p0 & p1] - 1))
x1 <- sum(xor(p0, p1))
cat(sprintf("V1: max relative gap %.3g on %d cells ; zero on one side only: %d\n", e1, sum(p0 & p1), x1))
if (!(e1 < 1e-9 && x1 == 0)) { cat("LOCK V1 FAILED: no index is read.\n"); quit(status = 1) }
# ---------- V2 ----------
use[, c_i := fifelse(f_i_full > 0, C_ri / f_i_full, 0)]
b <- use$f_i_full > 0 & use$f_i_conservative > 0 & !is.na(use$C_ri_conservative)
# Where c_i = 0 (zero stock), the relative gap is 0/0. The rule is then applied in its absolute
# form: c_i = 0 requires C_ri_conservative = 0. The 1e-9 threshold and the rest of the rule are unchanged.
# (Form fixed after a NaN stop at the first launch, before any reading of V2.)
bp <- b & use$c_i > 0; b0 <- b & use$c_i == 0
e2 <- max(abs((use$C_ri_conservative[bp] / use$f_i_conservative[bp]) / use$c_i[bp] - 1))
z2 <- sum(use$C_ri_conservative[b0] != 0)
s0 <- sum(use$c_i * f0); sR <- sum(d$C_ri, na.rm = TRUE)
cat(sprintf("V2: max relative gap %.3g on %d cells (c_i > 0) ; c_i = 0: %d cells, of which %d with C_ri_cons != 0 ; sum c_i f0 = %.6f, sum C_ri = %.6f (gap %.3g)\n",
            e2, sum(bp), sum(b0), z2, s0, sR, abs(s0 / sR - 1)))
if (!(e2 < 1e-9 && z2 == 0 && abs(s0 / sR - 1) < 1e-9)) { cat("LOCK V2 FAILED: no index is read.\n"); quit(status = 1) }
cat("LOCKS PASSED.\n\n")

# ---------- Sobol ----------
lo <- c(ff = 0.15, ks = 0.025, km = 0.5, pr = 0.435, Db = 12.15, om = 0.21)
hi <- c(ff = 0.45, ks = 0.075, km = 1.5, pr = 1.00,  Db = 36.45, om = 0.63)
set.seed(20260925); N <- 1000L
mk <- function() as.data.frame(mapply(function(a, b) runif(N, a, b), lo, hi))
X1 <- mk(); X2 <- mk(); names(X1) <- names(X2) <- names(central)
model_fn <- function(X) apply(X, 1, function(p) sum(use$c_i * f_cells(p[1], p[2], p[3], p[4], p[5], p[6])))
t0 <- Sys.time()
sa <- sobolSalt(model = model_fn, X1 = X1, X2 = X2, scheme = "A", conf = 0.95, nboot = 200)
cat("Sobol:", nrow(sa$X), "evaluations,", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
nm <- c(ff = "fast_fraction", ks = "slow_k", km = "k_mult", pr = "p_r", Db = "D_b", om = "omega")
out <- data.table(parameter = nm[names(central)],
                  S1 = sa$S[, 1], S1_lo = sa$S[, 4], S1_hi = sa$S[, 5],
                  ST = sa$T[, 1], ST_lo = sa$T[, 4], ST_hi = sa$T[, 5])
setorder(out, -ST); print(out, digits = 4)
cat(sprintf("sum S1 = %.4f ; sum ST = %.4f\n", sum(out$S1), sum(out$ST)))
top2 <- out$parameter[1:2]; s2 <- sum(out$ST[1:2])
cat(sprintf("H1: top two = %s ; sum ST = %.4f -> %s\n", paste(top2, collapse = " + "), s2,
            if (setequal(top2, c("fast_fraction", "p_r")) && s2 >= 0.75) "HOLDS" else "FAIL"))
stDb <- out[parameter == "D_b", ST]; stom <- out[parameter == "omega", ST]
cat(sprintf("H2: ST(D_b) = %.4f ; ST(omega) = %.4f -> %s\n", stDb, stom,
            if (max(stDb, stom) > 0.10) "FIRST ORDER" else "SECONDARY"))
cat(sprintf("\n[reported] C_ri central %.4f ; conservative %.4f (x%.4f, %+.1f %%) ; low stock %.4f (%+.2f %%) ; high %.4f (%+.2f %%)\n",
            sR, sum(d$C_ri_conservative, na.rm = TRUE), sum(d$C_ri_conservative, na.rm = TRUE) / sR,
            100 * (sum(d$C_ri_conservative, na.rm = TRUE) / sR - 1),
            sum(d$C_ri_lower, na.rm = TRUE), 100 * (sum(d$C_ri_lower, na.rm = TRUE) / sR - 1),
            sum(d$C_ri_upper, na.rm = TRUE), 100 * (sum(d$C_ri_upper, na.rm = TRUE) / sR - 1)))
od <- file.path(R, "sensitivity"); dir.create(od, showWarnings = FALSE)
fwrite(out, file.path(od, paste0("sobol_r6_", lab, ".csv")))
cat("written:", file.path(od, paste0("sobol_r6_", lab, ".csv")), "\n")
