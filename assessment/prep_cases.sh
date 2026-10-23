#!/bin/bash
# Prepares the inputs of tests R and Omega.
# Writes ONLY in $T (outside pipeline: the manifest is not touched).
#  - test YAML: copy of the YAML in service where ONLY the omega_cm_yr line of `default`
#    changes (values of the scenarios omega_0p1 / omega_1 / omega_5 of the same YAML). Guard: 1 line of diff.
#  - frac tables at f = 0.9 / 0.8 / 0.7: frac_f = max(V_balaye - f*V_retire, 0) / V_balaye, columns of the
#    original CSV. Lock R0: at f = 1 the formula gives back the frac column to 1e-9.
set -euo pipefail
T="${DFP_REPO:?set DFP_REPO}/assessment"
Y0="${DFP_REPO:?set DFP_REPO}/config/fi_parameters_with_freshness.yaml"
D="${DFP_DATA:?set DFP_DATA}/output_V6_rejeu/decharge"
mkdir -p "$T/yaml" "$T/frac"

[ "$(md5sum "$Y0" | cut -c1-8)" = "8619588d" ] || { echo "ABORT: YAML in service is not 8619588d"; exit 1; }
L=$(grep -n '^    omega_cm_yr: 0.42 ' "$Y0" | cut -d: -f1)
[ "$(echo "$L" | wc -w)" -eq 1 ] || { echo "ABORT: omega line of default not found or multiple"; exit 1; }
for v in 0.1 1.0 5.0; do
  out="$T/yaml/omega_${v}.yaml"
  sed "${L}s/omega_cm_yr: 0.42 /omega_cm_yr: ${v} /" "$Y0" > "$out"
  n=$(diff "$Y0" "$out" | grep -c '^[<>]' || true)
  [ "$n" -eq 2 ] || { echo "ABORT: $out differs by $n lines (1 changed line expected)"; exit 1; }
  echo "YAML omega=$v: line $L -> $(sed -n "${L}p" "$out" | cut -c1-40)"
done

export R_LIBS_USER="$HOME/R/library"
Rscript - "$D" "$T/frac" <<'REOF'
suppressMessages(library(data.table))
a <- commandArgs(TRUE); D <- a[1]; OUT <- a[2]
for (lab in c("smooth", "flag")) {
  x <- fread(file.path(D, sprintf("fraction_sur_place_%s.csv", lab)))
  f1 <- pmax(x$V_balaye - 1.0 * x$V_retire, 0) / x$V_balaye
  e <- max(abs(f1 - x$frac), na.rm = TRUE)
  cat(sprintf("%s: lock R0, max gap frac(f=1) reconstructed - CSV = %.3g (%d rows, %d NA)\n",
              lab, e, nrow(x), sum(is.na(f1) != is.na(x$frac))))
  if (!(e < 1e-9) || any(is.na(f1) != is.na(x$frac))) stop("LOCK R0 FAILED: wrong reconstruction")
  for (f in c(0.9, 0.8, 0.7)) {
    y <- copy(x)
    y[, frac := pmax(V_balaye - f * V_retire, 0) / V_balaye]
    y[, V_sur_place := pmax(V_balaye - f * V_retire, 0)]
    fwrite(y, file.path(OUT, sprintf("fraction_sur_place_%s_f%s.csv", lab, f)))
    cat(sprintf("  f = %.1f: mean frac %.4f (f=1: %.4f)\n", f, mean(y$frac, na.rm=TRUE), mean(x$frac, na.rm=TRUE)))
  }
}
REOF
ls -la "$T/yaml" "$T/frac"
