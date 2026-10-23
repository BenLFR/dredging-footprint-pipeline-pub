# -*- coding: utf-8 -*-
"""
extract_af_annual.py -- annual series of the airborne fraction (AF), 2015-2100, of the two reference runs
under the window (Label A = ocim_smooth_fen, Label B = ocim_flag_fen). Extraction, not a test: nothing is rerun.

End-of-year AF = atmospheric CO2 anomaly (as C mass) / cumulative injected carbon, in %; formula and
constants identical to read_remeasure.py (table `af`), 1x run (G4 passes for both).

LOCK, written before reading: the series must reproduce, to within 1e-4 point, the values printed by
read_remeasure.py:
  smooth: 2024 74.4196 ; 2033 74.6300 ; 2050 68.4704 ; 2100 63.2612 ; max 77.8578 in 2026
  flag   : 2024 62.9787 ; 2033 63.6180 ; 2050 59.1586 ; 2100 55.8185 ; max 65.8827 in 2027
If it fails for a label, the series of that label is not written.
Execution: sbatch (python3 + h5py), never on the interactive node.
Output: stdout (the log) and af_annuelle_20261013.csv in the folder passed as argument.
"""
import os, glob, sys
import numpy as np, h5py

RHO, MATM, GPMOL = 1025.0, 1.8e20, 12.0     # same values as read_remeasure.py (ns_step_eb.m)
WIN = os.environ["DFP_SCRATCH"] + "/ocim_fenetre_20261009"
EXPECTED = {"smooth": ({2024: 74.4196, 2033: 74.6300, 2050: 68.4704, 2100: 63.2612}, 77.8578, 2026),
           "flag": ({2024: 62.9787, 2033: 63.6180, 2050: 59.1586, 2100: 55.8185}, 65.8827, 2027)}
TOL = 1e-4


def T(a):
    return np.array(a).T


def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]


def vec(f, name):
    return np.asarray(f[name], dtype=float).ravel()


def volumes(f):
    M3d = T(f["M3d"]).ravel(order="F")
    g = f["grid"]
    V3 = (T(g["DXT3d"]) * T(g["DYT3d"]) * T(g["DZT3d"])).ravel(order="F")
    return V3[np.flatnonzero(M3d == 1)]


def series(lab):
    d = "%s/ocim_%s_fen" % (WIN, lab)
    fb = h5py.File(only_one(d + "/ocim_baseline_*.mat"), "r")
    fp = h5py.File(only_one(d + "/ocim_dredge_1x_*.mat"), "r")
    V = volumes(fp)
    J = np.asarray(fp["Jdredge_forcing"], dtype=float).ravel()
    rate = float(np.sum(J * V * RHO) * GPMOL * 1e-6)
    yr = vec(fp, "annee_pas").astype(int)
    act = vec(fp, "actif") != 0
    dp = vec(fp, "pco2a") - vec(fb, "pco2a")
    cum = np.cumsum(act) * rate
    i0 = int(np.flatnonzero(act)[0])
    af = np.full(dp.size, np.nan)
    af[i0:] = dp[i0:] * 1e-6 * MATM * GPMOL / cum[i0:] * 100.0
    return yr, af, cum, dp, i0


rows = {}
for lab in ("smooth", "flag"):
    yr, af, cum, dp, i0 = series(lab)
    steps = np.diff(yr[i0:])
    print("==== %s: %d steps from %d to %d ; annual step throughout: %s" % (lab, yr[i0:].size, yr[i0], yr[-1],
          bool(np.all(steps == 1))))
    milestones, amax, ymax = EXPECTED[lab]
    ok = bool(np.all(steps == 1))
    for y, v in milestones.items():
        k = int(np.flatnonzero(yr == y)[0])
        good = abs(af[k] - v) <= TOL
        ok &= good
        print("  end lock %d: %.4f against %.4f  %s" % (y, af[k], v, "PASS" if good else "FAIL"))
    kmax = int(np.nanargmax(af))
    good = abs(af[kmax] - amax) <= TOL and int(yr[kmax]) == ymax
    ok &= good
    print("  maximum lock: %.4f in %d against %.4f in %d  %s" % (af[kmax], yr[kmax], amax, ymax,
          "PASS" if good else "FAIL"))
    if not ok:
        print("  >>> LOCK FAILED: series %s not written" % lab)
        continue
    for k in range(i0, yr.size):
        rows.setdefault(int(yr[k]), {})[lab] = (af[k], cum[k] / 1e6, dp[k])

if len(rows) and all(len(v) == 2 for v in rows.values()):
    out_path = sys.argv[1] + "/af_annuelle_20261013.csv"
    with open(out_path, "w") as f:
        f.write("annee,af_A_smooth_pct,af_B_flag_pct,injecte_cumule_A_tC,injecte_cumule_B_tC\n")
        for y in sorted(rows):
            a, b = rows[y]["smooth"], rows[y]["flag"]
            f.write("%d,%.4f,%.4f,%.2f,%.2f\n" % (y, a[0], b[0], a[1], b[1]))
    print("WRITTEN: %s (%d rows)" % (out_path, len(rows)))
    for y in sorted(rows):
        print("  %d  A %.4f  B %.4f" % (y, rows[y]["smooth"][0], rows[y]["flag"][0]))
else:
    print("NOTHING WRITTEN")
