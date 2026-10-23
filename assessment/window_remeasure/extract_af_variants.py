# -*- coding: utf-8 -*-
"""
extract_af_variants.py -- annual series of the airborne fraction (AF), 2015-2100, of the three variants
under the window: injection at the bottom of the column (bottom_smooth, bottom_flag) and dominant column moved
by one step (shift_smooth). Extraction, not a test: nothing is rerun.

Copy of extract_af_annual.py: same constants, same formula, same lock; only
the folders read and the expected values change. Each variant has its own control and its own 1x run in
BASE/<subfolder> (as in read_remeasure.py).

LOCK, written before reading: the series must reproduce, to within 1e-4 point, the values printed by
read_remeasure.py:
  bottom_smooth: 2024 36.7873 ; 2033 44.9957 ; 2050 41.9086 ; 2100 38.9440 ; max 45.0721 in 2032
  bottom_flag   : 2024 34.1605 ; 2033 42.1476 ; 2050 39.3519 ; 2100 36.7088 ; max 42.2051 in 2032
  shift_smooth  : 2024 75.3481 ; 2033 74.5720 ; 2050 68.4275 ; 2100 63.2397 ; max 78.2571 in 2026
If it fails for a variant, the series of that variant is not written, and neither is the CSV.
Execution: sbatch (python3 + h5py), never on the interactive node.
Output: stdout (the log) and af_variantes_20261016.csv in the folder passed as argument.
"""
import os, glob, sys
import numpy as np, h5py

RHO, MATM, GPMOL = 1025.0, 1.8e20, 12.0     # same values as read_remeasure.py (ns_step_eb.m)
BASE = os.environ["DFP_SCRATCH"] + "/remesure_fen_20261010"
VARIANTS = (("bottom_smooth", "ocim_smooth_bot"), ("bottom_flag", "ocim_flag_bot"),
             ("shift_smooth", "ocim_smooth_shift"))          # as in read_remeasure.py
EXPECTED = {"bottom_smooth": ({2024: 36.7873, 2033: 44.9957, 2050: 41.9086, 2100: 38.9440}, 45.0721, 2032),
           "bottom_flag": ({2024: 34.1605, 2033: 42.1476, 2050: 39.3519, 2100: 36.7088}, 42.2051, 2032),
           "shift_smooth": ({2024: 75.3481, 2033: 74.5720, 2050: 68.4275, 2100: 63.2397}, 78.2571, 2026)}
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


def series(sd):
    d = "%s/%s" % (BASE, sd)
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


rows, all_ok = {}, True
for item_name, sd in VARIANTS:
    yr, af, cum, dp, i0 = series(sd)
    steps = np.diff(yr[i0:])
    print("==== %s (%s): %d steps from %d to %d ; annual step throughout: %s" % (item_name, sd, yr[i0:].size, yr[i0], yr[-1],
          bool(np.all(steps == 1))))
    milestones, amax, ymax = EXPECTED[item_name]
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
        print("  >>> LOCK FAILED: series %s not written" % item_name)
        all_ok = False
        continue
    for k in range(i0, yr.size):
        rows.setdefault(int(yr[k]), {})[item_name] = af[k]

nms = [n for n, _ in VARIANTS]
if all_ok and len(rows) and all(len(v) == len(nms) for v in rows.values()):
    out_path = sys.argv[1] + "/af_variantes_20261016.csv"
    with open(out_path, "w") as f:
        f.write("annee,af_A_bottom_pct,af_B_bottom_pct,af_A_shift_pct\n")
        for y in sorted(rows):
            f.write("%d,%.4f,%.4f,%.4f\n" % ((y,) + tuple(rows[y][n] for n in nms)))
    print("WRITTEN: %s (%d rows)" % (out_path, len(rows)))
else:
    print("NOTHING WRITTEN")
