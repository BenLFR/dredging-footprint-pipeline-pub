# -*- coding: utf-8 -*-
"""
lock_bottom.py -- test QF (reading rule written beforehand): locks of the INJECTION_MODE=bottom forcing
built on the capped source, against the installed forcing step7_<lab>_an (INJECTION_MODE=seafloor).
  B1  same mass: |q_bottom / q_annual - 1| <= 1e-9 (flux of Jdredge, 12 g/mol)
  B2  same horizontal distribution: for each column (i, j), |m_bottom - m_annual| <= 1e-9 x total
  B3  depth contrast: D_bottom > D_annual (mass-weighted depth of the layers); otherwise QF NULL
Geometry read from the OCIM 1x run under the cap (ocim_<lab>_an), as in d3_move.py.
Writes BASE/BOTTOM_<lab>_PASS if all three pass; exits with an error otherwise (the afterok OCIM does not start).
Usage: python3 lock_bottom.py <SVR> <BASE> <lab>
"""
import glob, os, sys
import numpy as np, h5py
from scipy.io import loadmat

SVR, BASE, LAB = sys.argv[1], sys.argv[2], sys.argv[3]
RHO = 1025.0                    # kg/m3: ns_step_eb.m, as in d3_move.py
TOL_MASS, TOL_COL = 1e-9, 1e-9  # B1: threshold of the Q0b test; B2: design choice, same order


def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]


O1 = only_one("%s/ocim_%s_an/ocim_dredge_1x_*.mat" % (SVR, LAB))
FA = only_one("%s/step7_%s_an/jdredge_ocim2_48l_*.mat" % (SVR, LAB))
FB = only_one("%s/step7_%s_bot/jdredge_ocim2_48l_*.mat" % (BASE, LAB))
print("geometry:", O1)
print("installed forcing (seafloor):", FA)
print("bottom forcing:", FB)

with h5py.File(O1, "r") as f:
    M3d = np.array(f["M3d"]).T
    g = f["grid"]
    V3 = np.array(g["DXT3d"]).T * np.array(g["DYT3d"]).T * np.array(g["DZT3d"]).T
    Z3 = np.array(g["ZT3d"]).T
shp = M3d.shape
wet = np.flatnonzero(M3d.ravel(order="F") == 1)
m = wet.size
Mw = V3.ravel(order="F")[wet] * RHO
Z = Z3.ravel(order="F")[wet]
ii, jj, kk = np.unravel_index(wet, shp, order="F")
col = ii * shp[1] + jj

JA = np.asarray(loadmat(FA)["Jdredge"], dtype=float).ravel()
JB = np.asarray(loadmat(FB)["Jdredge"], dtype=float).ravel()
if JA.size != m or JB.size != m:
    sys.exit("ABORT: lengths %d and %d != m %d" % (JA.size, JB.size, m))

mass = lambda x: float(np.sum(x * Mw) * 12.0 * 1e-6)   # gC/yr


def mean_depth(x):
    w = x * Mw
    return float(np.sum(Z * w) / np.sum(w))


qa, qb = mass(JA), mass(JB)
b1 = abs(qb / qa - 1) <= TOL_MASS
ma = np.bincount(col, weights=JA * Mw, minlength=shp[0] * shp[1])
mb = np.bincount(col, weights=JB * Mw, minlength=shp[0] * shp[1])
tot = float(ma.sum())
e2 = float(np.max(np.abs(mb - ma))) / tot
b2 = e2 <= TOL_COL
Da, Db = mean_depth(JA), mean_depth(JB)
b3 = Db > Da

print("\n===== %s =====" % LAB)
print("  B1 mass: q_annual %.6e g/yr ; q_bottom %.6e g/yr ; gap %.2e  %s" % (qa, qb, qb / qa - 1, "PASS" if b1 else "FAIL"))
print("  B2 columns: %d (yr) and %d (bottom) ; max |m_bottom - m_annual| / total = %.2e  %s"
      % (int(np.count_nonzero(ma)), int(np.count_nonzero(mb)), e2, "PASS" if b2 else "FAIL"))
print("  B3 mean depth: D_annual %.2f m ; D_bottom %.2f m  %s" % (Da, Db, "PASS" if b3 else "FAIL (QF NULL)"))
print("  (outside verdict) active cells (i, j, k): %d (yr) and %d (bottom)" % (np.count_nonzero(JA), np.count_nonzero(JB)))
if not (b1 and b2 and b3):
    sys.exit("\n>>> A B LOCK HAS FAILED: the bottom OCIM of %s does not start." % LAB)
open(os.path.join(BASE, "BOTTOM_%s_PASS" % LAB), "w").write(FB + "\n")
print(">>> BOTTOM FORCING %s READY (BOTTOM_%s_PASS written)." % (LAB, LAB))
