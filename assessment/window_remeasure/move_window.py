# -*- coding: utf-8 -*-
"""
move_window.py -- test DF (reading rule written beforehand): move by one grid cell the column that receives
the most mass (Label A), on the capped source (step7_smooth_an). Derived from d3_move.py:
same operator, same choice of neighbour, same locks C1-C6. Only the input and the R0 locks change:
  R0a  flux of 'Jdredge' (step7_smooth_an) = flux of 'Jdredge_forcing' of the 1x run under the cap, to 1e-9 relative
  R0d  the dominant column carries a mass share of 31.7 % determined to within 0.05 point (manuscript: "the largest cell held
       31.7%"); tolerance: half the printed precision
  (Unlike d3_move.py, the number of active cells and the mean depth are printed here, not locked.)
Operator: for each moved layer k, for each vector field of length m:
    y[nB] += y[nA] * Mw[nA] / Mw[nB] ; y[nA] = 0
Choice of neighbour: 4 neighbours (i+-1,j), (i,j+-1); eligible if wet at all layers k of the column;
    preference to untouched cells; then summed water mass closest to that of the origin; on a tie: north, east,
    south, west. None eligible: test NULL.
  C1  mass conserved per field to 1e-12 relative
  C2  active cells per field: n_after = n_before - (layers moved onto an already forced cell)
  C3  origin emptied, target filled (main field)
  C4  mean depth per field unchanged to 1e-6 m
  C5  ratios between fields identical to 1e-12
  C6  all targets are wet
Output: BASE/step7_smooth_an_shift/jdredge_ocim2_48l_anshift<dir>_20261010.mat and BASE/MOVE_PASS.
Usage: python3 move_window.py <SVR> <BASE>
"""
import glob, os, sys
import numpy as np, h5py
from scipy.io import loadmat, savemat

SVR, BASE = sys.argv[1], sys.argv[2]
RHO = 1025.0                       # kg/m3, as in ns_step_eb.m
PART_A, TOL_PART = 31.7, 0.05      # share of the largest cell, Label A (manuscript)


def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]


S7 = only_one(SVR + "/step7_smooth_an/jdredge_ocim2_48l_*.mat")
O1 = only_one(SVR + "/ocim_smooth_an/ocim_dredge_1x_*.mat")
print("capped forcing:", S7)
print("OCIM 1x run under the cap:", O1)

with h5py.File(O1, "r") as f:
    M3d = np.array(f["M3d"]).T
    g = f["grid"]
    V3 = np.array(g["DXT3d"]).T * np.array(g["DYT3d"]).T * np.array(g["DZT3d"]).T
    Z3 = np.array(g["ZT3d"]).T
    Jf = np.asarray(f["Jdredge_forcing"], dtype=float).ravel()
shp = M3d.shape
wet = np.flatnonzero(M3d.ravel(order="F") == 1)
m = wet.size
Mw = V3.ravel(order="F")[wet] * RHO
Z = Z3.ravel(order="F")[wet]
ii, jj, kk = np.unravel_index(wet, shp, order="F")
lut = {(int(a), int(b), int(c)): n for n, (a, b, c) in enumerate(zip(ii, jj, kk))}
print("shape M3d", shp, "| wet points", m)

mass = lambda x: float(np.sum(x * Mw) * 12.0 * 1e-6)          # gC/yr


def mean_depth(x):
    w = x * Mw
    return float(np.sum(Z * w) / np.sum(w)) if np.sum(w) != 0 else float("nan")


d = loadmat(S7)
var_keys = sorted(k for k in d.keys() if not k.startswith("__"))
vect = [k for k in var_keys if np.asarray(d[k]).size == m]
J = np.asarray(d["Jdredge"], dtype=float).ravel()

# ---- R0 ----------------------------------------------------------------------
print("\n===== R0: capped forcing =====")
if Jf.size != m:
    sys.exit("ABORT: Jdredge_forcing of length %d != m %d" % (Jf.size, m))
r0a = abs(mass(J) / mass(Jf) - 1) <= 1e-9
act = np.flatnonzero(J != 0)
mc = J * Mw
col = {}
for n in act:
    col[(int(ii[n]), int(jj[n]))] = col.get((int(ii[n]), int(jj[n])), 0.0) + mc[n]
tot = sum(col.values())
top = sorted(col.items(), key=lambda t: -t[1])
(i0, j0), m0 = top[0]
part0 = 100.0 * m0 / tot
r0d = abs(part0 - PART_A) <= TOL_PART
print("  R0a flux step7 / flux OCIM - 1 = %.2e  %s" % (mass(J) / mass(Jf) - 1, "PASS" if r0a else "FAIL"))
print("  R0d 1st column (%d,%d) in R indices, share %.3f %% (expected %.1f %% within %.2f)  %s"
      % (i0 + 1, j0 + 1, part0, PART_A, TOL_PART, "PASS" if r0d else "FAIL"))
print("  (outside verdict) active cells %d; columns %d; mean depth %.3f m" % (act.size, len(col), mean_depth(J)))
print("  first five columns (R indices, share %):",
      ", ".join("(%d,%d) %.2f" % (a + 1, b + 1, 100 * v / tot) for (a, b), v in top[:5]))
if not (r0a and r0d):
    sys.exit("\n>>> AN R0 LOCK HAS FAILED: nothing is built.")

# ---- survey ----------------------------------------------------------------------
kcol = sorted(int(kk[n]) for n in act if (int(ii[n]), int(jj[n])) == (i0, j0))
print("\n===== Dominant column (%d,%d) [R indices]: layers k (R indices) %s =====" % (i0 + 1, j0 + 1, [k + 1 for k in kcol]))
for k in kcol:
    n = lut[(i0, j0, k)]
    print("  k=%2d  z=%8.2f m  mass %.4e gC/yr  column share %.4f" % (k + 1, Z[n], mc[n], mc[n] / m0))
mA = sum(Mw[lut[(i0, j0, k)]] for k in kcol)
forced = set(act.tolist())
cand = []
print("\n  neighbour        wet at all these layers   untouched   water mass / origin")
for item_name, (di, dj) in (("north (j+1)", (0, 1)), ("east (i+1)", (1, 0)), ("south (j-1)", (0, -1)), ("west (i-1)", (-1, 0))):
    keys = [(i0 + di, j0 + dj, k) for k in kcol]
    is_wet = all(k in lut for k in keys)
    untouched = is_wet and all(lut[k] not in forced for k in keys)
    ratio = sum(Mw[lut[k]] for k in keys) / mA if is_wet else float("nan")
    print("  %-14s  %-28s  %-7s  %.6f" % (item_name, is_wet, untouched, ratio))
    if is_wet:
        cand.append((not untouched, abs(ratio - 1.0), len(cand), item_name, (di, dj)))
if not cand:
    open(os.path.join(BASE, "MOVE_NULL"), "w").write("no neighbour wet at all layers\n")
    sys.exit("\n>>> NO ELIGIBLE NEIGHBOUR: test DF NULL by construction. Nothing is written.")
cand.sort()
_, _, _, NAME, (DI, DJ) = cand[0]
print("\n  neighbour kept: %s -> column (%d,%d) in R indices" % (NAME, i0 + DI + 1, j0 + DJ + 1))
pairs = [(lut[(i0, j0, k)], lut[(i0 + DI, j0 + DJ, k)]) for k in kcol]


def move(x):
    y = x.copy()
    for nA, nB in pairs:
        y[nB] = y[nB] + y[nA] * Mw[nA] / Mw[nB]
        y[nA] = 0.0
    return y


# ---- build and locks ----------------------------------------------------------
out, bef, aft, ok = {}, {}, {}, True
print("\n===== C1 / C2 / C4, field by field =====")
for k in var_keys:
    v = np.asarray(d[k])
    if k not in vect:
        out[k] = v
        continue
    x = v.astype(float).ravel()
    y = move(x)
    out[k] = y.reshape(v.shape)
    m0_, m1_ = mass(x), mass(y)
    e = (m1_ / m0_ - 1) if m0_ != 0 else 0.0
    overlap = sum(1 for nA, nB in pairs if x[nA] != 0 and x[nB] != 0)
    a0, a1 = int(np.count_nonzero(x)), int(np.count_nonzero(y))
    p0, p1 = mean_depth(x), mean_depth(y)
    good = abs(e) <= 1e-12 and a1 == a0 - overlap and (np.isnan(p0) or abs(p1 - p0) <= 1e-6)
    ok &= good
    bef[k], aft[k] = m0_, m1_
    print("  %-28s mass %+.2e  active %d -> %d (overlaps %d)  depth %.4f -> %.4f  %s"
          % (k, e, a0, a1, overlap, p0, p1, "ok" if good else "!! FAILED"))
y0 = out["Jdredge"].ravel()
c3 = all(y0[nA] == 0.0 for nA, _ in pairs) and all(y0[nB] != 0.0 for _, nB in pairs)
c6 = all(0 <= nB < m for _, nB in pairs)
print("\n  C3 origin emptied / target filled: %s" % ("PASS" if c3 else "FAIL"))
print("  C6 wet targets: %s" % ("PASS" if c6 else "FAIL"))
ok &= c3 and c6
c5 = True
for k in vect:
    if bef["Jdredge"] == 0:
        continue
    ra, rb = bef[k] / bef["Jdredge"], aft[k] / aft["Jdredge"]
    c5 &= abs(rb - ra) <= 1e-12
print("  C5 ratios between fields: %s" % ("PASS" if c5 else "FAIL"))
ok &= c5
if not ok:
    sys.exit("\n>>> A C LOCK HAS FAILED: no file written.")

direction = {(0, 1): "N", (1, 0): "E", (0, -1): "S", (-1, 0): "O"}[(DI, DJ)]
OUTD = os.path.join(BASE, "step7_smooth_an_shift")
os.makedirs(OUTD, exist_ok=True)
if glob.glob(OUTD + "/*.mat"):
    sys.exit("ABORT: %s already contains a .mat" % OUTD)
OUT = os.path.join(OUTD, "jdredge_ocim2_48l_anshift%s_20261010.mat" % direction)
savemat(OUT, out, do_compression=True)
r = loadmat(OUT)
J1 = np.asarray(r["Jdredge"], dtype=float).ravel()
print("\n  written: %s (%.2f MB)" % (OUT, os.path.getsize(OUT) / 1e6))
print("  re-read: %d cells, flux / cap - 1 = %.2e, depth %.4f m" % (np.count_nonzero(J1), mass(J1) / mass(J) - 1, mean_depth(J1)))
open(os.path.join(BASE, "MOVE_PASS"), "w").write(OUT + "\n")
print(">>> MOVED FORCING READY (MOVE_PASS written).")
