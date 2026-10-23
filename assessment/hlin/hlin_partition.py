# -*- coding: utf-8 -*-
"""
hlin_partition.py -- new test of H-lin, set-up B (rule written beforehand).
Splits the smooth M4 forcing into two disjoint parts:
    F_A = the dominant column alone (all its layers),
    F_B = everything else,
so that F_A + F_B = F_ref, field by field, bit for bit. OCIM then runs on F_A and on F_B, and
read_hlin.py compares the sum of the responses with the response to F_ref.
Derived from d3_move.py: same reading of step7 and of the 1x OCIM run, same locks R0.
Operator: a mask (column / outside column), applied to EACH vector field of length m, with no rule by
field name. A mask splits any field exactly, additive or not.

ORDER IMPOSED (nothing is written if a lock fails):
  R0  reproduction of the forcing of test Q (identical to d3_move.py):
      R0a  flux of 'Jdredge' (step7) = flux of 'Jdredge_forcing' (OCIM 1x) to 1e-9 relative
      R0b  208 active cells
      R0c  mass-weighted mean depth = 26.0 m to within 0.1 m
      R0d  first column = (47,53) in R indices, share 32.6 % to within 0.05 point (test N)
  P1  for each vector field: F_A + F_B == F_ref, identically (equal array)
  P2  no non-zero cell in both F_A and F_B (main field Jdredge)
  P3  F_A carries only the active cells of the column; F_B carries 208 minus that number
  P4  share of F_A = 32.553 % to within 0.001 point (R0d of d3_move.py, printed to 3 decimals)
  Re-reading of the two .mat files written: P1-P3 redone on the re-read files.
Outputs: BASE/step7_smooth_colA/jdredge_ocim2_48l_colA_20260928.mat,
         BASE/step7_smooth_resteB/jdredge_ocim2_48l_resteB_20260928.mat, then BASE/PARTITION_PASS.
Usage: python3 hlin_partition.py <RUN> <BASE>
"""
import glob, os, sys
import numpy as np, h5py
from scipy.io import loadmat, savemat
from scipy.sparse import issparse

RUN, BASE = sys.argv[1], sys.argv[2]
RHO = 1025.0            # kg/m3, as in ns_step_eb.m
D_Q, TOL_D = 26.0, 0.1  # test Q, printed to 0.1 m; tolerance: design choice, as in d3_move.py
N_ACT_Q = 208           # test Q (summary_q.py)
COL_R, PART_N, TOL_PART = (47, 53), 32.6, 0.05     # test N, R indices, share printed to 0.1 point
PART_R0D, TOL_P4 = 32.553, 0.001                    # R0d measured by d3_move.py (logs/move.out), 3 decimals

def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]

for sd in ("step7_smooth_colA", "step7_smooth_resteB"):
    if glob.glob(os.path.join(BASE, sd, "*.mat")):
        sys.exit("ABORT: %s already contains a .mat" % os.path.join(BASE, sd))
if os.path.exists(os.path.join(BASE, "PARTITION_PASS")):
    sys.exit("ABORT: PARTITION_PASS already exists")

S7 = only_one(RUN + "/step7_smooth_m4/jdredge_ocim2_48l_*.mat")
O1 = only_one(RUN + "/ocim_smooth_m4/ocim_dredge_1x_*.mat")
print("M4 forcing:", S7)
print("run OCIM 1x:", O1)

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
print("shape M3d", shp, "| wet points", m)

mass = lambda x: float(np.sum(x * Mw) * 12.0 * 1e-6)          # gC/yr
def mean_depth(x):
    w = x * Mw
    return float(np.sum(Z * w) / np.sum(w)) if np.sum(w) != 0 else float("nan")

d = loadmat(S7)
var_keys = sorted(k for k in d.keys() if not k.startswith("__"))
vect = [k for k in var_keys if np.asarray(d[k]).size == m]
J = np.asarray(d["Jdredge"], dtype=float).ravel()
print("fields:", len(var_keys), "of which vectors (length m):", len(vect))
# Guard: a non-vector field is copied as is into F_A AND into F_B. If it carried forcing (sparse
# matrix, or array whose size is a multiple of m), it would be counted twice. We stop instead.
for k in var_keys:
    v = d[k]
    print("  %-28s %-10s shape %s" % (k, "vector" if k in vect else "copied", getattr(v, "shape", "?")))
    if issparse(v) or (k not in vect and np.asarray(v).size % m == 0 and np.asarray(v).size > 0):
        sys.exit("ABORT: field %s (sparse or of size multiple of m) would be copied into both parts" % k)

# ---- R0 (identical to d3_move.py) ------------------------------------------
print("\n===== R0: reproduction of the forcing of test Q =====")
if Jf.size != m:
    sys.exit("ABORT: Jdredge_forcing of length %d != m %d" % (Jf.size, m))
r0a = abs(mass(J) / mass(Jf) - 1) <= 1e-9
act = np.flatnonzero(J != 0)
r0b = act.size == N_ACT_Q
D0 = mean_depth(J)
r0c = abs(D0 - D_Q) <= TOL_D
mc = J * Mw
col = {}
for n in act:
    col[(int(ii[n]), int(jj[n]))] = col.get((int(ii[n]), int(jj[n])), 0.0) + mc[n]
tot = sum(col.values())
top = sorted(col.items(), key=lambda t: -t[1])
(i0, j0), m0 = top[0]
part0 = 100.0 * m0 / tot
r0d = (i0 + 1, j0 + 1) == COL_R and abs(part0 - PART_N) <= TOL_PART
print("  R0a flux step7 / flux OCIM - 1 = %.2e  %s" % (mass(J) / mass(Jf) - 1, "PASS" if r0a else "FAIL"))
print("  R0b active cells %d (expected %d)  %s" % (act.size, N_ACT_Q, "PASS" if r0b else "FAIL"))
print("  R0c mean depth %.3f m (expected %.1f)  %s" % (D0, D_Q, "PASS" if r0c else "FAIL"))
print("  R0d 1st column (%d,%d) in R indices, share %.3f %% (expected %s, %.1f %%)  %s"
      % (i0 + 1, j0 + 1, part0, COL_R, PART_N, "PASS" if r0d else "FAIL"))
if not (r0a and r0b and r0c and r0d):
    sys.exit("\n>>> AN R0 LOCK HAS FAILED: nothing is built.")

# ---- partition -----------------------------------------------------------------
inside = (ii == i0) & (jj == j0)                 # all wet cells of the column, all layers
print("\n===== Partition: column (%d,%d) in R indices, %d wet cells =====" % (i0 + 1, j0 + 1, int(inside.sum())))

def split_ab(dd):
    a, b = {}, {}
    for k in var_keys:
        v = np.asarray(dd[k])
        if k not in vect:
            a[k], b[k] = v, v
            continue
        x = v.ravel()
        ya, yb = np.where(inside, x, 0).astype(x.dtype), np.where(inside, 0, x).astype(x.dtype)
        a[k], b[k] = ya.reshape(v.shape), yb.reshape(v.shape)
    return a, b

def locks(ref, a, b, title):
    ok = True
    print("\n  --- %s ---" % title)
    for k in vect:
        # compared in float64: exact for double floats and integers < 2^53; types may
        # change on re-reading (logical -> uint8), values may not.
        x = np.asarray(ref[k]).ravel().astype(float)
        s = np.asarray(a[k]).ravel().astype(float) + np.asarray(b[k]).ravel().astype(float)
        p1 = np.array_equal(s, x, equal_nan=True)
        ok &= p1
        if not p1:
            print("  P1 %-28s FAIL (max |gap| %.3e)" % (k, float(np.max(np.abs(s - x)))))
    print("  P1 F_A + F_B == F_ref, %d fields: %s" % (len(vect), "PASS" if ok else "FAIL"))
    ja, jb = np.asarray(a["Jdredge"], dtype=float).ravel(), np.asarray(b["Jdredge"], dtype=float).ravel()
    p2 = not np.any((ja != 0) & (jb != 0))
    na, nb = int(np.count_nonzero(ja)), int(np.count_nonzero(jb))
    nc = int(np.count_nonzero(J[inside]))
    p3 = na == nc and nb == N_ACT_Q - nc and not np.any(ja[~inside] != 0)
    pa = 100.0 * mass(ja) / mass(J)
    p4 = abs(pa - PART_R0D) <= TOL_P4
    print("  P2 no active cell in both: %s" % ("PASS" if p2 else "FAIL"))
    print("  P3 active cells A %d (column: %d), B %d (expected %d): %s" % (na, nc, nb, N_ACT_Q - nc, "PASS" if p3 else "FAIL"))
    print("  P4 share of F_A %.4f %% (expected %.3f within %.3f point): %s" % (pa, PART_R0D, TOL_P4, "PASS" if p4 else "FAIL"))
    print("  (outside verdict) mean depth A %.3f m, B %.3f m, ref %.3f m" % (mean_depth(ja), mean_depth(jb), D0))
    return ok and p2 and p3 and p4

A, B = split_ab(d)
print("\n  layers of the column carrying forcing (R indices, depth, share of the column):")
for n in np.flatnonzero(inside & (J != 0)):
    print("    k=%2d  z=%8.2f m  %.4f" % (kk[n] + 1, Z[n], mc[n] / m0))
if not locks(d, A, B, "on the arrays built"):
    sys.exit("\n>>> A LOCK P FAILED: no file written.")

# ---- writing and re-reading ----------------------------------------------------
outputs = {}
for item_name, sd, dd in (("colA", "step7_smooth_colA", A), ("resteB", "step7_smooth_resteB", B)):
    od = os.path.join(BASE, sd)
    os.makedirs(od, exist_ok=True)
    out = os.path.join(od, "jdredge_ocim2_48l_%s_20260928.mat" % item_name)
    savemat(out, dd, do_compression=True)
    outputs[item_name] = out
    print("  written: %s (%.2f MB)" % (out, os.path.getsize(out) / 1e6))
rA, rB = loadmat(outputs["colA"]), loadmat(outputs["resteB"])
if not locks(d, rA, rB, "on the re-read files"):
    for f in outputs.values():
        os.rename(f, f + ".FAILED")
    sys.exit("\n>>> RE-READING MADE A LOCK FAIL: files renamed to .FAILED, no PARTITION_PASS.")
open(os.path.join(BASE, "PARTITION_PASS"), "w").write("%s\n%s\n" % (outputs["colA"], outputs["resteB"]))
print("\n>>> PARTITION READY (PARTITION_PASS written).")
