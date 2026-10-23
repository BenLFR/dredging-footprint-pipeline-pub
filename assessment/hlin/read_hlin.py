# read_hlin.py -- new test of H-lin: reading of set-ups A and B (rule written beforehand).
# Derived from read_d3.py: same calculation of the balance and of the AF (summary_q.py, test Q).
# Changes compared with read_d3.py, decided in the rule before any run:
#   - VS, spatial identity: Jdredge_forcing (OCIM 1x) == Jdredge (step7 input), cell by cell;
#     tested first on the M4 reference pair: if it fails there, its premise is false and we stop;
#   - L4 (homogeneity 100x / 1x) is reported and no longer masks the 1x AF;
#   - B: additive flux, share of F_A, PREC (balance <= 0.001 %) and E_B.
# ORDER: reference V0 and VS first; then, for each case, its locks, and the substance only if they pass.
# Usage: python3 read_hlin.py <RUN> <BASE_D3> <BASE>   (sbatch, never on the interactive node)
import csv, glob, os, sys
import numpy as np, h5py
from scipy.io import loadmat

RUN, BASE_D3, BASE = sys.argv[1], sys.argv[2], sys.argv[3]
RHO, MATM, GPMOL, NSTEP = 1025.0, 1.8e20, 12.0, 321   # ns_step_eb.m (OCIM code) / co2model_batch.m
TOL_BALANCE, TOL_SOURCE, TOL_CONTROL, TOL_LIN = 0.1, 1e-9, 1e-12, 5e-3   # L1 (Q0a), L2 (Q0b), L3 (Q0c), L4 (linearity)
TOL_VS = 1e-9          # VS: max |gap| <= 1e-9 x max |J|; design choice, same order as R0a and Q0b
TOL_V0 = 1e-4          # V0: AF points; the AFs of test Q are printed to 4 decimals
TOL_PREC = 0.001       # B, PREC: |balance| in % (1e-5 relative = 0.001 AF point, ten times below T1)
PART_A, TOL_PART = 32.553, 0.001                     # R0d of d3_move.py, in %
T1, G_AF = 0.01, 54.4364 - 49.8915                   # T1: design choice, from an independent review of the test; G: test Q
T2 = 0.01 * G_AF                                     # 0.045449 point: 0.01 G
AF_Q = {"smooth": 54.4364, "flag": 49.8915}          # test Q (summary_q.py)
AMPS = ("1x", "10x", "100x")

def T(a):
    return np.array(a).T

def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        return None, "%d files for %s" % (len(f), pat)
    return f[0], ""

def last_row(f, name, m):
    d = f[name]
    if d.shape[-1] == m:
        return np.asarray(d[-1, :], dtype=float)
    if d.shape[0] == m:
        return np.asarray(d[:, -1], dtype=float)
    sys.exit("ABORT: %s shape %s incompatible with m=%d" % (name, d.shape, m))

def geom(f):
    M3d = T(f["M3d"]).ravel(order="F")
    g = f["grid"]
    V3 = (T(g["DXT3d"]) * T(g["DYT3d"]) * T(g["DZT3d"])).ravel(order="F")
    Z3 = T(g["ZT3d"]).ravel(order="F")
    iocn = np.flatnonzero(M3d == 1)
    return V3[iocn], Z3[iocn]

def measure(base, pert):
    with h5py.File(base, "r") as fb, h5py.File(pert, "r") as fp:
        V, Z = geom(fp)
        m = V.size
        J = np.asarray(fp["Jdredge_forcing"], dtype=float).ravel()
        Mw = V * RHO
        w = J * Mw
        q = float(np.sum(w) * GPMOL * 1e-6)
        C_inj = q * NSTEP
        pa_b = float(np.asarray(fb["pco2a"]).ravel()[-1])
        pa_p = float(np.asarray(fp["pco2a"]).ravel()[-1])
        dC_atm = (pa_p - pa_b) * 1e-6 * MATM * GPMOL
        dC_oce = float(np.sum((last_row(fp, "DIC", m) - last_row(fb, "DIC", m)) * Mw) * GPMOL * 1e-6)
        act = J != 0
        return dict(q=q, C_inj=C_inj, dC_atm=dC_atm, pa_b=pa_b, J=J,
                    gap=(dC_atm + dC_oce - C_inj) / C_inj * 100.0, AF=dC_atm / C_inj * 100.0,
                    n_act=int(act.sum()), D=float(np.sum(w[act] * Z[act]) / np.sum(w[act])))

def run(d):
    b, e = only_one(d + "/ocim_baseline_*.mat")
    if b is None:
        return None, e
    out = {}
    for a in AMPS:
        p, e = only_one(d + "/ocim_dredge_%s_*.mat" % a)
        if p is None:
            return None, e
        out[a] = measure(b, p)
    return out, ""

def j_step7(d):
    f, e = only_one(d + "/jdredge_ocim2_48l_*.mat")
    if f is None:
        return None, e
    return np.asarray(loadmat(f)["Jdredge"], dtype=float).ravel(), f

def vs(Jf, J7):
    if Jf.size != J7.size:
        return False, float("nan")
    e = float(np.max(np.abs(Jf - J7)) / np.max(np.abs(J7)))
    return e <= TOL_VS, e

def ok(b):
    return "PASS" if b else "FAIL"

# ---- V0 and VS on the M4 references -----------------------------------------------
print("==================== REFERENCE LOCKS ====================")
refs, j7ref, all_ok = {}, {}, True
for lab in ("smooth", "flag"):
    R, e = run("%s/ocim_%s_m4" % (RUN, lab))
    if R is None:
        sys.exit("ABORT: reference %s unreadable (%s)" % (lab, e))
    J7, f7 = j_step7("%s/step7_%s_m4" % (RUN, lab))
    if J7 is None:
        sys.exit("ABORT: reference step7 %s (%s)" % (lab, f7))
    v0 = abs(R["1x"]["AF"] - AF_Q[lab]) <= TOL_V0
    p, e = vs(R["1x"]["J"], J7)
    print("  V0 %-6s AF 1x %.4f %% (test Q %.4f)  %s" % (lab, R["1x"]["AF"], AF_Q[lab], ok(v0)))
    print("  VS %-6s max|Jforcing - Jstep7| / max|J| = %.2e  %s   (%s)" % (lab, e, ok(p), os.path.basename(f7)))
    print("     (not part of the verdict) balance %.6f %% ; L4 %.2e" % (R["1x"]["gap"], R["100x"]["AF"] / R["1x"]["AF"] - 1))
    all_ok &= v0 and p
    refs[lab] = R
if not all_ok:
    sys.exit("\n>>> A REFERENCE LOCK FAILED. If it is VS, its premise is false (OCIM does not copy the "
             "vector as is): no substance is read.")

bounds = {}
for r in csv.DictReader(open(BASE_D3 + "/etape_S.csv")):
    if r.get("TV") not in (None, "", "NA"):
        bounds[(r["label"], r["variante"])] = 100.0 * float(r["TV"])

def locks_run(item_name, Rv, R0, J7, f7, rate=True):
    l1 = abs(Rv["1x"]["gap"]) <= TOL_BALANCE
    l3 = abs(Rv["1x"]["pa_b"] - R0["1x"]["pa_b"]) <= TOL_CONTROL
    p, e = vs(Rv["1x"]["J"], J7)
    print("  [%s] L1 balance %.6f %%  %s" % (item_name, Rv["1x"]["gap"], ok(l1)))
    print("  [%s] L3 control: gap %.3e ppm  %s" % (item_name, Rv["1x"]["pa_b"] - R0["1x"]["pa_b"], ok(l3)))
    print("  [%s] VS max|Jforcing - Jstep7| / max|J| = %.2e  %s   (%s)" % (item_name, e, ok(p), os.path.basename(f7)))
    res = l1 and l3 and p
    if rate:
        l2 = abs(Rv["1x"]["q"] / R0["1x"]["q"] - 1) <= TOL_SOURCE
        print("  [%s] L2 flux / M4 - 1 = %.2e  %s" % (item_name, Rv["1x"]["q"] / R0["1x"]["q"] - 1, ok(l2)))
        res &= l2
    print("  [%s] (reported, does not block) L4 AF_100x/AF_1x - 1 = %.2e  %s"
          % (item_name, Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1, "within 5e-3" if abs(Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1) <= TOL_LIN else "OUTSIDE 5e-3"))
    return res

# ---- A: K = 1 measured directly ------------------------------------------------
for lab in ("smooth", "flag"):
    print("\n==================== A: %s (100, 1) ====================" % lab)
    Rv, e = run("%s/ocim_%s_l100_k1" % (BASE, lab))
    J7, f7 = j_step7("%s/step7_%s_l100_k1" % (BASE_D3, lab))
    if Rv is None or J7 is None:
        print("  NOT READ: %s" % (e if Rv is None else f7))
        continue
    R0 = refs[lab]
    if not locks_run("A-" + lab, Rv, R0, J7, f7):
        print("  >>> A LOCK FAILED: the AF of this case is not read.")
        continue
    dAF = Rv["1x"]["AF"] - R0["1x"]["AF"]
    Rr = abs(dAF) / G_AF
    category = "negligible" if Rr <= 0.10 else "second order, reported" if Rr < 0.50 else "FIRST ORDER -> FAIL"
    print("  AF %.4f %% against M4 %.4f %%: dAF = %+.4f point ; R = %.4f -> %s" % (Rv["1x"]["AF"], R0["1x"]["AF"], dAF, Rr, category))
    b = bounds.get((lab, "l100_k1"))
    if b is None:
        print("  bound check: TV missing from etape_S.csv for (%s, l100_k1) -> NOT READ" % lab)
    else:
        print("  bound check: |dAF| %.4f <= 100 x TV = %.4f ? %s" % (abs(dAF), b,
              "YES" if abs(dAF) <= b else "NO -> the bound argument fails; bounds K = 4 and 64 withdrawn"))
        if lab == "flag" and 0.1 * G_AF < abs(dAF) <= b:
            print("  (note of the rule) effect between 0.1 G and the bound: negligibility rejected without contradicting the bound")
    print("  (not part of the verdict) D %.2f m against %.2f ; cells %d against %d ; dAF 10x %+.4f ; dAF 100x %+.4f"
          % (Rv["1x"]["D"], R0["1x"]["D"], Rv["1x"]["n_act"], R0["1x"]["n_act"],
             Rv["10x"]["AF"] - R0["10x"]["AF"], Rv["100x"]["AF"] - R0["100x"]["AF"]))

# ---- B: additivity of the column / rest partition (smooth) -----------------------
print("\n==================== B: smooth, dominant column + rest ====================")
if not os.path.exists(BASE + "/PARTITION_PASS"):
    sys.exit("  NOT READ: PARTITION_PASS missing (see logs/partition.out)\nEND")
RA, eA = run(BASE + "/ocim_smooth_colA")
RB, eB = run(BASE + "/ocim_smooth_resteB")
JA, fA = j_step7(BASE + "/step7_smooth_colA")
JB, fB = j_step7(BASE + "/step7_smooth_resteB")
if RA is None or RB is None or JA is None or JB is None:
    sys.exit("  NOT READ: %s %s %s %s\nEND" % (eA, eB, fA if JA is None else "", fB if JB is None else ""))
R0 = refs["smooth"]
va = locks_run("B-colA", RA, R0, JA, fA, rate=False)
vb = locks_run("B-resteB", RB, R0, JB, fB, rate=False)
jsum, e_add = RA["1x"]["J"] + RB["1x"]["J"], float("nan")
add = jsum.size == R0["1x"]["J"].size
if add:
    e_add = float(np.max(np.abs(jsum - R0["1x"]["J"])) / np.max(np.abs(R0["1x"]["J"])))
    add = e_add <= TOL_VS
print("  [B] VS+ Jforcing(A) + Jforcing(B) = Jforcing(ref), max|gap|/max|J| = %.2e  %s" % (e_add, ok(add)))
l2s = abs((RA["1x"]["q"] + RB["1x"]["q"]) / R0["1x"]["q"] - 1) <= TOL_SOURCE
pa = 100.0 * RA["1x"]["q"] / R0["1x"]["q"]
l2p = abs(pa - PART_A) <= TOL_PART
print("  [B] L2 (qA + qB) / q - 1 = %.2e  %s" % ((RA["1x"]["q"] + RB["1x"]["q"]) / R0["1x"]["q"] - 1, ok(l2s)))
print("  [B] L2 share qA / q = %.4f %% (expected %.3f within %.3f point)  %s" % (pa, PART_A, TOL_PART, ok(l2p)))
if not (va and vb and add and l2s and l2p):
    sys.exit("  >>> A LOCK FAILED: E_B is not read.\nEND")
prec = all(abs(R["1x"]["gap"]) <= TOL_PREC for R in (R0, RA, RB))
print("  [B] PREC |balance| <= %.3f %%: ref %.6f, A %.6f, B %.6f  %s" % (TOL_PREC, R0["1x"]["gap"], RA["1x"]["gap"], RB["1x"]["gap"], ok(prec)))
if not prec:
    sys.exit("  >>> PREC FAILS: B is NULL (instrument too imprecise for T1), not failed. E_B is not read.\nEND")
EB = {a: 100.0 * (R0[a]["dC_atm"] - RA[a]["dC_atm"] - RB[a]["dC_atm"]) / R0[a]["C_inj"] for a in AMPS}
e = abs(EB["1x"])
verdict = "PASS" if e <= T1 else "PASSED WITH ERROR (bounds + |E_B|)" if e <= T2 else "FAIL"
print("  AF ref %.4f %% ; AF colA %.4f %% ; AF resteB %.4f %%" % (R0["1x"]["AF"], RA["1x"]["AF"], RB["1x"]["AF"]))
print("  E_B (1x) = %+.5f point ; T1 = %.2f ; T2 = %.6f -> %s" % (EB["1x"], T1, T2, verdict))
print("  (not part of the verdict) E_B 10x %+.5f ; E_B 100x %+.5f point" % (EB["10x"], EB["100x"]))
print("  (not part of the verdict) D colA %.2f m, resteB %.2f m, ref %.2f m ; cells %d + %d = %d (ref %d)"
      % (RA["1x"]["D"], RB["1x"]["D"], R0["1x"]["D"], RA["1x"]["n_act"], RB["1x"]["n_act"],
         RA["1x"]["n_act"] + RB["1x"]["n_act"], R0["1x"]["n_act"]))
print("\nEND")
