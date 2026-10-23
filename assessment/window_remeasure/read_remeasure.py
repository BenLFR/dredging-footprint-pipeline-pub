# -*- coding: utf-8 -*-
"""
read_remeasure.py -- reading of the QF (injection depth) and DF (move of one column) tests, under the
window of Atwood et al. (2024) and on the capped source (reading rules written beforehand).
Derived from read_window.py: same balance, same AF (mole ratio, 12 g/mol on both sides), same thresholds.
ORDER IMPOSED: locks G0-G5 of each variant first; if one fails, the AF of the variant is not read.
  G0  same control: pCO2 2100 of the baseline = that of the window reference of the label, to 1e-12 ppm
  G1  window: 10 active steps, 2015-2024, scale 0.1; zero anomaly before 2015, non-zero in 2015 (= F1)
  G2  same source: mass injected in 1x = that of the window reference, to 1e-9 relative
  G3  mass balance in 2100 (1x) <= 0.1 %
  G4  linearity in 2100 <= 5e-3, with the F4 fallback (substance read on the 10x if only the 1x fails)
  G5  the forcing of the run differs from that of the reference (added before launch: without it, a
      FORCING_DIR error would read ΔAF = 0 as a result)
Substance, verdicts:
  QF  ΔAF_QF = AF_bottom - AF_ref, end 2100, per label. Predicted sign < 0. FAILS if >= 0 for at least one
      label. Material magnitude if |ΔAF_QF| >= 1 point for at least one label.
  DF  ΔAF_DF = AF_shift - AF_ref,smooth, end 2100. PASSES if |ΔAF_DF| <= 1 point, FAILS otherwise.
Outside verdict: AF at end 2024, 2033, 2050, 2100 and maximum; R = |ΔAF_QF| / G, G = AF_ref,smooth - AF_ref,flag;
  mean depths of the forcings; values under constant forcing (test Q, D3), without judgement.
Execution: sbatch (python3 + h5py), never on the interactive node.
"""
import os, glob, sys
import numpy as np, h5py

RHO, MATM, GPMOL = 1025.0, 1.8e20, 12.0     # ns_step_eb.m; co2model_batch.m (12 g/mol)
TOL_CONTROL, TOL_SOURCE, TOL_BALANCE, TOL_LIN = 1e-12, 1e-9, 0.1, 5e-3   # F0, F2, F3, F4 of read_window.py
THRESHOLD_PT = 1.0                               # AF time point: modelling choice, no source
YEARS = list(range(2015, 2025))
SCALE = 0.1
MILESTONES = (2024, 2033, 2050, 2100)            # as in read_window.py
AMPS = ("1x", "10x", "100x")

WIN = os.environ["DFP_SCRATCH"] + "/ocim_fenetre_20261009"
BASE = os.environ["DFP_SCRATCH"] + "/remesure_fen_20261010"
VARIANTS = (("bottom_smooth", "smooth", "ocim_smooth_bot"),
             ("bottom_flag", "flag", "ocim_flag_bot"),
             ("shift_smooth", "smooth", "ocim_smooth_shift"))
# Values under constant forcing 1780-2100, source before the cap (M4); descriptive comparison only
TEST_Q = {"smooth": (-18.62, 814.4, 26.0), "flag": (-15.65, 908.8, 24.3)}   # test Q: ΔAF, D_bottom, D_seafloor
D3_SHIFT, AUG_SHIFT = 0.0072, 0.0965        # D3 test (read_d3.py); earlier value
AF_WIN = {"smooth": 63.26, "flag": 55.82}    # rounded values of test B (read_window.py)


def T(a):
    return np.array(a).T


def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]


def vec(f, name):
    return np.asarray(f[name], dtype=float).ravel()


def line(f, name, m, k):
    d = f[name]
    if d.shape[-1] == m:
        return np.asarray(d[k, :], dtype=float)
    if d.shape[0] == m:
        return np.asarray(d[:, k], dtype=float)
    sys.exit("ABORT: %s shape %s incompatible with m=%d" % (name, d.shape, m))


def geom(f):
    M3d = T(f["M3d"]).ravel(order="F")
    g = f["grid"]
    V3 = (T(g["DXT3d"]) * T(g["DYT3d"]) * T(g["DZT3d"])).ravel(order="F")
    Z3 = T(g["ZT3d"]).ravel(order="F")
    iocn = np.flatnonzero(M3d == 1)
    return V3[iocn], Z3[iocn]


def rate(fp, V):
    J = np.asarray(fp["Jdredge_forcing"], dtype=float).ravel()
    return float(np.sum(J * V * RHO) * GPMOL * 1e-6)


def depth(fp, V, Z):
    w = np.asarray(fp["Jdredge_forcing"], dtype=float).ravel() * V * RHO
    return float(np.sum(Z * w) / np.sum(w))


def measure(fb, fp, V, k):
    m = V.size
    act = vec(fp, "actif") != 0
    C_inj = rate(fp, V) * int(act[:k + 1].sum())
    dpco2 = float(vec(fp, "pco2a")[k] - vec(fb, "pco2a")[k])
    dC_atm = dpco2 * 1e-6 * MATM * GPMOL
    dC_oce = float(np.sum((line(fp, "DIC", m, k) - line(fb, "DIC", m, k)) * V * RHO) * GPMOL * 1e-6)
    return dict(C_inj=C_inj, dpco2=dpco2, dC_atm=dC_atm, dC_oce=dC_oce, AF=dC_atm / C_inj * 100.0,
                gap=(dC_atm + dC_oce - C_inj) / C_inj * 100.0)


def open_pair(d):
    fb = h5py.File(only_one(d + "/ocim_baseline_*.mat"), "r")
    fp = {a: h5py.File(only_one(d + "/ocim_dredge_%s_*.mat" % a), "r") for a in AMPS}
    return fb, fp


def read_run(item_name, d, ref=None):
    """Locks then substance of a set of four scenarios. ref = (fb, fp, V) of the reference, or None for itself."""
    print("\n==================== %s ====================" % item_name)
    print("  folder: %s" % d)
    fb, fp = open_pair(d)
    V, Z = geom(fp["1x"])
    yr = vec(fp["1x"], "annee_pas")
    k_end = len(yr) - 1
    R = {a: measure(fb, fp[a], V, k_end) for a in AMPS}
    D = depth(fp["1x"], V, Z)
    ok = True
    if ref is not None:
        rfb, rfp, rV = ref
        # G0 -- same control
        pb, pr = float(vec(fb, "pco2a")[-1]), float(vec(rfb, "pco2a")[-1])
        g0 = abs(pb - pr) <= TOL_CONTROL
        print("  G0 control: pCO2 2100 %.10f against %.10f (reference) ; gap %.3e  %s" % (pb, pr, pb - pr,
              "PASS" if g0 else "FAIL"))
        ok &= g0
    # G1 -- window
    g1 = True
    for a in AMPS:
        act = vec(fp[a], "actif") != 0
        yrs = sorted(int(x) for x in vec(fp[a], "annee_pas")[act])
        scale = float(vec(fp[a], "dredge_scale")[0])
        i0 = int(np.flatnonzero(act)[0])
        dp = vec(fp[a], "pco2a") - vec(fb, "pco2a")
        before, start = float(np.max(np.abs(dp[:i0]))), float(dp[i0])
        good = yrs == YEARS and scale == SCALE and before == 0.0 and start != 0.0
        g1 &= good
        print("  G1 %-4s: active years %d-%d (%d) ; scale %g ; max|dpCO2| before 2015 = %.1e ; dpCO2 end 2015 = %.3e  %s"
              % (a, yrs[0], yrs[-1], len(yrs), scale, before, start, "PASS" if good else "FAIL"))
    ok &= g1
    if ref is not None:
        # G2 -- same source
        Rr = measure(rfb, rfp["1x"], rV, k_end)
        g2 = abs(R["1x"]["C_inj"] / Rr["C_inj"] - 1) <= TOL_SOURCE
        print("  G2 source: injected %.6e g against %.6e g (reference) ; gap %.2e  %s" % (
              R["1x"]["C_inj"], Rr["C_inj"], R["1x"]["C_inj"] / Rr["C_inj"] - 1, "PASS" if g2 else "FAIL"))
        ok &= g2
    # G3 -- mass balance
    g3 = abs(R["1x"]["gap"]) <= TOL_BALANCE
    print("  G3 balance 2100 (1x) %.6f %%  %s" % (R["1x"]["gap"], "PASS" if g3 else "FAIL"))
    ok &= g3
    # G4 -- linearity, and fallback
    l10 = R["10x"]["AF"] / R["1x"]["AF"] - 1
    l100 = R["100x"]["AF"] / R["1x"]["AF"] - 1
    l100_10 = R["100x"]["AF"] / R["10x"]["AF"] - 1
    g4 = abs(l10) <= TOL_LIN and abs(l100) <= TOL_LIN
    print("  G4 linearity: AF10/AF1 - 1 = %.2e ; AF100/AF1 - 1 = %.2e ; AF100/AF10 - 1 = %.2e  %s" % (
          l10, l100, l100_10, "PASS" if g4 else "FAIL"))
    if ref is not None:
        # G5 -- the forcing differs from that of the reference
        J = np.asarray(fp["1x"]["Jdredge_forcing"], dtype=float).ravel()
        Jr = np.asarray(rfp["1x"]["Jdredge_forcing"], dtype=float).ravel()
        ndiff = int(np.count_nonzero(J != Jr))
        g5 = ndiff > 0
        print("  G5 forcing different from the reference: %d different cells  %s" % (ndiff, "PASS" if g5 else "FAIL"))
        ok &= g5
    print("  (outside verdict) mean depth of the forcing %.2f m ; active cells %d" % (
          D, int(np.count_nonzero(np.asarray(fp["1x"]["Jdredge_forcing"])))))
    if not ok:
        print("  >>> A LOCK HAS FAILED: the AF is not read.")
        return None
    if g4:
        key, div = "1x", 1.0
    elif abs(l100_10) <= TOL_LIN:
        key, div = "10x", 10.0
        print("  >>> FALLBACK: the 1x is not resolved; substance read on the 10x, masses and dpCO2 divided by 10.")
    else:
        print("  >>> G4 FAILED, fallback impossible: the AF is not read.")
        return None
    print("  Substance (run %s):" % key)
    AF = {}
    for y in MILESTONES:
        k = int(np.flatnonzero(yr == y)[0])
        r = measure(fb, fp[key], V, k)
        AF[y] = r["AF"]
        print("    end %d: AF %.4f %% ; dpCO2 %.4e ppm ; atm %.4e Tg ; ocean %.4e Tg ; injected %.4e Tg ; balance %.4f %%" % (
              y, r["AF"], r["dpco2"] / div, r["dC_atm"] / div / 1e12, r["dC_oce"] / div / 1e12, r["C_inj"] / div / 1e12,
              r["gap"]))
    dp = vec(fp[key], "pco2a") - vec(fb, "pco2a")
    act = vec(fp[key], "actif") != 0
    cum = np.cumsum(act) * rate(fp[key], V)
    i0 = int(np.flatnonzero(act)[0])
    af = np.full(dp.size, np.nan)
    af[i0:] = dp[i0:] * 1e-6 * MATM * GPMOL / cum[i0:] * 100.0
    kmax = int(np.nanargmax(af))
    print("    AF maximum: %.4f %% at the end of %d" % (af[kmax], int(yr[kmax])))
    return dict(AF=AF, D=D, fb=fb, fp=fp, V=V)


# ---- references (test B passed): re-read with the same functions, G0 and G2 not applicable ----
REF = {}
for lab in ("smooth", "flag"):
    r = read_run("window reference %s" % lab, "%s/ocim_%s_fen" % (WIN, lab))
    if r is None:
        sys.exit("\n>>> THE REFERENCE %s CANNOT BE RE-READ: nothing is read." % lab)
    print("  (check) AF 2100 re-read %.4f %% against %.2f %% expected" % (r["AF"][2100], AF_WIN[lab]))
    REF[lab] = r
G = REF["smooth"]["AF"][2100] - REF["flag"]["AF"][2100]

# ---- variants ----
RES = {}
for item_name, lab, sd in VARIANTS:
    rf = REF[lab]
    RES[item_name] = read_run(item_name, "%s/%s" % (BASE, sd), ref=(rf["fb"], rf["fp"], rf["V"]))

# ---- verdicts ----
print("\n==================== VERDICTS (rule as written) ====================")
print("  G = AF_ref,smooth - AF_ref,flag = %.4f points (expected 7.44)" % G)
print("\n  QF -- injection depth (bottom - seafloor), end 2100")
dq = {}
for lab in ("smooth", "flag"):
    r = RES["bottom_" + lab]
    if r is None:
        print("    %-6s: NOT READ (lock failed)" % lab)
        continue
    dq[lab] = r["AF"][2100] - REF[lab]["AF"][2100]
    q, Db, Ds = TEST_Q[lab]
    print("    %-6s: ΔAF_QF = %+.4f points ; R = |ΔAF|/G = %.3f ; D %.2f m (bottom) against %.2f m (seafloor)"
          % (lab, dq[lab], abs(dq[lab]) / G, r["D"], REF[lab]["D"]))
    print("             milestones ΔAF: %s" % " ; ".join("%d %+.4f" % (y, r["AF"][y] - REF[lab]["AF"][y]) for y in MILESTONES))
    print("             (constant forcing, M4, test Q: %+.2f points ; D %.1f against %.1f m)" % (q, Db, Ds))
if any(v >= 0 for v in dq.values()):
    print("  >>> QF FAILS: ΔAF_QF >= 0 for %s. Stop and diagnose."
          % ", ".join(l for l, v in dq.items() if v >= 0))
elif len(dq) == 2:
    amp = "MATERIAL" if any(abs(v) >= THRESHOLD_PT for v in dq.values()) else "secondary"
    print("  >>> QF PASSES (predicted sign, ΔAF_QF < 0 for both labels); magnitude %s (threshold %.0f point)." % (amp, THRESHOLD_PT))
else:
    print("  >>> QF INCOMPLETE: one label is not read; no verdict.")

print("\n  DF -- move by one cell of the dominant column (Label A), end 2100")
r = RES["shift_smooth"]
if r is None:
    print("  >>> DF NOT READ (lock failed).")
else:
    dd = r["AF"][2100] - REF["smooth"]["AF"][2100]
    try:
        print("    forcing: %s" % open(BASE + "/MOVE_PASS").read().strip())
    except OSError:
        print("    forcing: MOVE_PASS unreadable")
    print("    ΔAF_DF = %+.4f points ; milestones: %s" % (dd, " ; ".join("%d %+.4f" % (y, r["AF"][y] - REF["smooth"]["AF"][y]) for y in MILESTONES)))
    print("    (constant forcing: D3 under M4 %+.4f ; earlier value %+.4f)" % (D3_SHIFT, AUG_SHIFT))
    print("  >>> DF %s (|ΔAF_DF| %s %.0f point)." % ("PASS" if abs(dd) <= THRESHOLD_PT else "FAIL",
          "<=" if abs(dd) <= THRESHOLD_PT else ">", THRESHOLD_PT))
print("\nEND")
