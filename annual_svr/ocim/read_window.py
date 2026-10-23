# read_window.py -- reading of the OCIM runs following the protocol of Atwood et al. (rule written beforehand)
# (2024, §2.2): dredging 2015-2024 at the annual mass Mijk, then zero until 2100; 1x, 10x, 100x.
# Derived from read_svr.py: same balance, same AF (mole ratio, 12 g/mol on both sides).
# ORDER IMPOSED: locks F0-F4 first; if one fails, the AF is not printed.
#   F0  control: final baseline pCO2 identical to that of the capped runs, to 1e-12 ppm
#   F1  window: 10 active steps, years 2015-2024, scale 0.1; zero anomaly before 2015, non-zero in 2015
#   F2  source: mass injected in 1x = annual flux of the capped 1x (decadal total), to 1e-9 relative
#   F3  mass balance in 2100 (1x) <= 0.1 %
#   F4  linearity in 2100: |AF_10x / AF_1x - 1| and |AF_100x / AF_1x - 1| <= 5e-3
#       Fallback written beforehand: if F4 fails only through the 1x (|AF_100x / AF_10x - 1| <= 5e-3), the 1x is declared not
#       resolved and the substance is read on the 10x (masses and dpCO2 divided by 10), saying so. Otherwise nothing is read.
# Substance: AF, dpCO2 (ppm), atmosphere / ocean masses at end of year 2024, 2033, 2050, 2100; maximum of the AF.
# Execution: sbatch (python3 + h5py), never on the interactive node.
import os, glob, sys
import numpy as np, h5py

RHO, MATM, GPMOL = 1025.0, 1.8e20, 12.0     # ns_step_eb.m; co2model_batch.m (12 g/mol)
TOL_CONTROL, TOL_SOURCE, TOL_BALANCE, TOL_LIN = 1e-12, 1e-9, 0.1, 5e-3   # same thresholds as read_svr.py (L3, L2, L1, L4)
YEARS = list(range(2015, 2025))
SCALE = 0.1
MILESTONES = (2024, 2033, 2050, 2100)   # 2033: ~9 years after the end, the horizon of Atwood et al. (2024, p. 6-7)

SVR = os.environ["DFP_SCRATCH"] + "/svr_annuel_ocim_20261005"
BASE = os.environ["DFP_SCRATCH"] + "/ocim_fenetre_20261009"
AMPS = ("1x", "10x", "100x")


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
    iocn = np.flatnonzero(M3d == 1)
    return V3[iocn]


def rate(fp, V):
    """Flux of the recorded forcing, g C/yr (12 g/mol, as in read_svr.py)."""
    J = np.asarray(fp["Jdredge_forcing"], dtype=float).ravel()
    return float(np.sum(J * V * RHO) * GPMOL * 1e-6)


def measure(fb, fp, V, k):
    """Balance at state k (end of year annee_pas[k])."""
    m = V.size
    act = vec(fp, "actif") != 0
    C_inj = rate(fp, V) * int(act[:k + 1].sum())
    dpco2 = float(vec(fp, "pco2a")[k] - vec(fb, "pco2a")[k])
    dC_atm = dpco2 * 1e-6 * MATM * GPMOL
    dC_oce = float(np.sum((line(fp, "DIC", m, k) - line(fb, "DIC", m, k)) * V * RHO) * GPMOL * 1e-6)
    return dict(C_inj=C_inj, dpco2=dpco2, dC_atm=dC_atm, dC_oce=dC_oce, AF=dC_atm / C_inj * 100.0,
                gap=(dC_atm + dC_oce - C_inj) / C_inj * 100.0)


for lab in ("smooth", "flag"):
    print("\n==================== %s ====================" % lab)
    d = "%s/ocim_%s_fen" % (BASE, lab)
    fb = h5py.File(only_one(d + "/ocim_baseline_*.mat"), "r")
    fp = {a: h5py.File(only_one(d + "/ocim_dredge_%s_*.mat" % a), "r") for a in AMPS}
    V = geom(fp["1x"])
    yr = vec(fp["1x"], "annee_pas")
    k_end = len(yr) - 1
    # F0 -- control
    with h5py.File(only_one("%s/ocim_%s_an/ocim_baseline_*.mat" % (SVR, lab)), "r") as fs:
        pb_svr = float(vec(fs, "pco2a")[-1])
    pb = float(vec(fb, "pco2a")[-1])
    f0 = abs(pb - pb_svr) <= TOL_CONTROL
    print("  F0 control: pCO2 2100 %.10f against %.10f (cap) ; gap %.3e  %s" % (pb, pb_svr, pb - pb_svr,
          "PASS" if f0 else "FAIL"))
    # F1 -- window
    f1 = True
    for a in AMPS:
        act = vec(fp[a], "actif") != 0
        yrs = sorted(int(x) for x in vec(fp[a], "annee_pas")[act])
        scale = float(vec(fp[a], "dredge_scale")[0])
        i0 = int(np.flatnonzero(act)[0])
        dp = vec(fp[a], "pco2a") - vec(fb, "pco2a")
        before, start = float(np.max(np.abs(dp[:i0]))), float(dp[i0])
        ok = yrs == YEARS and scale == SCALE and before == 0.0 and start != 0.0
        f1 = f1 and ok
        print("  F1 %-4s: active years %d-%d (%d) ; scale %g ; max|dpCO2| before 2015 = %.1e ; dpCO2 end 2015 = %.3e  %s"
              % (a, yrs[0], yrs[-1], len(yrs), scale, before, start, "PASS" if ok else "FAIL"))
    # F2 -- source
    with h5py.File(only_one("%s/ocim_%s_an/ocim_dredge_1x_*.mat" % (SVR, lab)), "r") as fs:
        q_svr = rate(fs, geom(fs))
    R = {a: measure(fb, fp[a], V, k_end) for a in AMPS}
    f2 = abs(R["1x"]["C_inj"] / q_svr - 1) <= TOL_SOURCE
    print("  F2 source: injected %.6e g against capped flux %.6e g/yr ; gap %.2e  %s" % (
          R["1x"]["C_inj"], q_svr, R["1x"]["C_inj"] / q_svr - 1, "PASS" if f2 else "FAIL"))
    # F3 -- mass balance
    f3 = abs(R["1x"]["gap"]) <= TOL_BALANCE
    print("  F3 balance 2100 (1x) %.6f %%  %s" % (R["1x"]["gap"], "PASS" if f3 else "FAIL"))
    # F4 -- linearity, and fallback written beforehand
    l10 = R["10x"]["AF"] / R["1x"]["AF"] - 1
    l100 = R["100x"]["AF"] / R["1x"]["AF"] - 1
    l100_10 = R["100x"]["AF"] / R["10x"]["AF"] - 1
    f4 = abs(l10) <= TOL_LIN and abs(l100) <= TOL_LIN
    print("  F4 linearity: AF10/AF1 - 1 = %.2e ; AF100/AF1 - 1 = %.2e ; AF100/AF10 - 1 = %.2e  %s" % (
          l10, l100, l100_10, "PASS" if f4 else "FAIL"))
    if not (f0 and f1 and f2 and f3):
        print("  >>> A LOCK HAS FAILED: the AF is not read.")
        continue
    if f4:
        ref, div = "1x", 1.0
    elif abs(l100_10) <= TOL_LIN:
        ref, div = "10x", 10.0
        print("  >>> FALLBACK: the 1x is not resolved; substance read on the 10x, masses and dpCO2 divided by 10.")
    else:
        print("  >>> F4 FAILED, fallback impossible: the AF is not read.")
        continue
    print("  Substance (run %s):" % ref)
    for y in MILESTONES:
        k = int(np.flatnonzero(yr == y)[0])
        r = measure(fb, fp[ref], V, k)
        print("    end %d: AF %.4f %% ; dpCO2 %.4e ppm ; atm %.4e Tg ; ocean %.4e Tg ; injected %.4e Tg ; balance %.4f %%" % (
              y, r["AF"], r["dpco2"] / div, r["dC_atm"] / div / 1e12, r["dC_oce"] / div / 1e12, r["C_inj"] / div / 1e12,
              r["gap"]))
    dp = vec(fp[ref], "pco2a") - vec(fb, "pco2a")
    act = vec(fp[ref], "actif") != 0
    cum = np.cumsum(act) * rate(fp[ref], V)
    i0 = int(np.flatnonzero(act)[0])
    af = np.full(dp.size, np.nan)
    af[i0:] = dp[i0:] * 1e-6 * MATM * GPMOL / cum[i0:] * 100.0
    kmax = int(np.nanargmax(af))
    print("    AF maximum: %.4f %% at the end of %d" % (af[kmax], int(yr[kmax])))
print("\nEND")
