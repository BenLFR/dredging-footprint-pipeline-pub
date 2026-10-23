# read_svr.py -- reading of the OCIM run dedicated to the annual SVR cap (rule written beforehand).
# Derived from read_d3.py, itself derived from summary_q.py (test Q): same balance, same AF.
# ORDER IMPOSED: locks L0-L4 first; if one fails, the AF is not printed.
#   L0  AF 1x re-read on the M4 runs = 54.4364 / 49.8915 to 1e-4 point (test Q printed to 4 decimals)
#   L1  mass balance <= 0.1 %
#   L2  injected flux / M4 flux = Q (BASE/Q_<lab>.txt, svr_forcing.R) to 1e-9 relative
#   L3  same control: final baseline pCO2 identical to that of M4, to 1e-12 ppm
#   L4  linearity: |AF_100x / AF_1x - 1| <= 5e-3
# Substance: dAF = AF_annual - AF_M4 (1x, 321 years); R = |dAF| / G, G = 4.5449 points (test Q).
#   R <= 0.10: does not depend on the accounting; 0.10 < R < 0.50: second order, reported; R >= 0.50: FIRST ORDER.
# Outside verdict: dpCO2_annual / dpCO2_M4 against Q; depth; cells. H-lin check: |dAF| <= TV x 100 (BASE/TV_<lab>.txt).
# Execution: sbatch (python3 + h5py), never on the interactive node.
import os, glob, sys
import numpy as np, h5py

RHO, MATM, GPMOL, NSTEP = 1025.0, 1.8e20, 12.0, 321   # ns_step_eb.m (OCIM code) / co2model_batch.m
TOL_REF, TOL_BALANCE, TOL_SOURCE, TOL_CONTROL, TOL_LIN = 1e-4, 0.1, 1e-9, 1e-12, 5e-3
AF_Q = {"smooth": 54.4364, "flag": 49.8915}            # test Q (summary_q.py)
G_AF = AF_Q["smooth"] - AF_Q["flag"]

RUN = os.environ["DFP_SCRATCH"] + "/output_V6_sala2021_r6_flux10"
BASE = os.environ["DFP_SCRATCH"] + "/svr_annuel_ocim_20261005"
AMPS = ("1x", "10x", "100x")

def T(a):
    return np.array(a).T

def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit("ABORT: %d files for %s" % (len(f), pat))
    return f[0]

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
        return dict(q=q, pa_b=pa_b, dpco2=pa_p - pa_b, gap=(dC_atm + dC_oce - C_inj) / C_inj * 100.0,
                    AF=dC_atm / C_inj * 100.0, n_act=int(act.sum()),
                    D=float(np.sum(w[act] * Z[act]) / np.sum(w[act])))

def run(d):
    b = only_one(d + "/ocim_baseline_*.mat")
    return {a: measure(b, only_one(d + "/ocim_dredge_%s_*.mat" % a)) for a in AMPS}

for lab in ("smooth", "flag"):
    print("\n==================== %s ====================" % lab)
    R0 = run("%s/ocim_%s_m4" % (RUN, lab))
    l0 = abs(R0["1x"]["AF"] - AF_Q[lab]) <= TOL_REF
    print("  L0 AF M4 re-read %.6f %% (test Q %.4f)  %s" % (R0["1x"]["AF"], AF_Q[lab], "PASS" if l0 else "FAIL"))
    if not l0:
        print("  >>> L0 FAILED: the reading does not reproduce the reference; nothing is read for %s." % lab)
        continue
    Q = float(open("%s/Q_%s.txt" % (BASE, lab)).read())
    TV = float(open("%s/TV_%s.txt" % (BASE, lab)).read())
    Rv = run("%s/ocim_%s_an" % (BASE, lab))
    l1 = abs(Rv["1x"]["gap"]) <= TOL_BALANCE
    l2 = abs(Rv["1x"]["q"] / R0["1x"]["q"] / Q - 1) <= TOL_SOURCE
    l3 = abs(Rv["1x"]["pa_b"] - R0["1x"]["pa_b"]) <= TOL_CONTROL
    l4 = abs(Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1) <= TOL_LIN
    print("  L1 balance %.6f %%  %s" % (Rv["1x"]["gap"], "PASS" if l1 else "FAIL"))
    print("  L2 flux / M4 = %.12f ; Q = %.12f ; gap %.2e  %s" % (Rv["1x"]["q"] / R0["1x"]["q"], Q,
          Rv["1x"]["q"] / R0["1x"]["q"] / Q - 1, "PASS" if l2 else "FAIL"))
    print("  L3 control: gap %.3e ppm  %s" % (Rv["1x"]["pa_b"] - R0["1x"]["pa_b"], "PASS" if l3 else "FAIL"))
    print("  L4 linearity AF_100x/AF_1x - 1 = %.2e  %s" % (Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1, "PASS" if l4 else "FAIL"))
    if not (l1 and l2 and l3 and l4):
        print("  >>> A LOCK HAS FAILED: the AF is not read.")
        continue
    dAF = Rv["1x"]["AF"] - R0["1x"]["AF"]
    Rr = abs(dAF) / G_AF
    category = ("does not depend on the accounting" if Rr <= 0.10 else "second order, reported" if Rr < 0.50
              else "FIRST ORDER")
    print("  AF %.4f %% against M4 %.4f %%: dAF = %+.4f point ; R = %.4f -> %s" % (Rv["1x"]["AF"], R0["1x"]["AF"], dAF, Rr, category))
    print("  (not part of the verdict) dpCO2 1x: %.6e ppm against %.6e ; ratio %.6f against Q %.6f" % (
          Rv["1x"]["dpco2"], R0["1x"]["dpco2"], Rv["1x"]["dpco2"] / R0["1x"]["dpco2"], Q))
    print("  (not part of the verdict) D %.2f m against %.2f ; cells %d against %d" % (Rv["1x"]["D"], R0["1x"]["D"], Rv["1x"]["n_act"], R0["1x"]["n_act"]))
    print("  H-lin check: |dAF| %.4f <= TV x 100 = %.4f ? %s" % (abs(dAF), 100 * TV,
          "YES" if abs(dAF) <= 100 * TV else "NO -> H-lin FAILS, bounds of step S of the D3 test withdrawn"))
print("\nEND")
