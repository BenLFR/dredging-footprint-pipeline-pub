# read_d3.py -- D3 under M4: reading of the OCIM runs (rule written beforehand).
# Derived from summary_q.py (test Q): same calculation of the balance and of the AF. Compares each variant with the
# M4 run of the same label (RUN/ocim_<lab>_m4).
# ORDER IMPOSED: locks L1-L4 first; if one fails for a variant, its AF is not printed.
#   L1  mass balance <= 0.1 % (Q0a)
#   L2  same injected flux as M4, to 1e-9 relative (Q0b)
#   L3  same control: final baseline pCO2 identical to that of M4, to 1e-12 ppm (Q0c)
#   L4  linearity: |AF_100x / AF_1x - 1| <= 5e-3 (interval [99.5; 100.5])
# Substance: dAF = AF_var - AF_M4 (1x, 321 years, 2.16 Pg C/ppm); R = |dAF| / G, G = 4.5449 points (test Q).
#   R <= 0.10 negligible; 0.10 < R < 0.50 reported; R >= 0.50 first order, FAILS.
# H-lin check (variants that have a bound at step S): |dAF| <= TV x 100.
# Outside verdict: rules of the pre-M4 D3 test (R1b 10 % on dpCO2 for lambda; 1 point for the displacement).
# Execution: sbatch (python3 + h5py), never on the interactive node.
import csv, glob, os, sys
import numpy as np, h5py

RHO, MATM, GPMOL, NSTEP = 1025.0, 1.8e20, 12.0, 321   # ns_step_eb.m (OCIM code) / co2model_batch.m
TOL_BALANCE, TOL_SOURCE, TOL_CONTROL, TOL_LIN = 0.1, 1e-9, 1e-12, 5e-3
G_AF = 54.4364 - 49.8915                              # points, test Q
AUG = {"l500_k16": {"smooth": 0.15, "flag": 0.17}, "shift": {"smooth": 0.0965}}   # results of the pre-M4 D3 test (points)

RUN = os.environ["DFP_SCRATCH"] + "/output_V6_sala2021_r6_flux10"
BASE = os.environ["DFP_SCRATCH"] + "/d3_m4_20260928"
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

bounds = {}
if os.path.exists(BASE + "/etape_S.csv"):
    for r in csv.DictReader(open(BASE + "/etape_S.csv")):
        if r.get("TV") not in (None, "", "NA"):
            bounds[(r["label"], r["variante"])] = 100.0 * float(r["TV"])

variants = sorted(os.path.basename(p)[len("ocim_"):] for p in glob.glob(BASE + "/ocim_*"))
print("variants found:", variants)
refs = {lab: run("%s/ocim_%s_m4" % (RUN, lab)) for lab in ("smooth", "flag")}
for lab, R in refs.items():
    print("reference M4 %s: AF 1x %.4f %% (test Q: %s)" % (lab, R["1x"]["AF"], {"smooth": "54,4364", "flag": "49,8915"}[lab]))

for v in variants:
    lab, tag = v.split("_", 1)
    tag = "shift" if tag.startswith("m4_shift") or tag == "m4shift" else tag
    print("\n==================== %s ====================" % v)
    Rv, R0 = run(BASE + "/ocim_" + v), refs[lab]
    l1 = abs(Rv["1x"]["gap"]) <= TOL_BALANCE
    l2 = abs(Rv["1x"]["q"] / R0["1x"]["q"] - 1) <= TOL_SOURCE
    l3 = abs(Rv["1x"]["pa_b"] - R0["1x"]["pa_b"]) <= TOL_CONTROL
    l4 = abs(Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1) <= TOL_LIN
    print("  L1 balance %.6f %%  %s" % (Rv["1x"]["gap"], "PASS" if l1 else "FAIL"))
    print("  L2 flux / M4 - 1 = %.2e  %s" % (Rv["1x"]["q"] / R0["1x"]["q"] - 1, "PASS" if l2 else "FAIL"))
    print("  L3 control: gap %.3e ppm  %s" % (Rv["1x"]["pa_b"] - R0["1x"]["pa_b"], "PASS" if l3 else "FAIL"))
    print("  L4 linearity AF_100x/AF_1x - 1 = %.2e  %s" % (Rv["100x"]["AF"] / Rv["1x"]["AF"] - 1, "PASS" if l4 else "FAIL"))
    if not (l1 and l2 and l3 and l4):
        print("  >>> A LOCK FAILED: the AF of this variant is not read.")
        continue
    dAF = Rv["1x"]["AF"] - R0["1x"]["AF"]
    Rr = abs(dAF) / G_AF
    category = ("negligible" if Rr <= 0.10 else "second order, reported" if Rr < 0.50 else "FIRST ORDER -> FAIL")
    print("  AF %.4f %% against M4 %.4f %%: dAF = %+.4f point ; R = %.4f -> %s" % (Rv["1x"]["AF"], R0["1x"]["AF"], dAF, Rr, category))
    print("  (not part of the verdict) D %.2f m against %.2f ; cells %d against %d ; relative dpCO2 %+.3f %%"
          % (Rv["1x"]["D"], R0["1x"]["D"], Rv["1x"]["n_act"], R0["1x"]["n_act"], 100 * (Rv["1x"]["dpco2"] / R0["1x"]["dpco2"] - 1)))
    if (lab, tag) in bounds:
        b = bounds[(lab, tag)]
        print("  H-lin check: |dAF| %.4f <= bound %.4f ? %s" % (abs(dAF), b, "YES" if abs(dAF) <= b else "NO -> H-lin FAILS, bounds of step S withdrawn"))
    if tag in AUG and lab in AUG[tag]:
        if tag == "l500_k16":
            r1b = abs(Rv["100x"]["dpco2"] / R0["100x"]["dpco2"] - 1)
            print("  (outside verdict) pre-M4 rule R1b: relative |dpCO2 100x| %.3f %% (<= 10 %%: robust) ; pre-M4 test: %+.2f point" % (100 * r1b, AUG[tag][lab]))
        else:
            print("  (outside verdict) pre-M4 rule: |dAF| <= 1 point ? %s ; pre-M4 test: %+.4f point" % (abs(dAF) <= 1.0, AUG[tag][lab]))
print("\nEND")
