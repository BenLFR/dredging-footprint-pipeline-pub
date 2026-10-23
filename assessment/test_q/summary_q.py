# Test Q (rule written beforehand): airborne fraction, injection at the real bottom (M4)
# against column bottom (_r6). Derived from an earlier OCIM mass-balance check (not in this repository): same calculation of the
# balance, plus the four runs, the check of the control and the mean depth of injection.
# Order imposed: locks Q0a-Q0c are printed first; if one fails, no AF is printed.
# Execution: sbatch (python3 + h5py), never on the interactive node.
import os, glob, sys
import numpy as np, h5py

RHO   = 1025.0   # kg/m3: ns_step_eb.m (read in the OCIM code)
MATM  = 1.8e20   # mol: ns_step_eb.m, same
GPMOL = 12.0     # g/mol C
NSTEP = 321      # co2model_batch.m: loop i = 2:322, dt = 1 yr
PGC_PPM_LIT = 2.124e15  # gC/ppm, literature constant, for comparison only
TOL_BALANCE = 0.1      # %: OCIM mass-balance rule, any gap above ~0.1 % fails
TOL_SOURCE = 1e-9    # relative: design choice, no source (Q0b)
TOL_CONTROL = 1e-12   # ppm: design choice, no source (Q0c)
THRESHOLD_MAT = 1.0      # AF point: design choice, no source

RUN = os.environ["DFP_SCRATCH"] + "/output_V6_sala2021_r6_flux10"
LABELS = ('smooth', 'flag')
MODES = (('r6', ''), ('m4', '_m4'))
AMPS = ('1x', '10x', '100x')

def T(a):
    return np.array(a).T   # h5py returns the axes reversed with respect to MATLAB

def only_one(pat):
    f = sorted(glob.glob(pat))
    if len(f) != 1:
        sys.exit(f"ABORT: {len(f)} files for {pat}")
    return f[0]

def last_row(f, name, m):
    d = f[name]
    if d.shape[-1] == m:
        return np.asarray(d[-1, :], dtype=float)
    if d.shape[0] == m:
        return np.asarray(d[:, -1], dtype=float)
    sys.exit(f"ABORT: {name} shape {d.shape} incompatible with m={m}")

def geom(f):
    M3d = T(f['M3d']).ravel(order='F')
    g = f['grid']
    V3 = (T(g['DXT3d']) * T(g['DYT3d']) * T(g['DZT3d'])).ravel(order='F')
    Z3 = T(g['ZT3d']).ravel(order='F')
    iocn = np.flatnonzero(M3d == 1)
    return V3[iocn], Z3[iocn]

def measure(base, pert):
    with h5py.File(base, 'r') as fb, h5py.File(pert, 'r') as fp:
        V, Z = geom(fp)
        m = V.size
        J = np.asarray(fp['Jdredge_forcing'], dtype=float).ravel()
        if J.size != m:
            sys.exit(f"ABORT: Jdredge_forcing {J.size} != m {m}")
        Mw = V * RHO
        w = J * Mw
        q = float(np.sum(w) * GPMOL * 1e-6)          # gC/yr, trend actually applied
        C_inj = q * NSTEP
        pa_b = float(np.asarray(fb['pco2a']).ravel()[-1])
        pa_p = float(np.asarray(fp['pco2a']).ravel()[-1])
        dC_atm = (pa_p - pa_b) * 1e-6 * MATM * GPMOL
        dC_oce = float(np.sum((last_row(fp, 'DIC', m) - last_row(fb, 'DIC', m)) * Mw) * GPMOL * 1e-6)
        yr = np.asarray(fp['year']).ravel()
        act = J != 0
        D = float(np.sum(w[act] * Z[act]) / np.sum(w[act]))   # m, ZT3d weighted by mass
        return dict(q=q, C_inj=C_inj, pa_b=pa_b, dpco2=pa_p - pa_b, dC_atm=dC_atm, dC_oce=dC_oce,
                    gap=(dC_atm + dC_oce - C_inj) / C_inj * 100.0, n_act=int(act.sum()), D=D,
                    yr0=float(yr[0]), yr1=float(yr[-1]),
                    AF=dC_atm / C_inj * 100.0, AF_lit=(pa_p - pa_b) * PGC_PPM_LIT / C_inj * 100.0)

R = {}
for lab in LABELS:
    for mode, suf in MODES:
        d = f"{RUN}/ocim_{lab}{suf}"
        base = only_one(f"{d}/ocim_baseline_*.mat")
        for amp in AMPS:
            R[(lab, mode, amp)] = measure(base, only_one(f"{d}/ocim_dredge_{amp}_*.mat"))
            R[(lab, mode, amp)]['file_order'] = only_one(f"{d}/ocim_dredge_{amp}_*.mat")

print("== LOCKS (read first) ==")
ok = True
for (lab, mode, amp), r in sorted(R.items()):
    if amp == '1x':
        p = abs(r['gap']) <= TOL_BALANCE
        ok &= p
        print(f"Q0a balance {lab:6s} {mode}: gap {r['gap']:+.6f} %  {'PASS' if p else 'FAIL'}")
for lab in LABELS:
    x = R[(lab, 'm4', '1x')]['q'] / R[(lab, 'r6', '1x')]['q'] - 1
    p = abs(x) <= TOL_SOURCE
    ok &= p
    print(f"Q0b source {lab:6s}: q_M4/q_r6 - 1 = {x:+.3e}  {'PASS' if p else 'FAIL'}")
pab = [r['pa_b'] for r in R.values()]
x = max(pab) - min(pab)
p = x <= TOL_CONTROL
ok &= p
print(f"Q0c control: range of final baseline pCO2a = {x:.3e} ppm  {'PASS' if p else 'FAIL'}")
if not ok:
    sys.exit("LOCK FAILED: no airborne fraction printed (rule written beforehand).")

print("\n== RAW RESULT ==")
print("label\tmode\tamp\tq_TgC_yr\tC_inj_TgC\tdpco2_ppm\tn_actives\tD_m\tyears\tAF_pct\tAF_2124_pct")
for (lab, mode, amp), r in sorted(R.items()):
    print(f"{lab}\t{mode}\t{amp}\t{r['q']/1e12:.6f}\t{r['C_inj']/1e12:.6f}\t{r['dpco2']:.9e}\t{r['n_act']}\t"
          f"{r['D']:.1f}\t{r['yr0']:.0f}-{r['yr1']:.0f}\t{r['AF']:.4f}\t{r['AF_lit']:.4f}")
print("\nlinearite (AF_amp / AF_1x - 1):")
for lab in LABELS:
    for mode, _ in MODES:
        a1 = R[(lab, mode, '1x')]['AF']
        print(f"  {lab} {mode}: 10x {R[(lab, mode, '10x')]['AF']/a1-1:+.2e}  100x {R[(lab, mode, '100x')]['AF']/a1-1:+.2e}")

print("\n== AGAINST THE RULE (1x) ==")
for lab in LABELS:
    r6, m4 = R[(lab, 'r6', '1x')], R[(lab, 'm4', '1x')]
    dAF, dD = m4['AF'] - r6['AF'], m4['D'] - r6['D']
    print(f"{lab}: D_r6 {r6['D']:.1f} m -> D_M4 {m4['D']:.1f} m ; AF_r6 {r6['AF']:.4f} % -> AF_M4 {m4['AF']:.4f} % ; "
          f"dAF {dAF:+.4f} pt")
