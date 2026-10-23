# -*- coding: utf-8 -*-
"""Out-of-verdict check: the 1x of the rerun of the OCIM resteB run must equal the 1x
of the incomplete folder, the model being deterministic. Read only.
Compares each numeric variable of the two .mat files (v7 via scipy, v7.3 via h5py): max |gap| and max |value|."""
import os, sys
import numpy as np

R = os.environ["DFP_SCRATCH"] + "/hlin_20260928/"
FILES = {
    "1x": (R + "ocim_smooth_resteB/ocim_dredge_1x_20261002_070118.mat",
           R + "ocim_smooth_resteB_incomplet_20260929/ocim_dredge_1x_20260928_192159.mat"),
    "baseline": (R + "ocim_smooth_resteB/ocim_baseline_20261002_070118.mat",
                 R + "ocim_smooth_resteB_incomplet_20260929/ocim_baseline_20260928_192159.mat"),
    "10x": (R + "ocim_smooth_resteB/ocim_dredge_10x_20261002_070118.mat",
            R + "ocim_smooth_resteB_incomplet_20260929/ocim_dredge_10x_20260928_192159.mat"),
}


def load(p):
    try:
        from scipy.io import loadmat
        d = loadmat(p)
        return {k: v for k, v in d.items() if not k.startswith("__")}
    except NotImplementedError:
        import h5py
        out = {}
        with h5py.File(p, "r") as f:
            def visit(item_name, obj):
                if isinstance(obj, h5py.Dataset):
                    out[item_name] = obj[()]
            f.visititems(visit)
        return out


for item_name, (a, b) in FILES.items():
    da, db = load(a), load(b)
    print("==", item_name, ":", len(da), "variables against", len(db))
    if set(da) != set(db):
        print("   differing names:", sorted(set(da) ^ set(db)))
    worst = 0.0
    for k in sorted(set(da) & set(db)):
        x, y = np.asarray(da[k]), np.asarray(db[k])
        if x.dtype.kind not in "fiuc" or y.dtype.kind not in "fiuc":
            continue
        if x.shape != y.shape:
            print("   %s: shapes %s against %s" % (k, x.shape, y.shape))
            continue
        e = float(np.nanmax(np.abs(x - y))) if x.size else 0.0
        m = float(np.nanmax(np.abs(x))) if x.size else 0.0
        if e > 0:
            print("   %s: max|gap| %.3e ; max|value| %.3e" % (k, e, m))
        worst = max(worst, e / m if m else e)
    print("   maximum relative gap over all variables:", "%.3e" % worst)
sys.stdout.flush()
