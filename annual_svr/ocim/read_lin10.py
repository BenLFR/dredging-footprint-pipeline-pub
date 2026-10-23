# read_lin10.py -- linearity at 10x of the OCIM run under the annual cap (rule written beforehand).
# read_svr.py computes the 1x, 10x and 100x amplitudes but prints only the 100x (L4).
# This script reuses ITS functions (measure, run), without copying them: it executes the part of read_svr.py that precedes
# its reading loop, after checking its fingerprint. No OCIM run is relaunched; the .mat files are only read.
# ORDER IMPOSED: lock first. V-lin: AF 1x and |AF_100x/AF_1x - 1| re-read = those printed by read_svr.py (54.4074 / 49.0434 % to
#   1e-4 point; 3.42e-5 / 8.63e-5 at the printed rounding). If it fails, nothing else is printed for the label.
# Substance: |AF_10x/AF_1x - 1| and |AF_100x/AF_1x - 1| <= 5e-3 (L4 threshold); outside verdict: dpCO2_10x/dpCO2_1x and
#   dpCO2_100x/dpCO2_1x (what the text calls "factors of 10.00 and 100.00").
# Execution: sbatch (python3 + h5py), never on the interactive node.
import hashlib, sys

SRC = "read_svr.py"
if hashlib.md5(open(SRC, "rb").read()).hexdigest() != "02239201d60dddf62a38f600c6eae7bf":
    sys.exit("ABORT: read_svr.py does not have the fingerprint 02239201")
code = open(SRC, encoding="utf-8").read()
cut = code.index('for lab in ("smooth", "flag"):')
exec(compile(code[:cut], SRC, "exec"))       # defines RUN, BASE, AMPS, measure, run, constants

AF_285 = {"smooth": 54.4074, "flag": 49.0434}   # values printed by read_svr.py
L4_285 = {"smooth": 3.42e-5, "flag": 8.63e-5}   # lock L4 of read_svr.py (printed to 3 digits)
for lab in ("smooth", "flag"):
    print("\n==================== %s ====================" % lab)
    Rv = run("%s/ocim_%s_an" % (BASE, lab))
    af1, af10, af100 = (Rv[a]["AF"] for a in AMPS)
    l4 = abs(af100 / af1 - 1)
    v = abs(af1 - AF_285[lab]) <= 1e-4 and float("%.2e" % l4) == L4_285[lab]
    print("  V-lin: AF 1x %.6f %% (earlier reference: %.4f) ; |AF100/AF1 - 1| %.2e (earlier reference: %.2e)  %s"
          % (af1, AF_285[lab], l4, L4_285[lab], "PASS" if v else "FAIL"))
    if not v:
        print("  >>> V-lin FAILED: the reading does not reproduce read_svr.py; nothing is read for %s." % lab)
        continue
    l10 = abs(af10 / af1 - 1)
    print("  AF 1x %.6f %% ; 10x %.6f %% ; 100x %.6f %%" % (af1, af10, af100))
    print("  |AF10/AF1 - 1| = %.2e  %s" % (l10, "PASS" if l10 <= 5e-3 else "FAIL"))
    print("  |AF100/AF1 - 1| = %.2e  %s" % (l4, "PASS" if l4 <= 5e-3 else "FAIL"))
    d1 = Rv["1x"]["dpco2"]
    print("  (not part of the verdict) dpCO2: 1x %.6e ; 10x %.6e ; 100x %.6e ppm" % (d1, Rv["10x"]["dpco2"], Rv["100x"]["dpco2"]))
    print("  (not part of the verdict) dpCO2_10x / dpCO2_1x = %.6f ; dpCO2_100x / dpCO2_1x = %.6f"
          % (Rv["10x"]["dpco2"] / d1, Rv["100x"]["dpco2"] / d1))
    print("  (not part of the verdict) balance 10x %.6f %% ; 100x %.6f %%" % (Rv["10x"]["gap"], Rv["100x"]["gap"]))
print("\nEND")
