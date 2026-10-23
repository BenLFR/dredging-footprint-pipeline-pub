#!/usr/bin/env python3
"""Inventory only. No R, no submission, no numerical sensitivity calculation."""
import os, argparse
import csv
import hashlib
import io
from pathlib import Path
import re
import sys

# Paths and line ranges from GRIT; expected public hashes: docs/CHECKSUMS.md.
ROOT = Path(os.environ["DFP_REPO"]) / "assessment/omega_pd"
PIPE = Path(os.environ["DFP_REPO"]) / "pipeline"
CENTRAL = Path(os.environ["DFP_SCRATCH"] + "/output_V6_sala2021_r6_flux10")
ARCHIVE = Path(os.environ["DFP_SCRATCH"] + "/output_V6_incert")
OUT = ARCHIVE / "omega_pd"
PINS = {
    "merge": (PIPE / "step5_merge_tiles_corr.R", "e40847579fa81eacadd1955e1013b42b", "135-154;296-320;495-605"),
    "step6": (PIPE / "step6_calculate_cri_corr.R", "b54cc42ebcc1d57e83aad28f42d82100", "97-160"),
    "yaml": (Path(os.environ["DFP_REPO"]) / "config/fi_parameters_with_freshness.yaml", "8619588df22088b94961f457437a5875", "1-51"),
    "launcher_merge": (PIPE / "merge_rerun.sbatch", "4dbc5b1d4229653698a83dcb7d5998dd", "48-51"),
    "launcher_step6": (PIPE / "step6_pr.sbatch", "fa0bacd5938d5ee9f809a6bef75ab7d3", "22-27"),
    "manifest_generator": (PIPE / "make_manifest.sh", "450594c105dbe26adcfe5fb1a33b6332", "whole file"),
    "pipeline_guard": (PIPE / "check_pipeline.sh", "0f520cc5694a16432e23a0dc784f6b66", "whole file"),
    "omega_prepare": (Path(os.environ["DFP_REPO"]) / "assessment/prep_cases.sh", "2fed504839e01e90176461a429807cfb", "14-21"),
    "omega_submit": (Path(os.environ["DFP_REPO"]) / "assessment/launch.sh", "c5c871f6d96d30ba28e14c19440daf8a", "27-37"),
    "omega_read": (Path(os.environ["DFP_REPO"]) / "assessment/read.R", "d2e4dedd1dbdeced55cb7f5dc7c2490f", "13-42"),
}
FILTER = re.compile("token|secret|key|password|client_id|auth", re.I)

def md5(path):
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""): # I/O block size only.
            digest.update(block)
    return digest.hexdigest()

def sole(directory, pattern):
    matches = sorted(directory.glob(pattern))
    # Exactly one file: design guard, not a scientific threshold.
    if len(matches) != 1:
        raise RuntimeError(f"Ambiguous/missing archive: {directory}/{pattern}")
    return matches[0]

def inventory():
    rows = []
    def add(name, path, expected=None, lines="binary; rows identified by grid_id"):
        if not path.is_file():
            raise RuntimeError(f"Missing input: {path}")
        digest = md5(path)
        if expected and digest != expected:
            raise RuntimeError(f"Fingerprint changed: {name}")
        rows.append(dict(id=name, path=str(path), md5=digest, bytes=path.stat().st_size, lines=lines))
    for name, (path, digest, lines) in PINS.items():
        add(name, path, digest, lines)
    for label in ("smooth", "flag"):
        # Exact central timestamps of the reference outputs.
        add(f"central_fi_{label}", CENTRAL / f"merge_{label}_pr/fi_grid_20260925_124425.parquet")
        add(f"central_cri_{label}", CENTRAL / f"cri_pr/cri_final_{label}_20260925_125232.parquet")
        add(f"frac_{label}", Path(os.environ["DFP_DATA"]) / f"output_V6_rejeu/decharge/fraction_sur_place_{label}.csv", lines="header and all data rows; grid_id")
        # Existing Omega scenarios, YAML in service.
        for case in ("ref", "om0.1", "om1.0", "om5.0"):
            add(f"{case}_fi_{label}", sole(ARCHIVE / case / f"merge_{label}_pr", "fi_grid_*.parquet"))
            add(f"{case}_cri_{label}", sole(ARCHIVE / case / "cri_pr", f"cri_final_{label}_*.parquet"))
    return rows

def tsv(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Freeze inputs.tsv, exclusively; no calculation")
    parser.add_argument("--check", action="store_true", help="Check every frozen input; no writes")
    args = parser.parse_args()
    if args.write and args.check:
        raise RuntimeError("Select write or check")
    if ROOT.resolve() != Path(__file__).resolve().parent:
        raise RuntimeError("Run only from the designated GRIT script directory")
    if OUT.is_symlink() or (OUT.exists() and any(OUT.iterdir())):
        raise RuntimeError("Output directory is a symlink or is not empty; no overwrite")
    rows = inventory()
    target = ROOT / "inputs.tsv"
    content = tsv(rows)
    if args.write:
        with target.open("x", encoding="utf-8", newline="") as stream:
            stream.write(content)
        print("Input manifest frozen: " + str(target))
    elif args.check:
        if target.read_text(encoding="utf-8") != content:
            raise RuntimeError("Input manifest changed; stop and review")
        print("Frozen input fingerprints unchanged")
    else:
        print("DRY preparation: would freeze " + str(target))
    for row in rows:
        print(row["md5"] + "  " + row["path"])

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        message = str(exc)
        print("ABORT: " + ("filtered diagnostic" if FILTER.search(message) else message), file=sys.stderr)
        sys.exit(1)
