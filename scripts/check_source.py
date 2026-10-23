#!/usr/bin/env python3
"""Check the archived source bytes; --syntax also parses code without executing the pipeline."""
import argparse
import csv
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys

SOURCE_DIRS = ('pipeline', 'config', 'ocim', 'annual_svr', 'assessment')
UNHASHED = {'config/paths.example.env'}


def check_manifest(root):
    root = Path(root).resolve()
    errors, declared = [], set()
    with (root / 'MANIFEST.tsv').open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if not rows:
        return ['Empty source manifest']
    for row in rows:
        rel = row.get('public_path', '')
        path = PurePosixPath(rel)
        if not rel or path.is_absolute() or '..' in path.parts or '\\' in rel or ':' in rel or path.as_posix() != rel:
            errors.append('Invalid source path: ' + rel)
            continue
        if rel in declared:
            errors.append('Duplicate source path: ' + rel)
        declared.add(rel)
        target = (root / rel).resolve()
        if not target.is_relative_to(root):
            errors.append('Source outside repository: ' + rel)
            continue
        if not target.is_file():
            errors.append('Missing source: ' + rel)
            continue
        expected = row.get('public_md5', '')
        if not re.fullmatch('[0-9a-f]{32}', expected) or hashlib.md5(target.read_bytes()).hexdigest() != expected:
            errors.append('Source fingerprint differs: ' + rel)
        if not re.fullmatch('[0-9a-f]{32}', row.get('grit_md5', '')):
            errors.append('Invalid original fingerprint: ' + rel)
    actual = {p.relative_to(root).as_posix() for name in SOURCE_DIRS
              for p in (root / name).rglob('*') if p.is_file()} - UNHASHED
    if actual != declared:
        errors.append('Source inventory differs: ' + repr(sorted(actual ^ declared)))
    return errors


def check_syntax(root):
    errors, counts = [], {'Python': 0, 'shell': 0}
    files = [p for name in (*SOURCE_DIRS, 'scripts', 'tests') for p in (root / name).rglob('*') if p.is_file()]
    bash = os.environ.get('CHECK_BASH') or shutil.which('bash')
    rscript = os.environ.get('CHECK_RSCRIPT') or shutil.which('Rscript')
    for p in sorted(files):
        if p.suffix == '.py':
            try:
                compile(p.read_bytes(), str(p), 'exec')
                counts['Python'] += 1
            except SyntaxError as e:
                errors.append(str(e))
        elif p.suffix in ('.sh', '.sbatch'):
            if not bash:
                errors.append('Bash not available')
                break
            result = subprocess.run([bash, '-n', str(p)], capture_output=True, text=True)
            if result.returncode:
                errors.append(str(p.relative_to(root)) + ': ' + result.stderr)
            else:
                counts['shell'] += 1
    if not rscript:
        errors.append('Rscript not available')
    else:
        result = subprocess.run([rscript, '--vanilla', str(root / 'scripts/check_r.R'), str(root)],
                                capture_output=True, text=True)
        if result.returncode:
            errors.append(result.stderr)
        else:
            print(result.stdout.strip())
    print('Syntax:', counts)
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--syntax', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    errors = check_manifest(root)
    if not errors and args.syntax:
        errors.extend(check_syntax(root))
    for error in errors:
        print('FAIL:', error)
    print('Source check:', 'FAIL' if errors else 'PASS')
    return bool(errors)


if __name__ == '__main__':
    sys.exit(main())
