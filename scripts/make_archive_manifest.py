"""Regenerate the source-archive inventory; --check compares it without writing."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def generate(root):
    metadata = json.loads((root / 'codemeta.json').read_text(encoding='utf-8'))
    tracked = subprocess.check_output(['git', '-C', str(root), 'ls-files', '-z']).decode('utf-8').split('\0')
    files = sorted(p for p in tracked if p and p != 'ARCHIVE_MANIFEST.yaml')
    rows = []
    for rel in files:
        path = (root / rel).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError('Archive path escapes repository: ' + rel)
        content = path.read_bytes()
        rows.append({'path': rel, 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()})
    data = {'schema_version': '1.0', 'generated_on': metadata['dateModified'],
            'repository': metadata['codeRepository'],
            'release': {'version': metadata['version'], 'status': 'unreleased', 'doi': None},
            'scope': 'Tracked source files; no study data or model matrices. This manifest excludes itself.',
            'source_provenance': 'MANIFEST.tsv', 'third_party_inputs': 'THIRD_PARTY_NOTICES.md',
            'artifacts': rows}
    # JSON is valid YAML 1.2; using it avoids an extra YAML dependency in source-only CI.
    return (json.dumps(data, ensure_ascii=False, indent=2) + '\n').encode('utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = root / 'ARCHIVE_MANIFEST.yaml'
    expected = generate(root)
    if args.check:
        if path.read_bytes() != expected:
            print('Archive manifest: FAIL; regenerate and review the changed files')
            return 1
        print('Archive manifest: PASS')
    else:
        path.write_bytes(expected)
        print('Archive manifest written; review before committing')
    return 0


if __name__ == '__main__':
    sys.exit(main())
