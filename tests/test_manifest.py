"""The source manifest rejects altered, missing and undeclared files."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from check_source import check_manifest


class ManifestTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'pipeline').mkdir()
        self.source = self.root / 'pipeline/example.py'
        self.source.write_bytes(b'x = 1\n')
        digest = hashlib.md5(self.source.read_bytes()).hexdigest()
        self.header = 'public_path\tgrit_path\tgrit_md5\tpublic_md5\tprovenance\n'
        self.row = f'pipeline/example.py\toriginal.py\t{digest}\t{digest}\tGRIT copy\n'
        self.manifest = self.root / 'MANIFEST.tsv'
        self.manifest.write_text(self.header + self.row, encoding='utf-8')

    def test_exact_copy_and_mutation(self):
        self.assertEqual(check_manifest(self.root), [])
        self.source.write_bytes(b'x = 2\n')
        self.assertTrue(check_manifest(self.root))

    def test_missing_extra_and_duplicate(self):
        self.source.unlink()
        self.assertTrue(check_manifest(self.root))
        self.source.write_bytes(b'x = 1\n')
        (self.root / 'pipeline/extra.py').write_bytes(b'x = 3\n')
        self.assertTrue(check_manifest(self.root))
        (self.root / 'pipeline/extra.py').unlink()
        self.manifest.write_text(self.header + self.row + self.row, encoding='utf-8')
        self.assertTrue(check_manifest(self.root))

    def test_path_outside_repository(self):
        self.manifest.write_text(self.header + self.row.replace('pipeline/example.py', '../outside.py'), encoding='utf-8')
        self.assertTrue(check_manifest(self.root))


if __name__ == '__main__':
    unittest.main()
