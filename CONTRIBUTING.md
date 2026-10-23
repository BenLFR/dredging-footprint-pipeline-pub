# Contributing

This is a research source archive with a reconstructed public layout. Read [README.md](README.md),
[INPUTS.md](docs/INPUTS.md) and [CHECKSUMS.md](docs/CHECKSUMS.md) before changing a script.

Run the source checks with Python, R and Bash installed:

```sh
python -m unittest discover -s tests -v
python scripts/check_source.py --syntax
python scripts/make_archive_manifest.py --check
```

These checks compare source bytes and parse syntax; they do not run the scientific analysis. Do not commit
AIS data, model matrices, credentials, local path settings or generated outputs. Use `config/paths.example.env`
to describe required paths, and keep actual values in the ignored `paths.env`.

For a source change, explain its scientific effect, the affected input/output contract and its validation.
Update the published checksum in `MANIFEST.tsv` only after reviewing the change. Keep the original GRIT
fingerprint as provenance; do not present changed code as the exact historical execution. Review checksum
guards that refer to changed source files, and preserve job-approval checks.

Record exact input versions, dependency versions, commands and output checksums when performing a new run.
Do not describe a syntax check as a reproduction of a reported result.

After reviewing changes, stage the intended file inventory and regenerate `ARCHIVE_MANIFEST.yaml` with
`python scripts/make_archive_manifest.py`, then stage the manifest and run its `--check` mode. This inventory
covers tracked source and documentation, not an execution archive or a set of scientific outputs.
