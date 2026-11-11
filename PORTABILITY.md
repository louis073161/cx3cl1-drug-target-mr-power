# Portability note

The 14-step reproduction chain was developed on a single Windows workstation,
and several scripts in `code/` carry an absolute `ROOT` constant pointing at
that machine. This does not affect any published number, but it does mean a
fresh copy must be told where it lives before the first run.

Do this once, from the root of the extracted package:

    python code/make_portable.py
    bash code/run_all.sh

`make_portable.py` only replaces two absolute directory strings in `code/*.py`
with the location of your copy. It is idempotent: running it again changes
nothing. It does not read or write any result value.

After this adaptation, the SHA256 hashes of the adapted `code/*.py` files in
your LOCAL copy differ from `results/a1_input_sha256.tsv`, which records the
unmodified files as deposited. This is expected: it is a path rewrite, not a
change to the analysis. To verify the deposited originals, check the unmodified
files in this ZIP before running the adapter, or re-download this archive.

The Word-layout scripts (`code/wb_format_*.py`) belong to the internal Chinese
report workflow and are not called by `run_all.sh`; they are kept only for
provenance.

`bash code/run_all.sh --with-fetch` additionally re-downloads public summary
statistics from FinnGen GCS, EBI FTP and Ensembl REST; the default offline run
needs no network at all.
