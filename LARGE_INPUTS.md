# Large input not shipped in this package

One original input is deliberately excluded because GitHub enforces a 100 MB
per-file limit and a web upload of ~500 MB is unreliable:

    data/33067605-GCST90012074-EFO_0009419.h.tsv.gz
    Fenland CX3CL1 plasma pQTL, GWAS Catalog accession GCST90012074
    size: 501,219,328 bytes

Its SHA256 is recorded in `results/a1_input_sha256.tsv`, so the checksum chain
is not broken: anyone can re-download the public file from the EBI FTP / GWAS
Catalog and verify it against the recorded hash. The download helper is
`code/a1_b3_parallel_dl.py` (used by `bash code/run_all.sh --with-fetch`).

The offline reproduction chain does NOT read this file. It uses the deposited
region-level extract `results/a1_b3_fenland_region.tsv` and
`results/a1_b3_coloc_full.json`, so `bash code/run_all.sh` completes without it.
`run_all.sh` regenerates `results/a1_input_sha256.tsv` from the files actually
present; a lightweight copy therefore lists one fewer row than the analysis-time
manifest. The analysis-time manifest and `results/a1_version_snapshot.json`
are retained unchanged as provenance records.

Also excluded, because they are download scratch files rather than inputs:
`data/_parts_gcst90012074/` (parallel-download chunks), `data/*.INCOMPLETE`
(aborted downloads), `data/tbi/` (tabix indexes used only by network re-fetch).

If a journal requires the full original input to sit inside the Zenodo deposit,
add `data/33067605-GCST90012074-EFO_0009419.h.tsv.gz` as a second file on the
Zenodo record. Zenodo accepts up to 100 files / 50 GB per record, and files can
still be added within 45 days of publication.
