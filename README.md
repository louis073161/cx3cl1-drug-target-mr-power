# A1 - Power and phenotype-definition heterogeneity in drug-target Mendelian randomization

Code, derived tables, vector figures and reporting materials for a retrospective
methodological study of circulating CX3CL1 and diabetic neuropathy across five
FinnGen releases (R9-R13), with cross-target replication, technical positive
controls and region-wide colocalisation.

## Identifiers and cross-links

| Artefact | Identifier | Role |
| --- | --- | --- |
| Retrospective OSF registration | https://osf.io/apwe2 ; DOI 10.17605/OSF.IO/APWE2 | Frozen analysis plan, Secondary Data Preregistration, deposited 11 November 2026, CC0 1.0 Universal |
| Zenodo archived snapshot | https://doi.org/10.5281/zenodo.24137602 ; DOI 10.5281/zenodo.24137602 | Immutable, versioned copy of this package |
| GitHub repository | https://github.com/louis073161/cx3cl1-drug-target-mr-power | Living code repository |

The OSF registration (10.17605/OSF.IO/APWE2, https://osf.io/apwe2) is frozen and is not modified by
this release. The Zenodo record links to the OSF registration and to this
GitHub repository (https://github.com/louis073161/cx3cl1-drug-target-mr-power); the manuscript, cover letter and declarations
cite all three identifiers.

## Reproduce offline (no network required)

On a computer other than the original analysis workstation, run
`python code/make_portable.py` once before the first reproduction. It rewrites
only the hard-coded directory constants inside `code/*.py` to the location of
your extracted copy; no statistic, figure or result value is changed. See
`PORTABILITY.md`. Then:

    bash code/run_all.sh

This runs the 14-step offline reproduction chain and recomputes every reported
statistic and figure from the deposited result tables. Network-dependent
re-fetching of public summary statistics is available through
`code/a1_multi_fetch.py`, `code/a1_tabix_fetch.py`, `code/a1_b1b2_multitarget.py`
and the B3 scripts; it is not required to verify the reported numbers.
`bash code/run_all.sh --with-fetch` runs the network-dependent steps as well.

## Layout

- `code/` - all analysis, figure and document-generation scripts; `run_all.sh`
  is the one-command reproduction entry point; `requirements.lock.txt` pins
  the Python dependencies.
- `results/` - result tables (TSV), JSON summaries, vector figures (PDF/SVG)
  and PNG previews, the input SHA256 manifest and the version snapshot.
- `data/` - the small public input snapshots and tabix indexes that are read by
  the offline chain.
- `results/a1_endpoint_dict.tsv` - one row per release-endpoint combination,
  with endpoint codes, public case/control counts and download URLs.
- `results/a1_input_sha256.tsv` - SHA256 of every input file recorded at
  analysis time. `run_all.sh` regenerates this manifest from the files actually
  present in a given copy, so a lightweight copy lists one fewer entry; this is
  expected and is explained in `LARGE_INPUTS.md`.

## Data sources

Open Targets Platform v26.9 (UKB-PPP CX3CL1 pQTL, study
UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1, n = 33,436), FinnGen public releases
R9-R13, the Fenland CX3CL1 plasma pQTL (GWAS Catalog GCST90012074, used for
region-wide colocalisation), and Ensembl REST 1000 Genomes phase 3 EUR LD.
All inputs are public; no individual-level data and no access tokens are
included. FinnGen and UKB-PPP participants are non-overlapping.

## Reporting restriction

The manuscript states "no supportive evidence; the available data are
insufficient to test an effect of moderate magnitude". It never states that
CX3CL1 has no effect on diabetic neuropathy, and it does not revise the
tissue-specific experimental hypothesis of the parent grant. Post-hoc power,
minimum detectable effects and the four-class reporting rule are reported for
every release-endpoint combination rather than only for the primary one.

## Licence

CC0 1.0 Universal (SPDX: CC0-1.0), matching the OSF registration. Third-party
summary statistics remain subject to their original data-use terms. See
`LICENSE.md`.
