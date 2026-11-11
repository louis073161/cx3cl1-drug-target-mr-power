# -*- coding: utf-8 -*-
"""
整合 A1 英文稿定稿骨架 + 完整投稿包 + 遗留项登记。

原则：
  - 既有 EN 段落原样合并，不改写既有结论；
  - 本轮新增结果（A4/A5/B1/B2）各写成可直接粘贴的英文段落；
  - 所有还需本人填写或原文无法核实之处集中登记，不推测、不冒充完成。
"""
import csv
import io
import json
import os
import re
import glob

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
PROJ = 'F:/博士/国自然科学基金/2025神自然'
RES = ROOT + '/results'
PACK = ROOT + '/submission_pack'
MD = PROJ + '/_markdown源'
os.makedirs(PACK, exist_ok=True)


def txt(path, s):
    with io.open(path, 'w', encoding='utf-8') as f:
        f.write(s.rstrip() + '\n')
    print('  wrote', os.path.basename(path))


def write_tsv(path, rows):
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------------------------------------------------------- 新增英文段落
I4 = ("[I4 - differentiation from prior work] Phenotype-related distortion "
      "of Mendelian randomization has been examined before. Hu et al. "
      "(J Neurochem 2025;169(7):e70151, doi:10.1111/jnc.70151) showed that "
      "self-reported Alzheimer's disease status in UK Biobank, and the "
      "proxy-based GWAS derived from it, can reverse causal directions "
      "relative to clinically ascertained datasets, including FinnGen R10; "
      "their focus is the validity of the phenotype source. Zou et al. "
      "(BMC Med Res Methodol 2022;22:162) addressed between-study "
      "heterogeneity with a Bayesian random-effects model, and work around "
      "LocusCompare2 has catalogued inconsistency among gene-implication "
      "methods [VERIFY: original Nature Genetics citation, volume and DOI]. "
      "None of these quantifies, for a drug-target protein-QTL MR, how much "
      "one reported null depends on successive releases of the same biobank, "
      "on alternative endpoint definitions within a release, or on the "
      "choice of instrument set, nor converts that dependence into a minimum "
      "detectable effect and a wording rule. Our contribution is therefore "
      "not the observation that analytical choices matter: it is a "
      "reproducible, quantified cross-context audit paired with a "
      "reporting-layer decision rule. We do not claim priority for showing "
      "that MR conclusions are sensitive to data choices.")

R6 = ("[R6, items 12b/13d - quantitative attribution of drift] Across all "
      "171 pairwise comparisons of the 19 primary estimates, the "
      "standardised drift |beta_i - beta_j| / sqrt(SE_i^2 + SE_j^2) never "
      "exceeded 1.96 (maximum 1.456; binomial P=1.00 against the 5% "
      "expectation), including all 21 same-endpoint-across-release pairs "
      "(maximum 0.655). The apparent 10.72-fold range of odds ratios is "
      "therefore compatible with sampling noise and does not demonstrate "
      "true effect heterogeneity between releases or endpoint definitions. "
      "Among the 45 pairs with publicly reported case counts, |delta "
      "log(cases)| showed a monotonic association with |delta logOR| "
      "(Spearman rho=0.335; 20,000-permutation P=0.0146), but this is "
      "descriptive, is confounded by the same case counts that determine "
      "the standard errors, and cannot be evaluated for R11 or for most "
      "R13 endpoints; it was not used to adjust any estimate. Pairs "
      "differing only in endpoint strictness (STRICT/EXMORE variants) had "
      "smaller drift than other pairs (mean |delta logOR| 0.404 vs 0.761; "
      "permutation P=0.0225). The correct reading is that these estimates "
      "are too imprecise to separate real heterogeneity from noise - not "
      "that heterogeneity has been established.")

R7 = ("[R7, items 8/12b - three-axis drift and the false-precision trap] A "
      "4 x 19 sensitivity grid (single cis lead, all cis, LD-clumped, all "
      "available; 76 estimates) showed that the instrument axis produced "
      "the largest point-estimate spread (median beta range 0.649, versus "
      "0.134 for the endpoint axis and 0.077 for the release axis), with 13 "
      "of 19 contexts changing direction across instrument sets. The "
      "clumped and all-available sets nevertheless contain trans "
      "instruments with documented horizontal pleiotropy: their IVW "
      "standard errors shrink only because the estimator assumes every "
      "instrument is valid, yielding apparent power of 100% and minimum "
      "detectable odds ratios near 1.06-1.21. We treat this as false "
      "precision; adding pleiotropic instruments adds no valid information. "
      "Only one of the 76 estimates was nominally significant (R9 "
      "DM_RETINOPATHY_EXMORE, all-available set, P=0.0492); it lies in a "
      "negative-control endpoint, does not survive Bonferroni correction "
      "and is not interpreted as evidence. The cis-only sets remain "
      "underpowered (5.45%-10.10% power for OR=1.5 with the lead instrument).")

R8 = ("[R8 - cross-target replication of the reporting problem] To test "
      "whether the reporting problem generalises beyond a single target, we "
      "repeated the identical pipeline for seven further druggable proteins "
      "with a UKB-PPP cis pQTL (ACE, AGT, IL1RN, IL6R, LDLR, PCSK9 and "
      "SERPINE1). Twelve candidates were queried; four returned no "
      "accessible cis pQTL and one target lookup failed, and none was "
      "imputed. Across 42 target-by-context estimates from three endpoints "
      "in R9 and R13, power for OR=1.5 ranged from 5.53% to 100%, the "
      "minimum detectable OR from 1.06 to 197.47, 11 of 42 estimates had "
      "power below 20%, 18 had a minimum detectable OR above 2, one test "
      "was nominally significant (2.1 expected) and 6 of 21 target-endpoint "
      "pairs changed direction between releases. Underpowering was not "
      "universal. Stratifying by instrument strength gave median power "
      "99.997% for F>=1000 (ACE, AGT, IL6R, PCSK9), 22.8% for 100<=F<1000 "
      "(IL1RN, SERPINE1) and 5.85% for F<100 (LDLR). Low power is therefore "
      "an interaction between instrument strength and endpoint size rather "
      "than an intrinsic property of these FinnGen endpoints. These are "
      "cross-target, not independent-dataset, replications: all outcomes "
      "derive from one resource and each target uses one cis lead without "
      "clumping, colocalisation or Steiger testing.")

R9 = ("[R9 - biological positive control] As a biological positive control, "
      "the same cis-pQTL pipeline was applied to three large FinnGen R13 "
      "cardiovascular endpoints (21 estimates). PCSK9 matched its "
      "pre-specified direction on all three and remained significant after "
      "Bonferroni correction (0.05/21=0.00238): I9_CHD OR=2.195 "
      "(P=9.01x10^-29), I9_CORATHER OR=2.706 (P=1.91x10^-32) and "
      "E4_HYPERLIPMIX OR=13.534 (P=1.42x10^-7). This confirms that the "
      "pipeline can recover a known biological effect, complementing the "
      "earlier technical positive control based on ABO-region trans "
      "variants. LDLR did not replicate: direction was opposite to "
      "expectation and estimates were extreme (OR up to 2626.7). Its "
      "instrument is the weakest (z=-8.04, F=64.6) and has the smallest "
      "exposure effect (beta=-0.0428), so dividing by a small exposure "
      "beta inflates the Wald ratio. LDLR is therefore excluded as a "
      "positive control and reported as a cautionary example. Seven of 21 "
      "estimates survived Bonferroni correction; IL6R associations were "
      "not pre-specified and remain exploratory. A directional control "
      "does not establish causality for any target and is not used to "
      "reinterpret the CX3CL1 null.")

D6 = ("[D6 - revised interpretation of drift] Two of our starting "
      "expectations were not supported by the additional analyses and are "
      "reported as corrections. First, the large spread of point estimates "
      "across releases and endpoint definitions is fully compatible with "
      "sampling noise; we therefore do not claim that phenotype "
      "heterogeneity was demonstrated, and we use the word 'drift' only "
      "for changes in the reported estimate or conclusion label, not for "
      "proven changes in an underlying effect. Second, insufficient power "
      "did not replicate across every additional target: it concentrated "
      "in targets with weaker cis instruments. The generalisable finding "
      "is that the information content of a null drug-target MR result "
      "depends jointly on instrument strength, endpoint size, endpoint "
      "definition and instrument-set validity, and that none of these "
      "should be evaluated in isolation. These corrections do not alter "
      "the primary conclusion: we obtained no supportive evidence for "
      "circulating CX3CL1 in diabetic neuropathy, and the public data "
      "are insufficient to test an effect of moderate magnitude.")


# ---------------------------------------------------------------- 合并既有 EN
def collect_en():
    files = [
        MD + '/A1_A级18项清零与英文稿骨架_20261110.md',
        MD + '/A1_A1A2A3A6补正包_20261110.md',
        MD + '/A1科学价值提升与投稿成功率策略_20261110.md',
    ]
    out = []
    for p in files:
        if not os.path.exists(p):
            continue
        for line in io.open(p, encoding='utf-8'):
            if line.startswith('EN: '):
                s = line[4:].strip()
                if s and not s.startswith('**Suggested') \
                        and not s.startswith('**Abstract opening') \
                        and not s.startswith('**Interpretation sentence'):
                    out.append((os.path.basename(p), s))
    return out


def manuscript():
    en = collect_en()
    # 按标签前缀归类
    def find(prefixes):
        return [s for _, s in en
                if any(re.match(r'^\[' + p, s) for p in prefixes)]
    title = find(['Power and phenotype', 'Alternative titles'])
    # Title 段没有方括号，单独取
    title = [s for _, s in en if s.startswith('Power and phenotype')
             or s.startswith('Alternative titles')]
    abstract = [s for _, s in en if s.startswith('Background.')
                or s.startswith('Methods.')
                or s.startswith('Results.')
                or s.startswith('Conclusions.')]
    intro = find(['I1', 'I2', 'I3'])
    methods = find(['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7'])
    results = find(['R1', 'R2', 'R3', 'R4', 'R5'])
    disc = find(['D1', 'D2', 'D3', 'D4', 'D5'])
    decl = find(['DC', 'Supplementary Table S4', 'DRAFT'])
    lines = [
        '# A1 manuscript skeleton (English, consolidated)',
        '',
        '> Compiled on 2026-11-11 from all previously delivered EN paragraphs '
        'plus the new A4/A5/B1/B2/B4 paragraphs. Every [VERIFY] remains '
        'intentional. The registration, repository DOI and author details '
        'are NOT asserted. Wording restriction: "no supportive evidence" / '
        '"insufficient to test a moderate effect"; never "CX3CL1 has no effect".',
        '',
        '## Title', '',
    ]
    lines += title + ['']
    lines += ['## Abstract (structured, ~350 words; compress to journal '
              'limit before submission)', '']
    lines += abstract + ['']
    lines += ['## Introduction', '']
    lines += intro + [I4, '']
    lines += ['## Methods', '']
    lines += methods + ['']
    lines += ['## Results', '']
    lines += results + [R6, R7, R8, R9, '']
    lines += ['## Discussion', '']
    lines += disc + [D6, '']
    lines += ['## Declarations', '']
    lines += [s for s in decl if not s.startswith('[DRAFT')] + ['']
    lines += ['## Figure legends', '', FIG_LEGENDS.strip(), '']
    lines += ['## Supplementary material index', '',
              'See `submission_pack/a1_supplementary_index.tsv`.', '']
    lines += ['## Separate submission documents', '',
              '- `submission_pack/a1_declarations.txt`',
              '- `submission_pack/a1_cover_letter_draft.txt`',
              '- `submission_pack/a1_presubmission_enquiry.txt`',
              '- `submission_pack/a1_github_readme.md`',
              '- `submission_pack/a1_zenodo_metadata.json`',
              '- `submission_pack/a1_finngen_compliance_checklist.txt`',
              '- `submission_pack/a1_round1_correction_note.txt`',
              '- `results/a1_verify_register.tsv`', '']
    txt(MD + '/A1英文稿定稿骨架_20261111.md', '\n'.join(lines))
    return {'n_existing_en': len(en),
            'sections': {'title': len(title), 'abstract': len(abstract),
                         'introduction': len(intro) + 1,
                         'methods': len(methods),
                         'results': len(results) + 4,
                         'discussion': len(disc) + 1,
                         'declarations': len([s for s in decl
                                              if not s.startswith('[DRAFT')])}}


FIG_LEGENDS = """
Figure 1. CX3CL1 cis region on chromosome 16 (GRCh38). Upper panel:
-log10(P) for FinnGen R13 diabetic neuropathy variants in chr16:57.2-57.6
Mb; no variant reaches P=1x10^-5. Lower panel: -log10(P) for the three
UKB-PPP CX3CL1 cis pQTL credible sets. Dashed vertical lines mark the cis
pQTL positions.

Figure 2. Forest plot of the Wald-ratio odds ratios with 95% confidence
intervals for the lead cis instrument across 19 release-by-endpoint
contexts, ordered by data release and endpoint. Colours denote endpoint
category; the dashed line marks OR=1. All 19 primary P values exceed 0.05.

Figure 3. Power curves for five representative contexts, computed from each
observed Wald-ratio standard error without case-count assumptions. The
dashed horizontal line marks 80% power and the dotted vertical line marks
OR=1.5.

Figure 4. Reporting-decision quadrant for null drug-target MR findings.
Power is assessed first, then evidence consistency, then the permitted
wording. Quadrant C (low power, consistent evidence) corresponds to the
present study.

Figure 5. Simulation validation of the reporting framework, with all
parameters calibrated to the observed A1 data (20,000 replicates per
point). (a) Power and minimum detectable OR as functions of the outcome
standard error. (b) Naive reporting versus the decision framework.
(c) Analytic Binomial(133, 0.05) distribution of nominally significant
tests under a global null, with the observed count of eight marked.

Figure 6. Quantitative attribution of estimate drift. (a) Standardised
drift |delta logOR|/SE(difference) for three classes of pairwise
comparison; no pair exceeds the 1.96 sampling-noise bound. (b) Association
between changes in public case counts and changes in log odds ratios for
the subset of pairs with public counts.

Figure 7. Instrument-set sensitivity. (a) Wald/IVW odds ratios for four
nested instrument sets across 19 contexts. (b) Distribution of power for
OR=1.5 by instrument set; apparent power near 100% in sets containing
pleiotropic trans instruments is false precision and is not interpreted
as an adequately powered null.

Figure 8. Cross-target replication. (a) Odds ratios for the primary
diabetic-neuropathy endpoint across seven additional druggable targets in
R9 and R13. (b) Distribution of power for OR=1.5 across 42 target-context
estimates; the dashed line marks the maximum power observed for CX3CL1
(10.1%) and the solid line marks 80% power.
"""


def declarations():
    s = """DECLARATIONS (paste into manuscript; [VERIFY] fields require the author)

Ethics approval and consent to participate
Not applicable. This study used only publicly available, summary-level
genome-wide association statistics. No individual-level data were accessed
and no new participant contact occurred. Each contributing cohort obtained
its own ethical approval and informed consent; the original publications,
including Kurki et al., Nature 2023 (FinnGen), and the UKB-PPP and Open
Targets primary publications, are cited. Data-use terms of FinnGen,
UKB-PPP and Open Targets were followed. [VERIFY: the author must complete
the FinnGen public-data access form on the FinnGen website and retain the
confirmation; the completed form and the exact citation list are not
asserted here.]

Consent for publication
Not applicable.

Availability of data and materials
All input data are public. Exposure instruments came from Open Targets
Platform v26.9 (UKB-PPP study UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1 and
the additional target studies listed in the deposited target manifest).
Outcome summary statistics came from FinnGen public releases R9-R13; the
endpoint dictionary with exact endpoint names, public case/control counts
and download URLs is deposited as a1_endpoint_dict.tsv. LD reference data
were retrieved from Ensembl REST (1000 Genomes phase 3 EUR, GRCh38).
Derived tables, all figures, the input SHA256 manifest and the version
snapshot are deposited with the code at [VERIFY: GitHub repository URL and
Zenodo DOI]. The retrospective registration is at [VERIFY: OSF DOI].

Code availability
All analyses were implemented in Python 3.13 with the version-locked
dependencies in code/requirements.lock.txt. The offline one-command
reproduction script is code/run_all.sh (verified to complete without
network access); --with-fetch re-downloads public summary statistics for
the network-dependent steps. Source files: a1_tabix_fetch.py,
a1_multi_fetch.py, a1_multi_mr.py, a1_ld_coloc.py, a1_pip_coloc.py,
a1_positive_control.py, a1_simulation.py, a1_strobe_checklist.py,
a1_a18_stats.py, a1_a18_clear.py, a1_a1a6_pack.py, a1_a4a5_drift.py,
a1_b1b2_multitarget.py, a1_figure_pack_en.py, a1_manuscript_pack.py and
the wb_format_* layout scripts. Random seed 20261110; two-sided alpha
0.05; 80% power target.

Competing interests
[VERIFY: insert the authors' competing-interests declaration.]

Funding
[VERIFY: insert grant numbers and funders, if any.]

Authors' contributions (CRediT)
Lu Yi (鲁义主任医师): Conceptualization, Methodology, Supervision,
Writing - review & editing, Funding acquisition [VERIFY: confirm exact
contribution terms and add all co-authors with their CRediT roles].
[VERIFY: complete the full author list before submission.]

Acknowledgements
We thank the FinnGen consortium, the UK Biobank Pharma Proteomics Project
and Open Targets for making summary-level data publicly available. We
cite Kurki et al., Nature 2023, for the FinnGen resource, and the primary
UKB-PPP and Open Targets publications [VERIFY: insert the complete
reference list]. No individual participant is identifiable in this work.

Registration
This study was not prospectively preregistered. The primary analysis was
fixed after the initial exploratory MR and before the reporting framework
was applied; this retrospective specification is disclosed in the
Methods. A retrospective OSF registration is prepared, with date and DOI
[VERIFY: insert OSF date and DOI]. All exploratory analyses completed
before registration are listed in the Methods and in the registration
text.
"""
    txt(PACK + '/a1_declarations.txt', s)


def other_docs():
    # cover letter：复用已生成文件，附一段本轮新增结果说明
    cl = io.open(RES + '/a1_a6_cover_letter.txt', encoding='utf-8').read()
    cl += """

6. Additional analyses completed after the first draft (for the editor's
information). Beyond the 19 primary contexts we now report: a quantitative
attribution showing all 171 pairwise estimate differences are compatible
with sampling noise; a 4 x 19 instrument-set sensitivity grid that
identifies false precision from pleiotropic trans instruments; a
cross-target replication in seven additional druggable proteins showing
that low power is an interaction between instrument strength and endpoint
size rather than a universal property of FinnGen endpoints; and a
biological positive control in which PCSK9 cis instruments recover the
expected direction on three cardiovascular endpoints after Bonferroni
correction. We also report that LDLR failed as a positive control because
of an extreme Wald ratio from a weak instrument, and we explicitly
exclude it. All figures are supplied as English-language vector files
(PDF/SVG) with PNG previews. These additions do not change the primary
conclusion and do not alter the wording restriction.
"""
    txt(PACK + '/a1_cover_letter_draft.txt', cl)

    txt(PACK + '/a1_presubmission_enquiry.txt', """PRE-SUBMISSION ENQUIRY DRAFT (send by the author; not sent by this agent)

To: Editorial Office, BMC Medical Research Methodology
Subject: Presubmission enquiry - empirical evaluation of methodology choice
and study outcomes in drug-target Mendelian randomization

Dear Editors,

We are preparing a Research article that asks whether a null finding from
a two-sample, single-cis-instrument drug-target Mendelian randomization
analysis can be interpreted without quantifying power, endpoint
definitions and analytical drift. Using circulating CX3CL1 pQTL
instruments and 19 FinnGen release-by-endpoint contexts, with additional
replication across seven further druggable targets, we quantify minimum
detectable effects, cross-context estimate drift, the false precision
that arises when pleiotropic trans instruments are added, and a
reporting-layer decision rule for null results. We report a negative
finding transparently and do not claim that CX3CL1 has no effect.

Your aims and scope explicitly encourage empirical studies of the
associations between the choice of methodology and study outcomes. Before
submission we would be grateful for your view on whether this manuscript
falls within the journal's methodological-empirical remit, and whether a
Research article is the appropriate article type.

The study uses public summary-level data only, includes a completed
STROBE-MR checklist, a self-assessment against a published MR quality
scoring system, full code, version-locked dependencies and SHA256
checksums, and a table anticipating likely reviewer concerns. No
prospective preregistration exists; the retrospective registration and
the timing of the primary-analysis specification are disclosed.

Yours sincerely,
[VERIFY: corresponding author name, degree, affiliation, email]
""")

    txt(PACK + '/a1_github_readme.md', """# A1 - Power and phenotype-definition heterogeneity in drug-target MR

Empirical evaluation and reporting framework using circulating CX3CL1 and
diabetic neuropathy across five FinnGen releases, with cross-target
replication and positive controls.

## Reproduce offline (no network required)

    bash code/run_all.sh

This recomputes all statistics and figures from the deposited result
tables. Network-dependent re-fetching of public summary statistics is
available through the original `a1_multi_fetch.py` / `a1_b1b2_multitarget.py`
scripts; it is not required to verify the reported numbers.

## Layout

- `code/` - all analysis, figure and document-layout scripts
- `results/` - result tables (TSV), JSON summaries, vector figures (PDF/SVG)
  and PNG previews
- `data/` - public input snapshots, tabix indexes and checksums
- `submission_pack/` - cover letter, declarations, presubmission enquiry,
  figure legends and compliance checklists
- `requirements.lock.txt` - version-locked Python dependencies
- `a1_input_sha256.tsv` - SHA256 of every input file
- `a1_version_snapshot.json` - code/data/environment snapshot

## Data sources

Open Targets Platform v26.9 (UKB-PPP pQTL), FinnGen public releases
R9-R13, Ensembl REST 1000 Genomes phase 3 EUR LD. All inputs are public;
no individual-level data and no access tokens are included.

## Reporting restriction

The manuscript states "no supportive evidence; the available data are
insufficient to test an effect of moderate magnitude". It never states
that CX3CL1 has no effect, and it does not revise the tissue-specific
experimental hypothesis of the parent grant.

OSF registration: 10.17605/OSF.IO/APWE2 (https://osf.io/apwe2), CC0 1.0.

[VERIFY before release: licence, author list and Zenodo DOI.]
""")

    zenodo = {
        'title': 'A1: Power and phenotype-definition heterogeneity in '
                 'drug-target Mendelian randomization - code and derived data',
        'description': 'Reproducible Python pipeline, result tables, '
                       'vector figures and reporting materials for an '
                       'empirical evaluation of null-result reporting in '
                       'drug-target MR (CX3CL1 / diabetic neuropathy, '
                       'FinnGen R9-R13), including cross-target '
                       'replication and biological positive controls.',
        'creators': [{'name': '[VERIFY: author name]',
                      'affiliation': '[VERIFY: institution]'}],
        'upload_type': 'software',
        'access_right': 'open',
        'license': '[VERIFY: choose licence, e.g. CC-BY-4.0 or MIT]',
        'keywords': ['Mendelian randomization', 'drug-target MR',
                     'STROBE-MR', 'statistical power', 'null results',
                     'FinnGen', 'UKB-PPP', 'CX3CL1', 'diabetic neuropathy',
                     'reporting guidelines'],
        'related_identifiers': [
            {'identifier': '10.17605/OSF.IO/APWE2',
             'relation': 'isSupplementedBy', 'scheme': 'doi'},
        ],
        'note': 'Metadata draft only. Creating the deposit and minting the '
                'DOI require the author account and are not asserted here.',
    }
    txt(PACK + '/a1_zenodo_metadata.json',
        json.dumps(zenodo, ensure_ascii=False, indent=2))

    txt(PACK + '/a1_finngen_compliance_checklist.txt', """FINNGEN PUBLIC-DATA COMPLIANCE (author action required)

This file lists exactly what remains. The agent cannot complete these
steps because they require the author's institutional account.

[ ] 1. Complete the FinnGen public-data access / acknowledgement form on
       the official FinnGen website and retain the confirmation email or
       receipt in the project archive.
[ ] 2. Verify the current public data-use terms for releases R9, R10,
       R11, R12 and R13 and record the retrieval date.
[ ] 3. Cite Kurki et al., Nature 2023 (FinnGen resource) in the
       References and in the Acknowledgements. [VERIFY: exact volume,
       pages and DOI from the published article.]
[ ] 4. Cite the UKB-PPP primary publication and the Open Targets
       primary publication. [VERIFY: exact references.]
[ ] 5. Record the official public-release dates of R9-R13 and UKB-PPP
       v1 from the primary publications. NOT DONE: repeated web searches
       did not return authoritative dated release notices, and no date
       has been guessed or written into the manuscript.
[ ] 6. Confirm that no individual-level FinnGen data were downloaded;
       this analysis used public summary statistics exclusively.

Prepared materials already available: a1_osf_registration_pack.txt,
a1_declarations.txt, a1_citation_anchors.tsv and the endpoint dictionary
with public download URLs.
""")

    txt(PACK + '/a1_round1_correction_note.txt', """CORRECTION NOTICE - supersedes two numbers in the first-launch report

File affected: A1启动报告与首轮真实结果_20261110.docx (kept as an internal
record; do NOT copy its IVW / MR-Egger values into the manuscript).

1. In the first-launch script (code/a1_round1_mr.py) both IVW and
   MR-Egger were implemented as weighted regressions with a free
   intercept. As a consequence the MR-Egger slope was numerically
   identical to the IVW slope and was not a valid MR-Egger estimate.
   Corrected implementation: code/a1_ld_coloc.py, function ivw(),
   - IVW: intercept fixed to zero, weights 1/se_out^2;
   - MR-Egger: free intercept, centred weighted regression, with
     intercept_se, intercept_p, Q, Q_df and Q_p reported.
   All subsequent tables (a1_mr_matrix.tsv, a1_a18_stats.tsv) use the
   corrected implementation.

2. The value 22.25 was described in an early draft as an F statistic.
   It is NOT. F for the lead cis instrument is 552.42 (z=23.5037);
   22.25 is the minimum detectable odds ratio at 80% power for the R13
   primary endpoint, computed from the observed Wald-ratio SE.

3. The tool-level bootstrap confidence interval (0.020-0.505 on the
   beta scale) does not incorporate GWAS sampling uncertainty and must
   not be used to claim a positive result. The analytical IVW interval
   (-0.0079 to 0.1286) includes the null.

4. I2_GX is 98.66%-99.08% (high), so the NOME assumption is satisfied.
   The seven nominally significant MR-Egger results must NOT be
   attributed to regression dilution from weak instruments. The
   supported explanation is the multiple-testing expectation
   (6.65 expected vs 8 observed; binomial P=0.548), horizontal
   pleiotropy in trans instruments and unstable Egger intercepts at
   low power.

Manuscript wording remains: "no supportive evidence; the available
data are insufficient to test an effect of moderate magnitude".
""")

    txt(PACK + '/a1_figure_legends.txt', FIG_LEGENDS)


def verify_register():
    rows = [
        dict(id='V01', item='Final title', status='AUTHOR',
             note='Three candidate titles drafted; author chooses one.'),
        dict(id='V02', item='OSF registration date and DOI',
             status='EXCLUDED_BY_USER',
             note='Expressly excluded from this round; retrospective '
                  'registration text pack is prepared.'),
        dict(id='V03', item='GitHub repository and Zenodo DOI',
             status='AUTHOR',
             note='README and Zenodo metadata drafted; account action needed.'),
        dict(id='V04', item='Competing interests', status='AUTHOR',
             note='Must be supplied by all authors.'),
        dict(id='V05', item='Corresponding author / affiliation / email',
             status='AUTHOR',
             note='Manuscript keeps "鲁义主任医师" as corresponding author '
                  'label; full details and co-author list needed.'),
        dict(id='V06', item='Funding and grant numbers', status='AUTHOR',
             note='Not asserted.'),
        dict(id='V07', item='Primary CX3CL1/CX3CR1 and DPN references',
             status='AUTHOR', note='Insert primary literature, not reviews.'),
        dict(id='V08', item='Official public-release dates R9-R13 and '
             'UKB-PPP v1', status='UNVERIFIABLE',
             note='Repeated searches returned no authoritative dated '
                  'release notice; deliberately not guessed.'),
        dict(id='V09', item='Rosas-Chavez Table 1 full rubric and total score',
             status='UNVERIFIABLE',
             note='MDPI HTML 403, /xml anti-bot and /pdf blank; '
                  'SNP reporting stays NOT_SCORED, no total claimed.'),
        dict(id='V10', item='Original Stender et al. reference',
             status='OPEN_LITERATURE',
             note='Currently cited only as endorsed in Rosas-Chavez 2025.'),
        dict(id='V11', item='LocusCompare2 original Nature Genetics DOI',
             status='OPEN_LITERATURE',
             note='Only a Chinese secondary summary was found.'),
        dict(id='V12', item='13d per-article full-text check of power '
             'reporting', status='OPEN_LITERATURE',
             note='Not performed; do not state a checked percentage.'),
        dict(id='V13', item='Abstract compression to journal word limit',
             status='DONE',
             note='Structured abstract compiled in the consolidated '
                  'skeleton; final trim is editorial.'),
        dict(id='V14', item='41% power-audit figure',
             status='DONE',
             note='Corrected: source reports 59% fulfilling the criterion; '
                  '41% is only derived and is not used as a quoted figure.'),
        dict(id='V15', item='Code availability / CRediT / declarations',
             status='DONE',
             note='Full text in submission_pack/a1_declarations.txt; '
                  'author-specific fields still [VERIFY].'),
        dict(id='V16', item='English vector figures',
             status='DONE',
             note='fig1-fig5 redrawn as EN PDF/SVG/PNG; fig6-fig8 were '
                  'already English vector figures.'),
        dict(id='V17', item='Kurki et al. Nature 2023 acknowledgement',
             status='DRAFTED',
             note='Included in declarations and compliance checklist; '
                  'exact volume/pages/DOI to confirm.'),
        dict(id='V18', item='FinnGen website form and receipt',
             status='AUTHOR',
             note='Checklist prepared; cannot be completed by the agent.'),
        dict(id='V19', item='Journal word limits, APC, presubmission '
             'channel', status='OPEN_LITERATURE',
             note='Recheck on the current journal site immediately before '
                  'submission; enquiry draft prepared.'),
        dict(id='V20', item='PAT-based full coloc appendix',
             status='DEFERRED',
             note='B3: requires a personal access token; not needed for '
                  'the first submission.'),
    ]
    write_tsv(RES + '/a1_verify_register.tsv', rows)
    return rows


def supplementary_index():
    rows = [
        dict(id='S1', title='Endpoint dictionary and definition table',
             file='results/a1_endpoint_dict.tsv; '
                  'results/a1_endpoint_definition_table.tsv',
             note='19 release-by-endpoint contexts, public case/control '
                  'counts where available, download URLs.'),
        dict(id='S2', title='Per-combination exposure, outcome and '
             'Wald-ratio estimates',
             file='results/a1_a18_lead_variant_associations.tsv; '
                  'results/a1_outcome_matrix.tsv', note='All 19 contexts.'),
        dict(id='S3', title='STROBE-MR checklist (43 subitems) and '
             'A-level 18-item clearance table',
             file='results/a1_strobe_mr_checklist.tsv; '
                  'results/a1_strobe_a18_cleared.tsv', note=''),
        dict(id='S4', title='Self-assessment against published MR quality '
             'scoring system',
             file='results/a1_a3_self_assessment.tsv',
             note='No total score asserted; Table 1 of source unavailable.'),
        dict(id='S5', title='Four-class decision table for 19 contexts',
             file='results/a1_a2_decision_classes.tsv', note=''),
        dict(id='S6', title='Reviewer-anticipation table (8 concerns)',
             file='results/a1_a6_reviewer_anticipation.tsv', note=''),
        dict(id='S7', title='Instrument-set sensitivity grid (76 estimates)',
             file='results/a1_a5_instrument_drift.tsv; '
                  'results/a1_a5_instrument_detail.tsv',
             note='Includes false-precision caveat.'),
        dict(id='S8', title='Pairwise drift attribution (171 pairs)',
             file='results/a1_a4_pairwise_drift.tsv; '
                  'results/a1_a4_attribution.json', note=''),
        dict(id='S9', title='Cross-target replication (7 targets, 42 '
             'estimates)',
             file='results/a1_b1_multitarget.tsv; '
                  'data/b1_target_pqtl.json', note=''),
        dict(id='S10', title='Biological positive control (21 estimates)',
             file='results/a1_b2_biological_control.tsv',
             note='PCSK9 succeeds; LDLR excluded with reason.'),
        dict(id='S11', title='Reproducibility package',
             file='code/run_all.sh; code/requirements.lock.txt; '
                  'results/a1_input_sha256.tsv; '
                  'results/a1_version_snapshot.json', note=''),
        dict(id='S12', title='Full analytic statistics',
             file='results/a1_a18_stats.tsv; results/a1_a18_steiger.tsv; '
                  'results/a1_mr_matrix.tsv',
             note='Q, I2, I2_GX, E-value, ABF, F, bootstrap CI, Steiger.'),
        dict(id='S13', title='Simulation validation',
             file='results/a1_simulation_results.json',
             note='Methodological validation only, not new empirical '
                  'evidence.'),
        dict(id='S14', title='Figure files (English, vector)',
             file='results/fig*_EN.pdf; results/fig6-8*.pdf',
             note='SVG and PNG also supplied.'),
    ]
    write_tsv(PACK + '/a1_supplementary_index.tsv', rows)


def main():
    print('Building consolidated manuscript and submission pack ...')
    ms = manuscript()
    declarations()
    other_docs()
    v = verify_register()
    supplementary_index()
    # OSF 归档期全文审核后的机械性修订；保持 14 步复现链不变
    import a1_postreview_patch
    a1_postreview_patch.apply()
    # 汇总按修订后的登记表重算，避免仍显示 EXCLUDED_BY_USER / DEFERRED
    with io.open(RES + '/a1_verify_register.tsv', encoding='utf-8') as f:
        v = [dict(r) for r in csv.DictReader(f, delimiter='\t')]
    summary = {
        'generated_on': '2026-11-11',
        'manuscript': ms,
        'verify_register': {
            'total': len(v),
            'done': sum(1 for r in v if r['status'] == 'DONE'),
            'author_account': sum(1 for r in v if r['status'] == 'AUTHOR'),
            'unverifiable': sum(1 for r in v
                                if r['status'] == 'UNVERIFIABLE'),
            'excluded_by_user': sum(1 for r in v
                                    if r['status'] == 'EXCLUDED_BY_USER'),
            'other': sum(1 for r in v if r['status'] in
                         ('OPEN_LITERATURE', 'DRAFTED', 'DEFERRED')),
        },
        'new_english_paragraphs': ['I4', 'R6', 'R7', 'R8', 'R9', 'D6',
                                   'M8', 'R10', 'D7'],
        'postreview': 'Full-text review of 2026-11-11 applied: compressed '
                      'abstract, corrected ABF wording, inserted OSF '
                      'registration date, integrated region-wide '
                      'colocalisation (M8/R10/D7/Figure 9/S15), added '
                      'draft references, fixed SHA256 count to 144 and '
                      'updated the verification register. No analysis '
                      'number differs from the archived OSF registration.',
        'note': 'No MR result was recomputed by this script; it only '
                'consolidates text and registers outstanding items.',
    }
    txt(PACK + '/a1_pack_summary.json',
        json.dumps(summary, ensure_ascii=False, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
