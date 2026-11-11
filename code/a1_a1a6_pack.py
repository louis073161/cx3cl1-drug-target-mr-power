# -*- coding: utf-8 -*-
"""A1-A2-A3-A6 补正包：按已落盘真实结果生成四分类判定、评分体系自评、
审稿意见预判表与 cover letter。全部为报告层材料，不重跑 MR、不产生新实证结论。

数据来源：a1_mr_matrix.tsv / a1_power_curve.tsv / a1_simulation_results.json /
a1_drift_summary.json / a1_strobe_a18_cleared.tsv 及既有 Word/Markdown 交付。
"""
import csv
import io
import json
import math
import os
from collections import defaultdict

RES = 'F:/博士/国自然科学基金/2025神自然/A1_启动/results'
TODAY = '2026-11-10'
Z_ALPHA = 1.959963984540054   # two-sided alpha=0.05
Z_POWER = 0.8416212335729143  # 80% power


def read_tsv(name):
    with io.open(os.path.join(RES, name), encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(name, rows):
    path = os.path.join(RES, name)
    if not rows:
        raise ValueError('empty rows: ' + name)
    cols = list(rows[0].keys())
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t',
                           lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def build_decision_classes():
    """A2：19 个 release×endpoint 组合逐行四分类，规则全部预设并写入 rationale。"""
    wald = [r for r in read_tsv('a1_mr_matrix.tsv')
            if r['method'] == 'Wald ratio (lead cis)']
    power = {(r['release'], r['endpoint']): float(r['power'])
             for r in read_tsv('a1_power_curve.tsv') if r['scenario'] == 'OR_1.5'}
    sim = json.load(io.open(os.path.join(RES, 'a1_simulation_results.json'),
                            encoding='utf-8'))
    bonf = sim['calibration']['bonferroni_threshold']

    by_endpoint = defaultdict(list)
    for r in wald:
        by_endpoint[r['endpoint']].append(float(r['or']))
    def direction_inconsistent(ep):
        vals = by_endpoint[ep]
        return len(vals) > 1 and any(v < 1 for v in vals) and any(v > 1 for v in vals)

    labels = {
        'I': 'I - supported association',
        'II': 'II - informative null (adequately powered)',
        'III': 'III - underpowered / not evaluable',
        'IV': 'IV - context-inconsistent',
    }
    rows = []
    for r in wald:
        key = (r['release'], r['endpoint'])
        se, p, orv = float(r['se']), float(r['p']), float(r['or'])
        ci_lo_or, ci_hi_or = math.exp(float(r['ci_lo'])), math.exp(float(r['ci_hi']))
        pw = power[key]
        mde = math.exp((Z_ALPHA + Z_POWER) * se)
        inc = direction_inconsistent(r['endpoint'])
        if p < bonf:
            cls = 'I'
            why = ('P=%.4g below Bonferroni threshold %.4g;'
                   ' would require replication before any causal claim' % (p, bonf))
        elif pw >= 0.8 and ci_hi_or < 1.5:
            cls = 'II'
            why = ('power at OR=1.5 %.1f%% >=80%% and CI excludes OR 1.5;'
                   ' an informative null could be reported' % (100 * pw))
        elif inc:
            cls = 'IV'
            why = ('point estimate changes sign across releases of the same endpoint;'
                   ' power at OR=1.5 only %.2f%% and MDE OR %.2f, so this is'
                   ' context-inconsistent AND underpowered, not a null'
                   % (100 * pw, mde))
        else:
            cls = 'III'
            why = ('no Bonferroni-significant result (P=%.3g), power at OR=1.5'
                   ' %.2f%% <80%%, MDE OR %.2f > any plausible moderate effect;'
                   ' the data cannot evaluate the hypothesis' % (p, 100 * pw, mde))
        rows.append({
            'release': r['release'], 'endpoint': r['endpoint'],
            'endpoint_type': r['endpoint_type'],
            'wald_beta': r['beta'], 'wald_se': r['se'], 'wald_p': r['p'],
            'or': '%.4f' % orv,
            'ci_lo_or': '%.4f' % ci_lo_or, 'ci_hi_or': '%.4f' % ci_hi_or,
            'power_or_1.5': '%.6f' % pw,
            'mde_or_80pct': '%.4f' % mde,
            'bonferroni_threshold': '%.8g' % bonf,
            'direction_inconsistent_within_endpoint': 'yes' if inc else 'no',
            'decision_class': cls, 'decision_class_label': labels[cls],
            'rationale': why,
            'evidence_files': 'a1_mr_matrix.tsv; a1_power_curve.tsv; a1_simulation_results.json',
            'generated': TODAY,
        })
    n = write_tsv('a1_a2_decision_classes.tsv', rows)
    counts = defaultdict(int)
    for r in rows:
        counts[r['decision_class']] += 1
    return rows, dict(counts), bonf


SELF_ASSESS = [
    # criterion, category, points_in_text, our_evidence, points, basis, note
    ('Study rationale / prior evidence', 'study design', '+2 / +1 / -1',
     'Intro cites STROBE-MR, Corbin 2025, Rosas-Chavez 2025 and prior CX3CL1/DPN MR;'
     ' the rationale is a reporting gap, not large observational evidence',
     '+1', 'explicit_in_text',
     'Conservative self-rating: mixed/small-study level; not +2 (no large observational study or trial)'),
    ('Bidirectional MR', 'study design', '+1',
     'Reverse MR not performed; Steiger orientation used with explicit biological reason (13c)',
     '0', 'explicit_in_text', 'Scored as not obtained; no deduction stated in source text'),
    ('Use of most recent GWAS data', 'study design', '0 or +1',
     'FinnGen R13 (latest release available at fetch) and UKB-PPP v1 used',
     '+1', 'explicit_in_text',
     'Official release dates still [VERIFY]; rating refers to version available at data fetch'),
    ('Matching ancestries', 'study design', '-1 to +1',
     'UKB-PPP European ancestry and FinnGen Finnish (European) ancestry',
     '+1', 'inferred_from_range',
     'Both European but not the same population; -1/+1 reconstructed from stated range'),
    ('Absence of participant overlap', 'study design', '-1 or +1',
     'UK (UKB-PPP) vs Finland (FinnGen); zero participant overlap reported',
     '+1', 'explicit_in_text', 'STROBE-MR item 10'),
    ('Replication in an independent dataset', 'study design', '+3 / -1',
     'Nineteen release-by-endpoint reruns use the SAME FinnGen resource;'
     ' no independent external outcome dataset was analysed',
     '-1', 'explicit_in_text',
     'HONEST DEDUCTION: cross-context comparison is NOT independent replication.'
     ' This is a stated limitation and motivates the multi-target replication queue (B1)'),
    ('Instrument strength (relevance)', 'statistical methods', '+1',
     'Lead cis F=552.42, exposure P=3.74e-122; threshold P<5e-8 or F>10 met',
     '+1', 'explicit_in_text', 'No weak-instrument claim permitted'),
    ('LD pruning', 'statistical methods', '+1',
     'Ensembl REST / 1000 Genomes EUR clumping at r2<0.001',
     '+1', 'explicit_in_text', 'Stricter than the source rubric threshold R2<0.1'),
    ('Confounders / independence assumption', 'statistical methods', '+1',
     'Trans instruments with extensive pleiotropic colocalisation excluded from primary analysis;'
     ' negative-control endpoints reported',
     '+1', 'explicit_in_text',
     'Based on precomputed colocalisation; no Phenoscanner-style full confounder lookup performed'),
    ('Mediator analysis beyond standard MR packages', 'statistical methods', '+1',
     'Not performed; this is a reporting study, not a mediation study',
     '0', 'explicit_in_text', 'Standard MR-Egger/PRESSO-type checks are not scored by the rubric'),
    ('SNP reporting', 'statistical methods', 'NOT STATED IN TEXT',
     'Lead variant beta/SE/P for exposure, outcome and Wald ratio reported for all 19 combos'
     ' (a1_a18_lead_variant_associations.tsv, Supplementary Table S2)',
     'NOT_SCORED', 'not_in_text',
     'Exact points/deductions unavailable because Table 1 body could not be retrieved'),
    ('Multiple-testing correction', 'statistical methods', 'up to +3',
     'Bonferroni threshold 3.7594e-4 applied across 133 tests; 8 nominal results'
     ' vs 6.65 expected (binomial P=0.548); none survived correction',
     '+3', 'explicit_in_text',
     'Upper bound of stated range; exact sub-tier needs Table 1'),
    ('Study power', 'statistical methods', '+2',
     'Post-hoc power and MDE derived from observed SEs for every combination;'
     ' absence of a priori calculation explicitly disclosed (STROBE-MR 4b)',
     '+2', 'explicit_in_text',
     'Rubric accepts power considered/calculated in Methods or Results'),
    ('Correct interpretation of results', 'interpretation', '+2 / -2',
     'No causal claim from uncorrected P values; pleiotropy, method concordance and'
     ' CI width discussed; wording restricted to "no supportive evidence"',
     '+2', 'explicit_in_text', 'STROBE-MR items 15/16a/16b'),
    ('STROBE-MR compliance table', 'STROBE-MR', '+1 (bonus)',
     '43-subitem STROBE-MR checklist deposited (a1_strobe_a18_cleared.tsv)',
     '+1', 'explicit_in_text', 'Bonus criterion added late by the source authors'),
]


def build_self_assessment():
    rows = []
    for crit, cat, pts, ev, our, basis, note in SELF_ASSESS:
        rows.append({
            'criterion': crit, 'category': cat,
            'points_stated_in_source_text': pts,
            'our_evidence': ev,
            'self_assessed_points': our,
            'point_basis': basis,
            'note': note,
            'total_score': 'NOT SUMMED',
            'total_score_reason':
                'Table 1 body could not be retrieved (publisher anti-bot on /xml);'
                ' SNP reporting points and correction sub-tiers are unstated in the'
                ' accessible text. Per-item ratings only; no total is asserted.',
            'source': 'Rosas-Chavez & Merriman 2025, doi:10.3390/gucdd3020008',
            'evidence_files': 'a1_strobe_mr_checklist.tsv; a1_a18_stats.tsv;'
                              ' a1_mr_matrix.tsv; a1_positive_control.tsv',
            'generated': TODAY,
        })
    return rows, write_tsv('a1_a3_self_assessment.tsv', rows)


REVIEWER = [
    ('R1', 'A single target cannot justify a general reporting framework.',
     'We explicitly frame this as a worked case study with 19 cross-context comparisons,'
     ' 20,000-replicate-per-point simulation and 56 calibration combinations, and we'
     ' list single-target generalisability as the main limitation; a multi-target'
     ' replication queue is stated as the next step, not claimed as completed.',
     'Discussion (Limitations); Methods (Simulation)',
     'a1_simulation_results.json; a1_mr_matrix.tsv; a1_positive_control.tsv'),
    ('R2', 'The null is just an after-the-fact "lack of power" excuse.',
     'No a priori power calculation was performed and this is disclosed under'
     ' STROBE-MR 4b. All power statements are quantified from observed SEs:'
     ' power at OR=1.5 5.45-10.10%, MDE OR 5.61-318.57, R13 primary CI OR 0.159-12.21.'
     ' We do not write "lack of power" without these numbers.',
     'Methods 4b; Results (Power); Supplementary decision-class table',
     'a1_power_curve.tsv; a1_a2_decision_classes.tsv'),
    ('R3', 'Eight nominally significant results look like selective reporting.',
     'Eight of 133 tests were nominally significant, versus 6.65 expected under a'
     ' global null (binomial P=0.548); none survived Bonferroni correction'
     ' (threshold 3.7594e-4). All 133 results are deposited, not selected.',
     'Results (Multiplicity)', 'a1_mr_matrix.tsv; a1_simulation_results.json'),
    ('R4', 'Weak instruments could explain the null.',
     'Lead cis F=552.42, I2_GX 98.66-99.08% (NOME satisfied), I2 0-70.58%;'
     ' the earlier "weak instrument" wording was corrected and is not used.'
     ' An earlier report mislabelled the MDE OR of 22.25 as an F statistic.',
     'Results 12b; Discussion (Limitations)', 'a1_a18_stats.tsv; a1_strobe_a18_cleared.tsv'),
    ('R5', 'Why are seven MR-Egger results nominally significant?',
     'Not regression dilution: I2_GX is high. The pattern is consistent with the'
     ' expected false-positive count under a global null, horizontal pleiotropy of'
     ' trans instruments (H4 up to 1.0) and unstable Egger intercepts under low power.',
     'Results 12b', 'a1_a18_stats.tsv; a1_robustness.json'),
    ('R6', 'The bootstrap CI excludes zero - why not report a positive effect?',
     'The instrument-level bootstrap interval (0.020-0.505) resamples instruments only'
     ' and omits GWAS sampling error. The analytic IVW 95% CI (-0.0079 to 0.1286)'
     ' includes the null and is the inferential interval; no positive claim is made.',
     'Results 12b', 'a1_a18_stats.tsv'),
    ('R7', 'How does this differ from Hu 2025, Zou 2024 and LocusCompare2?',
     'Hu 2025 addresses phenotype source (self-report vs clinical) across AD datasets;'
     ' we address within-database release version and endpoint definition for a'
     ' drug-target pQTL, quantified by power/MDE. Zou 2024 is a Bayesian posterior'
     ' decision rule (ROPE, ternary); ours is a reporting-layer classification and'
     ' does not propose a new estimator. LocusCompare2 addresses discordance between'
     ' gene-implication methods, not release/endpoint drift.',
     'Introduction (differentiation); Discussion',
     'a1_citation_anchors.tsv'),
    ('R8', 'This is only a re-analysis of public MR results with no independent evidence.',
     'We frame the 19 reruns as cross-context comparisons, add negative-control'
     ' endpoints, a technical positive control (P as small as 1.37e-19 in ABO trans'
     ' instruments, all 14 cis calibration combinations non-significant), a simulation'
     ' validation and external clinical/preclinical literature as context. We state'
     ' explicitly that these are not fully independent evidence sources and that'
     ' colocalisation and reverse MR are not counted as independent approaches.',
     'Methods (triangulation); Discussion',
     'a1_positive_control.tsv; a1_drift_summary.json; a1_endpoint_definition_table.tsv'),
]


def build_reviewer_table():
    rows = []
    for rid, comment, response, loc, ev in REVIEWER:
        rows.append({
            'id': rid, 'anticipated_reviewer_comment': comment,
            'our_response': response, 'manuscript_location': loc,
            'evidence_files': ev,
            'note': 'Prepared for transparent peer review; reviewer reports are'
                    ' published alongside the article by BMC Medical Research Methodology.',
            'generated': TODAY,
        })
    return rows, write_tsv('a1_a6_reviewer_anticipation.tsv', rows)


COVER_LETTER = """[DRAFT - replace every [VERIFY] field before submission]

To the Editors, BMC Medical Research Methodology

Dear Editors,

Please consider our manuscript, "[VERIFY: final title]", for publication as a
Research article.

1. Why this manuscript suits BMC Medical Research Methodology
Your aims and scope explicitly encourage empirical studies of the associations
between the choice of methodology and study outcomes. Our study is exactly such
an empirical evaluation: using one drug-target Mendelian randomization (MR) case
(circulating CX3CL1 and diabetic neuropathy), we quantified how endpoint
definitions, GWAS data releases and estimator choices change the reported
conclusion, and how much statistical power the publicly available data actually
provide. The contribution is a reporting framework for null MR findings, not a
claim about CX3CL1 biology. Because the journal does not make editorial decisions
on perceived interest or impact, we believe a scientifically valid negative
result with a reproducible methodological evaluation is an appropriate fit.

2. Journal-policy matters
We enclose a completed 43-subitem STROBE-MR checklist. Reporting guidelines,
data availability, ethics and competing interests are addressed in the
Declarations. The study uses only publicly available summary-level statistics
and required no additional ethical approval. Registration is retrospective;
the date and DOI are [VERIFY: insert OSF registration date and DOI], and all
analyses completed before registration are disclosed in the Methods. We do not
claim that CX3CL1 has no effect on diabetic neuropathy, and we do not claim to
propose a new statistical estimator. All analytical code, version-locked
dependencies, SHA256 checksums of 75 input files, a version snapshot and a
one-command offline reproduction script are available at
[VERIFY: insert GitHub repository and Zenodo DOI].

3. Competing interests
[VERIFY: insert the authors' competing-interests declaration].

4. Approval of all authors
All authors have read and approved the manuscript for submission. The
corresponding author is [VERIFY: insert name, degree and institutional
affiliation].

5. Originality and duplicate publication
The content has not been published previously and is not under consideration
for publication elsewhere. No part of the manuscript has been split from this
work for separate submission. The analyses form part of a broader
methodological programme; no salami publication is intended.

Additional points for the editor
- The primary analysis is designated explicitly in the Methods, with a
  transparent statement that the designation was fixed after the initial
  exploratory MR and before the reporting framework was applied.
- Triangulation is partial and is described as such: cross-context comparisons,
  negative-control endpoints, a technical positive control and simulation
  validation are presented, while colocalisation and reverse MR are not counted
  as independent approaches.
- A table anticipating eight likely reviewer concerns is supplied as
  Supplementary Material, in view of the journal's transparent peer-review
  policy.

Thank you for considering this submission.

Sincerely,
[VERIFY: corresponding author name and signature]
[VERIFY: institution, address, email]
"""


def main():
    rows, counts, bonf = build_decision_classes()
    sa_rows, n_sa = build_self_assessment()
    rv_rows, n_rv = build_reviewer_table()
    with io.open(os.path.join(RES, 'a1_a6_cover_letter.txt'), 'w',
                 encoding='utf-8') as f:
        f.write(COVER_LETTER)
    summary = {
        'generated': TODAY,
        'scope': 'A1 triangulation text; A2 decision classes; A3 self-assessment;'
                 ' A6 reviewer anticipation and cover letter',
        'no_new_empirical_analysis': True,
        'a2': {
            'n_combinations': len(rows),
            'bonferroni_threshold': bonf,
            'decision_class_counts': counts,
            'all_power_limited': all(float(r['power_or_1.5']) < 0.8 for r in rows),
            'n_direction_inconsistent': sum(
                1 for r in rows if r['direction_inconsistent_within_endpoint'] == 'yes'),
            'power_range': [min(float(r['power_or_1.5']) for r in rows),
                            max(float(r['power_or_1.5']) for r in rows)],
            'mde_range': [min(float(r['mde_or_80pct']) for r in rows),
                          max(float(r['mde_or_80pct']) for r in rows)],
            'mde_formula': 'exp((z_{1-alpha/2}+z_{power}) * SE_Wald)',
            'class_rule': 'I if P<Bonferroni; II if power>=0.8 and upper CI OR<1.5;'
                          ' IV if direction inconsistent within endpoint; else III',
        },
        'a3': {'n_items': n_sa,
               'total_score': 'NOT SUMMED',
               'reason': 'Table 1 body not retrievable (publisher anti-bot);'
                         ' SNP reporting points and correction sub-tiers unstated',
               'n_not_scored': sum(1 for r in sa_rows
                                   if r['self_assessed_points'] == 'NOT_SCORED')},
        'a6': {'n_anticipated_comments': n_rv,
               'cover_letter_items': 5,
               'cover_letter_verifies': COVER_LETTER.count('[VERIFY')},
        'outputs': ['a1_a2_decision_classes.tsv', 'a1_a3_self_assessment.tsv',
                    'a1_a6_reviewer_anticipation.tsv', 'a1_a6_cover_letter.txt',
                    'a1_a1a6_pack_summary.json'],
    }
    with io.open(os.path.join(RES, 'a1_a1a6_pack_summary.json'), 'w',
                 encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
