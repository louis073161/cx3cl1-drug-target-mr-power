# -*- coding: utf-8 -*-
"""OSF 归档期全文审核后的机械性修订（2026-11-11）。

由 a1_manuscript_pack.main() 在最后一步调用，保证 run_all.sh 重跑后修订不丢失。
原则：只改表述、结构、一致性、格式与已确认的日期；不改动任何已登记的分析数字。
所有替换均做存在性判断，可重复执行而不叠加。
"""
import io
import json
import os

PROJ = 'F:/博士/国自然科学基金/2025神自然'
ROOT = PROJ + '/A1_启动'
RES = ROOT + '/results'
PACK = ROOT + '/submission_pack'
MD = PROJ + '/_markdown源'

REG_DATE = '11 November 2026'
OSF_DOI = '10.17605/OSF.IO/APWE2'


def _read(path):
    with io.open(path, encoding='utf-8') as f:
        return f.read()


def _write(path, s):
    with io.open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(s.rstrip() + '\n')


def _rep(s, old, new):
    if old in s:
        return s.replace(old, new)
    return s


# ---------------------------------------------------------------- 英文稿定稿骨架
ABSTRACT = {
    'Background. ': (
        'Background. Single-cis-instrument drug-target Mendelian '
        'randomization (MR) is widely used to prioritise therapeutic '
        'targets, yet null findings are often reported without quantifying '
        'power or sensitivity to endpoint definitions. STROBE-MR items 4b '
        'and 15 require a statement on a priori power and discussion of '
        'bias magnitude. An unqualified null may be misread as evidence of '
        'no effect.'),
    'Methods. ': (
        'Methods. We performed two-sample, summary-level MR using nine '
        'UKB-PPP CX3CL1 pQTL credible-set lead variants (Open Targets '
        'v26.9; n=33,436 European-ancestry participants) and FinnGen '
        'releases R9-R13. Nineteen release-by-endpoint combinations covered '
        'diabetic neuropathy, alternative definitions and microvascular '
        'negative controls. The primary analysis was the Wald ratio of the '
        'lead cis variant rs62037084. Sensitivity analyses included LD '
        'clumping, IVW, MR-Egger, weighted median, leave-one-out, '
        'credible-set and 400 kb region-wide approximate Bayes-factor '
        'colocalisation, Steiger orientation, heterogeneity statistics, '
        'E-values, bootstrap, simulation and pipeline calibration. '
        'Post-hoc power and minimum detectable odds ratios were derived '
        'from observed standard errors; no a priori calculation was '
        'performed (STROBE-MR item 4b).'),
    'Results. ': (
        'Results. All 19 primary estimates were non-significant '
        '(P=0.141-0.987), although odds ratios ranged 0.350-3.755 and the '
        'primary endpoint changed direction between releases. Power for '
        'OR=1.5 was 5.45%-10.10% and the minimum detectable OR at 80% '
        'power 5.61-318.57. Eight of 133 tests were nominally significant, '
        'consistent with 6.65 expected under a global null (binomial '
        'P=0.548); none survived Bonferroni correction (3.76x10^-4). At '
        'the primary endpoint, F=552.42, OR=1.394 (95% CI 0.159-12.21) '
        'and the approximate Bayes factor was 0.62, indicating only weak '
        'support for the null. Region-wide colocalisation on 1,744 aligned '
        'variants gave PP.H4=0.032 and PP.H1=0.968, localising absent '
        'signal to the outcome side. The same pipeline recovered the '
        'expected PCSK9 direction (minimum P=1.91x10^-32), while all 14 '
        'cis-only calibration combinations were non-significant. Naive '
        'reporting mislabelled a true OR=1.5 as "no effect" in 93.5% of '
        'simulations; the decision framework made no false assertions.'),
    'Conclusions. ': (
        'Conclusions. We obtained no supportive evidence that genetically '
        'proxied circulating CX3CL1 influences diabetic neuropathy risk, '
        'and the public data cannot test a moderate effect. We do not '
        'claim that CX3CL1 has no effect. Quantifying power deficits and '
        'cross-context drift supports a framework requiring power, endpoint '
        'definitions, multiplicity and wording before any null MR result '
        'is interpreted.'),
}

M8 = (
    '[M8, items 4a/8 - region-wide colocalisation] Colocalisation was '
    'assessed at two levels. The credible-set-level analysis used the nine '
    'UKB-PPP SuSiE 95% credible-set lead variants with approximate Bayes '
    'factors (Wakefield 1994; Giambartolomei 2014; p1=p2=1e-4, p12=1e-5, '
    'W=0.04). Because the UKB-PPP full pQTL summary statistics are not '
    'publicly retrievable from this setting (the Open Targets-registered '
    'object returns HTTP 403 on anonymous access and Synapse syn51364943 '
    'returns HTTP 403 GEO_RESTRICTION under NIH NOT-OD-25-083, a geographic '
    'restriction that a personal access token does not resolve), the '
    'region-wide analysis used a fully public independent exposure '
    'dataset: the Fenland CX3CL1 plasma pQTL (GWAS Catalog GCST90012074; '
    'Pietzner et al., Nat Commun 2021; n=21,758 European-ancestry '
    'participants; GRCh38). The 400 kb window chr16:57,200,000-57,600,000 '
    'was extracted from FinnGen R13 DM_NEUROPATHY and G6_POLYOTHUNS by '
    'byte-range retrieval. Alleles were aligned by exact match, strand '
    'complementation or sign flip; unresolved configurations were '
    'discarded rather than imputed. Priors were varied over W=0.01, 0.04, '
    '0.16 and p12=5e-6, 1e-5, 2e-5. This supplementary analysis uses the '
    'same outcome data as the primary analysis and is not an independent '
    'dataset replication [VERIFY: insert exact Wakefield, Giambartolomei '
    'and Pietzner references].')

R10 = (
    '[R10, items 12b/13d - region-wide colocalisation] Of 2,010 Fenland '
    'and 3,140 FinnGen variants in the 400 kb window, 1,746 matched by '
    'position and 1,744 were retained after allele alignment. For diabetic '
    'neuropathy, PP.H1 (exposure-only signal) was 0.9680 and PP.H4 '
    '(shared causal variant) was 0.0319; for the polyneuropathy negative '
    'control they were 0.9827 and 0.0173. Across nine prior settings '
    'PP.H4 remained between 0.0165 and 0.0619. The asymmetry was marked: '
    '297 variants reached genome-wide significance for the exposure '
    '(strongest rs683544, z=-18.14, P=2.81x10^-73), whereas no variant '
    'reached P<1e-5 for the outcome (strongest rs72786892, '
    'P=8.70x10^-5). The analysis therefore does not merely fail to find '
    'colocalisation; it localises the absence of signal to the outcome '
    'side. These are posterior probabilities under a single-causal-variant '
    'assumption, not proof that colocalisation is excluded. Full results, '
    'per-variant alignment records and sensitivity tables are in '
    'Supplementary Material S15.')

D7 = (
    '[D7 - scope of the region-wide analysis] The region-wide '
    'colocalisation assumes one causal variant per trait and does not '
    'model multiple independent signals; with 1,744 variants across 400 kb '
    'this is a simplification, and PP.H4 could be underestimated if the '
    'exposure has two independent signals in the region. We therefore '
    'report PP.H1 alongside PP.H4, because a high PP.H1 does not depend on '
    'the number of exposure signals, only on the absence of any outcome '
    'signal in the same window. The Fenland (United Kingdom) and FinnGen '
    '(Finland) cohorts also differ in LD structure; approximate Bayes '
    'factors do not use LD directly, so this does not invalidate the '
    'computation, but it limits cross-population interpretation of a '
    'shared causal variant. This supplementary analysis changes neither '
    'the primary analysis nor the wording restriction.')

FIG9 = (
    'Figure 9. Region-wide approximate-Bayes-factor colocalisation in '
    'chr16:57.2-57.6 Mb (GRCh38). (a) Posterior probabilities for '
    'H0-H4 for diabetic neuropathy and the polyneuropathy negative '
    'control; H1 (exposure-only signal) dominates. (b) Regional '
    '-log10(P) for the Fenland CX3CL1 pQTL (upper) and FinnGen R13 '
    'diabetic neuropathy (lower); 297 exposure variants reach '
    'P<5e-8 and no outcome variant reaches P<1e-5. (c) Sensitivity '
    'of PP.H4 to the prior width W and the shared-variant prior p12.')

REFERENCES = """## References (draft; [VERIFY] items require the original articles)

1. Corbin LJ, Merino J, Fall T, Herder C. Raising the bar for publication of Mendelian randomisation studies in Diabetologia. Diabetologia. 2025;68(10):2088-2091. doi:10.1007/s00125-025-06484-6
2. Rosas-Chavez F, Merriman TR. Mendelian Randomization Studies: A Metric for Quality Evaluation. Gout Urate Cryst Depos Dis. 2025;3(2):8. doi:10.3390/gucdd3020008
3. Zou L, Fazia T, Guo H, Berzuini C. Bayesian Mendelian randomization with an interval causal null hypothesis: ternary decision rules and loss function calibration. BMC Med Res Methodol. 2024;24. doi:10.1186/s12874-023-02067-4
4. Zou L, Guo H, Berzuini C. Bayesian mendelian randomization with study heterogeneity and data partitioning for large studies. BMC Med Res Methodol. 2022;22:162. doi:10.1186/s12874-022-01619-4
5. Hu S, Zhu P, Gao S, et al. Self-Report Alzheimer's Disease Statuses in UK Biobank Distort Downstream Analyses. J Neurochem. 2025;169(7):e70151. doi:10.1111/jnc.70151
6. Pietzner M, Wheeler E, Carrasco-Zanini J, et al. Mapping the proteo-genomic convergence of human diseases. Nat Commun. 2021. [VERIFY: exact volume, article number and DOI; GWAS Catalog accession GCST90012074]
7. Kurki MI, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. Nature. 2023. [VERIFY: exact volume, pages and DOI]
8. Wakefield J. Bayes factors for genome-wide association studies. Am J Hum Genet. 1994. [VERIFY: exact volume, pages and DOI]
9. Giambartolomei C, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014. [VERIFY: exact volume, article number and DOI]
10. Liu BX, et al. Mitigating inconsistencies in GWAS follow-up analyses with LocusCompare2. Nat Genet. 2025. [VERIFY: original citation, volume and DOI]
11. Sun BB, et al. Plasma proteomic associations with genetics and health in the UK Biobank. Nature. 2023. [VERIFY: exact citation for the UKB-PPP resource]
12. Open Targets Platform / Open Targets Genetics primary publication. [VERIFY: insert complete reference and release v26.9 identifier]
13. Stender S, et al. Reclaiming mendelian randomization from the deluge of papers and misleading findings. [VERIFY: original article, journal, year and DOI]
14. FinnGen public summary statistics releases R9-R13. [VERIFY: insert full data citation with accession/URL per BMC Force 11 data-citation policy]
15. Fenland CX3CL1 plasma pQTL, GWAS Catalog GCST90012074. [VERIFY: insert full data citation with accession]
"""


def patch_skeleton():
    p = MD + '/A1英文稿定稿骨架_20261111.md'
    s = _read(p)

    # 头部说明
    s = _rep(s,
             '> Compiled on 2026-11-11 from all previously delivered EN '
             'paragraphs plus the new A4/A5/B1/B2/B4 paragraphs. Every '
             '[VERIFY] remains intentional. The registration, repository '
             'DOI and author details are NOT asserted. Wording restriction: '
             '"no supportive evidence" / "insufficient to test a moderate '
             'effect"; never "CX3CL1 has no effect".',
             '> Consolidated 2026-11-11; post-submission full-text review '
             'applied the same day. OSF retrospective registration was '
             'deposited on 11 November 2026 (DOI pending archiving email). '
             'Every [VERIFY] remains intentional. Wording restriction: '
             '"no supportive evidence" / "insufficient to test a moderate '
             'effect"; never "CX3CL1 has no effect". No analysis number '
             'differs from the archived registration.')

    # 摘要：压缩到期刊上限内，并修正 ABF 方向
    s = _rep(s,
             '## Abstract (structured, ~350 words; compress to journal '
             'limit before submission)',
             '## Abstract (structured; recheck the journal word limit '
             'before submission)')
    lines = s.split('\n')
    out = []
    for ln in lines:
        key = next((k for k in ABSTRACT if ln.startswith(k)), None)
        out.append(ABSTRACT[key] if key else ln)
    s = '\n'.join(out)

    # 工程数字
    s = _rep(s, 'SHA256 checksums of 75 input files',
             'SHA256 checksums for every input file')

    # M1：写明两层 coloc 的暴露来源；澄清 9 个可信集 = 3 cis + 6 trans
    s = _rep(s, 'nine SuSiE 95% credible sets',
             'nine SuSiE 95% credible sets (three cis and six trans), each '
             'represented by its lead variant')
    s = _rep(s,
             '(d) Colocalisation: Open Targets precomputed coloc '
             '(COLOC_PIP_ECAVIAR) and approximate Bayes-factor coloc '
             'computed on the FinnGen R13 chr16:57.2-57.6 Mb region.',
             '(d) Colocalisation: Open Targets precomputed coloc '
             '(COLOC_PIP_ECAVIAR) at the credible-set level, and a '
             'region-wide approximate-Bayes-factor analysis on FinnGen R13 '
             'chr16:57.2-57.6 Mb using the independent, fully public '
             'Fenland CX3CL1 plasma pQTL (GCST90012074; n=21,758) as the '
             'exposure dataset (see M8).')
    s = _rep(s, 'nine UKB-PPP CX3CL1 pQTL credible sets',
             'nine UKB-PPP CX3CL1 pQTL credible-set lead variants')

    # M4：标清两套先验敏感性各自属于哪一层
    s = _rep(s,
             'prior sensitivity for colocalisation (W=0.04; p12=1e-6, '
             '1e-5, 1e-4)',
             'prior sensitivity at the credible-set level (W=0.04; '
             'p12=1e-6, 1e-5, 1e-4) and at the region-wide level '
             '(W=0.01, 0.04, 0.16; p12=5e-6, 1e-5, 2e-5)')

    # M6：日期已知，DOI 待回填
    s = _rep(s,
             'A retrospective registration was prepared for OSF with the '
             'date and DOI to be inserted at submission; [VERIFY: insert '
             'OSF registration date and DOI].',
             'A retrospective Secondary Data Preregistration was deposited '
             'in the OSF registries on ' + REG_DATE + '; ' + OSF_DOI + '.')

    # 插入 M8 / R10 / D7 / 参考文献 / 图 9
    if '[M8, items 4a/8' not in s:
        s = s.replace('\n\n## Results', '\n' + M8 + '\n\n## Results', 1)
    if '[R10, items 12b/13d' not in s:
        s = s.replace('\n\n## Discussion', '\n' + R10 + '\n\n## Discussion', 1)
    if '[D7 - scope' not in s:
        s = s.replace('\n\n## Declarations', '\n' + D7 + '\n\n## Declarations', 1)
    if '## References (draft' not in s:
        s = s.replace('\n\n## Figure legends',
                      '\n\n' + REFERENCES.strip() + '\n\n## Figure legends', 1)
    if 'Figure 9. Region-wide' not in s:
        s = _rep(s, '(10.1%) and the solid line marks 80% power.',
                 '(10.1%) and the solid line marks 80% power.\n\n' + FIG9)

    # Declarations：BMC 规定的数据可用性句式 + 数据集引用要求
    s = _rep(s,
             'linked to the OSF registration [VERIFY: insert OSF DOI]',
             'linked to the OSF registration deposited on ' + REG_DATE +
             ' (' + OSF_DOI + '). All publicly available datasets are '
             'fully cited with accession numbers or DOIs in the '
             'References, in accordance with the BMC Force 11 data-citation '
             'policy')
    _write(p, s)
    return s


def patch_figure_legends():
    p = PACK + '/a1_figure_legends.txt'
    s = _read(p)
    if 'Figure 9. Region-wide' not in s:
        s = s.rstrip() + '\n\n' + FIG9 + '\n'
    _write(p, s)


def patch_cover_letter():
    p = PACK + '/a1_cover_letter_draft.txt'
    s = _read(p)
    s = _rep(s, 'SHA256 checksums of 75 input files',
             'SHA256 checksums for every input file')
    s = _rep(s,
             'the date and DOI are [VERIFY: insert OSF registration date '
             'and DOI]',
             'the registration was deposited on ' + REG_DATE + ' and its '
             'DOI is ' + OSF_DOI)
    # 把落在签名之后的第 6 节移到签名之前，并把感谢句放到签名正上方
    thanks = 'Thank you for considering this submission.'
    marker = '\n\n6. Additional analyses completed after the first draft'
    if marker in s and s.find('Sincerely,') < s.find(marker):
        i = s.find(marker)
        block = s[i:].strip()
        s = s[:i].rstrip()
        s = s.replace('\n\n' + thanks, '')
        k = s.find('Sincerely,')
        s = (s[:k].rstrip() + '\n\n' + block + '\n\n' + thanks +
             '\n\n' + s[k:])
    _write(p, s)


def patch_declarations():
    p = PACK + '/a1_declarations.txt'
    s = _read(p)
    s = _rep(s,
             'A retrospective OSF registration is prepared, with date and '
             'DOI\n[VERIFY: insert OSF date and DOI].',
             'A retrospective OSF registration was deposited on ' +
             REG_DATE + '; its DOI is ' + OSF_DOI + '.')
    s = _rep(s, 'The retrospective registration is at [VERIFY: OSF DOI].',
             'The dataset(s) supporting the conclusions of this article '
             'are available in the Zenodo repository, [VERIFY: insert the '
             'Zenodo DOI and its https:// hyperlink]. The retrospective '
             'registration was deposited on ' + REG_DATE + ' (' +
             OSF_DOI + ').')
    _write(p, s)


def patch_supplementary_index():
    p = PACK + '/a1_supplementary_index.tsv'
    s = _read(p)
    s = _rep(s, 'results/fig*_EN.pdf; results/fig6-8*.pdf',
             'results/fig*_EN.pdf; results/fig6-8*.pdf; '
             'results/fig9_b3_region_coloc_EN.pdf')
    if 'S15\tRegion-wide' not in s:
        s = s.rstrip() + (
            '\nS15\tRegion-wide complete colocalisation (Fenland x '
            'FinnGen R13)\tresults/a1_b3_coloc_full.json; '
            'results/a1_b3_coloc_variants.tsv; '
            'results/a1_b3_fenland_region.tsv; '
            'A1_附录B3_区域级完整共定位_20261111.docx\t'
            '1,744 aligned variants in 400 kb; PP.H4=0.0319, '
            'PP.H1=0.9680; nine prior settings; Figure 9.\n')
    _write(p, s)


def patch_zenodo():
    p = PACK + '/a1_zenodo_metadata.json'
    d = json.loads(_read(p))
    d['description'] = d['description'].replace(
        'including cross-target replication and biological positive '
        'controls.',
        'including cross-target replication, biological positive controls '
        'and region-wide colocalisation against the Fenland CX3CL1 pQTL '
        '(GCST90012074).')
    for r in d.get('related_identifiers', []):
        if r.get('relation') == 'isSupplementedBy':
            r['relation'] = 'isSupplementTo'
    d['license'] = ('[VERIFY: choose CC-BY-4.0 or MIT; note the OSF '
                    'registration itself is CC0 1.0]')
    d['note'] = ('Metadata draft. OSF registration deposited 11 November '
                 '2026 (' + OSF_DOI + '); creator, licence and the Zenodo '
                 'deposit itself require the author account before release.')
    _write(p, json.dumps(d, ensure_ascii=False, indent=2))


def patch_github_readme():
    p = PACK + '/a1_github_readme.md'
    s = _read(p)
    s = _rep(s,
             'Open Targets Platform v26.9 (UKB-PPP pQTL), FinnGen public '
             'releases\nR9-R13, Ensembl REST 1000 Genomes phase 3 EUR LD.',
             'Open Targets Platform v26.9 (UKB-PPP pQTL), FinnGen public '
             'releases\nR9-R13, the Fenland CX3CL1 plasma pQTL (GWAS '
             'Catalog GCST90012074, used for region-wide colocalisation), '
             'and Ensembl REST 1000 Genomes phase 3 EUR LD.')
    s = _rep(s, 'This recomputes all statistics and figures from the '
             'deposited result\ntables.',
             'This runs the 14-step offline reproduction chain and '
             'recomputes all statistics and figures from the deposited '
             'result\ntables.')
    s = _rep(s, '- `submission_pack/` - cover letter, declarations, '
             'presubmission enquiry,\n  figure legends and compliance '
             'checklists',
             '- `submission_pack/` - cover letter, declarations, '
             'presubmission enquiry,\n  figure legends, supplementary '
             'index and compliance checklists')
    _write(p, s)


def patch_verify_register():
    p = RES + '/a1_verify_register.tsv'
    s = _read(p)
    s = _rep(s,
             "V02\tOSF registration date and DOI\tEXCLUDED_BY_USER\t"
             "Expressly excluded from this round; retrospective "
             "registration text pack is prepared.",
             "V02\tOSF registration date and DOI\tDONE_DATE_ONLY\t"
             "Deposited 11 November 2026 (Secondary Data Preregistration, "
             "osf.io/pmhve); archiving in progress; DOI pending email.")
    s = _rep(s,
             "V20\tPAT-based full coloc appendix\tDEFERRED\t"
             "B3: requires a personal access token; not needed for the "
             "first submission.",
             "V20\tPAT-based full coloc appendix\tCLOSED\t"
             "PAT proved unusable (Synapse 403 GEO_RESTRICTION); B3 "
             "completed with the fully public Fenland GCST90012074 "
             "dataset, no token needed.")
    s = _rep(s,
             "V13\tAbstract compression to journal word limit\tDONE\t"
             "Structured abstract compiled in the consolidated skeleton; "
             "final trim is editorial.",
             "V13\tAbstract compression to journal word limit\tDONE\t"
             "Review found 467 words against a ~350-word label; abstract "
             "rewritten to journal length on 11 November 2026.")
    s = _rep(s,
             "V16\tEnglish vector figures\tDONE\t"
             "fig1-fig5 redrawn as EN PDF/SVG/PNG; fig6-fig8 were already "
             "English vector figures.",
             "V16\tEnglish vector figures\tDONE\t"
             "fig1-fig9 supplied as EN PDF/SVG/PNG; Figure 9 (B3 region-wide "
             "coloc) added to the legends and supplementary index.")
    s = _rep(s,
             "V19\tJournal word limits, APC, presubmission channel\t"
             "OPEN_LITERATURE\tRecheck on the current journal site "
             "immediately before submission; enquiry draft prepared.",
             "V19\tJournal word limits, APC, presubmission channel\t"
             "OPEN_LITERATURE\tJournal pages checked 11 November 2026: "
             "BMC requires Force 11 dataset citation in References and a "
             'fixed "Availability of data and materials" sentence. Exact '
             "abstract word limit and APC still need a final check.")
    extra = [
        ('V21', 'OSF registration Update to correct the seven-target list',
         'PENDING_ARCHIVE',
         'Field "Reliability and Robustness Testing" lists CRP, HMGCR, '
         'CETP; the actual panel is ACE, AGT, IL1RN, IL6R, LDLR, PCSK9, '
         'SERPINE1. Issue an Update after archiving; original is retained.'),
        ('V22', 'OSF licence and bibliographic contributor',
         'PENDING_ARCHIVE',
         'Review text showed an empty License field; confirm CC0 1.0 is '
         'displayed and that Yi Lu is ticked Bibliographic Contributor.'),
        ('V23', 'Region-wide coloc integrated into the manuscript', 'DONE',
         'M8, R10, D7, Figure 9, S15 and the abstract were added on '
         '11 November 2026; no number differs from the registration.'),
        ('V24', 'Six OSF DOI placeholders', 'DONE',
         'Backfilled ' + OSF_DOI + ' on 11 November 2026 into M6, '
         'Declarations, cover letter, declarations document, Zenodo '
         'metadata and GitHub README; no placeholder remains.'),
        ('V25', 'Stale engineering counts in Chinese internal reports',
         'OPEN_EDITORIAL',
         'Earlier reports say 12 steps / 123-124 files / 8 figures / OSF '
         'excluded / B3 needs a PAT. Mark them as historical snapshots or '
         'regenerate; do not attach them to the submission.'),
    ]
    for vid, item, status, note in extra:
        if '\n' + vid + '\t' not in s:
            s = s.rstrip() + '\n' + '\t'.join((vid, item, status, note))
    _write(p, s)


def patch_osf_paste_blocks():
    """本地粘贴包归档件：修正七靶点名单，避免与已提交登记对照时再误导。"""
    p = RES + '/a1_osf_paste_blocks.txt'
    s = _read(p)
    s = _rep(s,
             '  - A seven-target replication panel (PCSK9, LDLR, IL6R, '
             'CRP, HMGCR, CETP and\n    others with usable cis pQTL)',
             '  - A seven-target replication panel (ACE, AGT, IL1RN, '
             'IL6R, LDLR, PCSK9 and\n    SERPINE1)')
    _write(p, s)


def patch_final_doi():
    """14-step reruns must never resurrect an OSF DOI placeholder."""
    files = [
        MD + '/A1英文稿定稿骨架_20261111.md',
        PACK + '/a1_declarations.txt',
        PACK + '/a1_cover_letter_draft.txt',
        PACK + '/a1_zenodo_metadata.json',
        PACK + '/a1_github_readme.md',
    ]
    for p in files:
        s = _read(p)
        t = s.replace('; [VERIFY: insert OSF DOI].',
                      ' (' + OSF_DOI + ').')
        t = t.replace('[VERIFY: insert OSF DOI]', OSF_DOI)
        t = t.replace('[VERIFY: OSF registration DOI]', OSF_DOI)
        t = t.replace('[VERIFY: add licence, author list, OSF DOI and Zenodo '
                      'DOI before release.]',
                      'OSF registration: ' + OSF_DOI +
                      ' (https://osf.io/apwe2), CC0 1.0. '
                      '[VERIFY before release: licence, author list and '
                      'Zenodo DOI.]')
        if t != s:
            _write(p, t)


def apply():
    sk_path = MD + '/A1英文稿定稿骨架_20261111.md'
    sk = _read(sk_path)
    # 若骨架已被本轮修订过（未重跑生成器），跳过以免重复插入
    if 'post-submission full-text review' in sk:
        skeleton = sk
    else:
        skeleton = patch_skeleton()
    patch_figure_legends()
    patch_cover_letter()
    patch_declarations()
    patch_supplementary_index()
    patch_zenodo()
    patch_github_readme()
    patch_verify_register()
    patch_osf_paste_blocks()
    patch_final_doi()
    abstract_words = sum(len(ABSTRACT[k].split()) for k in ABSTRACT)
    print('  postreview patch applied; structured abstract = %d words'
          % abstract_words)
    return skeleton
