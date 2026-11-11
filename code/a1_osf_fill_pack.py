# -*- coding: utf-8 -*-
"""A1 OSF 回顾性登记：逐字段可粘贴包生成器

输入：OSF "Secondary Data Preregistration" 模板（schema 64775783798e08000a70407e）
的 27 个字段，已通过 api.osf.io 核实字段名与顺序。

输出：
  results/a1_osf_field_mapping.tsv   —— 字段 / 中文提示 / 状态 / 英文值（可复核）
  results/a1_osf_paste_blocks.txt    —— 纯文本粘贴块（避免从 Word 复制带入格式）

原则：只写已有真实结果与已核实事实；不确定处一律 [VERIFY]，不填推测数字。
"""
import csv
import io
import os

RES = os.path.join(os.path.dirname(__file__), '..', 'results')
TODAY = '2026-11-11'

# 字段顺序与名称均取自 OSF Secondary Data Preregistration 模板（v3, 5 页）
FIELDS = []


def f(page, name, ctype, cn_hint, en, status='READY'):
    FIELDS.append(dict(page=page, field=name, input_type=ctype,
                       cn_hint=cn_hint, en=en.strip(), status=status))


# ---------------------------------------------------------------- 第 1 页
f('Study Information', 'Research questions', 'long text',
  '直接粘贴。三个问题对应论文三条主线：能支持多强的结论、结论漂移多少、该怎样报告。',
  """
RQ1. In a drug-target two-sample Mendelian randomisation (MR) that has a single
independent cis instrument and no strong signal in the outcome region, what
strength of conclusion do publicly available summary statistics support?

RQ2. Within one exposure-outcome pair, how much do point estimates and the verbal
label attached to them drift across outcome-definition variants, across FinnGen
data releases, and across instrument-set choices?

RQ3. Which reporting rules prevent an underpowered null from being written as
evidence of no effect, and can such a rule set be stated as a transferable
decision procedure that other drug-target MR studies can apply?
""")

f('Study Information', 'Hypotheses', 'long text',
  '本页最敏感。必须写明这些是“初步结果已知之后形成的量化预判”，不是事前假设，否则构成 HARKing。照抄即可。',
  """
IMPORTANT FRAMING. The statements below are quantitative predictions, not
pre-data hypotheses. They were formulated after the exploratory analyses listed
under "Prior knowledge" had already been completed and after the first-round
result was known. We register them so that the confirmatory re-run is auditable,
and we explicitly do not claim that they were specified before the data were seen.
Any reader should treat them as falsifiable predictions whose evidential value is
limited by this ordering.

P1. The primary Wald ratio for the lead cis instrument will not reach the
Bonferroni-corrected threshold, and its 95% confidence interval will include the
null.

P2. Point estimates for the same exposure-outcome pair will vary across FinnGen
releases and across outcome-definition variants by substantially more than the
width of a single estimate implies, and at least one pair of releases will differ
in the sign of the point estimate.

P3. Post-hoc power at a true odds ratio of 1.5, derived from the observed standard
error, will be low in most release-endpoint combinations, and the minimum
detectable odds ratio at 80% power will exceed 2 in most of them.

P4. The number of nominally significant (P < 0.05) results across all
combination-by-method tests will not exceed the number expected under the global
null by more than chance would produce.

P5. A pre-specified reporting rule set applied to these results will classify no
combination as supporting a causal effect and will classify most combinations as
inconclusive rather than negative.
""")

# ---------------------------------------------------------------- 第 2 页
f('Data Description', 'Datasets used', 'long text',
  '逐个数据集写清来源、版本、人群、样本量口径。暴露侧 n=33,436 已核实（OT nSamples，2026-11-11 存档）。',
  """
EXPOSURE. Plasma/serum protein quantitative trait loci (pQTL) for CX3CL1
(UniProt P78423; Open Targets OID20976) from the UK Biobank Pharma Proteomics
Project (UKB-PPP), European ancestry, accessed through the Open Targets Platform
GraphQL API. Open Targets study identifier:
UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1. Nine variants were used from the SuSiE 95%
credible sets for this protein; eight of them could be retrieved in FinnGen and
one indel (16_57375150_AAAG_A) was absent from the FinnGen summary statistics and
is therefore excluded from all estimates rather than imputed. Per-protein sample
size: n = 33,436, all of European ancestry (Open Targets study record
UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1; nSamples = 33,436; ldPopulation = nfe with
relativeSampleSize = 1.0; verified through the Open Targets Platform GraphQL API on
2026-11-11 and archived as data/ot_cx3cl1_study_meta.json in the analysis
repository).

REGION-WIDE COLOCALISATION EXPOSURE. Because the UKB-PPP full pQTL summary
statistics are not retrievable from this setting (Open Targets-registered object
returns HTTP 403 on anonymous access; Synapse syn51364943 returns HTTP 403
GEO_RESTRICTION under NIH NOT-OD-25-083, a geographic restriction that a personal
access token does not resolve), region-wide colocalisation uses an independent
and fully public exposure dataset: the Fenland CX3CL1 plasma pQTL, GWAS Catalog
accession GCST90012074 (Pietzner et al., Nat Commun 2021, PMID 33067605), n =
21,758 European-ancestry participants, UK Biobank-free, GRCh38, downloaded from
EBI FTP without application. The analytical window is chr16:57,200,000-57,600,000
(GRCh38). Using a cohort from a different country than the outcome cohort also
makes the exposure and outcome samples structurally independent.

OUTCOMES. FinnGen public summary statistics, releases R9, R10, R11, R12 and R13.
Endpoints: DM_NEUROPATHY (primary), DM_POLYNEURO, G6_POLYOTHUNS, DM_RETINOPATHY
(STRICT and EXMORE definitions where available), DM_NEPHROPATHY (EXMORE where
available). Nineteen release-by-endpoint combinations in total. Case and control
counts were taken from the release manifests where these are public; releases R11
and R13 publish no machine-readable manifest, so counts for those combinations are
recorded as not available rather than estimated.

LD REFERENCE. 1000 Genomes Phase 3 European ancestry, accessed through the Ensembl
REST API, used for clumping and for the variant-level LD matrix.

COLOCALISATION INPUT. Open Targets credible-set and locus-level data for the
CX3CL1 region, used for ABF colocalisation between the pQTL and the outcome.

The exposure and outcome samples are non-overlapping: UKB-PPP participants are not
part of FinnGen.
""")

f('Data Description', 'Data availability', 'single select',
  '下拉选第一项：The dataset is publicly available。',
  """
The dataset is publicly available
""")

f('Data Description', 'Data access', 'short text',
  '写清“无需申请”与 FinnGen 条款要求，避免被质疑合规。',
  """
All summary statistics are downloadable without an application. FinnGen summary
statistics were obtained from the public Google Cloud Storage bucket and are used
under the FinnGen public data-use terms; no individual-level FinnGen data were
accessed. UKB-PPP pQTL results are obtained as aggregate results through the Open
Targets Platform; no individual-level UK Biobank data were accessed, so no UK
Biobank application was required for this work. 1000 Genomes and Ensembl
resources are open.
""")

f('Data Description', 'Data identifiers', 'short text',
  '这套标识符要能让审稿人不经联系你就复现取数。注意 R12 的 URL 路径与其他 release 不同。',
  """
Exposure study: UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1 (Open Targets Platform).
Lead cis variant: rs62037084 (GRCh38 16:57379748 C>T). Second cis variant:
rs223872. All nine credible-set variants and their identifiers are listed in the
endpoint/variant dictionary deposited with this registration.

Outcome URL pattern, FinnGen public bucket (all accessed 2026-11-10).
Replace {ENDPOINT} with the FinnGen endpoint code, for example DM_NEUROPATHY:
  R9   https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_{ENDPOINT}.gz
  R10  https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_{ENDPOINT}.gz
  R11  https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_{ENDPOINT}.gz
  R12  https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{ENDPOINT}.gz
  R13  https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{ENDPOINT}.gz
Note the additional /release/ path segment in R12. Tabix indices are fetched
from the same bucket.
""")

f('Data Description', 'Access date', 'short text',
  '填实际取数日。这个日期进审计轨迹，与文件时间戳一致（2026-11-10）。',
  """
2026-11-10 (all FinnGen releases, Open Targets and Ensembl queries). A second
retrieval for the multi-target replication panel was performed on 2026-11-10 as
well; the retrieval log with per-request status is deposited with this
registration.
""")

f('Data Description', 'Data collection procedures', 'long text',
  '说明定点取数方式——这是本文可复核性的技术要点。',
  """
Summary statistics were not downloaded in full. For each release-endpoint
combination, the tabix index was retrieved first and used to resolve the byte
offsets of the CX3CL1 region (GRCh38 chr16:57.35-57.42 Mb, +/-4 Mb window), and
only those byte ranges were requested over HTTP from the BGZF-compressed
summary-statistics file. A custom BGZF reader was used because the htslib binary
was not available in the analysis environment; the reader parses the tabix index
header directly. Each of the 133 region requests returned status OK.

Variants were harmonised on REF/ALT alleles against the exposure file. Exact
allele matches were used where available. One variant (10_63389329_G_A) was
present in FinnGen only at a +1 bp positional offset with reversed allele
orientation; it was matched at that offset, flagged as an approximate match and
reported as such rather than silently accepted. One indel (16_57375150_AAAG_A)
could not be retrieved and was excluded, never imputed. Variants absent from a
release were recorded as excluded, never imputed. The retrieval log records, per
combination, the number of raw lines, exact matches, approximate matches and
missing variants.
""")

f('Data Description', 'Data collection procedures documentation', 'file upload',
  '上传取数脚本两个（见附件策略），不要上传整包代码。',
  """
Upload: a1_tabix_fetch.py and a1_multi_fetch.py
""")

f('Data Description', 'Codebook', 'long text',
  '端点字典的字段释义。R11/R13 无公开 manifest 必须写明“记为空值”而非估算。',
  """
The deposited endpoint dictionary has one row per release-endpoint combination
with the following fields: release; endpoint (FinnGen endpoint code as published);
endpoint_type (primary, negative-control, or alternative-definition); n_cases;
n_controls; endpoint_name (plain-language label); case_control_source (manifest
file name, or NO_PUBLIC_MANIFEST); url.

Four endpoint roles are distinguished and stated in advance: DM_NEUROPATHY is the
primary outcome; DM_POLYNEURO and G6_POLYOTHUNS index the same clinical domain
under different coding; DM_RETINOPATHY and DM_NEPHROPATHY (STRICT / EXMORE
variants) serve as within-domain alternative definitions and as
negative-control outcomes, because CX3CL1 signalling is not expected to act
specifically on retinal or renal microvascular disease.

Where a release publishes no machine-readable manifest (R11, R13), n_cases and
n_controls are left empty and flagged, and no figure is inferred from press
releases or secondary sources.
""")

f('Data Description', 'Codebook documentation', 'file upload',
  '上传端点字典 TSV。',
  """
Upload: a1_endpoint_dict.tsv
""")

# ---------------------------------------------------------------- 第 3 页
f('Variables', 'Manipulated variables', 'long text',
  '本文无实验操纵。要写清“暴露由遗传变异代理”，否则审稿人会认为填错模板。',
  """
None. This is an analysis of observational summary statistics; no variable was
experimentally manipulated. The exposure is proxied by germline genetic variants
(one cis pQTL instrument in the primary analysis), which is the defining feature
of MR rather than an intervention. Four nested instrument sets are defined in
advance to separate instrument-strength effects from endpoint effects:
  S1  the single lead cis variant (primary analysis);
  S2  both cis variants in the credible set;
  S3  cis variants after LD clumping at r2 < 0.001;
  S4  all variants with a reported association, including trans variants
      (sensitivity only, see the assumption-violation field).
""")

f('Variables', 'Manipulated variables documentation', 'file upload',
  '本字段留空即可，不必上传文件。',
  """
Not applicable: no manipulated variables. No file uploaded.
""")

f('Variables', 'Measured variables', 'long text',
  '列出所有测量与派生量，并给出单位。审稿人据此判断 OR 的量纲。',
  """
Per variant and per release-endpoint combination: effect allele and other allele;
effect-allele frequency; beta and standard error of the variant on circulating
CX3CL1 (exposure, natural-log scale, per allele); beta and standard error of the
variant on the outcome (log odds ratio of the endpoint, per allele); the
corresponding P values; and, where public, case and control counts.

Derived quantities: Wald ratio beta and standard error; odds ratio with 95%
confidence interval; post-hoc power at a true odds ratio of 1.5; minimum
detectable odds ratio at 80% power; Cochran Q and I-squared across instruments;
I-squared for instrument strength (I2_GX); F statistic for the lead instrument;
MR-Egger intercept and its standard error; Steiger direction test; E-value;
approximate Bayes factor from colocalisation; and, for the replication panel,
the same quantities for seven additional drug targets.
""")

f('Variables', 'Missing data', 'long text',
  '必须写明“不插补”并点名缺失变异。这是可复核性的关键。',
  """
No missing data are imputed anywhere in the analysis. Three explicit cases:
  1. The indel 16_57375150_AAAG_A is present in the exposure credible set but
     absent from the FinnGen summary statistics; it is excluded from every
     estimate and reported as excluded.
  2. The variant 10_63389329_G_A is absent from FinnGen R9 and present in later
     releases, so R9 estimates use one fewer variant by construction.
  3. One variant (10_63389329_G_A) required a +1 bp approximate match and is
     flagged as such rather than treated as an exact match.
Case and control counts are missing for R11 and R13 because those releases
publish no machine-readable manifest; the affected combinations are analysed but
marked as not evaluable for any statistic that requires an outcome sample size
(Steiger test, sample-size-dependent approximations).
""")

f('Variables', 'Unit of analysis', 'long text',
  '一句话：变异-端点-release 组合。并声明不是独立样本。',
  """
One unit is a variant-by-endpoint-by-release combination, and the inferential unit
for the primary analysis is a single Wald ratio (one instrument, one outcome,
one release). The 19 combinations are not independent observations: they share the
same exposure, the same instrument, overlapping participants and overlapping
endpoint definitions. The 76 instrument-set estimates are a nested sensitivity
grid, not 76 independent tests. No result is selected on the basis of its P value.
""")

f('Variables', 'Statistical outliers', 'long text',
  '说明不剔离群值，并把 LDLR 极端估计作为反例保留——这一点会加分。',
  """
No estimate is removed as an outlier. Two deliberately retained extremes are
reported as findings rather than anomalies:
  1. In the multi-target replication panel, LDLR yielded an odds ratio of 2626.686
     with a direction opposite to the established biology. Its cis instrument is
     the weakest in the panel (F = 64.6, exposure beta = -0.0428), so the extreme
     value is a demonstration of what a weak cis instrument can produce, and is
     discussed as such rather than deleted.
  2. Individual Wald ratios with very wide confidence intervals are reported
     unchanged; their width is the object of study, and truncating them would
     inflate apparent precision.
No winsorisation, trimming or post-hoc exclusion is applied at any stage.
""")

f('Variables', 'Sampling weights', 'long text',
  '填 N/A 并说明原因，不要留空（必填项）。',
  """
Not applicable. The analysis uses published summary statistics, in which
individual-level sampling weights have already been applied by the source
consortium before effect estimation. No re-weighting is possible or appropriate at
the summary-statistics level. Effective sample size is used only where a
statistic formally requires it (Steiger test), and the definition used is stated
explicitly.
""")

# ---------------------------------------------------------------- 第 4 页
f('Knowledge of Data', 'Prior Publication/Dissemination', 'long text',
  '必须写“未发表”。如有预印本或会议摘要须如实列出，否则整个登记失效。',
  """
Nothing arising from this work has been published, posted as a preprint, presented
at a conference, or submitted to any journal as of the date of this registration.
Two internal technical reports dated 2026-11-10 (a start-up report with the
first-round MR result, and a robustness report covering LD clumping and
credible-set colocalisation) were circulated inside the research group only and
are not publicly accessible. They are cited in the "Prior knowledge" field so that
the sequence of work is auditable.
""")

f('Knowledge of Data', 'Prior knowledge', 'long text',
  '本页决定这次登记站不站得住。逐项列出登记前已做的分析和已知结果，一个字都不能省。',
  """
Substantial analysis had been completed before this registration, and we disclose
it in full.

Already completed before registration:
  1. A first-round two-sample MR of UKB-PPP CX3CL1 pQTL against FinnGen R13
     DM_NEUROPATHY, using the lead cis instrument. The result was null
     (beta = 0.3323, SE = 1.1073, P = 0.764, OR = 1.394).
  2. LD clumping of the credible-set variants against 1000 Genomes European
     reference, at r2 < 0.001; after clumping a single cis instrument remained.
  3. Credible-set-level colocalisation of the pQTL with the outcome in the CX3CL1
     region. The result argued against a shared causal variant (PP.H4 = 0.0584;
     CLPP = 3.56e-4; PP.H1 = 0.9415), and no variant in the outcome region reached
     P < 1e-5.
  4. Inspection of trans associations: ABO (9_133273983_A_G), FUT2, ECE1, GCNT1,
     CX3CR1 and JMJD1C loci were identified as associated with the exposure at
     genome-wide significance, and were recognised as likely horizontally
     pleiotropic.
  5. A region-wide colocalisation using the Fenland pQTL (GCST90012074) against
     FinnGen R13 DM_NEUROPATHY in chr16:57.2-57.6 Mb, on 1,744 aligned variants.
     The result again argued against a shared causal variant and, more
     informatively, attributed the absence of signal to the outcome side: PP.H4 =
     0.0319 with PP.H1 = 0.9680. Across nine prior settings PP.H4 remained
     between 0.0165 and 0.0619. In the same window 297 variants reached P < 5e-8
     for the exposure while none reached P < 1e-5 for the outcome.

Consequences that we state rather than hide:
  - This is a retrospective registration. The confirmatory re-run locks the
    endpoint set, the release set, the instrument sets, the power calculation and
    the reporting rules, but it does not create a pre-data record.
  - The choice to frame the paper as a methodological study about power and
    phenotype-definition heterogeneity was made after these exploratory results
    were known.
  - Any analysis added after this registration will be labelled exploratory and
    dated.
""")

# ---------------------------------------------------------------- 第 5 页
f('Analyses', 'Statistical models', 'long text',
  '主分析必须是单工具变量 Wald ratio，其余为敏感性分析。写清 IVW 过原点与 Egger 含截距——这是之前踩过的坑。',
  """
PRIMARY. Wald ratio for the single lead cis instrument: beta = beta_outcome /
beta_exposure, with the standard error obtained by first-order delta expansion.
This is the only estimate that carries inferential weight.

SENSITIVITY (reported in full, never substituted for the primary result):
  - Inverse-variance weighted estimate constrained through the origin (no
    intercept), for the multi-variant instrument sets.
  - MR-Egger with a free intercept; the intercept and its standard error are
    reported, and no Egger slope is interpreted when the instrument count is
    below three.
  - Weighted median.
  - Leave-one-out, applied to instrument sets only.
  - Steiger direction test, using an explicit outcome sample-size definition; the
    test is reported as not evaluable where no public case/control count exists.
  - ABF colocalisation between the pQTL and the outcome in the cis region.
  - A Monte Carlo simulation (20,000 replicates per scenario, fixed seed) used to
    validate the reporting rule set, not to generate empirical claims.

No model is fitted that the instrument count cannot support; where a method
requires three or more instruments it is reported as not applicable for
single-instrument combinations rather than run and discarded.
""")

f('Analyses', 'Effect size', 'long text',
  'OR 的量纲必须写清，否则审稿人会问“per 什么”。',
  """
The effect size is the odds ratio of the binary outcome per one-unit increase in
genetically predicted circulating CX3CL1 on the natural-log scale, as encoded by
the source pQTL betas. Because the exposure is expressed on the log scale, the
odds ratio is per log-unit and is not directly comparable to odds ratios reported
per standard deviation in studies that standardise the exposure differently; the
scale is stated wherever an estimate is reported.

Estimates are presented on the log odds ratio scale for meta-analytic comparison
and as odds ratios with 95% confidence intervals for interpretation. For
continuous secondary outcomes the same logic applies on the corresponding linear
scale, and this is stated where used.
""")

f('Analyses', 'Statistical power', 'long text',
  '核心字段。必须写“无事前效能计算”，并把事后 power 与 MDE 作为研究对象本身。这直接对应 STROBE-MR 4b。',
  """
No a priori power or sample-size calculation was performed, and none is claimed.
This is stated explicitly because the absence of an a priori calculation is itself
one of the objects of study.

Power is instead quantified post hoc, from the observed standard error of each
Wald ratio, and reported as an object of measurement rather than as a
justification. Across the 19 release-endpoint combinations, power at a true odds
ratio of 1.5 ranges from 5.45% to 10.10%, and the minimum detectable odds ratio at
80% power ranges from 5.61 to 318.57. These figures are computed with the same
formula applied uniformly to every combination, so that combinations can be
ranked against each other.

Two consequences are treated as results in their own right: (i) an interval that
excludes an odds ratio of 1.5 cannot be produced from these data, so the data
cannot distinguish a null from a moderate effect; (ii) any verbal label attached
to such an estimate is a statement about the data's resolution, not about the
biology.
""")

f('Analyses', 'Inference criteria', 'long text',
  '写清 Bonferroni 阈值、名义显著只作探索、四分类规则。这是本文的方法学产出。',
  """
The family of tests is defined as all combination-by-method estimates, giving 133
tests; the Bonferroni threshold is therefore 0.05/133 = 3.76e-4. Results below
that threshold would be called supportive; no result met it.

Uncorrected P < 0.05 is reported as nominal only, is never described as
significant in the inferential sense, and is accompanied by the number expected
under the global null (6.65 of 133) so that the observed count (8) can be judged
against chance (binomial P = 0.548).

Each combination is then assigned one of four labels by a rule fixed in advance:
  I    supportive evidence            (P below the Bonferroni threshold)
  II   adequately powered negative    (power >= 80% at OR = 1.5 and CI excluding 1)
  III  underpowered, inconclusive     (power < 80% and CI including 1)
  IV   cross-context conflict         (direction differs across releases or
                                       endpoint definitions for the same pair)
On the present data, 0 combinations are labelled I, 0 are labelled II, 15 are
labelled III and 4 are labelled IV. The label is reported for every combination
rather than only for the primary one.
""")

f('Analyses', 'Assumption Violation/Model Non-Convergence', 'long text',
  '逐条写出违背时的处置。这段能挡住“你怎么知道工具变量有效”的质疑。',
  """
Pre-specified handling of assumption violations:
  - Weak instrument: F statistics are reported for every instrument set. Where
    F >= 10 the instrument is treated as strong enough for a Wald ratio; where
    F < 10 the estimate is reported with an explicit warning and is excluded from
    any pooled interpretation. The lead cis instrument has F = 552.42.
  - Directionality: the Steiger test is run where an outcome sample size exists;
    where it does not, the test is recorded as not evaluable and no directional
    claim is made from it.
  - Horizontal pleiotropy: trans variants with known pleiotropic roles (ABO, FUT2,
    ECE1, GCNT1, CX3CR1, JMJD1C) are retained only in a sensitivity set, never in
    the primary analysis. MR-Egger intercepts are reported; where the intercept is
    unstable owing to few instruments, this is stated and the Egger estimate is
    not interpreted.
  - Heterogeneity: Cochran Q and I-squared are reported; I-squared for instrument
    strength (I2_GX) is reported so that regression dilution can be judged rather
    than assumed. I2_GX was 98.66%-99.08% across combinations, i.e. NOME is
    satisfied and a low I2_GX explanation for nominal findings is not available.
  - Colocalisation: where PP.H4 is low, a causal-sharing interpretation is
    withdrawn regardless of the MR P value.
  - Non-convergence: not applicable to Wald ratio. Where an iterative estimator
    fails to converge, the failure is reported and no estimate is substituted.
""")

f('Analyses', 'Reliability and Robustness Testing', 'long text',
  '把所有稳健性设计列全：19 组合、4 IV 集、7 靶点、仿真、阳性对照。',
  """
  - 19 release-endpoint combinations across five FinnGen releases, including
    within-domain alternative definitions and negative-control outcomes.
  - Four nested instrument sets, giving 76 estimates, treated as a sensitivity
    grid rather than as independent tests.
  - A seven-target replication panel (ACE, AGT, IL1RN, IL6R, LDLR, PCSK9 and
    SERPINE1) analysed with the identical code path, to test
    whether the power pattern is a property of these endpoints or of instrument
    strength interacting with endpoint size.
  - A technical positive control: 56 variant-by-endpoint tests, of which five
    reached P < 5e-8, all driven by the ABO trans variant, while all 14 cis-only
    tests were non-significant. This establishes that the pipeline can detect
    real signal and that the CX3CL1 null is not a pipeline artefact.
  - A region-wide approximate-Bayes-factor colocalisation in the full 400 kb
    window chr16:57.2-57.6 Mb, using 1,744 variants shared between the Fenland
    pQTL and FinnGen after allele alignment, reported alongside the
    credible-set-level analysis. Allele orientation is resolved by exact match,
    strand complementation or sign flip; variants whose configuration cannot be
    resolved are discarded rather than imputed. Posteriors are reported for
    prior widths W = 0.01, 0.04 and 0.16 and for shared-variant priors p12 =
    5e-6, 1e-5 and 2e-5.
  - A Monte Carlo simulation with 20,000 replicates per scenario and a fixed seed,
    used to quantify how often a naive reading of an underpowered null reaches the
    wrong conclusion.
  - A one-command offline reproduction script, version-locked dependencies, and
    SHA256 checksums for every input file, so that every reported number can be
    regenerated from the deposited inputs.
""")

f('Analyses', 'Exploratory analysis', 'long text',
  '明确把登记前完成的四类分析标为探索性，并声明新增分析的处理方式。',
  """
The following are exploratory and are labelled as such wherever reported:
  1. The first-round MR, LD clumping and credible-set colocalisation described
     under "Prior knowledge" - all completed before this registration.
  2. The seven-target replication panel, which was conceived as a follow-up to
     the single-target result rather than specified in advance.
  3. The region-wide colocalisation against the Fenland pQTL, which was added
     after the credible-set-level analysis had already been interpreted, and
     after it had become clear that the UKB-PPP full summary statistics were
     not retrievable from this setting.
  4. The technical positive control using the ABO trans variant.
  5. Any association observed in the trans sensitivity set.
  6. Any endpoint or release added after this registration.

Rule for additions: anything analysed after this registration is reported as
exploratory, dated, and listed in a registration update with a justification,
regardless of its result. No exploratory finding will be promoted to the primary
conclusion of the paper.
""")


def main():
    rows = []
    for i, d in enumerate(FIELDS, 1):
        rows.append({
            'order': i,
            'page': d['page'],
            'field': d['field'],
            'input_type': d['input_type'],
            'status': d['status'],
            'cn_hint': d['cn_hint'].replace('\n', ' '),
            'en_value': d['en'],
            'n_words': len(d['en'].split()),
        })
    cols = ['order', 'page', 'field', 'input_type', 'status', 'cn_hint',
            'en_value', 'n_words']
    with io.open(os.path.join(RES, 'a1_osf_field_mapping.tsv'), 'w',
                 encoding='utf-8-sig', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    lines = [
        'A1 - OSF RETROSPECTIVE REGISTRATION: PASTE BLOCKS',
        'Template: Secondary Data Preregistration (OSF), 27 fields, 5 pages',
        'Prepared %s. Copy each block into the matching OSF field.' % TODAY,
        'Do NOT paste this whole file into a single field.',
        '',
    ]
    for r in rows:
        lines += [
            '=' * 72,
            '[%02d] PAGE: %s' % (r['order'], r['page']),
            'FIELD: %s   (%s)' % (r['field'], r['input_type']),
            '-' * 72,
            r['en_value'].strip(),
            '',
        ]
    with io.open(os.path.join(RES, 'a1_osf_paste_blocks.txt'), 'w',
                 encoding='utf-8') as fh:
        fh.write('\n'.join(lines))

    print('fields', len(rows))
    print('total words', sum(r['n_words'] for r in rows))
    print('by page:')
    for p in ['Study Information', 'Data Description', 'Variables',
              'Knowledge of Data', 'Analyses']:
        n = [r for r in rows if r['page'] == p]
        print('  %-20s %2d fields  %5d words' % (p, len(n),
                                                 sum(x['n_words'] for x in n)))
    print('verify markers', sum(r['en_value'].count('[VERIFY') for r in rows))


if __name__ == '__main__':
    main()
