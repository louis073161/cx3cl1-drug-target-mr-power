# -*- coding: utf-8 -*-
"""
A 级 18 项清零 · 工程化与对照表更新（STROBE-MR 4a/4e/9a/9b/13d/19）

产出：
  code/requirements.lock.txt            真实运行环境的版本锁定
  code/run_all.sh                       一键复现入口（默认离线重算，--with-fetch 含联网取数）
  results/a1_input_sha256.tsv           输入文件校验和清单
  results/a1_version_snapshot.json      数据源版本、参数、随机种子、环境
  results/a1_osf_registration_pack.txt  OSF 回顾性登记可粘贴文本包
  results/a1_a18_prior_mr_comparison.tsv  13d 已发表 MR 研究口径对照表
  results/a1_strobe_a18_cleared.tsv     A 级 18 项逐条清零状态
  results/a1_a18_clearance_summary.json 清零统计与待外部动作
"""
import csv
import hashlib
import importlib.metadata as md
import io
import json
import os
import platform
import sys
from datetime import datetime

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
CODE = ROOT + '/code'
RES = ROOT + '/results'
DATA = ROOT + '/data'
TODAY = datetime.now().strftime('%Y-%m-%d')


def write_text(path, text):
    with io.open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(text)


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0].keys())
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def pkg_version(name):
    try:
        return md.version(name)
    except Exception:
        return 'NOT_INSTALLED_OR_UNUSED'


# ---------------------------------------------------------------- 9a 环境锁定
def requirements():
    lines = [
        '# A1 reproducibility lock — generated %s' % TODAY,
        '# Python %s (%s)' % (platform.python_version(), platform.system()),
        '# 分析链仅依赖下列包；python-docx 仅用于中文报告排版，不参与统计。',
    ]
    for p in ['numpy', 'scipy', 'matplotlib', 'pandas', 'pyarrow',
              'requests', 'python-docx']:
        lines.append('%s==%s' % (p, pkg_version(p)))
    lines += [
        '',
        '# 关键设置（与代码一致，见 results/a1_version_snapshot.json）',
        '# LD clumping: 1000 Genomes EUR via Ensembl REST, r2 < 0.001, unknown pairs filled 0',
        '# coloc: ABF prior variance W=0.04, p1=1e-4, p2=1e-4, p12=1e-5 (sensitivity 1e-6 / 1e-4)',
        '# weighted median: bootstrap 2000; A18 IVW bootstrap: 5000',
        '# random seeds: 20261110 (weighted median / simulation / A18 bootstrap)',
        '# two-sided alpha=0.05; CI=1.959964; power target=80%',
    ]
    write_text(CODE + '/requirements.lock.txt', '\n'.join(lines) + '\n')


def run_all():
    sh = '''#!/usr/bin/env bash
# A1 一键复现入口。默认只跑离线重算（不联网）；--with-fetch 才执行公开桶取数。
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONPATH="$ROOT/code${PYTHONPATH:+:$PYTHONPATH}"

if [[ "${1:-}" == "--with-fetch" ]]; then
  echo "[1/14] 公开桶取数（需联网；FinnGen GCS + tabix HTTP Range）"
  "$PY" code/a1_multi_fetch.py
  echo "[2/14] LD / clumping / coloc / 技术阳性对照 / 多靶点队列（需联网）"
  "$PY" code/a1_ld_coloc.py
  "$PY" code/a1_pip_coloc.py
  "$PY" code/a1_positive_control.py
  "$PY" code/a1_b1b2_multitarget.py
  echo "[3/14] B3 附录：Fenland pQTL 区域取数（需联网；EBI FTP，约 500 MB）"
  "$PY" code/a1_b3_parallel_dl.py
  "$PY" code/a1_b3_fenland_fetch.py
  echo "[4/14] B3 附录：区域级完整共定位（需联网；FinnGen GCS）"
  "$PY" code/a1_b3_coloc_full.py
else
  echo "[1/14] 跳过联网取数，使用已落盘的真实汇总统计（加 --with-fetch 可重取）"
  echo "[2/14] 跳过联网 LD / coloc / 多靶点取数；沿用已落盘结果与 data/b1_target_pqtl.json"
  echo "[3/14] 跳过 B3 的 EBI 下载；沿用 results/a1_b3_fenland_region.tsv"
  echo "[4/14] 跳过 B3 的 FinnGen 取数；沿用 results/a1_b3_coloc_full.json"
fi

echo "[5/14] 19 组合 MR 矩阵、效能曲线、结论漂移"
"$PY" code/a1_multi_mr.py
echo "[6/14] 图 2—4（内部留档中文版）"
"$PY" code/a1_multi_plots.py
echo "[7/14] 图 1 区段关联对照图"
"$PY" code/a1_region_plot.py
echo "[8/14] STROBE-MR 20 项逐条对照表"
"$PY" code/a1_strobe_checklist.py
echo "[9/14] 决策流程仿真验证（图 5）"
"$PY" code/a1_simulation.py
echo "[10/14] A 级 12b / 13c / 11a 实算"
"$PY" code/a1_a18_stats.py
echo "[11/14] A4 漂移归因 + A5 四类工具变量集 76 估计"
"$PY" code/a1_a4a5_drift.py
echo "[12/14] 图 1—5 英文矢量投稿版"
"$PY" code/a1_figure_pack_en.py
echo "[13/14] 图 9 B3 区域级共定位英文矢量版"
"$PY" code/a1_b3_figure.py
echo "[14/14] 英文稿定稿骨架、投稿包、可复核包、SHA256、版本快照与 A 级 18 项清零对照表"
"$PY" code/a1_manuscript_pack.py
"$PY" code/a1_a18_clear.py
echo "完成。全部结果写入 results/；投稿包写入 submission_pack/。"
'''
    write_text(CODE + '/run_all.sh', sh)
    try:
        os.chmod(CODE + '/run_all.sh', 0o755)
    except Exception:
        pass


# ---------------------------------------------------------------- 19 校验和与快照
SELF_GENERATED_LAST = ('a1_input_sha256.tsv', 'a1_version_snapshot.json')


def checksums():
    rows = []
    for folder, label in [(CODE, 'code'), (DATA, 'data'), (RES, 'results')]:
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            p = os.path.join(folder, name)
            if not os.path.isfile(p) or name.endswith('.INCOMPLETE'):
                continue
            if label == 'results' and name in SELF_GENERATED_LAST:
                continue   # 这两个文件在校验清单生成时才写出，避免自引用
            h = hashlib.sha256()
            with open(p, 'rb') as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b''):
                    h.update(chunk)
            rows.append({'folder': label, 'file': name,
                         'size_bytes': os.path.getsize(p),
                         'sha256': h.hexdigest(),
                         'local_mtime': datetime.fromtimestamp(
                             os.path.getmtime(p)).isoformat(timespec='seconds'),
                         'note': '本地文件修改时间，非数据源官方发布日期'})
    write_tsv(RES + '/a1_input_sha256.tsv', rows)
    return len(rows)


def version_snapshot(n_files):
    snap = {
        'generated': TODAY,
        'analysis': 'A1 drug-target MR power / phenotype-heterogeneity methodological paper',
        'data_sources': {
            'exposure': 'Open Targets Platform v26.9; UKB-PPP UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1; 9 SuSiE 95% credible sets; nSamples=33,436',
            'outcome': 'FinnGen public summary statistics R9, R10, R11, R12, R13 (DF13); 19 (release, endpoint) combinations',
            'ld_reference': '1000 Genomes phase 3 EUR via Ensembl REST (GRCh38)',
            'coloc': 'Open Targets precomputed coloc (COLOC_PIP_ECAVIAR); ABF coloc on FinnGen R13 chr16:57.2-57.6 Mb',
            'sample_overlap': 'None (UKB-PPP United Kingdom vs FinnGen Finland)',
        },
        'versions_explicitly_not_assumed': [
            'FinnGen R11 / R13 have no public manifest; case counts not imputed.',
            'Release public dates are NOT recorded here; cite primary FinnGen / UKB-PPP publications instead of guessing.',
            'OT 26.9 and FinnGen R9-R13 identifiers come from the fetched files, not from a guessed version list.',
        ],
        'parameters': {
            'clumping_r2_threshold': 0.001,
            'clumping_unknown_pair_fill': 0.0,
            'coloc_prior_variance_W': 0.04,
            'coloc_p1': 1e-4, 'coloc_p2': 1e-4, 'coloc_p12': 1e-5,
            'coloc_p12_sensitivity': [1e-6, 1e-4],
            'weighted_median_bootstrap': 2000,
            'a18_ivw_bootstrap': 5000,
            'simulation_replicates_per_point': 20000,
            'random_seed': 20261110,
            'alpha_two_sided': 0.05, 'ci_z': 1.959964, 'power_target': 0.80,
        },
        'environment': {
            'python': platform.python_version(),
            'os': platform.platform(),
            **{p: pkg_version(p) for p in
               ['numpy', 'scipy', 'matplotlib', 'pandas', 'pyarrow']},
        },
        'n_files_in_checksum_manifest': n_files,
        'api_log': 'results/a1_fetch_log.json: 133 region requests, 0 errors, 113 s; '
                   'results/a1_positive_control.json: 49 requests, 0 errors',
        'registration': {
            'type': 'retrospective (NOT prospective preregistration)',
            'osf_status': 'public and archived 11 November 2026 as a Secondary Data Preregistration (osf.io/apwe2) under project osf.io/pmhve; DOI 10.17605/OSF.IO/APWE2',
            'disclosure': 'exploratory MR, LD clumping and coloc were completed before registration and are disclosed',
        },
        'hosting': {
            'github': 'pending author action',
            'zenodo_doi': 'pending author action',
            'osf_doi': '10.17605/OSF.IO/APWE2',
        },
    }
    with io.open(RES + '/a1_version_snapshot.json', 'w', encoding='utf-8') as f:
        json.dump(snap, f, ensure_ascii=False, indent=2)


# ---------------------------------------------------------------- 9b OSF 登记文本包
def osf_pack():
    txt = '''A1 — OSF RETROSPECTIVE REGISTRATION (ready-to-paste text)
Prepared 2026-11-10. Registration itself, date and DOI require the author's OSF account.

TITLE
Power and phenotype-definition heterogeneity in drug-target Mendelian randomization:
a worked example using circulating CX3CL1 and diabetic neuropathy across five FinnGen releases

TYPE
Retrospective registration of an analysis plan for a methodological study using
publicly available summary-level data. This is NOT a prospective preregistration.

RETROSPECTIVE DISCLOSURE (must be pasted verbatim)
Before this registration, the authors had already completed exploratory analyses:
a first two-sample MR of UKB-PPP CX3CL1 pQTL against FinnGen R13 DM_NEUROPATHY,
LD clumping with 1000 Genomes EUR references, and credible-set-level colocalisation.
Those results are documented in two internal reports dated 2026-11-10. We do not
claim that the data were untouched before registration. This registration locks the
confirmatory re-run (nine instruments x multiple endpoints x multiple releases),
the power analysis, the reporting framework and the wording rules. Any deviation
will be logged with an OSF timestamp.

RESEARCH QUESTIONS
1. In a drug-target MR with a single independent cis instrument and no strong outcome
   signal in the region, what strength of conclusion can public data support?
2. How much do point estimates and conclusions drift across endpoint definitions,
   FinnGen releases and estimator choices?
3. Which reporting rules prevent an underpowered null from being written as "no effect"?

DESIGN AND DATA
Summary-level two-sample Mendelian randomisation. Exposure: Open Targets v26.9 / UKB-PPP
CX3CL1 pQTL, nine SuSiE 95% credible sets, n=33,436 European ancestry. Outcomes:
FinnGen public summary statistics R9-R13, endpoints DM_NEUROPATHY, DM_POLYNEURO,
G6_POLYOTHUNS and diabetic microvascular complications (including STRICT / EXMORE
variants) as negative controls, 19 combinations in total. UKB-PPP and FinnGen are
non-overlapping populations. Download URLs, tabix indices and manifests are recorded
in a1_endpoint_dict.tsv and a1_fetch_log.json.

PRIMARY OUTCOME MEASURES
Wald ratio beta, SE, P, OR and 95% CI for each combination; post-hoc power at OR=1.5
and minimum detectable OR at 80% power, both derived from the observed SE; range of
point estimates and number of sign flips; observed versus expected number of
nominally significant results across 133 combination x method tests.

ANALYSIS PLAN
Primary analysis is the Wald ratio of the single cis lead instrument rs62037084
(16_57379748_C_T). Sensitivity analyses: IVW, MR-Egger, weighted median and
leave-one-out after LD clumping at r2<0.001; Steiger direction test; I2, I2_GX,
Cochran Q, E-value and instrument-level bootstrap; ABF colocalisation; a simulation
validating the reporting decision framework; a pipeline / power calibration using
the same code on large FinnGen endpoints. No missing data are imputed; unmatched
variants are reported as excluded. Power is driven by observed SE, not by assumed
case counts.

HANDLING OF NULL RESULTS
All results are reported regardless of direction. An underpowered null must not be
written as "no effect"; nominally significant exploratory findings must not be
written as evidence for causality. Multiplicity and pleiotropy checks precede any
exploratory interpretation.

STATEMENT ON POWER (STROBE-MR item 4b)
No a priori power or sample-size calculation was performed. Post-hoc power and the
minimum detectable effect are the object of this study and will be reported as such.

DATA / CODE
All inputs, code, checksums and figures are local and will be deposited with DOI.
FinnGen, Open Targets and UKB-PPP terms and primary publications will be cited.

AUTHOR
Corresponding author: Dr Yi Lu (鲁义主任医师), Department of Anesthesiology,
Guangzhou Medical University Affiliated Traditional Chinese Medicine Hospital.
[Affiliation and ORCID to be completed at submission.]
'''
    write_text(RES + '/a1_osf_registration_pack.txt', txt)


# ---------------------------------------------------------------- 13d 文献对照表
def prior_mr_table():
    rows = [
        {'citation': 'Zhong J et al. J Pain Res 2026;19:1-9. doi:10.2147/JPR.S570828',
         'verified_content': '循环炎症蛋白 × 四类神经病理性疼痛双向两样本 MR；CX3CL1 与舌咽神经病风险关联显著（OR 0.182, 95%CI 0.048-0.691）',
         'outcome_definition': '四类神经病理性疼痛（非 DPN 特异）',
         'n_cases_controls': '未核实（本轮未取回全文）',
         'power_reported': '待全文核对',
         'coloc_or_fine_mapping': '无（已核实：未覆盖）',
         'difference_from_this_study': '本文为 DPN 特异结局、报告效能与端点定义漂移；不追求阳性，且做可信集级 coloc'},
        {'citation': 'Hu C et al. Front Immunol 2026. doi:10.3389/fimmu.2026.1681223',
         'verified_content': '731 免疫性状与 4,719/2,923 血浆蛋白对 DPN 的 MR + 中介 MR + 人腓肠神经单细胞与流式验证；全扫描设计',
         'outcome_definition': 'DPN（FinnGen，2,843 例 / 389,580 对照，以原文为准）',
         'n_cases_controls': '2843 / 389580（二手转述，须核对原文）',
         'power_reported': '待全文核对',
         'coloc_or_fine_mapping': '待全文核对',
         'difference_from_this_study': '本文为假说驱动的定向检验 + 多 release 稳健性 + 报告规范，不做全扫描'},
        {'citation': 'Huang K et al. Sci Rep 2025;15:23088. doi:10.1038/s41598-025-07178-0',
         'verified_content': '脑瘫的 SMR + MR + 共定位 + 可药性评估 + PheWAS，识别 CX3CL1、TYRO3 为候选靶点',
         'outcome_definition': '脑瘫（非 DPN）',
         'n_cases_controls': '未核实（本轮未取回全文）',
         'power_reported': '待全文核对',
         'coloc_or_fine_mapping': '有（已核实：含共定位）',
         'difference_from_this_study': '疾病领域不同；本文不移植其结论，只比较报告口径'},
        {'citation': 'Ji L et al. Diabetol Metab Syndr 2025;17. doi:10.1186/s13098-025-01696-7; Yu Y et al. Neurol Res 2025. PMID 40888768',
         'verified_content': '731 免疫细胞表型与糖尿病神经病变的双向 MR',
         'outcome_definition': '糖尿病神经病变',
         'n_cases_controls': '未核实（本轮未取回全文）',
         'power_reported': '待全文核对',
         'coloc_or_fine_mapping': '无（已核实：无分子 QTL 与共定位）',
         'difference_from_this_study': '本文用分子 pQTL 药物靶点框架并报效能下限'},
        {'citation': 'Wang Z et al. Front Endocrinol 2025',
         'verified_content': '179 种血浆脂质与糖尿病神经病变 MR，MVMR 校正 BMI 与 HbA1c',
         'outcome_definition': '糖尿病神经病变',
         'n_cases_controls': '未核实（本轮未取回全文）',
         'power_reported': '待全文核对',
         'coloc_or_fine_mapping': '待全文核对',
         'difference_from_this_study': '仅脂质组暴露；本文为蛋白药物靶点 + 报告规范'},
    ]
    write_tsv(RES + '/a1_a18_prior_mr_comparison.tsv', rows)


# ---------------------------------------------------------------- A 级 18 项状态
CLEAR = [
    ('1', '已清', '英文 Title / Abstract 成稿；2026-11-11 已并入定稿骨架并完成结构化压缩（四段式）',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; _markdown源/A1英文稿定稿骨架_20261111.md',
     '无', '题名与摘要首句均写明 two-sample summary-level MR；期刊字数确认后仅作编辑器删减'),
    ('2', '已清', 'Introduction §I1–I2 成稿：靶点背景 + 药物靶点 MR 价值 + 报告规范缺口；2026-11-11 补 I4 差异化段',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; _markdown源/A1英文稿定稿骨架_20261111.md; a1_citation_anchors.tsv',
     '无', '41% 已核实为推算值；原文直接报告 59% 满足效能条目，11% 分母为 2021 年后文章；正文不写 41%'),
    ('3', '已清', 'Introduction §I3 成稿：三个编号目标 + 条件性因果解释 + 回顾性无预设方向',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md', '无', '禁止写成验证性确证试验'),
    ('4a', '已清（结构）', 'Table 1 数据源表成稿：人群、祖源、release、端点定义、URL 全部来自已落盘文件',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; a1_endpoint_definition_table.tsv',
     '各 release 官方公开日期须引用原始文献填入，本轮未推测',
     'R11/R13 无 manifest 已显式标注 NO_PUBLIC_MANIFEST'),
    ('4b', '已清', 'Methods 明确声明未做事前效能/样本量计算，并说明事后效能与 MDE 为研究对象',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; a1_power_curve.tsv', '无',
     'OR=1.5 把握度 5.45%–10.10%；MDE OR 5.61–318.57'),
    ('4e', '已清（声明）/ 待用户动作', '伦理、知情同意、数据条款声明成稿；FinnGen 官网填表流程写入待办',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; a1_version_snapshot.json',
     '本人完成 FinnGen 官网填表并留存回执；补 Kurki et al. Nature 2023 致谢',
     '填表为账号/机构动作，不代做、不冒充完成'),
    ('5', '已清', 'Methods Assumptions 段逐条写三大 IV 假设 + InSIDE/NOME + 基因-环境等价',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md', '无',
     'I²_GX 实测 98.66%–99.08%，NOME 满足，已按实测写'),
    ('8', '已清', 'Methods Simulations 段引用仿真（20,000 次/点）与阳性对照（56 组合）作为 instrument validation',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md; a1_simulation_results.json; a1_positive_control.json',
     '无', '仿真明确标注不产生新实证结论'),
    ('9a', '已清', 'requirements.lock.txt + 版本快照记录全部设置；2026-11-11 在 Declarations 写入完整软件与版本说明',
     'requirements.lock.txt; a1_version_snapshot.json; submission_pack/a1_declarations.txt', '无',
     'Python %s；numpy/scipy/matplotlib 版本已锁定；代码清单已逐文件列明' % platform.python_version()),
    ('9b', '已清（文本）/ 待用户动作', 'OSF 回顾性登记可粘贴文本包成稿，含登记前探索性分析逐项披露',
     'a1_osf_registration_pack.txt; _markdown源/A1_A级18项清零与英文稿骨架_20261110.md',
     '本人注册 OSF 并把日期与 DOI 回填 Methods / Declarations / 版本快照',
     '不得写成前瞻性预注册'),
    ('11a', '已清', 'lead 位点三组关联表 19 组合全落盘；R13 主结局已写入定稿骨架 R1 段，其余入 Supplementary S2',
     'a1_a18_lead_variant_associations.tsv; _markdown源/A1英文稿定稿骨架_20261111.md; submission_pack/a1_supplementary_index.tsv',
     '无', 'R13 暴露 P=3.74×10⁻¹²²，结局 P=0.764，F=552.42'),
    ('12b', '已清（含口径更正）', '实算 I²、I²_GX、Cochran Q、E-value、Steiger、ABF、Bootstrap CI、F',
     'a1_a18_stats.tsv/json; a1_a18_steiger.tsv; _markdown源/A1_A级18项清零与英文稿骨架_20261110.md',
     '无', '原预设"I²_GX 低值解释名义显著"被实测推翻：I²_GX 98.66%–99.08%（高，NOME 满足）；改用多重检验期望 6.65 vs 实测 8（二项 P=0.548）、trans 水平多效性与 Egger 截距显著解释'),
    ('13c', '已清', 'Steiger 方向检验三种 N_out 口径实算；反向 MR 不适用理由成稿',
     'a1_a18_steiger.tsv; _markdown源/A1_A级18项清零与英文稿骨架_20261110.md',
     '无',
     'R13 主结局保守边界（N_out=8,290，来自 r13 专属端点页 4,145 例）P=2.49×10⁻²⁴，方向支持暴露→结局；'
     '8 个无公开病例数的组合如实标 NOT_EVALUABLE，不补数；近似实现已注明'),
    ('13d', '已清（表）/ 待全文核对', '已发表 MR 研究口径对照表成稿（数据源、端点、效能报告、结论措辞、差异）',
     'a1_a18_prior_mr_comparison.tsv; _markdown源/A1_A级18项清零与英文稿骨架_20261110.md',
     '逐篇取回全文核对"是否报效能"后把"待全文核对"替换为确切结论',
     '未核实字段一律留"待全文核对"，不编造'),
    ('15', '已清', '偏倚方向 + 量级表 6 条 + Discussion Limitations 段成稿',
     'a1_a18_bias_magnitude.tsv; _markdown源/A1_A级18项清零与英文稿骨架_20261110.md', '无',
     '明确写 F=552.4>>10，弱工具偏倚不是主因；主因为结局侧 SE 与暴露效应尺度'),
    ('16a', '已清', 'Discussion 总体解释段成稿，严格限定"未获得支持性证据、不足以检验中等效应"',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md', '无', '与撞车文献阳性口径切开'),
    ('16b', '已清', 'Discussion 单列基因-环境等价段：终身遗传暴露 ≠ 短期药物干预；循环 pQTL ≠ 脊髓局部信号',
     '_markdown源/A1_A级18项清零与英文稿骨架_20261110.md', '无', '不推翻标书假说，不改原三部分湿实验'),
    ('19', '已清（本地）/ 待用户动作',
     'run_all.sh（已扩展为 14 步，含 B3 附录）、依赖锁定、SHA256、版本快照、Code/Data availability 声明、GitHub README 与 Zenodo 元数据草稿均已成稿',
     'run_all.sh; requirements.lock.txt; a1_input_sha256.tsv; a1_version_snapshot.json; '
     'submission_pack/a1_declarations.txt; submission_pack/a1_github_readme.md; submission_pack/a1_zenodo_metadata.json',
     '本人建 GitHub 仓库 + Zenodo DOI + OSF 三者互链，并把 DOI 回填正文',
     '本地可复核已闭环（离线 14 步实跑零失败）；公开托管与 DOI 为账号动作，不代做'),
]


RESIDUAL_TYPE = {
    '4a': '投稿前引用核对（官方 release 日不可核实）',
    '4e': '需账号/机构动作',
    '9b': '需账号/机构动作（本轮用户已排除）',
    '13d': '投稿前全文核对',
    '19': '需账号/机构动作',
}


def clearance_table():
    src = list(csv.DictReader(io.open(RES + '/a1_strobe_mr_checklist.tsv',
                                     encoding='utf-8-sig'), delimiter='\t'))
    A = [r for r in src if r['priority'].strip() == 'A']
    by_id = {r['item'] + r['subitem']: r for r in A}
    rows = []
    for key, status, result, evidence, residual, note in CLEAR:
        base = by_id[key]
        rows.append({
            'item': base['item'], 'subitem': base['subitem'],
            'section': base['section'], 'verdict_before': base['verdict'],
            'official_item_zh': base['official_item_zh'],
            'action_original': base['action'],
            'cleared_status': status,
            'what_was_done': result,
            'evidence_files': evidence,
            'residual_action': residual,
            'residual_type': RESIDUAL_TYPE.get(key, '无'),
            'note': note,
            'assessed_on': TODAY,
        })
    write_tsv(RES + '/a1_strobe_a18_cleared.tsv', rows)

    counts, rtypes = {}, {}
    for r in rows:
        counts[r['cleared_status']] = counts.get(r['cleared_status'], 0) + 1
        rtypes[r['residual_type']] = rtypes.get(r['residual_type'], 0) + 1
    summary = {
        'generated': TODAY,
        'n_A_items': len(rows),
        'status_counts': counts,
        'residual_type_counts': rtypes,
        'fully_cleared_locally': sum(1 for r in rows
                                     if r['residual_type'] == '无'),
        'items_requiring_author_action': [
            {'item': r['item'], 'status': r['cleared_status'],
             'action': r['residual_action']}
            for r in rows if r['residual_type'].startswith('需账号')],
        'items_for_pre_submission_check': [
            {'item': r['item'], 'action': r['residual_action']}
            for r in rows if r['residual_type'].startswith('投稿前')],
        'corrections_to_original_checklist': [
            '12b: I2_GX measured 98.66%-99.08% (high, NOME satisfied); the original '
            'suggestion to use a low I2_GX to explain 7 nominally significant MR-Egger '
            'results is not supported by the data and has been rewritten.',
            '4b: F statistic is 552.42 (not 22.25); 22.25 is the 80%-power MDE OR.',
            '13c: Steiger implemented as a transparent approximation, not the official '
            'steiger package; 8 combinations without public case counts are NOT imputed.',
        ],
        'disclaimer': 'OSF registration, GitHub/Zenodo DOI, FinnGen web form and '
                      'full-text verification of prior MR papers require the author '
                      'account or full texts; materials are prepared, actions are not '
                      'claimed as completed.',
    }
    with io.open(RES + '/a1_a18_clearance_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    return summary


def main():
    requirements()
    run_all()
    osf_pack()
    prior_mr_table()
    s = clearance_table()
    n = checksums()          # 最后统一算校验和，纳入本轮全部新产出
    version_snapshot(n)
    print('checksum files', n)
    print('A items', s['n_A_items'], 'status', s['status_counts'])
    print('fully cleared locally', s['fully_cleared_locally'])
    print('author actions', len(s['items_requiring_author_action']))


if __name__ == '__main__':
    main()
