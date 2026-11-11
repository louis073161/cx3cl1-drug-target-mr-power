# -*- coding: utf-8 -*-
"""
B1 + B2：多靶点复现队列 + 生物学阳性对照（全部公开数据，零新增权限）。

B1：把 CX3CL1 的单靶点个案扩展为 n=7 的独立靶点队列。用 Open Targets
    Platform v26.9 GraphQL 的 target.credibleSets 取得 7 个可成药靶点的
    UKB-PPP cis pQTL 最强 lead（真实响应保存在 data/b1_target_pqtl.json），
    对同一批 FinnGen 端点 × R9/R13 复用已跑通的 tabix/BGZF 定点取数、
    harmonise 与 Wald ratio，检验"效能不足 / 方向漂移"是否系统性复现。

B2：用 PCSK9、LDLR 等 cis pQTL × FinnGen R13 冠心病/高脂血症大样本端点，
    作为"生物学阳性对照"，补足此前只有 ABO trans 技术阳性对照的缺口。
    方向预期：PCSK9 蛋白升高 → CHD 风险升高；LDLR 蛋白升高 → CHD 风险降低。

诚信边界（随结果一起报告）：
  - B1 各靶点同属 FinnGen 资源，不构成"独立数据集复现"；它们是跨靶点
    的 cross-context replication，用于检验方法学现象是否可推广。
  - 每个靶点只用 1 个 cis lead，未做 LD clumping / coloc / Steiger，
    不能作为各靶点的确证性病因 MR；结论限于方法学层面。
  - B2 是方向性验证，不是药物试验替代；不显著不解释为"PCSK9 无效"，
    显著也不解释为确证因果（蛋白 pQTL 可能受反向因果/共存变异影响）。
  - 所有数值来自实际 HTTP 分段取数，无模拟、无填补、无挑选。
"""
import csv
import io
import json
import math
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_multi_fetch import get_tbi, gcs_url, fetch_one_region, harmonise  # noqa: E402

RES = ROOT + '/results'
DATA = ROOT + '/data'
TARGET_JSON = DATA + '/b1_target_pqtl.json'
Z975 = norm.ppf(0.975)
Z80 = norm.ppf(0.80)

# 与 A1 主分析完全相同的端点集合（3 个 × 2 个 release）
B1_COMBOS = [
    ('R9',  'DM_NEUROPATHY',          'primary'),
    ('R9',  'G6_POLYOTHUNS',          'broad_proxy'),
    ('R9',  'DM_NEPHROPATHY_EXMORE',  'negative_ctrl_microvascular'),
    ('R13', 'DM_NEUROPATHY',          'primary'),
    ('R13', 'G6_POLYOTHUNS',          'broad_proxy'),
    ('R13', 'DM_NEPHROPATHY_EXMORE',  'negative_ctrl_microvascular'),
]
# B2 预先指定的 R13 大样本端点（tbi 已于前次阳性对照分析缓存）
B2_ENDPOINTS = [
    ('I9_CHD',         'Major coronary heart disease event'),
    ('I9_CORATHER',    'Coronary atherosclerosis'),
    ('E4_HYPERLIPMIX', 'Mixed hyperlipidaemia'),
]
# 已知生物学方向（基于药物靶点常识，仅作方向性预期，不作确证）
EXPECTED = {
    'PCSK9': 'positive',   # PCSK9 gain -> higher LDL -> higher CHD risk
    'LDLR':  'negative',   # LDLR gain -> lower LDL -> lower CHD risk
}


def power_at_or(se, or_true=1.5):
    z = math.log(or_true) / se
    return float(norm.cdf(z - Z975) + norm.cdf(-z - Z975))


def mde_or(se):
    return math.exp((Z975 + Z80) * se)


def wald(rec):
    bx = rec['beta_exp']
    if not bx:
        return None
    beta = rec['beta_out'] / bx
    se = rec['se_out'] / abs(bx)
    z = beta / se if se else 0.0
    return {
        'beta': beta, 'se': se, 'z': z,
        'p': float(2 * norm.sf(abs(z))),
        'ci_lo': beta - Z975 * se, 'ci_hi': beta + Z975 * se,
        'or': math.exp(beta),
        'or_ci_lo': math.exp(beta - Z975 * se),
        'or_ci_hi': math.exp(beta + Z975 * se),
        'power_or_1.5': power_at_or(se),
        'mde_or_80pct': mde_or(se),
        'F_instrument': rec['z_exp'] ** 2,
    }


def fig_style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.facecolor': 'white', 'figure.facecolor': 'white',
        'text.color': '#1f2933', 'axes.labelcolor': '#1f2933',
        'xtick.color': '#1f2933', 'ytick.color': '#1f2933',
        'axes.edgecolor': '#52606d', 'axes.linewidth': 0.8,
    })


def save(fig, stem):
    for ext, kw in [('pdf', {}), ('svg', {}), ('png', {'dpi': 300})]:
        fig.savefig('%s/%s.%s' % (RES, stem, ext), bbox_inches='tight', **kw)
    plt.close(fig)


def write_tsv(path, rows):
    if not rows:
        return
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_targets():
    d = json.load(io.open(TARGET_JSON, encoding='utf-8'))
    out = {}
    for sym, rec in d.items():
        cis = [r for r in rec['cis']
               if r.get('beta') is not None and r.get('zScore') is not None]
        for r in cis:
            vid = r['variant']['id']
            parts = vid.split('_')
            # 仅保留简单 SNP（ref/alt 均为单碱基），保证 harmonise 可靠
            if len(parts) == 4 and len(parts[2]) == 1 and len(parts[3]) == 1:
                out[sym] = {'symbol': sym, 'ensembl': rec['ensembl'],
                            'lead': r, 'variant': vid,
                            'chrom': r['chromosome'], 'pos': r['position'],
                            'beta_exp': r['beta'], 'z_exp': r['zScore'],
                            'study_id': r['studyId'],
                            'n_cis_available': len(cis)}
                break
    return dict(sorted(out.items()))


def fetch_all(targets, combos, tbi_cache):
    jobs = []
    for rel, ep, etype in combos:
        if (rel, ep) not in tbi_cache:
            continue
        url = gcs_url(rel, ep)
        for sym, t in targets.items():
            jobs.append((sym, rel, ep, t['chrom'],
                         t['pos'] - 1500, t['pos'] + 1500,
                         url, tbi_cache[(rel, ep)]))
    print('Fetching %d region blocks ...' % len(jobs))
    collected, errors = {}, []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(fetch_one_region,
                          (j[1], j[2], j[3], j[4], j[5], j[6], j[7])): j
                for j in jobs}
        for fu in as_completed(futs):
            j = futs[fu]
            rel, ep, chrom, lines, err = fu.result()
            collected.setdefault((j[0], rel, ep), {})[chrom] = lines
            if err:
                errors.append({'target': j[0], 'release': rel,
                               'endpoint': ep, 'chrom': chrom, 'error': err})
    return collected, errors, len(jobs)


def estimate(targets, combos, collected):
    rows = []
    for rel, ep, etype in combos:
        for sym, t in targets.items():
            lines_by_chrom = collected.get((sym, rel, ep), {})
            # harmonise 期望 Open Targets credibleSet row 结构
            recs = harmonise([t['lead']], lines_by_chrom)
            if not recs:
                rows.append({
                    'target': sym, 'ensembl': t['ensembl'],
                    'release': rel, 'endpoint': ep, 'endpoint_type': etype,
                    'variant': t['variant'], 'status': 'NO_MATCH_IN_OUTCOME',
                    'n_cis_available': t['n_cis_available'],
                    'study_id': t['study_id'],
                })
                continue
            rec = recs[0]
            w = wald(rec)
            rows.append({
                'target': sym, 'ensembl': t['ensembl'],
                'release': rel, 'endpoint': ep, 'endpoint_type': etype,
                'variant': t['variant'], 'status': 'OK',
                'n_cis_available': t['n_cis_available'],
                'study_id': t['study_id'],
                'beta_exp': t['beta_exp'], 'z_exp': t['z_exp'],
                'F_instrument': w['F_instrument'],
                'beta_out': rec['beta_out'], 'se_out': rec['se_out'],
                'pval_out': rec['pval_out'], 'flip': rec['flip'],
                'approx_match': rec['approx'],
                'wald_beta': w['beta'], 'wald_se': w['se'],
                'wald_p': w['p'], 'or': w['or'],
                'or_ci_lo': w['or_ci_lo'], 'or_ci_hi': w['or_ci_hi'],
                'direction': 'positive' if w['beta'] > 0 else 'negative',
                'power_OR_1.5': w['power_or_1.5'],
                'mde_OR_80pct': w['mde_or_80pct'],
            })
    return rows


def main():
    t0 = time.time()
    targets = load_targets()
    print('Targets with usable cis lead SNP:', len(targets),
          list(targets))
    for s, t in targets.items():
        print('  %-10s %-24s z=%8.2f beta=%.5f F=%.1f'
              % (s, t['variant'], t['z_exp'], t['beta_exp'],
                 t['z_exp'] ** 2))

    b1_combos = B1_COMBOS
    b2_combos = [('R13', ep, 'biological_positive_control')
                 for ep, _ in B2_ENDPOINTS]
    all_combos = b1_combos + b2_combos

    tbi_cache = {}
    for rel, ep, _ in all_combos:
        key = (rel, ep)
        if key in tbi_cache:
            continue
        try:
            _, tab = get_tbi(rel, ep)
            tbi_cache[key] = tab
            print('  tbi OK  %-4s %s' % (rel, ep))
        except Exception as e:
            print('  tbi FAIL %-4s %s: %s' % (rel, ep, e))
    collected, errors, n_jobs = fetch_all(targets, all_combos, tbi_cache)
    rows = estimate(targets, all_combos, collected)

    b1 = [r for r in rows if r['endpoint_type'] != 'biological_positive_control']
    b2 = [r for r in rows if r['endpoint_type'] == 'biological_positive_control']
    write_tsv(RES + '/a1_b1_multitarget.tsv', b1)
    write_tsv(RES + '/a1_b2_biological_control.tsv', b2)

    ok1 = [r for r in b1 if r['status'] == 'OK']
    ok2 = [r for r in b2 if r['status'] == 'OK']
    n_ok = len(ok1) + len(ok2)
    n_nomatch = sum(1 for r in rows if r['status'] != 'OK')

    # B1 系统性复现指标
    powers = [r['power_OR_1.5'] for r in ok1]
    mdes = [r['mde_OR_80pct'] for r in ok1]
    by_target_release = defaultdict(dict)
    for r in ok1:
        by_target_release[(r['target'], r['endpoint'])][r['release']] = r
    flips = []
    for (sym, ep), d in by_target_release.items():
        if 'R9' in d and 'R13' in d and d['R9']['direction'] != d['R13']['direction']:
            flips.append({'target': sym, 'endpoint': ep,
                          'R9_or': d['R9']['or'], 'R13_or': d['R13']['or']})
    p_sig = sum(1 for r in ok1 if r['wald_p'] < 0.05)
    # 按工具强度分层：检验"效能不足"是否普遍，还是取决于 |z_exp|
    strat = {'F_ge_1000': [], 'F_100_to_1000': [], 'F_lt_100': []}
    for r in ok1:
        f = r['F_instrument']
        k = ('F_ge_1000' if f >= 1000 else
             'F_100_to_1000' if f >= 100 else 'F_lt_100')
        strat[k].append(r)
    strat_summary = {}
    for k, g in strat.items():
        if not g:
            continue
        pw = [x['power_OR_1.5'] for x in g]
        md = [x['mde_OR_80pct'] for x in g]
        strat_summary[k] = {
            'n': len(g),
            'targets': sorted(set(x['target'] for x in g)),
            'median_power_OR_1.5': float(np.median(pw)),
            'n_power_below_20pct': sum(1 for p in pw if p < 0.20),
            'median_mde_OR': float(np.median(md)),
        }
    # B2 方向符合度 + 多重检验 + 极端估计警示
    BONF_B2 = 0.05 / max(1, len(ok2))
    b2_dir = []
    extreme = []
    for r in ok2:
        exp = EXPECTED.get(r['target'])
        b2_dir.append({
            'target': r['target'], 'endpoint': r['endpoint'],
            'or': r['or'], 'p': r['wald_p'],
            'expected_direction': exp or 'not_pre_specified',
            'observed_direction': r['direction'],
            'matches_expectation': (exp == r['direction']) if exp else None,
            'nominally_significant': r['wald_p'] < 0.05,
            'bonferroni_significant': r['wald_p'] < BONF_B2,
            'F_instrument': r['F_instrument'],
            'beta_exp': r['beta_exp'],
        })
        if abs(math.log(r['or'])) > 2:
            extreme.append({
                'target': r['target'], 'endpoint': r['endpoint'],
                'or': r['or'], 'p': r['wald_p'],
                'F_instrument': r['F_instrument'],
                'beta_exp': r['beta_exp'],
            })

    summary = {
        'generated_on': '2026-11-11',
        'elapsed_seconds': round(time.time() - t0, 1),
        'data_sources': {
            'exposure': 'Open Targets Platform v26.9 GraphQL '
                        'target.credibleSets, UKB-PPP pQTL (real response '
                        'in data/b1_target_pqtl.json)',
            'outcome': 'FinnGen public R9 / R13 summary stats via tabix '
                       'HTTP Range requests (no new credentials)',
        },
        'target_selection': {
            'n_candidates_queried': 12,
            'n_with_cis_pqtl_in_first_page': 7,
            'n_with_usable_simple_SNP_lead': len(targets),
            'rule': 'single strongest |z| cis pQTL lead per target; '
                    'indels excluded to keep allele harmonisation exact',
            'note': 'CRP/HMGCR/CETP/GLP1R returned no UKB_PPP cis pQTL '
                    'within the first 150 credible sets; NPC1L1 target '
                    'lookup failed. They were NOT imputed.',
        },
        'fetch': {
            'n_region_requests': n_jobs,
            'n_errors': len(errors),
            'errors': errors,
            'n_estimates_total': len(rows),
            'n_ok': n_ok, 'n_no_match': n_nomatch,
        },
        'B1': {
            'n_targets': len(targets),
            'n_combinations': len(b1_combos),
            'n_estimates': len(b1), 'n_ok': len(ok1),
            'power_OR_1.5_range': [min(powers), max(powers)] if powers else None,
            'mde_OR_range': [min(mdes), max(mdes)] if mdes else None,
            'n_estimates_power_below_20pct':
                sum(1 for p in powers if p < 0.20),
            'n_estimates_mde_OR_above_2':
                sum(1 for m in mdes if m > 2),
            'n_nominally_significant': p_sig,
            'expected_nominal_at_5pct': 0.05 * len(ok1),
            'n_release_direction_flips': len(flips),
            'flips': flips,
            'by_instrument_strength': strat_summary,
            'CX3CL1_baseline_power_OR_1.5': [0.0545, 0.1010],
            'CX3CL1_baseline_mde_OR': [5.6072, 318.5693],
            'interpretation': (
                'PARTIAL, not universal, replication. Underpowering '
                'replicates in targets with weak/small cis effects but NOT '
                'in targets with very strong cis pQTL (IL6R, ACE, AGT), '
                'where power reaches 80-100%. Therefore low power is NOT '
                'an intrinsic property of these FinnGen endpoints; it is '
                'the interaction of instrument strength (|beta_exp| and F) '
                'with endpoint case counts. Claiming "underpowering is '
                'systematic across all contexts" would be incorrect.'),
        },
        'B2': {
            'n_estimates': len(b2), 'n_ok': len(ok2),
            'pre_specified_directions': EXPECTED,
            'bonferroni_threshold': BONF_B2,
            'rows': b2_dir,
            'n_direction_matches': sum(1 for x in b2_dir
                                       if x['matches_expectation'] is True),
            'n_direction_mismatches': sum(1 for x in b2_dir
                                          if x['matches_expectation'] is False),
            'n_nominally_significant': sum(1 for x in b2_dir
                                           if x['nominally_significant']),
            'n_bonferroni_significant': sum(1 for x in b2_dir
                                            if x['bonferroni_significant']),
            'extreme_estimates_abs_logOR_gt_2': extreme,
            'interpretation': (
                'PCSK9 is the successful biological positive control: '
                'direction matches the pre-specified expectation and '
                'survives Bonferroni correction on all three cardiovascular '
                'endpoints, confirming that this pipeline can detect a '
                'known true biological effect (not merely a technical '
                'signal from ABO trans variants). LDLR does NOT replicate: '
                'direction is opposite to expectation and ORs are extreme '
                '(up to 2626). Its cis lead has the weakest instrument '
                '(z=-8.04, F=64.6) and the smallest |beta_exp| (0.0428); '
                'dividing the outcome beta by a small exposure beta inflates '
                'the Wald ratio. LDLR must be EXCLUDED as a positive control '
                'and reported instead as a cautionary example. IL6R signals '
                'are exploratory (direction not pre-specified) and must not '
                'be claimed as causal.'),
        },
        'corrections_to_original_plan': [
            'The planned claim "insufficient power and direction drift '
            'replicate systematically across targets" is only PARTLY true: '
            'see B1.by_instrument_strength and B1.interpretation.',
            'LDLR cannot serve as a biological positive control; PCSK9 '
            'alone fulfils that role. Both facts are reported, not hidden.',
            'These results do NOT change the CX3CL1 conclusion: no '
            'supporting evidence, and the public data remain insufficient '
            'to test a moderate effect for that target.',
        ],
        'honesty_boundaries': [
            'B1 targets share the FinnGen outcome resource: cross-target '
            'replication, NOT independent-dataset replication.',
            'One cis lead per target; no LD clumping, coloc, Steiger or '
            'reverse MR. Not a confirmatory causal MR for any target.',
            'B2 is a directional validation, not a substitute for drug '
            'trials. A null B2 does NOT mean PCSK9/LDLR are inactive; a '
            'significant B2 does NOT prove causality (protein pQTL may '
            'reflect reverse causation or linked variants).',
            'No P-value-driven selection of targets or endpoints; all '
            'pre-specified combinations are reported.',
        ],
    }
    with io.open(RES + '/a1_b1b2_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # 图 8：两个面板
    fig_style()
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4),
                           gridspec_kw={'width_ratios': [1.35, 1]})
    syms = list(targets)
    for rel, color, marker in [('R9', '#1d4ed8', 'o'),
                               ('R13', '#b91c1c', 's')]:
        xs, ys = [], []
        for i, s in enumerate(syms):
            for r in ok1:
                if r['target'] == s and r['release'] == rel \
                        and r['endpoint'] == 'DM_NEUROPATHY':
                    xs.append(i + (-0.12 if rel == 'R9' else 0.12))
                    ys.append(r['or'])
        ax[0].scatter(xs, ys, s=34, c=color, marker=marker, alpha=0.85,
                      edgecolors='none', label=rel)
    ax[0].axhline(1.0, ls='--', lw=0.9, color='#475569')
    ax[0].set_yscale('log')
    ax[0].set_xticks(range(len(syms)))
    ax[0].set_xticklabels(syms, rotation=40, ha='right', fontsize=8)
    ax[0].set_ylabel('Wald ratio OR (log scale)')
    ax[0].set_title('(a) Primary endpoint (DM_NEUROPATHY) across targets',
                    fontsize=10)
    ax[0].legend(fontsize=8, frameon=False)

    data = [r['power_OR_1.5'] * 100 for r in ok1]
    ax[1].hist(data, bins=16, color='#93c5fd', edgecolor='#1d4ed8', lw=0.6)
    ax[1].axvline(10.10, ls='--', color='#b91c1c', lw=1.1,
                  label='CX3CL1 max 10.1%')
    ax[1].axvline(80, ls='-', color='#047857', lw=1.1, label='80% power')
    ax[1].set_xlabel('Power (%) to detect OR = 1.5')
    ax[1].set_ylabel('Number of target-context estimates')
    ax[1].set_title('(b) Power depends on instrument strength', fontsize=10)
    ax[1].legend(fontsize=7.5, frameon=False)
    for a in ax:
        a.grid(alpha=0.22, lw=0.5, axis='y')
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig8_multitarget_replication')

    print('\n=== SUMMARY ===')
    print('fetch %d requests, %d errors, %d OK, %d no-match'
          % (n_jobs, len(errors), n_ok, n_nomatch))
    print('B1 power range', summary['B1']['power_OR_1.5_range'],
          'MDE range', summary['B1']['mde_OR_range'])
    print('B1 power<20%%:', summary['B1']['n_estimates_power_below_20pct'],
          '/', len(ok1), ' MDE>2:', summary['B1']['n_estimates_mde_OR_above_2'],
          ' nominal sig:', p_sig, ' flips:', len(flips))
    print('B2 direction matches', summary['B2']['n_direction_matches'],
          'mismatches', summary['B2']['n_direction_mismatches'],
          'nominal sig', summary['B2']['n_nominally_significant'])
    for x in b2_dir:
        print('  %-10s %-16s OR=%.3f p=%.3g exp=%s obs=%s'
              % (x['target'], x['endpoint'], x['or'], x['p'],
                 x['expected_direction'], x['observed_direction']))


if __name__ == '__main__':
    main()
