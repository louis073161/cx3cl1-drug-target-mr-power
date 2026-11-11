# -*- coding: utf-8 -*-
"""A1 阳性对照 / 效能基准：同一套 harmonise + Wald ratio 代码在大样本端点上的实跑。

目的（STROBE-MR 第 8 项 instrument validation、第 13b 项附加分析）：
回答审稿人"阴性是不是代码或 harmonise 失灵"。用与 A1 主分析完全相同
的取数、harmonise、Wald ratio 与效能计算，在 FinnGen R13 公开的大样本
端点上重跑同一批工具变量。

设计（预先指定候选集、全部如实披露、不做挑选）：
  工具变量：Open Targets v26.9 UKB-PPP CX3CL1 pQTL 9 个 SuSiE 可信集
            （3 cis / 6 trans），与 A1 主分析完全相同。
  结局：FinnGen R13 的 7 个公开大样本端点（含冠心病、冠状动脉粥样硬化、
        心绞痛、高脂血症、2 型糖尿病、硬心血管终点）。
  基准：A1 主分析 R13 DM_NEUROPATHY 的 lead cis Wald ratio SE=1.1073,
        OR=1.5 把握度 6.5%，80% 把握度最小可检出 OR=22.25。

诚信边界（必须随结果一起报告）：
  1. 这是"管线 / 工具变量验证"与"效能基准"，不是独立的阳性对照病因研究。
  2. 6 个 trans 位点（ABO / FUT2 等）存在已证实的水平多效性（Open Targets
     预计算共定位 H4 高达 1.0），任何显著结果均不得解释为 CX3CL1→疾病的
     因果效应；本分析只用来证明"管线在高效能情境下能检出真实关联信号"。
  3. 所有数值来自实际 HTTP 分段取数与运行，无模拟、无填补、无编造。
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_multi_fetch import (REGIONS, CS_JSON, get_tbi, gcs_url,   # noqa: E402
                            fetch_one_region, harmonise)

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
RES = ROOT + '/results'
Z975, Z80 = norm.ppf(0.975), norm.ppf(0.80)

# 预先指定的大样本端点（全部在 R13 公开桶中确认 .gz + .tbi 均可得）
ENDPOINTS = [
    ('I9_CHD', 'coronary heart disease'),
    ('I9_CORATHER', 'coronary atherosclerosis'),
    ('I9_ANGINA', 'angina pectoris'),
    ('I9_CVD_HARD', 'hard cardiovascular endpoint'),
    ('E4_HYPERLIPMIX', 'mixed hyperlipidaemia'),
    ('E4_HYPERLIPNAS', 'hypercholesterolaemia / hyperlipidaemia NOS'),
    ('T2D', 'type 2 diabetes'),
]

# A1 主分析基准（真实值，来自 a1_multiendpoint_mr.json）
BASELINE = {'release': 'R13', 'endpoint': 'DM_NEUROPATHY',
            'se': 1.1073, 'power_or_1.5': 0.065, 'mde_or_80pct': 22.25,
            'p': 0.764, 'or': 1.394}


def power_at_or(se, or_true):
    z = math.log(or_true) / se
    return float(norm.cdf(z - Z975) + norm.cdf(-z - Z975))


def mde_or(se):
    return math.exp((Z975 + Z80) * se)


def wald(r):
    be, se_o, bx = r['beta_out'], r['se_out'], r['beta_exp']
    if bx == 0:
        return None
    beta = be / bx
    se = se_o / abs(bx)
    z = beta / se
    p = float(2 * norm.sf(abs(z)))
    return {
        'beta': beta, 'se': se, 'z': z, 'p': p,
        'ci_lo': beta - Z975 * se, 'ci_hi': beta + Z975 * se,
        'or': math.exp(beta),
        'or_ci_lo': math.exp(beta - Z975 * se),
        'or_ci_hi': math.exp(beta + Z975 * se),
        'power_or_1.5': power_at_or(se, 1.5),
        'mde_or_80pct': mde_or(se),
        'F_instrument': r['z_exp'] ** 2,
        'se_improvement_vs_baseline': BASELINE['se'] / se,
    }


def main():
    t0 = time.time()
    cs = json.load(open(CS_JSON, encoding='utf-8'))['data']['credibleSets']['rows']
    print('Loading tabix indexes for %d endpoints ...' % len(ENDPOINTS))
    tbis = {}
    for ep, label in ENDPOINTS:
        try:
            _, t = get_tbi('R13', ep)
            tbis[ep] = t
            print('  OK  %-18s %s' % (ep, label))
        except Exception as e:
            print('  FAIL %-18s %s' % (ep, e))

    jobs = []
    for ep, _ in ENDPOINTS:
        if ep not in tbis:
            continue
        url = gcs_url('R13', ep)
        for chrom, s, e in REGIONS:
            jobs.append(('R13', ep, chrom, s, e, url, tbis[ep]))
    print('Fetching %d region blocks ...' % len(jobs))
    collected, errors = {}, []
    with ThreadPoolExecutor(max_workers=5) as ex:
        for fu in as_completed([ex.submit(fetch_one_region, j) for j in jobs]):
            rel, ep, chrom, lines, err = fu.result()
            collected.setdefault(ep, {})[chrom] = (
                collected.get(ep, {}).get(chrom, []) + lines)
            if err:
                errors.append({'endpoint': ep, 'chrom': chrom, 'error': err})

    rows = []
    for ep, label in ENDPOINTS:
        if ep not in collected:
            continue
        matched = harmonise(cs, collected[ep])
        for r in matched:
            w = wald(r)
            if not w:
                continue
            rows.append({
                'endpoint': ep, 'endpoint_label': label,
                'variant': r['variant'], 'rsids': r['rsids'],
                'chrom': r['chrom'], 'pos': r['pos'],
                'instrument_type': 'cis' if r['cis'] else 'trans',
                'nearest_gene': r['nearest_gene'],
                'exp_ref_alt': r['exp_ref_alt'], 'out_ref_alt': r['out_ref_alt'],
                'flip': r['flip'], 'approx': r['approx'],
                'column_check_ok': r['column_check_ok'],
                'beta_exp': r['beta_exp'], 'se_exp': r['se_exp'],
                'z_exp': r['z_exp'],
                'beta_out': r['beta_out'], 'se_out': r['se_out'],
                'pval_out': r['pval_out'], 'af_alt': r['af_alt'],
                'wald_beta': w['beta'], 'wald_se': w['se'],
                'wald_or': w['or'],
                'wald_ci_lo': w['or_ci_lo'], 'wald_ci_hi': w['or_ci_hi'],
                'wald_p': w['p'],
                'power_or_1.5': w['power_or_1.5'],
                'mde_or_80pct': w['mde_or_80pct'],
                'F_instrument': w['F_instrument'],
                'se_improvement_vs_baseline': w['se_improvement_vs_baseline'],
            })
        print('  %-18s matched=%d/9  raw_lines=%d' %
              (ep, len(matched),
               sum(len(v) for v in collected[ep].values())))

    # 汇总
    def summ(rs):
        if not rs:
            return {'n': 0}
        return {
            'n': len(rs),
            'n_p_lt_5e-8': sum(1 for r in rs if r['wald_p'] < 5e-8),
            'n_p_lt_0.05': sum(1 for r in rs if r['wald_p'] < 0.05),
            'min_p': min(r['wald_p'] for r in rs),
            'min_se': min(r['wald_se'] for r in rs),
            'max_power_or_1.5': max(r['power_or_1.5'] for r in rs),
            'min_mde_or': min(r['mde_or_80pct'] for r in rs),
        }

    cis = [r for r in rows if r['instrument_type'] == 'cis']
    trans = [r for r in rows if r['instrument_type'] == 'trans']
    lead = [r for r in cis if r['variant'] == '16_57379748_C_T']
    by_ep = {ep: summ([r for r in rows if r['endpoint'] == ep])
             for ep, _ in ENDPOINTS}

    out = {
        'purpose': 'pipeline / instrument validation + power benchmark',
        'design': 'pre-specified 7 R13 large-sample endpoints x 9 UKB-PPP '
                  'CX3CL1 pQTL instruments; all results reported, no selection',
        'baseline': BASELINE,
        'data_sources': {
            'exposure': 'Open Targets Platform v26.9, studyId '
                        'UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1 '
                        '(9 SuSiE credible sets, 3 cis / 6 trans)',
            'outcome': 'FinnGen DF13/R13 public summary statistics, '
                       'https://storage.googleapis.com/finngen-public-data-r13/'
                       'summary_stats/finngen_R13_<ENDPOINT>.gz (+.tbi)',
            'endpoint_manifest': 'no public manifest for R13; case counts NOT used',
        },
        'honesty_boundaries': [
            'This is a pipeline/instrument validation and power benchmark, '
            'NOT an independent positive-control aetiological study.',
            'The 6 trans instruments (ABO/FUT2 etc.) have documented horizontal '
            'pleiotropy (Open Targets precomputed coloc H4 up to 1.0); no '
            'significant result may be interpreted as a CX3CL1 causal effect.',
            'All numbers come from real HTTP range downloads and execution; '
            'no simulation, no imputation, no fabricated values.',
        ],
        'n_region_jobs': len(jobs), 'n_errors': len(errors), 'errors': errors,
        'elapsed_sec': round(time.time() - t0, 1),
        'n_rows': len(rows),
        'summary_all': summ(rows),
        'summary_cis': summ(cis),
        'summary_trans': summ(trans),
        'summary_lead_cis': summ(lead),
        'by_endpoint': by_ep,
        'rows': rows,
    }
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(RES + '/a1_positive_control.json', 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)

    keys = ['endpoint', 'endpoint_label', 'variant', 'rsids', 'chrom', 'pos',
            'instrument_type', 'nearest_gene', 'exp_ref_alt', 'out_ref_alt',
            'flip', 'approx', 'column_check_ok', 'beta_exp', 'se_exp', 'z_exp',
            'beta_out', 'se_out', 'pval_out', 'af_alt', 'wald_beta', 'wald_se',
            'wald_or', 'wald_ci_lo', 'wald_ci_hi', 'wald_p', 'power_or_1.5',
            'mde_or_80pct', 'F_instrument', 'se_improvement_vs_baseline']
    with open(RES + '/a1_positive_control.tsv', 'w', encoding='utf-8') as f:
        f.write('\t'.join(keys) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')

    print('\n=== 结果 ===')
    print('rows=%d  errors=%d  elapsed=%.1fs' %
          (len(rows), len(errors), time.time() - t0))
    for name, s in [('ALL', out['summary_all']), ('cis', out['summary_cis']),
                    ('trans', out['summary_trans']),
                    ('lead cis', out['summary_lead_cis'])]:
        print('  %-9s n=%2d  P<5e-8=%2d  P<0.05=%2d  minP=%.3g  '
              'minSE=%.4f  maxPower=%.3f  minMDE=%.2f'
              % (name, s.get('n', 0), s.get('n_p_lt_5e-8', 0),
                 s.get('n_p_lt_0.05', 0), s.get('min_p', float('nan')),
                 s.get('min_se', float('nan')),
                 s.get('max_power_or_1.5', float('nan')),
                 s.get('min_mde_or', float('nan'))))
    print('\n=== 各端点 ===')
    for ep, s in by_ep.items():
        print('  %-18s n=%2d  minP=%.3g  minSE=%.4f  maxPower=%.3f'
              % (ep, s.get('n', 0), s.get('min_p', float('nan')),
                 s.get('min_se', float('nan')),
                 s.get('max_power_or_1.5', float('nan'))))
    print('\nSaved:', RES + '/a1_positive_control.tsv',
          RES + '/a1_positive_control.json')


if __name__ == '__main__':
    main()
