# -*- coding: utf-8 -*-
"""
A1 方法学主线：跨 release × 跨端点 MR 矩阵、效能曲线与结论漂移量化。

输入全部为真实公开数据（a1_multi_fetch.py 产出），不做模拟、不填补缺失。
效能一律由各组合**实测** Wald ratio SE 驱动，不依赖任何病例数假设——
因为病例/对照口径本身正是本研究要揭示的异质性来源。

输出：
  results/a1_multiendpoint_mr.json
  results/a1_mr_matrix.tsv
  results/a1_power_curve.tsv
  results/a1_drift_summary.json
"""
import csv
import json
import math
import os
import sys
from collections import defaultdict

from scipy.stats import norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_ld_coloc import wald, ivw, weighted_median, mr_panel  # noqa: E402

MATRIX = ROOT + '/results/a1_outcome_matrix.tsv'
EDICT = ROOT + '/results/a1_endpoint_dict.tsv'
ROBUST = ROOT + '/results/a1_robustness.json'
RES = ROOT + '/results'

LEAD = '16_57379748_C_T'          # 唯一精确匹配的 cis lead 工具变量
ORS = [1.1, 1.2, 1.5, 2.0, 3.0]
Z80 = norm.ppf(0.80)               # 0.8416
Z975 = 1.959964


def read_tsv(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def power_at_or(se, or_true):
    """双侧 α=0.05 下检出 ln(OR) 真实效应的把握度。"""
    z = math.log(or_true) / se
    return float(norm.cdf(z - Z975) + norm.cdf(-z - Z975))


def min_detectable_or(se):
    return math.exp((Z975 + Z80) * se)


def analyse_combo(rel, ep, rows, clump_keep, meta):
    out = {'release': rel, 'endpoint': ep,
           'endpoint_type': meta.get('endpoint_type', ''),
           'n_cases': meta.get('n_cases', ''),
           'n_controls': meta.get('n_controls', ''),
           'n_instruments': len(rows),
           'n_approx': sum(1 for r in rows if r['approx'] == 'True'),
           'methods': {}, 'power': {}, 'panels': {}}

    lead = [r for r in rows if r['variant'] == LEAD]
    if lead:
        out['methods']['Wald ratio (lead cis)'] = wald(lead[0])

    clumped = [r for r in rows if r['variant'] in clump_keep]
    if clumped:
        out['panels']['clumped_r2_0.001'] = mr_panel(clumped, 'clumped')
    out['panels']['all_available'] = mr_panel(rows, 'all') if rows else None

    # 效能：完全由实测主分析 SE 驱动
    p = out['methods'].get('Wald ratio (lead cis)')
    if p:
        se = p['se']
        out['power'] = {
            'se': se,
            'power_by_OR': {'OR_%.1f' % o: power_at_or(se, o) for o in ORS},
            'min_detectable_OR_at_80pct': min_detectable_or(se),
            'se_fold_needed_for_OR1.5_at_80pct':
                se / (math.log(1.5) / (Z975 + Z80)),
            'or_range_of_95CI': [math.exp(p['ci_lo']), math.exp(p['ci_hi'])],
            'or_excluded_above': math.exp(p['ci_hi']),
            'or_excluded_below': math.exp(p['ci_lo']),
        }
    return out


def drift(group, label):
    """对一组组合的主分析 Wald ratio 做漂移量化。"""
    vals = []
    for c in group:
        p = c['methods'].get('Wald ratio (lead cis)')
        if p:
            vals.append({'release': c['release'], 'endpoint': c['endpoint'],
                         'beta': p['beta'], 'or': p['or'], 'p': p['p'],
                         'se': p['se']})
    if not vals:
        return {'label': label, 'n': 0}
    bs = [v['beta'] for v in vals]
    ps = [v['p'] for v in vals]
    os_ = [v['or'] for v in vals]
    signs = [1 if b > 0 else -1 for b in bs]
    flips = sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])
    n_sig = sum(1 for x in ps if x < 0.05)
    return {
        'label': label, 'n': len(vals),
        'beta_min': min(bs), 'beta_max': max(bs),
        'beta_range': max(bs) - min(bs),
        'or_min': min(os_), 'or_max': max(os_),
        'or_ratio_max_over_min': max(os_) / min(os_),
        'p_min': min(ps), 'p_max': max(ps),
        'n_nominal_p_lt_0.05': n_sig,
        'expected_nominal_at_0.05': 0.05 * len(vals),
        'n_sign_flips_in_listed_order': flips,
        'n_positive': sum(1 for s in signs if s > 0),
        'n_negative': sum(1 for s in signs if s < 0),
        'conclusion_classes': {
            'negative_and_nonsignificant': sum(1 for v in vals if v['p'] >= 0.05),
            'nominally_significant': n_sig,
        },
        'values': vals,
    }


def main():
    raw = read_tsv(MATRIX)
    robust = json.load(open(ROBUST, encoding='utf-8'))
    clump_keep = set(robust['clumping']['r2_0.001']['kept_ids'])
    ed = {(r['release'], r['endpoint']): r for r in read_tsv(EDICT)}

    # TSV 读回后把数值字段还原成 float / bool
    conv = {}
    for r in raw:
        key = (r['release'], r['endpoint'])
        x = dict(r)
        for k in ['beta_exp', 'se_exp', 'z_exp', 'beta_out', 'se_out',
                  'pval_out', 'af_alt', 'pos_offset']:
            x[k] = num(x[k])
        x['pos'] = int(x['pos'])
        x['cis'] = x['cis'] == 'True'
        x['flip'] = x['flip'] == 'True'
        x['approx'] = x['approx'] == 'True'
        conv.setdefault(key, []).append(x)

    combos = []
    for key in sorted(conv, key=lambda k: (int(k[0][1:]), k[1])):
        rel, ep = key
        combos.append(analyse_combo(rel, ep, conv[key], clump_keep,
                                    ed.get(key, {})))

    # ---- 长表：每组合 × 每方法 ----
    mrows = []
    for c in combos:
        for m, v in c['methods'].items():
            mrows.append({'release': c['release'], 'endpoint': c['endpoint'],
                          'endpoint_type': c['endpoint_type'],
                          'n_cases': c['n_cases'], 'n_controls': c['n_controls'],
                          'n_snp': v['n_snp'], 'method': m,
                          'beta': v['beta'], 'se': v['se'],
                          'ci_lo': v['ci_lo'], 'ci_hi': v['ci_hi'],
                          'p': v['p'], 'or': v['or']})
        for pname, panel in c['panels'].items():
            if not panel:
                continue
            for m, v in panel['methods'].items():
                mrows.append({'release': c['release'], 'endpoint': c['endpoint'],
                              'endpoint_type': c['endpoint_type'],
                              'n_cases': c['n_cases'], 'n_controls': c['n_controls'],
                              'n_snp': v['n_snp'],
                              'method': '%s | %s' % (pname, m),
                              'beta': v['beta'], 'se': v['se'],
                              'ci_lo': v['ci_lo'], 'ci_hi': v['ci_hi'],
                              'p': v['p'], 'or': v['or']})
    with open(RES + '/a1_mr_matrix.tsv', 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(mrows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(mrows)

    # ---- 效能曲线长表 ----
    pows = []
    for c in combos:
        if not c['power']:
            continue
        for k, v in c['power']['power_by_OR'].items():
            pows.append({'release': c['release'], 'endpoint': c['endpoint'],
                         'endpoint_type': c['endpoint_type'],
                         'se': c['power']['se'], 'scenario': k, 'power': v})
    with open(RES + '/a1_power_curve.tsv', 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(pows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(pows)

    # ---- 漂移量化 ----
    def pick(pred):
        return [c for c in combos if pred(c)]

    groups = {
        'primary_DM_NEUROPATHY_across_R9_R11_R12_R13':
            pick(lambda c: c['endpoint'] == 'DM_NEUROPATHY'),
        'broad_G6_POLYOTHUNS_across_R9_to_R13':
            pick(lambda c: c['endpoint'] == 'G6_POLYOTHUNS'),
        'negative_control_microvascular_all':
            pick(lambda c: c['endpoint_type'] == 'negative_ctrl_microvascular'),
        'R13_strictness_axis_retinopathy':
            pick(lambda c: c['release'] == 'R13' and 'RETINOPATHY' in c['endpoint']),
        'all_19_combos': combos,
    }
    drift_out = {k: drift(v, k) for k, v in groups.items()}
    json.dump(drift_out, open(RES + '/a1_drift_summary.json', 'w',
                              encoding='utf-8'), ensure_ascii=False, indent=2)

    result = {
        'design': {
            'exposure': 'Open Targets v26.9 UKB-PPP CX3CL1 pQTL '
                        'UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1 (9 credible sets)',
            'outcome': 'FinnGen public summary statistics, 5 releases '
                       '(R9/R10/R11/R12/R13) x multiple endpoints',
            'primary_instrument': LEAD,
            'clumping_threshold': 'r2<0.001, 1000 Genomes EUR (from a1_robustness.json)',
            'power_note': 'All power figures are driven by the OBSERVED '
                          'Wald ratio SE of each combo; no case-count assumption.',
            'independence_note': 'FinnGen releases are nested expansions of '
                                 'the SAME cohort, not independent studies; '
                                 'drift quantifies version sensitivity, '
                                 'not replicability.',
            'no_simulation': True,
        },
        'n_combos': len(combos),
        'combos': combos,
    }
    json.dump(result, open(RES + '/a1_multiendpoint_mr.json', 'w',
                           encoding='utf-8'), ensure_ascii=False, indent=2)

    # ---- 控制台摘要 ----
    print('=== 主分析 Wald ratio（lead cis）跨组合 ===')
    for c in combos:
        p = c['methods'].get('Wald ratio (lead cis)')
        pw = c['power']
        if p:
            print('  %-4s %-26s n=%d  beta=%+.4f SE=%.4f P=%.4f OR=%.3f '
                  '| pow(OR1.5)=%.3f  minOR80=%.2f'
                  % (c['release'], c['endpoint'], c['n_instruments'],
                     p['beta'], p['se'], p['p'], p['or'],
                     pw['power_by_OR']['OR_1.5'],
                     pw['min_detectable_OR_at_80pct']))
    print('\n=== 漂移 ===')
    for k, d in drift_out.items():
        if d.get('n'):
            print('  %-52s n=%2d OR=[%.3f,%.3f] ratio=%.2f '
                  'P=[%.3g,%.3g] nominal_sig=%d/%d'
                  % (k, d['n'], d['or_min'], d['or_max'],
                     d['or_ratio_max_over_min'], d['p_min'], d['p_max'],
                     d['n_nominal_p_lt_0.05'], d['n']))
    print('\nSaved:', RES + '/a1_multiendpoint_mr.json',
          RES + '/a1_mr_matrix.tsv', RES + '/a1_power_curve.tsv',
          RES + '/a1_drift_summary.json')


if __name__ == '__main__':
    main()
