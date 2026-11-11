# -*- coding: utf-8 -*-
"""基于 SuSiE 95% 可信集 PIP 的、先验无关的共定位核查（CLPP）。

数据源：
  * Open Targets Platform v26.9 GraphQL credibleSets.locus.rows
    （3 个 CX3CL1 cis 可信集的逐变异 PIP / beta / logBF，真实数据）
  * FinnGen R13 DM_NEUROPATHY chr16 区段 3,140 个变异（真实数据）

CLPP = Σ_i PIP_pQTL(i) × PIP_DPN(i)（eCAVIAR 风格）。
不依赖 p1/p2/p12 先验，用来检验 ABF coloc 对 p12 敏感这一局限。
"""
import json
import math
import os
import sys
from scipy.special import logsumexp

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_ld_coloc import (load_fg_region, log_abf, coloc_abf,
                         FG_TSV)  # noqa: E402

LOCUS = ROOT + '/data/ukbppp_cx3cl1_locus.json'
ROB = ROOT + '/results/a1_robustness.json'
OUT = ROOT + '/results/a1_pip_coloc.json'
W = 0.2 ** 2


def parse_vid(v):
    p = v.split('_')
    return {'chrom': p[0], 'pos': int(p[1]), 'ref': p[2], 'alt': p[3]}


def main():
    fg = load_fg_region()
    for r in fg:
        r['logABF'] = log_abf(r['beta'], r['se'], W)
    tot = logsumexp([r['logABF'] for r in fg])
    for r in fg:
        r['PIP'] = math.exp(r['logABF'] - tot)
    bypos = {}
    for r in fg:
        bypos.setdefault(r['pos'], []).append(r)

    cs = json.load(open(LOCUS, encoding='utf-8'))['data']['credibleSets']['rows']
    cis = [x for x in cs if not x['isTransQtl']]
    out = {'method': 'PIP-based colocalisation (CLPP), prior-independent',
           'prior_W_for_DPN_ABF': W,
           'n_FinnGen_region_variants': len(fg),
           'note': 'Open Targets locus.logBF 用其自身精细定位先验；'
                   '与 DPN 侧 Wakefield ABF 的先验强度不完全一致，'
                   '故 ABF 版本仅作敏感性分析',
           'credible_sets': []}
    all_cs = []
    for c in cis:
        rows = c['locus']['rows']
        rec = {'studyLocusId': c['studyLocusId'], 'lead_pos': c['position'],
               'cs_index': c['credibleSetIndex'],
               'n_cs_variants': len(rows), 'variants': []}
        for v in rows:
            q = parse_vid(v['variant']['id'])
            matches = bypos.get(q['pos'], [])
            hit = None
            for m in matches:
                if (m['ref'], m['alt']) == (q['ref'], q['alt']):
                    hit = m
                    break
            if hit is None and matches:
                hit = matches[0]
            item = {'variant': v['variant']['id'],
                    'pos': q['pos'], 'PIP_pQTL': v['posteriorProbability'],
                    'beta_pQTL': v['beta'], 'logBF_pQTL': v['logBF'],
                    'in_FinnGen': hit is not None,
                    'PIP_DPN': hit['PIP'] if hit else None,
                    'logABF_DPN': hit['logABF'] if hit else None,
                    'pval_DPN': hit['pval'] if hit else None}
            rec['variants'].append(item)
            all_cs.append(item)
        rec['CLPP'] = sum(x['PIP_pQTL'] * (x['PIP_DPN'] or 0)
                          for x in rec['variants'])
        rec['sum_PIP_DPN_over_CS'] = sum((x['PIP_DPN'] or 0)
                                         for x in rec['variants'])
        rec['max_PIP_DPN_in_CS'] = max((x['PIP_DPN'] or 0)
                                       for x in rec['variants'])
        out['credible_sets'].append(rec)

    out['CLPP_all_three_cis_CS'] = sum(r['CLPP'] for r in out['credible_sets'])
    out['sum_PIP_DPN_all_CS_variants'] = sum(
        (x['PIP_DPN'] or 0) for x in all_cs)

    # 敏感性：用 OT 精细定位 logBF 作 trait1 ABF，DPN Wakefield ABF 作 trait2
    shared = [x for x in all_cs if x['in_FinnGen']]
    la1 = [x['logBF_pQTL'] for x in shared]
    la2 = [x['logABF_DPN'] for x in shared]
    out['abf_coloc_over_CS_variants'] = coloc_abf(la1, la2)
    out['abf_coloc_n_variants'] = len(shared)
    out['abf_coloc_variants'] = [x['variant'] for x in shared]

    json.dump(out, open(OUT, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    print('=== PIP-based colocalisation ===')
    for r in out['credible_sets']:
        print('  CS idx=%s lead=%d n=%d  CLPP=%.3g  sumPIP_DPN=%.3g  '
              'maxPIP_DPN=%.3g' % (r['cs_index'], r['lead_pos'],
                                   r['n_cs_variants'], r['CLPP'],
                                   r['sum_PIP_DPN_over_CS'],
                                   r['max_PIP_DPN_in_CS']))
    print('  CLPP(all 3 CS) = %.3g' % out['CLPP_all_three_cis_CS'])
    print('  ABF coloc over %d CS variants: PP.H4=%.4g PP.H1=%.4g'
          % (len(shared), out['abf_coloc_over_CS_variants']['PP.H4'],
             out['abf_coloc_over_CS_variants']['PP.H1']))
    print('Saved:', OUT)


if __name__ == '__main__':
    main()
