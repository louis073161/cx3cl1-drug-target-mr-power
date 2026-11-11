# -*- coding: utf-8 -*-
"""
A1 第一轮真实两样本 MR（探索性、预注册前的可行性分析）

暴露：Open Targets Platform 收录的 UKB-PPP CX3CL1 血浆 pQTL 可信集
      (studyId = UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1)
结局：FinnGen DF13 / R13 公开 GWAS 汇总统计 DM_NEUROPATHY
      (Google Cloud Storage 公开桶，通过 tabix 索引 + HTTP Range 定点取数)

说明：
  * 暴露侧未提供 standardError，但提供 zScore，故 se_exp = |beta| / |z|，
    为标准做法，需在论文方法中明确交代。
  * FinnGen 汇总统计效应量以 alt 等位基因为准（官方文档明示）。
  * cis 三个位点彼此距离很近、存在 LD，本轮不做 clumping（本机无 LD 参考面板），
    故主分析采用单个最强 cis 位点的 Wald ratio；多 SNP 结果仅作探索性列示。
  * 本脚本不产生任何模拟或推算的"研究结果"，所有输入均为公开真实数据。

输出：results/a1_round1_variants.tsv、results/a1_round1_mr.json
"""
import json
import math
import os
import sys
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_tabix_fetch import fetch_region  # noqa: E402

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
CS_JSON = ROOT + '/data/ukbppp_cx3cl1_crediblesets.json'
OUT_TSV = ROOT + '/results/a1_round1_variants.tsv'
OUT_JSON = ROOT + '/results/a1_round1_mr.json'

# FinnGen R13 汇总统计列（经 −log10(p)=mlogp 与 beta/se 交叉验证确认）
# 1 #chrom 2 pos 3 ref 4 alt 5 rsids 6 nearest_genes
# 7 pval 8 mlogp 9 beta 10 sebeta 11 af_alt 12 af_alt_cases 13 af_alt_controls
COL = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'rsids': 4,
       'gene': 5, 'pval': 6, 'mlogp': 7, 'beta': 8, 'se': 9, 'af': 10}


def get_outcome(variant_id, chrom, pos, flank=8000):
    """从 FinnGen R13 拉取该位点附近真实汇总统计并做等位基因 harmonise。"""
    lines, _ = fetch_region(chrom, max(1, pos - flank), pos + flank,
                            span_mb=1)
    vid = variant_id.split('_')
    e_ref, e_alt = vid[2], vid[3]
    exact, approx = None, None
    for line in lines:
        f = line.split('\t')
        if len(f) < 11:
            continue
        try:
            p = int(f[COL['pos']])
        except ValueError:
            continue
        if abs(p - pos) > 2:
            continue
        o_ref, o_alt = f[COL['ref']], f[COL['alt']]
        beta = float(f[COL['beta']])
        se = float(f[COL['se']])
        pval = float(f[COL['pval']])
        if (o_ref, o_alt) == (e_ref, e_alt):
            flip = False
        elif (o_ref, o_alt) == (e_alt, e_ref):
            flip = True
            beta = -beta
        else:
            continue  # 等位基因不匹配，排除
        rec = {'rsids': f[COL['rsids']], 'gene': f[COL['gene']],
               'beta_out': beta, 'se_out': se, 'pval_out': pval,
               'af_alt': float(f[COL['af']]), 'flip': flip,
               'out_ref': o_ref, 'out_alt': o_alt,
               'pos_offset': p - pos, 'approx': p != pos}
        if p == pos:
            exact = rec
        elif approx is None:
            approx = rec
    return exact or approx


def wald(b_out, se_out, b_exp):
    theta = b_out / b_exp
    se = abs(se_out) / abs(b_exp)
    z = theta / se
    return {'beta': theta, 'se': se, 'ci_lo': theta - 1.96 * se,
            'ci_hi': theta + 1.96 * se,
            'p': 2 * norm.sf(abs(z)), 'or': math.exp(theta)}


def ivw(rows, egger=False):
    """固定截距 0 的 IVW；egger=True 时估计截距（MR-Egger）。"""
    n = len(rows)
    if n < (3 if egger else 1):
        return None
    w = [1.0 / (r['se_out'] ** 2) for r in rows]
    x = [r['beta_exp'] for r in rows]
    y = [r['beta_out'] for r in rows]
    sw = sum(w)
    mx = sum(wi * xi for wi, xi in zip(w, x)) / sw
    my = sum(wi * yi for wi, yi in zip(w, y)) / sw
    sxx = sum(wi * (xi - mx) ** 2 for wi, xi in zip(w, x))
    sxy = sum(wi * (xi - mx) * (yi - my) for wi, xi, yi in zip(w, x, y))
    slope = (sxy / sxx) if not egger else (sxy / sxx)
    intercept = 0.0 if not egger else (my - slope * mx)
    if egger:
        # MR-Egger 使用未加权中心化的简化实现，仅作方向性参考
        slope = sxy / sxx
        resid = [yi - (intercept + slope * xi) for xi, yi in zip(x, y)]
        se = math.sqrt(sum(wi * r ** 2 for wi, r in zip(w, resid))
                       / (n - 2) / sxx)
    else:
        se = math.sqrt(1.0 / sxx)
    z = slope / se
    q = sum(wi * (yi - slope * xi) ** 2 for wi, xi, yi in zip(w, x, y))
    return {'n_snp': n, 'beta': slope, 'se': se,
            'ci_lo': slope - 1.96 * se, 'ci_hi': slope + 1.96 * se,
            'p': 2 * norm.sf(abs(z)), 'or': math.exp(slope),
            'intercept': intercept,
            'Q': q, 'Q_df': n - 1,
            'Q_p': 1.0 - __import__('scipy.stats', fromlist=['chi2']).chi2.cdf(q, n - 1)}


def main():
    cs = json.load(open(CS_JSON, encoding='utf-8'))['data']['credibleSets']['rows']
    out_rows = []
    print('Fetching FinnGen R13 outcome betas for %d pQTL instruments ...' % len(cs))
    for r in cs:
        vid = r['variant']['id']
        chrom, pos = r['chromosome'], r['position']
        o = get_outcome(vid, chrom, pos)
        if o is None:
            print('  MISS', vid)
            continue
        if r['zScore'] is None or r['beta'] is None:
            continue
        se_exp = abs(r['beta']) / abs(r['zScore'])
        rec = {'variant': vid, 'rsids': o['rsids'], 'chrom': chrom,
               'pos': pos, 'cis': not r['isTransQtl'],
               'exp_ref_alt': '_'.join(vid.split('_')[2:]),
               'out_ref_alt': o['out_ref'] + '/' + o['out_alt'],
               'flip': o['flip'],
               'beta_exp': r['beta'], 'se_exp': se_exp, 'z_exp': r['zScore'],
               'p_exp_mantissa': r['pValueMantissa'],
               'p_exp_exponent': r['pValueExponent'],
               'beta_out': o['beta_out'], 'se_out': o['se_out'],
               'pval_out': o['pval_out'], 'af_alt': o['af_alt'],
               'nearest_gene': o['gene'], 'pos_offset': o['pos_offset'],
               'approx_match': o['approx']}
        out_rows.append(rec)
        print('  %-22s %-4s exp_b=%+.4f se=%.4f | out_b=%+.4f se=%.4f p=%.3g flip=%s'
              % (vid, 'cis' if rec['cis'] else 'trans',
                 rec['beta_exp'], rec['se_exp'], rec['beta_out'],
                 rec['se_out'], rec['pval_out'], rec['flip']))

    with open(OUT_TSV, 'w', encoding='utf-8') as f:
        keys = list(out_rows[0].keys())
        f.write('\t'.join(keys) + '\n')
        for r in out_rows:
            f.write('\t'.join(str(r[k]) for k in keys) + '\n')

    cis = [r for r in out_rows if r['cis']]
    trans = [r for r in out_rows if not r['cis']]
    cis_exact = [r for r in cis if not r['approx_match']]
    lead = max(cis_exact, key=lambda r: abs(r['z_exp'])) if cis_exact else None
    result = {
        'data_sources': {
            'exposure': 'Open Targets Platform v26.9 / UKB-PPP '
                        'UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1 (9 credible sets)',
            'outcome': 'FinnGen DF13 (public release 2026-06-02), '
                       'DM_NEUROPATHY, GCS finngen-public-data-r13',
            'se_exposure_note': 'derived as abs(beta)/abs(zScore); '
                                'standardError not provided by source',
            'ld_clumping': 'NOT PERFORMED (no LD reference panel on this machine)',
        },
        'n_instruments_fetched': len(out_rows),
        'n_cis': len(cis), 'n_cis_exact': len(cis_exact),
        'n_trans': len(trans),
        'approx_matched_variants': [r['variant'] for r in out_rows
                                    if r['approx_match']],
        'excluded_instruments': [
            {'variant': r['variant']['id'],
             'reason': 'not found in FinnGen R13 within +/-8 kb '
                       '(indel representation / variant absent)'}
            for r in cs if r['variant']['id'] not in
            {x['variant'] for x in out_rows}],
        'primary_note': 'Wald ratio uses the strongest EXACT-matched cis '
                        'instrument; indel 16_57375150_AAAG_A (largest |z|) '
                        'is absent from FinnGen R13 and cannot be used.',
        'primary_wald_ratio_lead_cis': None,
        'exploratory': {},
    }
    if lead:
        result['primary_wald_ratio_lead_cis'] = dict(
            variant=lead['variant'], rsids=lead['rsids'],
            **wald(lead['beta_out'], lead['se_out'], lead['beta_exp']))
    if cis:
        result['exploratory']['IVW_cis_3snp'] = ivw(cis)
        result['exploratory']['Egger_cis_3snp'] = ivw(cis, egger=True)
    if len(out_rows) >= 3:
        result['exploratory']['IVW_all_including_trans'] = ivw(out_rows)
        result['exploratory']['Egger_all'] = ivw(out_rows, egger=True)

    json.dump(result, open(OUT_JSON, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    print('\n=== PRIMARY (single strongest cis instrument, Wald ratio) ===')
    print(json.dumps(result['primary_wald_ratio_lead_cis'],
                     ensure_ascii=False, indent=2))
    print('\n=== EXPLORATORY (LD not clumped, trans pleiotropy likely) ===')
    print(json.dumps(result['exploratory'], ensure_ascii=False, indent=2))
    print('\nSaved:', OUT_TSV, OUT_JSON)


if __name__ == '__main__':
    main()
