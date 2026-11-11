# -*- coding: utf-8 -*-
"""
A1 稳健性核查：真实 EUR LD clumping + 共定位 + 重跑 MR

输入（全部为真实公开数据）：
  1) Open Targets Platform v26.9 / UKB-PPP CX3CL1 pQTL 9 个 SuSiE 可信集
  2) FinnGen DF13 / R13 DM_NEUROPATHY 汇总统计（chr16 CX3CL1 区段 3,140 行）
  3) Ensembl REST 1000 Genomes phase 3 EUR 群体 LD（r²）
  4) Open Targets 预计算共定位（COLOC_PIP_ECAVIAR）
  5) Risteys 端点页面报告的病例数（用于效能讨论，注明 release 未明示）

不做任何模拟；缺失数据一律标注，不填补。
"""
import csv
import json
import math
import os
import ssl
import sys
import urllib.request
from scipy.special import logsumexp
from scipy.stats import norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')

CS_JSON = ROOT + '/data/ukbppp_cx3cl1_crediblesets.json'
VAR_TSV = ROOT + '/results/a1_round1_variants.tsv'
FG_TSV = ROOT + '/results/fg_R13_DPN_chr16_CX3CL1.tsv'
ROUND1_JSON = ROOT + '/results/a1_round1_mr.json'
COLOC_ALL = ROOT + '/data/ukbppp_cx3cl1_coloc_all.json'
LD_JSON = ROOT + '/data/ld_ensembl_1000g_eur.json'
OUT_JSON = ROOT + '/results/a1_robustness.json'

COL = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'rsids': 4,
       'gene': 5, 'pval': 6, 'mlogp': 7, 'beta': 8, 'se': 9, 'af': 10}


# ---------------------------------------------------------------- LD 参考
def fetch_ld(rsids, chroms):
    """Ensembl REST 1000GENOMES:phase_3:EUR 单变异 LD 全表（真实数据）。"""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    out = {}
    for rs, ch in zip(rsids, chroms):
        url = ('https://rest.ensembl.org/ld/human/%s/'
               '1000GENOMES:phase_3:EUR' % rs)
        try:
            req = urllib.request.Request(
                url, headers={'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=60, context=ctx) as r:
                d = json.loads(r.read().decode('utf-8'))
            out[rs] = [{'v2': x['variation2'], 'r2': float(x['r2']),
                        'dprime': float(x['d_prime'])}
                       for x in d if isinstance(x, dict) and 'r2' in x]
            print('  LD %-16s %4d pairs' % (rs, len(out[rs])))
        except Exception as e:
            out[rs] = []
            print('  LD %-16s ERROR %s' % (rs, e))
    return out


def build_ld_matrix(cs, ld_raw):
    """构建 r² 矩阵：同染色体取 Ensembl 实测；未出现在 LD 表按 0 处理并标注。"""
    n = len(cs)
    ids = [r['variant']['id'] for r in cs]
    rss = [r['variant']['rsIds'][0] for r in cs]
    chs = [r['chromosome'] for r in cs]
    M = [[0.0] * n for _ in range(n)]
    unknown = []
    for i in range(n):
        for j in range(i + 1, n):
            if chs[i] != chs[j]:
                r2 = 0.0  # 不同染色体，结构上无 LD
                kind = 'different_chromosome'
            else:
                pairs = {x['v2']: x['r2'] for x in ld_raw.get(rss[i], [])}
                pairs2 = {x['v2']: x['r2'] for x in ld_raw.get(rss[j], [])}
                v = pairs.get(rss[j], pairs2.get(rss[i]))
                if v is None:
                    r2 = 0.0
                    kind = 'not_in_1000G_LD_panel__treated_as_0'
                    unknown.append([ids[i], ids[j]])
                else:
                    r2 = float(v)
                    kind = 'ensembl_1000G_EUR'
            M[i][j] = M[j][i] = r2
            if r2 > 0:
                print('  r2 %-22s x %-22s = %.6f (%s)'
                      % (ids[i], ids[j], r2, kind))
    return {'ids': ids, 'rsids': rss, 'chroms': chs, 'matrix': M,
            'unknown_pairs': unknown}


def clump(cs, M, threshold, unknown_fill=0.0, unknown_pairs=None):
    """按 |z| 降序贪心 clumping；保留 r² < threshold 的位点。"""
    ids, mat = M['ids'], M['matrix']
    idx = {v: i for i, v in enumerate(ids)}
    extra = {(a, b): unknown_fill for a, b in (unknown_pairs or [])}
    order = sorted(range(len(cs)),
                   key=lambda i: -abs(cs[i]['zScore'] or 0))
    kept, reason = [], {}
    for i in order:
        drop = None
        for k in kept:
            r2 = mat[i][k]
            a, b = ids[i], ids[k]
            if (a, b) in extra:
                r2 = extra[(a, b)]
            elif (b, a) in extra:
                r2 = extra[(b, a)]
            if r2 >= threshold:
                drop = (r2, ids[k])
                break
        if drop:
            reason[ids[i]] = 'r2=%.4f with %s' % drop
        else:
            kept.append(i)
    return {'kept_ids': [ids[i] for i in kept],
            'dropped': reason, 'threshold': threshold,
            'unknown_fill': unknown_fill}


# ---------------------------------------------------------------- MR 方法
def wald(row):
    b, se = row['beta_out'] / row['beta_exp'], abs(row['se_out'] / row['beta_exp'])
    z = b / se
    return {'method': 'Wald ratio', 'n_snp': 1, 'beta': b, 'se': se,
            'ci_lo': b - 1.96 * se, 'ci_hi': b + 1.96 * se,
            'p': float(2 * norm.sf(abs(z))), 'or': math.exp(b)}


def ivw(rows, egger=False):
    """IVW（截距固定为 0）与 MR-Egger（自由截距）各自的加权最小二乘解。

    注意：首轮 a1_round1_mr.py 把两者都实现成"自由截距"回归，导致斜率相同、
    MR-Egger 并非真正的 MR-Egger；此处按标准公式重算。
    """
    n = len(rows)
    if n < (3 if egger else 1):
        return None
    w = [1.0 / r['se_out'] ** 2 for r in rows]
    x = [r['beta_exp'] for r in rows]
    y = [r['beta_out'] for r in rows]
    from scipy.stats import chi2
    if not egger:
        # IVW: 过原点，权重 w = 1/se_y^2
        denom = sum(a * b * b for a, b in zip(w, x))
        if denom <= 0:
            return None
        slope = sum(a * b * c for a, b, c in zip(w, x, y)) / denom
        se = math.sqrt(1.0 / denom)
        intercept = 0.0
        resid = [yi - slope * xi for xi, yi in zip(x, y)]
        df = n - 1
    else:
        # MR-Egger: 加权回归含截距，中心化后求斜率
        sw = sum(w)
        mx = sum(a * b for a, b in zip(w, x)) / sw
        my = sum(a * b for a, b in zip(w, y)) / sw
        sxx = sum(a * (b - mx) ** 2 for a, b in zip(w, x))
        sxy = sum(a * (b - mx) * (c - my) for a, b, c in zip(w, x, y))
        if sxx <= 0:
            return None
        slope = sxy / sxx
        intercept = my - slope * mx
        resid = [yi - (intercept + slope * xi) for xi, yi in zip(x, y)]
        sigma2 = sum(a * r ** 2 for a, r in zip(w, resid)) / (n - 2)
        se = math.sqrt(sigma2 / sxx)
        df = n - 2
    q = sum(a * r ** 2 for a, r in zip(w, resid))
    out = {'method': 'MR-Egger' if egger else 'IVW (fixed)',
           'n_snp': n, 'beta': slope, 'se': se,
           'ci_lo': slope - 1.96 * se, 'ci_hi': slope + 1.96 * se,
           'p': float(2 * norm.sf(abs(slope / se))),
           'or': math.exp(slope), 'intercept': intercept,
           'Q': q, 'Q_df': df, 'Q_p': float(1 - chi2.cdf(q, df))}
    if egger:
        out['intercept_se'] = math.sqrt(
            sum(a * r ** 2 for a, r in zip(w, resid)) / (n - 2)
            * (1.0 / sum(w) + mx ** 2 / sxx))
        out['intercept_p'] = float(
            2 * norm.sf(abs(intercept / out['intercept_se'])))
    return out


def weighted_median(rows, n_boot=2000, seed=20261110):
    """加权中位数（Bowden 2016），SE 用 bootstrap。"""
    import random
    if len(rows) < 3:
        return None
    ratios = [r['beta_out'] / r['beta_exp'] for r in rows]
    w = [1.0 / (r['se_out'] / abs(r['beta_exp'])) ** 2 for r in rows]

    def wm(vals, wts):
        z = sorted(zip(vals, wts))
        tot = sum(wts)
        c = 0.0
        for v, wt in z:
            c += wt
            if c >= tot / 2:
                return v
        return z[-1][0]

    est = wm(ratios, w)
    rng = random.Random(seed)
    bs = []
    for _ in range(n_boot):
        ix = [rng.randrange(len(rows)) for _ in range(len(rows))]
        bs.append(wm([ratios[i] for i in ix], [w[i] for i in ix]))
    se = (sum((b - sum(bs) / len(bs)) ** 2 for b in bs)
          / (len(bs) - 1)) ** 0.5
    return {'method': 'Weighted median (bootstrap SE)', 'n_snp': len(rows),
            'beta': est, 'se': se, 'ci_lo': est - 1.96 * se,
            'ci_hi': est + 1.96 * se,
            'p': float(2 * norm.sf(abs(est / se))) if se > 0 else None,
            'or': math.exp(est), 'n_boot': n_boot}


def mr_panel(rows, label):
    out = {'label': label, 'n_snp': len(rows), 'variants': [r['variant'] for r in rows],
           'methods': {}}
    if len(rows) == 1:
        out['methods']['Wald ratio'] = wald(rows[0])
    a = ivw(rows, False)
    if a:
        out['methods']['IVW'] = a
        if a['Q_p'] < 0.05:
            a['random_effects_se'] = a['se'] * math.sqrt(max(1.0, a['Q'] / a['Q_df']))
            a['note'] = 'Q_p<0.05, 报告随机效应 SE'
    e = ivw(rows, True)
    if e:
        out['methods']['MR-Egger'] = e
    wm = weighted_median(rows)
    if wm:
        out['methods']['Weighted median'] = wm
    if len(rows) > 2:
        loo = []
        for k in range(len(rows)):
            sub = [r for i, r in enumerate(rows) if i != k]
            r2 = ivw(sub, False)
            loo.append({'removed': rows[k]['variant'],
                        'beta': r2['beta'], 'p': r2['p']})
        out['leave_one_out'] = loo
    out['F_per_instrument'] = {r['variant']: round(r['z_exp'] ** 2, 1)
                               for r in rows}
    out['direction_check_z2'] = {
        r['variant']: {'z_exp2': round(r['z_exp'] ** 2, 2),
                       'z_out2': round((r['beta_out'] / r['se_out']) ** 2, 3),
                       'exposure_gt_outcome':
                           r['z_exp'] ** 2 > (r['beta_out'] / r['se_out']) ** 2}
        for r in rows}
    return out


# ---------------------------------------------------------------- ABF / coloc
def log_abf(beta, se, W):
    V = se * se
    r = W / (V + W)
    z2 = (beta / se) ** 2
    return 0.5 * math.log(1.0 - r) + 0.5 * z2 * r


def coloc_abf(la1, la2, p1=1e-4, p2=1e-4, p12=1e-5):
    """Giambartolomei/Wallace 近似 Bayes factor 共定位（单因果变异假设）。"""
    n = len(la1)
    s1, s2, s4 = logsumexp(la1), logsumexp(la2), logsumexp(
        [a + b for a, b in zip(la1, la2)])
    d = s4 - s1 - s2
    log_s3 = s1 + s2 + math.log(max(1e-300, 1.0 - math.exp(min(0.0, d))))
    lh = [0.0,
          math.log(p1) + s1 - math.log(n),
          math.log(p2) + s2 - math.log(n),
          math.log(p1) + math.log(p2) + log_s3 - math.log(n * (n - 1)),
          math.log(p12) + s4 - math.log(n)]
    den = logsumexp(lh)
    pp = [math.exp(x - den) for x in lh]
    return {'n_snp': n, 'PP.H0': pp[0], 'PP.H1': pp[1], 'PP.H2': pp[2],
            'PP.H3': pp[3], 'PP.H4': pp[4],
            'logS1': s1, 'logS2': s2, 'logS4': s4,
            'p1': p1, 'p2': p2, 'p12': p12}


def load_fg_region():
    rows = []
    with open(FG_TSV, encoding='utf-8') as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if len(x) < 11:
                continue
            try:
                pos = int(x[COL['pos']])
                beta = float(x[COL['beta']])
                se = float(x[COL['se']])
                pval = float(x[COL['pval']])
            except ValueError:
                continue
            if se <= 0 or not (0 < pval <= 1):
                continue
            rows.append({'chrom': x[COL['chrom']], 'pos': pos,
                         'ref': x[COL['ref']], 'alt': x[COL['alt']],
                         'rsids': x[COL['rsids']], 'gene': x[COL['gene']],
                         'beta': beta, 'se': se, 'pval': pval,
                         'af_alt': float(x[COL['af']])})
    return rows


def main():
    cs = json.load(open(CS_JSON, encoding='utf-8'))['data']['credibleSets']['rows']
    with open(VAR_TSV, encoding='utf-8') as f:
        var = list(csv.DictReader(f, delimiter='\t'))
    for r in var:
        for k in ('beta_exp', 'se_exp', 'z_exp', 'beta_out', 'se_out',
                  'pval_out', 'af_alt', 'pos', 'pos_offset'):
            r[k] = float(r[k])
        r['cis'] = (r['cis'] == 'True')
    by_id = {r['variant']: r for r in var}
    fg = load_fg_region()
    print('Loaded: %d credible sets, %d harmonised instruments, '
          '%d FinnGen region variants' % (len(cs), len(var), len(fg)))

    # 1) 真实 LD
    print('\n[1] Fetching real 1000G EUR LD from Ensembl REST ...')
    ld_raw = fetch_ld([r['variant']['rsIds'][0] for r in cs],
                      [r['chromosome'] for r in cs])
    json.dump(ld_raw, open(LD_JSON, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    M = build_ld_matrix(cs, ld_raw)
    with open(ROOT + '/results/a1_ld_matrix.tsv', 'w', encoding='utf-8') as f:
        f.write('variant\t' + '\t'.join(M['ids']) + '\n')
        for i, v in enumerate(M['ids']):
            f.write(v + '\t' + '\t'.join('%.6f' % x for x in M['matrix'][i]) + '\n')

    # 2) clumping 多阈值 + 对未知 LD 的敏感性
    print('\n[2] LD clumping ...')
    clumps = {}
    for t in (0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8):
        c = clump(cs, M, t)
        c['available_after_harmonise'] = [v for v in c['kept_ids'] if v in by_id]
        c['cis_available'] = [v for v in c['available_after_harmonise']
                              if by_id[v]['cis']]
        clumps['r2_%.3g' % t] = c
        print('  r2<%-5g kept=%d  available=%d  cis_available=%d  %s'
              % (t, len(c['kept_ids']), len(c['available_after_harmonise']),
                 len(c['cis_available']), ','.join(c['cis_available'])))
    sens = clump(cs, M, 0.001, unknown_fill=1.0,
                 unknown_pairs=[tuple(x) for x in M['unknown_pairs']])
    sens['available_after_harmonise'] = [v for v in sens['kept_ids'] if v in by_id]
    sens['cis_available'] = [v for v in sens['available_after_harmonise']
                             if by_id[v]['cis']]
    clumps['r2_0.001__unknown_LD_assumed_1.0'] = sens
    with open(ROOT + '/results/a1_clumping.tsv', 'w', encoding='utf-8') as f:
        f.write('setting\tkept_exposure\tdropped_by_LD\tavailable_harmonised\tcis_available\n')
        for k, c in clumps.items():
            f.write('%s\t%s\t%s\t%s\t%s\n' %
                    (k, ';'.join(c['kept_ids']),
                     ';'.join('%s(%s)' % (a, b) for a, b in c['dropped'].items()),
                     ';'.join(c['available_after_harmonise']),
                     ';'.join(c['cis_available'])))

    # 3) 重跑 MR
    print('\n[3] Re-running MR on clumped sets ...')
    mr = {}
    sets = {
        'primary_cis_lead_only': [by_id['16_57379748_C_T']],
        'clumped_r2_0.001_cis_only': [by_id[v] for v in
                                      clumps['r2_0.001']['cis_available']],
        'clumped_r2_0.001_all': [by_id[v] for v in
                                 clumps['r2_0.001']['available_after_harmonise']],
        'clumped_r2_0.1_cis_only': [by_id[v] for v in
                                    clumps['r2_0.1']['cis_available']],
        'clumped_r2_0.2_cis_only': [by_id[v] for v in
                                    clumps['r2_0.2']['cis_available']],
        'all_8_unclumped': var,
    }
    for k, rows in sets.items():
        if rows:
            mr[k] = mr_panel(rows, k)
            m = mr[k]['methods'].get('IVW') or mr[k]['methods'].get('Wald ratio')
            print('  %-32s n=%d  beta=%+.4f  se=%.4f  p=%.4g'
                  % (k, len(rows), m['beta'], m['se'], m['p']))

    # 4) 区段 ABF / PIP（DPN 侧）
    print('\n[4] Region-wide ABF / PIP scan (FinnGen DPN) ...')
    W = 0.2 ** 2
    for r in fg:
        r['logABF'] = log_abf(r['beta'], r['se'], W)
    tot = logsumexp([r['logABF'] for r in fg])
    for r in fg:
        r['PIP'] = math.exp(r['logABF'] - tot)
    top = sorted(fg, key=lambda r: -r['PIP'])[:15]
    for r in top[:8]:
        print('  top DPN PIP  %d  %-14s p=%.3g  PIP=%.4f'
              % (r['pos'], r['rsids'], r['pval'], r['PIP']))
    lead_pos = [57379748, 57374992, 57375150]
    overlap = sum(r['PIP'] for r in fg if r['pos'] in lead_pos)
    srt = sorted(fg, key=lambda r: -r['PIP'])
    cum, cs95, i = 0.0, [], 0
    while i < len(srt) and cum < 0.95:
        cum += srt[i]['PIP']
        cs95.append(srt[i])
        i += 1
    minp = min(r['pval'] for r in fg)
    region = {'n_variants': len(fg), 'min_pval': minp,
              'top10': [{k: r[k] for k in ('pos', 'rsids', 'gene', 'pval',
                                           'beta', 'se', 'PIP')} for r in top[:10]],
              'PIP_sum_at_3_cis_pQTL_leads': overlap,
              'DPN_95pct_credible_set_size': len(cs95),
              'DPN_95pct_credible_set_positions': [r['pos'] for r in cs95[:30]],
              'any_pQTL_lead_in_DPN_95CS': any(r['pos'] in lead_pos for r in cs95),
              'n_nominal_p_lt_0.05': sum(1 for r in fg if r['pval'] < 0.05),
              'n_p_lt_1e-5': sum(1 for r in fg if r['pval'] < 1e-5)}
    print('  min p=%.3g  3-cis-lead PIP sum=%.6g  95%%CS size=%d  '
          'lead in CS=%s' % (minp, overlap, len(cs95),
                             region['any_pQTL_lead_in_DPN_95CS']))
    with open(ROOT + '/results/a1_dpn_region_abf.tsv', 'w', encoding='utf-8') as f:
        f.write('chrom\tpos\trsids\tgene\tbeta\tse\tpval\taf_alt\tlogABF\tPIP\n')
        for r in sorted(fg, key=lambda x: -x['PIP']):
            f.write('\t'.join(str(v) for v in
                             [r['chrom'], r['pos'], r['rsids'], r['gene'],
                              r['beta'], r['se'], r['pval'], r['af_alt'],
                              round(r['logABF'], 6), round(r['PIP'], 9)]) + '\n')

    # 5) coloc：A 对称（2 个共享 cis 位点）、B 全区段（保守假设）
    print('\n[5] Colocalisation (Wakefield ABF) ...')
    shared = [by_id[v] for v in ('16_57379748_C_T', '16_57374992_T_C')]
    coloc = {'method': 'Wakefield ABF / Giambartolomei coloc (single causal variant)',
             'prior_W': W, 'coloc_A_symmetric_shared_cis': {},
             'coloc_B_region_conservative': {}, 'sensitivity': []}
    la1 = [log_abf(r['beta_exp'], r['se_exp'], W) for r in shared]
    la2 = [log_abf(r['beta_out'], r['se_out'], W) for r in shared]
    coloc['coloc_A_symmetric_shared_cis'] = coloc_abf(la1, la2)
    # B：以 3,140 个 DPN 变异为集合，pQTL 仅在 2 个共享 lead 有信号
    #    其余位点 pQTL ABF 置为 1（无信息），为保守敏感性分析
    la1_b, la2_b = [], []
    for r in fg:
        v = next((x for x in shared if int(x['pos']) == r['pos']), None)
        la1_b.append(log_abf(v['beta_exp'], v['se_exp'], W) if v else 0.0)
        la2_b.append(r['logABF'])
    coloc['coloc_B_region_conservative'] = coloc_abf(la1_b, la2_b)
    print('  A (2 shared cis SNPs) PP.H4=%.4g  PP.H3=%.4g'
          % (coloc['coloc_A_symmetric_shared_cis']['PP.H4'],
             coloc['coloc_A_symmetric_shared_cis']['PP.H3']))
    print('  B (region, pQTL confined to leads) PP.H4=%.4g'
          % coloc['coloc_B_region_conservative']['PP.H4'])
    for w in (0.1 ** 2, 0.2 ** 2, 0.4 ** 2):
        for p12 in (1e-6, 1e-5, 1e-4):
            a1 = [log_abf(r['beta_exp'], r['se_exp'], w) for r in shared]
            a2 = [log_abf(r['beta_out'], r['se_out'], w) for r in shared]
            c = coloc_abf(a1, a2, p12=p12)
            coloc['sensitivity'].append({'W': w, 'p12': p12,
                                         'PP.H3': c['PP.H3'], 'PP.H4': c['PP.H4']})

    # 6) Open Targets 预计算共定位
    ot = json.load(open(COLOC_ALL, encoding='utf-8'))['data']['credibleSets']['rows']
    ot_sum = []
    for r in ot:
        c = r['colocalisation']
        ot_sum.append({'chrom': r['chromosome'], 'pos': r['position'],
                       'cis': not r['isTransQtl'], 'log10BF': r['credibleSetlog10BF'],
                       'n_coloc': c['count'],
                       'h4_max': max([x['h4'] for x in c['rows']], default=None),
                       'h4_min': min([x['h4'] for x in c['rows']], default=None)})
    print('\n[6] Open Targets precomputed coloc: cis loci counts =',
          [x['n_coloc'] for x in ot_sum if x['cis']])

    # 7) 效能
    r1 = json.load(open(ROUND1_JSON, encoding='utf-8'))
    se0 = r1['primary_wald_ratio_lead_cis']['se']
    power = {}
    for o in (1.1, 1.2, 1.5, 2.0, 3.0):
        b = math.log(o)
        power['OR_%.1f' % o] = float(norm.cdf(b / se0 - 1.96)
                                     + norm.cdf(-b / se0 - 1.96))
    power['min_detectable_OR_at_80pct'] = math.exp(2.8016 * se0)
    power['se_required_for_OR1.5_at_80pct'] = math.log(1.5) / 2.8016
    power['se_fold_reduction_needed'] = se0 / (math.log(1.5) / 2.8016)
    power['N_eff_fold_needed'] = power['se_fold_reduction_needed'] ** 2
    print('\n[7] Power with primary Wald se=%.4f: OR1.5 power=%.3f, '
          'MOR@80%%=%.1f, N_eff x%.0f needed'
          % (se0, power['OR_1.5'], power['min_detectable_OR_at_80pct'],
             power['N_eff_fold_needed']))

    result = {
        'generated': '2026-11-10',
        'data_sources': {
            'ld': 'Ensembl REST /ld/human/<rsid>/1000GENOMES:phase_3:EUR '
                  '(GRCh38); cross-chromosome pairs set to r2=0 structurally',
            'exposure': 'Open Targets Platform v26.9, UKB-PPP '
                        'UKB_PPP_EUR_CX3CL1_P78423_OID20976_v1, '
                        'nSamples=33436, 9 SuSiE credible sets',
            'outcome': 'FinnGen DF13/R13 DM_NEUROPATHY, '
                       'chr16:57,200,000-57,600,000, 3,140 variants',
            'outcome_cases': 'Risteys endpoint page reports 4,412 individuals '
                             'with DM_NEUROPATHY after genotype QC; '
                             'page does not state the release explicitly -> '
                             'used only for qualitative power discussion',
            'pqtL_full_sumstats': 'NOT publicly anonymous: OT registers '
                                  'gs://gentropy-vault/.../CX3CL1_P78423_'
                                  'OID20976_v1_Neurology.tar (403 anonymous); '
                                  'Synapse syn51364943 requires a free PAT',
            'ot_precomputed_coloc': 'Open Targets COLOC_PIP_ECAVIAR',
        },
        'ld_unknown_pairs': M['unknown_pairs'],
        'clumping': clumps,
        'mr': mr,
        'dpn_region_abf': region,
        'coloc': coloc,
        'ot_coloc': ot_sum,
        'power': power,
    }
    json.dump(result, open(OUT_JSON, 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    print('\nSaved:', OUT_JSON)


if __name__ == '__main__':
    main()
