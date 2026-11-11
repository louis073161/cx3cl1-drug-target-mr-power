# -*- coding: utf-8 -*-
"""
B3 步骤 2：真正的区域级共定位（region-wide ABF coloc）。

与已发表部分的口径差异（必须写进附录）：
  主分析：UKB-PPP 9 个 SuSiE 可信集变异 x FinnGen，coloc 为可信集级
          （Open Targets 预计算 + 本地 ABF），PP.H4=0.0584。
  本附录：Fenland (Pietzner 2021, GCST90012074, n=21,758 欧洲裔)
          x FinnGen，用区域全部可用变异做 ABF coloc，不再只用可信集变异。

为什么换成 Fenland 而不是 UKB-PPP 全量：
  见 a1_b3_fenland_fetch.py 头注：Synapse syn51364943 对本地区 403
  GEO_RESTRICTION，PAT 亦不可解；EBI FTP 无地理限制。
  副产物是方法学收益：暴露（英国 Fenland）与结局（芬兰 FinnGen）来自
  不同国家、不同队列，是 coloc 的理想独立样本结构。

统计：Wakefield ABF，单因果变异假设，先验 p1=p2=1e-4, p12=1e-5，
      与 a1_ld_coloc.py 的 coloc_abf() 完全同一实现（此处复制以免
      导入该脚本时触发其主流程）。
"""
import io
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from a1_multi_fetch import get_tbi, gcs_url          # noqa: E402
from a1_tabix_fetch import fetch_region               # noqa: E402

ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'results')
FEN_TSV = os.path.join(RES, 'a1_b3_fenland_region.tsv')

CHROM = '16'
START, END = 57200000, 57600000

# FinnGen 主结局 + 阴性对照端点
ENDPOINTS = [('R13', 'DM_NEUROPATHY'), ('R13', 'G6_POLYOTHUNS')]

# FinnGen sumstats 列序
COL = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'rsids': 4,
       'gene': 5, 'pval': 6, 'mlogp': 7, 'beta': 8, 'se': 9, 'af': 10}

W_MAIN = 0.2 ** 2
W_SENS = [0.1 ** 2, 0.2 ** 2, 0.4 ** 2]


# ---------------------------------------------------------------- 统计核心
# 与 a1_ld_coloc.py 逐字一致，保证与已发表部分同一口径
def log_abf(beta, se, W):
    V = se * se
    r = W / (V + W)
    z2 = (beta / se) ** 2
    return 0.5 * math.log(1.0 - r) + 0.5 * z2 * r


def logsumexp(xs):
    m = max(xs)
    if m == float('-inf'):
        return m
    return m + math.log(sum(math.exp(x - m) for x in xs))


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


# ---------------------------------------------------------------- 取数与对齐
def fetch_finngen(rel, ep):
    tbi_path, tbi = get_tbi(rel, ep)
    url = gcs_url(rel, ep)
    lines, _ = fetch_region(CHROM, START, END, url=url,
                            tbi_path=tbi_path, span_mb=6, tbi=tbi)
    rows = []
    for ln in lines:
        f = ln.split('\t')
        if len(f) < 11:
            continue
        try:
            pos = int(f[COL['pos']])
            beta = float(f[COL['beta']])
            se = float(f[COL['se']])
            pval = float(f[COL['pval']])
            af = f[COL['af']]
            af = float(af) if af not in ('NA', '', 'nan') else None
        except (ValueError, IndexError):
            continue
        if se <= 0 or not (0 < pval <= 1):
            continue
        if f[COL['chrom']].replace('chr', '') != CHROM:
            continue
        rows.append({'pos': pos, 'ref': f[COL['ref']].upper(),
                     'alt': f[COL['alt']].upper(), 'rsids': f[COL['rsids']],
                     'gene': f[COL['gene']], 'beta': beta, 'se': se,
                     'pval': pval, 'af': af})
    return rows


def load_fenland():
    rows = []
    with io.open(FEN_TSV, encoding='utf-8') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(head):
                continue
            try:
                beta = float(f[ix['beta']])
                se = float(f[ix['se']])
            except ValueError:
                continue
            if se <= 0:
                continue
            af = f[ix['eaf']]
            pv = f[ix['pval']]
            try:
                pv = float(pv) if pv not in ('NA', '', 'nan') else None
            except ValueError:
                pv = None
            rows.append({
                'pos': int(f[ix['pos']]),
                'effect_allele': f[ix['effect_allele']].upper(),
                'other_allele': f[ix['other_allele']].upper(),
                'beta': beta, 'se': se,
                'eaf': float(af) if af not in ('NA', '', 'nan') else None,
                'rsid': f[ix['rsid']], 'pval': pv,
            })
    return rows


COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def align(fen, fg):
    """把 Fenland 的 beta 对齐到 FinnGen 的 alt 等位基因方向。"""
    fe, fo = fen['effect_allele'], fen['other_allele']
    ga, gr = fg['alt'], fg['ref']
    if fe == ga and fo == gr:
        return fen['beta'], fen['eaf'], 'direct'
    if fe == gr and fo == ga:
        return -fen['beta'], (1 - fen['eaf']) if fen['eaf'] is not None else None, 'flipped'
    if COMP.get(fe) == ga and COMP.get(fo) == gr:
        return fen['beta'], fen['eaf'], 'complement'
    if COMP.get(fe) == gr and COMP.get(fo) == ga:
        return -fen['beta'], (1 - fen['eaf']) if fen['eaf'] is not None else None, 'complement_flipped'
    return None, None, None


def merge(fen_rows, fg_rows):
    by_pos_fg = {}
    for r in fg_rows:
        by_pos_fg.setdefault(r['pos'], []).append(r)
    shared = []
    stats = {'fen_in_region': len(fen_rows), 'fg_in_region': len(fg_rows),
             'pos_matched': 0, 'allele_resolved': 0, 'ambiguous_skipped': 0,
             'multi_allelic_positions': 0}
    for f in fen_rows:
        cands = by_pos_fg.get(f['pos'])
        if not cands:
            continue
        stats['pos_matched'] += 1
        if len(cands) > 1:
            stats['multi_allelic_positions'] += 1
        hit = None
        for g in cands:
            b, eaf, how = align(f, g)
            if b is not None:
                hit = (g, b, eaf, how)
                break
        if hit is None:
            stats['ambiguous_skipped'] += 1
            continue
        g, b, eaf, how = hit
        shared.append({
            'pos': f['pos'], 'rsid': f['rsid'],
            'ref': g['ref'], 'alt': g['alt'],
            'beta_exp': b, 'se_exp': f['se'], 'pval_exp': f['pval'],
            'eaf': eaf if eaf is not None else g['af'],
            'beta_out': g['beta'], 'se_out': g['se'], 'pval_out': g['pval'],
            'alignment': how,
        })
        stats['allele_resolved'] += 1
    stats['shared'] = len(shared)
    return shared, stats


def run_one(rel, ep, fen_rows):
    fg = fetch_finngen(rel, ep)
    shared, stats = merge(fen_rows, fg)
    out = {'release': rel, 'endpoint': ep,
           'n_fenland_region': stats['fen_in_region'],
           'n_finngen_region': stats['fg_in_region'],
           'n_shared': len(shared), 'merge_stats': stats}
    if len(shared) < 50:
        out['error'] = 'too few shared variants for stable ABF; reported as not evaluable'
        return out, []

    la1 = [log_abf(r['beta_exp'], r['se_exp'], W_MAIN) for r in shared]
    la2 = [log_abf(r['beta_out'], r['se_out'], W_MAIN) for r in shared]
    main = coloc_abf(la1, la2)
    out['coloc_W_0.04'] = main

    sens = []
    for W in W_SENS:
        a1 = [log_abf(r['beta_exp'], r['se_exp'], W) for r in shared]
        a2 = [log_abf(r['beta_out'], r['se_out'], W) for r in shared]
        s = coloc_abf(a1, a2)
        sens.append({'W': W, 'PP.H3': s['PP.H3'], 'PP.H4': s['PP.H4']})
    out['sensitivity_prior_W'] = sens

    # 先验敏感性
    pre = []
    for p12 in (1e-5, 5e-6, 2e-5):
        s = coloc_abf(la1, la2, p12=p12)
        pre.append({'p12': p12, 'PP.H4': s['PP.H4']})
    out['sensitivity_p12'] = pre

    # 两侧最强变异（用于报告，不用作结论）
    e_top = max(shared, key=lambda r: abs(r['beta_exp'] / r['se_exp']))
    o_top = max(shared, key=lambda r: abs(r['beta_out'] / r['se_out']))
    out['top_exposure_variant'] = {
        'rsid': e_top['rsid'], 'pos': e_top['pos'],
        'z': e_top['beta_exp'] / e_top['se_exp'], 'pval': e_top['pval_exp']}
    out['top_outcome_variant'] = {
        'rsid': o_top['rsid'], 'pos': o_top['pos'],
        'z': o_top['beta_out'] / o_top['se_out'], 'pval': o_top['pval_out']}
    out['n_outcome_p_below_1e5'] = sum(
        1 for r in shared if (r['pval_out'] is not None and r['pval_out'] < 1e-5))
    out['n_exposure_p_below_5e8'] = sum(
        1 for r in shared if (r['pval_exp'] is not None and r['pval_exp'] < 5e-8))
    return out, shared


def main():
    if not os.path.exists(FEN_TSV):
        print('缺少 %s，请先运行 a1_b3_fenland_fetch.py' % FEN_TSV)
        return 1
    fen_rows = load_fenland()
    print('Fenland 区域变异 %d' % len(fen_rows), flush=True)

    results = []
    all_shared = {}
    for rel, ep in ENDPOINTS:
        print('\n[%s %s] 取 FinnGen 区域 ...' % (rel, ep), flush=True)
        try:
            r, shared = run_one(rel, ep, fen_rows)
        except Exception as e:
            print('  失败 %s: %s' % (type(e).__name__, str(e)[:200]), flush=True)
            results.append({'release': rel, 'endpoint': ep,
                            'error': '%s: %s' % (type(e).__name__, str(e)[:200])})
            continue
        results.append(r)
        all_shared[(rel, ep)] = shared
        if 'coloc_W_0.04' in r:
            c = r['coloc_W_0.04']
            print('  共享变异 %d | PP.H3=%.4f PP.H4=%.4f | 结局最显著 P=%.3g'
                  % (r['n_shared'], c['PP.H3'], c['PP.H4'],
                     min(x['pval_out'] for x in shared)), flush=True)

    # 落盘
    out_json = {'method': 'Wakefield ABF region-wide coloc, single causal '
                          'variant, priors p1=p2=1e-4, p12=1e-5',
                'exposure': 'Fenland CX3CL1 pQTL (GCST90012074, n=21,758 '
                            'European, GRCh38)',
                'region': 'chr16:%d-%d' % (START, END),
                'W_main': W_MAIN,
                'results': results,
                'note': ('This is a true region-wide coloc using all variants '
                         'shared between the two studies, replacing the '
                         'credible-set-level approximation used in the main '
                         'analysis. Exposure and outcome are independent '
                         'cohorts from different countries (UK, Finland).')}
    with io.open(os.path.join(RES, 'a1_b3_coloc_full.json'), 'w',
                 encoding='utf-8') as fh:
        json.dump(out_json, fh, ensure_ascii=False, indent=1)

    # 变异级明细（取主结局）
    key = ('R13', 'DM_NEUROPATHY')
    if key in all_shared and all_shared[key]:
        cols = ['pos', 'rsid', 'ref', 'alt', 'beta_exp', 'se_exp', 'pval_exp',
                'eaf', 'beta_out', 'se_out', 'pval_out', 'alignment']
        with io.open(os.path.join(RES, 'a1_b3_coloc_variants.tsv'), 'w',
                     encoding='utf-8', newline='') as fh:
            fh.write('\t'.join(cols) + '\n')
            for r in sorted(all_shared[key], key=lambda x: x['pos']):
                fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
        print('\n写入 a1_b3_coloc_variants.tsv（%d 行）'
              % len(all_shared[key]), flush=True)

    print('\n写入 a1_b3_coloc_full.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
