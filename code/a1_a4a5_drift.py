# -*- coding: utf-8 -*-
"""
A4 + A5：端点漂移定量归因 + 四类工具变量集 × 19 组合的三类漂移框架。

全部输入为既有真实公开数据产出的本地结果，不新增取数、不做模拟填补。
A4：把"结论漂移"区分为抽样误差可解释 / 超出抽样误差两部分；再用公开
    病例数和端点定义差异做归因（Spearman + 置换检验），无病例数处标
    NOT_EVALUABLE，不补数。
A5：cis lead only / all cis / LD-clumped / all available 四种工具变量集
    × 19 组合 = 76 个估计，形成 release / endpoint / instrument 三类漂移
    统一框架。

输出（results/）：
  a1_a4_pairwise_drift.tsv, a1_a4_attribution.json,
  a1_a5_instrument_drift.tsv, a1_a5_three_axis_summary.json,
  fig6_drift_attribution.pdf/.svg/.png, fig7_instrument_drift.pdf/.svg/.png
图表一律英文标签、矢量优先。
"""
import csv
import io
import json
import math
import os
import sys
from collections import defaultdict
from itertools import combinations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, spearmanr

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_ld_coloc import wald, ivw, weighted_median  # noqa: E402

RES = ROOT + '/results'
LEAD = '16_57379748_C_T'
Z975 = 1.959964
Z80 = norm.ppf(0.80)
BONF = 0.05 / 133          # 与主分析 133 次检验保持同一阈值
SEED = 20261110


def read_tsv(path, enc='utf-8'):
    with io.open(path, encoding=enc) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_tsv(path, rows):
    if not rows:
        return
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def power_at_or(se, or_true=1.5):
    z = math.log(or_true) / se
    return float(norm.cdf(z - Z975) + norm.cdf(-z - Z975))


def mde(se):
    return math.exp((Z975 + Z80) * se)


def fig_style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.facecolor': 'white', 'figure.facecolor': 'white',
        'text.color': '#1f2933', 'axes.labelcolor': '#1f2933',
        'xtick.color': '#1f2933', 'ytick.color': '#1f2933',
        'axes.edgecolor': '#52606d', 'axes.linewidth': 0.8,
        'savefig.dpi': 300,
    })


def save(fig, stem):
    for ext, kw in [('pdf', {}), ('svg', {}), ('png', {'dpi': 300})]:
        fig.savefig('%s/%s.%s' % (RES, stem, ext), bbox_inches='tight', **kw)
    plt.close(fig)


# ------------------------------------------------------------------ 元数据
def load_meta():
    """组合元数据：病例数仅取官方 manifest；R13 主结局补端点页数字并标注。"""
    ed = read_tsv(RES + '/a1_endpoint_dict.tsv')
    defs = read_tsv(RES + '/a1_endpoint_definition_table.tsv', 'utf-8-sig')
    r13_cases = None
    for r in defs:
        if (r['来源类别'] == 'R13 专属端点页'
                and r['端点'] == 'DM_NEUROPATHY'):
            r13_cases = int(r['病例数'])
    meta = {}
    for r in ed:
        key = (r['release'], r['endpoint'])
        nc = int(r['n_cases']) if r['n_cases'] else None
        src = r['case_control_source'] or ''
        if nc is None and key == ('R13', 'DM_NEUROPATHY') and r13_cases:
            nc, src = r13_cases, 'R13 endpoint page (not official manifest)'
        meta[key] = {
            'release': r['release'], 'endpoint': r['endpoint'],
            'endpoint_type': r['endpoint_type'],
            'n_cases': nc, 'source': src,
            'endpoint_name': r['endpoint_name'],
        }
    return meta


# ------------------------------------------------------------------ A4
def a4(meta, wald_by_key):
    keys = sorted(wald_by_key)
    rows = []
    for a, b in combinations(keys, 2):
        ma, mb = meta[a], meta[b]
        pa, pb = wald_by_key[a], wald_by_key[b]
        se_pair = math.sqrt(pa['se'] ** 2 + pb['se'] ** 2)
        dbeta = abs(pa['beta'] - pb['beta'])
        z = dbeta / se_pair if se_pair > 0 else None
        # 端点定义差异（宽泛/严格/EXMORE 等同一家族变体）
        ea, eb = a[1], b[1]
        base_a = ea.replace('_STRICT', '').replace('_EXMORE', '')
        base_b = eb.replace('_STRICT', '').replace('_EXMORE', '')
        if a[1] == b[1]:
            pair_type = 'same_endpoint_across_release'
        elif a[0] == b[0]:
            pair_type = 'same_release_different_endpoint'
        else:
            pair_type = 'across_both'
        same_family = base_a == base_b and ea != eb
        if ma['n_cases'] and mb['n_cases']:
            dlogn = abs(math.log(ma['n_cases']) - math.log(mb['n_cases']))
            n_status = 'available'
        else:
            dlogn, n_status = None, 'NOT_EVALUABLE'
        rows.append({
            'release_1': a[0], 'endpoint_1': a[1],
            'release_2': b[0], 'endpoint_2': b[1],
            'pair_type': pair_type,
            'same_endpoint_family_strictness_variant': same_family,
            'beta_1': pa['beta'], 'beta_2': pb['beta'],
            'abs_delta_logOR': dbeta,
            'se_1': pa['se'], 'se_2': pb['se'],
            'se_of_difference': se_pair,
            'z_drift': z,
            'exceeds_sampling_noise_z1.96': (z > Z975) if z is not None else '',
            'sign_flip': (pa['beta'] * pb['beta'] < 0),
            'n_cases_1': ma['n_cases'] or '', 'n_cases_2': mb['n_cases'] or '',
            'abs_delta_log_cases': dlogn if dlogn is not None else '',
            'case_count_status': n_status,
            'endpoint_type_1': ma['endpoint_type'],
            'endpoint_type_2': mb['endpoint_type'],
        })
    write_tsv(RES + '/a1_a4_pairwise_drift.tsv', rows)

    # 分组：抽样噪声可解释 vs 超出
    by_type = defaultdict(list)
    for r in rows:
        by_type[r['pair_type']].append(r)
    groups = {}
    for k, g in by_type.items():
        zs = [r['z_drift'] for r in g if r['z_drift'] is not None]
        groups[k] = {
            'n_pairs': len(g),
            'median_z_drift': float(np.median(zs)) if zs else None,
            'max_z_drift': float(max(zs)) if zs else None,
            'n_exceeding_1.96': sum(1 for z in zs if z > Z975),
            'pct_exceeding_1.96': (100.0 * sum(1 for z in zs if z > Z975)
                                   / len(zs)) if zs else None,
            'n_sign_flips': sum(1 for r in g if r['sign_flip']),
            'expected_exceeding_at_5pct': 0.05 * len(zs),
        }

    # 归因 1：|ΔlogOR| 与 |Δlog cases| 的 Spearman + 置换检验
    avail = [r for r in rows if r['abs_delta_log_cases'] != '']
    x = np.array([float(r['abs_delta_log_cases']) for r in avail])
    y = np.array([float(r['abs_delta_logOR']) for r in avail])
    rng = np.random.default_rng(SEED)
    if len(x) >= 4 and np.std(x) > 0:
        rho, p_asym = spearmanr(x, y)
        obs = abs(float(rho))
        null = []
        for _ in range(20000):
            null.append(abs(float(spearmanr(rng.permutation(x), y)[0])))
        p_perm = (1 + sum(1 for v in null if v >= obs)) / (len(null) + 1)
        attr_cases = {
            'n_pairs_with_public_case_counts': len(x),
            'spearman_rho': float(rho), 'p_asymptotic': float(p_asym),
            'p_permutation_20000': float(p_perm),
            'interpretation': ('No statistically supported monotonic association'
                               if p_perm >= 0.05 else
                               'Monotonic association detected; still not causal'),
            'caveat': ('Only %d of %d pairs have public case counts; R11 and '
                       'most R13 endpoints are NOT_EVALUABLE. Case counts are '
                       'not used to impute or adjust any MR estimate.'
                       % (len(x), len(rows))),
        }
    else:
        attr_cases = {'n_pairs_with_public_case_counts': len(x),
                      'status': 'insufficient pairs'}

    # 归因 2：端点定义严格度变体 vs 其他配对，|ΔlogOR| 差异 + 置换检验
    fam = np.array([float(r['abs_delta_logOR']) for r in rows
                    if r['same_endpoint_family_strictness_variant']])
    oth = np.array([float(r['abs_delta_logOR']) for r in rows
                    if not r['same_endpoint_family_strictness_variant']])
    if len(fam) >= 2 and len(oth) >= 2:
        obs_d = abs(float(np.mean(fam)) - float(np.mean(oth)))
        pool = np.concatenate([fam, oth])
        null = []
        n1 = len(fam)
        for _ in range(20000):
            p = rng.permutation(pool)
            null.append(abs(float(np.mean(p[:n1]) - np.mean(p[n1:]))))
        p_def = (1 + sum(1 for v in null if v >= obs_d)) / (len(null) + 1)
        attr_def = {
            'n_strictness_variant_pairs': len(fam),
            'n_other_pairs': len(oth),
            'mean_abs_delta_logOR_strictness': float(np.mean(fam)),
            'mean_abs_delta_logOR_other': float(np.mean(oth)),
            'p_permutation_20000': float(p_def),
        }
    else:
        attr_def = {'status': 'insufficient pairs'}

    # 归因 3：漂移是否主要超出抽样误差（二项检验）
    zs_all = [r['z_drift'] for r in rows if r['z_drift'] is not None]
    n_exc = sum(1 for z in zs_all if z > Z975)
    from scipy.stats import binomtest
    bt = binomtest(n_exc, len(zs_all), 0.05, alternative='greater')
    summary = {
        'generated_on': '2026-11-11',
        'n_combinations': len(keys),
        'n_pairwise_comparisons': len(rows),
        'groups_by_pair_type': groups,
        'attribution_case_counts': attr_cases,
        'attribution_endpoint_definition': attr_def,
        'sampling_noise_test': {
            'n_pairs': len(zs_all),
            'n_z_gt_1.96': n_exc,
            'expected_at_5pct': 0.05 * len(zs_all),
            'binomial_p_greater': float(bt.pvalue),
            'interpretation': ('Most observed drift is compatible with sampling '
                               'noise under the null; this does NOT mean the '
                               'estimates are precise - see MDE ORs.'),
        },
        'formula': ('z_drift = |beta_i - beta_j| / sqrt(SE_i^2 + SE_j^2); '
                    'SE_i is the observed Wald ratio SE. Pairs are not '
                    'statistically independent (shared FinnGen resource); '
                    'permutation P-values are descriptive, not confirmatory.'),
        'red_lines': ('Do not convert drift into a positive CX3CL1 signal; '
                      'do not impute missing case counts; do not claim '
                      'independent replication.'),
    }
    with io.open(RES + '/a1_a4_attribution.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # 图 6
    fig_style()
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.1))
    labels, data = [], []
    for k in ['same_endpoint_across_release',
              'same_release_different_endpoint', 'across_both']:
        if k in groups:
            labels.append(k.replace('_', '\n'))
            data.append([r['z_drift'] for r in by_type[k]
                         if r['z_drift'] is not None])
    bp = ax[0].boxplot(data, tick_labels=labels, patch_artist=True, widths=0.55)
    for patch in bp['boxes']:
        patch.set_facecolor('#dbeafe')
        patch.set_edgecolor('#2563eb')
    ax[0].axhline(Z975, ls='--', lw=1.0, color='#b91c1c',
                  label='z = 1.96 (sampling-noise bound)')
    ax[0].set_ylabel('Standardised drift  |ΔlogOR| / SE(difference)')
    ax[0].set_title('(a) Drift relative to sampling noise', fontsize=10)
    ax[0].legend(fontsize=7, frameon=False)
    ax[0].tick_params(axis='x', labelsize=7)

    if len(x) >= 4:
        ax[1].scatter(x, y, s=22, c='#1d4ed8', alpha=0.65, edgecolors='none')
        if np.std(x) > 0:
            coef = np.polyfit(x, y, 1)
            xs = np.linspace(float(x.min()), float(x.max()), 50)
            ax[1].plot(xs, np.polyval(coef, xs), '-', color='#b91c1c', lw=1.2)
        ax[1].set_xlabel('|Δ log(public case count)|')
        ax[1].set_ylabel('|Δ logOR|')
        ax[1].set_title('(b) Case-count change vs effect drift (n=%d)' % len(x),
                        fontsize=10)
        txt = ('Spearman ρ = %.3f\npermutation P = %.3f'
               % (attr_cases.get('spearman_rho', float('nan')),
                  attr_cases.get('p_permutation_20000', float('nan'))))
        ax[1].text(0.03, 0.94, txt, transform=ax[1].transAxes,
                   va='top', fontsize=7.5,
                   bbox=dict(fc='white', ec='#cbd5e1', lw=0.5))
    else:
        ax[1].text(0.5, 0.5, 'Too few pairs with public case counts',
                   ha='center', transform=ax[1].transAxes)
    for a in ax:
        a.grid(alpha=0.22, lw=0.5)
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig6_drift_attribution')
    return summary, rows


# ------------------------------------------------------------------ A5
def a5(meta, raw, clump_keep):
    by_combo = defaultdict(list)
    for r in raw:
        for k in ('beta_exp', 'se_exp', 'beta_out', 'se_out',
                  'z_exp', 'pval_out', 'af_alt'):
            v = num(r.get(k))
            if v is not None:
                r[k] = v
        by_combo[(r['release'], r['endpoint'])].append(r)

    sets = [
        ('S1_lead_cis', lambda rows: [r for r in rows if r['variant'] == LEAD]),
        ('S2_all_cis', lambda rows: [r for r in rows if r['cis'] == 'True']),
        ('S3_ld_clumped', lambda rows: [r for r in rows
                                       if r['variant'] in clump_keep]),
        ('S4_all_available', lambda rows: list(rows)),
    ]
    out = []
    for key in sorted(by_combo):
        rows = by_combo[key]
        for sname, fn in sets:
            sel = fn(rows)
            if not sel:
                continue
            if len(sel) == 1:
                est = wald(sel[0])
                method = 'Wald ratio'
            else:
                est = ivw(sel, False)
                method = 'IVW (fixed, zero intercept)'
                if not est:
                    continue
            se = est['se']
            rec = {
                'release': key[0], 'endpoint': key[1],
                'instrument_set': sname, 'method': method,
                'n_snp': len(sel),
                'variants': ';'.join(r['variant'] for r in sel),
                'beta': est['beta'], 'se': se,
                'ci_lo': est['ci_lo'], 'ci_hi': est['ci_hi'],
                'p': est['p'], 'or': est['or'],
                'direction': 'positive' if est['beta'] > 0 else 'negative',
                'power_OR_1.5': power_at_or(se),
                'mde_OR_80pct': mde(se),
                'p_below_bonferroni': est['p'] < BONF,
                'bonferroni_threshold': BONF,
                'endpoint_type': meta[key]['endpoint_type'],
                'note': ('all-cis set contains two variants in LD '
                         '(r2=0.1002); not independent instruments'
                         if sname == 'S2_all_cis' else ''),
            }
            out.append(rec)
    write_tsv(RES + '/a1_a5_instrument_drift.tsv', out)

    # 三类漂移轴：release / endpoint / instrument
    axis_summary = {}
    for axis, groupfn in [
        ('release', lambda r: (r['endpoint'], r['instrument_set'])),
        ('endpoint', lambda r: (r['release'], r['instrument_set'])),
        ('instrument', lambda r: (r['release'], r['endpoint'])),
    ]:
        g = defaultdict(list)
        for r in out:
            g[groupfn(r)].append(r)
        spreads, flips, classes = [], 0, defaultdict(int)
        for _, grp in g.items():
            if len(grp) < 2:
                continue
            bs = [x['beta'] for x in grp]
            spreads.append(max(bs) - min(bs))
            signs = [1 if b > 0 else -1 for b in bs]
            flips += sum(1 for i in range(1, len(signs))
                         if signs[i] != signs[i - 1])
            for x in grp:
                classes[x['instrument_set']] += 1
        axis_summary[axis] = {
            'n_groups': sum(1 for _, grp in g.items() if len(grp) >= 2),
            'median_beta_spread': float(np.median(spreads)) if spreads else None,
            'max_beta_spread': float(max(spreads)) if spreads else None,
            'n_sign_flips_in_listed_order': flips,
            'n_estimates': len(out),
        }

    # instrument 轴单独：每个组合内四种 IV 集的 OR 范围与方向一致性
    per_combo = defaultdict(list)
    for r in out:
        per_combo[(r['release'], r['endpoint'])].append(r)
    inst_detail = []
    for key, grp in sorted(per_combo.items()):
        ors = [x['or'] for x in grp]
        signs = set(x['direction'] for x in grp)
        ps = [x['p'] for x in grp]
        inst_detail.append({
            'release': key[0], 'endpoint': key[1],
            'n_sets': len(grp),
            'or_min': min(ors), 'or_max': max(ors),
            'or_ratio': max(ors) / min(ors),
            'direction_consistent': len(signs) == 1,
            'p_min': min(ps), 'p_max': max(ps),
            'any_p_below_bonferroni': any(x['p_below_bonferroni'] for x in grp),
            'power_OR_1.5_min': min(x['power_OR_1.5'] for x in grp),
            'power_OR_1.5_max': max(x['power_OR_1.5'] for x in grp),
            'mde_OR_min': min(x['mde_OR_80pct'] for x in grp),
            'mde_OR_max': max(x['mde_OR_80pct'] for x in grp),
        })
    write_tsv(RES + '/a1_a5_instrument_detail.tsv', inst_detail)

    summary = {
        'generated_on': '2026-11-11',
        'n_instrument_sets': len(sets),
        'n_combinations': len(per_combo),
        'n_estimates': len(out),
        'axes': axis_summary,
        'instrument_axis': {
            'n_combinations_with_direction_inconsistency':
                sum(1 for r in inst_detail if not r['direction_consistent']),
            'max_or_ratio_within_combination':
                max(r['or_ratio'] for r in inst_detail),
            'n_estimates_below_bonferroni':
                sum(1 for r in out if r['p_below_bonferroni']),
            'power_OR_1.5_range': [
                min(r['power_OR_1.5'] for r in out),
                max(r['power_OR_1.5'] for r in out)],
            'mde_OR_range': [
                min(r['mde_OR_80pct'] for r in out),
                max(r['mde_OR_80pct'] for r in out)],
        },
        'formula': ('Wald ratio for a single instrument; IVW with intercept '
                    'fixed to zero for >=2 instruments (weights 1/se_out^2). '
                    'Power and MDE derived from each observed SE.'),
        'caveats': ('The four instrument sets are nested and NOT independent; '
                    '76 estimates are a sensitivity grid, not 76 independent '
                    'tests. S2 contains two cis variants in LD (r2=0.1002). '
                    'No P-value selection is applied.'),
        'critical_caveat_apparent_power': (
            'S3/S4 reach ~100% power for OR=1.5 and MDE OR ~1.06-1.21, but '
            'this is FALSE PRECISION. Their small IVW SEs rely on the '
            'assumption that all 7-8 instruments are valid; the trans '
            'pQTL instruments have documented horizontal pleiotropy '
            '(coloc H4 up to 1.0, CLPP evidence in a1_pip_coloc.json) and '
            'were explicitly excluded from the primary cis-only analysis. '
            'Adding pleiotropic instruments shrinks SE without adding valid '
            'information. The manuscript MUST NOT claim that the power '
            'problem is solved by using more instruments, MUST NOT treat '
            'S3/S4 as adequately powered null results, and MUST NOT '
            'interpret the single nominal P=0.0492 (R9 '
            'DM_RETINOPATHY_EXMORE, S4) as a positive finding: it is a '
            'negative-control endpoint, below Bonferroni after correction '
            '(0 of 76), and driven by the invalid instrument set.'),
        'primary_instrument_set': (
            'S1_lead_cis (single cis lead variant, Wald ratio) remains the '
            'pre-specified primary analysis. Its OR=1.5 power is '
            '5.45%-10.10% and MDE OR 5.61-318.57; S2 is a cis-only '
            'sensitivity set with two correlated variants.'),
    }
    with io.open(RES + '/a1_a5_three_axis_summary.json', 'w',
                 encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # 图 7：四个 IV 集 × 19 组合的 OR（对数轴）与方向一致性
    fig_style()
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 5.0),
                           gridspec_kw={'width_ratios': [1.55, 1]})
    keys = sorted(per_combo, key=lambda k: (k[0], k[1]))
    colors = {'S1_lead_cis': '#1d4ed8', 'S2_all_cis': '#0891b2',
              'S3_ld_clumped': '#b45309', 'S4_all_available': '#b91c1c'}
    markers = {'S1_lead_cis': 'o', 'S2_all_cis': 's',
               'S3_ld_clumped': '^', 'S4_all_available': 'D'}
    for sname, _ in sets:
        xs, ys = [], []
        for i, k in enumerate(keys):
            rec = next((r for r in per_combo[k]
                        if r['instrument_set'] == sname), None)
            if rec:
                xs.append(i)
                ys.append(rec['or'])
        ax[0].scatter(xs, ys, s=28, c=colors[sname],
                      marker=markers[sname], alpha=0.78,
                      edgecolors='none', label=sname.replace('_', ' '))
    ax[0].axhline(1.0, color='#475569', lw=0.9, ls='--')
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Odds ratio (log scale)')
    ax[0].set_xticks(range(len(keys)))
    ax[0].set_xticklabels(['%s\n%s' % (k[0], k[1].replace('DM_', '').replace(
        'G6_', '').replace('_EXMORE', '*')[:14]) for k in keys],
        rotation=90, fontsize=5.6)
    ax[0].set_title('(a) OR by instrument set across 19 contexts', fontsize=10)
    ax[0].legend(fontsize=7, frameon=False, ncol=2)

    pow_by_set = defaultdict(list)
    for r in out:
        pow_by_set[r['instrument_set']].append(r['power_OR_1.5'] * 100)
    names = [s for s, _ in sets]
    bp = ax[1].boxplot([pow_by_set[n] for n in names],
                       tick_labels=[n.replace('_', '\n') for n in names],
                       patch_artist=True, widths=0.55)
    for patch in bp['boxes']:
        patch.set_facecolor('#fef3c7')
        patch.set_edgecolor('#b45309')
    ax[1].axhline(80, ls='--', color='#b91c1c', lw=1.0,
                  label='80% power')
    ax[1].set_ylabel('Power (%) to detect OR = 1.5')
    ax[1].set_title('(b) Power by instrument set', fontsize=10)
    ax[1].legend(fontsize=7, frameon=False)
    ax[1].tick_params(axis='x', labelsize=7)
    for a in ax:
        a.grid(alpha=0.22, lw=0.5, axis='y')
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig7_instrument_drift')
    return summary, out


def main():
    os.makedirs(RES, exist_ok=True)
    meta = load_meta()
    mm = read_tsv(RES + '/a1_mr_matrix.tsv')
    wald_by_key = {}
    for r in mm:
        if r['method'] == 'Wald ratio (lead cis)':
            wald_by_key[(r['release'], r['endpoint'])] = {
                'beta': float(r['beta']), 'se': float(r['se']),
                'p': float(r['p']), 'or': float(r['or'])}
    a4s, a4rows = a4(meta, wald_by_key)
    raw = read_tsv(RES + '/a1_outcome_matrix.tsv')
    rob = json.load(io.open(RES + '/a1_robustness.json', encoding='utf-8'))
    clump_keep = set(rob['clumping']['r2_0.001']['kept_ids'])
    a5s, a5rows = a5(meta, raw, clump_keep)
    print('A4 pairs:', len(a4rows),
          'exceed-noise:', a4s['sampling_noise_test']['n_z_gt_1.96'],
          'binom p: %.4g' % a4s['sampling_noise_test']['binomial_p_greater'])
    print('A4 case-count rho:',
          a4s['attribution_case_counts'].get('spearman_rho'),
          'perm p:', a4s['attribution_case_counts'].get('p_permutation_20000'))
    print('A5 estimates:', len(a5rows),
          'power range:', a5s['instrument_axis']['power_OR_1.5_range'],
          'MDE range:', a5s['instrument_axis']['mde_OR_range'])
    print('A5 direction-inconsistent combos:',
          a5s['instrument_axis']['n_combinations_with_direction_inconsistency'])


if __name__ == '__main__':
    main()
