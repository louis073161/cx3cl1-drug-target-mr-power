# -*- coding: utf-8 -*-
"""
把 A1 的图 1—5 重绘为英文标签、矢量优先（PDF/SVG）+ PNG 的投稿版。
数据与既有图完全一致，仅替换文字语言、字体与输出格式；不改动任何数值。

输出前缀 *_EN，原中文 PNG 保留作为内部留档，不覆盖。
"""
import json
import io
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.stats import norm, binom

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_ld_coloc import load_fg_region  # noqa: E402

RES = ROOT + '/results'
Z975, Z80 = norm.ppf(0.975), norm.ppf(0.80)
COLORS = {'primary': '#c0392b', 'proxy_painful_DPN': '#8e44ad',
          'broad_polyneuropathy': '#2471a3',
          'negative_ctrl_microvascular': '#7f8c8d'}
LABELS = {'primary': 'Primary: diabetic neuropathy',
          'proxy_painful_DPN': 'Painful-neuropathy proxy',
          'broad_polyneuropathy': 'Broad polyneuropathy',
          'negative_ctrl_microvascular': 'Negative control: microvascular'}


def style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.facecolor': 'white', 'figure.facecolor': 'white',
        'text.color': '#1f2933', 'axes.labelcolor': '#1f2933',
        'xtick.color': '#1f2933', 'ytick.color': '#1f2933',
        'axes.edgecolor': '#52606d', 'axes.linewidth': 0.8,
    })


def save(fig, stem):
    for ext, kw in [('pdf', {}), ('svg', {}), ('png', {'dpi': 300})]:
        fig.savefig('%s/%s_EN.%s' % (RES, stem, ext),
                    bbox_inches='tight', **kw)
    plt.close(fig)


def fig1():
    fg = load_fg_region()
    cs = json.load(io.open(ROOT + '/data/ukbppp_cx3cl1_crediblesets.json',
                           encoding='utf-8'))['data']['credibleSets']['rows']
    cis = [r for r in cs if not r['isTransQtl']]
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True,
                           gridspec_kw={'height_ratios': [1.6, 1]})
    x = [r['pos'] / 1e6 for r in fg]
    y = [-math.log10(r['pval']) for r in fg]
    ax[0].scatter(x, y, s=6, c='#8c9bab', alpha=0.75, linewidths=0)
    for r in cis:
        ax[0].axvline(r['position'] / 1e6, color='#c0392b', lw=1.2,
                      ls='--', alpha=0.85)
    ax[0].axhline(5, color='#e67e22', lw=0.9, ls=':')
    ax[0].text(57.205, 5.12, 'P = 1e-5 (no variant reaches it)',
               fontsize=8, color='#e67e22')
    best = min(fg, key=lambda r: r['pval'])
    ax[0].annotate('%s  P=%.2g' % (best['rsids'], best['pval']),
                   xy=(best['pos'] / 1e6, -math.log10(best['pval'])),
                   xytext=(57.30, 3.3), fontsize=8,
                   arrowprops=dict(arrowstyle='->', lw=0.7, color='#555'))
    ax[0].set_ylabel('-log10(P)  FinnGen R13 diabetic neuropathy')
    ax[0].set_title('CX3CL1 cis region (chr16, GRCh38): no outcome '
                    'association signal at the pQTL locus', fontsize=10)
    ax[0].set_ylim(0, 5.6)
    ax[0].grid(axis='y', alpha=0.25, lw=0.5)
    for r in cis:
        p = r['pValueMantissa'] * (10 ** r['pValueExponent'])
        lp = -math.log10(p) if p > 0 else 0
        ax[1].vlines(r['position'] / 1e6, 0, lp, color='#c0392b', lw=3)
        label = r['variant']['rsIds'][0] if r['variant']['rsIds'] else ''
        ax[1].text(r['position'] / 1e6 + 0.004, lp + 6,
                   '%s  (P=%.1g)' % (label, p),
                   fontsize=7, ha='left', va='bottom', color='#c0392b')
    ax[1].set_ylabel('-log10(P)  UKB-PPP CX3CL1 pQTL')
    ax[1].set_xlabel('chr16 position (Mb, GRCh38)')
    ax[1].set_ylim(0, 215)
    ax[1].grid(axis='y', alpha=0.25, lw=0.5)
    for a in ax:
        a.set_xlim(57.2, 57.6)
        a.tick_params(labelsize=8)
    fig.tight_layout()
    save(fig, 'fig1_cx3cl1_region')


def fig2():
    m = json.load(io.open(RES + '/a1_multiendpoint_mr.json', encoding='utf-8'))
    rows = []
    for c in m['combos']:
        p = c['methods'].get('Wald ratio (lead cis)')
        if p:
            rows.append((c['release'], c['endpoint'], c['endpoint_type'], p))
    rows.sort(key=lambda x: (int(x[0][1:]), x[1]))
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(11, 8.8))
    y = list(range(len(rows)))
    for i, (rel, ep, et, p) in enumerate(rows):
        col = COLORS.get(et, '#555555')
        lo, hi = math.exp(p['ci_lo']), math.exp(p['ci_hi'])
        ax.plot([lo, hi], [i, i], color=col, lw=1.8, alpha=0.85, zorder=2)
        ax.scatter([p['or']], [i], s=42, color=col, zorder=3,
                   edgecolors='white', linewidths=0.7)
        ax.text(58, i, '%.2f (%.2f-%.2f)' % (p['or'], lo, hi),
                va='center', ha='left', fontsize=7.6, color='#34495e')
    ax.axvline(1, color='#2c3e50', lw=1.1, ls='--', zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(['%s  %s' % (r[0], r[1]) for r in rows], fontsize=8.2)
    ax.set_xscale('log')
    ax.set_xlim(0.004, 55)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel('Odds ratio (log scale; Wald ratio, single cis lead instrument)')
    ax.set_title('One instrument x 19 release-by-endpoint contexts: '
                 'point estimates and 95% confidence intervals',
                 fontsize=11, pad=12)
    ax.grid(axis='x', alpha=0.22, lw=0.5)
    handles = [plt.Line2D([0], [0], color=v, lw=3, label=LABELS[k])
               for k, v in COLORS.items()]
    ax.legend(handles=handles, loc='upper center',
              bbox_to_anchor=(0.42, -0.105), ncol=4, fontsize=8.0,
              frameon=True, facecolor='white', edgecolor='#cccccc')
    fig.text(0.012, 0.012,
             'Dashed line: OR=1. All 19 primary P values >0.05. '
             'Right column: OR (95% CI).', fontsize=8, color='#555555')
    fig.tight_layout(rect=(0, 0.075, 0.80, 1))
    save(fig, 'fig2_forest_multiendpoint')


def fig3():
    m = json.load(io.open(RES + '/a1_multiendpoint_mr.json', encoding='utf-8'))
    picks = [
        ('R13', 'DM_NEUROPATHY', '#c0392b', 'Primary R13 diabetic neuropathy'),
        ('R11', 'DM_NEUROPATHY', '#e67e22', 'Primary R11 diabetic neuropathy'),
        ('R13', 'G6_POLYOTHUNS', '#2471a3', 'Broad polyneuropathy R13'),
        ('R13', 'DM_RETINOPATHY', '#7f8c8d', 'Highest power: retinopathy R13'),
        ('R9', 'DM_POLYNEURO', '#8e44ad', 'Lowest power: proxy R9'),
    ]
    ses = {}
    for c in m['combos']:
        if c['power']:
            ses[(c['release'], c['endpoint'])] = c['power']['se']
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    xs = [1.0 + 0.02 * i for i in range(126)]
    for rel, ep, col, lab in picks:
        se = ses.get((rel, ep))
        if not se:
            continue
        z = [math.log(x) / se for x in xs]
        pw = [norm.cdf(v - Z975) + norm.cdf(-v - Z975) for v in z]
        ax.plot(xs, pw, color=col, lw=2, label=lab)
        p15 = (norm.cdf(math.log(1.5) / se - Z975)
               + norm.cdf(-math.log(1.5) / se - Z975))
        ax.scatter([1.5], [p15], s=32, color=col, zorder=4)
    ax.axhline(0.80, color='#2c3e50', lw=1.2, ls='--')
    ax.text(3.48, 0.815, '80% power', ha='right', fontsize=8.5,
            color='#2c3e50')
    ax.axvline(1.5, color='#95a5a6', lw=0.9, ls=':')
    ax.text(1.51, 0.015, 'OR=1.5', fontsize=8.2, color='#5d6d7e')
    ax.set_xlim(1.0, 3.5)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel('Assumed true odds ratio')
    ax.set_ylabel('Power at two-sided alpha = 0.05')
    ax.set_title('Power curves driven by each observed Wald-ratio SE '
                 '(no case-count assumptions)', fontsize=11, pad=10)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc='upper left', fontsize=8.4, frameon=True,
              facecolor='white', edgecolor='#cccccc')
    ax.text(2.05, 0.035,
            'Even the highest-power public context reaches only ~10% '
            'power for OR=1.5', fontsize=8.0, color='#34495e')
    fig.tight_layout()
    save(fig, 'fig3_power_curves')


def fig4():
    fig, ax = plt.subplots(figsize=(10.2, 6.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    ax.text(5, 9.55,
            'Reporting decisions for null findings in drug-target MR',
            ha='center', fontsize=12, weight='bold')
    ax.text(0.15, 6.8, 'H', ha='center', va='center', fontsize=11,
            color='#34495e')
    ax.text(0.15, 2.4, 'L', ha='center', va='center', fontsize=11,
            color='#34495e')
    ax.text(0.45, 8.0, 'power', ha='center', fontsize=8, color='#5d6d7e')
    ax.text(0.45, 3.6, 'power', ha='center', fontsize=8, color='#5d6d7e')
    ax.text(5, 0.28, 'Evidence consistency (colocalisation, pleiotropy, '
            'endpoint drift)  ->', ha='center', fontsize=9.5,
            color='#34495e')
    boxes = [
        (1.0, 5.35, 4.1, 3.5, '#fdecea', '#c0392b',
         'A  High power + consistent evidence',
         'Permitted: "no evidence of a causal effect; '
         'effects >= OR X can be excluded".\n'
         'Report the exclusion bound explicitly.'),
        (5.35, 5.35, 4.1, 3.5, '#fef5e7', '#b9770e',
         'B  High power + inconsistent evidence',
         'Permitted: "estimates are unstable; '
         'pleiotropy or phenotype heterogeneity suspected".\n'
         'No causal conclusion.'),
        (1.0, 1.15, 4.1, 3.5, '#eaf2f8', '#2471a3',
         'C  Low power + consistent evidence  <- this study',
         'Permitted: "no supportive evidence; available data '
         'insufficient to test a moderate effect".\n'
         'Forbidden: "the target has no effect".'),
        (5.35, 1.15, 4.1, 3.5, '#f4f6f7', '#5d6d7e',
         'D  Low power + inconsistent evidence',
         'Report only as "not interpretable".\n'
         'Redesign with a larger GWAS or independent '
         'instruments; no conclusion.'),
    ]
    for x, y, w, h, fc, ec, title, body in boxes:
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h,
            boxstyle='round,pad=0.10,rounding_size=0.12',
            facecolor=fc, edgecolor=ec, lw=1.4))
        ax.text(x + 0.22, y + h - 0.55, title, fontsize=9.4,
                color=ec, weight='bold', va='top')
        ax.text(x + 0.22, y + h - 1.15, body, fontsize=8.2,
                color='#1f2933', va='top', linespacing=1.5)
    ax.text(5, 9.05,
            'Assess power first, then evidence consistency, then constrain '
            'the wording. No quadrant permits a claim of "no effect" '
            'without a quantified exclusion bound.',
            ha='center', fontsize=8.4, color='#5d6d7e')
    fig.tight_layout()
    save(fig, 'fig4_decision_quadrant')


def fig5():
    d = json.load(io.open(RES + '/a1_simulation_results.json',
                          encoding='utf-8'))
    cal = d['calibration']
    pa, pb, pc = d['panel_A'], d['panel_B'], d['panel_C']
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.9))
    # Panel A
    ax = axes[0]
    x = np.array(pa['se_out'])
    ax.plot(x, np.array(pa['power_or_1.5']) * 100, color='#1f77b4', lw=2.2,
            label='Power for OR=1.5')
    ax.axhline(80, color='#e67e22', ls='--', lw=1.1)
    ax.text(0.0053, 82, '80% power', fontsize=8, color='#e67e22')
    ax.set_xscale('log')
    ax.set_xlabel('Outcome SE (log scale)')
    ax.set_ylabel('Power (%)', color='#1f77b4')
    ax.set_ylim(0, 105)
    ax2 = ax.twinx()
    ax2.plot(x, pa['mde_or_80pct'], color='#c0392b', lw=2.0,
             label='Minimum detectable OR')
    ax2.set_yscale('log')
    ax2.set_ylabel('Minimum detectable OR', color='#c0392b')
    base_so = cal['se_out_baseline']
    ax.plot([base_so], [cal['power_or_1.5_baseline'] * 100], 'o',
            color='#1f77b4', ms=7)
    ax.annotate('R13 primary\n%.2f%% / MDE %.2f'
                % (cal['power_or_1.5_baseline'] * 100,
                   cal['mde_or_80pct_baseline']),
                xy=(base_so, cal['power_or_1.5_baseline'] * 100),
                xytext=(0.045, 16), fontsize=7.5,
                arrowprops=dict(arrowstyle='->', lw=0.7))
    ax.set_title('(a) Power and minimum detectable effect', fontsize=10)
    ax.grid(alpha=0.2, lw=0.5)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=7.4,
              frameon=True)
    # Panel B
    ax = axes[1]
    xb = np.array([r['se_out'] for r in pb])
    ax.plot(xb, [r['naive_error_when_true_OR1.5'] * 100 for r in pb],
            color='#c0392b', lw=2.0, label='Naive: true OR=1.5 called "no effect"')
    ax.plot(xb, [r['naive_false_positive_when_null'] * 100 for r in pb],
            color='#e67e22', lw=1.6, ls=':', label='Naive: false positive under null')
    ax.plot(xb, [r['decision_error_when_true_OR1.5'] * 100 for r in pb],
            color='#1f77b4', lw=2.0, label='Framework: false assertion')
    ax.plot(xb, [100.0 if r['decision_inconclusive'] else 0.0 for r in pb],
            color='#2471a3', lw=1.4, ls='--',
            label='Framework: "not evaluable"')
    ax.set_xscale('log')
    ax.set_xlabel('Outcome SE (log scale)')
    ax.set_ylabel('Proportion (%)')
    ax.set_ylim(0, 102)
    ax.axvline(base_so, color='#95a5a6', ls='-.', lw=1.0)
    ax.text(base_so * 1.06, 4, 'Primary analysis', fontsize=7.5,
            color='#5d6d7e')
    ax.set_title('(b) Naive reporting vs decision framework '
                 '(%d replicates/point)' % d['n_sim_per_point'], fontsize=10)
    ax.legend(loc='center left', fontsize=7.0, frameon=True)
    ax.grid(alpha=0.2, lw=0.5)
    # Panel C: analytic binomial under the global null
    ax = axes[2]
    M, p0 = pc['M'], 0.05
    ks = np.arange(0, 26)
    pmf = binom.pmf(ks, M, p0)
    ax.bar(ks, pmf, color='#aec7e8', edgecolor='#1f77b4', lw=0.5,
           label='Binomial(133, 0.05) under global null')
    ax.axvline(pc['k_observed'], color='#c0392b', lw=2.2,
               label='Observed k = %d' % pc['k_observed'])
    ax.set_xlim(-0.5, 25.5)
    ax.set_xlabel('Number of P<0.05 results among 133 tests')
    ax.set_ylabel('Probability')
    ax.set_title('(c) Multiplicity creates "significant" findings',
                 fontsize=10)
    ax.text(0.98, 0.95,
            'Expected %.2f\nP(at least one nominal) = %.4f\n'
            'Binomial P (two-sided) = %.3f'
            % (pc['expected'], pc['prob_at_least_one_nominal'],
               pc['binom_p_two_sided']),
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(fc='white', ec='#cbd5e1', lw=0.5))
    ax.legend(loc='upper left', fontsize=7.4, frameon=True)
    ax.grid(alpha=0.2, lw=0.5, axis='y')
    fig.suptitle('Simulation validation of the null-reporting framework '
                 '(parameters calibrated to observed A1 data)',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, 'fig5_simulation_validation')


def main():
    style()
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    print('Saved English vector figures:')
    for f in sorted(os.listdir(RES)):
        if f.endswith('_EN.pdf'):
            print(' ', f)


if __name__ == '__main__':
    main()
