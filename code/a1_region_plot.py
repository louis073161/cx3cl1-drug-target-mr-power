# -*- coding: utf-8 -*-
"""CX3CL1 cis 区段：DPN 与 pQTL 关联信号对照图（全部真实数据）。"""
import json
import math
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
sys.path.insert(0, ROOT + '/code')
from a1_ld_coloc import load_fg_region  # noqa: E402

OUT = ROOT + '/results/fig1_cx3cl1_region.png'


def main():
    fg = load_fg_region()
    cs = json.load(open(ROOT + '/data/ukbppp_cx3cl1_crediblesets.json',
                        encoding='utf-8'))['data']['credibleSets']['rows']
    cis = [r for r in cs if not r['isTransQtl']]

    fig, ax = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True,
                           gridspec_kw={'height_ratios': [1.6, 1]})
    x = [r['pos'] / 1e6 for r in fg]
    y = [-math.log10(r['pval']) for r in fg]
    ax[0].scatter(x, y, s=6, c='#8c9bab', alpha=0.75, linewidths=0)
    for r in cis:
        ax[0].axvline(r['position'] / 1e6, color='#c0392b', lw=1.2, ls='--',
                      alpha=0.85)
    ax[0].axhline(5, color='#e67e22', lw=0.9, ls=':')
    ax[0].text(57.205, 5.12, 'P = 1e-5 (no variant reaches it)',
               fontsize=8, color='#e67e22')
    top = max(fg, key=lambda r: -r['pval'])
    best = min(fg, key=lambda r: r['pval'])
    ax[0].annotate('%s  P=%.2g' % (best['rsids'], best['pval']),
                   xy=(best['pos'] / 1e6, -math.log10(best['pval'])),
                   xytext=(57.30, 3.3), fontsize=8,
                   arrowprops=dict(arrowstyle='->', lw=0.7, color='#555'))
    ax[0].set_ylabel('-log10(P)  FinnGen R13 DPN', fontsize=9)
    ax[0].set_title('CX3CL1 cis region (chr16, GRCh38): no DPN association '
                    'signal at the pQTL locus', fontsize=10)
    ax[0].set_ylim(0, 5.6)
    ax[0].grid(axis='y', alpha=0.25, lw=0.5)

    for k, r in enumerate(cis):
        p = r['pValueMantissa'] * (10 ** r['pValueExponent'])
        lp = -math.log10(p) if p > 0 else 0
        ax[1].vlines(r['position'] / 1e6, 0, lp, color='#c0392b', lw=3)
        label = r['variant']['rsIds'][0] if r['variant']['rsIds'] else ''
        ax[1].text(r['position'] / 1e6 + 0.004, lp + 6,
                   '%s  (P=%.1g)' % (label, p),
                   fontsize=7, ha='left', va='bottom', color='#c0392b')
    ax[1].set_ylabel('-log10(P)  UKB-PPP CX3CL1 pQTL', fontsize=9)
    ax[1].set_xlabel('chr16 position (Mb, GRCh38)', fontsize=9)
    ax[1].set_ylim(0, 215)
    ax[1].grid(axis='y', alpha=0.25, lw=0.5)
    for a in ax:
        a.set_xlim(57.2, 57.6)
        a.tick_params(labelsize=8)
    plt.tight_layout()
    plt.savefig(OUT, dpi=220, facecolor='white')
    print('Saved', OUT)


if __name__ == '__main__':
    main()
