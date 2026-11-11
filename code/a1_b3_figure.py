# -*- coding: utf-8 -*-
"""B3 附录图：区域级共定位
fig9 (a) 区域轨迹：Fenland pQTL 与 FinnGen DPN 在同一 400 kb 区间的 -log10 P
fig9 (b) 单因果变异假设下五个后验概率（主结局 vs 阴性对照）
"""
import csv
import io
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, 'results')

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.linewidth': 0.8, 'pdf.fonttype': 42})

VARS = os.path.join(RES, 'a1_b3_coloc_variants.tsv')
CJ = os.path.join(RES, 'a1_b3_coloc_full.json')

rows = []
with io.open(VARS, encoding='utf-8') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        try:
            rows.append({'pos': int(r['pos']) / 1e6,
                         'pe': -np.log10(float(r['pval_exp'])),
                         'po': -np.log10(float(r['pval_out']))})
        except (ValueError, KeyError):
            continue

j = json.load(io.open(CJ, encoding='utf-8'))
main = [r for r in j['results'] if r['endpoint'] == 'DM_NEUROPATHY'][0]
neg = [r for r in j['results'] if r['endpoint'] == 'G6_POLYOTHUNS'][0]

fig, ax = plt.subplots(2, 1, figsize=(7.0, 5.6),
                       gridspec_kw={'height_ratios': [1.35, 1.0]})

# ---- (a) 区域轨迹
x = np.array([r['pos'] for r in rows])
pe = np.array([r['pe'] for r in rows])
po = np.array([r['po'] for r in rows])
ax[0].scatter(x, pe, s=6, c='#1f4e79', alpha=0.75, label='CX3CL1 pQTL (Fenland)', linewidths=0)
ax[0].scatter(x, po, s=6, c='#b03a2e', alpha=0.75, label='DM neuropathy (FinnGen R13)', linewidths=0)
ax[0].axhline(-np.log10(5e-8), color='#1f4e79', ls=':', lw=0.9)
ax[0].axhline(-np.log10(1e-5), color='#b03a2e', ls=':', lw=0.9)
ax[0].text(57.585, -np.log10(5e-8) + 1.6, 'pQTL 5×10$^{-8}$', fontsize=7,
           color='#1f4e79', ha='right')
ax[0].text(57.585, -np.log10(1e-5) - 3.0, 'outcome 10$^{-5}$', fontsize=7,
           color='#b03a2e', ha='right')
ax[0].axvspan(57.3725, 57.3851, color='#999999', alpha=0.18, lw=0)
ax[0].text(57.3788, pe.max() * 1.03, 'CX3CL1', ha='center',
           fontsize=8, color='#444444')
ax[0].set_ylim(0, pe.max() * 1.12)
ax[0].set_xlim(57.2, 57.6)
ax[0].set_xlabel('Position on chromosome 16 (Mb, GRCh38)')
ax[0].set_ylabel('$-\\log_{10}$ P')
ax[0].set_title('(a) Region-wide association, 1,744 shared variants', fontsize=10, loc='left')
ax[0].legend(frameon=False, fontsize=8, loc='upper right', markerscale=2.2)

# ---- (b) 后验概率
labels = ['H0', 'H1', 'H2', 'H3', 'H4']
keys = ['PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4']
mv = [main['coloc_W_0.04'][k] for k in keys]
nv = [neg['coloc_W_0.04'][k] for k in keys]
w = 0.36
xi = np.arange(len(labels))
ax[1].bar(xi - w / 2, mv, w, color='#b03a2e', label='DM neuropathy')
ax[1].bar(xi + w / 2, nv, w, color='#7f8c8d', label='Polyneuropathy (control)')
for i, v in enumerate(mv):
    if v >= 0.001:
        ax[1].text(i - w / 2, v + 0.02, '%.3f' % v, ha='center', fontsize=7)
for i, v in enumerate(nv):
    if v >= 0.20:
        ax[1].text(i + w / 2, v + 0.02, '%.3f' % v, ha='center', fontsize=7)
ax[1].set_xticks(xi)
ax[1].set_xticklabels(['H0\nno signal', 'H1\nexposure only', 'H2\noutcome only',
                       'H3\ntwo loci', 'H4\nshared'], fontsize=8)
ax[1].set_ylabel('Posterior probability')
ax[1].set_ylim(0, 1.42)
ax[1].set_title('(b) Colocalisation posteriors (single causal variant)', fontsize=10, loc='left')
ax[1].legend(frameon=False, fontsize=8, loc='upper center', ncol=2,
             bbox_to_anchor=(0.5, 1.02))
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)

fig.tight_layout()
for ext, kw in [('pdf', {}), ('svg', {}), ('png', {'dpi': 400})]:
    p = os.path.join(RES, 'fig9_b3_region_coloc_EN.%s' % ext)
    fig.savefig(p, bbox_inches='tight', **kw)
    print('saved', p)
print('PP.H1 main %.4f  PP.H4 main %.4f' % (main['coloc_W_0.04']['PP.H1'],
                                           main['coloc_W_0.04']['PP.H4']))
