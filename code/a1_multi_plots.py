# -*- coding: utf-8 -*-
"""A1 方法学主线三张图：跨组合森林图、效能曲线、阴性报告决策四象限。"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
RES = ROOT + '/results'
sys.path.insert(0, ROOT + '/code')

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['text.color'] = '#1f2933'
plt.rcParams['axes.labelcolor'] = '#1f2933'
plt.rcParams['xtick.color'] = '#1f2933'
plt.rcParams['ytick.color'] = '#1f2933'

Z975, Z80 = norm.ppf(0.975), norm.ppf(0.80)
COLORS = {'primary': '#c0392b', 'proxy_painful_DPN': '#8e44ad',
          'broad_polyneuropathy': '#2471a3',
          'negative_ctrl_microvascular': '#7f8c8d'}
LABELS = {'primary': '主结局 DPN', 'proxy_painful_DPN': '痛性代理',
          'broad_polyneuropathy': '泛化多发性神经病',
          'negative_ctrl_microvascular': '阴性对照（微血管并发症）'}


def forest():
    m = json.load(open(RES + '/a1_multiendpoint_mr.json', encoding='utf-8'))
    rows = []
    for c in m['combos']:
        p = c['methods'].get('Wald ratio (lead cis)')
        if p:
            rows.append((c['release'], c['endpoint'], c['endpoint_type'], p))
    rows.sort(key=lambda x: (int(x[0][1:]), x[1]))
    rows = rows[::-1]  # 顶部为 R9、底部为 R13，符合时间顺序阅读习惯

    fig, ax = plt.subplots(figsize=(11, 8.8))
    y = list(range(len(rows)))
    for i, (rel, ep, et, p) in enumerate(rows):
        col = COLORS.get(et, '#555555')
        lo, hi = math.exp(p['ci_lo']), math.exp(p['ci_hi'])
        ax.plot([lo, hi], [i, i], color=col, lw=1.8, alpha=0.85, zorder=2)
        ax.scatter([p['or']], [i], s=42, color=col, zorder=3,
                   edgecolors='white', linewidths=0.7)
        ax.text(58, i, '%.2f (%.2f–%.2f)' % (p['or'], lo, hi),
                va='center', ha='left', fontsize=7.6, color='#34495e')
    ax.axvline(1, color='#2c3e50', lw=1.1, ls='--', zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(['%s  %s' % (r[0], r[1]) for r in rows], fontsize=8.2)
    ax.set_xscale('log')
    ax.set_xlim(0.004, 55)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlabel('OR（对数刻度，Wald ratio，单一 cis lead 工具变量）', fontsize=9.5)
    ax.set_title('图1  同一工具变量 × 19 个（release, 端点）组合：点估计与 95% 置信区间',
                 fontsize=11, pad=12)
    ax.grid(axis='x', alpha=0.22, lw=0.5)
    handles = [plt.Line2D([0], [0], color=v, lw=3, label=LABELS[k])
               for k, v in COLORS.items()]
    ax.legend(handles=handles, loc='upper center',
              bbox_to_anchor=(0.42, -0.105), ncol=4, fontsize=8.2,
              frameon=True, facecolor='white', edgecolor='#cccccc')
    fig.text(0.012, 0.012, '虚线为 OR=1；19 个组合主分析 P 均 >0.05；'
             '右侧为 OR（95% CI）', fontsize=8, color='#555555')
    fig.tight_layout(rect=(0, 0.075, 0.80, 1))
    out = RES + '/fig2_forest_multiendpoint.png'
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print('saved', out)


def power_curve():
    m = json.load(open(RES + '/a1_multiendpoint_mr.json', encoding='utf-8'))
    picks = [
        ('R13', 'DM_NEUROPATHY', '#c0392b', '主结局 R13 DPN'),
        ('R11', 'DM_NEUROPATHY', '#e67e22', '主结局 R11 DPN'),
        ('R13', 'G6_POLYOTHUNS', '#2471a3', '泛化神经病变 R13'),
        ('R13', 'DM_RETINOPATHY', '#7f8c8d', '效能最高：视网膜病变 R13'),
        ('R9', 'DM_POLYNEURO', '#8e44ad', '效能最低：痛性代理 R9'),
    ]
    ses = {}
    for c in m['combos']:
        if c['power']:
            ses[(c['release'], c['endpoint'])] = c['power']['se']

    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    xs = [1.0 + 0.02 * i for i in range(126)]  # 1.0 – 3.5
    for rel, ep, col, lab in picks:
        se = ses.get((rel, ep))
        if not se:
            continue
        z = [math.log(x) / se for x in xs]
        pw = [norm.cdf(v - Z975) + norm.cdf(-v - Z975) for v in z]
        ax.plot(xs, pw, color=col, lw=2, label=lab)
        p15 = norm.cdf(math.log(1.5) / se - Z975) + norm.cdf(
            -math.log(1.5) / se - Z975)
        ax.scatter([1.5], [p15], s=32, color=col, zorder=4)
    ax.axhline(0.80, color='#2c3e50', lw=1.2, ls='--')
    ax.text(3.48, 0.815, '80% 把握度', ha='right', fontsize=8.5,
            color='#2c3e50')
    ax.axvline(1.5, color='#95a5a6', lw=0.9, ls=':')
    ax.text(1.51, 0.015, 'OR=1.5', fontsize=8.2, color='#5d6d7e')
    ax.set_xlim(1.0, 3.5)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel('假定的真实 OR', fontsize=9.5)
    ax.set_ylabel('双侧 α=0.05 的检验把握度', fontsize=9.5)
    ax.set_title('图2  效能曲线（由各组合实测 Wald ratio SE 驱动，非病例数假设）',
                 fontsize=11, pad=10)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc='upper left', fontsize=8.4, frameon=True,
              facecolor='white', edgecolor='#cccccc')
    ax.text(2.05, 0.035,
            '即使公开域中效能最高的组合，OR=1.5 的把握度也仅约 10%',
            fontsize=8.0, color='#34495e')
    fig.tight_layout()
    out = RES + '/fig3_power_curves.png'
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print('saved', out)


def decision_quadrant():
    fig, ax = plt.subplots(figsize=(10.2, 6.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    ax.text(5, 9.55, '图3  药物靶点 MR「阴性结果」报告决策四象限',
            ha='center', fontsize=12, weight='bold')
    ax.text(0.15, 5, '高\n检\n验\n效\n能', ha='center', va='center',
            fontsize=9.5, color='#34495e')
    ax.text(0.15, 2.2, '低\n检\n验\n效\n能', ha='center', va='center',
            fontsize=9.5, color='#34495e')
    ax.text(5, 0.28, '→ 证据一致性（共定位 / 多效性 / 端点定义漂移）',
            ha='center', fontsize=9.5, color='#34495e')

    boxes = [
        (1.0, 5.35, 4.1, 3.5, '#fdecea', '#c0392b',
         'A  高效能 + 证据一致',
         '允许：「未发现因果效应证据，\n可排除 OR≥X 的效应」\n必须报告可排除效应下界'),
        (5.35, 5.35, 4.1, 3.5, '#fef5e7', '#b9770e',
         'B  高效能 + 证据不一致',
         '允许：「结果不稳定，提示\n多效性或表型异质性」\n禁止下因果结论'),
        (1.0, 1.15, 4.1, 3.5, '#eaf2f8', '#2471a3',
         'C  低效能 + 证据一致   ← 本研究',
         '允许：「未获得支持性证据，\n现有数据不足以检验中等效应」\n禁止：「靶点无效」'),
        (5.35, 1.15, 4.1, 3.5, '#f4f6f7', '#5d6d7e',
         'D  低效能 + 证据不一致',
         '仅可报告为「不可解释」\n建议：换更大 GWAS / 更多独立\n工具变量后再设计，不作结论'),
    ]
    for x, y, w, h, fc, ec, title, body in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle='round,pad=0.10,rounding_size=0.12',
                                    facecolor=fc, edgecolor=ec, lw=1.4))
        ax.text(x + 0.22, y + h - 0.55, title, fontsize=9.6,
                color=ec, weight='bold', va='top')
        ax.text(x + 0.22, y + h - 1.15, body, fontsize=8.4,
                color='#1f2933', va='top', linespacing=1.5)
    ax.text(5, 9.05,
            '先判效能 → 再判证据一致性 → 最后限定措辞；任何象限都不得越过「未获得支持性证据」这一层。',
            ha='center', fontsize=8.6, color='#5d6d7e')
    fig.tight_layout()
    out = RES + '/fig4_decision_quadrant.png'
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print('saved', out)


if __name__ == '__main__':
    forest()
    power_curve()
    decision_quadrant()
