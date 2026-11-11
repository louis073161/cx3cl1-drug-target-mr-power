# -*- coding: utf-8 -*-
"""A1 仿真验证：低效能 / 端点异质性下，naive 报告与本决策流程的行为对比。

诚信声明：本脚本产生的是**方法学仿真**，不产生新的实证结论、不替代任何
真实数据分析。所有标定参数均来自 A1 已落盘的真实结果（见 CALIBRATION），
随机数种子固定，可一键复现。

三个面板：
  A  效能与最小可检出效应随结局侧 SE 的变化（解析式 + 实测点标定）
  B  蒙特卡洛对比：naive 报告 vs 决策流程的错误断言率与"不可判定率"
  C  蒙特卡洛：133 次（19 组合 x 7 方法）检验在全局 null 下的名义显著数分布，
     与实测 8 个、期望 6.65、二项 P=0.548 对照
"""
import json
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest, norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
RES = ROOT + '/results'
SEED = 20261110
N_SIM = 20000
Z975, Z80 = norm.ppf(0.975), norm.ppf(0.80)

# ---------------- 全部标定参数来自真实结果，不臆测 ----------------
CALIBRATION = {
    'beta_exp': 0.13660051070809745,       # UKB-PPP CX3CL1 cis lead, OT v26.9
    'se_exp': 0.005811879077927919,
    'z_exp': 23.50367391965749,
    'F_instrument': 552.42,                # = z_exp^2，见 a1_robustness.json
    'se_out_baseline': 0.1513,             # = 1.1072872 * beta_exp（R13 DPN）
    'se_wald_baseline': 1.1072872218114906,
    'power_or_1.5_baseline': 0.06549754388078935,
    'mde_or_80pct_baseline': 22.245940220211136,
    'se_wald_needed_80pct': math.log(1.5) / (Z975 + Z80),
    'n_combos': 19, 'n_methods': 133,
    'n_nominal_sig_observed': 8,
    'expected_nominal_sig': 133 * 0.05,
    'binom_p_observed': float(binomtest(8, 133, 0.05).pvalue),
    'bonferroni_threshold': 0.05 / 133,
    'or_range_observed': [0.350, 3.755],
    'or_ratio_observed': 10.72,
    'positive_control_best_cis_se_wald': 0.2683,   # I9_CVD_HARD, lead cis
    'positive_control_best_cis_power': 0.327,
    'positive_control_min_p': 1.37e-19,
}
BETA = CALIBRATION['beta_exp']
SE_OUT_NEEDED = CALIBRATION['se_wald_needed_80pct'] * BETA


def power_at_or(se, or_true):
    z = math.log(or_true) / se
    return float(norm.cdf(z - Z975) + norm.cdf(-z - Z975))


def mde_or(se):
    return math.exp((Z975 + Z80) * se)


def panel_a_data():
    xs = np.logspace(math.log10(0.005), math.log10(0.20), 160)
    se_wald = xs / BETA
    power = np.array([power_at_or(s, 1.5) for s in se_wald])
    mde = np.array([mde_or(s) for s in se_wald])
    return {'se_out': xs, 'se_wald': se_wald, 'power': power, 'mde': mde}


def panel_b_data(rng):
    """对每个 se_out 网格点，模拟 N_SIM 次两样本 Wald ratio 估计。"""
    xs = np.logspace(math.log10(0.005), math.log10(0.20), 22)
    rec = []
    for so in xs:
        se_w = so / BETA
        pw = power_at_or(se_w, 1.5)
        # 真效应情景：beta_true = ln(1.5)
        hat_true = rng.normal(math.log(1.5), se_w, N_SIM)
        p_true = 2 * norm.sf(np.abs(hat_true) / se_w)
        naive_err_true = float(np.mean(p_true >= 0.05))
        # 零效应情景
        hat_null = rng.normal(0.0, se_w, N_SIM)
        p_null = 2 * norm.sf(np.abs(hat_null) / se_w)
        naive_fp_null = float(np.mean(p_null < 0.05))
        # 决策流程：先判效能，效能不足则"不可判定"，不作确断断言
        decisive = pw >= 0.80
        decision_err_true = float(np.mean(p_true >= 0.05)) if decisive else 0.0
        decision_fp_null = float(np.mean(p_null < 0.05)) if decisive else 0.0
        rec.append({
            'se_out': float(so), 'se_wald': float(se_w),
            'power_or_1.5': pw, 'mde_or_80pct': mde_or(se_w),
            'naive_error_when_true_OR1.5': naive_err_true,
            'naive_false_positive_when_null': naive_fp_null,
            'decision_error_when_true_OR1.5': decision_err_true,
            'decision_false_positive_when_null': decision_fp_null,
            'decision_inconclusive': not decisive,
        })
    return rec


def panel_c_data(rng):
    M = CALIBRATION['n_methods']
    k = rng.binomial(M, 0.05, N_SIM)
    bt = binomtest(CALIBRATION['n_nominal_sig_observed'], M, 0.05)
    # Bonferroni 阈值下至少 1 个显著的概率
    p_bonf = CALIBRATION['bonferroni_threshold']
    prob_any_bonf = 1 - (1 - p_bonf) ** M
    prob_any_nominal = 1 - 0.95 ** M
    return {
        'M': M, 'N_sim': N_SIM, 'k_samples': k,
        'k_mean': float(k.mean()), 'k_sd': float(k.std()),
        'k_observed': CALIBRATION['n_nominal_sig_observed'],
        'expected': M * 0.05,
        'binom_p_two_sided': float(bt.pvalue),
        'prob_at_least_one_nominal': prob_any_nominal,
        'prob_at_least_one_bonferroni': prob_any_bonf,
        'p_gt_observed': float(np.mean(k >= CALIBRATION['n_nominal_sig_observed'])),
        'q95': int(np.quantile(k, 0.95)),
    }


def make_plot(a, b, c):
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams['text.color'] = '#1f2933'
    plt.rcParams['axes.labelcolor'] = '#1f2933'
    plt.rcParams['xtick.color'] = '#1f2933'
    plt.rcParams['ytick.color'] = '#1f2933'

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.9))

    # ---- A：效能 / MDE ----
    ax = axes[0]
    x = a['se_out']
    ax.plot(x, a['power'] * 100, color='#1f77b4', lw=2.2,
            label='OR=1.5 的把握度（左轴）')
    ax.axhline(80, color='#e67e22', ls='--', lw=1.1)
    ax.text(0.0053, 82, '80% 把握度', fontsize=8, color='#e67e22')
    ax.set_xscale('log')
    ax.set_xlabel('结局侧标准误 se_out（GRCh38 汇总统计，对数刻度）', fontsize=9)
    ax.set_ylabel('把握度（%）', fontsize=9, color='#1f77b4')
    ax.set_ylim(0, 105)
    ax2 = ax.twinx()
    ax2.plot(x, a['mde'], color='#c0392b', lw=2.0, ls='-',
             label='80% 把握度最小可检出 OR（右轴）')
    ax2.set_yscale('log')
    ax2.set_ylabel('最小可检出 OR', fontsize=9, color='#c0392b')
    ax2.tick_params(axis='y', colors='#c0392b')
    # 实测点
    base_so = CALIBRATION['se_out_baseline']
    ax.plot([base_so], [CALIBRATION['power_or_1.5_baseline'] * 100], 'o',
            color='#1f77b4', ms=7, zorder=5)
    ax.annotate('R13 DPN 主分析\n6.55% / MDE 22.25',
                xy=(base_so, 6.55), xytext=(0.045, 16),
                fontsize=7.8, color='#1f2933',
                arrowprops=dict(arrowstyle='->', color='#1f2933', lw=0.8))
    pc_so = 0.03664  # = 0.2683 * beta_exp，I9_CVD_HARD lead cis 实测
    ax.plot([pc_so], [CALIBRATION['positive_control_best_cis_power'] * 100], 's',
            color='#1f77b4', ms=7, zorder=5)
    ax.annotate('大样本端点最优\n32.7% / MDE 2.12',
                xy=(pc_so, 32.7), xytext=(0.0062, 62),
                fontsize=7.8, color='#1f2933',
                arrowprops=dict(arrowstyle='->', color='#1f2933', lw=0.8))
    ax.axvline(SE_OUT_NEEDED, color='#7f8c8d', ls=':', lw=1.2)
    ax.text(SE_OUT_NEEDED * 1.06, 3,
            '80%% 所需 se_out≈%.4f' % SE_OUT_NEEDED,
            fontsize=7, color='#7f8c8d', rotation=90, va='bottom')
    ax.set_title('A  效能与最小可检出效应（解析式，实测点标定）',
                 fontsize=10, pad=10)
    ax.grid(alpha=0.2, lw=0.5)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=7.5,
              frameon=True, facecolor='white', edgecolor='#cccccc')

    # ---- B：naive vs 决策流程 ----
    ax = axes[1]
    xb = [r['se_out'] for r in b]
    ax.plot(xb, [r['naive_error_when_true_OR1.5'] * 100 for r in b],
            color='#c0392b', lw=2.2, label='naive：真有 OR=1.5 却断言"无效"')
    ax.plot(xb, [r['naive_false_positive_when_null'] * 100 for r in b],
            color='#e67e22', lw=1.6, ls='--', label='naive：真无效应却断言"有效"')
    ax.plot(xb, [r['decision_error_when_true_OR1.5'] * 100 for r in b],
            color='#1f77b4', lw=2.2, label='决策流程：错误断言率')
    ax.plot(xb, [100.0 if r['decision_inconclusive'] else 0.0 for r in b],
            color='#7f8c8d', lw=1.5, ls=':', label='决策流程：判定"不可判定"')
    ax.set_xscale('log')
    ax.set_xlabel('结局侧标准误 se_out（对数刻度）', fontsize=9)
    ax.set_ylabel('比例（%）', fontsize=9)
    ax.set_ylim(0, 102)
    ax.axvline(base_so, color='#95a5a6', ls='-.', lw=1.0)
    ax.text(base_so * 1.06, 4, 'A1 主分析', fontsize=7.5,
            color='#7f8c8d', rotation=90, va='bottom')
    ax.set_title('B  naive 报告 vs 决策流程（%d 次/点仿真）' % N_SIM,
                 fontsize=10, pad=10)
    ax.legend(loc='center left', fontsize=7.3, frameon=True,
              facecolor='white', edgecolor='#cccccc')
    ax.grid(alpha=0.2, lw=0.5)

    # ---- C：多重检验噪声 ----
    ax = axes[2]
    k = c['k_samples']
    bins = np.arange(-0.5, 26.5, 1)
    ax.hist(k, bins=bins, density=True, color='#aec7e8',
            edgecolor='white', lw=0.6, label='仿真分布 Bin(133, 0.05)')
    ax.axvline(c['expected'], color='#1f77b4', ls='--', lw=1.3,
               label='期望 6.65')
    ax.axvline(c['k_observed'], color='#c0392b', lw=2.2,
               label='实测 8 个名义显著')
    ax.set_xlim(-0.5, 25.5)
    ax.set_xlabel('133 次检验中 P<0.05 的个数', fontsize=9)
    ax.set_ylabel('概率密度', fontsize=9)
    ax.set_title('C  方法自由度制造的"显著"（全局 null 仿真）',
                 fontsize=10, pad=10)
    ax.text(0.98, 0.95,
            '二项检验 P = %.3f\n至少 1 个名义显著的概率 = %.1f%%\n'
            'Bonferroni 阈值 %.2g 下实测 0 个'
            % (c['binom_p_two_sided'], c['prob_at_least_one_nominal'] * 100,
               CALIBRATION['bonferroni_threshold']),
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='#cccccc'))
    ax.legend(loc='upper left', fontsize=7.5, frameon=True,
              facecolor='white', edgecolor='#cccccc')
    ax.grid(alpha=0.2, lw=0.5)

    fig.suptitle('图5　药物靶点 MR 阴性报告决策流程的仿真验证（参数由 A1 实测标定；'
                 '仿真不产生新实证结论）', fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = RES + '/fig5_simulation_validation.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out


def main():
    os.makedirs(RES, exist_ok=True)
    rng = np.random.default_rng(SEED)
    a = panel_a_data()
    b = panel_b_data(rng)
    c = panel_c_data(rng)
    out_png = make_plot(a, b, c)

    payload = {
        'type': 'methodological simulation only — NOT new empirical evidence',
        'seed': SEED, 'n_sim_per_point': N_SIM,
        'calibration': CALIBRATION,
        'derived': {
            'se_out_needed_for_80pct_at_OR1.5': SE_OUT_NEEDED,
            'se_wald_needed_for_80pct_at_OR1.5':
                CALIBRATION['se_wald_needed_80pct'],
            'se_fold_needed': CALIBRATION['se_wald_baseline'] /
                              CALIBRATION['se_wald_needed_80pct'],
        },
        'panel_A': {'se_out': a['se_out'].tolist(),
                    'se_wald': a['se_wald'].tolist(),
                    'power_or_1.5': a['power'].tolist(),
                    'mde_or_80pct': a['mde'].tolist()},
        'panel_B': b,
        'panel_C': {k: v for k, v in c.items() if k != 'k_samples'},
        'honesty': [
            'All calibration parameters come from real downloaded/executed '
            'A1 results; no fabricated effect sizes or sample sizes.',
            'Panel A is analytic, Panels B and C are Monte Carlo with a '
            'fixed seed; they validate the decision workflow behaviour, '
            'not the CX3CL1 hypothesis.',
            'F=552.42 (lead cis z=23.50) is the instrument strength; '
            '22.25 is the minimum detectable OR, NOT the F statistic.',
        ],
    }
    with open(RES + '/a1_simulation_results.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    # 控制台摘要
    print('=== 标定参数（全部来自真实结果）===')
    print('  beta_exp=%.6f  F=%.2f  se_out(R13 DPN)=%.4f'
          % (BETA, CALIBRATION['F_instrument'],
             CALIBRATION['se_out_baseline']))
    print('  SE_wald=%.4f  power(OR1.5)=%.4f  MDE=%.2f'
          % (CALIBRATION['se_wald_baseline'],
             CALIBRATION['power_or_1.5_baseline'],
             CALIBRATION['mde_or_80pct_baseline']))
    print('  80%% 把握度需 se_out<=%.5f（当前需 %.2f 倍 SE 改善）'
          % (SE_OUT_NEEDED, payload['derived']['se_fold_needed']))
    print('\n=== B 面板（se_out=0.1513 处）===')
    near = min(b, key=lambda r: abs(r['se_out'] - CALIBRATION['se_out_baseline']))
    print('  se_out=%.4f power=%.3f' % (near['se_out'], near['power_or_1.5']))
    print('  naive 真效应误判"无效" = %.1f%%' %
          (near['naive_error_when_true_OR1.5'] * 100))
    print('  naive 零效应假阳性     = %.1f%%' %
          (near['naive_false_positive_when_null'] * 100))
    print('  决策流程错误断言       = %.1f%%（判定"不可判定"=%s）' %
          (near['decision_error_when_true_OR1.5'] * 100,
           near['decision_inconclusive']))
    print('\n=== C 面板 ===')
    print('  K 均值=%.3f  SD=%.3f  期望=%.2f  95%%分位=%d'
          % (c['k_mean'], c['k_sd'], c['expected'], c['q95']))
    print('  实测 K=8 的二项检验 P=%.4f' % c['binom_p_two_sided'])
    print('  M=133 时至少 1 个名义显著的概率 = %.2f%%'
          % (c['prob_at_least_one_nominal'] * 100))
    print('\nSaved:', out_png, RES + '/a1_simulation_results.json')


if __name__ == '__main__':
    main()
