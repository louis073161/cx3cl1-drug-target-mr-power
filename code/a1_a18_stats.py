# -*- coding: utf-8 -*-
"""
A 级 18 项清零 · 实算部分（对应 STROBE-MR 11a / 12b / 13c / 15）

只用本地已落盘的真实公开数据重算，不联网、不补缺失、不编数值：
  results/a1_outcome_matrix.tsv        19 个 (release, endpoint) × 工具变量
  results/a1_round1_variants.tsv       暴露侧 β/SE/z/P/EAF
  results/a1_robustness.json           LD clumping、coloc 先验、数据源元数据
  results/a1_endpoint_dict.tsv         各组合官方病例/对照数（无则为空）
  results/a1_mr_matrix.tsv             既有主分析 Wald ratio 结果

新增统计量：
  F = z_exp^2                                  工具变量强度
  Cochran Q / Q_df / Q_p / I^2                IVW 异质性
  I^2_GX (Bowden 2016)                        MR-Egger 的 NOME 诊断
  Steiger 方向检验（近似实现，见下方公式说明）  exp→out 方向
  E-value (VanderWeele & Ding 2017)           残余混杂敏感性
  Bootstrap CI（工具变量层重抽样）             IVW 百分位区间
  ABF (Wakefield 近似)                        数据对 H1/H0 的支持
  lead 位点三组关联表                        变异-暴露 / 变异-结局 / Wald ratio

重要口径：
  * Steiger 为本脚本基于 t^2 = R^2(N-2)/(1-R^2) 恒等式的透明近似，
    不是 Hemani 官方 steiger 包实现；二分类结局为观测尺度近似，未做
    liability 转换。为规避 R13 无官方 manifest，额外给出"最保守边界"
    N_out = 2*n_cases（对照不少于病例时有效样本量的下界）。
  * E-value 针对主分析点估计与其 95% CI；CI 跨 1 时 CI 侧 E-value = 1。
  * Bootstrap 只对工具变量重抽样，不反映 GWAS 抽样误差，不得当作
    主分析置信区间的替代。
"""
import csv
import io
import json
import math
import random
from collections import defaultdict

from scipy.stats import chi2, norm

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
RES = ROOT + '/results'
LEAD = '16_57379748_C_T'
N_EXP = 33436                 # UKB-PPP CX3CL1 pQTL 样本量（a1_robustness.json）
PRIOR_W = 0.04                # 与既有 coloc 同一先验方差
N_BOOT = 5000
SEED = 20261110
Z975 = 1.959964


def read_tsv(path):
    with io.open(path, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_tsv(path, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0].keys())
    with io.open(path, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in fields})


def r2_from_z(z, n):
    """由单变量回归 t/R^2 恒等式反推解释度：R^2 = t^2/(t^2 + N - 2)。"""
    if n is None or n <= 2:
        return None
    z2 = z * z
    return z2 / (z2 + n - 2.0)


def steiger(r2_exp, r2_out, n_exp, n_out):
    """独立样本相关系数差异的 Fisher z 检验（近似）。"""
    if None in (r2_exp, r2_out, n_out):
        return None
    r_exp, r_out = math.sqrt(max(r2_exp, 0)), math.sqrt(max(r2_out, 0))
    z1, z2 = math.atanh(min(r_exp, 0.999999)), math.atanh(min(r_out, 0.999999))
    se = math.sqrt(1.0 / (n_exp - 3) + 1.0 / (n_out - 3))
    z = (z1 - z2) / se
    return {'r_exp': r_exp, 'r_out': r_out,
            'z_diff': z1 - z2, 'z_stat': z,
            'p_two_sided': float(2 * norm.sf(abs(z))),
            'direction_supports_exp_to_out': bool(z > 0 and 2 * norm.sf(abs(z)) < 0.05)}


def e_value(orr):
    """E-value = RR + sqrt(RR(RR-1))；OR<1 先取倒数；RR<=1 时为 1。"""
    if orr is None or not math.isfinite(orr) or orr <= 0:
        return None
    rr = orr if orr >= 1 else 1.0 / orr
    if rr <= 1:
        return 1.0
    return rr + math.sqrt(rr * (rr - 1))


def ivw_slope(rows):
    w = [1.0 / r['se_out'] ** 2 for r in rows]
    denom = sum(a * r['beta_exp'] ** 2 for a, r in zip(w, rows))
    if denom <= 0:
        return None
    return sum(a * r['beta_exp'] * r['beta_out'] for a, r in zip(w, rows)) / denom


def ivw_q(rows):
    w = [1.0 / r['se_out'] ** 2 for r in rows]
    b = ivw_slope(rows)
    if b is None:
        return None
    q = sum(a * (r['beta_out'] - b * r['beta_exp']) ** 2
            for a, r in zip(w, rows))
    df = len(rows) - 1
    return {'Q': q, 'Q_df': df, 'Q_p': float(chi2.sf(q, df)),
            'I2_pct': max(0.0, (q - df) / q * 100.0) if q > 0 else 0.0,
            'beta_ivw': b}


def i2_gx(rows):
    """Bowden 2016 I^2_GX：以 1/se_exp^2 为权重的暴露效应异质性。"""
    w = [1.0 / r['se_exp'] ** 2 for r in rows]
    sw = sum(w)
    mu = sum(a * r['beta_exp'] for a, r in zip(w, rows)) / sw
    q = sum(a * (r['beta_exp'] - mu) ** 2 for a, r in zip(w, rows))
    df = len(rows) - 1
    if q <= 0:
        return {'Q_GX': q, 'df': df, 'I2_GX_pct': 0.0}
    return {'Q_GX': q, 'df': df,
            'I2_GX_pct': max(0.0, (q - df) / q * 100.0)}


def bootstrap_ivw(rows, n_boot=N_BOOT, seed=SEED):
    rng = random.Random(seed)
    vals = []
    for _ in range(n_boot):
        ix = [rng.randrange(len(rows)) for _ in range(len(rows))]
        v = ivw_slope([rows[i] for i in ix])
        if v is not None and math.isfinite(v):
            vals.append(v)
    if len(vals) < 100:
        return None
    vals.sort()
    return {'n_boot': len(vals),
            'beta_boot_mean': sum(vals) / len(vals),
            'ci_lo': vals[int(0.025 * len(vals))],
            'ci_hi': vals[int(0.975 * len(vals))]}


def log_abf(beta, se, W=PRIOR_W):
    V = se * se
    r = W / (V + W)
    return 0.5 * math.log(1.0 - r) + 0.5 * (beta / se) ** 2 * r


def main():
    outcome = read_tsv(RES + '/a1_outcome_matrix.tsv')
    variants = read_tsv(RES + '/a1_round1_variants.tsv')
    ed = {(r['release'], r['endpoint']): r
          for r in read_tsv(RES + '/a1_endpoint_dict.tsv')}
    robust = json.load(io.open(RES + '/a1_robustness.json', encoding='utf-8'))
    clump_keep = set(robust['clumping']['r2_0.001']['kept_ids'])
    vmap = {v['variant']: v for v in variants}

    # R13 无官方 manifest；主结局病例数仅见于 r13 专属端点页（4,145），
    # 且与通用页面的 4,412 不一致。此处只用于保守边界，并显式标注来源。
    endpoint_def = read_tsv(RES + '/a1_endpoint_definition_table.tsv')
    fallback_cases = {}
    for r in endpoint_def:
        if (r.get('release') == 'R13' and r.get('端点') == 'DM_NEUROPATHY'
                and r.get('病例数', '').isdigit()):
            fallback_cases[('R13', 'DM_NEUROPATHY')] = {
                'n_cases': int(r['病例数']),
                'source': 'r13.risteys.finngen.fi/endpoints/DM_NEUROPATHY '
                          '(非官方 manifest；通用页 4,412 例与之不一致)'}

    # 暴露侧 P 值需由 mantissa/exponent 还原（仅纳入 lead 关联表）
    def p_exp(v):
        m, e = num(v['p_exp_mantissa']), num(v['p_exp_exponent'])
        return m * (10.0 ** e) if (m is not None and e is not None) else None

    # 数值化并分组
    groups = defaultdict(list)
    for r in outcome:
        x = dict(r)
        for k in ['beta_exp', 'se_exp', 'z_exp', 'beta_out', 'se_out',
                  'pval_out', 'af_alt']:
            x[k] = num(x[k])
        x['flip'] = (x['flip'] == 'True')
        x['approx'] = (x['approx'] == 'True')
        groups[(x['release'], x['endpoint'])].append(x)

    stats_rows, steiger_rows, lead_rows, bias_rows = [], [], [], []

    for (rel, ep), rows in sorted(groups.items(), key=lambda kv: (int(kv[0][0][1:]), kv[0][1])):
        meta = ed.get((rel, ep), {})
        nc, nco = num(meta.get('n_cases')), num(meta.get('n_controls'))
        fb = fallback_cases.get((rel, ep))
        n_total = (nc + nco) if (nc and nco) else None
        n_eff = (4.0 / (1.0 / nc + 1.0 / nco)) if (nc and nco) else None
        nc_cons = nc if nc else (fb['n_cases'] if fb else None)
        n_conservative = (2.0 * nc_cons) if nc_cons else None  # 有效样本量下界

        clumped = [r for r in rows if r['variant'] in clump_keep]
        lead = [r for r in rows if r['variant'] == LEAD]

        rec = {'release': rel, 'endpoint': ep,
               'endpoint_type': meta.get('endpoint_type', ''),
               'n_cases_official': int(nc) if nc else '',
               'n_controls_official': int(nco) if nco else '',
               'n_snp_available': len(rows), 'n_snp_clumped': len(clumped)}
        rec.update({'Q': '', 'Q_df': '', 'Q_p': '', 'I2_pct': '',
                    'Q_GX': '', 'I2_GX_pct': '', 'F_lead': '',
                    'boot_ci_lo': '', 'boot_ci_hi': '', 'E_value_point': '',
                    'E_value_ci': '', 'ABF_log': '', 'steiger_p_observed': '',
                    'steiger_p_conservative': ''})

        if len(clumped) >= 2:
            q = ivw_q(clumped)
            gx = i2_gx(clumped)
            bs = bootstrap_ivw(clumped)
            if q:
                rec.update({'Q': q['Q'], 'Q_df': q['Q_df'], 'Q_p': q['Q_p'],
                            'I2_pct': q['I2_pct'], 'beta_ivw': q['beta_ivw']})
            rec.update({'Q_GX': gx['Q_GX'], 'I2_GX_pct': gx['I2_GX_pct']})
            if bs:
                rec.update({'boot_ci_lo': bs['ci_lo'], 'boot_ci_hi': bs['ci_hi'],
                            'n_boot': bs['n_boot']})

        if lead:
            L = lead[0]
            zexp = L['z_exp']
            rec['F_lead'] = zexp ** 2
            wr_beta = L['beta_out'] / L['beta_exp']
            wr_se = abs(L['se_out'] / L['beta_exp'])
            orr = math.exp(wr_beta)
            ci_lo_or = math.exp(wr_beta - Z975 * wr_se)
            ci_hi_or = math.exp(wr_beta + Z975 * wr_se)
            rec['OR_primary'] = orr
            rec['OR_ci_lo'] = ci_lo_or
            rec['OR_ci_hi'] = ci_hi_or
            rec['E_value_point'] = e_value(orr)
            rec['E_value_ci'] = 1.0 if (ci_lo_or <= 1 <= ci_hi_or) else e_value(
                ci_lo_or if abs(math.log(ci_lo_or)) < abs(math.log(ci_hi_or)) else ci_hi_or)
            rec['ABF_log'] = log_abf(L['beta_out'], L['se_out'])
            rec['ABF'] = math.exp(rec['ABF_log'])

            # Steiger：暴露侧 R^2 与结局侧 R^2（三种 N_out 口径）
            r2e = r2_from_z(zexp, N_EXP)
            zo = L['beta_out'] / L['se_out']
            specs = [('observed_N_total', n_total, 'official manifest'),
                     ('observed_N_eff', n_eff, 'official manifest, effective N'),
                     ('conservative_2x_cases', n_conservative,
                      fb['source'] if (fb and not nc) else
                      ('official manifest n_cases' if nc else 'unavailable'))]
            for label, nout, src in specs:
                r2o = r2_from_z(zo, nout)
                s = steiger(r2e, r2o, N_EXP, nout)
                steiger_rows.append({
                    'release': rel, 'endpoint': ep, 'variant': LEAD,
                    'n_out_definition': label,
                    'n_out': int(nout) if nout else '',
                    'n_out_source': src,
                    'r2_exposure': r2e,
                    'r2_outcome': r2o if r2o is not None else '',
                    'z_stat': s['z_stat'] if s else '',
                    'p_two_sided': s['p_two_sided'] if s else '',
                    'direction_supports_exp_to_out':
                        s['direction_supports_exp_to_out'] if s else 'NOT_EVALUABLE',
                    'note': 'approximate t^2/R^2 identity; binary outcome on observed scale; not Hemani steiger implementation'
                    if s else 'no public case counts for this release; not imputed'})
                if s and label == 'observed_N_eff':
                    rec['steiger_p_observed'] = s['p_two_sided']
                if s and label == 'conservative_2x_cases':
                    rec['steiger_p_conservative'] = s['p_two_sided']
                    rec['r2_outcome_conservative'] = r2o

            # 11a：lead 位点三组关联表
            v = vmap.get(LEAD, {})
            lead_rows.append({
                'release': rel, 'endpoint': ep,
                'endpoint_type': meta.get('endpoint_type', ''),
                'variant': LEAD, 'rsid': v.get('rsids', L.get('rsids', '')),
                'cis': 'True', 'eaf': L['af_alt'],
                'exp_beta': L['beta_exp'], 'exp_se': L['se_exp'],
                'exp_z': L['z_exp'], 'exp_p': p_exp(v),
                'out_beta': L['beta_out'], 'out_se': L['se_out'],
                'out_p': L['pval_out'],
                'wald_beta': wr_beta, 'wald_se': wr_se,
                'wald_p': float(2 * norm.sf(abs(wr_beta / wr_se))),
                'wald_or': orr, 'wald_ci_lo_or': ci_lo_or,
                'wald_ci_hi_or': ci_hi_or,
                'F': zexp ** 2,
                'source': 'a1_outcome_matrix.tsv / a1_round1_variants.tsv'})

        stats_rows.append(rec)

    # 15：偏倚方向与量级（全部取自实测）
    from a1_multi_mr import power_at_or, min_detectable_or  # noqa: E402
    # 直接由 lead 位点实测 Wald ratio SE 计算，避免引用未定义字段
    pw, mde = [], []
    for r in lead_rows:
        pw.append(power_at_or(r['wald_se'], 1.5))
        mde.append(min_detectable_or(r['wald_se']))
    bias_rows = [
        {'bias': '不精确 / 低效能（主因）', 'direction': '偏向 null、并低估不确定性',
         'magnitude': 'OR=1.5 把握度 %.2f%%–%.2f%%；80%% 把握度 MDE OR %.2f–%.2f；主分析 95%%CI 的 OR 上界 %.2f'
                      % (min(pw) * 100, max(pw) * 100, min(mde), max(mde),
                         max(r['wald_ci_hi_or'] for r in lead_rows)),
         'mitigation': '报告事后效能、MDE 与 CI 宽度；不写"无效"'},
        {'bias': '水平多效性（trans 工具变量）', 'direction': '可制造伪阳性或方向相反的名义显著',
         'magnitude': '6 个 trans 位点在 Open Targets 预计算 coloc 中有 4–2,481 条共定位，H4 最高 1.0（ABO 位点）',
         'mitigation': '主分析仅用单一 cis lead；trans 只作敏感性并明确不作因果解释'},
        {'bias': '端点错分与对照定义漂移', 'direction': '可双向漂移，甚至翻转方向',
         'magnitude': '19 组合 OR 0.350–3.755（10.72 倍）；主结局 R9/R11 与 R12/R13 方向翻转；R9 对照 271,817 → R12 92,612',
         'mitigation': '多端点 × 多 release 全披露；端点定义对照表；不挑单一组合报告'},
        {'bias': '弱工具变量偏倚', 'direction': '若存在则偏向 null',
         'magnitude': 'lead cis F=%.2f >> 10；I²_GX 实测 %.2f%%–%.2f%%（NOME 满足）→ 不是本研究主要问题'
                      % (next(r['F'] for r in lead_rows),
                         min(float(r['I2_GX_pct']) for r in stats_rows if r['I2_GX_pct'] != ''),
                         max(float(r['I2_GX_pct']) for r in stats_rows if r['I2_GX_pct'] != '')),
         'mitigation': '正式稿不得写"弱工具导致阴性"；只写结局侧 SE 与暴露效应尺度'},
        {'bias': '样本重叠', 'direction': '若存在则偏向 confounded 估计',
         'magnitude': 'UKB-PPP（英国）与 FinnGen（芬兰）为互不重叠队列，重叠人数 0',
         'mitigation': '按 STROBE-MR 第 10 项明确写"无样本重叠"'},
        {'bias': '基因-环境等价 / 终身暴露', 'direction': '外推到药物干预时方向不可判',
         'magnitude': '循环 pQTL 代表终身遗传差异，不等于成年期短期 CX3CL1 干预；不等于脊髓局部配受体信号',
         'mitigation': 'Discussion 单列一段；因果语言仅在 IV 假设成立时使用'},
    ]

    write_tsv(RES + '/a1_a18_stats.tsv', stats_rows)
    write_tsv(RES + '/a1_a18_steiger.tsv', steiger_rows)
    write_tsv(RES + '/a1_a18_lead_variant_associations.tsv', lead_rows)
    write_tsv(RES + '/a1_a18_bias_magnitude.tsv', bias_rows)

    payload = {
        'type': 'A-level STROBE-MR clear-up: real recalculation from local public summary statistics',
        'generated': '2026-11-10',
        'inputs': ['a1_outcome_matrix.tsv', 'a1_round1_variants.tsv',
                   'a1_robustness.json', 'a1_endpoint_dict.tsv'],
        'formulae': {
            'F': 'z_exp^2',
            'I2': 'max(0,(Q-df)/Q)*100, Q = sum(w*(beta_out - beta_IVW*beta_exp)^2), w=1/se_out^2',
            'I2_GX': 'max(0,(Q_GX-(k-1))/Q_GX)*100, Q_GX = sum((beta_exp-mu)^2/se_exp^2)',
            'Steiger': 'R2 = z^2/(z^2+N-2); Fisher z of sqrt(R2); independent-sample z test; approximate, NOT official steiger package',
            'E_value': 'RR + sqrt(RR(RR-1)), OR<1 inverted; CI side = 1 when CI crosses 1',
            'bootstrap': 'instrument-level resampling only, B=%d, seed=%d; NOT GWAS sampling uncertainty' % (N_BOOT, SEED),
            'ABF': '0.5*log(1-W/(V+W)) + 0.5*z^2*W/(V+W), W=%.2f (same prior as coloc)' % PRIOR_W,
        },
        'n_combos': len(stats_rows),
        'stats': stats_rows,
        'steiger': steiger_rows,
        'lead_variant': lead_rows,
        'bias_magnitude': bias_rows,
        'honesty': [
            'Steiger 采用透明近似实现，并给出 N_out 三种口径；R13 无官方 manifest，仅用 2*n_cases 保守边界，不推测对照数。',
            'E-value 源自观察性研究的混杂敏感性框架，在 MR 中只能作为残余多效性的辅助量化，不能替代 IV 假设检验。',
            'bootstrap CI 仅反映工具变量构成的不确定性，不替代主分析 Wald ratio 置信区间。',
            'I²_GX 实测结果若与原清单预设不一致，一律以实测为准，见汇总说明。',
        ],
    }
    with io.open(RES + '/a1_a18_stats.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print('combos', len(stats_rows), 'steiger rows', len(steiger_rows),
          'lead rows', len(lead_rows))
    print('I2_GX range', min(float(r['I2_GX_pct']) for r in stats_rows if r['I2_GX_pct'] != ''),
          max(float(r['I2_GX_pct']) for r in stats_rows if r['I2_GX_pct'] != ''))
    print('I2 range', min(float(r['I2_pct']) for r in stats_rows if r['I2_pct'] != ''),
          max(float(r['I2_pct']) for r in stats_rows if r['I2_pct'] != ''))
    print('steiger conservative p range',
          min(float(r['p_two_sided']) for r in steiger_rows if r['p_two_sided'] != ''),
          max(float(r['p_two_sided']) for r in steiger_rows if r['p_two_sided'] != ''))
    print('power OR1.5 range %.2f%% - %.2f%%' % (min(pw) * 100, max(pw) * 100))
    print('MDE range', min(mde), max(mde))


if __name__ == '__main__':
    import sys
    sys.path.insert(0, ROOT + '/code')
    main()
