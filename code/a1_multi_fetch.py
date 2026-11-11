# -*- coding: utf-8 -*-
"""
A1 方法学主线：跨 release × 跨端点批量定点取数。

对 (release, endpoint) 组合，仅通过 tabix 索引 + HTTP Range 抓取 9 个
UKB-PPP CX3CL1 pQTL 工具变量所在区间的真实 FinnGen 汇总统计，做等位基因
harmonise，产出可复核的长表。不下载约 800 MB/端点 的全文件，不做任何模拟。

输出：
  results/a1_endpoint_dict.tsv   端点字典（病例/对照数来自官方 manifest）
  results/a1_outcome_matrix.tsv  每个组合 × 每个变异的真实结局效应
  results/a1_fetch_log.json      取数与匹配日志（含失败原因）
  data/tbi/*.tbi                 缓存的 tabix 索引
"""
import json
import math
import os
import ssl
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_tabix_fetch import Tabix, fetch_region  # noqa: E402

ROOT = 'F:/博士/国自然科学基金/2025神自然/A1_启动'
DATA = ROOT + '/data'
RES = ROOT + '/results'
TBI_DIR = DATA + '/tbi'
CS_JSON = DATA + '/ukbppp_cx3cl1_crediblesets.json'

COL = {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'rsids': 4,
       'gene': 5, 'pval': 6, 'mlogp': 7, 'beta': 8, 'se': 9, 'af': 10}

# 9 个工具变量所在的 7 个区间（合并同一染色体的邻近位点，减少请求数）
REGIONS = [
    ('16', 57374000, 57381000),      # 3 个 cis: 57374992 / 57375150 / 57379748
    ('1',  21345000, 21346000),      # 21345480
    ('3',  39265000, 39266500),      # 39265671
    ('9',  133272000, 133275000),    # 133273983
    ('9',  76489400, 76490500),      # 76489957
    ('10', 63388800, 63389800),      # 63389329
    ('19', 48703200, 48704200),      # 48703728
]

# (release, endpoint) 清单。R9/R10/R12 的病例对照数取自各自官方 manifest；
# R11/R13 无公开 manifest，标注为 None，效能分析一律使用实测 SE，不依赖病例数假设。
COMBOS = [
    ('R9',  'DM_NEUROPATHY'),
    ('R9',  'DM_POLYNEURO'),
    ('R9',  'G6_POLYOTHUNS'),
    ('R9',  'DM_RETINOPATHY_EXMORE'),
    ('R9',  'DM_NEPHROPATHY_EXMORE'),
    ('R10', 'G6_POLYOTHUNS'),
    ('R11', 'DM_NEUROPATHY'),
    ('R11', 'G6_POLYOTHUNS'),
    ('R12', 'DM_NEUROPATHY'),
    ('R12', 'G6_POLYOTHUNS'),
    ('R12', 'DM_RETINOPATHY_STRICT'),
    ('R12', 'DM_NEPHROPATHY_EXMORE'),
    ('R13', 'DM_NEUROPATHY'),
    ('R13', 'G6_POLYOTHUNS'),
    ('R13', 'DM_RETINOPATHY'),
    ('R13', 'DM_RETINOPATHY_STRICT'),
    ('R13', 'DM_RETINOPATHY_EXMORE'),
    ('R13', 'DM_NEPHROPATHY'),
    ('R13', 'DM_NEPHROPATHY_EXMORE'),
]

# 端点类型：primary=主结局，proxy=痛性代理，broad=泛化神经病变，
# negative_ctrl=糖尿病微血管并发症（阴性对照），strictness=定义严格度对照
def etype(ep):
    if ep == 'DM_NEUROPATHY':
        return 'primary'
    if ep == 'DM_POLYNEURO':
        return 'proxy_painful_DPN'
    if ep == 'G6_POLYOTHUNS':
        return 'broad_polyneuropathy'
    if 'RETINOPATHY' in ep or 'NEPHROPATHY' in ep:
        return 'negative_ctrl_microvascular'
    return 'other'


def gcs_url(rel, ep):
    b = 'finngen-public-data-' + rel.lower()
    if rel == 'R12':
        p = 'summary_stats/release/finngen_R12_%s.gz' % ep
    else:
        p = 'summary_stats/finngen_%s_%s.gz' % (rel, ep)
    return 'https://storage.googleapis.com/%s/%s' % (b, p)


def load_manifests():
    """从三个官方 manifest 提取端点病例/对照数与路径来源。"""
    out = {}
    for rel, fn in [('R9', 'r9_manifest.tsv'), ('R10', 'r10_manifest.tsv'),
                    ('R12', 'r12_manifest.tsv')]:
        path = os.path.join(DATA, fn)
        if not os.path.exists(path):
            continue
        with open(path, encoding='utf-8') as f:
            head = f.readline().rstrip('\n').split('\t')
            for line in f:
                c = line.rstrip('\n').split('\t')
                if len(c) < 5:
                    continue
                rec = dict(zip(head, c))
                out[(rel, c[0])] = {
                    'n_cases': int(c[3]), 'n_controls': int(c[4]),
                    'name': rec.get('name') or rec.get('phenotype') or '',
                    'source': fn,
                }
    return out


def get_tbi(rel, ep):
    """下载并缓存 tabix 索引，返回 (本地路径, Tabix 对象)。"""
    os.makedirs(TBI_DIR, exist_ok=True)
    local = os.path.join(TBI_DIR, 'finngen_%s_%s.gz.tbi' % (rel, ep))
    if not os.path.exists(local):
        url = gcs_url(rel, ep) + '.tbi'
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'A1-MR/1.0'})
                with urllib.request.urlopen(req, timeout=90, context=ctx) as r:
                    b = r.read()
                with open(local, 'wb') as f:
                    f.write(b)
                break
            except Exception as e:
                if attempt == 2:
                    raise RuntimeError('tbi download failed %s %s: %s' % (rel, ep, e))
                time.sleep(1)
    return local, Tabix(local)


def fetch_one_region(args):
    """单个 (combo, region) 取数；返回原始行列表。"""
    rel, ep, chrom, s, e, url, tbi = args
    last = None
    for attempt in range(3):
        try:
            lines, _ = fetch_region(chrom, s, e, url=url, tbi=tbi, span_mb=1)
            return rel, ep, chrom, lines, None
        except Exception as ex:
            last = str(ex)
            time.sleep(1)
    return rel, ep, chrom, [], last


def harmonise(cs, lines_by_chrom):
    """在已抓取的区间行中，为 9 个工具变量做等位基因 harmonise。"""
    out = []
    for r in cs:
        vid = r['variant']['id']
        parts = vid.split('_')
        chrom, pos = r['chromosome'], int(r['position'])
        e_ref, e_alt = parts[2], parts[3]
        exact, approx = None, None
        for line in lines_by_chrom.get(chrom, []):
            f = line.split('\t')
            if len(f) < 11:
                continue
            try:
                p = int(f[COL['pos']])
            except (ValueError, IndexError):
                continue
            if abs(p - pos) > 2:
                continue
            o_ref, o_alt = f[COL['ref']], f[COL['alt']]
            try:
                beta = float(f[COL['beta']])
                se = float(f[COL['se']])
                pval = float(f[COL['pval']])
                mlogp = float(f[COL['mlogp']])
            except (ValueError, IndexError):
                continue
            # 列定义交叉验证：mlogp 应约等于 -log10(pval)
            col_ok = (pval > 0 and abs(mlogp + math.log10(pval)) < 0.05)
            if (o_ref, o_alt) == (e_ref, e_alt):
                flip = False
            elif (o_ref, o_alt) == (e_alt, e_ref):
                flip = True
                beta = -beta
            else:
                continue
            rec = {'variant': vid, 'chrom': chrom, 'pos': pos,
                   'cis': not r['isTransQtl'],
                   'exp_ref_alt': e_ref + '/' + e_alt,
                   'out_ref_alt': o_ref + '/' + o_alt, 'flip': flip,
                   'beta_exp': r['beta'],
                   'se_exp': abs(r['beta']) / abs(r['zScore']),
                   'z_exp': r['zScore'],
                   'beta_out': beta, 'se_out': se, 'pval_out': pval,
                   'af_alt': float(f[COL['af']]),
                   'rsids': f[COL['rsids']], 'nearest_gene': f[COL['gene']],
                   'pos_offset': p - pos, 'approx': p != pos,
                   'column_check_ok': col_ok}
            if p == pos:
                exact = rec
            elif approx is None:
                approx = rec
        rec = exact or approx
        if rec:
            out.append(rec)
    return out


def main():
    os.makedirs(RES, exist_ok=True)
    cs = json.load(open(CS_JSON, encoding='utf-8'))['data']['credibleSets']['rows']
    man = load_manifests()

    # 端点字典
    edict = []
    for rel, ep in COMBOS:
        m = man.get((rel, ep))
        edict.append({
            'release': rel, 'endpoint': ep, 'endpoint_type': etype(ep),
            'n_cases': m['n_cases'] if m else '',
            'n_controls': m['n_controls'] if m else '',
            'endpoint_name': m['name'] if m else '',
            'case_control_source': m['source'] if m else 'NO_PUBLIC_MANIFEST',
            'url': gcs_url(rel, ep),
        })
    with open(RES + '/a1_endpoint_dict.tsv', 'w', encoding='utf-8') as f:
        keys = ['release', 'endpoint', 'endpoint_type', 'n_cases',
                'n_controls', 'endpoint_name', 'case_control_source', 'url']
        f.write('\t'.join(keys) + '\n')
        for d in edict:
            f.write('\t'.join(str(d[k]) for k in keys) + '\n')

    print('Loading %d tabix indexes ...' % len(COMBOS))
    tbis = {}
    for rel, ep in COMBOS:
        try:
            _, t = get_tbi(rel, ep)
            tbis[(rel, ep)] = t
            print('  OK  %s %-26s refs=%s' % (rel, ep, list(t.refs)[:2]))
        except Exception as e:
            print('  FAIL %s %-26s %s' % (rel, ep, e))

    # 构建 (combo, region) 任务
    jobs = []
    for rel, ep in COMBOS:
        if (rel, ep) not in tbis:
            continue
        url = gcs_url(rel, ep)
        for chrom, s, e in REGIONS:
            jobs.append((rel, ep, chrom, s, e, url, tbis[(rel, ep)]))
    print('Fetching %d region blocks with 5 threads ...' % len(jobs))

    collected = {}
    errors = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=5) as ex:
        futs = [ex.submit(fetch_one_region, j) for j in jobs]
        done = 0
        for fu in as_completed(futs):
            rel, ep, chrom, lines, err = fu.result()
            key = (rel, ep)
            collected.setdefault(key, {})[chrom] = (
                collected.get(key, {}).get(chrom, []) + lines)
            if err:
                errors.append({'release': rel, 'endpoint': ep,
                               'chrom': chrom, 'error': err})
            done += 1
            if done % 25 == 0:
                print('  %d/%d  %.0fs' % (done, len(jobs), time.time() - t0))

    # harmonise 并输出长表
    matrix = []
    log = {'n_combos': len(COMBOS), 'n_region_jobs': len(jobs),
           'elapsed_sec': round(time.time() - t0, 1),
           'errors': errors, 'combos': []}
    for d in edict:
        rel, ep = d['release'], d['endpoint']
        key = (rel, ep)
        if key not in collected:
            log['combos'].append({'release': rel, 'endpoint': ep,
                                  'status': 'NO_INDEX', 'n_matched': 0})
            continue
        lbc = collected[key]
        rows = harmonise(cs, lbc)
        for r in rows:
            r = dict(r)
            r['release'] = rel
            r['endpoint'] = ep
            r['endpoint_type'] = d['endpoint_type']
            r['n_cases'] = d['n_cases']
            r['n_controls'] = d['n_controls']
            matrix.append(r)
        n_lines = sum(len(v) for v in lbc.values())
        log['combos'].append({
            'release': rel, 'endpoint': ep, 'status': 'OK',
            'n_raw_lines': n_lines,
            'n_matched': len(rows),
            'n_approx': sum(1 for r in rows if r['approx']),
            'n_column_check_failed': sum(1 for r in rows
                                         if not r['column_check_ok']),
            'missing_variants': [c['variant']['id'] for c in cs
                                 if c['variant']['id'] not in
                                 {r['variant'] for r in rows}],
        })
        print('  %s %-26s raw=%5d matched=%d/9' %
              (rel, ep, n_lines, len(rows)))

    keys = ['release', 'endpoint', 'endpoint_type', 'n_cases', 'n_controls',
            'variant', 'rsids', 'chrom', 'pos', 'cis', 'exp_ref_alt',
            'out_ref_alt', 'flip', 'beta_exp', 'se_exp', 'z_exp',
            'beta_out', 'se_out', 'pval_out', 'af_alt', 'nearest_gene',
            'pos_offset', 'approx', 'column_check_ok']
    with open(RES + '/a1_outcome_matrix.tsv', 'w', encoding='utf-8') as f:
        f.write('\t'.join(keys) + '\n')
        for r in matrix:
            f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
    json.dump(log, open(RES + '/a1_fetch_log.json', 'w', encoding='utf-8'),
              ensure_ascii=False, indent=2)
    print('\nSaved matrix rows=%d  errors=%d  elapsed=%.1fs' %
          (len(matrix), len(errors), time.time() - t0))
    print('Outputs:', RES + '/a1_outcome_matrix.tsv',
          RES + '/a1_endpoint_dict.tsv', RES + '/a1_fetch_log.json')


if __name__ == '__main__':
    main()
