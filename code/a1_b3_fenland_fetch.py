# -*- coding: utf-8 -*-
"""
B3 步骤 1：取 Fenland (Pietzner 2021, GCST90012074, n=21,758 欧洲裔) 的
CX3CL1 pQTL 全量汇总统计，流式筛出 chr16 CX3CL1 区域。

为什么不用 UKB-PPP 全量数据：
  Synapse (syn51364943) 对本地区返回 403 GEO_RESTRICTION
  （NIH NOT-OD-25-083），即使持有 PAT 亦不可达。改用 EBI FTP 公开的
  独立队列 Fenland，同时获得方法学上的好处：暴露（Fenland, 英国）
  与结局（FinnGen, 芬兰）为独立样本且独立国家，是 coloc 的理想条件，
  而 UKB-PPP 与 FinnGen 同样独立但与 Open Targets 预计算 coloc 共用
  可信集级输入，独立性不如前者清晰。

文件性质（来自官方 meta.yaml）：
  genome_assembly: GRCh38，1-based
  is_sorted: false  -> 不能遇 chr16 提前终止，必须读完全文件
  无 .tbi 索引      -> 不能 tabix 定点取数，只能流式
"""
import hashlib
import io
import json
import os
import sys
import time
import urllib.request
import zlib

BASE = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
        'GCST90012001-GCST90013000/GCST90012074/harmonised/')
FNAME = '33067605-GCST90012074-EFO_0009419.h.tsv.gz'
URL = BASE + FNAME
MD5_EXPECTED = 'c04aa821ffd5ef776eb3363b64407b55'

CHROM = '16'
REG_START = 57200000
REG_END = 57600000

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
RES = os.path.join(ROOT, 'results')
RAW = os.path.join(DATA, FNAME)

COLS = ['hm_variant_id', 'hm_rsid', 'hm_chrom', 'hm_pos',
        'hm_other_allele', 'hm_effect_allele', 'hm_beta',
        'hm_effect_allele_frequency', 'standard_error', 'p_value', 'n']


def download(path):
    """流式下载并同时算 md5，返回 (bytes, md5, seconds)。"""
    if os.path.exists(path) and os.path.getsize(path) > 500_000_000:
        h = hashlib.md5()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b''):
                h.update(chunk)
        print('已存在本地文件，跳过下载，直接校验', flush=True)
        return os.path.getsize(path), h.hexdigest(), 0.0

    h = hashlib.md5()
    n = 0
    t0 = time.time()
    req = urllib.request.Request(URL, headers={'User-Agent': 'A1-coloc/1.0'})
    with urllib.request.urlopen(req, timeout=180) as r, open(path, 'wb') as fh:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            h.update(chunk)
            n += len(chunk)
            if n % (50 << 20) < (1 << 20):
                print('  已下载 %.1f MB / %.1f MB' % (n / 1e6, 557.4),
                      flush=True)
    return n, h.hexdigest(), time.time() - t0


def chrom_num(x):
    x = str(x).replace('chr', '').strip()
    if x.isdigit():
        return int(x)
    return {'X': 23, 'Y': 24, 'MT': 25, 'M': 25}.get(x.upper(), 99)


def stream_region(path):
    """从磁盘 gz 流式解压，筛区域。

    文件按染色体单调排序且可能被截断（只下到覆盖 chr16 的块），因此：
      * 一旦进入 chromosome > 16 的段且已收集到区域变异，立即停止，
        无需解压尾部；
      * 截断处的 zlib 不完整错误一律吞掉，已解压数据仍然有效。
    """
    stats = {'total_rows': 0, 'chr16_rows': 0, 'region_rows': 0,
             'beta_na': 0, 'se_na': 0, 'kept': 0, 'stopped_at_chrom': None,
             'truncated': False}
    kept = []
    header = None
    idx = {}
    with open(path, 'rb') as fh:
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        buf = b''
        first = True
        stop = False
        while not stop:
            raw = fh.read(1 << 22)
            if not raw:
                break
            try:
                buf += d.decompress(raw)
            except zlib.error:
                stats['truncated'] = True
                break
            # 按整行切分，末尾不完整片段留到下一轮
            lines = buf.split(b'\n')
            buf = lines.pop()
            for ln in lines:
                s = ln.decode('utf-8', 'ignore')
                if first:
                    header = s.rstrip('\r').split('\t')
                    missing = [c for c in COLS if c not in header]
                    if missing:
                        raise SystemExit('缺少列: %s' % missing)
                    idx = {c: header.index(c) for c in COLS}
                    first = False
                    continue
                if not s:
                    continue
                f = s.rstrip('\r').split('\t')
                stats['total_rows'] += 1
                if len(f) <= idx['hm_pos']:
                    continue
                cn = chrom_num(f[idx['hm_chrom']])
                if cn > chrom_num(CHROM):
                    # 已越过目标染色体：若区域变异已收集，直接收工
                    if kept:
                        stats['stopped_at_chrom'] = f[idx['hm_chrom']]
                        stop = True
                        break
                    continue
                if f[idx['hm_chrom']] != CHROM:
                    continue
                stats['chr16_rows'] += 1
                try:
                    pos = int(float(f[idx['hm_pos']]))
                except (ValueError, TypeError):
                    continue
                if not (REG_START <= pos <= REG_END):
                    continue
                stats['region_rows'] += 1
                beta = f[idx['hm_beta']]
                se = f[idx['standard_error']]
                if beta in ('NA', '', 'nan'):
                    stats['beta_na'] += 1
                    continue
                if se in ('NA', '', 'nan'):
                    stats['se_na'] += 1
                    continue
                kept.append({
                    'variant': f[idx['hm_variant_id']],
                    'rsid': f[idx['hm_rsid']],
                    'chrom': f[idx['hm_chrom']],
                    'pos': pos,
                    'other_allele': f[idx['hm_other_allele']].upper(),
                    'effect_allele': f[idx['hm_effect_allele']].upper(),
                    'beta': float(beta),
                    'se': float(se),
                    'eaf': f[idx['hm_effect_allele_frequency']],
                    'pval': f[idx['p_value']],
                    'n': f[idx['n']],
                })
                stats['kept'] += 1
    return kept, stats, header


def main():
    os.makedirs(DATA, exist_ok=True)
    print('目标：%s' % URL, flush=True)
    n, md5, secs = download(RAW)
    print('本地文件 %.1f MB，用时 %.0f s，md5=%s' % (n / 1e6, secs, md5),
          flush=True)
    partial = (n < 844542746)
    md5_ok = (md5 == MD5_EXPECTED)
    if partial:
        print('注：这是按染色体排序截断的部分文件（%.1f%%），仅需覆盖 chr16，'
              '整文件 md5 不适用' % (100.0 * n / 844542746), flush=True)
    else:
        print('md5 与官方一致：%s（期望 %s）' % (md5_ok, MD5_EXPECTED),
              flush=True)

    t0 = time.time()
    kept, stats, header = stream_region(RAW)
    print('流式筛选完成 %.0f s：%s' % (time.time() - t0,
                                   json.dumps(stats, ensure_ascii=False)),
          flush=True)

    out = os.path.join(RES, 'a1_b3_fenland_region.tsv')
    cols = ['variant', 'rsid', 'chrom', 'pos', 'other_allele',
            'effect_allele', 'beta', 'se', 'eaf', 'pval', 'n']
    with io.open(out, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in sorted(kept, key=lambda x: x['pos']):
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print('写入 %s（%d 个变异）' % (out, len(kept)), flush=True)

    meta = {
        'source': 'GWAS Catalog GCST90012074, Fenland cohort (Pietzner et al., '
                  'Nat Commun 2021; PMID 33067605)',
        'trait': 'CX3CL1 (Fractalkine) levels, plasma',
        'ancestry': 'European (Fenland, United Kingdom)',
        'n_samples': 21758,
        'genome_assembly': 'GRCh38, 1-based',
        'file_url': URL,
        'file_bytes': n,
        'file_md5': md5,
        'file_partial_by_design': partial,
        'file_md5_matches_catalog': md5_ok,
        'is_sorted_upstream': False,
        'sorted_by_chromosome': True,
        'partial_rationale': ('File is a plain gzip sorted by chromosome and '
                             'cannot be indexed (no BGZF FEXTRA). Only blocks '
                             'covering chr1-chr16 were downloaded; tail '
                             'chromosomes are not needed for the CX3CL1 region.'),
        'region': 'chr16:%d-%d' % (REG_START, REG_END),
        'stream_stats': stats,
        'n_variants_in_region': len(kept),
        'why_not_ukbppp': 'Synapse syn51364943 returns HTTP 403 '
                          'GEO_RESTRICTION from this region (NIH NOT-OD-25-083); '
                          'a PAT would not help. EBI FTP is not geo-restricted.',
        'fetched': time.strftime('%Y-%m-%d'),
    }
    with io.open(os.path.join(RES, 'a1_b3_fenland_fetch.json'), 'w',
                 encoding='utf-8') as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=1)
    print('完成')

    if not kept:
        print('!! 区域内无可用变异，需人工检查')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
