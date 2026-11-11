# -*- coding: utf-8 -*-
"""
A1 定点取数工具：通过 tabix (.tbi) 索引 + HTTP Range，仅下载 BGZF 文件中
指定基因组区间的块，避免下载整个 FinnGen GWAS 汇总统计（约 801 MB）。

用法：
  python a1_tabix_fetch.py --inspect
  python a1_tabix_fetch.py --region chr16:57300000-57500000 --out out.tsv

依赖：仅标准库（urllib、zlib、struct）。无需 R、pysam 或 htslib。
"""
import argparse
import gzip
import struct
import sys
import zlib
from urllib.request import Request, urlopen

BASE = ('https://storage.googleapis.com/finngen-public-data-r13/'
        'summary_stats/finngen_R13_DM_NEUROPATHY.gz')
DEFAULT_TBI = ('F:/博士/国自然科学基金/2025神自然/A1_启动/data/'
               'finngen_R13_DM_NEUROPATHY.gz.tbi')


def u32(b, o):
    return struct.unpack_from('<I', b, o)[0]


def i32(b, o):
    return struct.unpack_from('<i', b, o)[0]


def u64(b, o):
    return struct.unpack_from('<Q', b, o)[0]


class Tabix:
    def __init__(self, path):
        with gzip.open(path, 'rb') as f:
            self.b = f.read()
        b = self.b
        assert b[:4] == b'TBI\x01', 'not a tabix index: %r' % b[:4]
        o = 4
        self.n_ref = i32(b, o); o += 4
        self.format = i32(b, o); o += 4
        # htslib 新版 TBI 头：依次写出 sc, bc, ec, meta_char, line_skip, l_nm
        self.sc = i32(b, o); o += 4
        self.bc = i32(b, o); o += 4
        self.ec = i32(b, o); o += 4
        self.meta_char = i32(b, o); o += 4
        self.line_skip = i32(b, o); o += 4
        self.l_nm = i32(b, o); o += 4
        self.names_raw = b[o:o + self.l_nm]; o += self.l_nm
        self.refs = {}
        for _ in range(self.n_ref):
            n_bin = i32(b, o); o += 4
            bins = {}
            for _ in range(n_bin):
                bin_id = u32(b, o); o += 4
                n_chunk = i32(b, o); o += 4
                chunks = []
                for _ in range(n_chunk):
                    chunks.append((u64(b, o), u64(b, o + 8)))
                    o += 16
                bins[bin_id] = chunks
            n_intv = i32(b, o); o += 4
            ioff = [u64(b, o + 8 * k) for k in range(n_intv)]
            o += 8 * n_intv
            # names 是连续 null-terminated 字符串，按 ref 顺序分配
            self.refs.setdefault('_order', []).append((bins, ioff))
        names = [x for x in self.names_raw.split(b'\x00') if x]
        self.refs = {n.decode(): self.refs['_order'][i]
                     for i, n in enumerate(names)}
        self.n_no_coor = u64(b, o) if o + 8 <= len(b) else None

    def lower_bound_voff(self, chrom, start):
        """线性索引给出覆盖 start 的最小虚拟偏移。"""
        bins, ioff = self.refs[chrom]
        idx = max(0, (start - 1) >> 14)
        idx = min(idx, len(ioff) - 1)
        voff = ioff[idx]
        # 同时纳入 bin 索引中可能更早的候选块
        reg_bins = reg2bins(max(0, start - 1), start)
        cand = [voff]
        for bid in reg_bins:
            for cb, ce in bins.get(bid, []):
                cand.append(cb)
        return min(cand)


def reg2bins(beg, end):
    """htslib 的 reg2bins：返回与 [beg,end) 重叠的 bin id 列表。"""
    out = [0]
    end -= 1
    for shift, offset in ((14, 4681), (17, 585), (20, 73), (23, 9),
                          (26, 1)):
        lo = offset + (beg >> shift)
        hi = offset + (end >> shift)
        out.extend(range(lo, hi + 1))
    return out


def http_range(url, start, length):
    end = start + length - 1
    req = Request(url, headers={'Range': 'bytes=%d-%d' % (start, end),
                                'User-Agent': 'A1-MR/1.0'})
    with urlopen(req, timeout=90) as r:
        return r.read()


def bgzf_decompress_from(buf, start_uoffset=0):
    """从一段 BGZF 字节流中逐个解压 gzip member，返回拼接的明文。"""
    out = []
    pos = 0
    while pos + 2 <= len(buf) and buf[pos:pos + 2] == b'\x1f\x8b':
        try:
            d = zlib.decompressobj(31)
            chunk = d.decompress(buf[pos:])
        except zlib.error:
            break
        if not d.eof:
            break  # 末块不完整，忽略
        if not out and start_uoffset:
            chunk = chunk[start_uoffset:]
        out.append(chunk)
        consumed = len(buf[pos:]) - len(d.unused_data)
        if consumed <= 0:
            break
        pos += consumed
    return b''.join(out)


def fetch_region(chrom, start, end, url=BASE, tbi_path=DEFAULT_TBI,
                 span_mb=4, tbi=None):
    t = tbi or Tabix(tbi_path)
    # 不同 FinnGen release 的 tabix 参考序列命名可能是 "16" 或 "chr16"
    if chrom not in t.refs:
        c2 = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
        if c2 in t.refs:
            chrom = c2
        else:
            raise KeyError('chrom not in index: %s; available=%s'
                           % (chrom, list(t.refs)[:30]))
    voff = t.lower_bound_voff(chrom, start)
    coff = voff >> 16
    uoff = voff & 0xFFFF
    step = span_mb * 1024 * 1024
    lines = []
    cursor = coff
    seen_text = b''
    for _ in range(12):
        raw = http_range(url, cursor, step)
        if not raw:
            break
        txt = bgzf_decompress_from(raw, uoff if cursor == coff else 0)
        seen_text += txt
        # 判断是否已读过 end
        rows = txt.decode('utf-8', 'replace').splitlines()
        maxpos = 0
        for r in rows:
            if r.startswith('#') or not r.strip():
                continue
            f = r.split('\t')
            if len(f) < 2:
                f = r.split()
            try:
                p = int(f[1])
            except (ValueError, IndexError):
                continue
            maxpos = max(maxpos, p)
            if start <= p <= end:
                lines.append(r)
        if maxpos > end:
            break
        cursor += step
        if len(raw) < step:
            break
    return lines, seen_text[:2000]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inspect', action='store_true')
    ap.add_argument('--region', help='e.g. chr16:57300000-57500000')
    ap.add_argument('--out')
    ap.add_argument('--tbi', default=DEFAULT_TBI)
    a = ap.parse_args()
    t = Tabix(a.tbi)
    if a.inspect:
        print('n_ref=%d format=%d sc=%d bc=%d ec=%d meta=%d skip=%d l_nm=%d'
              % (t.n_ref, t.format, t.sc, t.bc, t.ec, t.meta_char,
                 t.line_skip, t.l_nm))
        print('refs(%d)=' % len(t.refs), list(t.refs)[:30])
        for c in list(t.refs)[:3]:
            bins, ioff = t.refs[c]
            print('  %s: bins=%d linear_bins=%d first_ioff=%d'
                  % (c, len(bins), len(ioff), ioff[0] if ioff else -1))
        return
    if not a.region:
        ap.error('--region required (or use --inspect)')
    chrom, rng = a.region.split(':')
    s, e = map(int, rng.split('-'))
    lines, head = fetch_region(chrom, s, e, tbi_path=a.tbi)
    print('region=%s:%d-%d rows=%d' % (chrom, s, e, len(lines)))
    if lines:
        print('HEADER/FIRST:', lines[0][:300])
        print('LAST:', lines[-1][:300])
    else:
        print('decompressed head:', head[:500])
    if a.out:
        with open(a.out, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        print('written', a.out)


if __name__ == '__main__':
    main()
