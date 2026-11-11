# -*- coding: utf-8 -*-
"""
B3 步骤 1a：EBI FTP 单连接限速（约 15 KB/s），改用并行 Range 分块下载。

实测（2026-11-11）：
  单连接      约 15 KB/s
  8 并行      约 134 KB/s
  24 并行     约 490 KB/s
说明限速接近"按连接"而非"按 IP"，故并行可接近线性扩展。

实现：每块写独立临时文件，全部完成后按序拼接，最后校验 md5。
不使用多线程写同一句柄，避免 Windows 下的并发 seek 竞态。
"""
import hashlib
import io
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

URL = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
       'GCST90012001-GCST90013000/GCST90012074/harmonised/'
       '33067605-GCST90012074-EFO_0009419.h.tsv.gz')
TOTAL = 844542746          # 官方 Content-Length
MD5_EXPECTED = 'c04aa821ffd5ef776eb3363b64407b55'
CHUNK = 2 * 1024 * 1024
# 32 并发时约 44% 的块返回 502 Bad Gateway。根因不是 EBI 拒绝，而是本机
# 出网统一走本地代理（http://127.0.0.1:57141），并发过高会把代理压垮。
# 取 16 并发 + 8 次退避重试；已成功的块由 cached 分支跳过。
WORKERS = 16

# 关键优化：文件按染色体单调排序（实测 138 MB 前缀覆盖 chr1→chr4），
# 且是普通 gzip 而非 BGZF（无 FEXTRA，不能建 tabix 索引），故只能顺序读到
# chr16。chr16 约在 65% 处，无需下载尾部 chr17-Y。按块号上限截断，
# 拼接出的是合法但截断的 gz，流式解压到目标区域即止。
MAX_BLOCK = 340          # 340 x 2 MB ≈ 680 MB，足以覆盖 chr16 及其后
 
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(DATA, '33067605-GCST90012074-EFO_0009419.h.tsv.gz')
PARTS = os.path.join(DATA, '_parts_gcst90012074')


def fetch_range(args):
    idx, start, end = args
    path = os.path.join(PARTS, '%05d.part' % idx)
    need = end - start + 1
    if os.path.exists(path) and os.path.getsize(path) == need:
        return idx, need, 'cached'
    last = None
    for attempt in range(8):
        try:
            req = urllib.request.Request(
                URL, headers={'User-Agent': 'A1-coloc/1.0',
                              'Range': 'bytes=%d-%d' % (start, end)})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) != need:
                raise IOError('short read %d/%d' % (len(data), need))
            with open(path, 'wb') as fh:
                fh.write(data)
            return idx, need, 'ok'
        except Exception as e:      # noqa: BLE001
            last = '%s: %s' % (type(e).__name__, str(e)[:80])
            time.sleep(2 + attempt * 3)
    return idx, 0, 'FAIL %s' % last


def main():
    os.makedirs(PARTS, exist_ok=True)
    jobs = []
    idx = 0
    for s in range(0, TOTAL, CHUNK):
        if idx > MAX_BLOCK:
            break
        e = min(s + CHUNK - 1, TOTAL - 1)
        jobs.append((idx, s, e))
        idx += 1
    complete = (jobs[-1][2] == TOTAL - 1)
    print('总块数 %d（%s），并发 %d'
          % (len(jobs), '完整文件' if complete else '部分文件：仅到块 %d'
             % MAX_BLOCK, WORKERS), flush=True)

    t0 = time.time()
    done = 0
    failed = []
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(fetch_range, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            jidx, n, st = f.result()
            if st.startswith('FAIL'):
                failed.append((jidx, st))
            done += n
            if i % 20 == 0 or i == len(jobs):
                el = time.time() - t0
                print('  %d/%d 块，%.1f MB，%.1f s，%.0f KB/s'
                      % (i, len(jobs), done / 1e6, el,
                         done / 1024 / max(el, 1e-9)), flush=True)

    if failed:
        print('失败块 %d 个：%s' % (len(failed), failed[:5]), flush=True)
        return 2

    print('拼接 ...', flush=True)
    h = hashlib.md5()
    with open(OUT, 'wb') as out:
        for jidx, _, _ in jobs:
            p = os.path.join(PARTS, '%05d.part' % jidx)
            with open(p, 'rb') as fh:
                while True:
                    b = fh.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
                    h.update(b)
    got = h.hexdigest()
    size = os.path.getsize(OUT)
    print('拼接完成 %d 字节，md5=%s' % (size, got), flush=True)
    print('总用时 %.1f 分钟' % ((time.time() - t0) / 60), flush=True)
    if complete:
        print('md5 与官方一致：%s（期望 %s）'
              % (got == MD5_EXPECTED, MD5_EXPECTED), flush=True)
        if got != MD5_EXPECTED:
            print('!! md5 不符，需重下')
            return 3
    else:
        print('部分文件（截断于块 %d，%d 字节）；不校验整文件 md5，'
              '按设计只需覆盖到 chr16' % (MAX_BLOCK, size), flush=True)
    # 清理分块
    for jidx, _, _ in jobs:
        try:
            os.remove(os.path.join(PARTS, '%05d.part' % jidx))
        except OSError:
            pass
    try:
        os.rmdir(PARTS)
    except OSError:
        pass
    return 0


if __name__ == '__main__':
    sys.exit(main())
