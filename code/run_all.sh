#!/usr/bin/env bash
# A1 一键复现入口。默认只跑离线重算（不联网）；--with-fetch 才执行公开桶取数。
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python}"
export PYTHONPATH="$ROOT/code${PYTHONPATH:+:$PYTHONPATH}"

if [[ "${1:-}" == "--with-fetch" ]]; then
  echo "[1/14] 公开桶取数（需联网；FinnGen GCS + tabix HTTP Range）"
  "$PY" code/a1_multi_fetch.py
  echo "[2/14] LD / clumping / coloc / 技术阳性对照 / 多靶点队列（需联网）"
  "$PY" code/a1_ld_coloc.py
  "$PY" code/a1_pip_coloc.py
  "$PY" code/a1_positive_control.py
  "$PY" code/a1_b1b2_multitarget.py
  echo "[3/14] B3 附录：Fenland pQTL 区域取数（需联网；EBI FTP，约 500 MB）"
  "$PY" code/a1_b3_parallel_dl.py
  "$PY" code/a1_b3_fenland_fetch.py
  echo "[4/14] B3 附录：区域级完整共定位（需联网；FinnGen GCS）"
  "$PY" code/a1_b3_coloc_full.py
else
  echo "[1/14] 跳过联网取数，使用已落盘的真实汇总统计（加 --with-fetch 可重取）"
  echo "[2/14] 跳过联网 LD / coloc / 多靶点取数；沿用已落盘结果与 data/b1_target_pqtl.json"
  echo "[3/14] 跳过 B3 的 EBI 下载；沿用 results/a1_b3_fenland_region.tsv"
  echo "[4/14] 跳过 B3 的 FinnGen 取数；沿用 results/a1_b3_coloc_full.json"
fi

echo "[5/14] 19 组合 MR 矩阵、效能曲线、结论漂移"
"$PY" code/a1_multi_mr.py
echo "[6/14] 图 2—4（内部留档中文版）"
"$PY" code/a1_multi_plots.py
echo "[7/14] 图 1 区段关联对照图"
"$PY" code/a1_region_plot.py
echo "[8/14] STROBE-MR 20 项逐条对照表"
"$PY" code/a1_strobe_checklist.py
echo "[9/14] 决策流程仿真验证（图 5）"
"$PY" code/a1_simulation.py
echo "[10/14] A 级 12b / 13c / 11a 实算"
"$PY" code/a1_a18_stats.py
echo "[11/14] A4 漂移归因 + A5 四类工具变量集 76 估计"
"$PY" code/a1_a4a5_drift.py
echo "[12/14] 图 1—5 英文矢量投稿版"
"$PY" code/a1_figure_pack_en.py
echo "[13/14] 图 9 B3 区域级共定位英文矢量版"
"$PY" code/a1_b3_figure.py
echo "[14/14] 英文稿定稿骨架、投稿包、可复核包、SHA256、版本快照与 A 级 18 项清零对照表"
"$PY" code/a1_manuscript_pack.py
"$PY" code/a1_a18_clear.py
echo "完成。全部结果写入 results/；投稿包写入 submission_pack/。"
