# -*- coding: utf-8 -*-
"""复用既有固定版式生成器，把 A1-A3+A6 补正包 Markdown 源转为正式 .docx。

版式：题名方正小标宋简体二号；正文仿宋_GB2312 三号、1.5 倍行距、首行缩进 2 字；
一级标题黑体；二级标题楷体_GB2312；三线表、表头黑体；A4；页脚居中页码；
英文段落 Times New Roman 小四。
"""
import importlib.util
import os

from docx import Document

ROOT = 'F:/博士/国自然科学基金/2025神自然'
FMT = ROOT + '/A1_启动/code/wb_format_a1_strategy_1129.py'
SRC = ROOT + '/_markdown源/A1_A1A2A3A6补正包_20261110.md'
OUT = ROOT + '/A1_A1A2A3A6补正包_20261110.docx'

spec = importlib.util.spec_from_file_location('wb_fmt', FMT)
fmt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fmt)

fmt.SRC = SRC
fmt.OUT_REPORT = OUT
fmt.build_report()

d = Document(OUT)
d.core_properties.title = 'A1 论文 A1-A3 + A6 补正包'
d.core_properties.author = '鲁义主任医师'
d.save(OUT)
print('final', OUT, os.path.getsize(OUT), 'bytes',
      'paragraphs', len(d.paragraphs), 'tables', len(d.tables))
