# -*- coding: utf-8 -*-
"""按用户固定公文版式生成 OSF 归档期全文审核报告 .docx。"""
from docx import Document
import wb_format_final_1129 as F

F.SRC = 'F:/博士/国自然科学基金/2025神自然/_markdown源/A1_OSF归档期全文审核报告_20261111.md'
F.OUT = 'F:/博士/国自然科学基金/2025神自然/A1_OSF归档期全文审核报告_20261111.docx'
F.main()

d = Document(F.OUT)
d.core_properties.title = 'A1 OSF 归档期全文审核报告（2026-11-11）'
d.core_properties.author = '鲁义主任医师'
d.save(F.OUT)
print('metadata fixed:', F.OUT)
