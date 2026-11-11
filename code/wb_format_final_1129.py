# -*- coding: utf-8 -*-
"""把 A1 A 级 18 项清零 Markdown 源转为正式 .docx（用户固定版式）。

版式：题名方正小标宋简体二号；正文仿宋_GB2312 三号、1.5 倍行距、首行缩进 2 字；
一级标题黑体；二级标题楷体_GB2312；表格三线表、表头黑体；A4；页脚居中页码。
英文稿段落用 Times New Roman 小四，与中文正文区分但不破坏整体版式。
"""
import io
import os
import re

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

ROOT = 'F:/博士/国自然科学基金/2025神自然'
SRC = ROOT + '/_markdown源/A1剩余任务提前清零报告_20261111.md'
OUT = ROOT + '/A1剩余任务提前清零报告_20261111.docx'

TITLE_FONT, TITLE_SIZE = '方正小标宋简体', 22      # 二号
BODY_FONT, BODY_SIZE = '仿宋_GB2312', 16          # 三号
H1_FONT, H1_SIZE = '黑体', 16
H2_FONT, H2_SIZE = '楷体_GB2312', 16
TABLE_FONT, TABLE_SIZE = '仿宋_GB2312', 10.5       # 五号
EN_FONT, EN_SIZE = 'Times New Roman', 12           # 小四


def set_run_font(run, font, size, bold=False, eastasia=None):
    run.font.name = font
    run.font.size = Pt(size)
    run.font.bold = bold
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.rFonts
    if rfonts is None:
        rfonts = OxmlElement('w:rFonts')
        rpr.insert(0, rfonts)
    rfonts.set(qn('w:eastAsia'), eastasia or font)
    rfonts.set(qn('w:ascii'), font)
    rfonts.set(qn('w:hAnsi'), font)


def add_rich_text(p, text, font, size, bold=False, eastasia=None):
    """处理 **加粗** 与 `行内代码`，其余按原样输出。"""
    text = text.replace('`', '')
    parts = re.split(r'(\*\*[^*]+\*\*)', text)
    for part in parts:
        if not part:
            continue
        is_bold = part.startswith('**') and part.endswith('**')
        clean = part[2:-2] if is_bold else part
        r = p.add_run(clean)
        set_run_font(r, font, size, bold or is_bold, eastasia)


def set_three_line_borders(table):
    tblpr = table._tbl.tblPr
    old = tblpr.find(qn('w:tblBorders'))
    if old is not None:
        tblpr.remove(old)
    borders = OxmlElement('w:tblBorders')
    for edge in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        el = OxmlElement('w:' + edge)
        if edge in ('top', 'bottom'):
            el.set(qn('w:val'), 'single')
            el.set(qn('w:sz'), '8')
            el.set(qn('w:color'), '000000')
        else:
            el.set(qn('w:val'), 'nil')
        borders.append(el)
    tblpr.append(borders)

    # 表头行下横线
    for cell in table.rows[0].cells:
        tcpr = cell._tc.get_or_add_tcPr()
        tcb = OxmlElement('w:tcBorders')
        b = OxmlElement('w:bottom')
        b.set(qn('w:val'), 'single')
        b.set(qn('w:sz'), '8')
        b.set(qn('w:color'), '000000')
        tcb.append(b)
        tcpr.append(tcb)
    # 表头跨页重复
    trpr = table.rows[0]._tr.get_or_add_trPr()
    th = OxmlElement('w:tblHeader')
    th.set(qn('w:val'), 'true')
    trpr.append(th)


def add_table(doc, rows):
    ncols = max(len(r) for r in rows)
    t = doc.add_table(rows=len(rows), cols=ncols)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for i, row in enumerate(rows):
        for j in range(ncols):
            cell = t.cell(i, j)
            cell.text = ''
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.15
            txt = row[j] if j < len(row) else ''
            r = p.add_run(txt.replace('`', '').replace('**', ''))
            set_run_font(r, '黑体' if i == 0 else TABLE_FONT,
                         TABLE_SIZE, bold=(i == 0))
    set_three_line_borders(t)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)
    return t


def add_page_number(doc):
    section = doc.sections[0]
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run()
    set_run_font(r, BODY_FONT, 10.5)
    fld = OxmlElement('w:fldSimple')
    fld.set(qn('w:instr'), 'PAGE')
    r._element.addnext(fld)


def is_table_sep(line):
    return bool(re.match(r'^\s*\|[\s:|-]+\|\s*$', line))


def split_row(line):
    return [c.strip() for c in line.strip().strip('|').split('|')]


def main():
    lines = io.open(SRC, encoding='utf-8').read().splitlines()
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    sec.top_margin, sec.bottom_margin = Cm(3.0), Cm(2.5)
    sec.left_margin, sec.right_margin = Cm(2.5), Cm(2.5)
    add_page_number(doc)

    # 默认正文样式
    normal = doc.styles['Normal']
    normal.font.name = BODY_FONT
    normal.font.size = Pt(BODY_SIZE)
    normal._element.rPr.rFonts.set(qn('w:eastAsia'), BODY_FONT)

    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue

        # 表格块
        if line.lstrip().startswith('|') and i + 1 < len(lines) and is_table_sep(lines[i + 1]):
            rows = [split_row(line)]
            i += 2
            while i < len(lines) and lines[i].lstrip().startswith('|'):
                rows.append(split_row(lines[i]))
                i += 1
            add_table(doc, rows)
            continue

        if line.startswith('# '):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_after = Pt(18)
            p.paragraph_format.line_spacing = 1.5
            r = p.add_run(line[2:].strip())
            set_run_font(r, TITLE_FONT, TITLE_SIZE)
        elif line.startswith('## '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(line[3:].strip())
            set_run_font(r, H1_FONT, H1_SIZE, bold=True)
        elif line.startswith('### '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(line[4:].strip())
            set_run_font(r, H2_FONT, H2_SIZE, bold=True)
        elif line.startswith('EN: '):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 1.3
            p.paragraph_format.space_after = Pt(6)
            add_rich_text(p, line[4:].strip(), EN_FONT, EN_SIZE)
        elif line.startswith('- '):
            p = doc.add_paragraph(style='List Bullet')
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.left_indent = Pt(32)
            add_rich_text(p, line[2:].strip(), BODY_FONT, BODY_SIZE)
        else:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.first_line_indent = Pt(BODY_SIZE * 2)
            p.paragraph_format.space_after = Pt(2)
            add_rich_text(p, line.strip(), BODY_FONT, BODY_SIZE)
        i += 1

    doc.core_properties.title = 'A1 剩余任务提前清零报告（2026-11-11）'
    doc.core_properties.author = '鲁义主任医师'
    doc.save(OUT)
    print('saved', OUT, os.path.getsize(OUT), 'bytes')
    print('paragraphs', len(doc.paragraphs), 'tables', len(doc.tables))


if __name__ == '__main__':
    main()
