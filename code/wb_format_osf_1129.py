# -*- coding: utf-8 -*-
"""把 OSF 登记操作指引 Markdown 转为正式 .docx（用户固定版式）。

相对 wb_format_a18_1129.py 的两处扩展：
  1. 支持 #### 四级标题（字段名级），用黑体小四加粗；
  2. 支持 ENBLOCK: ... ENDENBLOCK 英文块，整块排为一个段落、内部保留软换行，
     避免被拆成十几个独立段落而破坏可读性。
其余版式不变：题名方正小标宋简体二号；正文仿宋_GB2312 三号、1.5 倍行距、
首行缩进 2 字；一级标题黑体；二级标题楷体_GB2312；三线表；A4；页脚居中页码。
"""
import io
import os
import re
import sys

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

ROOT = 'F:/博士/国自然科学基金/2025神自然'
SRC = ROOT + '/_markdown源/A1_OSF登记操作指引_20261111.md'
OUT = ROOT + '/A1_OSF回顾性登记操作指引_20261111.docx'
# 可选：python wb_format_osf_1129.py <源.md> <输出.docx>
if len(sys.argv) >= 3:
    SRC, OUT = sys.argv[1], sys.argv[2]

TITLE_FONT, TITLE_SIZE = '方正小标宋简体', 22      # 二号
BODY_FONT, BODY_SIZE = '仿宋_GB2312', 16          # 三号
H1_FONT, H1_SIZE = '黑体', 16
H2_FONT, H2_SIZE = '楷体_GB2312', 16
H3_FONT, H3_SIZE = '黑体', 12.5
TABLE_FONT, TABLE_SIZE = '仿宋_GB2312', 10.5       # 五号
EN_FONT, EN_SIZE = 'Times New Roman', 11


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

    for cell in table.rows[0].cells:
        tcpr = cell._tc.get_or_add_tcPr()
        tcb = OxmlElement('w:tcBorders')
        b = OxmlElement('w:bottom')
        b.set(qn('w:val'), 'single')
        b.set(qn('w:sz'), '8')
        b.set(qn('w:color'), '000000')
        tcb.append(b)
        tcpr.append(tcb)
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


def add_en_block(doc, block_lines):
    """英文可粘贴块：一个段落，内部软换行，左缩进与浅灰底纹便于识别。"""
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    pf = p.paragraph_format
    pf.line_spacing = 1.25
    pf.space_before = Pt(3)
    pf.space_after = Pt(8)
    pf.left_indent = Pt(18)
    pf.right_indent = Pt(6)
    ppr = p._p.get_or_add_pPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), 'F4F4F4')
    ppr.append(shd)

    first = True
    for ln in block_lines:
        if not first:
            r0 = p.add_run()
            set_run_font(r0, EN_FONT, EN_SIZE)
            r0.add_break()
        r = p.add_run(ln)
        set_run_font(r, EN_FONT, EN_SIZE)
        first = False
    return p


def main():
    lines = io.open(SRC, encoding='utf-8').read().splitlines()
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    sec.top_margin, sec.bottom_margin = Cm(3.0), Cm(2.5)
    sec.left_margin, sec.right_margin = Cm(2.5), Cm(2.5)
    add_page_number(doc)

    normal = doc.styles['Normal']
    normal.font.name = BODY_FONT
    normal.font.size = Pt(BODY_SIZE)
    normal._element.rPr.rFonts.set(qn('w:eastAsia'), BODY_FONT)

    i = 0
    n_blocks = 0
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue

        # 英文块
        if line.strip() == 'ENBLOCK:':
            i += 1
            buf = []
            while i < len(lines) and lines[i].strip() != 'ENDENBLOCK':
                buf.append(lines[i].rstrip())
                i += 1
            i += 1
            while buf and not buf[-1].strip():
                buf.pop()
            while buf and not buf[0].strip():
                buf.pop(0)
            add_en_block(doc, buf)
            n_blocks += 1
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
        elif line.startswith('#### '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(10)
            p.paragraph_format.space_after = Pt(3)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(line[5:].strip())
            set_run_font(r, H3_FONT, H3_SIZE, bold=True)
        elif line.startswith('## '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(14)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(line[3:].strip())
            set_run_font(r, H1_FONT, H1_SIZE, bold=True)
        elif line.startswith('### '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(10)
            p.paragraph_format.space_after = Pt(4)
            p.paragraph_format.keep_with_next = True
            r = p.add_run(line[4:].strip())
            set_run_font(r, H2_FONT, H2_SIZE, bold=True)
        elif line.startswith('EN: '):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 1.3
            p.paragraph_format.space_after = Pt(6)
            add_rich_text(p, line[4:].strip(), 'Times New Roman', 12)
        elif line.startswith('- '):
            p = doc.add_paragraph(style='List Bullet')
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.left_indent = Pt(32)
            add_rich_text(p, line[2:].strip(), BODY_FONT, BODY_SIZE)
        elif re.match(r'^\d+\.\s', line):
            p = doc.add_paragraph(style='List Number')
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.left_indent = Pt(32)
            add_rich_text(p, re.sub(r'^\d+\.\s', '', line), BODY_FONT, BODY_SIZE)
        else:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.first_line_indent = Pt(BODY_SIZE * 2)
            p.paragraph_format.space_after = Pt(2)
            add_rich_text(p, line.strip(), BODY_FONT, BODY_SIZE)
        i += 1

    doc.core_properties.title = 'A1 论文 OSF 回顾性登记操作指引'
    doc.core_properties.author = '鲁义主任医师'
    doc.save(OUT)
    print('saved', OUT, os.path.getsize(OUT), 'bytes')
    print('paragraphs', len(doc.paragraphs), 'tables', len(doc.tables),
          'en_blocks', n_blocks)


if __name__ == '__main__':
    main()
