# -*- coding: utf-8 -*-
"""生成 A1 策略报告（正式版式）与一页可勾选清单。

正式报告：题名方正小标宋简体二号；正文仿宋_GB2312 三号、1.5 倍行距、首行缩进 2 字；
一级标题黑体；二级标题楷体_GB2312；三线表、表头黑体；A4；页脚居中页码。
一页清单：无页码、紧凑三号题名 + 五号表格，可打印勾选。
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
SRC = ROOT + '/_markdown源/A1科学价值提升与投稿成功率策略_20261110.md'
OUT_REPORT = ROOT + '/A1科学价值提升与投稿接受率策略_20261110.docx'
OUT_CHECK = ROOT + '/A1科学价值提升与投稿勾选清单_20261110.docx'

TITLE_FONT, TITLE_SIZE = '方正小标宋简体', 22
BODY_FONT, BODY_SIZE = '仿宋_GB2312', 16
H1_FONT, H1_SIZE = '黑体', 16
H2_FONT, H2_SIZE = '楷体_GB2312', 16
TABLE_FONT, TABLE_SIZE = '仿宋_GB2312', 10.5
EN_FONT, EN_SIZE = 'Times New Roman', 12


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


def set_three_line_borders(table, header_font_size=None):
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


def add_table(doc, rows, size=TABLE_SIZE, line_spacing=1.15):
    ncols = max(len(r) for r in rows)
    t = doc.add_table(rows=len(rows), cols=ncols)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for i, row in enumerate(rows):
        for j in range(ncols):
            cell = t.cell(i, j)
            cell.text = ''
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(1)
            p.paragraph_format.space_after = Pt(1)
            p.paragraph_format.line_spacing = line_spacing
            txt = row[j] if j < len(row) else ''
            r = p.add_run(txt.replace('`', '').replace('**', ''))
            set_run_font(r, '黑体' if i == 0 else TABLE_FONT,
                         size, bold=(i == 0))
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


def build_report():
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
    while i < len(lines):
        line = lines[i].rstrip()
        if not line.strip():
            i += 1
            continue
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

    doc.core_properties.title = 'A1 论文科学价值提升与投稿接受率策略'
    doc.core_properties.author = '鲁义主任医师'
    doc.save(OUT_REPORT)
    print('saved report', OUT_REPORT, os.path.getsize(OUT_REPORT), 'bytes',
          'paragraphs', len(doc.paragraphs), 'tables', len(doc.tables))


CHECK_ITEMS = [
    ('A1', '证据三角显式化，声明 a priori 主要分析', '本期'),
    ('A2', '效能量化语句 + 结论四分类判定表', '本期'),
    ('A3', '按已发表评分体系逐项自评（Supplementary）', '本期'),
    ('A4', '端点漂移定量归因（病例数 ↔ 效应漂移）', '本期'),
    ('A5', 'IV 集维度 4×19=76，三类漂移统一框架', '本期'),
    ('A6', '审稿意见预判表 + cover letter 五项', '本期'),
    ('A7', '图矢量英文化、代码声明、CRediT、清单附件', '本期'),
    ('B1', '多靶点复现队列（n≥5），最高价值增益', '备胎'),
    ('B2', 'PCSK9/LDLR 生物学阳性对照', '备胎'),
    ('B3', 'PAT + 完整 coloc 附录', '备胎'),
    ('B4', '与 Hu/Zou/LocusCompare2 差异化段落', '本期'),
    ('C1', 'OSF 登记 + GitHub + Zenodo DOI 互链', '需本人'),
    ('C2', 'FinnGen 官网填表留回执 + 致谢', '需本人'),
    ('C3', '预印本 + 编辑部预投稿咨询', '需本人'),
    ('C4', '方法学共同作者 / 同行预读', '需本人'),
]


def build_checklist():
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    sec.top_margin, sec.bottom_margin = Cm(1.6), Cm(1.4)
    sec.left_margin, sec.right_margin = Cm(1.8), Cm(1.8)
    normal = doc.styles['Normal']
    normal.font.name = BODY_FONT
    normal.font.size = Pt(10.5)
    normal._element.rPr.rFonts.set(qn('w:eastAsia'), BODY_FONT)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(6)
    r = p.add_run('A1 科学价值提升与投稿勾选清单')
    set_run_font(r, TITLE_FONT, 18)

    for txt in (
        '2026-11-10｜通讯作者：鲁义主任医师｜本清单为策略核对用，未改动任何分析结果',
        '结论：① 瓶颈是稿件归类（纯 MR 结果论文易被直接拒稿）；② 最高杠杆是多靶点复现队列 B1；③ 首选 BMC Med Res Methodol 方法学栏目。',
        '口径：只写"未获得支持性证据、现有公开数据不足以检验中等效应"；不写"CX3CL1 无效"，不推翻标书假说。',
    ):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing = 1.05
        r = p.add_run(txt)
        set_run_font(r, BODY_FONT, 10.5)

    rows = [['□', '编号', '动作', '归属']]
    rows += [['□', k, v, g] for k, v, g in CHECK_ITEMS]
    t = doc.add_table(rows=len(rows), cols=4)
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    widths = [Cm(1.0), Cm(1.4), Cm(12.4), Cm(2.0)]
    for i, row in enumerate(rows):
        for j in range(4):
            cell = t.cell(i, j)
            cell.width = widths[j]
            cell.text = ''
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.0
            r = p.add_run(row[j])
            set_run_font(r, '黑体' if i == 0 else BODY_FONT,
                         10.5, bold=(i == 0))
    set_three_line_borders(t)

    for txt in (
        '红线：不挑 P<0.05 组合翻阳性；不用 bootstrap CI 报阳性；41% 须标"由 59% 推算"；11% 分母是 2021 年后文章；不拆稿。',
        '建议顺序：A1 → A2 → A3/A6 → A4/A5 → 英文定稿 + A7 → C1/C2/C3 → 投稿；B1 不阻塞首次投稿。',
    ):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(4)
        p.paragraph_format.space_after = Pt(0)
        p.paragraph_format.line_spacing = 1.05
        r = p.add_run(txt)
        set_run_font(r, BODY_FONT, 10)

    doc.core_properties.title = 'A1 科学价值提升与投稿勾选清单'
    doc.core_properties.author = '鲁义主任医师'
    doc.save(OUT_CHECK)
    print('saved checklist', OUT_CHECK, os.path.getsize(OUT_CHECK), 'bytes',
          'paragraphs', len(doc.paragraphs), 'tables', len(doc.tables))


if __name__ == '__main__':
    build_report()
    build_checklist()
