# -*- coding: utf-8 -*-
"""把审核修订后的英文稿定稿骨架输出为可编辑 Word（英文期刊式排版）。"""
import io
import re
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, Cm

SRC = ('F:/博士/国自然科学基金/2025神自然/_markdown源/'
       'A1英文稿定稿骨架_20261111.md')
OUT = ('F:/博士/国自然科学基金/2025神自然/'
       'A1英文稿定稿骨架_20261111.docx')

FONT = 'Times New Roman'


def run(p, text, size=12, bold=False, italic=False):
    r = p.add_run(text)
    r.font.name = FONT
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic
    return r


def main():
    lines = io.open(SRC, encoding='utf-8').read().splitlines()
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21), Cm(29.7)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    sec.left_margin = sec.right_margin = Cm(2.5)
    n = doc.styles['Normal']
    n.font.name = FONT
    n.font.size = Pt(12)

    for line in lines:
        s = line.rstrip()
        if not s.strip():
            continue
        if s.startswith('# '):
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_after = Pt(14)
            run(p, s[2:].strip(), 16, bold=True)
        elif s.startswith('## '):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(12)
            p.paragraph_format.space_after = Pt(6)
            p.paragraph_format.keep_with_next = True
            run(p, s[3:].strip(), 13, bold=True)
        elif s.startswith('> '):
            p = doc.add_paragraph()
            p.paragraph_format.space_after = Pt(8)
            run(p, s[2:].strip(), 10, italic=True)
        elif s.startswith('- ') or s.startswith('  S') or s.startswith('  I ') \
                or s.startswith('  II') or s.startswith('  III') \
                or s.startswith('  IV'):
            p = doc.add_paragraph(style='List Bullet')
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.3
            run(p, s.strip().lstrip('- ').strip(), 12)
        elif s.startswith('EN: '):
            p = doc.add_paragraph()
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.space_after = Pt(6)
            run(p, s[4:].strip(), 12)
        else:
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            p.paragraph_format.line_spacing = 1.5
            p.paragraph_format.space_after = Pt(6)
            run(p, s.strip(), 12)

    doc.core_properties.title = 'A1 consolidated English manuscript skeleton'
    doc.core_properties.author = 'Yi Lu'
    doc.save(OUT)
    print('saved', OUT, 'paragraphs', len(doc.paragraphs))


if __name__ == '__main__':
    main()
