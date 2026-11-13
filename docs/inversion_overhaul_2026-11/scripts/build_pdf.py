"""Build a PDF from a Markdown report (headings, paragraphs, lists, tables, figures).

A small converter on reportlab, enough for REPORT.md: no pandoc or LaTeX needed.

Run: python build_pdf.py REPORT.md OUT.pdf
"""

import re
import sys
from pathlib import Path

from PIL import Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image as RLImage
from reportlab.platypus import (
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

WIDTH = A4[0] - 3.6 * cm
INK, MUTED, RULE = (
    colors.HexColor("#0b0b0b"),
    colors.HexColor("#52514e"),
    colors.HexColor("#c8c7c2"),
)


def styles():
    base = getSampleStyleSheet()
    body = ParagraphStyle(
        "body",
        parent=base["BodyText"],
        fontName="Helvetica",
        fontSize=9.5,
        leading=13.2,
        textColor=INK,
        alignment=TA_LEFT,
        spaceAfter=5,
    )
    return dict(
        body=body,
        h1=ParagraphStyle(
            "h1",
            parent=body,
            fontName="Helvetica-Bold",
            fontSize=17,
            leading=21,
            spaceBefore=2,
            spaceAfter=9,
        ),
        h2=ParagraphStyle(
            "h2",
            parent=body,
            fontName="Helvetica-Bold",
            fontSize=13,
            leading=16,
            spaceBefore=11,
            spaceAfter=5,
        ),
        h3=ParagraphStyle(
            "h3",
            parent=body,
            fontName="Helvetica-Bold",
            fontSize=10.5,
            leading=13.5,
            spaceBefore=8,
            spaceAfter=3,
        ),
        bullet=ParagraphStyle(
            "bullet", parent=body, leftIndent=13, bulletIndent=3, spaceAfter=2.5
        ),
        cell=ParagraphStyle(
            "cell", parent=body, fontSize=7.8, leading=9.8, spaceAfter=0
        ),
        head=ParagraphStyle(
            "head",
            parent=body,
            fontName="Helvetica-Bold",
            fontSize=7.8,
            leading=9.8,
            spaceAfter=0,
        ),
        caption=ParagraphStyle(
            "caption",
            parent=body,
            fontSize=8.3,
            leading=11,
            textColor=MUTED,
            spaceAfter=8,
        ),
        code=ParagraphStyle(
            "code",
            parent=body,
            fontName="Courier",
            fontSize=7.8,
            leading=10,
            leftIndent=8,
            backColor=colors.HexColor("#f1f0ec"),
            spaceAfter=6,
        ),
    )


def inline(text):
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"!\[[^\]]*\]\([^)]*\)", "", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)
    text = re.sub(
        r"`([^`]+)`",
        lambda m: (
            '<font face="Courier" size="8.3">'
            + m.group(1).replace("*", "&#42;").replace("_", "&#95;")
            + "</font>"
        ),
        text,
    )
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])", r"<i>\1</i>", text)
    return text


def table(rows, st):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [c for i, c in enumerate(cells) if i != 1]  # drop the separator row
    ncol = max(len(c) for c in cells)
    cells = [c + [""] * (ncol - len(c)) for c in cells]
    weight = [max(9, min(60, max(len(r[j]) for r in cells))) for j in range(ncol)]
    widths = [WIDTH * w / sum(weight) for w in weight]
    data = [
        [Paragraph(inline(c), st["head"] if i == 0 else st["cell"]) for c in r]
        for i, r in enumerate(cells)
    ]
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(
        TableStyle(
            [
                ("LINEBELOW", (0, 0), (-1, 0), 0.8, INK),
                ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
            ]
        )
    )
    return t


def build(source, target):
    source = Path(source)
    st = styles()
    story, para, lines = [], [], source.read_text().splitlines()

    def flush():
        if para:
            story.append(Paragraph(inline(" ".join(para)), st["body"]))
            para.clear()

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            flush()
            block = []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                block.append(
                    lines[i]
                    .replace("&", "&amp;")
                    .replace("<", "&lt;")
                    .replace(">", "&gt;")
                    .replace(" ", "&nbsp;")
                )
                i += 1
            story.append(Paragraph("<br/>".join(block), st["code"]))
        elif line.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(lines[i])
                i += 1
            story += [table(rows, st), Spacer(1, 7)]
            continue
        elif m := re.match(r"!\[[^\]]*\]\(([^)]+)\)", line.strip()):
            flush()
            path = source.parent / m.group(1)
            w, h = Image.open(path).size
            width = min(WIDTH, 17.4 * cm)
            img = RLImage(str(path), width=width, height=width * h / w)
            caption = None
            if (
                i + 2 < len(lines)
                and lines[i + 1].strip() == ""
                and lines[i + 2].startswith("*")
            ):
                caption = Paragraph(
                    inline(lines[i + 2].strip().strip("*")), st["caption"]
                )
                i += 2
            story.append(
                KeepTogether([img, Spacer(1, 3)] + ([caption] if caption else []))
            )
        elif m := re.match(r"(#{1,3}) (.*)", line):
            flush()
            story.append(Paragraph(inline(m.group(2)), st[f"h{len(m.group(1))}"]))
        elif m := re.match(r"\s*(-|\d+\.) (.*)", line):
            flush()
            item = [m.group(2)]
            while (
                i + 1 < len(lines)
                and lines[i + 1].startswith("  ")
                and lines[i + 1].strip()
                and not re.match(r"\s*(-|\d+\.) ", lines[i + 1])
            ):
                i += 1
                item.append(lines[i].strip())
            mark = "•" if m.group(1) == "-" else m.group(1)
            story.append(
                Paragraph(inline(" ".join(item)), st["bullet"], bulletText=mark)
            )
        elif not line.strip():
            flush()
        else:
            para.append(line.strip())
        i += 1
    flush()

    def footer(canvas, doc):
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(1.8 * cm, 1.1 * cm, source.stem.replace("_", " "))
        canvas.drawRightString(A4[0] - 1.8 * cm, 1.1 * cm, f"{doc.page}")

    doc = SimpleDocTemplate(
        str(target),
        pagesize=A4,
        leftMargin=1.8 * cm,
        rightMargin=1.8 * cm,
        topMargin=1.7 * cm,
        bottomMargin=1.8 * cm,
        title=lines[0].lstrip("# "),
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
    print("written", sys.argv[2])
