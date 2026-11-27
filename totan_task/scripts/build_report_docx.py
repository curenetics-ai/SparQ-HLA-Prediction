#!/usr/bin/env python3
"""
Render docs/week2_processing_report.md to a Word (.docx) file.

Built locally with python-docx rather than transcribing a base64 export, because
hand-copied base64 corrupts the archive (a single transposed character produces
a file with the right length and a failing CRC).

Supports the subset of Markdown used by the report: ATX headings, pipe tables,
fenced code blocks, bullet and ordered lists, and inline **bold**, *italic* and
`code`.

Run from the project root:
    python ./scripts/build_report_docx.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

TASK_ROOT = Path(__file__).resolve().parents[1]   # totan_task/
REPO_ROOT = Path(__file__).resolve().parents[2]   # repository root (Week 1 inputs)
ROOT = TASK_ROOT                                  # this task's outputs
SRC = ROOT / "reports" / "week2_processing_report.md"
OUT = ROOT / "reports" / "Week2_Antigen_Processing_Benchmark_Report.docx"

INLINE = re.compile(r"(\*\*.+?\*\*|(?<!\*)\*[^*]+?\*(?!\*)|`[^`]+?`)")


def add_runs(paragraph, text: str) -> None:
    """Write text into a paragraph, honouring **bold**, *italic* and `code`."""
    for part in INLINE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            paragraph.add_run(part[2:-2]).bold = True
        elif part.startswith("`") and part.endswith("`") and len(part) > 2:
            r = paragraph.add_run(part[1:-1])
            r.font.name = "Consolas"
            r.font.size = Pt(9.5)
            r.font.color.rgb = RGBColor(0xA0, 0x30, 0x30)
        elif part.startswith("*") and part.endswith("*") and len(part) > 2:
            paragraph.add_run(part[1:-1]).italic = True
        else:
            paragraph.add_run(part)


def split_row(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def is_separator(line: str) -> bool:
    return bool(re.fullmatch(r"\|[\s:|-]+\|", line.strip()))


def build() -> int:
    if not SRC.exists():
        sys.exit(f"missing source: {SRC}")
    lines = SRC.read_text(encoding="utf-8").splitlines()

    doc = Document()
    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(10.5)

    i = 0
    n_tables = n_heads = n_code = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        # fenced code block
        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Pt(18)
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(8)
            r = p.add_run("\n".join(buf))
            r.font.name = "Consolas"
            r.font.size = Pt(9)
            n_code += 1
            continue

        # table
        if stripped.startswith("|") and i + 1 < len(lines) and is_separator(lines[i + 1]):
            header = split_row(stripped)
            i += 2
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            t = doc.add_table(rows=1, cols=len(header))
            t.style = "Light Grid Accent 1"
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for c, text in enumerate(header):
                cell = t.rows[0].cells[c]
                cell.text = ""
                add_runs(cell.paragraphs[0], text)
                for run in cell.paragraphs[0].runs:
                    run.bold = True
                    run.font.size = Pt(9)
            for row in rows:
                cells = t.add_row().cells
                for c, text in enumerate(row[:len(header)]):
                    cells[c].text = ""
                    add_runs(cells[c].paragraphs[0], text)
                    for run in cells[c].paragraphs[0].runs:
                        run.font.size = Pt(9)
            doc.add_paragraph()
            n_tables += 1
            continue

        # headings
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            if level == 1:
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
                r = p.add_run(m.group(2))
                r.bold = True
                r.font.size = Pt(20)
            else:
                h = doc.add_heading(level=min(level, 4))
                for run in list(h.runs):
                    run.text = ""
                add_runs(h, m.group(2))
            n_heads += 1
            i += 1
            continue

        # bullet
        if stripped.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, stripped[2:])
            i += 1
            continue

        # ordered
        m = re.match(r"^(\d+)\.\s+(.*)$", stripped)
        if m:
            p = doc.add_paragraph(style="List Number")
            add_runs(p, m.group(2))
            i += 1
            continue

        p = doc.add_paragraph()
        add_runs(p, stripped)
        i += 1

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"  headings {n_heads}, tables {n_tables}, code blocks {n_code}")
    return 0


if __name__ == "__main__":
    sys.exit(build())
