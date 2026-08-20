# -*- coding: utf-8 -*-
"""83_composite_figures.py

Merge multi-panel main figures into single composite PDFs as required by the
BMC Bioinformatics submission guidelines ("multi-panel figures should be
submitted as a single composite file"):

  Fig3.pdf  = Fig3A (volcano) + Fig3B (module-trait) + Fig3C (candidate
              logFC) + Fig3D (permutation background), arranged 2x2.
  Fig4.pdf  = Fig4A (GO bubble, full width) + Fig4B (AhR coverage) +
              Fig4C (sensitivity heatmap), arranged A over B/C.

Panel letters are added in a narrow strip above each panel. Vector content is
preserved (source pages are embedded as PDF objects, not rasterised).
"""

import argparse
import os
from pathlib import Path

import pymupdf


def fit_rect(panel_w: float, panel_h: float, cell: pymupdf.Rect) -> pymupdf.Rect:
    """Return a rect inside `cell` that preserves the source aspect ratio."""
    scale = min(cell.width / panel_w, cell.height / panel_h)
    w = panel_w * scale
    h = panel_h * scale
    x0 = cell.x0 + (cell.width - w) / 2
    y0 = cell.y0 + (cell.height - h) / 2
    return pymupdf.Rect(x0, y0, x0 + w, y0 + h)


def build_fig3(panels, out):
    doc = pymupdf.open()
    page_w, margin, gap, strip = 1000.0, 20.0, 20.0, 30.0
    col_w = (page_w - 2 * margin - gap) / 2
    src = [pymupdf.open(p) for p in panels]
    sizes = [(p[0].rect.width, p[0].rect.height) for p in src]
    # row heights = max(panel height at col width) + letter strip
    row_h = []
    for pair in [(0, 1), (2, 3)]:
        hh = []
        for i in pair:
            w, h = sizes[i]
            hh.append(h * (col_w / w))
        row_h.append(max(hh) + strip)
    page_h = 2 * margin + row_h[0] + gap + row_h[1]
    page = doc.new_page(width=page_w, height=page_h)
    cells = []
    y = margin
    for r, pair in enumerate([(0, 1), (2, 3)]):
        for k, i in enumerate(pair):
            x = margin + k * (col_w + gap)
            cells.append(pymupdf.Rect(x, y + strip, x + col_w, y + row_h[r]))
        y += row_h[r] + gap
    for i, (s, rect) in enumerate(zip(src, cells)):
        w, h = sizes[i]
        page.show_pdf_page(fit_rect(w, h, rect), s, 0)
        page.insert_text((rect.x0 + 4, rect.y0 - 8),
                        "abcd"[i], fontsize=11, fontname="hebo",
                         color=(0, 0, 0))
    doc.save(out, deflate=True)
    for s in src:
        s.close()
    print(f"wrote {out} ({page_w:.0f} x {page_h:.0f} pt)")


def build_fig4(panels, out):
    doc = pymupdf.open()
    page_w, margin, gap, strip = 1000.0, 20.0, 20.0, 30.0
    src = [pymupdf.open(p) for p in panels]
    sizes = [(p[0].rect.width, p[0].rect.height) for p in src]
    wA, hA = sizes[0]
    wB, hB = sizes[1]
    wC, hC = sizes[2]
    row1_w = page_w - 2 * margin
    row1_h = hA * (row1_w / wA) + strip
    col_w = (page_w - 2 * margin - gap) / 2
    row2_h = max(hB * (col_w / wB), hC * (col_w / wC)) + strip
    page_h = 2 * margin + row1_h + gap + row2_h
    page = doc.new_page(width=page_w, height=page_h)

    rects = [
        pymupdf.Rect(margin, margin + strip, margin + row1_w, margin + row1_h),
        pymupdf.Rect(margin, margin + row1_h + gap + strip,
                     margin + col_w, margin + row1_h + gap + row2_h),
        pymupdf.Rect(margin + col_w + gap, margin + row1_h + gap + strip,
                     margin + col_w + gap + col_w,
                     margin + row1_h + gap + row2_h),
    ]
    for i, (s, rect) in enumerate(zip(src, rects)):
        w, h = sizes[i]
        page.show_pdf_page(fit_rect(w, h, rect), s, 0)
        page.insert_text((rect.x0 + 4, rect.y0 - 8),
                        "abc"[i], fontsize=11, fontname="hebo",
                         color=(0, 0, 0))
    doc.save(out, deflate=True)
    for s in src:
        s.close()
    print(f"wrote {out} ({page_w:.0f} x {page_h:.0f} pt)")


def build_fig1_from_png(png, out, page_w_pt=595.0):
    """Single-page PDF from the canonical Fig1 PNG (~360 dpi at A4 width)."""
    from PIL import Image
    im = Image.open(png)
    page_h = page_w_pt * im.height / im.width
    doc = pymupdf.open()
    page = doc.new_page(width=page_w_pt, height=page_h)
    page.insert_image(page.rect, filename=png)
    doc.save(out, deflate=True)
    print(f"wrote {out} ({page_w_pt:.0f} x {page_h:.0f} pt)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figdir", default=r"E:\sheng xin\tbi\figures")
    args = ap.parse_args()
    figdir = Path(args.figdir)

    f3 = figdir / "Figure_03_Discovery"
    build_fig3(
        [f3 / "volcano_deg_GSE58485_only.pdf",
         f3 / "module_trait_GSE58485_only.pdf",
         f3 / "figure3E_candidate_logFC_GSE58485_only.pdf",
         f3 / "permutation_background.pdf"],
        figdir / "Figure_03_Discovery" / "Fig3.pdf")

    f4 = figdir / "Figure_04_Enrichment"
    build_fig4(
        [f4 / "GO_BP_bubble.pdf",
         f4 / "figure4B_ahr_v6.pdf",
         f3 / "sensitivity_heatmap.pdf"],
        figdir / "Figure_04_Enrichment" / "Fig4.pdf")

    build_fig1_from_png(
        r"E:\sheng xin\tbi\figures\Figure_01_Workflow\Figure1_workflow_v5.png",
        figdir / "Figure_01_Workflow" / "Fig1_StudyDesign_singlepage.pdf")


if __name__ == "__main__":
    main()
