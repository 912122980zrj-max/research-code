# -*- coding: utf-8 -*-
"""84_build_sci_rep_supplements.py

Build the Scientific Reports submission supplements:
  1. Supplementary_Information.pdf  - methods + figure/table legends text,
     followed by one page per supplementary figure (S1a-c, S2, S3a-d, S4a-d,
     S5a-o), each page carrying its legend above the figure.
  2. Supplementary_Data.zip        - raw supplementary-table CSVs, processed
     source matrices and data-preparation notes.
"""

import argparse
import os
import re
import zipfile
from pathlib import Path

import pymupdf


def md_to_blocks(path):
    blocks = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        level = 0
        while line.startswith("#"):
            level += 1
            line = line[1:].strip()
        if line.startswith(">"):
            line = line[1:].strip()
        line = re.sub(r"\*\*(.+?)\*\*", r"\1", line)
        line = re.sub(r"`(.+?)`", r"\1", line)
        line = line.replace("_", " ")
        if line.startswith(("- ", "* ")):
            line = "\u2022 " + line[2:]
        blocks.append((level, line))
    return blocks


def wrap_text(page, text, font, size, rect, bold=False):
    """Insert text with basic word wrapping; return the y position after."""
    words = text.split()
    line = ""
    y = rect.y0
    for w in words:
        trial = (line + " " + w).strip()
        wlen = pymupdf.get_text_length(trial, fontname=font, fontsize=size)
        if wlen > rect.width and line:
            page.insert_text((rect.x0, y), line, fontname=font, fontsize=size)
            y += size + 3.5
            line = w
            if y > rect.y1 - size:
                return y, True
        else:
            line = trial
    if line:
        page.insert_text((rect.x0, y), line, fontname=font, fontsize=size)
        y += size + 3.5
    return y, False


def render_text(pdf_path, md_paths, page_w=595.0, page_h=842.0):
    doc = pymupdf.open()
    page = doc.new_page(width=page_w, height=page_h)
    y = 48.0
    for md in md_paths:
        for level, text in md_to_blocks(md):
            if y > page_h - 80:
                page = doc.new_page(width=page_w, height=page_h)
                y = 48.0
            if level == 1:
                page.insert_text((50, y), text, fontname="hebo", fontsize=15)
                y += 22
            elif level == 2:
                page.insert_text((50, y), text, fontname="hebo", fontsize=12.5)
                y += 19
            else:
                rect = pymupdf.Rect(50, y, page_w - 50, page_h - 50)
                y, overflow = wrap_text(page, text, "helv", 10.5, rect)
                if overflow:
                    page = doc.new_page(width=page_w, height=page_h)
                    y = 48.0
                    y, _ = wrap_text(page, text, "helv", 10.5,
                                     pymupdf.Rect(50, y, page_w - 50, page_h - 50))
                y += 4
    doc.save(pdf_path, garbage=4, deflate=True, clean=True)
    doc.close()
    print("wrote", pdf_path)


def parse_legends(path):
    text = Path(path).read_text(encoding="utf-8")
    legends = {}
    pattern = re.compile(r"\*\*Fig\. S(\d+)\.[^*]+\*\*(.*?)(?=\*\*Fig\. S|\n## |\Z)", re.S)
    for m in pattern.finditer(text):
        num = m.group(1)
        body = re.sub(r"\*\*|`", "", m.group(2))
        body = " ".join(body.split())
        legends[num] = body
    return legends


def figure_order(figdir):
    names = [f"S1a_algorithm_benchmark.pdf", "S1b_oof_roc.pdf", "S1c_calibration_DCA.pdf",
             "S2_direction_replication.pdf",
             "S3a_QC_violin.pdf", "S3b_marker_dotplot.pdf", "S3c_marker_heatmap.pdf",
             "S3d_sample_composition.pdf",
             "S4a_MuSiC_proportions.pdf", "S4b_BisqueRNA_proportions.pdf",
             "S4c_CIBERSORTx_proportions.pdf", "S4d_LOO_validation.pdf",
             "S5a_TGFb_network_Sham.pdf", "S5b_TGFb_network_TBI.pdf",
             "S5c_KO_volcano_P2RY6.pdf", "S5d_KO_volcano_P2RX7.pdf",
             "S5e_docking_energy_summary.pdf"]
    for g in ["AHR", "CRABP2", "CYP1A2", "CYP1B1", "ESR1", "HSD11B1",
              "MAOA", "NFE2L2", "P2RY6", "PGF"]:
        letter = chr(102 + ["AHR", "CRABP2", "CYP1A2", "CYP1B1", "ESR1",
                            "HSD11B1", "MAOA", "NFE2L2", "P2RY6", "PGF"].index(g))
        names.append(f"S5{letter}_pose_{g}.png")
    return names


def append_figures(pdf_path, figdir, legends, page_w=595.0, page_h=842.0):
    doc = pymupdf.open(pdf_path)
    files = figure_order(figdir)
    for name in files:
        key = re.match(r"S(\d)", name).group(1)
        legend = legends.get(key, "")
        page = doc.new_page(width=page_w, height=page_h)
        label = re.sub(r"_.*$", "", name)
        header = f"Supplementary Figure {label}"
        page.insert_text((50, 40), header, fontname="hebo", fontsize=12)
        rect = pymupdf.Rect(50, 50, page_w - 50, 210)
        y_end, _ = wrap_text(page, legend, "helv", 9.5, rect)
        fig_top = max(y_end + 10, 126)
        img_path = os.path.join(figdir, name)
        target = pymupdf.Rect(50, fig_top, page_w - 50, page_h - 40)
        if name.endswith(".pdf"):
            src = pymupdf.open(img_path)
            p = src[0]
            pw, ph = p.rect.width, p.rect.height
            scale = min(target.width / pw, target.height / ph)
            w, h = pw * scale, ph * scale
            x0 = target.x0 + (target.width - w) / 2
            y0 = target.y0 + (target.height - h) / 2
            page.show_pdf_page(pymupdf.Rect(x0, y0, x0 + w, y0 + h), src, 0)
            src.close()
        else:
            page.insert_image(target, filename=img_path, keep_proportion=True)
    tmp = pdf_path + ".tmp"
    doc.save(tmp, garbage=4, deflate=True, clean=True)
    doc.close()
    os.replace(tmp, pdf_path)
    print("figures appended to", pdf_path)


def build_data_zip(zip_path, roots):
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for root in roots:
            for dirpath, _, files in os.walk(root):
                for f in files:
                    full = os.path.join(dirpath, f)
                    rel = os.path.relpath(full, os.path.dirname(root))
                    z.write(full, arcname=rel)
    print("wrote", zip_path, os.path.getsize(zip_path), "bytes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir",
                    default=r"E:\sheng xin\tbi\submission\03_Supplementary\SciRep")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    base = Path(r"E:\sheng xin\tbi\submission\03_Supplementary")
    figdir = base / "Figures"
    methods = base / "Text" / "SUPPLEMENTARY_METHODS_RESULTS_EN.md"
    legends_md = base / "Text" / "SUPPLEMENTARY_LEGENDS_EN.md"

    si_pdf = outdir / "Supplementary_Information.pdf"
    render_text(str(si_pdf), [str(methods), str(legends_md)])
    legends = parse_legends(str(legends_md))
    append_figures(str(si_pdf), str(figdir), legends)

    data_zip = outdir / "Supplementary_Data.zip"
    build_data_zip(str(data_zip), [
        str(base / "Tables" / "Raw_CSV"),
        str(Path(r"E:\sheng xin\tbi\submission\04_Source_Data")),
    ])


if __name__ == "__main__":
    main()
