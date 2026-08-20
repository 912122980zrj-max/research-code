# -*- coding: utf-8 -*-
"""
76_compound_control.py —— v7 补充：非 BaP 化合物对照（ChEMBL 单源）
依据：2026-08-11 终审意见（建议补做芘/菲对照）。
说明：SwissTargetPrediction/PharmMapper/SEA 在本会话中无法程序化提交
（SwissTargetPrediction 为 JS 表单；SEA 403；CTD 为 JS 路由），
因此对照使用 ChEMBL pChEMBL>=5 的人源单蛋白靶点，并如实标注为单源对照。
流程与主分析一致：Ensembl 1:1 直系同源 -> 5,043 矩阵覆盖 -> 364 关键基因交集
-> 置换检验（seed=2026, 10,000 次）。
输出：results/16_Orthology/v7_compound_control/
"""
import csv, json, random, sys, time, urllib.parse
from pathlib import Path

import requests

ROOT = Path(r"E:\sheng xin\tbi")
OUT = ROOT / "results" / "16_Orthology" / "v7_compound_control"
OUT.mkdir(parents=True, exist_ok=True)

BASE = "https://rest.ensembl.org"
H = {"Content-Type": "application/json", "Accept": "application/json"}

# ChEMBL pChEMBL>=5 人源单蛋白靶点（assay type 任意；B/F 且 conf>=7 为 0 条）
COMPOUNDS = {
    "pyrene": {
        "chembl_id": "CHEMBL279564",
        "cid": 31423,
        "targets": [
            ("CYP1A1", "Cytochrome P450 1A1", 7.39, "A"),
            ("CYP1A2", "Cytochrome P450 1A2", 8.15, "A"),
            ("CYP1B1", "Cytochrome P450 1B1", 8.70, "B"),
            ("CYP3A4", "Cytochrome P450 3A4", 5.00, "A"),
        ],
    },
    "phenanthrene": {
        "chembl_id": "CHEMBL46730",
        "cid": 705,
        "targets": [
            ("TDP1", "Tyrosyl-DNA phosphodiesterase 1", 5.70, "B"),
            ("TSHR", "Thyrotropin receptor", 6.00, "F"),
        ],
    },
}


def fetch_json(url, params=None):
    for attempt in range(3):
        try:
            r = requests.get(url, params=params, headers=H, timeout=60)
            if r.status_code == 200:
                return r.json()
            print("HTTP", r.status_code, url, file=sys.stderr)
        except Exception as exc:
            print("retry", url, repr(exc)[:100], file=sys.stderr)
        time.sleep(2)
    return None


def gene_to_id(symbol):
    url = f"{BASE}/xrefs/symbol/homo_sapiens/{urllib.parse.quote(symbol)}"
    d = fetch_json(url)
    if d and d:
        return d[0].get("id")
    return None


def orthologs(gene_id):
    url = f"{BASE}/homology/id/homo_sapiens/{gene_id}"
    d = fetch_json(url, params={"type": "orthologues",
                                "target_species": "mus_musculus"})
    out = []
    if not d or "data" not in d:
        return out
    mouse_ids = set()
    for item in d["data"]:
        for hom in item.get("homologies", []):
            mid = hom.get("target", {}).get("id")
            if mid:
                mouse_ids.add(mid)
            out.append({
                "mouse_id": mid,
                "mouse_label": None,
                "orthology_type": hom.get("type"),
                "orthology_class": hom.get("homology_type"),
            })
    for mid in mouse_ids:
        lu = fetch_json(f"{BASE}/lookup/id/{mid}")
        label = None
        if lu:
            label = lu.get("display_name")
        for o in out:
            if o["mouse_id"] == mid:
                o["mouse_label"] = label
    return out


def load_gene_list(path):
    with open(path, encoding="utf-8") as f:
        return [r["gene"] for r in csv.DictReader(f)]


matrix_genes = set(load_gene_list(ROOT / "results" / "16_Orthology" /
                                  "v7_single_cohort" /
                                  "deg_all_GSE58485_only.csv"))
key_genes = set(load_gene_list(ROOT / "results" / "16_Orthology" /
                               "v7_single_cohort" /
                               "key_genes_GSE58485_only.csv"))
universe = sorted(matrix_genes)
N = len(universe)
print("matrix genes:", N, "key genes:", len(key_genes))

rows = []
for compound, info in COMPOUNDS.items():
    mapped = []
    for symbol, pref, pchembl, atype in info["targets"]:
        gid = gene_to_id(symbol)
        if not gid:
            mapped.append({"human": symbol, "gene_id": None,
                           "mouse": None, "type": None})
            continue
        homs = orthologs(gid)
        one2one = [h for h in homs if h["orthology_type"] == "ortholog_one2one"]
        if one2one:
            h = one2one[0]
            mapped.append({"human": symbol, "gene_id": gid,
                           "mouse": h["mouse_label"], "type": "1:1"})
        else:
            labels = [h["mouse_label"] for h in homs if h.get("mouse_label")]
            mapped.append({"human": symbol, "gene_id": gid,
                           "mouse": ",".join(labels) if labels else None,
                           "type": "1:many" if homs else "none"})
    eval_genes = sorted({m["mouse"] for m in mapped
                         if m.get("type") == "1:1" and m["mouse"] in matrix_genes})
    cand = sorted(set(eval_genes) & key_genes)
    # 置换（K = 候选数；无候选则记 0）
    random.seed(2026)
    k = len(cand)
    obs = k
    expected = k * len(key_genes) / N
    emp_p = None
    ci = None
    if k > 0:
        key_upper = {g.upper() for g in key_genes}
        null = []
        uni = universe
        for _ in range(10000):
            picked = random.sample(uni, k)
            null.append(sum(1 for g in picked if g.upper() in key_upper))
        emp_p = (1 + sum(x >= obs for x in null)) / 10001
        ci = (min(null), max(null))
        fold = obs / expected if expected > 0 else float("nan")
    else:
        fold = 0.0
    rows.append({
        "compound": compound, "chembl_id": info["chembl_id"],
        "pubchem_cid": info["cid"],
        "human_targets": ";".join(s for s, *_ in info["targets"]),
        "mapped_1to1": ";".join(m["mouse"] for m in mapped
                                if m.get("type") == "1:1" and m.get("mouse")),
        "evaluable": ";".join(eval_genes),
        "n_evaluable": len(eval_genes),
        "candidates": ";".join(cand),
        "n_candidates": len(cand),
        "expected": round(expected, 6),
        "enrichment_fold": round(fold, 3) if fold == fold else "NA",
        "empirical_p": emp_p,
        "null_range": f"{ci}" if ci else "NA",
    })
    print(compound, "| mapped:", [(m["human"], m["mouse"], m["type"])
                                  for m in mapped],
          "| evaluable:", eval_genes, "| candidates:", cand,
          "| expected:", round(expected, 6),
          "| fold:", round(fold, 3) if fold == fold else "NA")

with open(OUT / "compound_control_results.csv", "w", newline="",
          encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

with open(OUT / "COMPOUND_CONTROL_REPORT.txt", "w", encoding="utf-8") as f:
    f.write("v7 non-BaP compound control (ChEMBL single-source)\n")
    f.write(f"generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"universe: {N} genes; key genes: {len(key_genes)}; "
            f"seed=2026; 10,000 draws\n")
    f.write("Note: SwissTargetPrediction/PharmMapper/SEA/CTD not "
            "programmatically accessible in this session; control is "
            "ChEMBL-only and is a partial control.\n")
    for row in rows:
        f.write(json.dumps(row, ensure_ascii=False) + "\n")
print("done")
