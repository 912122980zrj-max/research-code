# -*- coding: utf-8 -*-
"""
02b_tier_targets.py —— BaP 靶点证据分级（审查交接 T3 / P6）

对 data/04_BaP_targets/ 下各来源靶点按证据强度分级：
  - ChEMBL：按 BaP 分子（CHEMBL312）的活性记录重拉 API；
    Tier A = 存在 pChEMBL>=5 且 assay_type in {B,F} 且 confidence>=7 的活性；
    其余 ChEMBL 靶点一律 Tier C。
  - SwissTargetPrediction：Probability* >=0.5 -> A；0.1-0.5 -> B；<0.1 -> C。
  - PharmMapper：z>=1.0 -> B；0.5-1.0 -> C。
  - SEA：保留但标注 Tanimoto<0.5 弱相似，定级 C。

输出：
  data/04_BaP_targets/union_targets_tiered.tsv
  data/04_BaP_targets/target_summary_tiered.json
  data/04_BaP_targets/chembl_activity_cache.tsv（首次 API 拉取后的缓存）

运行：python scripts/02b_tier_targets.py [--refresh]
"""

import argparse
import csv
import glob
import io
import json
import os
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (TBI-pipeline/1.0; research; contact: local)"}
BA_CHEMBL_ID = "CHEMBL31184"  # PubChem RegistryID 交叉引用确认（SEA 文件名中的 CHEMBL312 非 ChEMBL 分子号）

# 官方符号别名映射（第二轮审查 §4.1-1：Tier A 含 BHLHE76(=AHR) 与 CYPIA2(=CYP1A2)）
ALIAS_MAP = {
    "BHLHE76": "AHR",
    "CYPIA2": "CYP1A2",
}


def log(msg):
    print(msg, flush=True)


def http_json(url, timeout=60, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"HTTP failed: {url} -> {last}")


def detect_delimiter(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except Exception:  # noqa: BLE001
        return "\t"


# ---------------- ChEMBL 活性记录 ----------------
def chembl_activities(mol_id, organism="Homo sapiens"):
    hits = []
    offset = 0
    limit = 100
    while True:
        q = urllib.parse.urlencode({
            "molecule_chembl_id": mol_id,
            "target_organism": organism,
            "limit": limit,
            "offset": offset,
        })
        d = http_json("https://www.ebi.ac.uk/chembl/api/data/activity.json?" + q,
                      timeout=90)
        batch = d.get("activities", [])
        hits.extend(batch)
        if len(batch) < limit:
            break
        offset += limit
    return hits


def load_chembl_cache(out_dir):
    p = os.path.join(out_dir, "chembl_cache.tsv")
    mapping = {}  # target_id -> set(gene)
    if not os.path.exists(p):
        return mapping
    with open(p, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = (row.get("gene_symbol") or "").strip()
            tid = (row.get("target_id") or "").strip()
            if g and tid:
                mapping.setdefault(tid.upper(), set()).add(g.upper())
    return mapping


def chembl_target_genes_api(tid):
    """单个 target 的基因名提取（与 02 脚本一致，仅取直接符号与括号符号）。"""
    d = http_json(f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json",
                  timeout=90)
    genes = set()
    for tc in d.get("target_components", []):
        comp = tc.get("component", tc)
        for syn in comp.get("target_component_synonyms", []):
            s = str(syn.get("component_synonym", ""))
            st = str(syn.get("syn_type", ""))
            if st == "GENE_SYMBOL" or re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", s):
                genes.add(s)
        desc = str(comp.get("description", ""))
        for m in re.finditer(r"\(([A-Z][A-Z0-9]{1,10})\)", desc):
            genes.add(m.group(1))
    return sorted(genes)


def uniprot_map(accs, to="GENE_NAME"):
    """UniProt accession -> 基因名（批量查询，仅少量回退时使用）。"""
    if not accs:
        return {}
    mapping = {}
    accs = sorted(set(accs))
    for i in range(0, len(accs), 80):
        chunk = accs[i:i + 80]
        query = " OR ".join("accession:" + urllib.parse.quote(a) for a in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + urllib.parse.quote(query)
               + "&fields=accession,gene_primary&format=json&size=500")
        try:
            d = http_json(url, timeout=90)
            for r in d.get("results", []):
                gene = None
                for g in r.get("genes", []):
                    gn = g.get("geneName") or {}
                    if gn.get("value"):
                        gene = gn["value"]
                        break
                if gene:
                    mapping.setdefault(r["primaryAccession"], gene)
        except Exception as e:  # noqa: BLE001
            log(f"  [warn] UniProt 批量查询失败（{e}），继续")
    return mapping


def uniprot_mnemonic_map(mnemonics):
    mapping = {}
    mnemonics = sorted(set(mnemonics))
    for i in range(0, len(mnemonics), 80):
        chunk = mnemonics[i:i + 80]
        query = " OR ".join("id:" + urllib.parse.quote(m) for m in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + urllib.parse.quote(query)
               + "&fields=accession,gene_primary,id&format=json&size=500")
        try:
            d = http_json(url, timeout=90)
            for r in d.get("results", []):
                gene = None
                for g in r.get("genes", []):
                    gn = g.get("geneName") or {}
                    if gn.get("value"):
                        gene = gn["value"]
                        break
                entry = r.get("uniProtkbId", "")
                if gene and entry:
                    mapping.setdefault(entry, gene)
                    mapping.setdefault(entry.split("_")[0], gene)
        except Exception as e:  # noqa: BLE001
            log(f"  [warn] UniProt 条目名批量查询失败（{e}），继续")
    return mapping


def load_activity_cache(out_dir):
    p = os.path.join(out_dir, "chembl_activity_cache.tsv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save_activity_cache(out_dir, records):
    p = os.path.join(out_dir, "chembl_activity_cache.tsv")
    fieldnames = ["target_chembl_id", "assay_chembl_id", "pchembl_value",
                  "assay_type",
                  "confidence_score", "standard_type", "standard_value"]
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(records)


def load_assay_cache(out_dir):
    p = os.path.join(out_dir, "chembl_assay_cache.tsv")
    if not os.path.exists(p):
        return {}
    out = {}
    with open(p, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            aid = (row.get("assay_chembl_id") or "").strip()
            if aid:
                out[aid] = row
    return out


def save_assay_cache(out_dir, records):
    p = os.path.join(out_dir, "chembl_assay_cache.tsv")
    fieldnames = ["assay_chembl_id", "confidence_score",
                  "confidence_description"]
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(records)


def fetch_assay_confidence(assay_ids, out_dir):
    cache = load_assay_cache(out_dir)
    missing = sorted(set(aid for aid in assay_ids if aid) - set(cache))
    for aid in missing:
        try:
            a = http_json(f"https://www.ebi.ac.uk/chembl/api/data/assay/{aid}.json",
                          timeout=60)
            cache[aid] = {
                "assay_chembl_id": aid,
                "confidence_score": str(a.get("confidence_score") or "").strip(),
                "confidence_description": str(
                    a.get("confidence_description") or "").strip(),
            }
            time.sleep(0.15)
        except Exception as e:  # noqa: BLE001
            log(f"  [warn] assay {aid} 抓取失败: {e}")
            cache[aid] = {"assay_chembl_id": aid, "confidence_score": "",
                          "confidence_description": "fetch_failed"}
    save_assay_cache(out_dir, [cache[k] for k in sorted(cache)])
    return cache


def tier_chembl(out_dir, union_genes, refresh=False):
    log("== ChEMBL 分级 ==")
    cache = load_activity_cache(out_dir)
    if cache is not None and not refresh:
        log(f"  使用 ChEMBL 活性缓存：{len(cache)} 条")
    else:
        log(f"  拉取 {BA_CHEMBL_ID} 的 Homo sapiens 活性记录（分页，可能需要数分钟）...")
        acts = chembl_activities(BA_CHEMBL_ID)
        log(f"  活性记录：{len(acts)} 条")
        cache = []
        for a in acts:
            cache.append({
                "target_chembl_id": (a.get("target_chembl_id") or "").strip(),
                "assay_chembl_id": (a.get("assay_chembl_id") or "").strip(),
                "pchembl_value": str(a.get("pchembl_value") or "").strip(),
                "assay_type": str(a.get("assay_type") or "").strip(),
                "confidence_score": str(a.get("confidence_score") or "").strip(),
                "standard_type": str(a.get("standard_type") or "").strip(),
                "standard_value": str(a.get("standard_value") or "").strip(),
            })
        save_activity_cache(out_dir, cache)

    assay_conf = fetch_assay_confidence(
        [r["assay_chembl_id"] for r in cache], out_dir)

    target_genes = load_chembl_cache(out_dir)
    missing_tids = sorted({r["target_chembl_id"] for r in cache
                           if r["target_chembl_id"]} - set(target_genes))
    for tid in missing_tids:
        try:
            gs = chembl_target_genes_api(tid)
            target_genes.setdefault(tid, set()).update(g.upper() for g in gs)
            time.sleep(0.2)
        except Exception as e:  # noqa: BLE001
            log(f"  [warn] target {tid} 基因名解析失败: {e}")

    per_gene = {}  # gene -> list of qualifying info
    all_genes = set()
    for r in cache:
        tid = r["target_chembl_id"]
        if not tid:
            continue
        genes = target_genes.get(tid, set())
        if not genes:
            continue
        try:
            pchembl = float(r["pchembl_value"])
        except (TypeError, ValueError):
            pchembl = float("nan")
        try:
            conf_rec = assay_conf.get(r["assay_chembl_id"], {})
            conf = int(conf_rec.get("confidence_score", "") or -1)
        except (TypeError, ValueError):
            conf = -1
        assay = r["assay_type"].upper()
        qualifies = (pchembl >= 5 and assay in ("B", "F") and conf >= 7)
        for g in genes:
            if g not in union_genes:
                continue
            all_genes.add(g)
            info = per_gene.setdefault(g, {
                "max_pchembl": float("-inf"),
                "n_qualifying": 0,
                "n_activities": 0,
                "assay_types": set(),
                "max_confidence": -1,
            })
            info["n_activities"] += 1
            info["assay_types"].add(assay or "NA")
            info["max_confidence"] = max(info["max_confidence"], conf)
            if pchembl == pchembl:  # not NaN
                info["max_pchembl"] = max(info["max_pchembl"], pchembl)
            if qualifies:
                info["n_qualifying"] += 1

    rows = []
    for g in sorted(all_genes):
        info = per_gene[g]
        maxp = info["max_pchembl"]
        maxp_s = (f"{maxp:.2f}" if math.isfinite(maxp) else "NA")
        tier = "A" if info["n_qualifying"] > 0 else "C"
        detail = (f"n_qualifying={info['n_qualifying']}; n_activities="
                  f"{info['n_activities']}; max_pchembl={maxp_s}; "
                  f"assay_types={','.join(sorted(info['assay_types']))}; "
                  f"max_assay_confidence={info['max_confidence']}")
        rows.append({"gene_symbol": g, "source": "ChEMBL", "tier": tier,
                     "score": maxp_s, "detail": detail,
                     "target_id": ""})
    log(f"  ChEMBL 分级基因（在并集内）：{len(rows)}；Tier A："
        f"{sum(1 for r in rows if r['tier'] == 'A')}")
    return rows


# ---------------- SwissTargetPrediction ----------------
def parse_swisstarget(path, union_genes):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return rows
        fld = [x.strip() for x in reader.fieldnames]
        name_col = None
        prob_col = None
        for fn in fld:
            if "common name" in fn.lower() or "gene" in fn.lower():
                name_col = fn
            if "probability" in fn.lower():
                prob_col = fn
        if name_col is None:
            name_col = fld[1] if len(fld) > 1 else fld[0]
        for row in reader:
            r = {k.strip(): (v or "").strip() for k, v in row.items()}
            val = r.get(name_col, "")
            genes = set()
            if re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", val):
                genes.add(val)
            for m in re.finditer(r"\(([A-Z][A-Z0-9]{1,10})\)", val):
                genes.add(m.group(1))
            prob = None
            if prob_col:
                pv = r.get(prob_col, "").replace("%", "").strip()
                try:
                    prob = float(pv)
                except ValueError:
                    prob = None
            for g in genes:
                if g.upper() not in union_genes:
                    continue
                tier = "A" if prob is not None and prob >= 0.5 else (
                    "B" if prob is not None and prob >= 0.1 else "C")
                rows.append({"gene_symbol": g.upper(), "source":
                             "SwissTargetPrediction", "tier": tier,
                             "score": "NA" if prob is None else f"{prob:.4g}",
                             "detail": "probability", "target_id": ""})
    return rows


# ---------------- PharmMapper ----------------
def parse_pharmmapper(path, union_genes):
    lines = open(path, encoding="utf-8", errors="replace").readlines()
    start = None
    for i, line in enumerate(lines):
        if "Pharma Model" in line and ("Uniplot" in line or "UniProt" in line):
            start = i
            break
    if start is None:
        return []
    reader = csv.DictReader(io.StringIO("".join(lines[start:])))
    rows_raw = []
    mnemonics = set()
    accs = set()
    for row in reader:
        r = {k.strip(): (v or "").strip() for k, v in row.items()}
        try:
            z = float(r.get("zscore", "0") or 0)
        except ValueError:
            z = 0.0
        uni = r.get("Uniplot", "") or r.get("Target Name", "") or ""
        tok = uni.split("_")[0] if "_" in uni else uni
        rows_raw.append({"token": tok, "z": z, "uni": uni})
        if re.fullmatch(r"[A-Z0-9]{2,10}_HUMAN", uni):
            mnemonics.add(uni)
        elif re.fullmatch(r"[A-Z0-9]{2,10}", tok):
            mnemonics.add(tok)
        if re.fullmatch(r"[A-NP-RT-Z][0-9][A-Z0-9]{4}", tok):
            accs.add(tok)
    # 带 "_HUMAN" 的条目一律用 UniProt 官方符号（避免 ADA17/ADAM17 这类别名陷阱）
    need_map = sorted(m for m in mnemonics
                      if "_" in m or m.upper() not in union_genes)
    mapped = uniprot_mnemonic_map(need_map) if need_map else {}
    acc_map = uniprot_map(accs) if accs else {}
    gene_for = {}
    for m in mnemonics:
        if m in mapped:
            gene_for[m] = mapped[m].upper()
        else:
            short = m.split("_")[0].upper()
            if short in union_genes:
                gene_for[m] = short
    for acc, gene in acc_map.items():
        if gene.upper() in union_genes:
            gene_for[acc] = gene.upper()
    best = {}  # gene -> max z
    for rr in rows_raw:
        g = gene_for.get(rr["uni"]) or gene_for.get(rr["token"])
        if not g or g not in union_genes:
            continue
        best[g] = max(best.get(g, -1), rr["z"])
    out = []
    for g, z in sorted(best.items()):
        tier = "B" if z >= 1.0 else "C"
        out.append({"gene_symbol": g, "source": "PharmMapper", "tier": tier,
                    "score": f"{z:.3f}", "detail": "zscore", "target_id": ""})
    return out


# ---------------- SEA ----------------
def parse_sea(out_dir, union_genes):
    rows = []
    sea_files = glob.glob(os.path.join(out_dir, "sea_result*.tsv")) + \
        glob.glob(os.path.join(out_dir, "sea_result*.csv"))
    tanimoto = None
    top_files = glob.glob(os.path.join(out_dir, "sea_result_*", "top10_*.tsv"))
    if top_files:
        with open(top_files[0], encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    tanimoto = float(row.get("Tanimoto_Score", ""))
                    break
                except (TypeError, ValueError):
                    continue
    seen = set()
    for sf in sea_files:
        with open(sf, encoding="utf-8", errors="replace") as f:
            delim = detect_delimiter(sf)
            reader = csv.DictReader(f, delimiter=delim)
            for row in reader:
                up = (row.get("Uniprot_ID") or "").strip()
                try:
                    z = float(row.get("Z_score", "0") or 0)
                except ValueError:
                    z = float("nan")
                for acc, g in uniprot_map([up]).items() if up else []:
                    g = g.upper()
                    if g in union_genes and g not in seen:
                        seen.add(g)
                        tani_s = "NA" if tanimoto is None else f"{tanimoto:.3f}"
                        rows.append({
                            "gene_symbol": g, "source": "SEA", "tier": "C",
                            "score": tani_s,
                            "detail": (f"SEA_Z={z:.2f}; tanimoto={tani_s}"),
                            "target_id": ""})
    return rows


# ---------------- 主流程 ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--refresh", action="store_true",
                    help="忽略 ChEMBL 活性缓存重新拉取")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    out_dir = args.out or os.path.join(root, "data", "04_BaP_targets")
    os.makedirs(out_dir, exist_ok=True)

    union_path = os.path.join(out_dir, "union_targets.tsv")
    union_genes = set()
    with open(union_path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            union_genes.add((row.get("gene_symbol") or "").strip().upper())
    log(f"并集靶点：{len(union_genes)}")

    all_rows = []
    all_rows += tier_chembl(out_dir, union_genes, refresh=args.refresh)

    stp_files = glob.glob(os.path.join(out_dir, "*wiss*arget*.csv")) + \
        glob.glob(os.path.join(out_dir, "swiss*.csv"))
    stp_path = None
    for p in stp_files:
        if "SwissTargetPrediction" in os.path.basename(p) or "swisstarget" in os.path.basename(p).lower():
            stp_path = p
            break
    if stp_path is None and stp_files:
        stp_path = stp_files[0]
    if stp_path:
        stp_rows = parse_swisstarget(stp_path, union_genes)
        log(f"SwissTargetPrediction 分级基因：{len(stp_rows)}；"
            f"Tier A/B/C = "
            f"{sum(1 for r in stp_rows if r['tier']=='A')}/"
            f"{sum(1 for r in stp_rows if r['tier']=='B')}/"
            f"{sum(1 for r in stp_rows if r['tier']=='C')}")
        all_rows += stp_rows
    else:
        log("  [warn] 未找到 SwissTargetPrediction 结果文件")

    pm_files = glob.glob(os.path.join(out_dir, "pharmmapper*"))
    if pm_files:
        pm_rows = parse_pharmmapper(pm_files[0], union_genes)
        log(f"PharmMapper 分级基因：{len(pm_rows)}；Tier B/C = "
            f"{sum(1 for r in pm_rows if r['tier']=='B')}/"
            f"{sum(1 for r in pm_rows if r['tier']=='C')}")
        all_rows += pm_rows
    else:
        log("  [warn] 未找到 PharmMapper 结果文件")

    all_rows += parse_sea(out_dir, union_genes)

    # 去重（同一 gene+source 取最高 tier）
    dedup = {}
    for r in all_rows:
        key = (r["gene_symbol"], r["source"])
        if key not in dedup:
            dedup[key] = r
        else:
            old = dedup[key]
            order = {"A": 0, "B": 1, "C": 2}
            if order.get(r["tier"], 9) < order.get(old["tier"], 9):
                dedup[key] = r
    rows = sorted(dedup.values(), key=lambda r: (r["gene_symbol"], r["source"]))

    # 第二轮审查 §4.1-1：附加官方符号列（别名→官方符号）
    for r in rows:
        r["official_symbol"] = ALIAS_MAP.get(r["gene_symbol"], r["gene_symbol"])

    out_tsv = os.path.join(out_dir, "union_targets_tiered.tsv")
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_symbol", "source", "tier",
                                          "score", "detail", "target_id",
                                          "official_symbol"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"写出：{out_tsv}（{len(rows)} 行）")

    # 去别名后的唯一基因汇总（第二轮审查通过条件 1）
    tier_order = {"A": 0, "B": 1, "C": 2}
    dedup_by_symbol = {}
    for r in rows:
        sym = r["official_symbol"]
        key = (sym, r["source"])
        if key not in dedup_by_symbol:
            dedup_by_symbol[key] = dict(r, gene_symbol=sym)
        else:
            old = dedup_by_symbol[key]
            if tier_order.get(r["tier"], 9) < tier_order.get(old["tier"], 9):
                dedup_by_symbol[key] = dict(r, gene_symbol=sym)
    dedup_rows = sorted(dedup_by_symbol.values(),
                        key=lambda r: (r["gene_symbol"], r["source"]))
    dedup_tsv = os.path.join(out_dir, "union_targets_tiered_dedup.tsv")
    with open(dedup_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_symbol", "source", "tier",
                                          "score", "detail", "target_id",
                                          "official_symbol"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(dedup_rows)
    log(f"写出（去别名）：{dedup_tsv}（{len(dedup_rows)} 行）")

    best_tier = {}
    for r in rows:
        g = r["official_symbol"]
        cur = best_tier.get(g, "Z")
        if tier_order.get(r["tier"], 9) < tier_order.get(cur, 9):
            best_tier[g] = r["tier"]
    tier_a_genes = sorted(g for g, t in best_tier.items() if t == "A")
    tier_a_rows = [r for r in rows if r["tier"] == "A"]
    summary = {
        "n_union": len(union_genes),
        "n_tiered_rows": len(rows),
        "n_tiered_rows_dedup": len(dedup_rows),
        "n_genes_with_tier": len(best_tier),
        "tier_counts": {t: sum(1 for v in best_tier.values() if v == t)
                        for t in ("A", "B", "C")},
        "tierA_unique_genes": tier_a_genes,
        "tierA_rows": [r["gene_symbol"] for r in tier_a_rows],
        "alias_pairs_resolved": [{"alias": a, "official": o}
                                 for a, o in sorted(ALIAS_MAP.items())],
        "signature_genes": {g: best_tier.get(g, "MISSING")
                            for g in ("P2RY6", "CCNA2", "HSD11B1")},
    }
    with open(os.path.join(out_dir, "target_summary_tiered.json"), "w",
              encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    log("分级汇总：")
    log(json.dumps(summary, ensure_ascii=False, indent=2))
    log("完成。")


if __name__ == "__main__":
    main()
    sys.exit(0)
