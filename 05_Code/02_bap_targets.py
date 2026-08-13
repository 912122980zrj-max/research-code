# -*- coding: utf-8 -*-
"""
02_bap_targets.py —— 苯并[a]芘（BaP）靶点预测（复现模板论文的"三库并集"思路）

自动完成：
  1) PubChem 获取 BaP 的 SMILES / InChIKey
  2) ChEMBL API：按 InChIKey 查化合物 -> 拉取 Homo sapiens 生物活性靶点 -> 基因名
  3) CTD：下载 BaP 的化学-基因互作表（若 CTD 下载被反爬，则提示手动下载）
  4) 若存在手动提交的 SEA / PharmMapper 结果文件，一并解析并入并集
  5) 合并、去重、输出 data/04_BaP_targets/union_targets.tsv

运行：python scripts/02_bap_targets.py
要求：Python 3.9+，仅标准库。
"""

import argparse
import csv
import glob
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import http.cookiejar

UA = {"User-Agent": "Mozilla/5.0 (TBI-pipeline/1.0; research; contact: local)"}
BA_NAME = "benzo[a]pyrene"


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


def http_text(url, timeout=90, retries=3):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"HTTP failed: {url} -> {last}")


# ---------------- 1. PubChem ----------------
def pubchem_props():
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(BA_NAME)
           + "/property/SMILES,ConnectivitySMILES,InChIKey,IUPACName/JSON")
    d = http_json(url)
    props = d["PropertyTable"]["Properties"][0]
    return {
        "smiles": props.get("SMILES") or props.get("ConnectivitySMILES", ""),
        "isomeric_smiles": props.get("SMILES", ""),
        "inchikey": props.get("InChIKey", ""),
        "iupac": props.get("IUPACName", ""),
    }


# ---------------- 2. ChEMBL ----------------
def chembl_id_from_pubchem(inchikey):
    """ChEMBL 收录但 inspect 接口查不到时，用 PubChem 的 RegistryID 交叉引用解析。"""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
           + urllib.parse.quote(inchikey) + "/xrefs/RegistryID/JSON")
    d = http_json(url)
    ids = []
    for info in d.get("InformationList", {}).get("Information", []):
        ids.extend(info.get("RegistryID", []))
    chembl = [x for x in ids if re.fullmatch(r"CHEMBL\d+", x)]
    return chembl[0] if chembl else None


def chembl_compound(inchikey):
    try:
        d = http_json(f"https://www.ebi.ac.uk/chembl/api/data/molecule/inspect/{inchikey}.json")
        if d.get("molecule_chembl_id"):
            return d
    except Exception:  # noqa: BLE001
        pass
    # 回退：PubChem 交叉引用 -> ChEMBL ID
    mol_id = chembl_id_from_pubchem(inchikey)
    if mol_id:
        d = http_json(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{mol_id}.json")
        return d
    return {}


def chembl_bioactivities(mol_id):
    hits = []
    offset = 0
    limit = 100
    while True:
        q = urllib.parse.urlencode({
            "molecule_chembl_id": mol_id,
            "target_organism": "Homo sapiens",
            "limit": limit,
            "offset": offset,
        })
        d = http_json("https://www.ebi.ac.uk/chembl/api/data/activity.json?" + q)
        batch = d.get("activities", [])
        hits.extend(batch)
        if len(batch) < limit:
            break
        offset += limit
    return hits


def chembl_target_genes(target_id):
    """从靶点 JSON 中提取基因名/UniProt 号。"""
    d = http_json(f"https://www.ebi.ac.uk/chembl/api/data/target/{target_id}.json")
    genes, accs = set(), set()
    for tc in d.get("target_components", []):
        comp = tc.get("component", tc)  # 不同版本结构略有差异
        acc = comp.get("accession")
        if acc:
            accs.add(acc)
        for syn in comp.get("target_component_synonyms", []):
            s = str(syn.get("component_synonym", ""))
            st = str(syn.get("syn_type", ""))
            if st == "GENE_SYMBOL" or re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", s):
                genes.add(s)
        desc = str(comp.get("description", ""))
        for m in re.finditer(r"\(([A-Z][A-Z0-9]{1,10})\)", desc):
            genes.add(m.group(1))
    return sorted(genes), sorted(accs)


def uniprot_map(accs, to="GENE_NAME"):
    """UniProt accession -> 基因名（用 search API 批量查，稳定且快）。"""
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
            print(f"  [warn] UniProt 批量查询失败（{e}），继续")
    return mapping


# ---------------- ChEMBL 缓存（加快重复运行） ----------------
def load_chembl_cache(out_dir):
    p = os.path.join(out_dir, "chembl_cache.tsv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save_chembl_cache(out_dir, records):
    p = os.path.join(out_dir, "chembl_cache.tsv")
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_symbol", "source", "target_id"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(records)


# ---------------- 3. CTD ----------------
def parse_ctd_file(path):
    genes = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = next(csv.reader(io.StringIO(line)))
            if len(parts) < 2:
                continue
            if parts[0].lower() in ("gene symbol", "input", "input id"):
                continue
            try:
                n = int(parts[2]) if len(parts) >= 3 else 1
            except ValueError:
                n = 1
            genes.append((parts[0].strip(), n))
    return [(g, n) for g, n in genes if re.fullmatch(r"[A-Za-z][A-Za-z0-9-]+", g)]


def ctd_download(out_path):
    """CTD 化学-基因互作表。BaP 的 MeSH accession 为 D001564。
    流程：带会话打开详情页 -> 从 HTML 中解析导出(CSV)链接 -> 下载。"""
    if os.path.exists(out_path):
        try:
            with open(out_path, encoding="utf-8", errors="replace") as f:
                head = f.read(3000)
            if "Gene Symbol" in head:
                print(f"  [CTD] 使用已下载文件（{os.path.basename(out_path)}）")
                return parse_ctd_file(out_path)
        except Exception:  # noqa: BLE001
            pass
    cj = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
    detail_url = "https://ctdbase.org/detail.go?type=chem&acc=D001564&view=gene"
    try:
        req = urllib.request.Request(detail_url, headers=UA)
        html = opener.open(req, timeout=120).read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        print(f"  [CTD] 自动下载失败（{e}）")
        print("  [CTD] 请手动打开 https://ctdbase.org/detail.go?type=chem&acc=D001564&view=gene")
        print("  [CTD] 点击页面下方 Download -> CSV，另存为 data/04_BaP_targets/CTD_BaP_gene_ixns.tsv 后重跑本脚本")
        return []
    links = []
    for link in re.findall(r'href="([^"]*6578706f7274=1[^"]*)"', html):
        links.append(link.replace("&amp;", "&"))
    csv_link = None
    for link in links:
        if re.search(r"[?&]d-\d+-e=1($|&)", link):
            csv_link = urllib.parse.urljoin(detail_url, link)
            break
    if not csv_link:
        print("  [CTD] 未从页面解析到 CSV 导出链接，请手动下载")
        return []
    try:
        req = urllib.request.Request(csv_link, headers=UA)
        text = opener.open(req, timeout=180).read().decode("utf-8", "replace")
    except Exception as e:  # noqa: BLE001
        print(f"  [CTD] 下载 CSV 失败（{e}），请手动下载")
        return []
    if "Gene Symbol" not in text[:3000] and "Input" not in text[:3000]:
        print("  [CTD] 返回内容不是 CSV 表格（可能仍需登录或反爬），请手动下载")
        return []
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    return parse_ctd_file(out_path)


# ---------------- 4. SEA / PharmMapper / SwissTargetPrediction 解析 ----------------
def discover_files(out_dir, patterns):
    """按多个文件名模式扫描目录，返回去重后的文件列表。"""
    found = []
    for pat in patterns:
        found.extend(glob.glob(os.path.join(out_dir, pat)))
    return sorted(set(found))


def detect_delimiter(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except Exception:  # noqa: BLE001
        return ","


def parse_sea(path):
    """解析 SEA 结果（兼容两种格式）：
    - SEA v3 网页下载：Target_Chembl_ID, Target_Species, Z_score, p(P_value),
      Target_Name, Uniprot_ID, Sequence
    - 旧格式：Target|GeneName|Description 的 CSV
    返回 (直接基因名 set, ChEMBL ID set, UniProt ID set)。
    """
    genes, chembl_ids, uniprot_ids = set(), set(), set()
    if not os.path.exists(path):
        return genes, chembl_ids, uniprot_ids
    delim = detect_delimiter(path)
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            return genes, chembl_ids, uniprot_ids
        fld = [x.strip() for x in reader.fieldnames]
        if any("Target_Chembl_ID" in x or "Uniprot_ID" in x for x in fld):
            for row in reader:
                r = {k.strip(): (v or "").strip() for k, v in row.items()}
                tid = r.get("Target_Chembl_ID", "")
                if re.fullmatch(r"CHEMBL\d+", tid):
                    chembl_ids.add(tid)
                up = r.get("Uniprot_ID", "")
                if re.fullmatch(r"[A-NP-RT-Z][0-9][A-Z0-9]{4}", up):
                    uniprot_ids.add(up)
                for m in re.finditer(r"\(([A-Z][A-Z0-9]{1,10})\)",
                                     r.get("Target_Name", "")):
                    genes.add(m.group(1))
        else:
            # 旧格式：首列可能形如 P02768|ALB|Serum albumin
            for row in reader:
                first = list(row.values())[0] if row else ""
                for token in re.split(r"[\|,;]", first):
                    token = token.strip()
                    if re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", token):
                        genes.add(token)
    return genes, chembl_ids, uniprot_ids


def parse_pharmmapper(path, z_threshold=0.5):
    """解析 PharmMapper 结果 CSV。
    输出文件常见列：Pharma Model,...,zscore,...,Name,Class,Uniplot,...
    Uniplot 形如 'BMP2_HUMAN'（基因_物种）。
    返回 (直接基因名 set, UniProt accession set)。
    """
    genes, uniprot_ids, mnemonics = set(), set(), set()
    if not os.path.exists(path):
        return genes, uniprot_ids, mnemonics
    lines = open(path, encoding="utf-8", errors="replace").readlines()
    start = None
    for i, line in enumerate(lines):
        if "Pharma Model" in line and ("Uniplot" in line or "UniProt" in line):
            start = i
            break
    if start is None:
        return genes, uniprot_ids, mnemonics
    reader = csv.DictReader(io.StringIO("".join(lines[start:])))
    for row in reader:
        r = {k.strip(): (v or "").strip() for k, v in row.items()}
        try:
            z = float(r.get("zscore", "0") or 0)
        except ValueError:
            z = 0.0
        if z < z_threshold:
            continue
        uni = r.get("Uniplot", "") or r.get("Target Name", "") or ""
        if re.fullmatch(r"[A-Z0-9]{2,10}_HUMAN", uni):
            mnemonics.add(uni)
        tok = uni.split("_")[0] if "_" in uni else uni
        if re.fullmatch(r"[A-NP-RT-Z][0-9][A-Z0-9]{4}", tok):
            uniprot_ids.add(tok)
        elif re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", tok):
            genes.add(tok)
    return genes, uniprot_ids, mnemonics


def uniprot_mnemonic_map(mnemonics):
    """UniProt 条目名（如 AK1C2_HUMAN）-> 正式基因名（AKR1C2）。"""
    mapping = {}
    mnemonics = sorted(set(mnemonics))
    if not mnemonics:
        return mapping
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
        except Exception as e:  # noqa: BLE001
            print(f"  [warn] UniProt 条目名批量查询失败（{e}），继续")
    return mapping


def parse_swisstarget(path):
    """解析 SwissTargetPrediction 的结果 CSV（替代/补充 SEA）。
    常见表头：Target, Common name, UniProt ID, ChEMBL ID, Probability*
    网站：http://www.swisstargetprediction.ch
    """
    genes = set()
    if not os.path.exists(path):
        return genes
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return genes
        name_col = None
        for fn in reader.fieldnames:
            if "common name" in fn.lower() or "gene" in fn.lower():
                name_col = fn
                break
        if name_col is None:
            name_col = reader.fieldnames[1] if len(reader.fieldnames) > 1 else reader.fieldnames[0]
        for row in reader:
            val = (row.get(name_col) or "").strip()
            # 单令牌基因符号，如 CYP1A1 / ABCB1 / MMP9
            if re.fullmatch(r"[A-Z][A-Z0-9]{1,10}", val):
                genes.add(val)
            # 名称里的括号符号，如 "Cytochrome P450 1A1 (CYP1A1)"
            for m in re.finditer(r"\(([A-Z][A-Z0-9]{1,10})\)", val):
                genes.add(m.group(1))
    return genes


def sea_ids_to_genes(chembl_ids, uniprot_ids):
    """把 SEA 返回的 CHEMBL/UniProt ID 映射为基因名。"""
    genes = set()
    for tid in sorted(chembl_ids):
        try:
            g, accs = chembl_target_genes(tid)
            genes.update(g)
        except Exception as e:  # noqa: BLE001
            print(f"  [warn] SEA target {tid} 映射失败: {e}")
    if uniprot_ids:
        mapping = uniprot_map(sorted(uniprot_ids))
        genes.update(v for v in mapping.values() if v)
    return genes


# ---------------- 主流程 ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None,
                    help="输出目录（默认 <pipeline>/data/04_BaP_targets）")
    ap.add_argument("--with-ctd", action="store_true",
                    help="是否纳入 CTD 化学-基因互作作为靶点来源（默认关闭，避免稀释）")
    ap.add_argument("--ctd-min-count", type=int, default=5,
                    help="CTD 基因的最低 Interaction Count（默认 5）")
    ap.add_argument("--pm-zscore", type=float, default=0.5,
                    help="PharmMapper 结果的 z-score 阈值（默认 0.5）")
    ap.add_argument("--refresh", action="store_true",
                    help="忽略 ChEMBL 缓存，重新从 API 拉取")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    out_dir = args.out or os.path.join(root, "data", "04_BaP_targets")
    os.makedirs(out_dir, exist_ok=True)

    print("== 1/5 PubChem ==")
    props = pubchem_props()
    print("  SMILES:", props["smiles"])
    print("  InChIKey:", props["inchikey"])

    print("== 2/5 ChEMBL ==")
    mol = chembl_compound(props["inchikey"])
    mol_id = mol.get("molecule_chembl_id")
    chembl_records = []
    cache = load_chembl_cache(out_dir) if not args.refresh else None
    if cache is not None:
        chembl_records = cache
        print("  使用 ChEMBL 缓存:", len(chembl_records), "条")
    elif mol_id:
        print("  molecule_chembl_id:", mol_id)
        acts = chembl_bioactivities(mol_id)
        print("  bioactivity records:", len(acts))
        tids = sorted({a["target_chembl_id"] for a in acts if a.get("target_chembl_id")})
        for tid in tids:
            try:
                genes, accs = chembl_target_genes(tid)
            except Exception as e:  # noqa: BLE001
                print(f"  [warn] target {tid} 解析失败: {e}")
                continue
            mapping = uniprot_map(accs) if accs else {}
            for acc, gene in mapping.items():
                if gene:
                    genes.append(gene)
        for g in sorted(set(genes)):
            chembl_records.append({"gene_symbol": g, "source": "ChEMBL",
                                   "target_id": tid})
        print("  ChEMBL unique genes:", len({r["gene_symbol"] for r in chembl_records}))
        save_chembl_cache(out_dir, chembl_records)
    else:
        print("  [warn] ChEMBL 未收录该 InChIKey，跳过 ChEMBL 靶点")

    print("== 3/5 CTD（可选，默认关闭）==")
    ctd_records = []
    if args.with_ctd:
        ctd_path = os.path.join(out_dir, "CTD_BaP_gene_ixns.tsv")
        ctd_gene_n = ctd_download(ctd_path)
        ctd_genes = sorted({g for g, n in ctd_gene_n if n >= args.ctd_min_count})
        ctd_records = [{"gene_symbol": g, "source": "CTD", "target_id": ""}
                       for g in ctd_genes]
        print(f"  CTD genes (Interaction Count >= {args.ctd_min_count}):", len(ctd_records))
    else:
        print("  未启用。如需纳入 CTD，加 --with-ctd（并可用 --ctd-min-count 控制置信度）。")

    print("== 4/5 SEA / PharmMapper / SwissTargetPrediction（自动扫描 data/04_BaP_targets/）==")
    sea_files = discover_files(out_dir, [
        "sea_result*.csv", "sea_result*.tsv",
        "sea_results.csv", "sea_results.tsv"])
    sea_genes = set()
    sea_chembl, sea_uniprot = set(), set()
    for sf in sea_files:
        g, c, u = parse_sea(sf)
        sea_genes |= g
        sea_chembl |= c
        sea_uniprot |= u
    sea_genes |= sea_ids_to_genes(sea_chembl, sea_uniprot)
    if sea_files:
        print("  找到 SEA 结果文件:", [os.path.basename(x) for x in sea_files],
              "| SEA genes:", len(sea_genes))
    else:
        print("  未找到 SEA 结果文件（sea_result*.csv/tsv）")

    pm_files = discover_files(out_dir, ["pharmmapper*"])
    pm_genes, pm_accs, pm_mnemonics = set(), set(), set()
    for pf in pm_files:
        g, a, m = parse_pharmmapper(pf, z_threshold=args.pm_zscore)
        pm_genes |= g
        pm_accs |= a
        pm_mnemonics |= m
    if pm_mnemonics:
        pm_genes |= set(v for v in uniprot_mnemonic_map(pm_mnemonics).values() if v)
    if pm_accs:
        pm_map = uniprot_map(sorted(pm_accs))
        pm_genes |= set(v for v in pm_map.values() if v)
    if pm_files:
        print("  找到 PharmMapper 结果文件:", [os.path.basename(x) for x in pm_files],
              "| PharmMapper genes:", len(pm_genes))
    else:
        print("  未找到 PharmMapper 结果文件（pharmmapper*）")

    stp_files = discover_files(out_dir, [
        "swisstarget_results.csv", "*wiss*arget*.csv", "swiss*.csv"])
    stp_genes = set()
    for tf in stp_files:
        stp_genes |= parse_swisstarget(tf)
    if stp_files:
        print("  找到 SwissTargetPrediction 文件:",
              [os.path.basename(x) for x in stp_files],
              "| genes:", len(stp_genes))
    else:
        print("  未找到 SwissTargetPrediction 文件（swisstarget_results.csv 或 SwissTargetPrediction .csv）")
    sea_records = [{"gene_symbol": g, "source": "SEA", "target_id": ""}
                   for g in sorted(sea_genes)]
    pm_records = [{"gene_symbol": g, "source": "PharmMapper", "target_id": ""}
                  for g in sorted(pm_genes)]
    stp_records = [{"gene_symbol": g, "source": "SwissTargetPrediction", "target_id": ""}
                   for g in sorted(stp_genes)]
    print("  SEA genes:", len(sea_genes), "| PharmMapper genes:", len(pm_genes),
          "| SwissTargetPrediction genes:", len(stp_genes))

    print("== 5/5 合并并集 ==")
    all_records = chembl_records + ctd_records + sea_records + pm_records + stp_records
    union = {}
    for rec in all_records:
        g = rec["gene_symbol"].upper()
        if g in union:
            union[g]["source"] = union[g]["source"] + "," + rec["source"]
            if rec["target_id"] and rec["target_id"] not in union[g]["target_id"]:
                union[g]["target_id"] = (union[g]["target_id"] + ";" + rec["target_id"])
        else:
            union[g] = dict(rec, gene_symbol=g)
    rows = sorted(union.values(), key=lambda r: r["gene_symbol"])
    out_tsv = os.path.join(out_dir, "union_targets.tsv")
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_symbol", "source", "target_id"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "smiles": props["smiles"],
        "inchikey": props["inchikey"],
        "n_chembl": len({r["gene_symbol"] for r in chembl_records}),
        "n_ctd": len({r["gene_symbol"] for r in ctd_records}),
        "n_sea": len(sea_genes),
        "n_pharmmapper": len(pm_genes),
        "n_swisstarget": len(stp_genes),
        "n_union": len(rows),
    }
    with open(os.path.join(out_dir, "target_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("  union targets:", len(rows))
    print("  ->", out_tsv)
    print("  提示：SEA 对 BaP 这类平面多环芳烃可能返回空结果（相似配体少），属正常现象；"
          "可用 SwissTargetPrediction 结果（swisstarget_results.csv）替代或补充。"
          "默认并集 = ChEMBL + SEA + PharmMapper + SwissTargetPrediction（与模板论文思路一致）；"
          "若开启 --with-ctd，会按 Interaction Count 阈值加入 CTD 基因。")


if __name__ == "__main__":
    main()
    sys.exit(0)
