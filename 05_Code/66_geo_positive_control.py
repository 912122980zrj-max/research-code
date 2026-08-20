# ============================================================
# 66_geo_positive_control.py —— B-6：GEO 正对照检索
# 依据：修改报告_BaP_TBI_最终版.md 阶段 B-6
#   4h 检索（GEO "benzo[a]pyrene" + Mus musculus + 表达芯片/RNA-seq，
#   脑相关优先，每组 >=3 例）；超时转 CTD 弱替代。
# 输出：results/16_Orthology/geo_positive_control_search.csv
# ============================================================
import csv, json, os, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET

OUT = os.path.join(os.getcwd(), "results", "16_Orthology")
os.makedirs(OUT, exist_ok=True)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "codex-research@example.invalid"

def eutil(endpoint, params, retries=5):
    params["tool"] = "codex-tbi"
    params["email"] = EMAIL
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "codex/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(1.5 + i)

def esearch(term, retmax=200):
    xml = eutil("esearch.fcgi", {"db": "gds", "term": term, "retmax": retmax,
                                 "retmode": "xml"})
    root = ET.fromstring(xml)
    ids = [x.text for x in root.findall(".//IdList/Id")]
    count = root.findtext(".//Count")
    return ids, int(count) if count else 0

def esummary(ids):
    xml = eutil("esummary.fcgi", {"db": "gds", "id": ",".join(ids),
                                  "retmode": "xml"})
    root = ET.fromstring(xml)
    out = []
    for doc in root.findall(".//DocSum"):
        rec = {"uid": doc.findtext("Id")}
        for item in doc.findall("Item"):
            rec[item.get("Name")] = "".join(item.itertext())
        out.append(rec)
    return out

def run():
    log = os.path.join(OUT, "geo_positive_control.log")
    with open(log, "w", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] search started\n")

    queries = {
        "mouse_bap_all": 'benzo[a]pyrene[All Fields] AND Mus musculus[Organism] AND gse[ETYP]',
        "mouse_bap_expr": ('benzo[a]pyrene[All Fields] AND Mus musculus[Organism] AND gse[ETYP] '
                           'AND ("expression profiling by array"[Filter] OR '
                           '"expression profiling by high throughput sequencing"[Filter])'),
        "mouse_bap_brain": 'benzo[a]pyrene[All Fields] AND Mus musculus[Organism] AND gse[ETYP] AND brain[All Fields]',
    }
    seen = {}
    summary = {}
    for name, term in queries.items():
        try:
            ids, count = esearch(term)
        except Exception as e:
            with open(log, "a", encoding="utf-8") as f:
                f.write(f"[{time.strftime('%H:%M:%S')}] {name} FAIL {e!r}\n")
            continue
        with open(log, "a", encoding="utf-8") as f:
            f.write(f"[{time.strftime('%H:%M:%S')}] {name}: count={count} ids={len(ids)}\n")
        if not ids:
            continue
        # 分批 esummary（每批 100）
        for i in range(0, len(ids), 100):
            batch = ids[i:i+100]
            try:
                recs = esummary(batch)
            except Exception as e:
                with open(log, "a", encoding="utf-8") as f:
                    f.write(f"  batch {i} FAIL {e!r}\n")
                continue
            for r in recs:
                uid = r.get("uid")
                summary[uid] = r
                seen.setdefault(uid, []).append(name)
        time.sleep(0.5)

    rows = []
    for uid, r in summary.items():
        title = r.get("title", "")
        acc = r.get("accession", "")
        taxon = r.get("taxon", "")
        gds_type = r.get("gdsType", "")
        n_samples = r.get("n_samples", "")
        platform = r.get("platform", "")
        gse = ""
        for item in r.get("GDS", "").replace(",", " ").split():
            if item.startswith("GSE"):
                gse = item
                break
        text = (title + " " + r.get("summary", "")).lower()
        brain = any(k in text for k in ["brain", "cortex", "hippocamp", "cerebr",
                                        "striatum", "neuron", "neuro", "cerebellum"])
        rows.append({
            "uid": uid, "accession": acc, "gse": gse, "title": title,
            "taxon": taxon, "type": gds_type, "n_samples": n_samples,
            "platform": platform, "brain_related": "yes" if brain else "",
            "query_hits": ";".join(seen[uid]),
            "summary": r.get("summary", "")
        })

    rows.sort(key=lambda x: (x["brain_related"] != "yes", x["accession"]))
    with open(os.path.join(OUT, "geo_positive_control_search.csv"), "w",
              encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["uid", "accession", "gse", "title", "taxon", "type",
                            "n_samples", "platform", "brain_related", "query_hits",
                            "summary"])
        w.writeheader()
        w.writerows(rows)
    with open(log, "a", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] done; unique series={len(rows)}\n")
    print("unique hits:", len(rows), "->", os.path.join(OUT, "geo_positive_control_search.csv"))

if __name__ == "__main__":
    run()
