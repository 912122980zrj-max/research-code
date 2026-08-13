# ============================================================
# 61b_fetch_ensembl_orthologs.py —— B-1 数据抓取（Ensembl REST）
# 依据：修改报告_BaP_TBI_最终版.md 阶段 B-1
# 用法：python scripts/61b_fetch_ensembl_orthologs.py
# 说明：网络不佳时可反复重跑；已抓取条目写入缓存，断点续传。
# 输出：
#   results/16_Orthology/ensembl_lookup_cache.json   （symbol->ENSG）
#   results/16_Orthology/ensembl_homology_cache.json （ENSG->orthologs）
#   results/16_Orthology/ensembl_gene_cache.json     （ENSMUSG->symbol）
#   results/16_Orthology/orthology_raw_ensembl.tsv
# ============================================================
import csv, json, os, sys, threading, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.getcwd()
OUT = os.path.join(ROOT, "results", "16_Orthology")
os.makedirs(OUT, exist_ok=True)

TIER_CSV = os.path.join(OUT, "bap_targets_364_tier.csv")
LOOKUP_CACHE = os.path.join(OUT, "ensembl_lookup_cache.json")
HOM_CACHE = os.path.join(OUT, "ensembl_homology_cache.json")
GENE_CACHE = os.path.join(OUT, "ensembl_gene_cache.json")
RAW_TSV = os.path.join(OUT, "orthology_raw_ensembl.tsv")

BASE = "https://rest.ensembl.org"
MAX_RETRY = 5
TIMEOUT = 45
WORKERS = 8

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(OUT, "ensembl_fetch.log"), "a", encoding="utf-8") as f:
        f.write(line + "\n")

def http_json(path, data=None, retries=MAX_RETRY):
    url = BASE + path
    last = None
    for i in range(retries):
        try:
            headers = {"Content-Type": "application/json"}
            body = json.dumps(data).encode() if data is not None else None
            req = urllib.request.Request(url, data=body, headers=headers)
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return json.loads(r.read().decode("utf-8", errors="replace"))
        except Exception as e:
            last = e
            time.sleep(0.5 + i * 0.7)
    raise last

def load_cache(path):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {}

def save_cache(path, obj):
    with _lock:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(obj, f, ensure_ascii=False, indent=1)

_lock = threading.Lock()

def main():
    if os.path.exists(os.path.join(OUT, "ensembl_fetch.log")):
        os.remove(os.path.join(OUT, "ensembl_fetch.log"))
    log("Ensembl fetch started")

    # release 信息
    try:
        rel = http_json("/info/software/version/?content-type=application/json",
                        retries=3)
        release = str(rel.get("release", ""))
        log("Ensembl release: " + release)
    except Exception as e:
        release = "unknown"
        log("release fetch failed: " + repr(e))

    symbols = []
    with open(TIER_CSV, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            symbols.append(row["gene_symbol"])
    log(f"target genes: {len(symbols)}")

    lookup_cache = load_cache(LOOKUP_CACHE)
    hom_cache = load_cache(HOM_CACHE)
    gene_cache = load_cache(GENE_CACHE)

    def lookup_symbol(sym):
        if sym in lookup_cache:
            return lookup_cache[sym]
        try:
            obj = http_json("/lookup/symbol/homo_sapiens/"
                            + urllib.parse.quote(sym) + "?content-type=application/json")
            gid = obj.get("id") if isinstance(obj, dict) else None
        except Exception:
            gid = None
        lookup_cache[sym] = gid
        save_cache(LOOKUP_CACHE, lookup_cache)
        return gid

    def fetch_homology(ensg):
        if ensg in hom_cache:
            return hom_cache[ensg]
        try:
            obj = http_json("/homology/id/human/" + ensg
                            + "?target_species=mouse&type=orthologues&content-type=application/json")
            data = obj.get("data", []) if isinstance(obj, dict) else []
            homs = []
            for d in data:
                for h in d.get("homologies", []):
                    tgt = h.get("target", {})
                    homs.append({
                        "type": h.get("type"),
                        "mouse_id": tgt.get("id"),
                        "perc_id": tgt.get("perc_id"),
                        "method_link_type": h.get("method_link_type")
                    })
            hom_cache[ensg] = homs
        except Exception:
            hom_cache[ensg] = None
        save_cache(HOM_CACHE, hom_cache)
        return hom_cache[ensg]

    def lookup_gene(ensmusg):
        if ensmusg in gene_cache:
            return gene_cache[ensmusg]
        try:
            obj = http_json("/lookup/id/" + ensmusg + "?content-type=application/json")
            gene_cache[ensmusg] = obj.get("display_name") if isinstance(obj, dict) else None
        except Exception:
            gene_cache[ensmusg] = None
        save_cache(GENE_CACHE, gene_cache)
        return gene_cache[ensmusg]

    def process_symbol(sym):
        ensg = lookup_symbol(sym)
        if not ensg:
            return [{"human_symbol": sym, "human_ensembl": "",
                     "mouse_ensembl": "", "mouse_symbol": "",
                     "orthology_type": "", "perc_id": "",
                     "ensembl_release": release}]
        homs = fetch_homology(ensg)
        if not homs:
            return [{"human_symbol": sym, "human_ensembl": ensg,
                     "mouse_ensembl": "", "mouse_symbol": "",
                     "orthology_type": "", "perc_id": "",
                     "ensembl_release": release}]
        out = []
        for h in homs:
            msym = lookup_gene(h["mouse_id"]) if h["mouse_id"] else None
            out.append({"human_symbol": sym, "human_ensembl": ensg,
                        "mouse_ensembl": h["mouse_id"] or "",
                        "mouse_symbol": msym or "",
                        "orthology_type": h["type"] or "",
                        "perc_id": h["perc_id"] or "",
                        "ensembl_release": release})
        return out

    rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futures = {ex.submit(process_symbol, s): s for s in symbols}
        for fut in as_completed(futures):
            rows.extend(fut.result())
            done += 1
            if done % 20 == 0:
                log(f"processed {done}/{len(symbols)}")

    with open(RAW_TSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "human_symbol", "human_ensembl", "mouse_ensembl",
            "mouse_symbol", "orthology_type", "perc_id", "ensembl_release"],
            delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"done; rows={len(rows)} -> {RAW_TSV}")

if __name__ == "__main__":
    main()
