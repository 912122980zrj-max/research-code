# ============================================================
# 82_reference_doi_authors.py —— m10 应对方针
# 1) BMC 稿 38 条参考文献按 REFERENCES_VERIFIED.md 追加 DOI；
# 2) 三份手稿中 [8]（Goal 等）作者列表补全为 4 人。
# 映射：reference_renumber_map.csv（旧->新）。
# ============================================================
import csv
import os
import re

ROOT = r"E:\sheng xin\tbi"
OLD_DOI = {
    1: "10.1016/j.envint.2026.110390",
    2: "10.1007/s10565-025-10109-4",
    3: "10.1016/j.celrep.2023.113395",
    4: "10.1016/j.neuron.2024.06.021",
    5: "10.1038/s41586-025-09449-2",
    6: "10.3389/fneur.2023.1087767",
    7: "10.1016/j.intimp.2019.105980",
    8: "10.1016/j.crneur.2023.100122",
    9: "10.1016/j.etp.2015.11.007",
    10: "10.4103/1673-5374.284996",
    11: "10.1186/1471-2105-9-559",
    12: "10.1093/nar/gkv007",
    13: "10.1038/nmeth.3337",
    14: "10.1038/s41467-021-21246-9",
    15: "10.1016/j.patter.2022.100434",
    16: "10.1016/j.cell.2021.04.048",
    17: "10.1038/s41592-019-0619-0",
    18: "10.1093/nar/gkx374",
    19: "10.1038/nbt1284",
    20: "10.1038/s41586-021-03819-2",
    21: "10.1021/acs.jcim.1c00203",
    22: "10.1093/nar/gkac833",
    23: "10.7326/M18-1376",
    25: "10.1093/nar/gky1075",
    26: "10.1093/nar/gkz382",
    27: "10.1038/s41587-019-0114-2",
    28: "10.1038/s41467-018-08023-x",
    29: "10.1038/s41467-020-15816-6",
    30: "10.1016/j.xinn.2021.100141",
    31: "10.1093/biostatistics/kxj037",
    32: "10.12688/f1000research.73600.2",
    33: "10.1136/bmj-2024-082505",
    34: "10.1186/1471-2105-12-77",
    35: "10.18637/jss.v028.i05",
    36: "10.18637/jss.v033.i01",
    37: "10.1145/2939672.2939785",
    38: "10.1007/978-0-387-98141-3",
}

old2new = {}
with open(os.path.join(ROOT, "reports", "reference_renumber_map.csv"),
          encoding="utf-8") as f:
    r = csv.reader(f)
    next(r)
    for old, new in r:
        old2new[int(old)] = int(new)
new2old = {v: k for k, v in old2new.items()}

BMC = os.path.join(ROOT, "submission", "01_Manuscript", "MANUSCRIPT_v7_BMC_EN.md")
lines = open(BMC, encoding="utf-8").read().splitlines()
out = []
in_refs = False
for ln in lines:
    if ln.startswith("## References"):
        in_refs = True
        out.append(ln)
        continue
    if in_refs:
        m = re.match(r"^(\d+)\.\s(.*)$", ln)
        if m:
            n = int(m.group(1))
            old = new2old[n]
            doi = OLD_DOI.get(old)
            body = m.group(2)
            if "https://doi.org/" not in body and doi:
                body = body.rstrip() + " https://doi.org/" + doi
            if old == 8:
                body = re.sub(r"Goal A, Raj K, Singh S, et al\.",
                              "Goal A, Raj K, Singh S, Arora R.", body)
            out.append(f"{n}. {body}")
            continue
        if ln.startswith("##"):
            in_refs = False
    out.append(ln)
open(BMC, "w", encoding="utf-8", newline="\n").write("\n".join(out) + "\n")

# EN/ZH masters: fix Goal authors
for path, pat in [
    (os.path.join(ROOT, "reports", "manuscript", "MANUSCRIPT_FINAL_v7_EN.md"),
     r"Goal A, Raj K, Singh S, et al\."),
    (os.path.join(ROOT, "reports", "manuscript", "MANUSCRIPT_FINAL_v7_ZH.md"),
     r"Goal A, Raj K, Singh S, et al\."),
]:
    txt = open(path, encoding="utf-8").read()
    txt2 = re.sub(pat, "Goal A, Raj K, Singh S, Arora R.", txt)
    open(path, "w", encoding="utf-8", newline="\n").write(txt2)
    print("fixed authors in", os.path.basename(path),
          "changed:", txt != txt2)
print("DOI entries written to BMC manuscript")
