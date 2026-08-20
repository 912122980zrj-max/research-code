# ============================================================
# 78_renumber_references.py
# 按 BMC 投稿稿（MANUSCRIPT_v7_BMC_EN.md）中参考文献的首次出现顺序
# 重新编号（旧编号 1-39 去掉未使用的 24 后 → 新编号 1-38），并同步：
#   - EN/ZH v7 主稿（正文引用 + 文末参考文献列表）
#   - 补充方法、补充图例（02_Main_Figures/FIGURE_LEGENDS.md）
#   - 回复审稿人信（RESPONSE_TO_REVIEWERS.md）
# 输出：编号映射表（控制台 + reports/reference_renumber_map.csv）
# ============================================================
import csv
import os
import re

ROOT = r"E:\sheng xin\tbi"

BMC = os.path.join(ROOT, "submission", "01_Manuscript", "MANUSCRIPT_v7_BMC_EN.md")
EN = os.path.join(ROOT, "reports", "manuscript", "MANUSCRIPT_FINAL_v7_EN.md")
ZH = os.path.join(ROOT, "reports", "manuscript", "MANUSCRIPT_FINAL_v7_ZH.md")
SUPP_METHODS = os.path.join(
    ROOT, "submission", "03_Supplementary", "Text", "SUPPLEMENTARY_METHODS_RESULTS_EN.md"
)
FIG_LEGENDS = os.path.join(ROOT, "submission", "02_Main_Figures", "FIGURE_LEGENDS.md")
RESPONSE = os.path.join(ROOT, "submission", "01_Manuscript", "RESPONSE_TO_REVIEWERS.md")

CIT_RE = re.compile(r"\[([0-9]+(?:[-,][0-9]+)*)\]")
REF_NUM_RE = re.compile(r"^(\d+)\.\s")          # EN/BMC style: "12. ..."
REF_BRACKET_RE = re.compile(r"^\[(\d+)\]\s*")   # ZH style: "[12] ..."


def expand(token):
    nums = []
    for part in token.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            nums.extend(range(int(a), int(b) + 1))
        else:
            nums.append(int(part))
    return nums


def fmt(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 1:
            out.append(f"{nums[i]}-{nums[j]}")
        else:
            out.append(str(nums[i]))
        i = j + 1
    return "[" + ",".join(out) + "]"


def rewrite_citations(text, mapping):
    def repl(m):
        return fmt([mapping[n] for n in expand(m.group(1))])
    return CIT_RE.sub(repl, text)


def reorder_references(text, mapping, bracket_style=False):
    """Reorder the numbered reference-list entries under References heading."""
    lines = text.splitlines(keepends=True)
    in_refs = False
    entries = []          # (old_no, line)
    body, tail = [], []
    for ln in lines:
        if ln.startswith("## References") or ln.startswith("## 参考文献"):
            in_refs = True
            body.append(ln)
            continue
        if in_refs:
            m = (REF_BRACKET_RE if bracket_style else REF_NUM_RE).match(ln)
            if m:
                entries.append((int(m.group(1)), ln))
            elif ln.startswith("#") or ln.startswith("##"):
                in_refs = False
                tail.append(ln)
            else:
                tail.append(ln)   # software-citation note and blank lines
        else:
            body.append(ln)
    if not entries:
        raise RuntimeError("No reference entries found")
    new_entries = sorted(
        ((mapping[old], ln) for old, ln in entries if old in mapping),
        key=lambda x: x[0],
    )
    for new_no, ln in new_entries:
        body.append(f"{new_no}. {ln[len(re.match(r'^(\d+)\.\s*', ln).group(0)):]}"
                    if not bracket_style else
                    f"[{new_no}] {ln[len(re.match(r'^\[\d+\]\s*', ln).group(0)):]}")
    return "".join(body + tail)


def main():
    bmc = open(BMC, encoding="utf-8").read()

    # ---- build mapping from BMC first-appearance order ----
    mapping = {}
    for m in CIT_RE.finditer(bmc):
        for n in expand(m.group(1)):
            if n not in mapping:
                mapping[n] = len(mapping) + 1
    expected_old = set(range(1, 40)) - {24}
    if set(mapping) != expected_old:
        raise SystemExit(
            "Unexpected old-number set: missing=%s extra=%s"
            % (sorted(expected_old - set(mapping)), sorted(set(mapping) - expected_old))
        )
    if sorted(mapping.values()) != list(range(1, 39)):
        raise SystemExit("New numbering is not 1..38")

    print("old -> new")
    for old in sorted(mapping):
        print(f"  {old:>2} -> {mapping[old]:>2}")
    with open(os.path.join(ROOT, "reports", "reference_renumber_map.csv"),
              "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["old_number", "new_number"])
        for old in sorted(mapping):
            w.writerow([old, mapping[old]])

    # ---- BMC manuscript: rewrite citations + reorder list ----
    bmc_new = reorder_references(rewrite_citations(bmc, mapping), mapping)
    open(BMC, "w", encoding="utf-8").write(bmc_new)

    # ---- EN master ----
    en = open(EN, encoding="utf-8").read()
    en_new = reorder_references(rewrite_citations(en, mapping), mapping)
    open(EN, "w", encoding="utf-8").write(en_new)

    # ---- ZH master (bracketed list style) ----
    zh = open(ZH, encoding="utf-8").read()
    zh_new = reorder_references(rewrite_citations(zh, mapping), mapping,
                                bracket_style=True)
    open(ZH, "w", encoding="utf-8").write(zh_new)

    # ---- supplementary methods / figure legends / response letter ----
    for path in (SUPP_METHODS, FIG_LEGENDS, RESPONSE):
        txt = open(path, encoding="utf-8").read()
        open(path, "w", encoding="utf-8").write(rewrite_citations(txt, mapping))

    print("\nDone. Files updated.")


if __name__ == "__main__":
    main()
