#!/usr/bin/env python3
"""Gate G11 of TASK_vignette_revision.md: every numeral in the prose lines added or changed between
ac76f82 and HEAD, each with a class.

Prose = lines outside fenced code blocks. Before the numerals are read off a line:
  - inline R (`r ...`) is replaced by a marker: its value is computed when the file is knitted;
  - citation keys (@key2004) are replaced by a marker: they are not numerals.
What is left is listed, with where it stands: in plain text, in a code span (`...`) or in math ($...$ or a
$$ block). Comment lines added inside R chunks that contain a numeral are listed in a second block.
Exit status 1 when a numeral in plain text has no class.

Usage: python3 dev/scratch/g11_numerals.py > dev/audit/vignette_revision/g11_numerals.txt
"""
import re
import subprocess
import sys

BASE = "ac76f82"
FILES = ["vignettes/weightedsurv_examples.Rmd",
         "vignettes/articles/martingale_resampling_applications.Rmd",
         "vignettes/articles/weightedcox_methodology.Rmd",
         "vignettes/articles/weightedkm_methodology.Rmd",
         "vignettes/articles/weightedrmst_methodology.Rmd",
         "README.md"]
NUM = re.compile(r"(?<![A-Za-z_])\d+(?:[.,]\d+)*%?")
INLINE = "<inline R>"
CITE = "<citation>"

# (pattern matched against a window of text around the numeral, class). First match wins.
RULES = [
    (r"FH\(\d, ?\d\)|MB\(12\)|zero_one\(6\)", "design constant: name of a weighting scheme"),
    (r"first 20 resampled", "design constant: number of curves drawn by show_resamples (fixed in KM_diff())"),
    (r"95% quantile|pointwise 95%|0\.95 quantile", "design constant: confidence level"),
    (r"\(1000 draws\)", "design constant: draws = 1000 in the chunk that follows"),
    (r"draws = 5000|around 500 for", "design constant: draws argument (sentence of ac76f82; one clause added)"),
    (r"had 10 events", "design constant: the starting rule of item 4 (10th event in each arm)"),
    (r"at 48 months|at 72 months", "design constant: time points chosen for the at-risk counts (the counts are inline-computed)"),
    (r"coded 0/1", "design constant: coding of the treatment"),
    (r"ratio near 1|close to 1", "design constant: reference value of the ESS ratio"),
    (r"Step 4|steps 1--3|step 2|step 3", "reference to a numbered step of the text"),
    (r"must be 0\)", "design constant: argument value (draws must be 0)"),
    (r"Grade 3 tumours", "design constant: name of the grade category (the percentages are inline-computed)"),
    (r"a term is 0 when", "formula constant: the zero term of the formula (wording of ac76f82)"),
    (r"a weight of 2 is not two copies", "illustration, not a result (sentence of ac76f82)"),
]


def added_lines(path):
    out = subprocess.run(["git", "diff", "-U0", BASE, "HEAD", "--", path], capture_output=True, text=True).stdout
    res, n = set(), None
    for ln in out.splitlines():
        m = re.match(r"@@ -\S+ \+(\d+)(?:,(\d+))? @@", ln)
        if m:
            n = int(m.group(1))
            continue
        if n is None or ln.startswith("+++") or ln.startswith("---"):
            continue
        if ln.startswith("+"):
            res.add(n)
            n += 1
    return res


def classify(window):
    for pat, cls in RULES:
        if re.search(pat, window):
            return cls
    return "UNCLASSIFIED"


def scan(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    add = added_lines(path)
    # lines inside fenced code, and lines inside a display-math block that spans several lines
    fenced_line, math_line = set(), set()
    fenced = open_math = False
    for i, ln in enumerate(lines, 1):
        if ln.lstrip().startswith("```"):
            fenced = not fenced
            fenced_line.add(i)
            continue
        if fenced:
            fenced_line.add(i)
            continue
        k = ln.count("$$")
        if open_math or k % 2 == 1:
            math_line.add(i)
        if k % 2 == 1:
            open_math = not open_math
    prose, comments = [], []
    for i, ln in enumerate(lines, 1):
        if i not in add:
            continue
        if i in fenced_line:
            if ln.lstrip().startswith("#") and NUM.search(ln):
                comments.append((i, ", ".join(NUM.findall(ln)), ln.strip()))
            continue
        txt = re.sub(r"`r [^`]*`", INLINE, ln)
        txt = re.sub(r"@[A-Za-z][A-Za-z0-9_]*", CITE, txt)
        if i in math_line:
            segs = [("math", 0, txt)]
        else:
            segs, pos = [], 0
            for m in re.finditer(r"`[^`]*`|\$\$.*?\$\$|\$[^$]*\$", txt):
                segs.append(("plain", pos, txt[pos:m.start()]))
                segs.append(("code" if m.group().startswith("`") else "math", m.start(), m.group()))
                pos = m.end()
            segs.append(("plain", pos, txt[pos:]))
        for where, off, seg in segs:
            for m in NUM.finditer(seg):
                tok = m.group()
                a, b = off + m.start(), off + m.end()
                window = txt[max(0, a - 45):b + 45]
                if where == "math":
                    cls = "formula constant: index, exponent or tolerance of a formula"
                elif where == "code":
                    cls = "design constant: argument value or name in a code span"
                elif re.match(r"\s*\d+\.\s", txt) and a < 6:
                    cls = "list-item number"
                else:
                    cls = classify(window)
                prose.append((i, where, tok, cls, window.strip()))
    return len(add), prose, comments


def main():
    bad = 0
    for path in FILES:
        n_add, prose, comments = scan(path)
        plain = [r for r in prose if r[1] == "plain"]
        print(f"\n== {path}: {n_add} added or changed lines; numerals literally in the added prose: "
              f"{len(prose)} ({len(plain)} in plain text, {sum(r[1] == 'code' for r in prose)} in code spans, "
              f"{sum(r[1] == 'math' for r in prose)} in math)")
        for i, where, tok, cls, window in prose:
            if cls == "UNCLASSIFIED":
                bad += 1
            print(f"  L{i:<4} {tok:<7} [{where}] {cls}  | ...{window}...")
        if comments:
            print("  -- comment lines added inside R chunks that contain a numeral (not prose):")
            for i, toks, ln in comments:
                print(f"  L{i:<4} {toks:<10} | {ln[:130]}")
    print(f"\nUNCLASSIFIED numerals in plain text: {bad}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
