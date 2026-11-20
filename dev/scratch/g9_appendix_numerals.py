#!/usr/bin/env python3
"""Gate G9 of TASK_estimated_weights_appendix.md: every numeral in the prose of the new article and in the
pointer sentences added to four files (lines added or changed between a15ff33 and HEAD), each with a class.

Prose = lines outside fenced code blocks. Before the numerals are read off a line:
  - inline R (`r ...`) is replaced by a marker: its value is computed when the file is knitted;
  - citation keys (@key2004) are replaced by a marker: they are not numerals.
Lines of the copied evaluation table (table rows between "**Evaluation evidence.**" and "Source:") are counted,
not listed numeral by numeral. Exit status 1 when a numeral in plain text has no class.

Usage: python3 dev/scratch/g9_appendix_numerals.py > dev/audit/estimated_weights_appendix/g9_numerals.txt
"""
import re
import subprocess
import sys

BASE = "a15ff33"
FILES = ["vignettes/articles/estimated_weights_methodology.Rmd",
         "vignettes/weightedsurv_examples.Rmd",
         "vignettes/articles/weightedkm_methodology.Rmd",
         "vignettes/articles/weightedcox_methodology.Rmd",
         "vignettes/articles/weightedrmst_methodology.Rmd"]
NUM = re.compile(r"(?<![A-Za-z_])\d+(?:[.,]\d+)*%?")
INLINE = "<inline R>"
CITE = "<citation>"

# (pattern matched against a window of text around the numeral, class). First match wins.
RULES = [
    (r"toc_depth|VignetteEncoding", "YAML header of the file, not prose"),
    (r"FH\(\d, ?\d\)", "design constant: name of a weighting scheme"),
    (r"section \d|sections \d|Section \d|Lemma \d|Lemmas \d|Stage 2e|study M7|Table M7", "reference to a numbered section, lemma, stage or table"),
    (r"nominal 0\.95|two-sided 5%|month 18", "copied evaluation table: its design (nominal level, test level, time point), as in the source"),
    (r"2,000 replicates|have 500|first 500|same first 500", "copied evaluation table: number of replicates, as in the source"),
    (r"unless `draws` is 0|a term is 0 when|numerically 0 or 1|away from 0 and 1|other subjects|\(0 for the others\)", "formula constant: the value 0 or 1 in a condition"),
    (r"ratio near 1|variance 1", "formula constant: reference value"),
    (r"Three facts|at one point only", "count in the text, not a result"),
]

TABLE_START = "**Evaluation evidence.**"
TABLE_END = "Source: `dev/audit/REPORT_km_fixes_stage2e.md`"

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
    in_tab = False
    scan.table_rows = getattr(scan, "table_rows", 0)
    scan.table_nums = getattr(scan, "table_nums", 0)
    for i, ln in enumerate(lines, 1):
        if ln.startswith(TABLE_START):
            in_tab = True
        if ln.startswith(TABLE_END):
            in_tab = False
        if i not in add:
            continue
        if in_tab and ln.startswith("|"):
            scan.table_rows += 1
            scan.table_nums += len(NUM.findall(ln))
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
    print(f"\nCopied evaluation table: {scan.table_rows} table lines holding {scan.table_nums} numerals (source: "
          f"dev/audit/REPORT_km_fixes_stage2e.md, section C); not listed one by one.")
    print(f"UNCLASSIFIED numerals in plain text: {bad}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
