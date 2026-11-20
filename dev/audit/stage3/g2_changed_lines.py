# G2: classify every changed line of the two Rmd files (e95728a..HEAD) as prose (outside a chunk),
# chunk header, or chunk code. usage: python3 g2_changed_lines.py <base-sha>
import subprocess, sys, re
base = sys.argv[1]
for f in ["vignettes/weightedsurv_examples.Rmd", "vignettes/articles/martingale_resampling_applications.Rmd"]:
    d = subprocess.run(["git", "diff", "-U0", base, "--", f], capture_output=True, text=True).stdout
    def classify(text):
        out, inchunk = [], False
        for ln in text.split("\n"):
            if re.match(r"^```\{r", ln): out.append("chunk header"); inchunk = True
            elif re.match(r"^```\s*$", ln) and inchunk: out.append("chunk end"); inchunk = False
            else: out.append("chunk code" if inchunk else "prose")
        return out
    old = classify(subprocess.run(["git", "show", base + ":" + f], capture_output=True, text=True).stdout)
    new = classify(open(f).read())
    counts = {}
    o = n = 0
    for ln in d.split("\n"):
        m = re.match(r"^@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@", ln)
        if m: o, n = int(m.group(1)), int(m.group(2)); continue
        if ln.startswith("---") or ln.startswith("+++"): continue
        if ln.startswith("-"):
            k = ("removed", old[o - 1]); o += 1
        elif ln.startswith("+"):
            k = ("added", new[n - 1]); n += 1
        else: continue
        counts[k] = counts.get(k, 0) + 1
        if k[1] != "prose": print("  ", f, k, ln[:100])
    print(f, counts)
