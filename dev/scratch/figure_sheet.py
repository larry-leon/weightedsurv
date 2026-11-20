#!/usr/bin/env python3
"""Before/after figure sheet for TASK_vignette_revision.md (Gate G12).
Pairs the figures of the knit at ac76f82 (dev/scratch/vignette_before/) with those of the knit at HEAD
(dev/scratch/vignette_after/), by order of appearance in the knitted file, and writes one self-contained HTML
file (images embedded).   python3 dev/scratch/figure_sheet.py <out.html>"""
import base64, hashlib, html, os, re, sys

FILES = [("examples", "vignettes/weightedsurv_examples.Rmd"),
         ("martingale", "vignettes/articles/martingale_resampling_applications.Rmd"),
         ("coxmeth", "vignettes/articles/weightedcox_methodology.Rmd"),
         ("kmmeth", "vignettes/articles/weightedkm_methodology.Rmd"),
         ("rmstmeth", "vignettes/articles/weightedrmst_methodology.Rmd"),
         ("sims", "vignettes/articles/weightedcox_wald-vs-logrank_simulations.Rmd")]
WHAT = {
    "examples": ["Flowchart of the weighting schemes (wt.rg.S)", "Weight functions (plot_weight_schemes)",
                 "GBSG K-M curves, plot_weighted_km and survfit", "GBSG K-M curves with 95% CIs",
                 "Rotterdam K-M curves, un-weighted and weighted", "GBSG and Rotterdam weighted K-M curves",
                 "Rotterdam weighted survival difference: resampled curves and band",
                 "Rotterdam weighted and GBSG: survival-difference bands and cumulative RMST",
                 "GBSG and cross-trial (CTC) K-M curves", "Cumulative weighted log-rank z (cross-trial)",
                 "GBSG and cross-trial: survival-difference bands and cumulative RMST",
                 "GBSG subgroups with ITT band", "GBSG er > 0: resampled curves and band",
                 "check.seKM: weighted K-M standard errors against survfit (draws = 0)",
                 "check.seKM: resampling standard errors (draws = 5000) against survfit",
                 "Rotterdam weighted K-M curves (for the adjustedCurves comparison)",
                 "adjustedCurves iptw_km curves", "Standard errors: weightedsurv (grey) and adjustedCurves (red)"],
    "martingale": ["GBSG K-M curves", "Rotterdam K-M curves, un-weighted and weighted", "GBSG and Rotterdam weighted K-M curves",
                   "Rotterdam weighted survival difference: resampled curves", "Rotterdam weighted survival difference: band",
                   "GBSG and Rotterdam weighted: survival-difference bands and cumulative RMST"],
    "coxmeth": ["Weight functions (plot_weight_schemes)"],
}

def figs(which, nm):
    d = os.path.join("dev/scratch", "vignette_" + which, nm)
    md = os.path.join(d, nm + ".md")
    if not os.path.exists(md):
        return None
    seen, out = set(), []
    for f in re.findall(r"fig-[A-Za-z0-9_.-]*?\.png", open(md, encoding="utf-8").read()):
        if f not in seen and os.path.exists(os.path.join(d, f)):
            seen.add(f); out.append(os.path.join(d, f))
    return out

def chunk(path):
    return re.sub(r"-(\d+)\.png$", r" (figure \1)", os.path.basename(path)[4:])

def img(path):
    b = open(path, "rb").read()
    return '<img src="data:image/png;base64,%s" alt="%s">' % (base64.b64encode(b).decode(), html.escape(path)), hashlib.sha256(b).hexdigest()

rows, summary = [], []
for nm, rmd in FILES:
    b, a = figs("before", nm), figs("after", nm)
    if nm == "sims":
        # the file is not changed by this task: one knit serves as before and after
        a = b
    if b is None and a is None:
        summary.append((rmd, "not knitted", 0, 0)); continue
    b, a = b or [], a or []
    if not b and not a:
        summary.append((rmd, "no figures (no evaluated chunk draws one)", 0, 0)); continue
    assert len(b) == len(a), (nm, len(b), len(a))
    same = 0
    rows.append('<h2>%s</h2>' % html.escape(rmd))
    if nm == "sims":
        rows.append('<p class="note">This file is byte-identical at ac76f82 and at HEAD; it was knitted once, and the same figure is shown on both sides.</p>')
    for k, (fb, fa) in enumerate(zip(b, a)):
        ib, hb = img(fb); ia, ha = img(fa)
        ident = hb == ha; same += ident
        what = WHAT.get(nm, [""] * len(b))[k] if k < len(WHAT.get(nm, [])) else ""
        rows.append('<section><h3>%d. %s</h3><p class="lab">%s &middot; chunk before: <code>%s</code> &middot; chunk after: <code>%s</code> &middot; <span class="%s">%s</span></p>'
                    '<div class="pair"><figure><figcaption>Before (ac76f82)</figcaption>%s</figure><figure><figcaption>After (HEAD)</figcaption>%s</figure></div></section>'
                    % (k + 1, html.escape(what), html.escape(rmd), html.escape(chunk(fb)), html.escape(chunk(fa)),
                       "same" if ident else "diff", "image files byte-identical" if ident else "image files differ", ib, ia))
    summary.append((rmd, "%d figure pair(s)" % len(b), len(b), same))

tab = "".join("<tr><td><code>%s</code></td><td>%s</td><td>%s</td></tr>" % (html.escape(r), html.escape(s), ("%d of %d" % (i, n)) if n else "&ndash;") for r, s, n, i in summary)
doc = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>weightedsurv vignette figures, before and after</title><style>
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:24px;color:#1a1a1a;background:#fff;line-height:1.45}
h1{font-size:22px}h2{font-size:18px;margin-top:40px;border-top:2px solid #444;padding-top:14px}h3{font-size:15px;margin:26px 0 2px}
.lab{font-size:12.5px;color:#444;margin:0 0 8px}.note{font-size:13px;color:#444}
.pair{display:grid;grid-template-columns:1fr 1fr;gap:14px}figure{margin:0;border:1px solid #ccc;padding:6px}
figcaption{font-size:12.5px;font-weight:600;margin-bottom:4px}img{width:100%%;height:auto;display:block}
.same{color:#555}.diff{color:#8a4b00;font-weight:600}table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #bbb;padding:4px 8px;text-align:left}
@media(max-width:800px){.pair{grid-template-columns:1fr}}
</style></head><body><h1>weightedsurv vignettes: every figure before (ac76f82) and after (HEAD) the vignette revision</h1>
<p class="note">Figures are paired by order of appearance in the knitted file and labelled by file and chunk. Both sets were knitted with <code>knitr::knit</code> (png, 96 dpi) against the same build of the package (<code>R/</code> is identical at the two commits).</p>
<table><tr><th>File</th><th>Figures</th><th>Pairs with byte-identical image files</th></tr>%s</table>%s</body></html>""" % (tab, "".join(rows))
open(sys.argv[1], "w", encoding="utf-8").write(doc)
for r, s, n, i in summary:
    print("%-70s %-45s identical %d of %d" % (r, s, i, n))
print("wrote", sys.argv[1], "%.1f MB" % (os.path.getsize(sys.argv[1]) / 1e6))
