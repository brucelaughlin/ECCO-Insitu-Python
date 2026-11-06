"""
The report site's questions page (``questions.html``): common questions as links into a preset view of the report,
the cost map explorer or the in situ data viewer, each with a choice of obs class (and cost type, variable).

Written by ``profile-report --site`` next to ``index.html``. Each question names a page and its settings; the
pages keep their whole view in the address, so a question is just a link. The same list can later be the
actions of a question-answering assistant (docs/planning/llm_assistant.md).
"""
import html as _html
import json

# (group, question, what you'll see, link template). Templates: {src} {ct} {V} are filled from the choices on
# the page; 'report:' = index.html, 'explorer:' = viewer.html, 'data:' = the in situ data viewer.
QUESTIONS = [
    ('The data themselves', 'What do the observations look like with depth (mean T and S)?',
     'analysis 01: global mean T and S against depth, and how many values there are at each depth',
     'report:?src={src}&ct={ct}&an=01'),
    ('The data themselves', 'Where are the observations, and how far are they from the climatology or the model?',
     'the in situ data viewer: anomaly per geodesic bin at the surface; zoom in for single profiles',
     'data:?src={src}&ref={ref}&V={V}&k=0'),
    ('How good is the fit, and at what depth?', 'At what depths do the observations disagree most with the reference?',
     'analysis 02: cost per datum against depth, unsplit and the two terms of the split',
     'report:?src={src}&ct={ct}&an=02'),
    ('How good is the fit, and at what depth?', 'Is the misfit a steady offset or variability in time?',
     'analysis 02 with the split: the mean term is the steady offset, the anomaly term the variability',
     'report:?src={src}&ct={ct}&an=02'),
    ('How good is the fit, and at what depth?', 'How good is the fit over the upper ocean, to 100, 250, 500, 700 or 2000 m?',
     'analysis 10: cost per datum over each depth range from the surface (0–100 m, 0–250 m, … 0–2000 m)',
     'report:?src={src}&ct={ct}&an=10'),
    ('Where?', 'Where is the cost high, all depths together?',
     'the cost map explorer: unsplit cost per datum per geodesic bin, all depths, all years',
     'explorer:?src={src}&ct={ct}&V={V}&part=unsplit&depth=all&years=all'),
    ('Where?', 'Where is the cost high at one depth?',
     'the cost map explorer at 1000 m (change the depth in the Depth row)',
     'explorer:?src={src}&ct={ct}&V={V}&part=unsplit&depth=k{k1000}&years=all'),
    ('Where?', 'Where are the steady offsets (the mean term)?',
     'the cost map explorer: the split mean term per bin at 1000 m',
     'explorer:?src={src}&ct={ct}&V={V}&part=mean&depth=k{k1000}'),
    ('Where?', 'Is the reference warmer or colder (saltier or fresher) than the observations, and where?',
     'the cost map explorer: mean reference − obs per bin (red: the reference is warmer/saltier, blue: colder/fresher), 1000 m',
     'explorer:?src={src}&ct={ct}&V={V}&part=bias&depth=k{k1000}&years=all'),
    ('Where?', 'How does the fit vary with latitude and depth?',
     'analysis 07: zonal mean cost per datum against latitude and depth',
     'report:?src={src}&ct={ct}&an=07'),
    ('When?', 'How does the fit change through time?',
     'analysis 03: cost per datum per month; analysis 04 is the total per month',
     'report:?src={src}&ct={ct}&an=03'),
    ('When?', 'Is the newer data fit worse than the data inside the climatology\'s years?',
     'supplemental analysis 11: cost against depth, inside vs after the climatology\'s years',
     'report:?src={src}&ct={ct}&an=11'),
    ('When?', 'How did one region change over the years?',
     'the cost map explorer by year: use the time slider, or select a region for its per-year bars',
     'explorer:?src={src}&ct={ct}&V={V}&part=bias&depth=k{k1000}&years=2015-2015'),
    ('Outliers and peaks', 'Are there outliers?',
     'analysis 08: histogram of the cost per datum (a long right tail or a second peak means outliers)',
     'report:?src={src}&ct={ct}&an=08'),
    ('Outliers and peaks', 'Where does a peak in the cost against depth come from?',
     'the cost map explorer opens the peak tool on 2000–2600 m (change the range there): sampling vs signal, '
     'the bins and areas that carry it, and their platforms',
     'explorer:?src={src}&ct={ct}&V={V}&part=mean&peak=2000-2600'),
    ('Outliers and peaks', 'Which instruments (floats, ships) fit worst in an area?',
     'the cost map explorer: select a region on the map; its panel lists the platforms and their shares',
     'explorer:?src={src}&ct={ct}&V={V}&part=unsplit&depth=k{k1000}&years=all'),
    ('Outliers and peaks', 'What does one profile look like against the climatology and the model?',
     'the in situ data viewer: zoom in until the red dots show, click one',
     'data:?src={src}&ref={ref}&V={V}&k={k1000}'),
    ('How good is the fit, and at what depth?', 'How far are the observations from the reference, in units of their assumed error?',
     'analysis 09: cost per datum (ΔV/σ)² against depth, all factors removed: 1 means one assumed error',
     'report:?src={src}&ct={ct}&an=09'),
    ('Comparing', 'How do the obs classes compare?',
     'analysis 02 for the chosen class; switch class with the Source tabs (cost type and analysis stay)',
     'report:?src={src}&ct={ct}&an=02'),
    ('Comparing', 'Is the fit better against the climatology or the model?',
     'analysis 02; switch between vs_clim and vs_model with the Cost type tabs (axes are shared)',
     'report:?src={src}&ct={ct}&an=02'),
    ('Comparing', 'Does the split depend on pooling the bin means over all classes (as ECCO does)?',
     'supplemental analysis 12: the anomaly term with own-class vs pooled bin means',
     'report:?src={src}&ct={ct}&an=12'),
]

CSS = """
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #222; font-size: 15px; }
header { padding: 0.8em 1.5em; border-bottom: 2px solid #d0d0d0; display: flex; gap: 1.2em; align-items: center; flex-wrap: wrap; }
header h1 { font-size: 1.15em; margin: 0; }
header a { color: #1a5ca8; font-size: 0.9em; }
#choose { background: #f6f7f9; border-bottom: 1px solid #e3e5e8; padding: 0.6em 1.5em; font-size: 0.92em; display: flex; gap: 1.4em;
  flex-wrap: wrap; align-items: center; position: sticky; top: 0; }
#choose label { display: flex; gap: 0.4em; align-items: center; }
select { font-size: 0.95em; padding: 0.15em; }
main { max-width: 62em; padding: 0.5em 1.5em 3em; }
h2 { font-size: 1.02em; color: #345; margin: 1.4em 0 0.4em; text-transform: uppercase; letter-spacing: 0.04em; }
a.q { display: block; border: 1px solid #dde3ea; border-radius: 6px; padding: 0.55em 0.8em; margin: 0.35em 0; text-decoration: none;
  color: #1a3d6b; background: #fff; }
a.q:hover { border-color: #1a5ca8; background: #f2f7fd; }
a.q b { font-weight: 600; }
a.q span { display: block; color: #667; font-size: 0.86em; margin-top: 0.15em; }
a.q i { float: right; font-style: normal; font-size: 0.78em; color: #889; }
p.lede { color: #555; }
"""


def page(title, sources, cost_types, depth, bin_viewer=None):
    """The questions page. ``sources`` / ``cost_types``: what the report has; ``depth``: the 97 levels (m)."""
    k1000 = min(range(len(depth)), key=lambda k: abs(depth[k] - 1000)) if depth else 0
    groups, items = [], {}
    for g, q, what, link in QUESTIONS:
        if g not in items:
            groups.append(g); items[g] = []
        items[g].append((q, what, link))
    body = ''.join(f'<h2>{_html.escape(g)}</h2>' + ''.join(
        f'<a class="q" data-t="{_html.escape(link)}" href="#"><i>{link.split(":")[0].replace("explorer", "cost map explorer").replace("data", "data viewer")}</i>'
        f'<b>{_html.escape(q)}</b><span>{_html.escape(what)}</span></a>' for q, what, link in items[g]) for g in groups)
    opt = lambda vals, sel: ''.join(f'<option{" selected" if v == sel else ""}>{_html.escape(v)}</option>' for v in vals)
    src0 = 'PFL' if 'PFL' in sources else sources[0]
    ct0 = 'vs_model' if 'vs_model' in cost_types else cost_types[0]
    js = f"""
const BINV = {json.dumps((bin_viewer or '').rstrip('/') + '/' if bin_viewer else '')}, K1000 = {k1000};
function fill() {{
  const src = document.getElementById('src').value, ct = document.getElementById('ct').value, V = document.getElementById('V').value;
  const ref = ct === 'vs_model' ? 'model' : 'clim';
  document.querySelectorAll('a.q').forEach(a => {{
    let t = a.dataset.t.replace(/{{src}}/g, src).replace(/{{ct}}/g, ct).replace(/{{V}}/g, V).replace(/{{ref}}/g, ref).replace(/{{k1000}}/g, K1000);
    const [kind, rest] = [t.slice(0, t.indexOf(':')), t.slice(t.indexOf(':') + 1)];
    if (kind === 'report') a.href = 'index.html' + rest;
    else if (kind === 'explorer') a.href = 'viewer.html' + rest;
    else if (kind === 'data') {{ if (BINV) a.href = BINV + rest + '&back=' + encodeURIComponent(location.href.split('?')[0]); else a.style.display = 'none'; }}
  }});
}}
['src', 'ct', 'V'].forEach(id => document.getElementById(id).addEventListener('change', fill)); fill();
"""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Questions — {_html.escape(title)}</title><style>{CSS}</style></head>
<body>
<header><h1>What do you want to know?</h1><a href="index.html">&larr; Profile report</a><a href="viewer.html">Cost map explorer</a>
  <a href="guide.html">Guide</a></header>
<div id="choose">For: <label>obs class <select id="src">{opt(sources, src0)}</select></label>
  <label>compared with <select id="ct">{opt(cost_types, ct0)}</select></label>
  <label>variable <select id="V">{opt(['T', 'S'], 'T')}</select></label>
  <span style="color:#777">(each question opens with these; you can change everything there)</span></div>
<main><p class="lede">Each question opens a page already set up to answer it. The analyses are those of the report
({_html.escape(title)}).</p>{body}</main>
<script>{js}</script>
</body></html>
"""
