"""
The report site's user guide (``guide.html``): how to use the report, the cost map explorer and the in situ data
viewer, and how the attribution tools work, step by step, with the formulas.

Written by ``profile-report --site`` next to ``index.html``. The report's Notes & errata page keeps the run-specific
findings, decisions and corrections; this page is the manual and does not depend on the run.
"""
import html as _html

# (id, title, body). Order = the table of contents.
SECTIONS = [
    ('start', 'Start here', """
<p>The site has three parts, linked from the top right of every page:</p>
<table class="g">
<tr><th>Page</th><th>What it shows</th><th>Use it to</th></tr>
<tr><td><b>Profile report</b> (the main page)</td><td>Fixed figures of the analyses (01&ndash;10 core, 11&ndash;12 supplemental) for each
  obs class and cost type</td><td>see the overall picture: where and at what depths the observations disagree with the
  climatology or the model</td></tr>
<tr><td><b>Cost map explorer</b></td><td>The same costs drawn live: maps for any depth, years and cost part; the depth profile;
  regions; the peak tool</td><td>find <em>where</em> a feature in the report comes from</td></tr>
<tr><td><b>In situ data viewer</b></td><td>The observations themselves: anomalies per geodesic bin, every profile when zoomed in, the
  profile panel</td><td>look at the profiles and platforms behind a feature</td></tr>
</table>
<p>A typical path: a bump in an analysis&nbsp;02 curve &rarr; &ldquo;investigate&rdquo; (below the figure) &rarr; the explorer finds the peak,
separates sampling from signal, outlines the bins that carry it and lists their platforms &rarr; &ldquo;show profiles&rdquo; opens the
data viewer on one platform &rarr; click a profile to see it whole. Section&nbsp;<a href="#attribution">Attribution</a> explains each step.</p>
<p class="note">This page is the manual. Run-specific findings, decisions and corrections are on the report&rsquo;s
<b>Notes &amp; errata</b> page.</p>
"""),
    ('terms', 'Terms', """
<table class="g">
<tr><th>Term</th><th>Meaning</th></tr>
<tr><td>obs class</td><td>A kind of instrument: CTD/bottle casts (CTD_WOD), gliders (GLD_WOD), ice-tethered profilers (ITP), seal tags
  (MEOP), moorings (MRB_WOD), Argo floats (PFL), XBTs (XBT_WOD). Each is analyzed on its own.</td></tr>
<tr><td>reference, V<sub>ref</sub></td><td>What an observation V is compared with: the WOA23 climatology (<code>vs_clim</code>) or the ECCO
  v4r6 model (<code>vs_model</code>; monthly means sampled at the profile&rsquo;s grid cell and month).</td></tr>
<tr><td>misfit &Delta;V</td><td>V &minus; V<sub>ref</sub>, in &deg;C or psu.</td></tr>
<tr><td>&sigma;</td><td>The assumed error of an observation (its uncertainty), from the processing chain.</td></tr>
<tr><td>cost per datum</td><td>(&Delta;V/&sigma;)&sup2; for one value (&ldquo;true&rdquo; weights). 1 means one assumed error away; 4 means two.</td></tr>
<tr><td>geodesic bin</td><td>One of 10,242 cells of roughly equal area (about 220&nbsp;km across) covering the globe.</td></tr>
<tr><td>depth level</td><td>One of the 97 standard depths (2&nbsp;m to 6000&nbsp;m): 2&ndash;5&nbsp;m apart in the top 100&nbsp;m (23 levels), wider below
  (50&nbsp;m around 1000&nbsp;m, 100&nbsp;m from 1100 to 2000&nbsp;m), 200&nbsp;m apart below 2000&nbsp;m.</td></tr>
<tr><td>unsplit cost</td><td>w&middot;(&Delta;V)&sup2; per value: the whole misfit.</td></tr>
<tr><td>split: anomaly and mean term</td><td>ECCO&rsquo;s split of the cost per geodesic bin <i>i</i> and level <i>k</i>:
  <b>anomaly term</b> per value w&middot;(&Delta;V &minus; &Delta;V&#772;<sub>i,k</sub>)&sup2;, the misfit after removing the bin&rsquo;s
  long-term mean misfit &Delta;V&#772;<sub>i,k</sub>; and <b>mean term</b> per bin and level
  J<sub>i,k</sub> = w&#772;&middot;(&Delta;V&#772;<sub>i,k</sub>)&sup2;, the steady offset itself, with
  w&#772; = 1/max(s/n, s/9), s = mean(1/w) over the bin&rsquo;s n values. (The &ldquo;9&rdquo; caps how much n independent values can
  sharpen the bin mean.)</td></tr>
<tr><td>weights final / nosf / true</td><td>w = &gamma;&middot;sf/&sigma;&sup2; as ECCO uses it (final); without the spatial scaling sf (nosf);
  without &gamma; and sf, i.e. 1/&sigma;&sup2; (true). The explorer and the attribution use <b>true</b> weights.</td></tr>
<tr><td>kept / rejected</td><td>The processing chain&rsquo;s quality control (step 06) sets a value&rsquo;s weight to 0 when it fails a test; such a
  value is <em>rejected</em> and enters no cost. Everything here uses kept values only.</td></tr>
<tr><td>platform</td><td>The instrument a profile comes from: an Argo float (WMO number), an ITP, a seal tag, a mooring; for the WOD ship
  classes, a ship or program code (coarser).</td></tr>
</table>
"""),
    ('report', 'The profile report', """
<p>Choose the obs class, cost type and analysis with the tabs. Click a figure to see it full width. The side panel describes
the figure; <b>Terminology</b> and <b>Metadata</b> open in it. Within a class, a figure uses the same axes in every cost type,
so the tabs can be compared by eye.</p>
<p>Analysis&nbsp;02 (cost per datum against depth) is where most investigations start. Below it, the
<b>&ldquo;Investigate it&rdquo;</b> box opens the cost map explorer on the same class and cost type with a chosen cost part, variable
and depth range, and runs the peak tool there (section <a href="#peak">The peak tool</a>).</p>
"""),
    ('explorer', 'The cost map explorer', """
<p><b>Rows of buttons</b> (top): class, cost type, variable, cost part (unsplit, anomaly, mean term, number of points), the mean-bias
diagnostic (obs &minus; reference, not a cost), depth (all, one level, 100&nbsp;m bands, groups of levels, flagged spikes) and time
(all years, or any range; the mean term has no years).</p>
<p><b>Maps:</b> the global map (equal area) and the two polar maps; each bin colored by the selection. Hover a bin for its value;
<b>click a bin</b> for its breakdown (the &ldquo;Clicked bin&rdquo; panel: per level, per year, and a link to its profiles).</p>
<p><b>Depth profile of the selection:</b> the curve analysis&nbsp;02 draws, for the current class, cost type, variable, part and years:
for costs, at each level the mean over all its values (cost per datum); for the mean term, the mean over bins; for &ldquo;# points&rdquo;,
the number of values. On the right, the number of bins per level (log scale).</p>
<p><b>Depth range</b> (the shaded band on the profile): set it with the from/to menus or on the plot &mdash; drag the band&rsquo;s top or
bottom edge (the handles), drag inside to move it, drag across an empty part to draw a new one. Choosing a range does nothing by
itself. It is used by:</p>
<ul><li><b>Investigate the peak in this range</b> (the peak tool), and</li>
<li><b>also draw the curve over the same bins at every depth of the range</b>: a second, blue curve over only the bins that have data
at every level of the range (section <a href="#samebins">Same bins</a>).</li></ul>
<p><b>Region</b>: &ldquo;select region&rdquo;, then drag a box on any map (shift-drag removes bins, shift-click toggles one bin). The region is
outlined in purple and the region panel opens (section <a href="#region">Region vs the rest</a>).</p>
"""),
    ('viewer', 'The in situ data viewer', """
<p>The observations themselves: per geodesic bin, the mean anomaly (face color) and its standard deviation (edge color) at any depth
or depth range and any years, against the climatology or the model. Zoomed in past the threshold, each bin splits into one cell per
profile, with the true locations as red dots. Polar maps are below the world map.</p>
<p><b>Platform box:</b> type or pick a platform to ring all its profiles in cyan, at any zoom (its track).</p>
<p><b>Click a profile</b> (a red dot or a ring, zoomed in) for the <b>profile panel</b>:</p>
<ul>
<li><b>this profile</b>: observation (black), climatology (blue), model (red; monthly mean), and up to 60 neighbors in gray (same bin,
  &plusmn;2 years, other platforms); on the right, obs &minus; reference with &plusmn;1&sigma; and &plusmn;2&sigma; bands; &times; marks values the
  chain rejected; the values at the map&rsquo;s depth are listed, with the rejection reason if any.</li>
<li><b>this platform over time</b>: all the platform&rsquo;s profiles as a depth &times; time picture of obs &minus; reference.</li>
<li><b>largest misfits here</b>: the bin&rsquo;s profiles ranked by |&Delta;V|/&sigma; at the map&rsquo;s depth; click a row to open one.</li>
</ul>
"""),
    ('attribution', 'Attribution: from a cost peak to the profiles behind it', """
<p>&ldquo;Attribution&rdquo; here means: for a feature in a cost curve (a peak at some depth), find how much of it is due to <em>which data are
counted</em> at that depth, <em>which bins</em> carry it, and <em>which platforms and profiles</em> are in them. Every step is a fixed
calculation on the saved results; nothing is fitted or judged. The page reports the numbers; reading them is up to you.</p>
<p>The steps, in the order the tools apply them:</p>
<ol>
<li><a href="#peak">Find the peak</a> in the chosen depth range.</li>
<li><a href="#samebins">Same bins</a>: is the peak in the water, or in the sampling?</li>
<li><a href="#contrib">Bin contributions</a>: how much each bin adds to the curve at the peak.</li>
<li><a href="#clusters">Clusters</a>: group the largest contributors into areas, and name them.</li>
<li><a href="#region">Region vs the rest</a>: does the feature disappear without the area?</li>
<li><a href="#platforms">Platforms</a>: which instruments, profile by profile.</li>
<li><a href="#profiles">Profiles</a>: look at them.</li>
</ol>
<p>A worked example (PFL, vs_model, T, mean term, 2000&ndash;2600&nbsp;m) runs through all of them at the <a href="#example">end</a>.</p>
"""),
    ('peak', '1. The peak', """
<p>Let c<sub>k</sub> be the curve at level k (what the depth profile draws). For costs it is the cost per datum over all values at the
level: c<sub>k</sub> = &Sigma;<sub>i</sub> C<sub>i,k</sub> / &Sigma;<sub>i</sub> N<sub>i,k</sub>, with C<sub>i,k</sub> the summed cost and
N<sub>i,k</sub> the number of values in bin i. For the mean term it is the mean over the B<sub>k</sub> bins with data:
c<sub>k</sub> = &Sigma;<sub>i</sub> J<sub>i,k</sub> / B<sub>k</sub>. These are exactly analysis&nbsp;02&rsquo;s curves.</p>
<p>The <b>peak</b> is the level p with the largest c<sub>k</sub> in the depth range. Its <b>neighbors</b> are the nearest levels with data
above (p&minus;) and below (p+); the <b>excess</b> is E = c<sub>p</sub> &minus; (c<sub>p&minus;</sub> + c<sub>p+</sub>)/2.</p>
"""),
    ('samebins', '2. Same bins: sampling or signal?', """
<p>A curve like c<sub>k</sub> can jump because the water changes, or because <em>different bins</em> have data at different depths
(e.g. most Argo floats stop at 2000&nbsp;m; below it only a few regions are sampled). To tell them apart, the same calculation is
repeated over only the bins S that have data at p and at both neighbors:</p>
<p class="f">c<sub>k</sub><sup>S</sup> = the curve at level k over the bins in S only, &nbsp;
E<sup>S</sup> = c<sub>p</sub><sup>S</sup> &minus; (c<sub>p&minus;</sub><sup>S</sup> + c<sub>p+</sub><sup>S</sup>)/2,</p>
<p class="f">share of the excess due to sampling = 1 &minus; E<sup>S</sup>/E (clipped to 0&ndash;100&nbsp;%).</p>
<p>If the peak survives in a fixed set of bins (E<sup>S</sup> &asymp; E), it is in the water. If it disappears (E<sup>S</sup> &asymp; 0), it comes
from the change in which bins are counted: the bins that have data only at the peak&rsquo;s depth are unusual. That does <em>not</em> make
the peak wrong: those bins really are badly fit; but it is a regional feature, not a feature of the global water column.</p>
<p>The checkbox &ldquo;also draw the curve over the same bins at every depth of the range&rdquo; draws c<sup>S</sup> over a whole range
(S = bins with data at every level of the range), so you can see it by eye.</p>
"""),
    ('contrib', '3. Bin contributions', """
<p>Each bin&rsquo;s part of the curve at the peak, chosen so the parts add up:</p>
<ul>
<li>costs: bin i contributes its summed cost C<sub>i,p</sub>; the curve is &Sigma;C/&Sigma;N, so the bin&rsquo;s <b>share</b> of the level is
  C<sub>i,p</sub>/&Sigma;<sub>j</sub>C<sub>j,p</sub>.</li>
<li>mean term: bin i contributes J<sub>i,p</sub>; share J<sub>i,p</sub>/&Sigma;<sub>j</sub>J<sub>j,p</sub>.</li>
</ul>
<p>A share always comes with the bin&rsquo;s share of the <em>bins</em> (or values): a region holding 4&nbsp;% of the bins and 34&nbsp;% of the term
carries 8&times; its share. A large share from many bins is ordinary; a large share from few bins is the signal.</p>
"""),
    ('clusters', '4. Clusters and their names', """
<p>The bins at the peak whose contribution is in the top 10&nbsp;% are grouped into <b>clusters</b>: two high bins belong to the same
cluster if their centers are within 600&nbsp;km (chained, so a cluster can be long). Each cluster gets its center (the mean
direction of its bins&rsquo; centers), its number of bins, its share of the bins and its share of the summed value at the peak. The
largest cluster is outlined on the maps and becomes the region; click another row to switch.</p>
<p><b>Names</b> are descriptive only: the smallest Natural Earth sea, gulf or bay containing most of the cluster&rsquo;s bin centers;
otherwise a coarse hand-drawn box (e.g. &ldquo;subpolar North Atlantic&rdquo;, &ldquo;Brazil Basin&rdquo;, &ldquo;Antarctic Circumpolar Current
(south of Australia)&rdquo;); otherwise the ocean. The coordinates are always shown next to the name.</p>
<p>A cluster holds only the high bins. To include their surroundings (to compare with nearby, well-fit bins), widen the region by
dragging on the maps.</p>
"""),
    ('region', '5. Region vs the rest', """
<p>For any set of bins R (a cluster, or one you drag):</p>
<ul>
<li>the depth profile draws the curve over R (purple), over the rest of the world (gray) and over all bins (light red), each formed
  as in step&nbsp;1. If the region&rsquo;s curve runs past the plot, it gets its own axis along the top. Levels where R has fewer than 5 bins
  are drawn pale. <b>The decisive check:</b> if the feature is gone from the gray curve, R accounts for it.</li>
<li>a table: R&rsquo;s bins with data, values (points), the curve over R and over the rest (per datum, or per bin for the mean term),
  and R&rsquo;s share of the summed value, with how many times its share of bins that is;</li>
<li>per year: R (bars) against the rest (line); for the mean term, which has no years, obs &minus; reference per year;</li>
<li>R&rsquo;s ten largest bins (click one for its breakdown).</li>
</ul>
<p><b>Caution:</b> a region drawn around the darkest bins will always carry a large share. Judge it by its share relative to its share
of bins and values, and by the rest-of-the-world curve, not by the share alone.</p>
"""),
    ('platforms', '6. Platforms, profile by profile', """
<p>The region panel lists the platforms in R with their profiles, years, mean obs &minus; reference and <b>share of the region&rsquo;s
value</b>. The shares are computed profile by profile from every profile&rsquo;s own values (the data viewer&rsquo;s profile details), kept
values at the selected levels only, true weights (w = 1/&sigma;&sup2;):</p>
<ul>
<li><b>unsplit cost</b>: each value&rsquo;s own (&Delta;V/&sigma;)&sup2;, summed per platform;</li>
<li><b>anomaly term</b>: each value&rsquo;s (&Delta;V &minus; &Delta;V&#772;<sub>i,k</sub>)&sup2;/&sigma;&sup2;, with &Delta;V&#772;<sub>i,k</sub> the bin&rsquo;s mean
  misfit over all years;</li>
<li><b>mean term</b>: J<sub>i,k</sub> = w&#772;&middot;(&Delta;V&#772;)&sup2; is shared among the bin&rsquo;s n values as
  w&#772;&middot;&Delta;V&#772;&middot;&Delta;V<sub>v</sub>/n for value v. The parts add up to J exactly (since &Sigma;<sub>v</sub>&Delta;V<sub>v</sub> =
  n&middot;&Delta;V&#772;); a value pulling against its bin&rsquo;s mean gets a negative part.</li>
</ul>
<p>Shares are of the region&rsquo;s total over its platforms. With a time window, values outside it are left out of the shares, but the
bin means of the split use all years, as the analysis does.</p>
<p><b>Checked</b> (2026-11-04): the mean term rebuilt this way from the profile details agrees with the analysis for every bin
(PFL vs model, T, 2400&nbsp;m: 1,046 bins, identical counts) to 0.14&nbsp;% of its total; the difference is the details&rsquo; rounding
to 0.001. If the details can&rsquo;t be loaded, the table falls back to splitting each (bin, year)&rsquo;s value among its platforms by
their numbers of profiles there, and says so.</p>
"""),
    ('profiles', '7. Profiles', """
<p>&ldquo;show profiles&rdquo; opens the in situ data viewer with the platform&rsquo;s profiles ringed, at the same class, reference and depth.
Zoom in and click one for the profile panel (section <a href="#viewer">The in situ data viewer</a>). What to look for:</p>
<ul>
<li><b>against its neighbors</b> (gray): if the platform agrees with other platforms nearby and the reference doesn&rsquo;t, the
  difference is the reference&rsquo;s (or the water&rsquo;s); if it alone is off, suspect the platform;</li>
<li><b>over time</b> (second tab): a difference that persists from profile to profile is a water-mass (or model) difference; one
  that creeps is a likely sensor drift; a single bad profile stands alone;</li>
<li><b>against &sigma;</b>: how many assumed errors the misfit is; above about 3 is unusual;</li>
<li><b>rejected values</b> (&times;): the chain already excluded them.</li>
</ul>
"""),
    ('example', 'Worked example', """
<p>PFL (Argo), vs_model, T, mean term, depth range 2000&ndash;2600&nbsp;m (run 20261027_142247):</p>
<ol>
<li><b>Peak:</b> 2400&nbsp;m, 5.57 over 1,046 bins; neighbors 2200&nbsp;m 4.90 (1,181 bins) and 2600&nbsp;m 2.62 (410 bins).</li>
<li><b>Same bins:</b> over the 410 bins with data at all three levels, 2.83 vs 2.79: 98&nbsp;% of the excess is from which bins are
  counted. Most floats stop at 2000&nbsp;m; the few deeper ones are in a few regions.</li>
<li><b>Clusters:</b> the largest is the subpolar North Atlantic (56&deg;N 38&deg;W): 4.1&nbsp;% of the bins, 33.8&nbsp;% of the term.</li>
<li><b>Region vs the rest:</b> the region&rsquo;s mean term per bin is 45.9 against 3.84 for the rest of the world.</li>
<li><b>Platforms:</b> four floats deployed there since 2021 (6903082, 6903871, 6903085, 1902621) carry about 40&nbsp;% of the region&rsquo;s term,
  each with obs &minus; model about &minus;0.37 to &minus;0.40&nbsp;&deg;C.</li>
<li><b>Profiles:</b> float 6903082 at 2400&nbsp;m in March 2025: obs 2.810&nbsp;&deg;C, model 3.269&nbsp;&deg;C (&minus;0.46&nbsp;&deg;C, 3.4&sigma;); it agrees with
  the other floats nearby; the cold layer at 2000&ndash;2500&nbsp;m persists over its life.</li>
</ol>
<p><i>Reading (not a test the tools run):</i> the model is warmer than the floats in the subpolar North Atlantic at 2000&ndash;2500&nbsp;m; the
global mean term peaks at 2400&nbsp;m because few other regions are sampled there. The model values are monthly means, not the model&rsquo;s
values at each profile&rsquo;s time.</p>
"""),
    ('cautions', 'Cautions', """
<ul>
<li><b>Few bins or values</b> make any number noisy: look at the counts shown next to everything (pale = under 5 bins).</li>
<li><b>Regions drawn after looking</b> always look important; compare shares with shares of bins and values.</li>
<li><b>vs_model is approximate:</b> monthly-mean model fields, so sub-monthly variability counts as misfit; no model values after
  December 2025.</li>
<li><b>Precision:</b> the explorer stores 3 significant digits per value; the profile details 0.001 (&sigma; 0.0001).</li>
<li><b>Platform IDs</b> for the WOD ship classes are ship or program codes, not single instruments.</li>
</ul>
"""),
]

CSS = """
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #222; font-size: 15px; display: flex; }
#toc { position: sticky; top: 0; height: 100vh; overflow-y: auto; width: 17em; flex-shrink: 0; background: #f6f7f9;
  border-right: 1px solid #ddd; padding: 1em 0.9em; font-size: 0.9em; }
#toc h1 { font-size: 1.05em; margin: 0 0 0.6em; }
#toc a { display: block; color: #1a5ca8; text-decoration: none; padding: 0.18em 0.3em; border-radius: 3px; }
#toc a:hover { background: #e8f0fb; }
#toc a.sub { padding-left: 1.3em; font-size: 0.95em; color: #345; }
#toc a.on { background: #1a5ca8; color: #fff; }
#toc .links { margin-top: 1.2em; border-top: 1px solid #ddd; padding-top: 0.6em; }
#toc input { width: 100%; padding: 0.3em; margin-bottom: 0.6em; border: 1px solid #ccc; border-radius: 4px; }
main { max-width: 52em; padding: 1.2em 2em 4em; line-height: 1.6; }
main h2 { font-size: 1.3em; margin: 1.8em 0 0.4em; padding-top: 0.4em; border-top: 1px solid #eee; }
main h2:first-child { margin-top: 0.2em; border-top: none; }
table.g { border-collapse: collapse; margin: 0.6em 0; font-size: 0.92em; }
table.g td, table.g th { border: 1px solid #ddd; padding: 4px 8px; text-align: left; vertical-align: top; }
table.g th { background: #f2f2f2; }
p.f { background: #f8f6fb; border-left: 3px solid #6c3483; padding: 0.4em 0.8em; font-family: Georgia, serif; }
p.note { color: #666; font-style: italic; }
code { background: #f0f0f0; padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }
.hide { display: none; }
"""

JS = """
const links = [...document.querySelectorAll('#toc a[href^="#"]')], secs = links.map(a => document.getElementById(a.getAttribute('href').slice(1)));
const onScroll = () => { let cur = secs[0]; for (const s of secs) if (s && s.getBoundingClientRect().top < 120) cur = s;
  links.forEach(a => a.classList.toggle('on', a.getAttribute('href') === '#' + cur.id)); };
document.addEventListener('scroll', onScroll); onScroll();
document.getElementById('q').addEventListener('input', e => { const q = e.target.value.trim().toLowerCase();
  secs.forEach((s, i) => { const hit = !q || s.textContent.toLowerCase().includes(q); s.classList.toggle('hide', !hit);
    links[i].classList.toggle('hide', !hit); }); });
"""

SUBS = {'peak', 'samebins', 'contrib', 'clusters', 'region', 'platforms', 'profiles', 'example'}


def page(title, bin_viewer=None):
    """The guide as one HTML page."""
    toc = ''.join(f'<a href="#{i}" class="{"sub" if i in SUBS else ""}">{t}</a>' for i, t, _ in SECTIONS)
    body = ''.join(f'<section id="{i}"><h2>{t}</h2>{b}</section>' for i, t, b in SECTIONS)
    extra = (f'<a href="{_html.escape(bin_viewer)}">In situ data viewer &rarr;</a>' if bin_viewer else '')
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Guide — {_html.escape(title)}</title><style>{CSS}</style></head>
<body>
<nav id="toc"><h1>Guide</h1><input id="q" placeholder="search this guide" aria-label="search this guide">{toc}
  <div class="links"><a href="index.html">&larr; Profile report</a><a href="viewer.html">Cost map explorer &rarr;</a>{extra}</div></nav>
<main>{body}</main>
<script>{JS}</script>
</body></html>
"""
