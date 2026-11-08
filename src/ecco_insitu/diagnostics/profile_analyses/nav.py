"""
The site's navigation bar: the same links, in the same place, on every page (the analyses, the cost map explorer,
the in situ data viewer, the questions page and the guide), with the current page marked. Replaces the old
"back to the report" links, so any page is one click from any other.

The in situ data viewer is a separate site; it builds the same bar in its own page (bin_viewer/viewer.html,
``siteNav``) from the address of the analyses site it was opened from.
"""
import html as _html

# (key, file, label, tooltip); the data viewer's address is filled in per site
PAGES = [
    ('analyses', 'index.html', 'Analyses', 'The figures of the analyses, per obs class and cost type'),
    ('explorer', 'viewer.html', 'Cost map explorer', 'Cost maps for any depth, years and region; the peak tool'),
    ('data', None, 'In situ data viewer', 'The observations themselves: every profile, per geodesic bin'),
    ('questions', 'questions.html', 'What do you want to know?', 'Common questions, each opening a page set up to answer it'),
    ('guide', 'guide.html', 'Guide', 'How to use every page; what the terms mean; how the attribution works'),
]

CSS = """
.navrow { display: flex; align-items: center; gap: 1em; padding: 0.35em 1em; background: #f2f6fb; border-bottom: 1px solid #d6e0ec; }
.navrow .site { font: 600 13px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #345; white-space: nowrap; }
.navrow .run { font: 12px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #778; margin-left: auto; }
.sitenav { display: flex; gap: 0.25em; align-items: center; flex-wrap: wrap; font-family: -apple-system, BlinkMacSystemFont,
  'Segoe UI', sans-serif; font-size: 13px; }
.sitenav a { text-decoration: none; color: #1a5ca8; padding: 0.22em 0.65em; border: 1px solid #c9d7ea; border-radius: 4px;
  background: #fff; white-space: nowrap; }
.sitenav a:hover { background: #eef4fc; border-color: #1a5ca8; }
.sitenav a.here { background: #1a5ca8; color: #fff; border-color: #1a5ca8; cursor: default; }
.sitenav a.navq { color: #6c3483; border-color: #d9c8e6; }
.sitenav a.navq.here { background: #6c3483; color: #fff; border-color: #6c3483; }
"""


def row(current, title, bin_viewer=None, ids=None):
    """The navigation row at the very top of a page: the site's name, the bar, and the run's title at the right."""
    return (f'<div class="navrow"><span class="site">ECCO in situ profiles</span>{bar(current, bin_viewer, ids)}'
            f'<span class="run" title="the run these pages are made from">{_html.escape(title)}</span></div>')


def bar(current, bin_viewer=None, ids=None):
    """The bar's HTML. ``current``: this page's key; ``bin_viewer``: the data viewer's address (its link is left
    out without one); ``ids``: {key: element id} for links a page updates itself (e.g. to carry the obs class)."""
    out = []
    for key, href, label, tip in PAGES:
        if key == 'data':
            if not bin_viewer:
                continue
            href = bin_viewer.rstrip('/') + '/'
        cls = ' '.join(c for c in ('here' if key == current else '', 'navq' if key == 'questions' else '') if c)
        idattr = f' id="{ids[key]}"' if ids and key in ids else ''
        out.append(f'<a href="{_html.escape(href)}"{idattr}{f" class={chr(34)}{cls}{chr(34)}" if cls else ""} '
                   f'title="{_html.escape(tip)}">{_html.escape(label)}</a>')
    return '<nav class="sitenav" aria-label="pages of this site">' + ''.join(out) + '</nav>'
