// The ask-a-question box (prototype). Added to the report site's pages only by the local assistant server
// (ecco-insitu assistant), so the published site never shows it. It sends the question and the page's current
// view to the server, and shows the answer, the required AI notice and links to the views behind the numbers.
(() => {
  'use strict';
  const PAGES = {'index.html': 'report', '': 'report', 'viewer.html': 'explorer', 'questions.html': 'questions', 'guide.html': 'guide'};
  const here = PAGES[location.pathname.split('/').pop()] || 'other';
  let STATUS = null;
  const css = `
#asst-btn { position: fixed; right: 18px; bottom: 18px; z-index: 50; background: #6c3483; color: #fff; border: none;
  border-radius: 22px; padding: 0.6em 1.1em; font: 600 14px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
  box-shadow: 0 2px 8px rgba(0,0,0,0.25); cursor: pointer; }
#asst { position: fixed; right: 18px; bottom: 68px; z-index: 50; width: min(30em, calc(100vw - 36px)); max-height: 45vh;
  display: none; flex-direction: column; background: #fff; border: 1px solid #cdbfdc; border-radius: 8px;
  box-shadow: 0 4px 18px rgba(0,0,0,0.2); font: 14px/1.45 -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: #222; }
#asst.on { display: flex; }
#asst .hd { padding: 0.5em 0.8em; background: #f4effa; border-bottom: 1px solid #e3d9ee; border-radius: 8px 8px 0 0;
  display: flex; justify-content: space-between; align-items: center; }
#asst .hd b { color: #4a235a; }
#asst .hd span { font-size: 0.8em; color: #777; flex: 1; margin-left: 0.6em; }
#asst .hd button { background: none; border: 1px solid #cdbfdc; border-radius: 4px; color: #4a235a; cursor: pointer;
  font-size: 12px; padding: 0 0.45em; margin-left: 0.3em; line-height: 1.6; }
#asst .hd button:hover { background: #e9def3; }
#asst .log { overflow-y: auto; padding: 0.5em 0.8em; flex: 1; }
#asst .q { margin: 0.5em 0 0.2em; font-weight: 600; color: #333; }
#asst .a { white-space: pre-wrap; margin: 0.1em 0 0.3em; }
#asst .a.err { color: #b03a2e; }
#asst .links a { display: block; color: #1a5ca8; font-size: 0.9em; margin: 0.15em 0; }
#asst .notice { font-size: 0.78em; color: #8a6d00; background: #fff8e1; border-left: 3px solid #e0b100; padding: 0.2em 0.5em; margin: 0.3em 0 0.6em; }
#asst form { display: flex; gap: 0.4em; padding: 0.5em 0.8em; border-top: 1px solid #eee; }
#asst textarea { flex: 1; resize: none; height: 3.2em; font: inherit; padding: 0.3em; border: 1px solid #ccc; border-radius: 4px; }
#asst button.go { background: #6c3483; color: #fff; border: none; border-radius: 4px; padding: 0 1em; font-weight: 600; cursor: pointer; }
#asst button.go[disabled] { opacity: 0.5; }
#asst .hint { font-size: 0.8em; color: #888; padding: 0 0.8em 0.5em; }`;

  function currentView() {
    const q = Object.fromEntries(new URLSearchParams(location.search));
    return Object.assign({page: here}, q);
  }
  function href(v) {
    const q = new URLSearchParams(v.params).toString();
    if (v.page === 'report') return 'index.html' + (q ? '?' + q : '');
    if (v.page === 'explorer') return 'viewer.html' + (q ? '?' + q : '');
    if (v.page === 'questions') return 'questions.html';
    if (v.page === 'guide') return 'guide.html';
    if (v.page === 'data' && STATUS && STATUS.bin_viewer) return STATUS.bin_viewer + '?' + q + '&back=' + encodeURIComponent(location.href);
    return null;
  }
  function el(tag, cls, text) { const e = document.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; }

  async function ask(question, log, btn) {
    log.appendChild(el('div', 'q', question));
    const a = el('div', 'a', 'Working…'); log.appendChild(a); log.scrollTop = log.scrollHeight;
    btn.disabled = true;
    try {
      const r = await fetch('/assistant/ask', {method: 'POST', headers: {'Content-Type': 'application/json'},
                                               body: JSON.stringify({question, view: currentView()})});
      const out = await r.json();
      if (!r.ok) throw new Error(out.error || r.status);
      a.textContent = out.answer || '(no answer)';
      if (out.views && out.views.length) {
        const ls = el('div', 'links');
        for (const v of out.views) { const h = href(v); if (!h) continue; const l = el('a', null, '→ ' + (v.caption || v.page)); l.href = h; ls.appendChild(l); }
        log.appendChild(ls);
      }
      log.appendChild(el('div', 'notice', out.notice || STATUS.notice));
    } catch (e) { a.textContent = 'Error: ' + e.message; a.classList.add('err'); }
    btn.disabled = false; log.scrollTop = log.scrollHeight;
  }

  async function init() {
    try { const r = await fetch('/assistant/status'); if (!r.ok) return; STATUS = await r.json(); } catch (e) { return; }
    const st = el('style'); st.textContent = css; document.head.appendChild(st);
    const btn = el('button', null, 'Ask a question'); btn.id = 'asst-btn'; document.body.appendChild(btn);
    const box = el('div'); box.id = 'asst';
    box.innerHTML = `<div class="hd"><b>Ask a question</b><span>${STATUS.offline ? 'offline test mode (no AI)' : 'prototype · ' + STATUS.model.replace('us.anthropic.', '')}</span>
      <button class="clr" title="Clear the questions and answers">clear</button><button class="x" title="Close (the button at the bottom right opens it again; Esc also closes)">✕</button></div>
      <div class="log"></div><div class="hint">E.g. “Where does the peak near 2400 m in Argo T against the model come from?” ·
      “Which obs class fits the model worst above 2000 m?” · “Is the XBT misfit a steady offset?” The answer uses this page's
      settings when the question doesn't name them.</div>
      <form><textarea placeholder="Your question about the report's data…"></textarea><button class="go" type="submit">Ask</button></form>`;
    document.body.appendChild(box);
    const log = box.querySelector('.log'), ta = box.querySelector('textarea'), go = box.querySelector('button.go');
    const show = on => { box.classList.toggle('on', on); btn.textContent = on ? 'Close the questions' : 'Ask a question'; if (on) ta.focus(); };
    btn.addEventListener('click', () => show(!box.classList.contains('on')));
    box.querySelector('.hd .x').addEventListener('click', () => show(false));
    box.querySelector('.hd .clr').addEventListener('click', () => { log.innerHTML = ''; });
    document.addEventListener('keydown', e => { if (e.key === 'Escape' && box.classList.contains('on')) show(false); });
    box.querySelector('form').addEventListener('submit', e => { e.preventDefault(); const q = ta.value.trim(); if (!q || go.disabled) return; ta.value = '';
      box.querySelector('.hint').style.display = 'none'; ask(q, log, go); });
    ta.addEventListener('keydown', e => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); box.querySelector('form').requestSubmit(); } });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
