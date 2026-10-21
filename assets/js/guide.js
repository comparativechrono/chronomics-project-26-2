/* =====================================================================
   The project guide in the left pane.

   - The tabs at the top of the pane (#guideTabs, written by
     build/make_page.py) are the stages of the project. The tab of the
     page you are reading is lit, and each tab fills in as you tick the
     checklist of its stage.
   - Checklists: <ul class="checklist"><li data-check="ID">…</li></ul>
     get a box to tick; what you tick is kept in this browser.
   - Buttons of the guide's pages (they work like the other ▶ buttons):
       data-stage-nb="r|stage05-structure|Stage 5 – …"   open, or make,
                                 ~/notebooks/stage05-structure.ipynb
       data-logbook="1"          open ~/project/logbook.md (made the
                                 first time from a template)
   ===================================================================== */
(function () {
  'use strict';
  const MG = window.MG;
  const { h, esc, bus, toast, store, $, $$ } = MG;
  MG.actions = MG.actions || {};
  const KEY = 'guideChecks';
  let checks = store.get(KEY, {}) || {};

  /* ---------------- checklists ---------------- */
  function initChecklists() {
    $$('.help-page ul.checklist > li[data-check]').forEach((li) => {
      if (li.querySelector(':scope > input.ck')) return;
      const id = li.dataset.check;
      const box = h('input.ck', { type: 'checkbox', 'aria-label': 'Done' });
      box.checked = !!checks[id];
      li.classList.toggle('done', box.checked);
      box.addEventListener('change', () => {
        if (box.checked) checks[id] = Date.now();
        else delete checks[id];
        store.set(KEY, checks);
        li.classList.toggle('done', box.checked);
        update();
      });
      li.prepend(box);
    });
  }
  function progressOf(stage) {
    const items = $$(`.help-page[data-stage="${stage}"] ul.checklist > li[data-check]`);
    return { done: items.filter((li) => checks[li.dataset.check]).length, total: items.length };
  }
  function update() {
    $$('#guideTabs .gt[data-stage]').forEach((b) => {
      const p = progressOf(b.dataset.stage);
      const pct = p.total ? Math.round((100 * p.done) / p.total) : 0;
      b.classList.toggle('complete', p.total > 0 && p.done === p.total);
      b.classList.toggle('started', p.done > 0 && p.done < p.total);
      b.title = b.dataset.title + (p.total ? ` – ${p.done} of ${p.total} ticked` : '');
      const bar = b.querySelector('.gt-bar i');
      if (bar) bar.style.width = pct + '%';
    });
    $$('.help-page[data-stage] .stage-progress').forEach((el) => {
      const page = el.closest('.help-page');
      const p = progressOf(page.dataset.stage);
      el.textContent = p.total ? `${p.done} of ${p.total} ticked` : '';
    });
  }

  /* ---------------- the tab of the page that is shown ---------------- */
  function light(id) {
    const page = id && document.getElementById(id);
    if (!page) return;
    const stage = page.dataset.stage;
    const inGuide = page.dataset.group === 'The project';
    let on = null;
    $$('#guideTabs .gt').forEach((b) => {
      const hit = inGuide ? b.dataset.stage === (stage || '0') : !!b.dataset.tools;
      b.classList.toggle('on', hit);
      b.setAttribute('aria-current', hit ? 'page' : 'false');
      if (hit) on = b;
    });
    // keep the lit tab in view when the row of tabs scrolls
    if (on && on.scrollIntoView) on.scrollIntoView({ block: 'nearest', inline: 'nearest' });
  }
  bus.on('help:page', (d) => light(d.id));

  /* ---------------- buttons ---------------- */
  MG.actions.stageNb = async (spec) => {
    const [lang, stem, title] = String(spec).split('|');
    const fs = MG.app.fs;
    const path = `${MG.app.HOME}/notebooks/${stem}.ipynb`;
    await MG.app.showWorkbench('notebook');
    const app = MG.app.nbApp;
    if (fs.exists(path)) return app.open(path);
    const cells = [
      { type: 'markdown', source: `# ${title || stem}\n\n*Your notebook for this stage of the natural short sleep project. Add the guide's code to it with* Add to notebook*, run each cell with Shift+Enter, and write what you find in the text cells – this notebook is part of your record.*` },
      { type: 'code', source: '' }
    ];
    const v = await app.newNotebook(lang === 'python' ? 'python' : 'r', { path, cells });
    toast(`Made <b>${esc(fs.pretty(path))}</b> for this stage.`, 'ok', 3500);
    return v;
  };

  const LOGBOOK = `# My project logbook – short-sleep mutations across the mammals

Keep this file up to date: what you did, what you found, what puzzled you, and what to ask at the next meeting.
Write the date at the top of each entry. Your supervisor may ask to see it, and it is the raw material of your report.

## Questions for my next meeting

-

## Entries

### (date) – Stage 0
- What I read:
- What I did:
- What I found:
- What I do not understand yet:
`;
  MG.actions.logbook = async () => {
    const fs = MG.app.fs;
    const path = `${MG.app.HOME}/project/logbook.md`;
    if (!fs.exists(path)) {
      fs.mkdirp(`${MG.app.HOME}/project`);
      fs.writeText(path, LOGBOOK);
      toast('Made your logbook: <b>~/project/logbook.md</b>.', 'ok', 3500);
    }
    await MG.app.editFile(path);
  };

  /* ---------------- start ---------------- */
  function start() {
    initChecklists();
    update();
    const shown = $('.help-page:not([hidden])');
    if (shown) light(shown.id);
  }
  if (document.body.classList.contains('ready')) start();
  else bus.on('app:ready', start);
})();
