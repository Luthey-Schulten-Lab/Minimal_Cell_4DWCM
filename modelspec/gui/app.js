'use strict';
/* 4DWCM model browser. Data: the ModelSpec JSON embedded in #spec (python -m modelspec build).
   evaluate() is a port of modelspec/impact.py; on load it is checked against the Python results in S.single_knockouts.
   dumpYaml()/loadYaml() write and read the same subset as modelspec/perturbation.py. */

const S = JSON.parse(document.getElementById('spec').textContent);
const $ = (sel, el = document) => el.querySelector(sel);
const store = {
  get(k, d) { try { const v = localStorage.getItem('wcm.' + k); return v === null ? d : JSON.parse(v); } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem('wcm.' + k, JSON.stringify(v)); } catch (e) { /* storage unavailable */ } },
};

function app(el, ...kids) {  // Node.append(null) would insert the text "null"
  for (const c of kids.flat(Infinity)) if (c !== null && c !== undefined && c !== false) el.appendChild(c instanceof Node ? c : document.createTextNode(String(c)));
  return el;
}
function h(tag, attrs, ...kids) {
  const el = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') el.className = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (k === 'html') el.innerHTML = v;
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const c of kids.flat(Infinity)) {
    if (c === null || c === undefined || c === false) continue;
    el.appendChild(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return el;
}
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = v => {
  if (v === null || v === undefined || v === '') return '—';
  if (typeof v !== 'number') return String(v);
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (Number.isInteger(v) && a < 1e12) return v.toLocaleString('en-US');
  if (a >= 1e5 || a < 1e-3) return v.toExponential(3);
  return (+v.toPrecision(5)).toString();
};
function toast(msg) { const t = $('#toast'); t.textContent = msg; t.classList.add('show'); clearTimeout(toast._t); toast._t = setTimeout(() => t.classList.remove('show'), 1800); }
async function copy(text, what) {
  try { await navigator.clipboard.writeText(text); toast('Copied ' + (what || text)); }
  catch (e) { const ta = h('textarea', null, text); document.body.appendChild(ta); ta.select(); document.execCommand('copy'); ta.remove(); toast('Copied ' + (what || '')); }
}

/* ------------------------------------------------------------------ LaTeX (KaTeX from cdnjs; plain text when offline) */
function texSpecies(id) {
  const m = id.match(/^(?:M_)?(.*?)(?:_(c|e))?$/);
  const base = m[1].replace(/__/g, '-').replace(/_/g, '\\_');
  return '[\\mathrm{' + base + '}' + (m[2] === 'e' ? '_{ext}' : '') + ']';
}
function texSym(name, ctx) {  // ctx: {Sub1: species, KmSub1: species, ...}
  if (ctx && ctx[name]) return ctx[name];
  const m = name.match(/^(Km)?(Sub|Prod)(\d+)$/);
  if (m) return (m[1] ? 'K_{' : '') + (m[2] === 'Sub' ? 'S' : 'P') + '_{' + m[3] + '}' + (m[1] ? '}' : '');
  return { kcatF: 'k_{cat}^{F}', kcatR: 'k_{cat}^{R}', Enzyme: '[E]', onoff: '\\mathrm{on/off}', K_uptake: 'k_{uptake}', P_R: 'P_{R}',
    Radius: 'r_{cell}', k_atp: 'k_{ATP}', k_aa: 'k_{aa}', k_tRNA: 'k_{tRNA}', k_cat: 'k_{cat}' }[name] || ('\\mathrm{' + name.replace(/_/g, '\\_') + '}');
}
function toTex(expr, ctx) {
  const toks = expr.replace(/\$/g, '').match(/\d+\.?\d*(?:e[-+]?\d+)?|[A-Za-z_][\w]*|[()+\-*/^]/g) || [];
  let i = 0;
  const peek = () => toks[i], next = () => toks[i++];
  function atom() {
    const t = next();
    if (t === '(') { const e = sum(); next(); return '\\left(' + e + '\\right)'; }
    if (t === '-') return '-' + atom();
    if (/^[\d.]/.test(t)) return t;
    return texSym(t, ctx);
  }
  function pow() { let a = atom(); while (peek() === '^') { next(); a = a + '^{' + atom() + '}'; } return a; }
  function prod() {
    let a = pow();
    while (peek() === '*' || peek() === '/') { const op = next(); const b = pow(); a = op === '*' ? a + ' \\cdot ' + b : '\\frac{' + a + '}{' + b + '}'; }
    return a;
  }
  function sum() { let a = prod(); while (peek() === '+' || peek() === '-') { const op = next(); a = a + ' ' + op + ' ' + prod(); } return a; }
  return sum();
}
function mmTex(r) {  // Rxns_ODE.Enzymatic with the real species substituted
  const S = [], P = [];
  r.substrates.forEach(([id, n]) => { for (let k = 0; k < n; k++) S.push(texSpecies(id)); });
  r.products.forEach(([id, n]) => { for (let k = 0; k < n; k++) P.push(texSpecies(id)); });
  const km = (arr, tag) => arr.map((sp, k) => '\\frac{' + sp + '}{K_{' + tag + ',' + (k + 1) + '}}');
  const num = 'k_{cat}^{F} ' + km(S, 'S').join(' ') + (P.length ? ' - k_{cat}^{R} ' + km(P, 'P').join(' ') : '');
  const den = km(S, 'S').map(x => '\\left(1 + ' + x + '\\right)').join('') + (P.length ? ' + ' + km(P, 'P').map(x => '\\left(1 + ' + x + '\\right)').join('') + ' - 1' : '');
  return 'v = \\mathrm{on/off} \\cdot [E] \\cdot \\frac{' + num + '}{' + den + '}';
}
function customTex(r) {
  const ctx = {};
  r.params.forEach(p => { if (/^(Sub|Prod)\d+$/.test(p.name)) ctx[p.name] = texSpecies(p.value); });
  r.substrates.forEach(([id], k) => { ctx['Sub' + (k + 1)] = texSpecies(id); });
  r.products.forEach(([id], k) => { ctx['Prod' + (k + 1)] = texSpecies(id); });
  return 'v = ' + toTex(r.rate_law, ctx);
}
function trnaTex(r) {
  const E = '\\mathrm{' + r.synthetase.replace('_', '\\_') + '}', aa = texSpecies(r.amino_acid);
  const c = (...parts) => E + '{\\cdot}\\mathrm{' + parts.join('{\\cdot}') + '}';
  return [E + ' + \\mathrm{ATP} \\xrightarrow{k_{ATP}} ' + c('ATP'),
    c('ATP') + ' + ' + aa + ' \\xrightarrow{k_{aa}} ' + c('aa'),
    c('aa') + ' + \\mathrm{tRNA} \\xrightarrow{k_{tRNA}} ' + c('aa', 'tRNA'),
    c('aa', 'tRNA') + ' \\xrightarrow{k_{cat}} ' + E + ' + \\mathrm{aa\\text{-}tRNA} + \\mathrm{AMP} + \\mathrm{PP_i}'].join(' \\\\ ');
}
function texBlock(tex, plain) {
  const el = h('div', { class: 'tex' });
  if (window.katex) { try { katex.render(tex, el, { displayMode: true, throwOnError: true }); return el; } catch (e) { console.warn('katex', e.message); } }
  el.appendChild(h('div', { class: 'plain' }, plain || tex));
  return el;
}
function legendTable(rows) {  // rows: [symbolTex, meaning, linkId|null, valueText|null]
  return h('div', { class: 'card', style: 'margin-top:8px' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Symbol'), h('th', null, 'Meaning'), h('th', null, 'Value / link'))),
    h('tbody', null, rows.map(([sym, meaning, link, val]) => h('tr', null, h('td', null, texInline(sym)), h('td', { class: 'muted' }, meaning),
      h('td', null, val !== null && val !== undefined ? [fmt(val), ' '] : null, link ? (link.startsWith('gip:') ? h('a', { class: 'link', onclick: () => go('constants', link.slice(4)) }, link.slice(4) + (RD.gip_constants[link.slice(4)] ? ' = ' + fmt(RD.gip_constants[link.slice(4)].value) + ' ' + (RD.gip_constants[link.slice(4)].unit || '') : '')) : lnk(link)) : null))))));
}
function texInline(tex) { const el = h('span'); if (window.katex) { try { katex.render(tex, el, { throwOnError: true }); return el; } catch (e) { /* fall through */ } } el.textContent = tex; return el; }
function legendFor(r) {
  const rows = [];
  const spName = id => (SP[id] && SP[id].name) || id;
  if (r.kind === 'ode_mm') {
    let k = 0; r.substrates.forEach(([id, n]) => { for (let j = 0; j < n; j++) { k++; rows.push([texSpecies(id), 'concentration of ' + spName(id) + ' (mM)', id, null]); rows.push(['K_{S,' + k + '}', 'Michaelis constant of ' + spName(id) + ' (mM)', id, (r.params.find(p => p.key === 'Km:' + id) || {}).value]); } });
    k = 0; r.products.forEach(([id, n]) => { for (let j = 0; j < n; j++) { k++; rows.push([texSpecies(id), 'concentration of ' + spName(id) + ' (mM)', id, null]); rows.push(['K_{P,' + k + '}', 'Michaelis constant of ' + spName(id) + ' (mM)', id, (r.params.find(p => p.key === 'Km:' + id) || {}).value]); } });
    rows.push(['k_{cat}^{F}', 'forward catalytic rate constant (1/s)', null, r.kcatF], ['k_{cat}^{R}', 'reverse catalytic rate constant (1/s)', null, r.kcatR]);
    const enz = r.requires.flatMap(g => g.species);
    rows.push(['[E]', enz.length ? 'enzyme concentration from the protein count (' + (r.gpr === 'or' ? 'sum of isozymes' : r.gpr === 'and' ? 'minimum over subunits' : 'single enzyme') + ')' : 'fixed 0.001 mM ("default" enzyme)', enz[0] || null, null]);
    enz.slice(1).forEach(e => rows.push(['[E]', 'also', e, null]));
    rows.push(['\\mathrm{on/off}', '1, or 0 when the reaction is disabled by a perturbation', null, 1]);
  } else if (r.kind === 'ode_custom') {
    r.substrates.forEach(([id]) => rows.push([texSpecies(id), 'concentration of ' + spName(id) + (id.endsWith('_e') ? ' in the medium, fixed' : '') + ' (mM)', id, null]));
    r.products.forEach(([id]) => rows.push([texSpecies(id), 'concentration of ' + spName(id) + ' (mM)', id, null]));
    r.params.forEach(p => rows.push([texSym(p.key), { kcatF: 'forward rate constant', kcatR: 'reverse rate constant', K_uptake: 'fixed uptake rate (mM/s)', P_R: 'membrane permeability (m/s)', Radius: 'cell radius (m)' }[p.key] || p.name, null, p.value]));
  } else if (r.kind === 'cme_trna') {
    rows.push(['E', 'the synthetase ' + spName(r.synthetase), r.synthetase, null], [texSpecies(r.amino_acid), spName(r.amino_acid) + ' (count)', r.amino_acid, null]);
    (r.trnas || []).forEach(t => rows.push(['\\mathrm{tRNA}', 'free tRNA ' + t + '; charged form ' + t + '_ch', t, null]));
    r.params.forEach(p => rows.push([texSym(p.key), { k_atp: 'ATP binding to the synthetase', k_aa: 'amino-acid binding', k_tRNA: 'tRNA binding', k_cat: 'transfer and release' }[p.key] + ' (' + (p.unit || '') + ')', null, p.value]));
  }
  return rows.length ? legendTable(rows) : null;
}
function rateLawBlock(r) {
  if (r.kind === 'ode_mm') return r.rate_law ? texBlock(mmTex(r), r.rate_law) : null;
  if (r.kind === 'ode_custom') return texBlock(customTex(r), r.rate_law);
  if (r.kind === 'cme_trna') return texBlock('\\begin{aligned}' + trnaTex(r) + '\\end{aligned}', r.rate_law);
  return null;
}

/* ------------------------------------------------------------------ indexes */
const GENES = S.genes, SP = S.species, RX = S.reactions, PR = S.processes;
const RD = Object.assign({ families: {}, cme_families: {}, rate_constants: {}, per_gene_constants: {}, diffusion_constants: {}, diffusion_profiles: {}, species_diffusion: {}, gip_constants: {}, regions: [] }, S.rdme || {});
const PGC_OF = {};  // per-gene constant name -> [template, instance]
for (const [t, pg] of Object.entries(RD.per_gene_constants)) for (const i of pg.instances) PGC_OF[i.name] = [t, i];
const FAM = Object.assign({}, RD.families, RD.cme_families);
const GENE_RATE_KEYS = { translation: ['rdme:translation'], rnap_binding: Object.keys(FAM).filter(k => FAM[k].family === 'rnap_binding' || FAM[k].family === 'rnap_binding_next'),
  mrna_degradation: ['rdme:mrna_degradation'], secy_insertion: ['rdme:secy_insertion'], transcription: ['cme:transcription', 'cme:transcription_long'] };
const FAM_RATE_KEY = {}; for (const [k, fs] of Object.entries(GENE_RATE_KEYS)) for (const f of fs) FAM_RATE_KEY[f] = k;
const GIP_EDITABLE = ['rnaPolKcat', 'rnaPolKd', 'riboKcat', 'riboKd', 'ctRNAconc', 'ATPconc', 'UTPconc', 'CTPconc', 'GTPconc'];
function famName(id) { return FAM[id] ? FAM[id].name : id; }
function subst(t, vars) { return t.replace(/<n>/g, vars.n || '<n>').replace(/<i\+1>/g, vars.i !== undefined ? String(vars.i + 1) : '<i+1>').replace(/<i>/g, vars.i !== undefined ? String(vars.i) : '<i>'); }
function instEq(f, inst) {
  const subs = inst.subs || f.template_subs.map(t => subst(t, inst.vars)), prods = inst.prods || f.template_prods.map(t => subst(t, inst.vars));
  const one = x => resolveId(x) ? lnk(x) : x;
  return [subs.map((x, i) => [i ? ' + ' : '', one(x)]), ' → ', prods.map((x, i) => [i ? ' + ' : '', one(x)])];
}
function rateUnit(order) { return order === 2 ? '1/M/s' : '1/s'; }
const geneList = Object.values(GENES).sort((a, b) => a.num.localeCompare(b.num));
const byNum = {}; geneList.forEach(g => { byNum[g.num] = g; });
function geneSpecies(locus) { const g = GENES[locus]; return g.type === 'protein' ? ['R_' + g.num, 'P_' + g.num] : ['R_' + g.num]; }
function speciesGene(id) { const s = SP[id]; return s && s.gene ? s.gene : null; }
function rname(id) { return RX[id] ? (RX[id].name || id) : id; }
function label(id) {
  const s = SP[id];
  if (s && s.gene) { const g = GENES[s.gene]; return (g.gene_name || g.symbol || '') ; }
  return s && s.name && s.name !== id ? s.name : '';
}
const PRODUCERS = {}, CONSUMERS = {};
for (const r of Object.values(RX)) {
  if (!r.active) continue;
  for (const [s] of r.products) (PRODUCERS[s] = PRODUCERS[s] || []).push(r.id);
  for (const [s] of r.substrates) (CONSUMERS[s] = CONSUMERS[s] || []).push(r.id);
  if (r.reversible) {
    for (const [s] of r.substrates) (PRODUCERS[s] = PRODUCERS[s] || []).push(r.id);
    for (const [s] of r.products) (CONSUMERS[s] = CONSUMERS[s] || []).push(r.id);
  }
}
function producers(sid) {  // impact.producers
  return Object.values(RX).filter(r => r.active && (
    (r.products.some(([s]) => s === sid) && r.kcatF !== 0) || (r.substrates.some(([s]) => s === sid) && r.reversible))).map(r => r.id);
}
const PRODUCER_CACHE = {};
for (const [sid, sp] of Object.entries(SP)) if (sp.kind === 'metabolite' || sp.kind === 'protein_form') PRODUCER_CACHE[sid] = producers(sid);

/* ------------------------------------------------------------------ impact (port of modelspec/impact.py) */
const RANK = { ok: 0, reduced: 1, blocked: 2 };
function groups(reqs, lost, red) {
  let status = 'ok'; const why = [];
  for (const g of reqs) {
    const m = g.species, gone = m.filter(s => lost.has(s)), weak = m.filter(s => red.has(s) && !lost.has(s));
    let st;
    if (g.rule === 'all' && gone.length) { st = g.on_loss || 'blocked'; why.push(`needs all of ${m.join(', ')}; lost ${gone.join(', ')}`); }
    else if (g.rule === 'any' && gone.length && gone.length === m.length) { st = g.on_loss || 'blocked'; why.push(`needs one of ${m.join(', ')}; all lost`); }
    else if (gone.length || weak.length) {
      st = 'reduced';
      why.push(gone.length ? `${gone.join(', ')} lost, ${m.filter(s => !lost.has(s)).join(', ')} remain` : `${weak.join(', ')} reduced`);
    } else continue;
    if (RANK[st] > RANK[status]) status = st;
  }
  return [status, why];
}
function evaluate(lost, red, disabled) {
  lost = new Set(lost); red = new Set(red || []);
  const rx = {}, pr = {};
  for (const id of disabled || []) rx[id] = { status: 'blocked', why: ['disabled'] };
  for (const [id, r] of Object.entries(RX)) { if (!r.active || rx[id]) continue; const [st, why] = groups(r.requires, lost, red); if (st !== 'ok') rx[id] = { status: st, why }; }
  for (const [id, p] of Object.entries(PR)) { const [st, why] = groups(p.requires, lost, red); if (st !== 'ok') pr[id] = { status: st, why }; }
  let changed = true;
  while (changed) {
    changed = false;
    for (const [id, p] of Object.entries(PR)) for (const dep of (p.depends_on || [])) {
      let src, lbl;
      if (dep.process) { src = pr[dep.process]; lbl = 'process ' + dep.process; }
      else {
        const hit = Object.entries(rx).filter(([rid, v]) => v.status === 'blocked' && RX[rid].kind === dep.reaction_kind).map(([rid]) => rid).sort();
        src = hit.length ? { status: 'blocked' } : null; lbl = 'reaction ' + hit.join(', ');
      }
      if (!src) continue;
      const st = src.status === 'blocked' ? dep.effect : 'reduced';
      const cur = pr[id], reason = `depends on ${lbl} (${src.status})`;
      if (!cur || RANK[st] > RANK[cur.status]) { pr[id] = { status: st, why: (cur ? cur.why : []).concat([reason]) }; changed = true; }
      else if (!cur.why.includes(reason)) cur.why.push(reason);
    }
  }
  const blocked = new Set(Object.entries(rx).filter(([, v]) => v.status === 'blocked').map(([k]) => k));
  const mets = {};
  if (blocked.size) for (const [sid, prods] of Object.entries(PRODUCER_CACHE))
    if (prods.length && prods.every(p => blocked.has(p))) mets[sid] = { status: 'blocked', why: ['every producing reaction is blocked: ' + prods.join(', ')] };
  return { reactions: rx, processes: pr, metabolites: mets };
}
const KO1 = {};
function singleKO(locus) { return KO1[locus] || (KO1[locus] = evaluate(geneSpecies(locus))); }
function impactClass(locus) {
  const r = singleKO(locus), all = [...Object.values(r.reactions), ...Object.values(r.processes)];
  if (all.some(v => v.status === 'blocked')) return 'blocks';
  if (all.length) return 'reduces';
  return 'none';
}
const SELFCHECK = (() => {
  if (!S.single_knockouts) return { ok: 0, bad: [], n: 0 };
  const bad = [];
  for (const [locus, py] of Object.entries(S.single_knockouts)) {
    const js = singleKO(locus);
    for (const k of ['reactions', 'processes', 'metabolites']) {
      const a = Object.fromEntries(Object.entries(js[k]).map(([i, v]) => [i, v.status]));
      if (JSON.stringify(Object.entries(a).sort()) !== JSON.stringify(Object.entries(py[k]).sort())) { bad.push(locus + ' ' + k); break; }
    }
  }
  return { n: Object.keys(S.single_knockouts).length, bad };
})();

/* ------------------------------------------------------------------ perturbation state */
const MODES = {
  full: 'Full knockout: no transcription, no protein or RNA at t=0',
  expression_only: 'Block expression: t=0 protein and mRNA remain and are only diluted (no protein degradation in the model)',
  initial_only: 'Deplete at t=0 only: the gene is still transcribed, so the protein is made again',
};
const MODE_HELP = {
  full: ['Full knockout', 'The gene is switched off and its product is gone from the start.',
    ['Promoter strength set to 0, so RNAP never transcribes it.', 'Protein and mRNA (or tRNA) start at 0.', 'Hardcoded roles that need the protein (e.g. DnaA at oriC) never fire.'],
    'The standard in-silico knockout.'],
  expression_only: ['Expression only', 'The gene is switched off at t=0, but the cell starts with its normal protein and mRNA.',
    ['Promoter strength set to 0: no new mRNA, so no new protein.', 'The existing mRNA is degraded normally.', 'The model has no protein degradation, so the existing protein only dilutes as the cell grows and divides.'],
    'Models a cell that loses the gene after it has already made the product; within one cycle the effect can be small.'],
  initial_only: ['Initial only', 'The product is removed at t=0 but the gene is left on.',
    ['Protein and mRNA start at 0.', 'The promoter is untouched (with its floor of 45), so the gene is transcribed and the protein is made again.'],
    'Not a knockout: a depletion-and-recovery experiment.'],
};
function infoTip(mode) {  // hover or focus shows the floating bubble (#tip-layer), so panel scroll areas can't clip it
  const [title, lead] = MODE_HELP[mode];
  return h('span', { class: 'tip', tabindex: 0, 'data-mode': mode, 'aria-label': title + ': ' + lead, 'aria-describedby': 'tip-layer' }, h('span', { class: 'tip-i', 'aria-hidden': 'true' }, 'i'));
}
function tipLayer() {
  let el = document.getElementById('tip-layer');
  if (!el) { el = h('div', { id: 'tip-layer', class: 'tip-bubble', role: 'tooltip' }); document.body.appendChild(el); }
  return el;
}
function showTip(anchor) {
  const [title, lead, points, foot] = MODE_HELP[anchor.dataset.mode];
  const el = tipLayer(); el.innerHTML = '';
  app(el, h('strong', null, title), h('span', { class: 'tip-lead' }, lead), h('ul', null, points.map(x => h('li', null, x))), h('span', { class: 'tip-foot' }, foot));
  el.style.display = 'block';
  const a = anchor.getBoundingClientRect(), b = el.getBoundingClientRect(), pad = 8;
  let left = a.left + a.width / 2 - b.width / 2;
  left = Math.max(pad, Math.min(window.innerWidth - b.width - pad, left));
  let top = a.bottom + 6;
  if (top + b.height > window.innerHeight - pad) top = Math.max(pad, a.top - b.height - 6);   // flip above near the bottom
  el.style.left = left + 'px'; el.style.top = top + 'px';
}
function hideTip() { const el = document.getElementById('tip-layer'); if (el) el.style.display = 'none'; }
document.addEventListener('mouseover', e => { const t = e.target.closest && e.target.closest('.tip'); if (t) showTip(t); });
document.addEventListener('mouseout', e => { const t = e.target.closest && e.target.closest('.tip'); if (t && !t.contains(e.relatedTarget)) hideTip(); });
document.addEventListener('focusin', e => { if (e.target.classList && e.target.classList.contains('tip')) showTip(e.target); });
document.addEventListener('focusout', e => { if (e.target.classList && e.target.classList.contains('tip')) hideTip(); });
document.addEventListener('keydown', e => { if (e.key === 'Escape') hideTip(); });
document.addEventListener('scroll', hideTip, true);
function modeTips() {  // one tip per mode, for a legend
  return h('span', { class: 'mode-legend' }, Object.keys(MODE_HELP).map(m => h('span', { class: 'mode-item' }, m.replace('_', ' '), ' ', infoTip(m))));
}
function emptyPert() {
  return { schema: S.schema, name: 'perturbation', description: '', knockouts: [], knockdowns: [], initial_protein_counts: {},
    initial_mrna_means: {}, initial_metabolites_mM: {}, medium_mM: {}, reaction_parameters: [], disabled_reactions: [], rdme_rate_constants: [], gene_rate_scales: [], diffusion_scales: [], gip_constants: [] };
}
let PERT = Object.assign(emptyPert(), store.get('pert', {}));
function savePert() { store.set('pert', PERT); renderCart(); refreshDetail(); }
function setListEntry(key, match, entry) {  // replace-or-remove an entry in one of the list sections
  PERT[key] = (PERT[key] || []).filter(e => !match(e));
  if (entry) PERT[key].push(entry);
  savePert();
}
function koOf(locus) { return PERT.knockouts.find(k => k.gene === locus); }
function kdOf(locus) { return PERT.knockdowns.find(k => k.gene === locus); }
function setKO(locus, mode) {
  PERT.knockdowns = PERT.knockdowns.filter(k => k.gene !== locus);
  const k = koOf(locus);
  if (!mode) PERT.knockouts = PERT.knockouts.filter(x => x.gene !== locus);
  else if (k) k.mode = mode; else PERT.knockouts.push({ gene: locus, mode });
  if (mode === 'full') { const g = GENES[locus]; delete PERT.initial_protein_counts['P_' + g.num]; delete PERT.initial_mrna_means['R_' + g.num]; }
  if (!PERT.name || PERT.name === 'perturbation' || /^ko_/.test(PERT.name)) PERT.name = autoName();
  savePert();
}
function setKD(locus, scale) {
  PERT.knockouts = PERT.knockouts.filter(k => k.gene !== locus);
  PERT.knockdowns = PERT.knockdowns.filter(k => k.gene !== locus);
  if (scale !== null) PERT.knockdowns.push({ gene: locus, promoter_scale: scale });
  savePert();
}
function setMap(key, id, v) { if (v === null) delete PERT[key][id]; else PERT[key][id] = v; savePert(); }
function setParam(rid, key, mode, v) {
  PERT.reaction_parameters = PERT.reaction_parameters.filter(p => !(p.reaction === rid && p.parameter === key));
  if (v !== null) PERT.reaction_parameters.push(mode === 'scale' ? { reaction: rid, parameter: key, scale: v } : { reaction: rid, parameter: key, value: v });
  savePert();
}
function toggleDisabled(rid) {
  const i = PERT.disabled_reactions.indexOf(rid);
  if (i >= 0) PERT.disabled_reactions.splice(i, 1); else PERT.disabled_reactions.push(rid);
  savePert();
}
function autoName() {
  const k = PERT.knockouts.map(x => x.gene.replace('JCVISYN3A_', ''));
  return k.length ? ('ko_' + k.join('_')).slice(0, 80) : 'perturbation';
}
function pertLostReduced() {  // perturbation.resolve: lost / reduced species
  const lost = new Set(), red = new Set();
  for (const k of PERT.knockouts) {
    const sp = geneSpecies(k.gene);
    if (k.mode !== 'initial_only') sp.forEach(s => lost.add(s)); else sp.forEach(s => red.add(s));
  }
  for (const k of PERT.knockdowns) if (k.promoter_scale < 1) geneSpecies(k.gene).forEach(s => red.add(s));
  for (const [sid, v] of Object.entries(PERT.initial_protein_counts)) if (SP[sid] && SP[sid].initial_count && v < SP[sid].initial_count) red.add(sid);
  for (const e of PERT.gene_rate_scales || []) if (e.gene !== '*' && GENES[e.gene] && e.scale < 1) geneSpecies(e.gene).forEach(x => red.add(x));
  return [lost, red];
}
function pertImpact() {
  const [lost, red] = pertLostReduced();
  return evaluate(lost, red, PERT.disabled_reactions);
}
function pertSize() {
  return PERT.knockouts.length + PERT.knockdowns.length + PERT.reaction_parameters.length + PERT.disabled_reactions.length +
    ['rdme_rate_constants', 'gene_rate_scales', 'diffusion_scales', 'gip_constants'].reduce((a, k) => a + (PERT[k] || []).length, 0) +
    ['initial_protein_counts', 'initial_mrna_means', 'initial_metabolites_mM', 'medium_mM'].reduce((a, k) => a + Object.keys(PERT[k]).length, 0);
}
function pertForExport() {
  const p = { schema: S.schema, name: PERT.name || 'perturbation', description: PERT.description || '',
    created: new Date().toISOString().slice(0, 16) + 'Z', model: { fingerprint: S.meta.fingerprint, git_commit: S.meta.git_commit } };
  for (const k of ['knockouts', 'knockdowns', 'initial_protein_counts', 'initial_mrna_means', 'initial_metabolites_mM', 'medium_mM', 'reaction_parameters', 'disabled_reactions',
    'rdme_rate_constants', 'gene_rate_scales', 'diffusion_scales', 'gip_constants'])
    p[k] = PERT[k] || [];
  return p;
}

/* ------------------------------------------------------------------ YAML (same subset as perturbation.dump / loads) */
const TOP = ['schema', 'name', 'description', 'created', 'model', 'knockouts', 'knockdowns', 'reaction_parameters', 'disabled_reactions',
  'rdme_rate_constants', 'gene_rate_scales', 'diffusion_scales', 'gip_constants', 'initial_protein_counts', 'initial_mrna_means', 'initial_metabolites_mM', 'medium_mM'];
function ys(v) {
  if (v === null || v === undefined) return 'null';
  if (typeof v === 'boolean') return v ? 'true' : 'false';
  if (typeof v === 'number') return Number.isInteger(v) && !Object.is(v, -0) ? String(v) : (String(v).includes('.') || String(v).includes('e') ? String(v) : v.toFixed(1));
  return JSON.stringify(String(v));
}
function dumpYaml(p) {
  const out = [`# 4DWCM perturbation, modelspec schema ${S.schema}. Check with: python -m modelspec check <this file>`];
  for (const k of TOP) {
    if (!(k in p)) continue;
    const v = p[k];
    if (Array.isArray(v)) {
      if (!v.length) { out.push(`${k}: []`); continue; }
      out.push(`${k}:`);
      for (const it of v) {
        if (it && typeof it === 'object') Object.entries(it).forEach(([kk, vv], i) => out.push(`${i ? '    ' : '  - '}${kk}: ${ys(vv)}`));
        else out.push(`  - ${ys(it)}`);
      }
    } else if (v && typeof v === 'object') {
      if (!Object.keys(v).length) { out.push(`${k}: {}`); continue; }
      out.push(`${k}:`);
      for (const [kk, vv] of Object.entries(v)) out.push(`  ${kk}: ${ys(vv)}`);
    } else out.push(`${k}: ${ys(v)}`);
  }
  return out.join('\n') + '\n';
}
function yScalar(s) {
  s = s.trim();
  if (s === '' || s === 'null' || s === '~') return null;
  if (s === '{}') return {};
  if (s === '[]') return [];
  if (s === 'true' || s === 'false') return s === 'true';
  if (s[0] === '"') return JSON.parse(s);
  if (s[0] === "'") return s.slice(1, -1).replace(/''/g, "'");
  if (/^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(s)) return Number(s);
  return s;
}
function stripComment(line) {
  let q = null, out = '';
  for (const ch of line) { if (q) { if (ch === q) q = null; } else if (ch === '"' || ch === "'") q = ch; else if (ch === '#') break; out += ch; }
  return out.replace(/\s+$/, '');
}
function loadYaml(text) {
  const root = {}; let key = null, cur = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = stripComment(raw);
    if (!line.trim()) continue;
    const ind = line.length - line.trimStart().length, body = line.trim();
    if (ind === 0) {
      const i = body.indexOf(':'); if (i < 0) throw new Error('cannot read line: ' + raw);
      key = body.slice(0, i).trim(); const v = body.slice(i + 1); root[key] = v.trim() ? yScalar(v) : null; cur = null;
      if (root[key] && typeof root[key] === 'string' && root[key].startsWith('{') && root[key] !== '{}') throw new Error('flow maps are not supported here; use block style (python -m modelspec check reads both)');
    } else if (body.startsWith('- ')) {
      if (root[key] === null) root[key] = [];
      const item = body.slice(2);
      if (/^[A-Za-z_]\w*\s*:/.test(item)) { const i = item.indexOf(':'); cur = { [item.slice(0, i).trim()]: yScalar(item.slice(i + 1)) }; root[key].push(cur); }
      else { root[key].push(yScalar(item)); cur = null; }
    } else {
      const i = body.indexOf(':'); const k = body.slice(0, i).trim(), v = yScalar(body.slice(i + 1));
      if (cur && ind >= 4) cur[k] = v; else { if (root[key] === null) root[key] = {}; root[key][k] = v; }
    }
  }
  return root;
}
function validateLocal(p) {  // quick checks; the server (python -m modelspec serve / check) is authoritative
  const e = [];
  for (const k of p.knockouts || []) { if (!GENES[k.gene]) e.push('unknown gene ' + k.gene); if (!MODES[k.mode || 'full']) e.push('bad mode ' + k.mode); }
  for (const k of p.knockdowns || []) if (!GENES[k.gene]) e.push('unknown gene ' + k.gene);
  for (const [key, kind] of [['initial_protein_counts', 'protein'], ['initial_mrna_means', 'mRNA'], ['initial_metabolites_mM', 'metabolite'], ['medium_mM', 'medium']])
    for (const id of Object.keys(p[key] || {})) if (!SP[id] || SP[id].kind !== kind) e.push(`${key}: ${id} is not a ${kind}`);
  for (const r of p.reaction_parameters || []) if (!RX[r.reaction] || !RX[r.reaction].params.some(x => x.key === r.parameter)) e.push(`unknown parameter ${r.reaction}.${r.parameter}`);
  for (const r of p.disabled_reactions || []) if (!RX[r]) e.push('unknown reaction ' + r);
  for (const x of p.rdme_rate_constants || []) if (!RD.rate_constants[x.constant]) e.push('unknown RDME rate constant ' + x.constant);
  for (const x of p.gene_rate_scales || []) { if (x.gene !== '*' && !GENES[x.gene]) e.push('unknown gene ' + x.gene); if (!GENE_RATE_KEYS[x.rate]) e.push('unknown gene rate ' + x.rate); }
  for (const x of p.diffusion_scales || []) if (!RD.diffusion_constants[x.constant] && x.constant !== 'rna_diff') e.push('unknown diffusion constant ' + x.constant);
  for (const x of p.gip_constants || []) if (!GIP_EDITABLE.includes(x.constant)) e.push('not an editable GIP constant: ' + x.constant);
  if (!/^[A-Za-z0-9_.\-]{1,80}$/.test(p.name || '')) e.push('name must be letters, digits, _ . - (max 80)');
  return e;
}

/* ------------------------------------------------------------------ routing and lists */
const TABS = [['overview', 'Overview'], ['genes', 'Genes'], ['rna', 'RNA'], ['proteins', 'Proteins'], ['metabolites', 'Metabolites'], ['reactions', 'Reactions'], ['processes', 'Processes'], ['constants', 'Constants'], ['inputs', 'Inputs']];
const PROTEIN_KINDS = ['protein', 'protein_form', 'complex'], RNA_KINDS = ['mRNA', 'tRNA', 'rRNA'], MET_KINDS = ['metabolite', 'medium'];
let TAB = store.get('tab', 'overview'), SEL = store.get('sel', null), Q = '', FILT = store.get('filt', {});
const HIST = [];
function go(type, id, push = true) {
  if (push && SEL) HIST.push(SEL);
  SEL = { type, id }; store.set('sel', SEL);
  // the list pane stays on the tab the user chose; links only change the detail pane
  renderList(); renderDetail(); $('#detail').scrollTop = 0; $('#detail').focus({ preventScroll: true });
}
function openId(id) {
  const r = resolveId(id);
  if (r) return go(r[0], r[1]);
}
const NO_AUTO = new Set(['ATPase', 'decay', 'conversion', 'zero']);  // ids that are also ordinary words in product names / prose
function resolveId(id) {  // -> [pageType, pageId] or null
  if (GENES[id]) return ['gene', id];
  if (RX[id]) return ['reaction', id];
  if (FAM[id]) return ['family', id];
  if (SP[id]) return ['species', id];
  if (PR[id]) return ['process', id];
  let m;
  if ((m = id.match(/^(?:G|RP|RPM|DM)_(\d{4})(?:_.*)?$/)) && byNum[m[1]]) return ['gene', byNum[m[1]].locus];        // gene copies, RNAP on gene, counters
  if ((m = id.match(/^(?:RB|D|DT)_(\d{4})(?:_.*)?$/)) && SP['R_' + m[1]]) return ['species', 'R_' + m[1]];             // ribosome- or degradosome-bound mRNA
  if ((m = id.match(/^R_(\d{4})_(?:d|ch)$/)) && SP['R_' + m[1]]) return ['species', 'R_' + m[1]];                       // read mRNA, charged tRNA
  if ((m = id.match(/^(?:C_P|S|PM)_(\d{4})$/)) && SP['P_' + m[1]]) return ['species', 'P_' + m[1]];                    // membrane precursor, SecY-bound, made counter
  if ((m = id.match(/^P_(\d{4})_(?:TC|atp|atp_aa.*)$/)) && SP['P_' + m[1]]) return ['species', 'P_' + m[1]];           // translation cost token, synthetase states
  if (id === 'oriC' || id === 'replisome' || /^ori_(?:HA|LA\d|ss|rep2?)_DnaA(?:_\d*)?$/.test(id)) return ['process', 'replication_initiation'];
  if (id === 'decay') return ['family', 'rdme:decay_removal'];
  if (/^profile_\d+$/.test(id) && RD.diffusion_profiles[id]) return ['constants', id];
  if (RD.rate_constants[id] || RD.gip_constants[id] || (RD.diffusion_constants[id] && id !== 'zero') || id === 'rna_diff') return ['constants', id];
  if (PGC_OF[id] || RD.per_gene_constants[id]) return ['constants', id];                                                 // per-gene rate constants
  return null;
}
function lnk(id, text) { return h('a', { class: 'link', onclick: () => openId(id), title: id }, text || id); }
function badge(text, cls) { return h('span', { class: 'badge ' + (cls || '') }, text); }
function statusBadge(st) { return badge(st, 'b-' + st); }
function essBadge(e) { return e ? badge(e === 'Quasiessential' ? 'Quasi-essential' : e, 'b-' + e) : null; }

function renderTabs() {
  const nav = $('#tabs'); nav.innerHTML = '';
  for (const [id, name] of TABS) nav.appendChild(h('button', { role: 'tab', 'aria-selected': TAB === id ? 'true' : 'false', onclick: () => {
    TAB = id; store.set('tab', TAB); if (id === 'overview') SEL = null; if (id === 'inputs') SEL = { type: 'inputs' }; if (id === 'constants') SEL = { type: 'constants', id: 'all' };
    renderTabs(); renderFilters(); renderList(); renderDetail();
  } }, name));
}
function opts(values) { return values.map(v => Array.isArray(v) ? h('option', { value: v[0] }, v[1]) : h('option', { value: v }, v)); }
function renderFilters() {
  const f = $('#filters'); f.innerHTML = '';
  const lt = TAB === 'overview' ? 'genes' : TAB;   // the overview lists genes, so it shares their filters
  const sel = (key, values) => {
    const s = h('select', { onchange: e => { FILT[lt + '.' + key] = e.target.value; store.set('filt', FILT); renderList(); if (lt === 'constants' && SEL && SEL.type === 'constants') renderDetail(); } }, opts(values));
    s.value = FILT[lt + '.' + key] || (Array.isArray(values[0]) ? values[0][0] : values[0]); f.appendChild(s);
  };
  if (lt === 'genes') {
    sel('type', [['', 'All types'], ['protein', 'Protein'], ['tRNA', 'tRNA'], ['rRNA', 'rRNA']]);
    sel('ess', [['', 'Any essentiality'], ['Essential', 'Essential'], ['Quasiessential', 'Quasi-essential'], ['Nonessential', 'Nonessential']]);
    sel('impact', [['', 'Any knockout impact'], ['blocks', 'KO blocks something'], ['reduces', 'KO only reduces'], ['none', 'KO: no modeled role']]);
    sel('loc', [['', 'Any localization'], ...[...new Set(geneList.map(g => g.localization).filter(Boolean))].sort().map(v => [v, v])]);
  } else if (lt === 'rna') {
    sel('type', [['', 'All RNA'], ['mRNA', 'mRNA'], ['tRNA', 'tRNA'], ['rRNA', 'rRNA']]);
    sel('ess', [['', 'Any essentiality (gene)'], ['Essential', 'Essential'], ['Quasiessential', 'Quasi-essential'], ['Nonessential', 'Nonessential']]);
  } else if (lt === 'proteins') {
    sel('kind', [['', 'All'], ['protein', 'Proteins'], ['protein_form', 'Protein forms (ODE)'], ['complex', 'Complexes and intermediates']]);
    sel('ess', [['', 'Any essentiality'], ['Essential', 'Essential'], ['Quasiessential', 'Quasi-essential'], ['Nonessential', 'Nonessential']]);
    sel('loc', [['', 'Any localization'], ...[...new Set(geneList.map(g => g.localization).filter(Boolean))].sort().map(v => [v, v])]);
    sel('impact', [['', 'Any knockout impact'], ['blocks', 'KO blocks something'], ['reduces', 'KO only reduces'], ['none', 'KO: no modeled role']]);
  } else if (lt === 'metabolites') {
    sel('kind', [['', 'All'], ['metabolite', 'Cytoplasm'], ['medium', 'Medium']]);
    sel('use', [['', 'Used and unused'], ['used', 'Used by a reaction'], ['unused', 'Unused']]);
  } else if (lt === 'constants') {
    const cat = constCatalogue(), n = pred => cat.filter(pred).length;
    sel('cat', [['', `All categories (${cat.length})`], ...CONST_CATS.map(([k, v]) => [k, `${v} (${n(c => c.cat === k)})`])]);
    const procs = [...new Set(cat.flatMap(c => c.procs))].filter(k => PR[k]);
    sel('proc', [['', 'Any process'], ...Object.keys(PR).filter(k => procs.includes(k)).map(k => [k, `${PR[k].name} (${n(c => c.procs.includes(k))})`])]);
    sel('state', [['', 'Any status'], ['used', 'Used by a reaction'], ['unused', 'Unused'], ['edited', 'Edited in this perturbation']]);
  } else if (lt === 'reactions') {
    sel('layer', [['', 'All layers'], ['ODE', 'ODE metabolism'], ['CME', 'CME tRNA charging'], ['CMEt', 'CME transcription'], ['RDME', 'RDME (spatial)']]);
    // grouped: spreadsheet subsystems, then processes (which also cover the RDME / CME families)
    const subCount = {}, procCount = {};
    for (const r of Object.values(RX)) { if (r.subsystem) subCount[r.subsystem] = (subCount[r.subsystem] || 0) + 1; if (r.process) procCount[r.process] = (procCount[r.process] || 0) + 1; }
    for (const f of Object.values(FAM)) if (f.process) procCount[f.process] = (procCount[f.process] || 0) + 1;
    const s = h('select', { 'aria-label': 'Subsystem or process', onchange: e => { FILT['reactions.sub'] = e.target.value; store.set('filt', FILT); renderList(); } },
      h('option', { value: '' }, 'Any subsystem or process'),
      h('optgroup', { label: 'Subsystem (spreadsheet reactions)' }, Object.keys(subCount).sort().map(v => h('option', { value: 'sub:' + v }, `${v} (${subCount[v]})`))),
      h('optgroup', { label: 'Process (RDME, CME and spreadsheet)' }, Object.keys(PR).filter(k => procCount[k]).map(k => h('option', { value: 'proc:' + k }, `${PR[k].name} (${procCount[k]})`))));
    const cur = FILT['reactions.sub'] || '';
    s.value = cur && !cur.includes(':') ? 'sub:' + cur : cur;   // older saved filters stored the bare subsystem
    f.appendChild(s);
    sel('act', [['', 'Active and inactive'], ['1', 'Active only'], ['0', 'Not in the model']]);
  }
}
function listItems() {
  const lt = TAB === 'overview' ? 'genes' : TAB, q = Q.toLowerCase(), F = k => FILT[lt + '.' + k] || '';
  const m = (...xs) => !q || xs.some(x => x && String(x).toLowerCase().includes(q));
  if (lt === 'genes') return geneList.filter(g => (!F('type') || g.type === F('type')) && (!F('ess') || g.essentiality === F('ess')) &&
      (!F('loc') || g.localization === F('loc')) && (!F('impact') || impactClass(g.locus) === F('impact')) &&
      m(g.locus, g.num, g.gene_name, g.symbol, g.product, g.product_proteomics, 'P_' + g.num, 'R_' + g.num))
    .map(g => ({ key: 'gene:' + g.locus, on: () => go('gene', g.locus), id: g.locus.replace('JCVISYN3A_', '') + (g.gene_name || g.symbol ? '  ' + (g.gene_name || g.symbol) : ''),
      right: [koOf(g.locus) ? badge('KO', 'b-ko') : kdOf(g.locus) ? badge('KD', 'b-reduced') : null, g.type !== 'protein' ? badge(g.type) : essBadge(g.essentiality)],
      sub: g.product_proteomics || g.product }));
  if (lt === 'rna') return Object.values(SP).filter(s => RNA_KINDS.includes(s.kind) && (!F('type') || s.kind === F('type')) &&
      (!F('ess') || (GENES[s.gene] || {}).essentiality === F('ess')) && m(s.id, s.name, s.gene, (GENES[s.gene] || {}).gene_name, (GENES[s.gene] || {}).product_proteomics, (GENES[s.gene] || {}).product))
    .sort((a, b) => a.id.localeCompare(b.id))
    .map(s => { const g = GENES[s.gene] || {}; return { key: 'species:' + s.id, on: () => go('species', s.id), id: s.id + (g.gene_name || g.symbol ? '  ' + (g.gene_name || g.symbol) : ''),
      right: [koOf(s.gene) ? badge('KO', 'b-ko') : null, badge(s.kind)], sub: s.kind === 'mRNA' ? (g.product_proteomics || g.product) : s.name }; });
  if (lt === 'proteins') return Object.values(SP).filter(s => PROTEIN_KINDS.includes(s.kind) && (!F('kind') || s.kind === F('kind')) &&
      (!F('ess') || (s.kind === 'protein' && (GENES[s.gene] || {}).essentiality === F('ess'))) &&
      (!F('loc') || (s.kind === 'protein' && (GENES[s.gene] || {}).localization === F('loc'))) &&
      (!F('impact') || (s.kind === 'protein' && impactClass(s.gene) === F('impact'))) &&
      m(s.id, s.name, s.gene, s.carrier, (GENES[s.gene] || {}).gene_name, (GENES[s.gene] || {}).product_proteomics))
    .sort((a, b) => PROTEIN_KINDS.indexOf(a.kind) - PROTEIN_KINDS.indexOf(b.kind) || a.id.localeCompare(b.id))
    .map(s => { const g = GENES[s.gene] || {};
      if (s.kind === 'protein') return { key: 'species:' + s.id, on: () => go('species', s.id), id: s.id + (g.gene_name ? '  ' + g.gene_name : ''),
        right: [koOf(s.gene) ? badge('KO', 'b-ko') : kdOf(s.gene) ? badge('KD', 'b-reduced') : null, essBadge(g.essentiality)], sub: g.product_proteomics || g.product };
      if (s.kind === 'protein_form') return { key: 'species:' + s.id, on: () => go('species', s.id), id: s.id, right: [badge('form')],
        sub: 'form of ' + s.carrier + ' ' + ((GENES[speciesGene(s.carrier)] || {}).gene_name || '') };
      const nm = (s.name || s.id).replace(/ assembly intermediate: (16S|23S) rRNA \+ /, ' intermediate: +');
      return { key: 'species:' + s.id, on: () => go('species', s.id), id: nm.length > 40 ? nm.slice(0, 39) + '…' : nm, right: [badge(s.assembly ? s.assembly : 'complex')], sub: s.id }; });
  if (lt === 'metabolites') return Object.values(SP).filter(s => MET_KINDS.includes(s.kind) &&
      (!F('kind') || s.kind === F('kind')) && (!F('use') || (F('use') === 'unused') === !!s.unused) && m(s.id, s.name, s.kegg))
    .sort((a, b) => a.id.localeCompare(b.id))
    .map(s => ({ key: 'species:' + s.id, on: () => go('species', s.id), id: s.id.length > 34 ? s.id.slice(0, 32) + '…' : s.id,
      right: [s.kind === 'medium' ? badge('medium') : s.kind === 'protein_form' ? badge('protein form') : s.kind === 'complex' ? badge('complex') : null,
        (PERT.initial_metabolites_mM[s.id] !== undefined || PERT.medium_mM[s.id] !== undefined) ? badge('edited', 'b-reduced') : null, s.unused ? badge('unused', 'b-unused') : null],
      sub: s.name + (s.initial_mM !== undefined && s.initial_mM !== null ? ` · ${fmt(s.initial_mM)} mM` : s.medium_mM !== undefined ? ` · ${fmt(s.medium_mM)} mM` : '') }));
  if (lt === 'reactions') {
    const subF = F('sub') && !F('sub').includes(':') ? 'sub:' + F('sub') : F('sub');
    const procF = subF.startsWith('proc:') ? subF.slice(5) : '', sysF = subF.startsWith('sub:') ? subF.slice(4) : '';
    const fams = Object.values(FAM).filter(f => (!F('layer') || (F('layer') === 'RDME' && f.layer === 'RDME') || (F('layer') === 'CMEt' && f.layer === 'CME')) && !sysF && (!procF || f.process === procF) && F('act') !== '0' &&
        m(f.id, f.name, f.family, ...f.template_subs, ...f.template_prods, ...(f.instances.length === 1 ? [...(f.instances[0].subs || []), ...(f.instances[0].prods || [])] : [])))
      .map(f => ({ key: 'family:' + f.id, on: () => go('family', f.id), id: f.name.length > 46 ? f.name.slice(0, 45) + '…' : f.name,
        right: [(PERT.rdme_rate_constants || []).some(x => FAM[f.id].instances[0].rate_name === x.constant) || (PERT.gene_rate_scales || []).some(x => x.gene === '*' && GENE_RATE_KEYS[x.rate].includes(f.id)) ? badge('edited', 'b-reduced') : null,
          badge(f.layer + (f.count > 1 ? ' ×' + f.count : ''))],
        sub: (f.process && PR[f.process] ? PR[f.process].name + ' · ' : '') + f.template_subs.join(' + ') + ' → ' + f.template_prods.join(' + ') }));
    if (F('layer') === 'RDME' || F('layer') === 'CMEt') return fams;
    return fams.concat(Object.values(RX).filter(r => (!F('layer') || r.layer === F('layer')) && (!sysF || r.subsystem === sysF) && (!procF || r.process === procF) &&
      (!F('act') || String(+r.active) === F('act')) &&
      m(r.id, r.name, r.subsystem, r.enzyme_str, r.synthetase, r.formula, ...r.substrates.map(x => x[0]), ...r.products.map(x => x[0])))
    .sort((a, b) => a.id.localeCompare(b.id))
    .map(r => ({ key: 'reaction:' + r.id, on: () => go('reaction', r.id), id: r.id, dim: !r.active,
      right: [PERT.disabled_reactions.includes(r.id) ? badge('off', 'b-blocked') : PERT.reaction_parameters.some(p => p.reaction === r.id) ? badge('edited', 'b-reduced') : null,
        !r.active ? badge('not in model', 'b-inactive') : badge(r.layer)],
      sub: r.name || (r.subsystem || '') })));
  }
  if (lt === 'constants') {
    const cat = constCatalogue().filter(c => (!F('cat') || c.cat === F('cat')) && (!F('proc') || c.procs.includes(F('proc'))) &&
      (!F('state') || (F('state') === 'used' ? c.used : F('state') === 'unused' ? !c.used : constEdited(c))) && m(c.key, c.sub));
    const allLabel = F('cat') ? CONST_CATS.find(x => x[0] === F('cat'))[1] : 'All constants';
    return [{ key: 'constants:all', on: () => go('constants', 'all'), id: allLabel, right: [], sub: 'the filtered constants as tables, with perturbation controls' }]
      .concat(cat.map(c => ({ key: 'constants:' + c.key, on: () => go('constants', c.key), id: c.key, dim: !c.used,
        right: [constEdited(c) ? badge('edited', 'b-reduced') : null, badge(c.badge)], sub: c.sub })));
  }
  if (lt === 'processes') return Object.values(PR).filter(p => m(p.id, p.name, p.description, p.layer))
    .map(p => ({ key: 'process:' + p.id, on: () => go('process', p.id), id: p.name, right: [badge(p.layer)], sub: p.scope === 'per_gene' ? 'every gene' : p.description }));
  if (lt === 'inputs') return S.inputs.filter(r => m(r.file, r.sheet))
    .map(r => ({ key: 'inputs', on: () => go('inputs', r.file + '|' + (r.sheet || '')), id: r.sheet || r.file, dim: !r.used,
      right: [r.used ? badge('read', 'b-ok') : badge('unused', 'b-inactive')], sub: r.sheet ? r.file : '' }));
  return [];
}
function renderList() {
  const items = listItems(), ul = $('#list'); ul.innerHTML = '';
  const active = Object.entries(FILT).filter(([k, v]) => v && k.startsWith((TAB === 'overview' ? 'genes' : TAB) + '.')).length;
  $('#count').textContent = items.length + ' shown' + (active ? ` · ${active} filter${active > 1 ? 's' : ''} on` : '');
  if (!items.length) ul.appendChild(h('li', { class: 'dim', style: 'cursor:default' }, h('span', { class: 'sub' }, 'Nothing matches these filters.')));
  const cur = SEL ? SEL.type + ':' + SEL.id : '';
  const frag = document.createDocumentFragment();
  for (const it of items) frag.appendChild(h('li', { onclick: it.on, 'aria-current': it.key === cur ? 'true' : null, class: it.dim ? 'dim' : null },
    h('span', { class: 'id' }, it.id), h('span', null, it.right), h('span', { class: 'sub' }, it.sub)));
  ul.appendChild(frag);
}

/* ------------------------------------------------------------------ developer: code references */
function editorHref(f, l) {
  const tpl = store.get('editor', '');
  if (!tpl) return null;
  return tpl.replace('{abs}', (S.meta.head || '') + '/' + f).replace('{file}', f).replace('{line}', l);
}
function refRow(r) {
  const href = editorHref(r.f, r.l);
  const loc = href ? h('a', { class: 'loc', href }, `${r.f}:${r.l}`) : h('span', { class: 'loc', title: 'Copy path:line', onclick: () => copy(`${r.f}:${r.l}`) }, `${r.f}:${r.l}`);
  return h('div', { class: 'ref' + (r.c ? ' c' : '') }, loc, h('span', { class: 'fn' }, r.fn), r.t ? h('pre', null, r.t) : null);
}
function refList(refs, limit = 40) {
  if (!refs || !refs.length) return h('div', { class: 'faint' }, 'none');
  const live = refs.filter(r => !r.c), dead = refs.filter(r => r.c);
  const box = h('div', null, live.slice(0, limit).map(refRow));
  if (live.length > limit) box.appendChild(h('details', null, h('summary', null, `${live.length - limit} more`), live.slice(limit).map(refRow)));
  if (dead.length) box.appendChild(h('details', null, h('summary', null, `${dead.length} commented-out`), dead.map(refRow)));
  return box;
}
const HARDCODED_NOTE = () => h('div', { class: 'faint', style: 'font-size:12px;margin-bottom:6px' },
  'Lines where this exact id is written as a string literal, so the code gives it a special role that no spreadsheet edit can change. Lines that build names from a variable (for every gene at once) are listed under generic references.');
function devBox(title, ...kids) { return h('div', { class: 'devbox dev-only' }, h('h4', null, title), kids); }
function srcLine(src) {
  if (!src) return null;
  if (src.sheet) return h('div', { class: 'src' }, `${src.file} › ${src.sheet} › row ${src.row}, column "${src.col}"`);
  return h('div', { class: 'src' }, `${src.file}${src.locus_tag ? ' › ' + src.feature + ' ' + src.locus_tag : ''}`);
}
function anchorRow(a) { const r = S.code && S.code.anchors[a]; return r ? refRow({ f: r.f, l: r.l, fn: r.fn + `  (lines ${r.l}–${r.end})`, t: '' }) : h('div', { class: 'faint' }, a); }
function patternRefs(kind) {
  if (!S.code) return null;
  const pre = (S.code.prefix_kinds[kind] || []);
  return pre.map(p => h('details', null, h('summary', null, `"${p}" + id — ${S.code.prefixes[p]} (${(S.code.patterns[p] || []).filter(r => !r.c).length} places build this name for every ${kind === 'metabolite' || kind === 'medium' ? 'metabolite' : 'gene'})`),
    refList(S.code.patterns[p] || [], 200)));
}

/* ------------------------------------------------------------------ network sketch */
function ego(left, center, right, opt = {}) {
  const W = 800, bw = 250, bh = 36, gap = 8, cx = W / 2, maxN = 12;
  const cut = a => a.length > maxN ? a.slice(0, maxN - 1).concat([{ label: `+${a.length - maxN + 1} more`, more: true }]) : a;
  left = cut(left); right = cut(right);
  const n = Math.max(left.length, right.length, 1), H = Math.max(n * (bh + gap) + 40, 120 + (opt.top ? 60 : 0));
  const top0 = opt.top ? 60 : 0;
  const col = (arr, x) => arr.map((d, i) => ({ ...d, x, y: 20 + top0 / 2 + i * (bh + gap) + ((n - arr.length) * (bh + gap)) / 2 }));
  const L = col(left, 10), R = col(right, W - bw - 10), C = { ...center, x: cx - bw / 2, y: H / 2 - bh / 2 + top0 / 4 };
  const svg = [`<svg viewBox="0 0 ${W} ${H + top0 / 2}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="${esc(opt.aria || 'neighborhood')}">`];
  const edge = (x1, y1, x2, y2, cls = '') => svg.push(`<path class="ed ${cls}" d="M${x1},${y1} C${(x1 + x2) / 2},${y1} ${(x1 + x2) / 2},${y2} ${x2},${y2}"/>`);
  for (const d of L) edge(d.x + bw, d.y + bh / 2, C.x, C.y + bh / 2, d.cat ? 'cat' : '');
  for (const d of R) edge(C.x + bw, C.y + bh / 2, d.x, d.y + bh / 2, d.cat ? 'cat' : '');
  if (opt.top) opt.top.forEach((d, i) => { d.x = cx - bw / 2 + (i - (opt.top.length - 1) / 2) * (bw + 10); d.y = 6; edge(d.x + bw / 2, d.y + bh, C.x + bw / 2, C.y, 'cat'); });
  const node = (d, cls) => {
    const t = (s, n) => esc(s.length > n ? s.slice(0, n - 1) + '…' : s);
    svg.push(`<g class="nd ${cls} ${d.status || ''} ${d.id && !d.more ? 'click' : ''}" data-id="${esc(d.id || '')}" transform="translate(${d.x},${d.y})">` +
      `<title>${esc(d.title || d.id || d.label)}</title><rect width="${bw}" height="${bh}" rx="6"/>` +
      `<text x="10" y="${d.sub ? 15 : 22}">${t(d.label, 32)}</text>${d.sub ? `<text class="s" x="10" y="29">${t(d.sub, 42)}</text>` : ''}</g>`);
  };
  L.forEach(d => node(d, '')); R.forEach(d => node(d, '')); (opt.top || []).forEach(d => node(d, '')); node(C, 'center');
  if (opt.lt) svg.push(`<text class="lbl" x="10" y="12">${esc(opt.lt)}</text>`);
  if (opt.rt) svg.push(`<text class="lbl" x="${W - 10}" y="12" text-anchor="end">${esc(opt.rt)}</text>`);
  svg.push('</svg>');
  const wrap = h('div', { class: 'svgwrap', html: svg.join('') });
  wrap.addEventListener('click', e => { const g = e.target.closest('.nd.click'); if (g && g.dataset.id) openId(g.dataset.id); });
  return wrap;
}
function rxNode(rid, imp) {
  const r = RX[rid];
  return { id: rid, label: rid, sub: r.name || r.subsystem || '', title: rid + ' — ' + (r.name || rid) + (r.enzyme_str ? '\nenzyme ' + r.enzyme_str : ''), status: imp && imp.reactions[rid] ? imp.reactions[rid].status : '' };
}

/* ------------------------------------------------------------------ auto-link ids in prose */
const TOKEN_RE = /[A-Za-z0-9_]+(?::[A-Za-z][A-Za-z0-9_]*)?/g;
function linkable(tok) { return tok.length >= 3 && !NO_AUTO.has(tok) && (/[0-9_A-Z]/.test(tok.slice(1)) || (SP[tok] && SP[tok].kind === 'protein_form')) && resolveId(tok) !== null; }
function linkifyTree(root) {
  const skip = new Set(['A', 'INPUT', 'TEXTAREA', 'SELECT', 'SCRIPT', 'STYLE', 'BUTTON', 'OPTION', 'H2']);
  const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, { acceptNode: n => {
    for (let e = n.parentElement; e && e !== root; e = e.parentElement)
      if (skip.has(e.tagName) || e.namespaceURI === 'http://www.w3.org/2000/svg' || (typeof e.className === 'string' && /\b(tex|yaml|katex|loc|ref|mono|tip-bubble|svgwrap)\b/.test(e.className))) return NodeFilter.FILTER_REJECT;
    TOKEN_RE.lastIndex = 0; let m;
    while ((m = TOKEN_RE.exec(n.textContent))) if (linkable(m[0])) return NodeFilter.FILTER_ACCEPT;
    return NodeFilter.FILTER_SKIP; } });
  const nodes = []; let n; while ((n = w.nextNode())) nodes.push(n);
  for (const t of nodes) {
    const frag = document.createDocumentFragment(), text = t.textContent; let last = 0, m; TOKEN_RE.lastIndex = 0;
    while ((m = TOKEN_RE.exec(text))) {
      if (!linkable(m[0])) continue;
      frag.appendChild(document.createTextNode(text.slice(last, m.index))); frag.appendChild(lnk(m[0])); last = m.index + m[0].length;
    }
    if (!last) continue;
    frag.appendChild(document.createTextNode(text.slice(last))); t.parentNode.replaceChild(frag, t);
  }
}

/* ------------------------------------------------------------------ detail pages */
let CURRENT_IMPACT = null;
function refreshDetail() { renderList(); renderDetail(); }
function renderDetail() {
  const d = $('#detail'); d.innerHTML = '';
  CURRENT_IMPACT = pertSize() ? pertImpact() : null;
  if (HIST.length) d.appendChild(h('div', { class: 'crumb' }, h('a', { class: 'link', onclick: () => { const p = HIST.pop(); go(p.type, p.id, false); } }, '← back')));
  if (!SEL || TAB === 'overview' && !SEL) { d.appendChild(pageOverview()); linkifyTree(d); return; }
  const f = { gene: pageGene, species: pageSpecies, reaction: pageReaction, process: pageProcess, inputs: pageInputs, family: pageFamily, constants: pageConstants }[SEL.type];
  try { d.appendChild(f ? f(SEL.id) : pageOverview()); linkifyTree(d); }
  catch (e) { console.error(e); d.appendChild(h('div', { class: 'note bad' }, 'Could not render ' + SEL.id + ': ' + e.message)); }
}
function tile(k, v, n) { return h('div', { class: 'tile' }, h('div', { class: 'k' }, k), h('div', { class: 'v' }, v), n ? h('div', { class: 'n' }, n) : null); }
function impactBlock(imp, heading) {
  const rows = [];
  const add = (kind, id, v) => rows.push(h('tr', null, h('td', null, statusBadge(v.status)), h('td', null, kind),
    h('td', null, kind === 'process' ? lnk(id, PR[id].name) : [lnk(id), kind === 'reaction' ? h('div', { class: 'faint', style: 'font-size:12px' }, rname(id)) : null]),
    h('td', { class: 'muted' }, v.why.join('; '), kind === 'process' && PR[id].consequence ? h('div', { style: 'color:var(--ink);margin-top:2px' }, PR[id].consequence) : null)));
  const order = o => Object.entries(o).sort((a, b) => RANK[b[1].status] - RANK[a[1].status] || a[0].localeCompare(b[0]));
  order(imp.processes).forEach(([id, v]) => add('process', id, v));
  order(imp.reactions).forEach(([id, v]) => add('reaction', id, v));
  order(imp.metabolites).forEach(([id, v]) => add('metabolite', id, v));
  const box = h('div', null, heading ? h('h3', null, heading) : null);
  if (!rows.length) box.appendChild(h('div', { class: 'card muted' }, 'No modeled reaction or process depends on it. Removing it only frees the NTPs and amino acids spent making it.'));
  else box.appendChild(h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Effect'), h('th', null, 'Kind'), h('th', null, 'What'), h('th', null, 'Why'))), h('tbody', null, rows))));
  return box;
}
function rolesTable(sid) {
  const roles = (SP[sid] && SP[sid].roles) || [];
  if (!roles.length) return null;
  return h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Role'), h('th', null, 'In'), h('th', null, 'Layer'), h('th', null, 'Detail'))),
    h('tbody', null, roles.map(r => {
      const t = r.ref.slice(0, r.ref.indexOf(':')), id = r.ref.slice(r.ref.indexOf(':') + 1); const o = t === 'rxn' ? RX[id] : t === 'fam' ? FAM[id] : PR[id];
      if (!o) return null;
      const st = CURRENT_IMPACT && (t === 'rxn' ? CURRENT_IMPACT.reactions[id] : t === 'proc' ? CURRENT_IMPACT.processes[id] : null);
      return h('tr', null, h('td', null, r.role), h('td', null, lnk(id, t === 'rxn' ? id : o.name), ' ', st ? statusBadge(st.status) : null, !o.active && t === 'rxn' ? badge('not in model', 'b-inactive') : null),
        h('td', null, o.layer), h('td', { class: 'muted' }, t === 'rxn' ? (o.name || o.subsystem || '') : t === 'fam' ? (o.regions || []).join(', ') : (o.scope === 'per_gene' ? 'every gene' : '')));
    }))));
}

function pageOverview() {
  const box = h('div');
  const counts = { protein: 0, tRNA: 0, rRNA: 0 }; geneList.forEach(g => counts[g.type]++);
  const act = Object.values(RX).filter(r => r.active);
  app(box, h('div', { class: 'title' }, h('h2', null, 'JCVI-syn3A in the 4DWCM')),
    h('p', { class: 'lede' }, 'What the simulator is built from: every gene, species, reaction and hardcoded role, read from the same inputs the run reads. ' +
      'Pick genes, metabolites or reactions to build a perturbation in the panel on the right; it exports a YAML file that ',
      h('code', null, 'python -m modelspec check'), ' validates.'),
    h('div', { class: 'tiles' }, tile('Genes', geneList.length, `${counts.protein} protein · ${counts.tRNA} tRNA · ${counts.rRNA} rRNA`),
      tile('Metabolites', Object.values(SP).filter(s => s.kind === 'metabolite').length, Object.values(SP).filter(s => s.kind === 'medium').length + ' more held fixed in the medium'),
      tile('Reactions in the model', act.length + (RD.n_reactions || 0) + (RD.n_cme_transcription || 0), `${act.filter(r => r.layer === 'ODE').length} ODE · ${RD.n_reactions || 0} RDME · ${(RD.n_cme_transcription || 0) + act.filter(r => r.layer === 'CME').length} CME`),
      tile('Hardcoded processes', Object.keys(PR).length, 'roles written in Python'),
      tile('Cell volume', fmt(S.meta.volume_L * 1e15) + ' fL', `1 particle = ${fmt(S.meta.mM_per_particle * 1000)} µM`)));

  // essentiality x static knockout impact
  const cls = ['blocks', 'reduces', 'none'], names = { blocks: 'KO blocks a reaction or process', reduces: 'KO only reduces', none: 'no modeled role' };
  const ess = ['Essential', 'Quasiessential', 'Nonessential'], M = {};
  geneList.filter(g => g.type === 'protein').forEach(g => { const k = (g.essentiality || '?') + '|' + impactClass(g.locus); M[k] = (M[k] || 0) + 1; });
  const cell = (e, c) => h('td', { class: 'num' }, h('a', { class: 'link', onclick: () => { TAB = 'genes'; store.set('tab', TAB); FILT['genes.ess'] = e; FILT['genes.impact'] = c; FILT['genes.type'] = 'protein'; store.set('filt', FILT); renderTabs(); renderFilters(); renderList(); } }, M[e + '|' + c] || 0));
  app(box, h('h3', null, 'Single-gene knockouts: static impact vs. experimental essentiality'),
    h('div', { class: 'card' }, h('table', { class: 'matrix' },
      h('thead', null, h('tr', null, h('th', null, 'Essentiality (experiment)'), cls.map(c => h('th', null, names[c])))),
      h('tbody', null, ess.map(e => h('tr', null, h('td', { class: 'rowh' }, essBadge(e)), cls.map(c => cell(e, c))))))),
    h('p', { class: 'muted' }, 'Static impact is structural: which reactions lose their only enzyme and which processes lose a part. ' +
      'Essential genes with no modeled role are the ones the model cannot yet call essential; nonessential genes that block something are candidates for a model disagreement. Click a count to list the genes.'));

  const caveats = [
    ['Promoter floor', 'promoter strength = min(765, max(45, count)): setting a protein\'s initial count to 0 does not stop its gene being expressed. Use a knockout, which sets the promoter to 0.'],
    ['No protein degradation', 'the RDME has mRNA degradation but no protein degradation, so blocking expression only removes a protein by dilution at division.'],
    ['RNAP starts at zero', 'every polymerase is assembled from P_0645 ×2 + P_0804 + P_0803; knocking out any of them stops all transcription.'],
    ['Ribosomes are placed at t=0', 'ribosomal protein and rRNA knockouts stop new ribosomes only; the initial ones keep translating.'],
    ['SMC floor', 'numSmc = max(1, …): knocking out P_0415 leaves one loop extruder.'],
    ['Enzyme-less reactions', `${act.filter(r => r.gpr === 'default').map(r => r.id).join(', ')} use a fixed 0.001 mM enzyme and cannot be knocked out.`],
    ['Other-Random-Binding', 'the sheet is parsed by a function the model never calls (commented out in defineRxns), so its reactions are not in the model.'],
    ['Medium is fixed', 'medium species are constant parameters of the ODE, not state; changing one changes the whole run.'],
  ];
  app(box, h('h3', null, 'Before you knock something out'), h('div', { class: 'card' }, h('table', null, h('tbody', null, caveats.map(([a, b]) => h('tr', null, h('td', { style: 'width:190px;font-weight:600' }, a), h('td', { class: 'muted' }, b)))))));

  // top hub enzymes
  const hub = {};
  for (const r of act) for (const g of r.requires) for (const s of g.species) if (s.startsWith('P_')) hub[s] = (hub[s] || 0) + 1;
  const top = Object.entries(hub).sort((a, b) => b[1] - a[1]).slice(0, 10), mx = top.length ? top[0][1] : 1;
  app(box, h('h3', null, 'Proteins needed by the most reactions'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Protein'), h('th', null, 'Gene'), h('th', null, 'Product'), h('th', null, 'Reactions'))), h('tbody', null, top.map(([p, n]) => {
    const g = GENES[speciesGene(p)];
    return h('tr', null, h('td', { style: 'width:90px' }, lnk(p)), h('td', null, g ? (g.gene_name || '') : ''), h('td', { class: 'muted' }, g ? g.product_proteomics : ''),
      h('td', { style: 'width:170px' }, h('span', { class: 'bar', style: `width:${Math.round(120 * n / mx)}px` }), ' ', n));
  })))));

  app(box, devBox('Build', h('div', { class: 'src' }, `commit ${S.meta.git_commit} on ${S.meta.git_branch} · inputs ${S.meta.fingerprint} · built ${S.meta.built}`),
    h('div', { class: 'src' }, SELFCHECK.bad.length ? `impact engine: ${SELFCHECK.bad.length} of ${SELFCHECK.n} genes differ from Python (${SELFCHECK.bad.slice(0, 4).join(', ')})`
      : `impact engine matches modelspec/impact.py for all ${SELFCHECK.n} single knockouts`),
    h('div', { style: 'margin-top:8px' }, 'Editor link template (optional): ',
      h('input', { style: 'width:420px;max-width:100%', placeholder: 'vscode://vscode-remote/ssh-remote+HOST{abs}:{line}', value: store.get('editor', ''),
        onchange: e => { store.set('editor', e.target.value.trim()); renderDetail(); } })),
    h('div', { class: 'faint', style: 'font-size:12px;margin-top:4px' }, 'With it set, every file:line opens in your editor; without it, clicking copies file:line. {abs} = absolute path, {file} = repo path.'),
    h('h4', { style: 'margin-top:12px' }, 'Reaction id conventions'),
    h('table', null, h('thead', null, h('tr', null, h('th', null, 'Id part'), h('th', null, 'Meaning'))), h('tbody', null, S.id_parts.map(([a, b]) => h('tr', null, h('td', { class: 'mono', style: 'width:200px' }, a), h('td', { class: 'muted' }, b))))),
    h('h4', { style: 'margin-top:12px' }, 'Species naming in the code'),
    h('table', null, h('thead', null, h('tr', null, h('th', null, 'Pattern'), h('th', null, 'Meaning'))), h('tbody', null, S.naming.map(([a, b]) => h('tr', null, h('td', { class: 'mono', style: 'width:200px' }, a), h('td', { class: 'muted' }, b)))))));
  return box;
}

function koControls(g) {
  const ko = koOf(g.locus), kd = kdOf(g.locus);
  const modeSel = h('select', { class: 'btn', onchange: e => setKO(g.locus, e.target.value) }, opts(Object.keys(MODES).map(m => [m, MODES[m].split(':')[0]])));
  if (ko) modeSel.value = ko.mode;
  const kdIn = h('input', { type: 'number', min: 0, step: 0.05, value: kd ? kd.promoter_scale : 0.25, title: 'promoter strength multiplier' });
  const icIn = h('input', { type: 'number', min: 0, step: 1, value: PERT.initial_protein_counts['P_' + g.num] ?? g.initial_protein ?? 0 });
  return h('div', { class: 'actions' },
    ko ? [h('button', { class: 'btn danger on', onclick: () => setKO(g.locus, null) }, 'Knocked out ✕'), modeSel, infoTip(ko.mode)]
      : [h('button', { class: 'btn danger', onclick: () => setKO(g.locus, 'full') }, 'Knock out'), h('span', { class: 'faint', style: 'font-size:12px;align-self:center' }, 'modes: '), modeTips()],
    h('span', { class: 'inline-form' }, h('button', { class: 'btn', onclick: () => setKD(g.locus, +kdIn.value) }, kd ? 'Update knockdown' : 'Knock down ×'), kdIn,
      kd ? h('button', { class: 'btn sm', onclick: () => setKD(g.locus, null) }, 'remove') : null),
    g.type === 'protein' && !(ko && ko.mode === 'full') ? h('span', { class: 'inline-form' },
      h('button', { class: 'btn', onclick: () => { const v = parseInt(icIn.value, 10); setMap('initial_protein_counts', 'P_' + g.num, Number.isFinite(v) && v !== g.initial_protein ? v : null); } }, 'Set initial protein'), icIn) : null);
}
function roleNetwork(g, main) {
  const roles = (SP[main] && SP[main].roles) || [];
  const right = roles.map(r => { const tp = r.ref.slice(0, r.ref.indexOf(':')), id = r.ref.slice(r.ref.indexOf(':') + 1);
    if (tp === 'rxn') return { ...rxNode(id, CURRENT_IMPACT), sub: r.role, title: `${id} — ${rname(id)}\n${r.role}`, cat: /enzyme|isozyme|ODE form|synthetase|subunit/.test(r.role) };
    if (tp === 'fam') return { id, label: FAM[id] ? FAM[id].name : id, sub: r.role + ' · RDME', title: (FAM[id] ? FAM[id].name : id) + '\n' + r.role };
    return { id, label: PR[id].name, sub: r.role, title: PR[id].name + '\n' + r.role, status: CURRENT_IMPACT && CURRENT_IMPACT.processes[id] ? CURRENT_IMPACT.processes[id].status : '' }; });
  const rid = 'R_' + g.num;
  if (main === rid && g.type === 'protein')        // an mRNA's uses are per-gene RDME reactions, not roles
    right.unshift({ id: 'rdme:ribosome_binding', label: 'Translation', sub: 'ribosome binds ' + rid }, { id: 'rdme:degradosome_binding', label: 'mRNA degradation', sub: 'degradosome binds ' + rid });
  const left = [{ id: g.locus, label: g.locus.replace('JCVISYN3A_', 'gene '), sub: 'G_' + g.num + '_C1/C2' }];
  if (main.startsWith('P_')) left.push({ id: rid, label: rid, sub: 'mRNA' }, { label: 'ribosomeP', id: 'ribosomeP', sub: 'translation' });
  else left.push({ label: 'RNAP', id: 'RNAP', sub: 'transcription' });
  return ego(left, { label: main, sub: g.gene_name || g.type, id: null }, right.length ? right : [{ label: 'no modeled role' }], { lt: 'made from', rt: 'used by', aria: 'roles of ' + main });
}
function pageProtein(pid) {
  const s = SP[pid], g = GENES[s.gene], box = h('div'), ko = koOf(g.locus);
  app(box, h('div', { class: 'crumb' }, 'protein · ', lnk(g.locus, g.locus + ' gene page')),
    h('div', { class: 'title' }, h('h2', null, g.gene_name || g.symbol || pid), h('span', { class: 'mono' }, pid),
      essBadge(g.essentiality), g.localization ? badge(g.localization) : null, g.function ? badge(g.function) : null, ko ? badge('knocked out: ' + ko.mode, 'b-ko') : null),
    h('p', { class: 'lede' }, g.product_proteomics || g.product), koControls(g));
  if (ko) app(box, h('div', { class: 'note' + (ko.mode === 'initial_only' ? '' : ' info') }, MODES[ko.mode]));
  const rd = RD.species_diffusion[pid], dc = rd ? RD.diffusion_profiles[rd.profile] : null;
  app(box, h('div', { class: 'tiles' },
    tile('Initial count', fmt(g.initial_protein), g.initial_rule || 'placed at t=0' + (g.sim_count !== g.initial_protein ? ` (${fmt(g.sim_count)} in the sheet)` : '')),
    tile('Measured count', fmt(g.exp_count), 'proteomics (Breuer et al.)'),
    tile('Length', fmt(g.aa_len) + ' aa', ''),
    tile('Made by', 'translation', 'of ' + 'R_' + g.num + (g.membrane_insertion ? ', then SecY insertion' : '')),
    tile('Degradation', 'none', 'no protein degradation in the model; diluted at division')));
  if (g.membrane_insertion) app(box, h('div', { class: 'note info' }, `Trans-membrane: translated as the precursor C_P_${g.num} in the cytoplasm, bound by SecY as S_${g.num}, then inserted as ${pid} in the membrane.`));
  if (s.carries) app(box, h('div', { class: 'note info' }, 'Also exists in the ODE as ', s.carries.map((f, i) => [i ? ', ' : '', lnk(f)]), ' (protein_metabolites.xlsx); reactions on those forms need this protein.'));
  app(box, h('h3', null, 'Where it acts'), roleNetwork(g, pid));
  const rt = rolesTable(pid); if (rt) app(box, rt);
  app(box, geneReactionsTable(g, f => /translation|ribosome_binding|secy/.test(f.family)));
  app(box, diffusionBlock(pid));
  app(box, impactBlock(singleKO(g.locus), 'If its gene is knocked out (full)'));
  if (S.code) app(box, devBox('Hardcoded in the source', HARDCODED_NOTE(), refList([pid, 'C_P_' + g.num].flatMap(i => S.code.refs[i] || []))),
    devBox('Generic references (built for every gene)', patternRefs('protein')));
  return box;
}
function pageRNA(rid) {
  const s = SP[rid], g = GENES[s.gene], box = h('div'), ko = koOf(g.locus);
  const title = s.kind === 'mRNA' ? 'mRNA of ' + (g.gene_name || g.symbol || g.locus) : (g.product || rid);
  app(box, h('div', { class: 'crumb' }, s.kind + ' · ', lnk(g.locus, g.locus + ' gene page')),
    h('div', { class: 'title' }, h('h2', null, title), h('span', { class: 'mono' }, rid), badge(s.kind), essBadge(g.essentiality), ko ? badge('knocked out: ' + ko.mode, 'b-ko') : null),
    h('p', { class: 'lede' }, s.kind === 'mRNA' ? 'Encodes ' + (g.product_proteomics || g.product) + ' (' + lnk('P_' + g.num).textContent + ').' : s.kind === 'tRNA' ? `Carries ${g.trna_aa} to the ribosome; charged in the global CME by ${g.trna_aa}TRS.` : 'Structural RNA of the ribosome.'),
    koControls(g));
  if (ko) app(box, h('div', { class: 'note' + (ko.mode === 'initial_only' ? '' : ' info') }, MODES[ko.mode]));
  const d = RD.species_diffusion[rid], own = d ? Object.entries(d.own).find(([k]) => !k.endsWith('Dna')) : null;
  const t = [];
  if (s.kind === 'mRNA') t.push(tile('At t=0', fmt(g.mrna_initial_mean), 'Poisson mean (2 × average count)'), tile('Average count', fmt(g.mrna_avg), 'Diffusion Coefficient sheet'));
  else t.push(tile('At t=0', fmt(g.initial_rna ?? 0), s.kind === 'tRNA' ? 'initializeTRNA places 199' : 'made only by transcription'));
  t.push(tile('Length', fmt(g.rna_len) + ' nt', g.long_transcription ? 'several RNAP at once' : ''), tile('Promoter strength', fmt(g.promoter), g.promoter_rule));
  if (own) t.push(tile('Diffusion', fmt(own[1]) + ' m²/s', 'from its length (Stokes–Einstein); half in the DNA region'));
  if (s.kind === 'mRNA') { const deg = (FAM['rdme:mrna_degradation'] || { instances: [] }).instances.find(i => i.gene === g.locus); if (deg) t.push(tile('Degradation', fmt(deg.rate) + ' 1/s', 'once bound by a degradosome (88 nt/s ÷ length)')); }
  if (s.kind === 'tRNA') t.push(tile('Amino acid', g.trna_aa, 'charged form ' + rid + '_ch'));
  app(box, h('div', { class: 'tiles' }, t));
  app(box, h('h3', null, 'Where it acts'), roleNetwork(g, rid));
  const rt = rolesTable(rid); if (rt) app(box, rt);
  app(box, geneReactionsTable(g, f => !/translation$|secy/.test(f.family) || f.family === 'ribosome_binding' || f.family === 'mrna_recycle'));
  app(box, diffusionBlock(rid));
  if (s.kind !== 'mRNA') app(box, impactBlock(singleKO(g.locus), 'If its gene is knocked out (full)'));
  if (S.code) app(box, devBox('Hardcoded in the source', HARDCODED_NOTE(), refList([rid, 'G_' + g.num + '_C1'].flatMap(i => S.code.refs[i] || []))),
    devBox('Generic references (built for every gene)', patternRefs(s.kind)));
  return box;
}
function pageGene(locus) {
  const g = GENES[locus], box = h('div'), pid = 'P_' + g.num, rid = 'R_' + g.num;
  const ko = koOf(locus);
  app(box, h('div', { class: 'crumb' }, `${g.type} gene · ${g.locus}`),
    h('div', { class: 'title' }, h('h2', null, g.gene_name || g.symbol || g.locus), h('span', { class: 'mono' }, g.type === 'protein' ? `${pid} · ${rid}` : rid),
      essBadge(g.essentiality), g.localization ? badge(g.localization) : null, g.function ? badge(g.function) : null, ko ? badge('knocked out: ' + ko.mode, 'b-ko') : null),
    h('p', { class: 'lede' }, g.product_proteomics || g.product),
    koControls(g));
  if (ko) box.appendChild(h('div', { class: 'note' + (ko.mode === 'initial_only' ? '' : ' info') }, MODES[ko.mode]));
  const t = [];
  if (g.type === 'protein') {
    t.push(tile('Initial protein', fmt(g.initial_protein), g.initial_rule || (g.sim_count !== g.initial_protein ? `${fmt(g.sim_count)} in the sheet; integer split over compartments` : 'placed at t=0')));
    t.push(tile('Measured count', fmt(g.exp_count), 'proteomics (Breuer et al.)'));
    t.push(tile('Promoter strength', fmt(g.promoter), g.promoter_rule));
    t.push(tile('mRNA at t=0', fmt(g.mrna_initial_mean), 'Poisson mean (2 × average count)'));
    t.push(tile('Length', fmt(g.aa_len) + ' aa', `${fmt(g.rna_len)} nt`));
  } else {
    t.push(tile('Initial RNA', fmt(g.initial_rna ?? 0), g.type === 'tRNA' ? 'initializeTRNA places 199' : 'made only by transcription'));
    t.push(tile('Promoter strength', fmt(g.promoter), g.promoter_rule));
    t.push(tile('Length', fmt(g.rna_len) + ' nt', g.long_transcription ? 'multi-RNAP transcription' : ''));
    if (g.trna_aa) t.push(tile('Amino acid', g.trna_aa, 'charged by ' + g.trna_aa + 'TRS'));
  }
  t.push(tile('Location', `${fmt(g.start)}–${fmt(g.end)}`, g.strand === 1 ? 'forward strand' : 'reverse strand'));
  box.appendChild(h('div', { class: 'tiles' }, t));
  if (g.type === 'protein' && g.sim_count < 45) box.appendChild(h('div', { class: 'note' }, `Promoter floor applies: the sheet count is ${g.sim_count}, but the promoter strength is 45. Setting the initial count to 0 does not stop expression.`));

  // network: gene -> RNA -> protein -> roles
  const main = g.type === 'protein' ? pid : rid;
  app(box, h('div', { class: 'actions', style: 'margin-top:0' }, g.type === 'protein' ? h('button', { class: 'btn sm', onclick: () => go('species', pid) }, 'Protein page ' + pid) : null,
    h('button', { class: 'btn sm', onclick: () => go('species', rid) }, (g.type === 'protein' ? 'mRNA' : g.type) + ' page ' + rid)));
  app(box, h('h3', null, 'Where it acts'), roleNetwork(g, main));
  const rt = rolesTable(main); if (rt) box.appendChild(rt);
  if (g.carries || (SP[pid] && SP[pid].carries)) box.appendChild(h('div', { class: 'note info' }, `This protein also exists as ODE species ${SP[pid].carries.join(', ')} (protein_metabolites.xlsx); reactions on those forms need it.`));

  app(box, geneReactionsTable(g));
  app(box, diffusionBlock(main));
  app(box, impactBlock(singleKO(locus), 'If knocked out (full)'));
  if (g.orthologs) app(box, h('h3', null, 'Counts in other organisms'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Organism'), h('th', null, 'Protein count'))), h('tbody', null,
    Object.entries(g.orthologs).map(([k, v]) => h('tr', null, h('td', null, k), h('td', { class: 'num' }, fmt(v))))))));

  if (S.code) {
    const lit = [pid, rid, 'G_' + g.num + '_C1', 'C_P_' + g.num].flatMap(i => (S.code.refs[i] || []).map(r => ({ ...r, id: i })));
    const procs = Object.values(PR).filter(p => p.scope === 'per_gene' && (p.id !== 'membrane_insertion' || g.membrane_insertion) && (g.type === 'protein' || p.id === 'transcription'));
    app(box, devBox('Inputs', srcLine(g.src), srcLine(g.src_counts), srcLine(g.src_mrna)),
      devBox('Species built for this gene', h('div', { class: 'mono', style: 'font-size:12px' },
        (g.type === 'protein' ? [`G_${g.num}_C1`, `G_${g.num}_C2`, `RP_${g.num}_C1`, `RP_${g.num}_C2`, rid, rid + '_d', `RB_${g.num}`, `D_${g.num}`, `DT_${g.num}`, g.membrane_insertion ? `C_P_${g.num}` : null, g.membrane_insertion ? `S_${g.num}` : null, pid, pid + '_TC', `PM_${g.num}`, `RPM_${g.num}`, `DM_${g.num}`]
          : [`G_${g.num}_C1`, `G_${g.num}_C2`, `RP_${g.num}_C1`, `RP_${g.num}_C2`, rid, g.type === 'tRNA' ? rid + '_ch' : null, `RPM_${g.num}`]).filter(Boolean).join('  '))),
      devBox('Hardcoded in the source', HARDCODED_NOTE(), refList(lit)),
      devBox('Per-gene processes it goes through', procs.map(p => h('details', null, h('summary', null, p.name), (p.anchors || []).map(anchorRow)))),
      devBox('Generic references (built for every gene)', patternRefs(g.type === 'protein' ? 'protein' : g.type)));
  }
  return box;
}

function pageSpecies(id) {
  const s = SP[id], box = h('div');
  if (!s) return h('div', null, 'unknown species ' + id);
  if (s.kind === 'protein' && GENES[s.gene]) return pageProtein(id);
  if (RNA_KINDS.includes(s.kind) && GENES[s.gene]) return pageRNA(id);
  app(box, h('div', { class: 'crumb' }, s.kind.replace('_', ' ')), h('div', { class: 'title' }, h('h2', null, s.name || id), h('span', { class: 'mono' }, id),
    s.kegg ? badge('KEGG ' + s.kegg) : null, s.carrier ? h('span', null, 'form of ', lnk(s.carrier)) : null));
  if (s.note) box.appendChild(h('p', { class: 'lede' }, s.note));
  if (s.unused) box.appendChild(h('div', { class: 'note' }, 'Not used: ' + s.unused));
  const t = [];
  if (s.kind === 'metabolite' && s.initial_mM !== undefined) {
    const edited = PERT.initial_metabolites_mM[id];
    t.push(tile('Initial concentration', fmt(s.initial_mM) + ' mM', edited !== undefined ? `perturbed to ${fmt(edited)} mM` : 'Intracellular Metabolites'));
    t.push(tile('Initial particles', fmt(s.initial_count), 'round(mM × Nₐ × V)'));
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: edited ?? s.initial_mM });
    box.appendChild(h('div', { class: 'actions' }, h('span', { class: 'inline-form' }, h('button', { class: 'btn', onclick: () => setMap('initial_metabolites_mM', id, +inp.value === s.initial_mM ? null : +inp.value) }, 'Set initial mM'), inp),
      edited !== undefined ? h('button', { class: 'btn sm', onclick: () => setMap('initial_metabolites_mM', id, null) }, 'reset') : null));
  }
  if (s.kind === 'medium') {
    const edited = PERT.medium_mM[id];
    t.push(tile('Medium concentration', fmt(s.medium_mM) + ' mM', edited !== undefined ? `perturbed to ${fmt(edited)} mM` : 'Simulation Medium, held fixed'));
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: edited ?? s.medium_mM });
    box.appendChild(h('div', { class: 'actions' }, h('span', { class: 'inline-form' }, h('button', { class: 'btn', onclick: () => setMap('medium_mM', id, +inp.value === s.medium_mM ? null : +inp.value) }, 'Set medium mM'), inp),
      h('button', { class: 'btn danger', onclick: () => setMap('medium_mM', id, 0) }, 'Remove from medium'),
      edited !== undefined ? h('button', { class: 'btn sm', onclick: () => setMap('medium_mM', id, null) }, 'reset') : null));
  }
  if (s.kind === 'complex') t.push(tile('At t=0', fmt(s.initial_count), ''));
  if (t.length) box.appendChild(h('div', { class: 'tiles' }, t));
  if (s.kind === 'complex' && s.note && !s.assembly) box.appendChild(h('div', { class: 'note info' }, s.note));
  if (s.kind === 'complex' && s.assembly) app(box, h('p', { class: 'muted' }, 'Part of ', lnk(s.assembly === 'SSU' ? 'ssu_assembly' : 'lsu_assembly', s.assembly === 'SSU' ? '30S subunit assembly' : '50S subunit assembly'), '.'));
  const prods = [...new Set(PRODUCERS[id] || [])], cons = [...new Set(CONSUMERS[id] || [])];
  if (prods.length || cons.length) app(box, h('h3', null, 'Reactions'), ego(prods.map(r => rxNode(r, CURRENT_IMPACT)), { label: id, sub: s.name, id: null }, cons.map(r => rxNode(r, CURRENT_IMPACT)),
    { lt: 'produced by', rt: 'consumed by', aria: 'reactions of ' + id }), h('p', { class: 'faint', style: 'font-size:12px' }, 'Reversible reactions appear on both sides. Enzymes are under each reaction name.'));
  const rt = rolesTable(id); if (rt) app(box, h('h3', null, 'Roles'), rt);
  app(box, diffusionBlock(id));
  if (CURRENT_IMPACT && CURRENT_IMPACT.metabolites[id]) box.appendChild(h('div', { class: 'note bad' }, 'With the current perturbation: ' + CURRENT_IMPACT.metabolites[id].why[0]));
  if (S.code) app(box, devBox('Input', srcLine(s.src)), devBox('Hardcoded in the source', HARDCODED_NOTE(), refList(S.code.refs[id])),
    ['metabolite', 'medium'].includes(s.kind) ? devBox('Generic references (built for every metabolite)', patternRefs('metabolite')) : null);
  return box;
}

function side(list) { return list.map(([s, st], i) => [i ? ' + ' : '', st !== 1 ? `${st} ` : '', lnk(s)]); }
function pageReaction(rid) {
  const r = RX[rid], box = h('div'), off = PERT.disabled_reactions.includes(rid);
  const st = CURRENT_IMPACT && CURRENT_IMPACT.reactions[rid];
  app(box, h('div', { class: 'crumb' }, `${r.layer} reaction · ${r.sheet}`),
    h('div', { class: 'title' }, h('h2', null, r.name || rid), h('span', { class: 'mono' }, rid), r.subsystem ? badge(r.subsystem) : null, !r.active ? badge('not in model', 'b-inactive') : null, st ? statusBadge(st.status) : null),
    h('div', { class: 'eq' }, side(r.substrates), r.reversible ? '  ⇌  ' : '  →  ', side(r.products)));
  if (r.process && PR[r.process]) app(box, h('p', { class: 'muted' }, 'Part of the process ', lnk(r.process, PR[r.process].name), '.'));
  if (r.note) box.appendChild(h('div', { class: 'note' }, r.note));
  if (r.warning) box.appendChild(h('div', { class: 'note' }, r.warning));
  if (st) box.appendChild(h('div', { class: 'note bad' }, 'With the current perturbation: ' + st.why.join('; ')));
  if (r.active) box.appendChild(h('div', { class: 'actions' }, h('button', { class: 'btn danger' + (off ? ' on' : ''), onclick: () => toggleDisabled(rid) }, off ? 'Disabled ✕' : 'Disable reaction')));
  // enzyme
  const enz = r.requires.length ? r.requires.map(g => h('div', null, g.via ? 'carrier (protein metabolite): ' : g.rule === 'any' ? 'any of (isozymes, counts summed): ' : g.species.length > 1 ? 'all of (complex, minimum count): ' : 'enzyme: ',
    g.species.map((s, i) => [i ? ', ' : '', lnk(s), label(s) ? ` ${label(s)}` : '']))) : h('div', { class: 'muted' }, r.gpr === 'default' ? 'no enzyme: a fixed 0.001 mM ("default") — cannot be knocked out' : 'none');
  app(box, h('h3', null, 'Needs'), h('div', { class: 'card' }, enz));
  if (r.kind === 'cme_trna') app(box, h('h3', null, 'Charges'), h('div', { class: 'card' }, (r.trnas || []).map((t, i) => [i ? ', ' : '', lnk(t)]), ' with ', lnk(r.amino_acid)));
  // network
  const [lostNow] = pertLostReduced();
  const top = r.requires.flatMap(g => g.species).slice(0, 3).map(s => ({ id: s, label: s, sub: label(s), status: lostNow.has(s) ? 'blocked' : '' }));
  app(box, h('h3', null, 'Species'), ego(r.substrates.map(([s, n]) => ({ id: s, label: s, sub: (SP[s] && SP[s].name) || '' })), { label: rid, sub: r.subsystem },
    r.products.map(([s, n]) => ({ id: s, label: s, sub: (SP[s] && SP[s].name) || '' })), { top, lt: 'substrates', rt: 'products', aria: 'species of ' + rid }));
  const law = rateLawBlock(r);
  if (law) app(box, h('h3', null, 'Rate law'), law, legendFor(r), r.kind === 'ode_mm' ? h('p', { class: 'faint', style: 'font-size:12px' }, 'Random-binding reversible Michaelis–Menten (Rxns_ODE.Enzymatic); [E] is the enzyme concentration from the protein count; on/off is 1 unless the reaction is disabled.')
    : r.kind === 'cme_trna' ? h('p', { class: 'faint', style: 'font-size:12px' }, 'Four mass-action steps in the global CME (Rxns_CME.tRNAcharging), one chain per tRNA of this amino acid.') : null);
  // params
  const rows = r.params.map(p => {
    const cur = PERT.reaction_parameters.find(x => x.reaction === rid && x.parameter === p.key);
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? (cur.scale ?? cur.value) : 1, style: 'width:80px' });
    const mode = h('select', null, opts([['scale', '×'], ['value', '=']])); if (cur && cur.value !== undefined) mode.value = 'value';
    return h('tr', null, h('td', null, p.name, p.species ? [' ', lnk(p.species)] : null), h('td', { class: 'num' }, fmt(p.value)), h('td', { class: 'muted' }, p.unit || ''),
      h('td', null, p.key && r.active && typeof p.value === 'number' ? h('span', { class: 'inline-form' }, mode, inp,
        h('button', { class: 'btn sm', onclick: () => { const v = +inp.value; setParam(rid, p.key, mode.value, (mode.value === 'scale' && v === 1) || (mode.value === 'value' && v === p.value) ? null : v); } }, cur ? 'update' : 'set'),
        cur ? h('button', { class: 'btn sm', onclick: () => setParam(rid, p.key, 'scale', null) }, 'reset') : null) : null),
      h('td', { class: 'dev-only src' }, p.src ? `row ${p.src.row}` : '', p.stored_as_text ? h('div', { class: 'faint' }, 'text cell') : null));
  });
  app(box, h('h3', null, 'Parameters'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Parameter'), h('th', null, 'Value'), h('th', null, 'Unit'), h('th', null, 'Perturb'), h('th', { class: 'dev-only' }, 'Sheet row'))), h('tbody', null, rows))));
  if (S.code) {
    const builder = { ode_mm: r.active ? 'processes/Rxns_ODE.py::defineRandomBindingRxns' : null, ode_custom: 'processes/Rxns_ODE.py::defineNonRandomBindingRxns', cme_trna: 'processes/Rxns_CME.py::tRNAcharging' }[r.kind];
    app(box, devBox('Input', h('div', { class: 'src' }, `input_data/kinetic_params.xlsx › ${r.sheet}` + (r.params[0] && r.params[0].src ? ` › rows ${Math.min(...r.params.map(p => p.src.row))}–${Math.max(...r.params.map(p => p.src.row))}` : '')),
      r.sbml_id ? h('div', { class: 'src' }, `input_data/Syn3A_updated.xml › reaction ${r.sbml_id} (matched by name)`) : null),
      devBox('Built by', builder ? anchorRow(builder) : h('div', { class: 'faint' }, 'no builder: not in the model'), r.kind === 'ode_mm' ? anchorRow('processes/Rxns_ODE.py::getEnzymeConc') : null),
      devBox('Hardcoded in the source', HARDCODED_NOTE(), refList(S.code.refs[rid])));
  }
  return box;
}

function pageProcess(pid) {
  const p = PR[pid], box = h('div'), st = CURRENT_IMPACT && CURRENT_IMPACT.processes[pid];
  app(box, h('div', { class: 'crumb' }, `${p.layer} · ${p.scope === 'per_gene' ? 'runs for every gene' : 'global'}`), h('div', { class: 'title' }, h('h2', null, p.name), st ? statusBadge(st.status) : null),
    h('p', { class: 'lede' }, p.description));
  if (p.consequence) box.appendChild(h('div', { class: 'note' }, 'If it fails: ' + p.consequence));
  if (st) box.appendChild(h('div', { class: 'note bad' }, 'With the current perturbation: ' + st.why.join('; ')));
  if (p.requires.length) app(box, h('h3', null, 'Needs'), h('div', { class: 'card' }, p.requires.map(g => h('div', { style: 'margin:3px 0' }, g.rule === 'any' ? 'one of: ' : 'all of: ',
    g.species.map((s, i) => [i ? ', ' : '', lnk(s), label(s) ? h('span', { class: 'faint' }, ' ' + label(s)) : null]), g.on_loss ? h('span', { class: 'faint' }, `  (loss only ${g.on_loss})`) : null))));
  if ((p.depends_on || []).length) app(box, h('h3', null, 'Depends on'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Needs'), h('th', null, 'If that is blocked, this is'), h('th', null, 'Basis'))), h('tbody', null, p.depends_on.map(d => h('tr', null,
    h('td', null, d.process ? lnk(d.process, PR[d.process].name) : `any ${d.reaction_kind} reaction`), h('td', null, statusBadge(d.effect)),
    h('td', null, d.verified ? h('span', { class: 'verified' }, '✓ checked in code: ' + d.verified) : h('span', { class: 'unverified' }, 'inferred, not checked in code'))))))));
  const roles = Object.entries(p.species_roles || {});
  if (roles.length) app(box, h('h3', null, 'Species'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Species'), h('th', null, 'Gene / name'), h('th', null, 'Role'))), h('tbody', null, roles.map(([s, role]) => h('tr', null, h('td', { style: 'width:30%' }, lnk(s)), h('td', { class: 'muted' }, label(s)), h('td', null, role)))))));
  const hits = geneList.filter(g => { const v = singleKO(g.locus).processes[pid]; return v; });
  if (hits.length) app(box, h('h3', null, 'Single knockouts that hit it'), h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Locus'), h('th', null, 'Gene'), h('th', null, 'Product'), h('th', null, 'Effect on this process'))), h('tbody', null, hits.map(g => h('tr', null,
    h('td', { style: 'width:120px' }, lnk(g.locus, g.locus.replace('JCVISYN3A_', ''))), h('td', null, g.gene_name || ''), h('td', { class: 'muted' }, g.product_proteomics || g.product), h('td', null, statusBadge(singleKO(g.locus).processes[pid].status))))))));
  if (p.applies_to) {
    const rows = p.applies_to.map(l => GENES[l]).map(g => h('tr', null, h('td', { style: 'width:120px' }, lnk(g.locus, g.locus.replace('JCVISYN3A_', ''))),
      h('td', null, g.gene_name || ''), h('td', { class: 'muted' }, g.product_proteomics || g.product), h('td', null, essBadge(g.essentiality))));
    app(box, h('h3', null, `${p.applies_to.length} genes go through this process`), h('div', { class: 'card' }, h('table', null,
      h('thead', null, h('tr', null, h('th', null, 'Locus'), h('th', null, 'Gene'), h('th', null, 'Product'), h('th', null, 'Essentiality'))), h('tbody', null, rows))));
  }
  if ((p.families || []).length) app(box, h('h3', null, `${p.families.length} RDME / CME reaction ${p.families.length > 1 ? 'types' : 'type'}`), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Reaction'), h('th', null, 'Layer'), h('th', null, 'Instances'), h('th', null, 'Template'))),
    h('tbody', null, p.families.map(fid => { const f = FAM[fid]; return h('tr', null, h('td', null, lnk(fid, f.name)), h('td', null, f.layer), h('td', { class: 'num' }, f.count),
      h('td', { class: 'muted mono', style: 'font-size:12px' }, f.template_subs.join(' + ') + ' → ' + f.template_prods.join(' + '))); })))));
  if (p.reactions) app(box, h('h3', null, `${p.reactions.length} reactions`), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Reaction'), h('th', null, 'Name'), h('th', null, 'Enzyme'))), h('tbody', null, p.reactions.map(rid => h('tr', null,
      h('td', null, lnk(rid)), h('td', { class: 'muted' }, rname(rid)), h('td', null, RX[rid].requires.flatMap(g => g.species).map((x, i) => [i ? ', ' : '', lnk(x)]))))))));
  box.appendChild(devBox('Implemented in', (p.anchors || []).map(anchorRow)));
  if ((p.code_species || []).length) box.appendChild(devBox('Species named inside these functions', h('div', null, p.code_species.map((x, i) => [i ? ', ' : '', lnk(x)]))));
  return box;
}

const CONST_CATS = [['shared', 'Shared rate constants'], ['assembly', 'Ribosome assembly rate constants'], ['per_gene', 'Per-gene rate constants'],
  ['diffusion', 'Diffusion coefficients'], ['profile', 'Diffusion profiles'], ['gip', 'GIP rate-formula constants']];
const GIP_PROCESS = { riboKcat: ['translation'], riboKd: ['translation'], ctRNAconc: ['translation'], ribosomeConc: [],
  rnaPolKcat: ['transcription'], rnaPolKd: ['transcription'], rrnaPolKcat: ['transcription'], ATPconc: ['transcription'], UTPconc: ['transcription'],
  CTPconc: ['transcription'], GTPconc: ['transcription'], Ecoli_V: ['transcription', 'translation', 'mrna_degradation'], avgdr: ['transcription', 'translation', 'mrna_degradation'], countToMiliMol: [] };
function famProcs(fids) { return [...new Set(fids.map(f => FAM[f] && FAM[f].process).filter(Boolean))]; }
function speciesProcs(prof) {  // processes whose RDME families use a species with this diffusion profile
  const out = new Set();
  for (const f of Object.values(RD.families)) for (const d of f.diffusion || []) if ((d.profiles || []).includes(prof) && f.process) out.add(f.process);
  return [...out];
}
let CONST_CAT_CACHE = null;
function constCatalogue() {
  if (CONST_CAT_CACHE) return CONST_CAT_CACHE;
  const out = [];
  for (const [k, v] of Object.entries(RD.rate_constants)) out.push({ key: k, cat: v.group === 'assembly' ? 'assembly' : 'shared', procs: famProcs(v.used_by), used: v.used_by.length > 0,
    sub: fmt(v.value) + ' ' + v.unit + (v.used_by.length ? '' : ' · unused'), badge: v.group === 'assembly' ? 'assembly' : 'rate' });
  for (const [k, v] of Object.entries(RD.per_gene_constants)) out.push({ key: k, cat: 'per_gene', procs: famProcs(v.used_by), used: true, rate: v.gene_rate,
    sub: v.what + ' · ' + fmt(v.min) + '–' + fmt(v.max) + ' ' + v.unit, badge: 'per gene ×' + v.count });
  const dcUsers = {};
  for (const [pid, p] of Object.entries(RD.diffusion_profiles)) for (const [, , dc] of p.transitions) (dcUsers[dc] = dcUsers[dc] || new Set()).add(pid);
  for (const [k, v] of Object.entries(RD.diffusion_constants)) { if (k === 'zero') continue;
    const profs = [...(dcUsers[k] || [])]; out.push({ key: k, cat: 'diffusion', procs: [...new Set(profs.flatMap(speciesProcs))], used: profs.length > 0, sub: fmt(v.value) + ' m²/s', badge: 'diffusion' }); }
  out.push({ key: 'rna_diff', cat: 'diffusion', procs: [...new Set(Object.keys(RD.diffusion_profiles).filter(p => RD.diffusion_profiles[p].transitions.some(t => t[2].startsWith('rna_diff'))).flatMap(speciesProcs))],
    used: true, sub: `${RD.rna_diff_count || ''} RNAs and ribosome intermediates, ${fmt((RD.rna_diff_range || [0])[0])}–${fmt((RD.rna_diff_range || [0, 0])[1])} m²/s`, badge: 'diffusion' });
  for (const [k, p] of Object.entries(RD.diffusion_profiles)) out.push({ key: k, cat: 'profile', procs: speciesProcs(k), used: p.members > 0 && k !== 'profile_1',
    sub: `${p.members} species · ${p.kinds.join(', ')}`, badge: 'profile' });
  for (const [k, v] of Object.entries(RD.gip_constants)) out.push({ key: k, cat: 'gip', procs: GIP_PROCESS[k] || [], used: (GIP_PROCESS[k] || []).length > 0,
    sub: fmt(v.value) + ' ' + (v.unit || '') + (GIP_EDITABLE.includes(k) ? '' : ' · derived'), badge: 'GIP' });
  return (CONST_CAT_CACHE = out);
}
function constEdited(c) {
  if (c.cat === 'shared' || c.cat === 'assembly') return (PERT.rdme_rate_constants || []).some(x => x.constant === c.key);
  if (c.cat === 'per_gene') return (PERT.gene_rate_scales || []).some(x => x.rate === c.rate);
  if (c.cat === 'diffusion') return (PERT.diffusion_scales || []).some(x => x.constant === c.key);
  if (c.cat === 'gip') return (PERT.gip_constants || []).some(x => x.constant === c.key);
  return false;
}
function rcLink(name) {  // any rate constant name -> link to its row on the Constants tab
  return (RD.rate_constants[name] || PGC_OF[name]) ? h('a', { class: 'link mono', style: 'font-size:12px', onclick: () => go('constants', name) }, name) : h('span', { class: 'mono' }, name);
}
function diffusionTable(f) {
  if (!(f.diffusion || []).length) return null;
  return h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Species'), h('th', null, 'Diffusion'), h('th', null, 'Profile'))),
    h('tbody', null, f.diffusion.map(d => h('tr', null, h('td', null, resolveId(d.example) && d.species === d.example ? lnk(d.example) : [h('span', { class: 'mono' }, d.species), d.example !== d.species && resolveId(d.example) ? [' e.g. ', lnk(d.example)] : null]),
      h('td', null, d.immobile ? h('span', { class: 'muted' }, 'immobile — ' + d.immobile.split(': ').slice(1).join(': ') || d.immobile)
        : d.constants.map((c, i) => h('div', null, c.name.startsWith('RNA length') ? c.name : h('a', { class: 'link mono', style: 'font-size:12px', onclick: () => go('constants', c.name) }, c.name), ' ',
          c.min === c.max ? fmt(c.min) : fmt(c.min) + '–' + fmt(c.max), ' m²/s'))),
      h('td', null, (d.profiles || []).map((p, i) => [i ? ', ' : '', lnk(p)]), d.immobile ? h('span', { class: 'faint' }, 'default') : null)))),
    h('div', { class: 'faint', style: 'font-size:12px;margin-top:6px' }, 'Each species moves between lattice regions with the listed coefficients; the full region-by-region table is on its profile. ', RD.immobile_note || '')));
}
function famRateControl(f) {
  const inst0 = f.instances[0], rc = RD.rate_constants[inst0.rate_name];
  if (rc) {
    const cur = (PERT.rdme_rate_constants || []).find(x => x.constant === inst0.rate_name);
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? cur.scale : 1, style: 'width:80px' });
    return h('div', { class: 'actions' }, h('span', { class: 'inline-form' }, 'Rate constant ', h('a', { class: 'link', onclick: () => go('constants', inst0.rate_name) }, inst0.rate_name), ' × ', inp,
      h('button', { class: 'btn sm', onclick: () => setListEntry('rdme_rate_constants', x => x.constant === inst0.rate_name, +inp.value === 1 ? null : { constant: inst0.rate_name, scale: +inp.value }) }, cur ? 'update' : 'set'),
      cur ? h('button', { class: 'btn sm', onclick: () => setListEntry('rdme_rate_constants', x => x.constant === inst0.rate_name, null) }, 'reset') : null),
      h('span', { class: 'faint', style: 'font-size:12px' }, 'shared by every reaction that uses this constant'));
  }
  const key = FAM_RATE_KEY[f.id];
  if (key) {
    const cur = (PERT.gene_rate_scales || []).find(x => x.gene === '*' && x.rate === key);
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? cur.scale : 1, style: 'width:80px' });
    return h('div', { class: 'actions' }, h('span', { class: 'inline-form' }, 'Scale the ' + key.replace('_', ' ') + ' rate of every gene × ', inp,
      h('button', { class: 'btn sm', onclick: () => setListEntry('gene_rate_scales', x => x.gene === '*' && x.rate === key, +inp.value === 1 ? null : { gene: '*', rate: key, scale: +inp.value }) }, cur ? 'update' : 'set'),
      cur ? h('button', { class: 'btn sm', onclick: () => setListEntry('gene_rate_scales', x => x.gene === '*' && x.rate === key, null) }, 'reset') : null),
      h('span', { class: 'faint', style: 'font-size:12px' }, 'single genes: use the table on the gene page'));
  }
  return null;
}
function pageFamily(id) {
  const f = FAM[id], box = h('div');
  if (!f) return h('div', null, 'unknown reaction family ' + id);
  const single = f.count === 1;
  app(box, h('div', { class: 'crumb' }, `${f.layer} reaction${single ? '' : ' family'} · ${f.regions.join(', ')}`),
    h('div', { class: 'title' }, h('h2', null, f.name), h('span', { class: 'mono' }, f.id), badge(f.layer), !single ? badge(f.count + (f.per_gene ? ' genes' : ' reactions')) : null),
    h('div', { class: 'eq' }, single ? instEq(f, f.instances[0]) : [f.template_subs.join(' + '), ' → ', f.template_prods.join(' + ')]),
    h('p', { class: 'lede' }, f.description, f.per_gene ? ' <n> stands for the locus number: one such reaction exists for each gene.' : ''),
    f.process && PR[f.process] ? h('p', { class: 'muted' }, 'Part of the process ', lnk(f.process, PR[f.process].name), '.') : null,
    famRateControl(f));
  app(box, h('h3', null, 'Where it runs'), h('div', { class: 'card' }, 'Lattice regions: ', f.regions.map(r => badge(r)), f.regions.length > 1 ? h('span', { class: 'faint' }, '  (the same reaction is defined in each; a particle reacts in whichever region it is in)') : null));
  const inst0 = f.instances[0], rc = RD.rate_constants[inst0.rate_name];
  app(box, h('h3', null, 'Rate'));
  if (f.rate_tex) app(box, texBlock(f.rate_tex), legendTable((f.legend || []).map(([sym, mean, link]) => [sym, mean, link, null])), h('p', { class: 'faint', style: 'font-size:12px' }, f.rate_note || ''));
  const pgc = PGC_OF[inst0.rate_name];
  if (rc) app(box, h('div', { class: 'card' }, 'Rate constant ', rcLink(inst0.rate_name), ` = ${fmt(rc.value)} ${rc.unit}`, rc.note ? h('span', { class: 'muted' }, ' — ' + rc.note) : null));
  else if (pgc) app(box, h('div', { class: 'card' }, 'Rate constant: one per gene, named ', h('a', { class: 'link mono', style: 'font-size:12px', onclick: () => go('constants', pgc[0]) }, pgc[0]), ' (e.g. ', rcLink(inst0.rate_name), ')'));
  else if (single) app(box, h('div', { class: 'card' }, 'Rate constant ', rcLink(inst0.rate_name), ` = ${fmt(inst0.rate)} ${rateUnit(inst0.order)}`));
  else app(box, h('div', { class: 'tiles' }, tile('Lowest', fmt(f.rate_min) + ' ' + rateUnit(f.order), ''), tile('Median', fmt(f.rate_median) + ' ' + rateUnit(f.order), ''), tile('Highest', fmt(f.rate_max) + ' ' + rateUnit(f.order), '')));
  if (!single) {
    const insts = f.instances.slice().sort((a, b) => b.rate - a.rate);
    const row = i => h('tr', null, h('td', null, i.gene ? lnk(i.gene, i.gene.replace('JCVISYN3A_', '')) : Object.entries(i.vars).map(([k, v]) => `${k}=${v}`).join(' ')), h('td', null, i.gene && GENES[i.gene] ? (GENES[i.gene].gene_name || '') : ''),
      h('td', { class: 'muted' }, instEq(f, i)), h('td', null, f.layer === 'CME' ? h('span', { class: 'faint' }, 'computed') : rcLink(i.rate_name)), h('td', { class: 'num' }, fmt(i.rate)));
    const head = h('thead', null, h('tr', null, h('th', null, f.per_gene ? 'Gene' : 'Index'), h('th', null, 'Name'), h('th', null, 'Reaction'), h('th', null, 'Rate constant'), h('th', null, 'Rate (' + rateUnit(f.order) + ')')));
    app(box, h('h3', null, `All ${f.count} instances`), h('div', { class: 'card' }, h('table', null, head, h('tbody', null, insts.slice(0, 25).map(row))),
      insts.length > 25 ? h('details', null, h('summary', null, `${insts.length - 25} more`), h('table', null, h('tbody', null, insts.slice(25).map(row)))) : null));
  }
  if (f.layer === 'RDME') app(box, h('h3', null, 'Species and diffusion'), diffusionTable(f));
  else app(box, h('h3', null, 'Diffusion'), h('div', { class: 'card muted' }, 'The global CME is well mixed: its species have no position, so no diffusion applies.'));
  if (S.code) app(box, devBox('Built by', ({ rdme: 'processes/Rxns_RDME.py', cme: 'processes/Rxns_CME.py' }[f.id.split(':')[0]]) ? h('div', { class: 'src' }, { rdme: 'processes/Rxns_RDME.py / ImportInitialConditions.py (recorded from the builders)', cme: 'processes/Rxns_CME.py (recorded from the builders)' }[f.id.split(':')[0]]) : null));
  return box;
}
function pageConstants(which) {
  const box = h('div');
  const fc = FILT['constants.cat'] || '', fp = FILT['constants.proc'] || '', fs = FILT['constants.state'] || '';
  const byKey = Object.fromEntries(constCatalogue().map(c => [c.key, c]));
  const keep = k => { const c = byKey[k]; if (!c) return true; if (which && which !== 'all' && k === which) return true;
    return (!fp || c.procs.includes(fp)) && (!fs || (fs === 'used' ? c.used : fs === 'unused' ? !c.used : constEdited(c))); };
  const show = cat => !fc || fc === cat || (which && which !== 'all' && byKey[which] && byKey[which].cat === cat);
  const filtered = [fc && CONST_CATS.find(x => x[0] === fc)[1], fp && PR[fp] && 'process: ' + PR[fp].name, fs && { used: 'used', unused: 'unused', edited: 'edited in this perturbation' }[fs]].filter(Boolean);
  const rcRows = Object.entries(RD.rate_constants).filter(([k]) => keep(k)), dcRows = Object.entries(RD.diffusion_constants).filter(([k]) => k !== 'zero' && keep(k)),
    gRows = Object.entries(RD.gip_constants).filter(([k]) => keep(k));
  const scaleCtl = (key, name, curPred, mk) => {
    const cur = (PERT[key] || []).find(curPred);
    const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? (cur.scale ?? cur.value) : (mk === 'value' ? (RD.gip_constants[name] || {}).value : 1), style: 'width:80px' });
    return h('span', { class: 'inline-form' }, mk === 'value' ? '= ' : '× ', inp,
      h('button', { class: 'btn sm', onclick: () => { const v = +inp.value; const noop = mk === 'value' ? v === (RD.gip_constants[name] || {}).value : v === 1;
        setListEntry(key, curPred, noop ? null : (mk === 'value' ? { constant: name, value: v } : { constant: name, scale: v })); } }, cur ? 'update' : 'set'),
      cur ? h('button', { class: 'btn sm', onclick: () => setListEntry(key, curPred, null) }, 'reset') : null);
  };
  const hl = k => which === k ? 'background:var(--accent-soft)' : null;
  app(box, h('div', { class: 'title' }, h('h2', null, 'Constants'), filtered.length ? badge('filtered: ' + filtered.join(' · '), 'b-reduced') : null,
    filtered.length ? h('button', { class: 'btn sm', onclick: () => { delete FILT['constants.cat']; delete FILT['constants.proc']; delete FILT['constants.state']; store.set('filt', FILT); renderFilters(); renderList(); renderDetail(); } }, 'clear filters') : null), h('p', { class: 'lede' }, 'Numbers written in the Python rather than in a spreadsheet: jLM rate constants shared by whole reaction families, diffusion coefficients, and the constants inside the GIP rate formulas. Each can be scaled or set from a perturbation file.'));
  const rcTable = rows => h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Name'), h('th', null, 'Value'), h('th', null, 'Unit'), h('th', null, 'Used by'), h('th', null, 'Note'), h('th', null, 'Perturb'))),
    h('tbody', null, rows.map(([k, v]) => h('tr', { style: hl(k), 'data-k': k }, h('td', { class: 'mono' }, k), h('td', { class: 'num' }, fmt(v.value)), h('td', null, v.unit),
      h('td', null, v.used_by.length ? v.used_by.map((f, i) => [i ? ', ' : '', lnk(f, famName(f))]) : h('span', { class: 'faint' }, 'unused')), h('td', { class: 'muted' }, v.note || ''),
      h('td', null, v.used_by.length ? scaleCtl('rdme_rate_constants', k, x => x.constant === k, 'scale') : null))))));
  const section = (cat, title, rows, table) => { if (show(cat) && rows.length) app(box, h('h3', null, title + ` (${rows.length})`), table(rows)); };
  section('shared', 'Shared rate constants (Rxns_RDME, ImportInitialConditions)', rcRows.filter(([, v]) => v.group !== 'assembly'), rcTable);
  section('assembly', 'Ribosome assembly rate constants (addRibosomeBiogenesis)', rcRows.filter(([, v]) => v.group === 'assembly'), rcTable);
  // per-gene constants: one row per template, every gene's value underneath
  const openT = PGC_OF[which] ? PGC_OF[which][0] : (RD.per_gene_constants[which] ? which : null);
  const pgRows = Object.entries(RD.per_gene_constants).filter(([t]) => keep(t));
  if (show('per_gene') && pgRows.length) app(box, h('h3', null, `Per-gene rate constants (one per gene, computed by GIP_rates) (${pgRows.length})`), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Name'), h('th', null, 'Rate'), h('th', null, 'Genes'), h('th', null, 'Range'), h('th', null, 'Formula'), h('th', null, 'Perturb'))),
    h('tbody', null, pgRows.map(([t, pg]) => {
      const cur = (PERT.gene_rate_scales || []).find(x => x.gene === '*' && x.rate === pg.gene_rate);
      const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? cur.scale : 1, style: 'width:70px' });
      const insts = pg.instances.slice().sort((x, y) => x.name.localeCompare(y.name));
      return h('tr', { style: hl(t) || (openT === t ? 'background:var(--accent-soft)' : null), 'data-k': t }, h('td', { class: 'mono' }, t), h('td', null, pg.what, h('div', null, pg.used_by.map((f, i) => [i ? ', ' : '', lnk(f, famName(f))]))),
        h('td', { class: 'num' }, pg.count), h('td', { class: 'num' }, `${fmt(pg.min)}–${fmt(pg.max)} ${pg.unit}`, h('div', { class: 'faint' }, 'median ' + fmt(pg.median))), h('td', { class: 'muted' }, pg.formula),
        h('td', null, h('span', { class: 'inline-form' }, 'all × ', inp, h('button', { class: 'btn sm', onclick: () => setListEntry('gene_rate_scales', x => x.gene === '*' && x.rate === pg.gene_rate, +inp.value === 1 ? null : { gene: '*', rate: pg.gene_rate, scale: +inp.value }) }, cur ? 'update' : 'set')),
          h('details', PGC_OF[which] && openT === t ? { open: true } : null, h('summary', null, 'each gene'), h('table', null, h('tbody', null, insts.map(i => h('tr', { style: which === i.name ? 'background:var(--accent-soft)' : null, 'data-k': i.name },
            h('td', { class: 'mono', style: 'font-size:12px' }, i.name), h('td', null, lnk(i.gene, (GENES[i.gene] || {}).gene_name || i.gene.replace('JCVISYN3A_', ''))), h('td', { class: 'num' }, fmt(i.value)))))))));
    })))));
  if (show('diffusion') && (dcRows.length || keep('rna_diff'))) app(box, h('h3', null, 'Diffusion constants (Diffusion.py)'), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Name'), h('th', null, 'Value (m²/s)'), h('th', null, 'Applies to'), h('th', null, 'Perturb'))),
    h('tbody', null, dcRows.map(([k, v]) => h('tr', { style: hl(k), 'data-k': k }, h('td', { class: 'mono' }, k), h('td', { class: 'num' }, fmt(v.value)),
      h('td', { class: 'muted' }, { diffPtn: 'proteins in cytoplasm', diffPtnDna: 'proteins entering / inside the DNA region', diffMemPtn: 'membrane proteins in the membrane', diffDeg: 'degradosome and bound mRNA', RNAP_diff: 'RNA polymerase', ribo_diff: 'ribosomes', tribo_diff: 'translating ribosomes (RB_ states)' }[k] || ''),
      h('td', null, scaleCtl('diffusion_scales', k, x => x.constant === k, 'scale')))),
      !keep('rna_diff') ? null : h('tr', { style: hl('rna_diff'), 'data-k': 'rna_diff' }, h('td', { class: 'mono' }, 'rna_diff'), h('td', { class: 'num' }, 'per RNA'), h('td', { class: 'muted' }, `every RNA and ribosome assembly intermediate (${RD.rna_diff_count || '?'} species, ${fmt((RD.rna_diff_range || [0])[0])}–${fmt((RD.rna_diff_range || [0, 0])[1])} m²/s): Stokes–Einstein from its length (Diffusion.RNA_diff_coeff); halved inside the DNA region`), h('td', null, scaleCtl('diffusion_scales', 'rna_diff', x => x.constant === 'rna_diff', 'scale')))))));
  if (show('gip') && gRows.length) app(box, h('h3', null, `GIP rate-formula constants (utility/GIP_rates.py) (${gRows.length})`), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Name'), h('th', null, 'Value'), h('th', null, 'Unit'), h('th', null, 'Note'), h('th', null, 'Perturb'))),
    h('tbody', null, gRows.map(([k, v]) => h('tr', { style: hl(k), 'data-k': k }, h('td', { class: 'mono' }, k), h('td', { class: 'num' }, fmt(v.value)), h('td', null, v.unit || ''), h('td', { class: 'muted' }, v.note || ''),
      h('td', null, GIP_EDITABLE.includes(k) ? scaleCtl('gip_constants', k, x => x.constant === k, 'value') : h('span', { class: 'faint' }, 'derived'))))))));
  const profRows = Object.values(RD.diffusion_profiles).filter(p => keep(p.id));
  if (show('profile') && profRows.length) app(box, h('h3', null, `Diffusion profiles (${profRows.length})`), h('p', { class: 'muted' }, (RD.immobile_note || '') + ' ', Object.entries(RD.immobile_patterns || {}).map(([k, v]) => h('div', null, h('span', { class: 'mono' }, k), ': ' + v))),
    h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Profile'), h('th', null, 'Species'), h('th', null, 'Kinds'), h('th', null, 'Example'), h('th', null, 'Moves'))),
    h('tbody', null, profRows.map(p => h('tr', { style: hl(p.id), 'data-k': p.id }, h('td', { class: 'mono' }, p.id), h('td', { class: 'num' }, p.members), h('td', { class: 'muted' }, p.kinds.join(', ')), h('td', null, p.example === '*' ? 'default (all species, zero)' : lnk(p.example)),
      h('td', null, h('details', null, h('summary', null, p.transitions.length + ' transitions'), h('table', null, h('tbody', null, p.transitions.map(([a, b, d]) => h('tr', null, h('td', null, a), h('td', null, '→ ' + b), h('td', { class: 'mono' }, d)))))))))))));
  if (box.querySelectorAll('h3').length === 0) app(box, h('div', { class: 'card muted' }, 'No constants match these filters.'));
  if (which && which !== 'all') setTimeout(() => { const row = box.querySelector(`tr[data-k="${CSS.escape(which)}"]`) || [...box.querySelectorAll('tr')].find(r => r.style.background); if (row) row.scrollIntoView({ block: row.getBoundingClientRect().height > innerHeight / 2 ? 'start' : 'center' }); }, 0);
  return box;
}
function geneReactionsTable(g, keep) {
  const rows = [];
  for (const f of Object.values(FAM)) {
    if (!f.per_gene || (keep && !keep(f))) continue;
    for (const i of f.instances) if (i.gene === g.locus) rows.push([f, i]);
  }
  if (!rows.length) return null;
  const key = f => FAM_RATE_KEY[f.id];
  return h('div', null, h('h3', null, `Reactions for this gene (${rows.length})`), h('div', { class: 'card' }, h('table', null,
    h('thead', null, h('tr', null, h('th', null, 'Type'), h('th', null, 'Reaction'), h('th', null, 'Layer'), h('th', null, 'Rate'), h('th', null, 'Perturb'))),
    h('tbody', null, rows.map(([f, i]) => {
      const k = key(f), cur = k && (PERT.gene_rate_scales || []).find(x => x.gene === g.locus && x.rate === k);
      const inp = h('input', { type: 'number', min: 0, step: 'any', value: cur ? cur.scale : 1, style: 'width:70px' });
      return h('tr', null, h('td', null, lnk(f.id, f.name)), h('td', { class: 'muted' }, instEq(f, i)), h('td', null, f.layer),
        h('td', { class: 'num' }, fmt(i.rate) + ' ' + rateUnit(i.order), f.layer === 'RDME' ? h('div', { style: 'font-size:11px' }, RD.rate_constants[i.rate_name] ? 'shared ' : '', rcLink(i.rate_name)) : null),
        h('td', null, k ? h('span', { class: 'inline-form' }, '× ', inp, h('button', { class: 'btn sm', onclick: () => setListEntry('gene_rate_scales', x => x.gene === g.locus && x.rate === k, +inp.value === 1 ? null : { gene: g.locus, rate: k, scale: +inp.value }) }, cur ? 'update' : 'set'),
          cur ? h('button', { class: 'btn sm', onclick: () => setListEntry('gene_rate_scales', x => x.gene === g.locus && x.rate === k, null) }, 'reset') : null) : h('span', { class: 'faint' }, 'via the constant')));
    })))));
}
function diffusionBlock(id) {
  const d = RD.species_diffusion[id];
  if (!d) return null;
  const p = RD.diffusion_profiles[d.profile];
  const val = dc => d.own[dc] !== undefined ? d.own[dc] : (RD.diffusion_constants[dc] || {}).value;
  return h('div', null, h('h3', null, 'Diffusion'), h('div', { class: 'card' }, h('div', { class: 'muted', style: 'margin-bottom:6px' }, 'Profile ', lnk(p.id), `, shared with ${p.members - 1} other species (${p.kinds.join(', ')}).`, Object.keys(d.own).length ? ' Own coefficients from the RNA length: ' + Object.entries(d.own).map(([k, v]) => `${k.endsWith('Dna') ? 'in DNA' : 'free'} ${fmt(v)} m²/s`).join(', ') : ''),
    h('details', null, h('summary', null, p.transitions.length + ' region transitions'), h('table', null, h('thead', null, h('tr', null, h('th', null, 'From'), h('th', null, 'To'), h('th', null, 'Constant'), h('th', null, 'm²/s'))),
      h('tbody', null, p.transitions.map(([a, b, dc]) => { const real = dc.startsWith('rna_diff') ? Object.keys(d.own).find(k => k.endsWith('Dna') === dc.endsWith('Dna')) || dc : dc; return h('tr', null, h('td', null, a), h('td', null, b), h('td', { class: 'mono' }, real), h('td', { class: 'num' }, fmt(val(real)))); }))))));
}
function pageInputs(key) {
  const box = h('div');
  app(box, h('div', { class: 'title' }, h('h2', null, 'Input files')), h('p', { class: 'lede' }, 'Every file and sheet in input_data/, and where the code reads it. Unused sheets are safe to edit but change nothing.'));
  const files = [...new Set(S.inputs.map(r => r.file))];
  for (const f of files) {
    const rows = S.inputs.filter(r => r.file === f);
    app(box, h('h3', null, f, ' ', S.meta.input_sha256 && S.meta.input_sha256[f] ? h('span', { class: 'faint mono', style: 'text-transform:none' }, S.meta.input_sha256[f].slice(0, 12)) : null),
      h('div', { class: 'card' }, h('table', null, h('thead', null, h('tr', null, h('th', null, 'Sheet'), h('th', null, 'Rows'), h('th', null, 'Read by the code'), h('th', null, 'Functions'))), h('tbody', null, rows.map(r => h('tr', { style: key === r.file + '|' + (r.sheet || '') ? 'background:var(--accent-soft)' : null },
        h('td', { style: 'width:30%' }, r.sheet || h('span', { class: 'faint' }, 'whole file')), h('td', { class: 'num', style: 'width:70px' }, r.rows !== undefined ? r.rows + ' rows' : ''),
        h('td', { style: 'width:80px' }, r.used ? badge('read', 'b-ok') : badge('unused', 'b-inactive')),
        h('td', null, h('div', { class: 'dev-only' }, refList(r.reads, 6)), h('span', { class: 'muted' }, r.reads.filter(x => !x.c).length ? r.reads.filter(x => !x.c).map(x => x.fn).filter((v, i, a) => a.indexOf(v) === i).join(', ') : ''))))))));
  }
  return box;
}

/* ------------------------------------------------------------------ cart */
let SERVED = false, SERVER_MSG = null;
function renderCart() {
  const c = $('#cart'); c.innerHTML = '';
  const nameIn = h('input', { value: PERT.name, placeholder: 'name (file name)', 'aria-label': 'perturbation name', onchange: e => { PERT.name = e.target.value.trim(); store.set('pert', PERT); renderCart(); } });
  const descIn = h('textarea', { placeholder: 'description (optional)', 'aria-label': 'description', onchange: e => { PERT.description = e.target.value; store.set('pert', PERT); renderCart(); } }, PERT.description || '');
  c.appendChild(h('div', { class: 'cart-h' }, h('h2', null, 'Perturbation'), nameIn, descIn));
  const b = h('div', { class: 'cart-b' }); c.appendChild(b);
  const items = [];
  const rm = fn => h('button', { title: 'remove', onclick: fn }, '×');
  for (const k of PERT.knockouts) {
    const g = GENES[k.gene];
    const sel = h('select', { onchange: e => setKO(k.gene, e.target.value) }, opts(Object.keys(MODES).map(m => [m, m.replace('_', ' ')]))); sel.value = k.mode;
    items.push(h('div', { class: 'edit' }, h('div', null, badge('KO', 'b-ko'), ' ', lnk(k.gene, `${k.gene.replace('JCVISYN3A_', '')} ${g.gene_name || g.symbol || ''}`), ' ', sel, ' ', infoTip(k.mode)), rm(() => setKO(k.gene, null)), h('div', { class: 'w' }, g.product_proteomics || g.product)));
  }
  for (const k of PERT.knockdowns) { const g = GENES[k.gene]; items.push(h('div', { class: 'edit' }, h('div', null, badge('KD', 'b-reduced'), ' ', lnk(k.gene, `${k.gene.replace('JCVISYN3A_', '')} ${g.gene_name || ''}`), ` promoter × ${k.promoter_scale}`), rm(() => setKD(k.gene, null)), h('div', { class: 'w' }, `${fmt(g.promoter)} → ${fmt(g.promoter * k.promoter_scale)}`))); }
  const mapLabels = { initial_protein_counts: ['initial', 'initial_count', ''], initial_mrna_means: ['initial mRNA', 'initial_mean', ''], initial_metabolites_mM: ['initial', 'initial_mM', ' mM'], medium_mM: ['medium', 'medium_mM', ' mM'] };
  for (const [key, [word, field, unit]] of Object.entries(mapLabels)) for (const [id, v] of Object.entries(PERT[key]))
    items.push(h('div', { class: 'edit' }, h('div', null, `${word} `, lnk(id), ` = ${fmt(v)}${unit}`), rm(() => setMap(key, id, null)), h('div', { class: 'w' }, `was ${fmt(SP[id] ? SP[id][field] : null)}${unit}`)));
  for (const p of PERT.reaction_parameters) { const par = RX[p.reaction] && RX[p.reaction].params.find(x => x.key === p.parameter); items.push(h('div', { class: 'edit' }, h('div', null, lnk(p.reaction), `.${p.parameter} ${p.scale !== undefined ? '× ' + p.scale : '= ' + fmt(p.value)}`), rm(() => setParam(p.reaction, p.parameter, 'scale', null)), h('div', { class: 'w' }, par ? `${fmt(par.value)} → ${fmt(p.scale !== undefined ? par.value * p.scale : p.value)}` : ''))); }
  for (const x of PERT.rdme_rate_constants || []) items.push(h('div', { class: 'edit' }, h('div', null, 'rate constant ', h('a', { class: 'link', onclick: () => go('constants', x.constant) }, x.constant), x.scale !== undefined ? ` × ${x.scale}` : ` = ${fmt(x.value)}`), rm(() => setListEntry('rdme_rate_constants', y => y.constant === x.constant, null)), h('div', { class: 'w' }, RD.rate_constants[x.constant] ? `was ${fmt(RD.rate_constants[x.constant].value)} ${RD.rate_constants[x.constant].unit}` : '')));
  for (const x of PERT.gene_rate_scales || []) items.push(h('div', { class: 'edit' }, h('div', null, x.rate.replace('_', ' ') + ' rate of ', x.gene === '*' ? 'every gene' : lnk(x.gene, x.gene.replace('JCVISYN3A_', '')), ` × ${x.scale}`), rm(() => setListEntry('gene_rate_scales', y => y.gene === x.gene && y.rate === x.rate, null)), h('div', { class: 'w' }, 'computed per gene, then scaled')));
  for (const x of PERT.diffusion_scales || []) items.push(h('div', { class: 'edit' }, h('div', null, 'diffusion ', h('a', { class: 'link', onclick: () => go('constants', x.constant) }, x.constant), ` × ${x.scale}`), rm(() => setListEntry('diffusion_scales', y => y.constant === x.constant, null)), h('div', { class: 'w' }, RD.diffusion_constants[x.constant] ? `was ${fmt(RD.diffusion_constants[x.constant].value)} m²/s` : 'every RNA')));
  for (const x of PERT.gip_constants || []) items.push(h('div', { class: 'edit' }, h('div', null, 'GIP constant ', h('a', { class: 'link', onclick: () => go('constants', x.constant) }, x.constant), ` = ${fmt(x.value)}`), rm(() => setListEntry('gip_constants', y => y.constant === x.constant, null)), h('div', { class: 'w' }, RD.gip_constants[x.constant] ? `was ${fmt(RD.gip_constants[x.constant].value)} ${RD.gip_constants[x.constant].unit || ''}` : '')));
  for (const r of PERT.disabled_reactions) items.push(h('div', { class: 'edit' }, h('div', null, badge('off', 'b-blocked'), ' ', lnk(r), ' disabled'), rm(() => toggleDisabled(r)), h('div', { class: 'w' }, RX[r] ? RX[r].subsystem : '')));
  if (PERT.knockouts.length) b.appendChild(h('div', { class: 'faint', style: 'font-size:12px;margin:4px 0 2px' }, 'Knockout modes: ', modeTips()));
  b.appendChild(items.length ? h('div', null, items) : h('div', { class: 'empty' }, 'Nothing yet. Open a gene and press Knock out, or set a concentration or rate constant.'));

  if (items.length) {
    const imp = pertImpact(), cnt = o => Object.values(o).reduce((a, v) => (a[v.status] = (a[v.status] || 0) + 1, a), {});
    const pr = cnt(imp.processes), rx = cnt(imp.reactions), me = cnt(imp.metabolites);
    b.append(h('h3', null, 'Static impact'), h('div', { class: 'sumline' },
      pr.blocked ? badge(`${pr.blocked} processes blocked`, 'b-blocked') : null, pr.reduced ? badge(`${pr.reduced} processes reduced`, 'b-reduced') : null,
      rx.blocked ? badge(`${rx.blocked} reactions blocked`, 'b-blocked') : null, rx.reduced ? badge(`${rx.reduced} reactions reduced`, 'b-reduced') : null,
      me.blocked ? badge(`${me.blocked} metabolites without a source`, 'b-blocked') : null,
      !Object.keys(imp.processes).length && !Object.keys(imp.reactions).length ? h('span', { class: 'faint' }, 'no reaction or process loses a part') : null),
      h('details', null, h('summary', null, 'details'), impactBlock(imp)));
  }
  const errs = validateLocal(pertForExport());
  const yaml = dumpYaml(pertForExport());
  b.append(h('h3', null, 'YAML'), h('div', { class: 'yaml' }, yaml));
  const fileIn = h('input', { type: 'file', accept: '.yaml,.yml,.json', style: 'display:none', onchange: e => { const f = e.target.files[0]; if (f) f.text().then(t => importText(t, f.name)); e.target.value = ''; } });
  b.append(h('div', { class: 'row-btns' },
    h('button', { class: 'btn primary', disabled: !!errs.length || !items.length, onclick: () => download(yaml) }, DOWNLOADS ? 'Save .json' : 'Download .yaml'),
    h('button', { class: 'btn', onclick: () => copy(yaml, 'YAML') }, 'Copy'),
    SERVED ? h('button', { class: 'btn', disabled: !!errs.length || !items.length, onclick: () => serverSave(false) }, 'Save to repo') : null,
    SERVED ? h('button', { class: 'btn', disabled: !items.length, onclick: serverValidate }, 'Validate') : null,
    h('button', { class: 'btn', onclick: () => fileIn.click() }, 'Open…'), fileIn,
    h('button', { class: 'btn', id: 'clear-pert', onclick: e => {
      const b = e.currentTarget;
      if (items.length && b.dataset.armed !== '1') { b.dataset.armed = '1'; b.textContent = 'Click again to clear'; b.classList.add('danger'); setTimeout(() => { if (b.isConnected) { b.dataset.armed = ''; b.textContent = 'Clear'; b.classList.remove('danger'); } }, 4000); return; }
      PERT = emptyPert(); savePert(); SERVER_MSG = null; } }, 'Clear')));
  if (SERVED && items.length && PERT.knockouts.length > 1) b.appendChild(h('div', { class: 'row-btns' }, h('button', { class: 'btn', onclick: saveEachKO }, `Save ${PERT.knockouts.length} single-knockout files`)));
  const msgs = h('div', { class: 'msgs' }, errs.map(e => h('div', { class: 'e' }, e)));
  if (SERVER_MSG) msgs.append(...SERVER_MSG.map(([k, t]) => h('div', { class: k }, t)));
  b.appendChild(msgs);
  b.appendChild(h('p', { class: 'faint', style: 'font-size:12px' }, SERVED ? 'Served by python -m modelspec serve: Save writes perturbations/<name>.yaml after the Python validator accepts it.'
    : DOWNLOADS ? 'Save .json writes the perturbation as JSON, which the simulator reads as is: python Whole_Cell_Minimal_Cell.py ... -p <name>.json. Check it first with python3 -m modelspec check <name>.json.'
    : 'Opened as a file: download the YAML and check it with python -m modelspec check. Run python -m modelspec serve to save straight into the repo.'));
  if (SERVED) serverList(b);
  linkifyTree(b);
}
// In a claude.ai artifact the page cannot start downloads itself; the `downloads` capability asks the viewer instead.
// Its allowlist has no .yaml, so the shared page saves .json (JSON is valid YAML: `-p file.json` works as is).
let DOWNLOADS = null;
if (window.claude && typeof window.claude.use === 'function') window.claude.use('downloads').then(d => { DOWNLOADS = d; renderCart(); }).catch(() => {});
async function download(text) {
  if (DOWNLOADS) {
    const name = (PERT.name || 'perturbation') + '.json';
    try { await DOWNLOADS.save({ filename: name, data: JSON.stringify(pertForExport(), null, 1) + '\n' }); toast('Saved ' + name); }
    catch (e) { toast(e && e.code === 'declined' ? 'Save cancelled' : 'Could not save: ' + ((e && e.message) || e)); }
    return;
  }
  const a = h('a', { href: URL.createObjectURL(new Blob([text], { type: 'text/yaml' })), download: (PERT.name || 'perturbation') + '.yaml' });
  document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}
function importText(text, name) {
  try {
    const p = name && name.endsWith('.json') ? JSON.parse(text) : loadYaml(text);
    const e = validateLocal(Object.assign(emptyPert(), p));
    PERT = Object.assign(emptyPert(), p); delete PERT.created; delete PERT.model;
    for (const k of ['knockouts', 'knockdowns', 'reaction_parameters', 'disabled_reactions', 'rdme_rate_constants', 'gene_rate_scales', 'diffusion_scales', 'gip_constants']) PERT[k] = PERT[k] || [];
    for (const k of ['initial_protein_counts', 'initial_mrna_means', 'initial_metabolites_mM', 'medium_mM']) PERT[k] = PERT[k] || {};
    SERVER_MSG = [[e.length ? 'w' : 'o', `Opened ${name || 'file'}` + (p.model && p.model.fingerprint && p.model.fingerprint !== S.meta.fingerprint ? ' — built against other inputs (' + p.model.fingerprint + ')' : '')]];
    savePert();
  } catch (err) { SERVER_MSG = [['e', 'Could not read ' + (name || 'file') + ': ' + err.message]]; renderCart(); }
}
async function api(path, body) {
  const r = await fetch(path, body ? { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) } : undefined);
  return [r.status, await r.json()];
}
async function serverValidate() {
  const [code, res] = await api('/api/validate', { perturbation: pertForExport() });
  SERVER_MSG = [...(res.errors || []).map(e => ['e', e]), ...(res.warnings || []).map(w => ['w', w])];
  if (code === 200) SERVER_MSG.push(['o', `valid: ${res.resolution.edits.length} edits`], ...res.resolution.edits.map(e => ['', `${e.target}: ${fmt(e.from)} → ${fmt(e.to)}  [${e.where}]${e.note ? ' — ' + e.note : ''}`]));
  renderCart();
}
async function serverSave(overwrite, pert) {
  pert = pert || pertForExport();
  let [code, res] = await api('/api/save', { perturbation: pert, overwrite });
  if (code === 409 && confirm(`${pert.name}.yaml exists in perturbations/. Overwrite?`)) [code, res] = await api('/api/save', { perturbation: pert, overwrite: true });
  SERVER_MSG = code === 200 ? [['o', 'Saved ' + res.saved], ...(res.warnings || []).map(w => ['w', w])] : (res.errors || ['save failed']).map(e => ['e', e]);
  renderCart();
}
async function saveEachKO() {
  const base = pertForExport(), saved = [];
  for (const k of PERT.knockouts) {
    const p = Object.assign({}, base, { name: 'ko_' + k.gene.replace('JCVISYN3A_', ''), knockouts: [k], knockdowns: [], initial_protein_counts: {}, initial_mrna_means: {},
      initial_metabolites_mM: {}, medium_mM: {}, reaction_parameters: [], disabled_reactions: [], description: `single knockout of ${k.gene} (${k.mode})` });
    const [code, res] = await api('/api/save', { perturbation: p, overwrite: true });
    saved.push(code === 200 ? ['o', 'Saved ' + res.saved] : ['e', p.name + ': ' + (res.errors || []).join('; ')]);
  }
  SERVER_MSG = saved; renderCart();
}
async function serverList(b) {
  try {
    const [, res] = await api('/api/list');
    if (!res.files.length) return;
    b.append(h('h3', null, 'In perturbations/'), h('div', null, res.files.map(f => h('div', null, h('a', { class: 'link', onclick: async () => { const [, r] = await api('/api/file/' + encodeURIComponent(f)); importText(JSON.stringify(r.perturbation), f.replace('.yaml', '.json')); } }, f)))));
  } catch (e) { /* not served */ }
}

/* ------------------------------------------------------------------ resizable panes */
const PANE = { list: { v: '--w-list', def: 320, min: 200, max: 640 }, cart: { v: '--w-cart', def: 360, min: 260, max: 760 } };
function setPane(which, px, save = true) {
  const p = PANE[which], detailMin = 360;
  const other = which === 'list' ? paneWidth('cart') : paneWidth('list');
  const room = window.innerWidth - other - 12 - detailMin;          // keep the middle pane usable
  px = Math.round(Math.max(p.min, Math.min(p.max, room, px)));
  document.documentElement.style.setProperty(p.v, px + 'px');
  const sep = document.getElementById(which === 'list' ? 'split-list' : 'split-cart');
  sep.setAttribute('aria-valuenow', px); sep.setAttribute('aria-valuemin', p.min); sep.setAttribute('aria-valuemax', p.max);
  if (save) store.set('pane.' + which, px);
}
function paneWidth(which) { return parseFloat(getComputedStyle(document.documentElement).getPropertyValue(PANE[which].v)) || PANE[which].def; }
function initSplitter(id, which) {
  const sep = document.getElementById(id);
  // the list grows when its right edge moves right; the cart grows when its left edge moves left
  const sign = which === 'list' ? 1 : -1;
  sep.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    e.preventDefault(); sep.setPointerCapture(e.pointerId);
    const x0 = e.clientX, w0 = paneWidth(which);
    sep.classList.add('dragging'); document.body.classList.add('resizing');
    const move = ev => setPane(which, w0 + sign * (ev.clientX - x0), false);
    const up = ev => { sep.releasePointerCapture(ev.pointerId); sep.removeEventListener('pointermove', move); sep.removeEventListener('pointerup', up); sep.removeEventListener('pointercancel', up);
      sep.classList.remove('dragging'); document.body.classList.remove('resizing'); store.set('pane.' + which, paneWidth(which)); };
    sep.addEventListener('pointermove', move); sep.addEventListener('pointerup', up); sep.addEventListener('pointercancel', up);
  });
  sep.addEventListener('dblclick', () => setPane(which, PANE[which].def));
  sep.addEventListener('keydown', e => {
    const step = e.shiftKey ? 60 : 20;
    if (e.key === 'ArrowLeft') { e.preventDefault(); setPane(which, paneWidth(which) - sign * step); }
    else if (e.key === 'ArrowRight') { e.preventDefault(); setPane(which, paneWidth(which) + sign * step); }
    else if (e.key === 'Home' || e.key === 'Enter') { e.preventDefault(); setPane(which, PANE[which].def); }
  });
}
for (const w of ['list', 'cart']) setPane(w, store.get('pane.' + w, PANE[w].def), false);
initSplitter('split-list', 'list'); initSplitter('split-cart', 'cart');
window.addEventListener('resize', () => { for (const w of ['list', 'cart']) setPane(w, paneWidth(w), false); });

/* ------------------------------------------------------------------ boot */
function applyTheme(t) { if (t) document.documentElement.setAttribute('data-theme', t); else document.documentElement.removeAttribute('data-theme'); }
applyTheme(store.get('theme', null));
$('#theme').onclick = () => {
  const dark = matchMedia('(prefers-color-scheme: dark)').matches, cur = document.documentElement.getAttribute('data-theme') || (dark ? 'dark' : 'light');
  const next = cur === 'dark' ? 'light' : 'dark'; applyTheme(next); store.set('theme', next);
};
const devBoxEl = $('#dev'); devBoxEl.checked = store.get('dev', false); document.body.classList.toggle('dev', devBoxEl.checked);
devBoxEl.onchange = () => { document.body.classList.toggle('dev', devBoxEl.checked); store.set('dev', devBoxEl.checked); };
$('#meta').textContent = `commit ${S.meta.git_commit || '?'} · inputs ${S.meta.fingerprint} · built ${S.meta.built}`;
$('#q').addEventListener('input', e => { Q = e.target.value; renderList(); });
document.addEventListener('keydown', e => { if (e.key === '/' && document.activeElement.tagName !== 'INPUT' && document.activeElement.tagName !== 'TEXTAREA') { e.preventDefault(); $('#q').focus(); } });
if (SEL && !['inputs', 'constants'].includes(SEL.type) && !({ gene: GENES, species: SP, reaction: RX, process: PR, family: FAM }[SEL.type] || {})[SEL.id]) SEL = null;
renderTabs(); renderFilters(); renderList(); renderDetail(); renderCart();
if (location.protocol === 'http:') fetch('/api/list').then(r => { if (r.ok) { SERVED = true; renderCart(); } }).catch(() => {});
if (SELFCHECK.bad.length) console.warn('impact engine differs from Python for', SELFCHECK.bad);
