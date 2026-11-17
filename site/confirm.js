// The allele confirmation page: map the reads to the genome plus the introduced cassette (in a worker, in WebAssembly) and say where the cassette
// joins the genome and whether it is intact. The files stay in this browser; web addresses are fetched by this page.
import { MAX_DOWNLOAD, mb, checkUrl, download as downloadFile, createReferencePanel } from './reference-panel.js?v=__BUILD__';

const $ = (id) => document.getElementById(id);
const show = (id, text) => { $(id).textContent = text; };

// How far a file may unpack when it is loaded whole into memory. `?limit_mb=N` can lower it (the tests use that). The download limit is MAX_DOWNLOAD (reference-panel.js).
const LIMIT_BYTES = Math.min(parseFloat(new URLSearchParams(location.search).get('limit_mb')) || 1024, 1024) * 1e6;

// A gene that the insertion hits links to its FungiDB gene page. The link works for gene IDs that FungiDB uses.
const GENE_LINK = 'https://fungidb.org/gene/';

let activeWorker = null;
let activeAbort = null;
let running = false;
let cancelled = false;
const links = {}; // the object addresses of the download links, to release them

/** A download link for `text`. */
function offer(id, filename, text, type) {
  if (links[id]) URL.revokeObjectURL(links[id]);
  links[id] = URL.createObjectURL(new Blob([text], { type }));
  const a = $(id);
  a.href = links[id];
  a.download = filename;
  a.hidden = false;
}

// ---- The target choice shows its own fields ------------------------------------------------------------------------------------------------
function showTarget() {
  const kind = $('confirm-target-kind').value;
  for (const group of document.querySelectorAll('[data-target-group]')) group.hidden = group.dataset.targetGroup !== kind;
}
$('confirm-target-kind').addEventListener('change', showTarget);
showTarget();

// ---- The construct ----------------------------------------------------------------------------------------------------------------------
/** FASTA text of the construct: the file or the pasted text. A text with no `>` line gets the name "construct". Returns '' when there is none. */
async function constructText() {
  const file = $('confirm-construct-file').files[0];
  let text = file ? await file.text() : $('confirm-construct-text').value;
  text = text.trim();
  if (!text) return '';
  return /^>/m.test(text) ? `${text}\n` : `>construct\n${text}\n`;
}
async function describeConstruct() {
  try {
    const text = await constructText();
    if (!text) { show('confirm-construct-info', ''); return; }
    const records = text.split(/^>/m).slice(1);
    const bases = records.reduce((n, r) => n + r.split('\n').slice(1).join('').replace(/\s/g, '').length, 0);
    show('confirm-construct-info', `${records.length} sequence${records.length === 1 ? '' : 's'}, ${bases.toLocaleString()} bases.`);
  } catch (error) {
    show('confirm-construct-info', `The file cannot be read: ${error.message}`);
  }
}
$('confirm-construct-file').addEventListener('change', describeConstruct);
$('confirm-construct-text').addEventListener('input', describeConstruct);

// ---- Downloads --------------------------------------------------------------------------------------------------------------------------
/** Run one message in a worker of its own. The worker is stopped when the answer comes, or when the run is cancelled. */
function callWorker(message, onProgress) {
  return new Promise((resolve, reject) => {
    const w = new Worker('worker.js?v=__BUILD__', { type: 'module' });
    activeWorker = { w, reject };
    const finish = () => { w.terminate(); if (activeWorker && activeWorker.w === w) activeWorker = null; };
    w.onmessage = (event) => {
      const m = event.data;
      if (m.type === 'status' && message.step === 'confirm') show('confirm-status', m.text);
      if (m.type === 'status' && (message.step === 'breakpoints' || message.step === 'assembly')) show('bp-status', m.text);
      if (m.type === 'progress' && onProgress) onProgress(m);
      if (m.type === 'error') { finish(); reject(new Error(m.message)); }
      if (m.type === 'done') { finish(); resolve(m); }
    };
    w.onerror = (event) => { finish(); reject(new Error(event.message || 'the worker stopped')); };
    w.postMessage({ limit: LIMIT_BYTES, ...message });
  });
}

/** The first reads of a download, kept in a worker (see `fetchReads` in worker.js). */
async function downloadFirstReads(address, label, maxReads, budget) {
  const url = checkUrl(address, label);
  const m = await callWorker({ step: 'fetchreads', url: url.href, maxReads, maxKept: budget.left, maxDownload: budget.left, label },
    (p) => show('confirm-status', `Reading ${label}: ${mb(p.downloaded)} downloaded, ${p.reads.toLocaleString()} reads kept`));
  budget.left -= m.bytes.length;
  const name = `${decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'reads')}`.replace(/\.gz$/, '');
  return new File([m.bytes], `${name}.gz`);
}

// ---- The reference genome: the same panel as the main page (reference-panel.js) -------------------------------------------------------------
let loading = false;
// What the menu or the addresses loaded. A file chosen in an input always wins over it.
const loaded = { reference: null, gff: null, source: { reference: '', gff: '' } };
const refFile = () => $('reference').files[0] || loaded.reference;
const gffFile = () => $('gff').files[0] || loaded.gff;
createReferencePanel({
  $, show, state: loaded,
  canLoad: () => !loading && !running,
  setLoading: (flag) => { loading = flag; },
  onLoaded: () => {},
  texts: {
    noGenes: ' NCBI has no gene models for this genome, so a gene cannot be found by name unless you upload a GFF3 file; a region still works.',
    next: ' Now give the cassette, the target and the reads.',
  },
});

// ---- The run ----------------------------------------------------------------------------------------------------------------------------
const number = (id, fallback) => { const v = Number($(id).value); return Number.isFinite(v) && v >= 0 ? v : fallback; };

/** The target of the form: { target } for a region, { gene } for a gene name, {} for none. Throws an Error that starts with "Type a" for a bad region. */
function targetOfForm() {
  const kind = $('confirm-target-kind').value;
  if (kind === 'region') {
    const m = $('confirm-region').value.trim().match(/^(\S+):([\d,]+)-([\d,]+)$/);
    if (!m) throw new Error('Type a region such as Chr1:50,001-51,500 (the contig, a colon, the start, a dash and the end).');
    const [start, end] = [Number(m[2].replace(/,/g, '')), Number(m[3].replace(/,/g, ''))];
    if (!(start >= 1 && end >= start)) throw new Error('Type a region whose start is 1 or more and not above its end.');
    return { target: { contig: m[1], start, end } };
  }
  if (kind === 'gene') {
    const gene = $('confirm-gene').value.trim();
    if (!gene) throw new Error('Type the name of the gene.');
    return { gene };
  }
  return {};
}

$('confirm-run').addEventListener('click', async () => {
  if (running) return;
  if (loading) { show('confirm-status', 'Wait until the genome has loaded.'); return; }
  const reference = refFile();
  if (!reference) {
    show('confirm-status', 'Choose the reference genome first (section 1): load one from the NCBI list, upload a file, or download it from an address.');
    return;
  }
  let construct;
  try { construct = await constructText(); } catch (error) { show('confirm-status', `Give the construct: the file cannot be read (${error.message}).`); return; }
  if (!construct) { show('confirm-status', 'Give the construct: a FASTA file or the pasted sequence of the cassette (section 2).'); return; }
  let wanted;
  try { wanted = targetOfForm(); } catch (error) { show('confirm-status', error.message); return; }
  const gff = gffFile();
  if (wanted.gene && !gff) {
    show('confirm-status', loaded.reference && !loaded.gff && !$('gff').files[0]
      ? 'Choose a GFF3 file or give a region: this genome has no gene models, and a gene name needs them (section 1).'
      : 'Choose the GFF3 file that holds the gene: a gene name needs it (section 1).');
    return;
  }
  const localReads = [...$('confirm-reads').files];
  const addresses = $('confirm-reads-urls').value.split('\n').map((x) => x.trim()).filter(Boolean);
  if (!localReads.length && !addresses.length) { show('confirm-status', 'Choose the reads: files from your computer or https addresses (section 3).'); return; }
  let urls;
  try {
    urls = addresses.map((a) => checkUrl(a, 'the reads').href); // all checked before any download
  } catch (error) { show('confirm-status', `Error: ${error.message}`); return; }

  running = true;
  cancelled = false;
  $('confirm-run').disabled = true;
  $('confirm-cancel').hidden = false;
  $('confirm-result').hidden = true;
  $('locus-section').hidden = true;
  activeAbort = new AbortController();
  try {
    const budget = { left: MAX_DOWNLOAD };
    const limitThousands = number('confirm-read-limit', 0);
    const maxReads = Math.round(limitThousands * 1000);
    const reads = [...localReads];
    for (const [i, address] of urls.entries()) {
      const label = `reads file ${i + 1} of ${urls.length}`;
      reads.push(maxReads ? await downloadFirstReads(address, label, maxReads, budget)
        : await downloadFile({ url: address }, { label, signal: activeAbort.signal, budget, onProgress: (text) => show('confirm-status', text) }));
    }
    show('confirm-status', 'Starting');
    const m = await callWorker({
      step: 'confirm', reference, construct: new TextEncoder().encode(construct), reads, kind: $('confirm-kind').value, gff: gff || null,
      ...wanted, maxReads, geneLink: GENE_LINK,
      options: { minReads: Math.max(1, Math.round(number('confirm-min-reads', 3))), minDepth: number('confirm-min-depth', 10), tolerance: $('confirm-tolerance').value === '' ? -1 : Math.round(number('confirm-tolerance', 20)) },
    }, (p) => show('confirm-status', `Mapping: ${p.reads.toLocaleString()} reads`));
    render(m);
    renderLocus(m.locus, 'The main insertion site found from the reads.');
    show('confirm-status', `Done: ${m.reads.toLocaleString()} reads in ${m.seconds.toFixed(1)} s.`);
  } catch (error) {
    show('confirm-status', cancelled || error.name === 'AbortError' ? 'Cancelled.' : `Error: ${error.message}`);
  } finally {
    running = false;
    activeAbort = null;
    $('confirm-run').disabled = false;
    $('confirm-cancel').hidden = true;
  }
});

$('confirm-cancel').addEventListener('click', () => {
  cancelled = true;
  if (activeAbort) activeAbort.abort();
  if (activeWorker) { activeWorker.w.terminate(); activeWorker.reject(new DOMException('cancelled', 'AbortError')); activeWorker = null; }
});

// ---- The result -------------------------------------------------------------------------------------------------------------------------
const CHECK_LABELS = {
  integration_found: 'Cassette joined to the genome',
  both_sides_joined: 'Both ends joined',
  single_site: 'Single insertion site',
  at_target: 'Place of insertion',
  target_replaced: 'Target replaced',
  cassette_coverage: 'Cassette coverage',
  cassette_sequence: 'Cassette sequence',
  cassette_structure: 'Cassette structure',
  copy_number: 'Copy number',
  whole_insertion_in_one_read: 'Whole insertion in one read',
};
const CHIPS = { pass: ['✓ Pass', 'pass'], fail: ['✗ Fail', 'fail'], inconclusive: ['? Unclear', 'unclear'], not_applicable: ['– n/a', 'na'] };

function el(tag, text, className) {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  if (className) e.className = className;
  return e;
}

function verdictSentence(data) {
  const label = (c) => CHECK_LABELS[c.name] || c.name;
  const bad = data.checks.filter((c) => c.status === 'fail').map(label);
  const unclear = data.checks.filter((c) => c.status === 'inconclusive').map(label);
  const s = data.sites[0];
  if (data.verdict === 'confirmed') {
    const where = s ? `${s.construct} joins ${s.contig} after base ${s.left} and before base ${s.right} (${s.orientation} orientation), replacing ${s.replaced_bases} genome bases` : 'the cassette was found';
    return `Confirmed. ${where.charAt(0).toUpperCase()}${where.slice(1)}; every check passed.`;
  }
  if (data.verdict === 'not_confirmed') return `Not confirmed. ${bad.length === 1 ? 'This check failed' : 'These checks failed'}: ${bad.join('; ')}. The details are below.`;
  return `Inconclusive. There is not enough coverage to decide: ${unclear.join('; ')}.${bad.length ? ` These checks failed: ${bad.join('; ')}.` : ''}`;
}

function siteTable(data) {
  const wrap = el('div');
  if (!data.sites.length) { wrap.appendChild(el('p', 'No read joins the genome to the cassette.', 'hint')); return wrap; }
  const table = el('table', undefined, 'grid');
  const head = table.createTHead().insertRow();
  for (const h of ['Genome contig', 'Cassette', 'Left flank', 'Right flank', 'Orientation', 'Genome bases replaced', 'Shared bases at the junctions']) head.appendChild(el('th', h));
  const body = table.createTBody();
  for (const s of data.sites) {
    const row = body.insertRow();
    const flank = (pos, reads, at) => (pos === null ? 'not found' : `${pos} (${reads} reads; cassette position ${at === null ? '?' : at})`);
    const cells = [s.contig, s.construct, flank(s.left, s.reads_left, s.construct_at_left), flank(s.right, s.reads_right, s.construct_at_right), s.orientation,
      s.replaced_bases === null ? 'unknown' : String(s.replaced_bases), `${Math.max(0, s.homology_left)} left, ${Math.max(0, s.homology_right)} right`];
    for (const c of cells) row.insertCell().textContent = c;
  }
  wrap.appendChild(table);
  return wrap;
}

/** The depth along the first cassette as an SVG: the depth, the genome median, the gaps, the junctions and the differences from the given sequence. */
function depthPlot(data, depthTsv) {
  const cassette = data.cassettes[0];
  const depth = depthTsv.split('\n').slice(1).filter(Boolean).map((l) => Number(l.split('\t')[2]));
  if (!cassette || !depth.length) return null;
  const NS = 'http://www.w3.org/2000/svg';
  const make = (tag, attrs, text) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {})) e.setAttribute(k, String(v));
    if (text !== undefined) e.textContent = text;
    return e;
  };
  const W = 900, H = 200, L = 54, R = 16, T = 18, B = 34;
  const n = depth.length;
  const top = Math.max(10, data.genome_depth * 1.3, ...depth);
  const x = (pos) => L + ((pos - 1) / Math.max(1, n - 1)) * (W - L - R);
  const y = (d) => T + (1 - Math.min(d, top) / top) * (H - T - B);
  const svg = make('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', class: 'depth-plot' });
  const label = `Depth along ${cassette.name} (${n} bases): median ${cassette.median_depth}, genome median ${data.genome_depth}; ${cassette.gaps.length} gap(s); ${cassette.variants.length} difference(s) from the given sequence.`;
  svg.setAttribute('aria-label', label);
  svg.appendChild(make('title', {}, label));
  for (const [a, b] of cassette.gaps) svg.appendChild(make('rect', { class: 'gap', x: x(a), y: T, width: Math.max(2, x(b) - x(a)), height: H - T - B }));
  // The depth, as the mean of buckets so that a long cassette has at most one point per pixel.
  const buckets = Math.min(n, W - L - R);
  let d = '';
  for (let i = 0; i < buckets; i++) {
    const from = Math.floor((i * n) / buckets), to = Math.max(from + 1, Math.floor(((i + 1) * n) / buckets));
    let sum = 0;
    for (let k = from; k < to; k++) sum += depth[k];
    d += `${i === 0 ? 'M' : 'L'}${x(from + 1 + (to - from - 1) / 2).toFixed(1)},${y(sum / (to - from)).toFixed(1)}`;
  }
  svg.appendChild(make('path', { class: 'area', d: `${d}L${x(n).toFixed(1)},${y(0)}L${x(1).toFixed(1)},${y(0)}Z` }));
  svg.appendChild(make('path', { class: 'line', d }));
  svg.appendChild(make('line', { class: 'median', x1: L, x2: W - R, y1: y(data.genome_depth), y2: y(data.genome_depth) }));
  svg.appendChild(make('text', { x: W - R, y: y(data.genome_depth) - 4, 'text-anchor': 'end' }, `genome median ${data.genome_depth}x`));
  const s = data.sites[0];
  for (const [pos, name] of [[s && s.construct_at_left, 'left flank'], [s && s.construct_at_right, 'right flank']]) {
    if (!pos) continue;
    svg.appendChild(make('line', { class: 'junction', x1: x(pos), x2: x(pos), y1: T, y2: H - B }));
    svg.appendChild(make('text', { x: x(pos), y: T - 5, 'text-anchor': pos > n / 2 ? 'end' : 'start' }, name));
  }
  for (const v of cassette.variants) svg.appendChild(make('path', { class: 'variant', d: `M${x(v.pos) - 4},${H - B + 12}L${x(v.pos) + 4},${H - B + 12}L${x(v.pos)},${H - B + 3}Z` }));
  svg.appendChild(make('line', { class: 'axis', x1: L, x2: W - R, y1: H - B, y2: H - B }));
  svg.appendChild(make('line', { class: 'axis', x1: L, x2: L, y1: T, y2: H - B }));
  for (const [pos, anchor] of [[1, 'start'], [Math.round(n / 2), 'middle'], [n, 'end']]) svg.appendChild(make('text', { x: x(pos), y: H - 6, 'text-anchor': anchor }, String(pos)));
  svg.appendChild(make('text', { x: L - 6, y: T + 8, 'text-anchor': 'end' }, `${Math.round(top)}x`));
  svg.appendChild(make('text', { x: L - 6, y: H - B, 'text-anchor': 'end' }, '0'));
  return svg;
}

function render(m) {
  const data = JSON.parse(m.json);
  $('confirm-result').hidden = false;
  const verdict = $('confirm-verdict');
  verdict.dataset.verdict = data.verdict;
  verdict.textContent = verdictSentence(data);
  const warnings = $('confirm-warnings');
  warnings.replaceChildren(...data.warnings.map((w) => el('p', `Warning: ${w}`)));
  warnings.hidden = !data.warnings.length;
  const notes = $('confirm-notes');
  notes.replaceChildren(...data.notes.map((n) => el('p', `Note: ${n}`)));
  notes.hidden = !data.notes.length;
  const target = data.target ? `The target is ${data.target.contig}:${data.target.start}-${data.target.end}${data.target_inferred ? ' (inferred from the homology arms of the construct)' : ''}. ` : 'No target was given. ';
  $('confirm-sites').replaceChildren(el('p', `${target}${m.reads.toLocaleString()} reads; genome median depth ${data.genome_depth}x.`, 'hint'), siteTable(data));
  const plot = depthPlot(data, m.depth);
  $('confirm-plot').replaceChildren(...(plot ? [plot] : [el('p', 'No depth to show.', 'hint')]));
  show('confirm-plot-note', 'The shaded stretches have too little coverage. A line marks each junction. A small triangle marks a difference from the given cassette sequence.');
  const body = $('confirm-checks').tBodies[0];
  body.replaceChildren();
  for (const c of data.checks) {
    const row = body.insertRow();
    row.dataset.check = c.name;
    row.dataset.status = c.status === 'inconclusive' ? 'inconclusive' : c.status;
    row.insertCell().textContent = CHECK_LABELS[c.name] || c.name;
    const [text, klass] = CHIPS[c.status] || [c.status, 'na'];
    row.insertCell().appendChild(el('span', text, `chip ${klass}`));
    row.insertCell().textContent = c.detail;
  }
  offer('confirm-download-text', 'confirm.txt', m.text, 'text/plain');
  offer('confirm-download-json', 'confirm.json', m.json, 'application/json');
  offer('confirm-download-junctions', 'junctions.tsv', m.junctions, 'text/tab-separated-values');
  offer('confirm-download-depth', 'cassette_depth.tsv', m.depth, 'text/tab-separated-values');
  $('confirm-result').scrollIntoView({ block: 'nearest' });
}

// ---- The locus drawing: from the reads (above) or from a breakpoint table -------------------------------------------------------------------
/** Put `text` on the clipboard. Falls back to selecting nothing when the browser refuses (the text is shown in the page anyway). */
async function copyText(text, button) {
  try {
    await navigator.clipboard.writeText(text);
    button.textContent = 'Copied';
  } catch (error) {
    button.textContent = 'Select and copy';
  }
  setTimeout(() => { button.textContent = 'Copy'; }, 1500);
}

/** A cell with a region as selectable text and a Copy button. */
function regionCell(row, text) {
  const cell = row.insertCell();
  cell.appendChild(el('code', text, 'region'));
  const button = el('button', 'Copy', 'copy');
  button.type = 'button';
  button.setAttribute('aria-label', `Copy ${text}`);
  button.addEventListener('click', () => copyText(text, button));
  cell.appendChild(button);
}

/** Draw the result of `locus_json` (a JSON text, or "null"). */
function renderLocus(json, caption) {
  const section = $('locus-section');
  const data = json ? JSON.parse(json) : null;
  if (!data) { section.hidden = true; return; }
  section.hidden = false;
  show('locus-caption', caption);
  // The drawing comes from the WebAssembly package as SVG text; it is parsed as SVG, so no script can run in it.
  const doc = new DOMParser().parseFromString(data.svg, 'image/svg+xml');
  $('locus-svg').replaceChildren(document.importNode(doc.documentElement, true));
  const regions = $('locus-regions').tBodies[0];
  regions.replaceChildren();
  for (const r of data.regions) {
    const row = regions.insertRow();
    row.insertCell().textContent = r.name;
    regionCell(row, r.region);
    regionCell(row, r.region_dots);
  }
  const genes = $('locus-genes');
  const body = genes.tBodies[0];
  body.replaceChildren();
  genes.hidden = !data.genes.length;
  show('locus-genes-note', data.genes.length
    ? `The links go to fungidb.org. They work for gene IDs that FungiDB uses (the gene IDs of a FungiDB genome); for other genomes FungiDB shows its not-found page.`
    : 'No gene overlaps the replaced bases, so there is no gene to link.');
  for (const g of data.genes) {
    const row = body.insertRow();
    const cell = row.insertCell();
    const a = el('a', g.id);
    a.href = `${GENE_LINK}${encodeURIComponent(g.id)}`;
    a.target = '_blank';
    a.rel = 'noopener noreferrer';
    cell.appendChild(a);
    row.insertCell().textContent = g.name || '';
    row.insertCell().textContent = g.product || '';
    regionCell(row, g.region);
    row.insertCell().textContent = g.strand;
    row.insertCell().textContent = `${g.removed_bases.toLocaleString()} of ${g.gene_bases.toLocaleString()}`;
    row.insertCell().textContent = `${g.coding_removed.toLocaleString()} of ${g.coding_bases.toLocaleString()}`;
    row.insertCell().textContent = g.effect;
  }
  const tsv = (header, rows) => `${header.join('\t')}\n${rows.map((r) => r.join('\t')).join('\n')}\n`;
  offer('locus-download-svg', 'insertion_locus.svg', data.svg, 'image/svg+xml');
  offer('locus-download-genes', 'insertion_genes.tsv', tsv(['gene', 'name', 'product', 'region', 'strand', 'gene_bases', 'removed_bases', 'coding_bases', 'coding_removed', 'effect'],
    data.genes.map((g) => [g.id, g.name || '', g.product || '', g.region, g.strand, g.gene_bases, g.removed_bases, g.coding_bases, g.coding_removed, g.effect])), 'text/tab-separated-values');
  offer('locus-download-regions', 'insertion_regions.tsv', tsv(['name', 'region', 'region_dots'], data.regions.map((r) => [r.name, r.region, r.region_dots])), 'text/tab-separated-values');
}

// ---- Draw from a breakpoint table or an assembly (no reads) ------------------------------------------------------------------------------------
let breakpointTable = '';

/** Fill the row menu from the table text. Chooses the first row that overlaps the locus. */
async function loadRows(table) {
  breakpointTable = table;
  $('bp-row').replaceChildren();
  $('bp-row').disabled = true;
  $('bp-draw').disabled = true;
  const m = await callWorker({ step: 'breakpoints', table });
  const rows = JSON.parse(m.rows);
  if (!rows.length) { show('bp-status', 'No gap between alignments was found in this table.'); return rows; }
  $('bp-row').replaceChildren(...rows.map((r) => {
    const option = el('option', `${r.sample || r.contig} · ${r.ref_scaffold}:${r.left}-${r.right} · ${r.query_gap_bp.toLocaleString()} bp in place of ${r.ref_gap_bp.toLocaleString()} bp${r.overlaps_locus ? ' · overlaps the locus' : ''}`);
    option.value = String(r.index);
    return option;
  }));
  const first = rows.find((r) => r.overlaps_locus) || rows[0];
  $('bp-row').value = String(first.index);
  $('bp-row').disabled = false;
  $('bp-draw').disabled = false;
  return rows;
}

async function drawRow() {
  const gff = gffFile();
  const m = await callWorker({ step: 'breakpoints', table: breakpointTable, row: Number($('bp-row').value), gff: gff || null, geneLink: GENE_LINK });
  renderLocus(m.locus, gff
    ? 'The row of the breakpoint table, drawn in the genes of the reference.'
    : 'The row of the breakpoint table. No genes file was chosen (section 1), so no genes are drawn.');
  $('locus-section').scrollIntoView({ block: 'nearest' });
}

$('bp-file').addEventListener('change', async () => {
  const file = $('bp-file').files[0];
  breakpointTable = '';
  $('bp-row').replaceChildren();
  $('bp-row').disabled = true;
  $('bp-draw').disabled = true;
  show('bp-status', '');
  if (!file) return;
  try {
    if (file.size > 50e6) throw new Error('the file is larger than 50 MB');
    const rows = await loadRows(await file.text());
    if (rows.length) show('bp-status', `${rows.length} row${rows.length === 1 ? '' : 's'} read.`);
  } catch (error) {
    breakpointTable = '';
    show('bp-status', `Error: ${error.message}`);
  }
});

$('bp-draw').addEventListener('click', async () => {
  if (!breakpointTable || running || loading) { if (loading) show('bp-status', 'Wait until the genome has loaded.'); return; }
  show('bp-status', 'Drawing');
  try {
    await drawRow();
    show('bp-status', 'Done.');
  } catch (error) {
    show('bp-status', `Error: ${error.message}`);
  }
});

/** The locus box: `SCAFFOLD:START-END` or `SCAFFOLD:START..END` is a region; anything else is a gene name. */
function locusOfForm() {
  const text = $('bp-locus').value.trim();
  if (!text) throw new Error('Type the locus: a region such as scf_000011:282,234-285,136, or a gene name.');
  const m = text.match(/^(\S+):([\d,]+)(?:-|\.\.)([\d,]+)$/);
  if (!m) return { gene: text };
  const [start, end] = [Number(m[2].replace(/,/g, '')), Number(m[3].replace(/,/g, ''))];
  if (!(start >= 1 && end >= start)) throw new Error('Type a region whose start is 1 or more and not above its end.');
  return { locus: { contig: m[1], start, end } };
}

$('bp-find').addEventListener('click', async () => {
  if (running || loading) { show('bp-status', loading ? 'Wait until the genome has loaded.' : 'Wait until the other run is done.'); return; }
  const assembly = $('bp-assembly').files[0];
  const reference = refFile();
  if (!reference) { show('bp-status', 'Choose the reference genome first (section 1).'); return; }
  if (!assembly) { show('bp-status', 'Choose the assembly (a FASTA file).'); return; }
  let wanted;
  try { wanted = locusOfForm(); } catch (error) { show('bp-status', error.message); return; }
  const gff = gffFile();
  if (wanted.gene && !gff) { show('bp-status', 'Choose the GFF3 file that holds the gene (section 1), or type a region.'); return; }
  running = true;
  $('bp-find').disabled = true;
  try {
    const sample = assembly.name.replace(/\.gz$/, '').replace(/\.(fasta|fa|fna)$/, '').replace(/\.sorted$/, '') || 'assembly';
    const m = await callWorker({ step: 'assembly', reference, assembly, gff: gff || null, sample, preset: $('bp-preset').value, ...wanted },
      (p) => show('bp-status', p.text || ''));
    const result = JSON.parse(m.result);
    offer('bp-download-table', 'breakpoint_summary.tsv', result.table, 'text/tab-separated-values');
    offer('bp-download-fasta', 'inserted_sequences.fasta', result.fasta, 'text/plain');
    const where = `${m.locus.contig}:${m.locus.start.toLocaleString()}-${m.locus.end.toLocaleString()}`;
    if (result.contigs_near_locus === 0) { show('bp-status', `No contig of the assembly aligns near ${where}: the assembly does not reach this locus (or it is too different; try a wider setting).`); return; }
    const rows = await loadRows(result.table);
    const hits = rows.filter((r) => r.overlaps_locus);
    if (!hits.length) { show('bp-status', `${result.contigs_near_locus} contig(s) align near ${where}, but no insertion or deletion overlaps it. Other gaps nearby: ${rows.length}.`); return; }
    await drawRow();
    show('bp-status', `${hits.length} gap${hits.length === 1 ? '' : 's'} overlap${hits.length === 1 ? 's' : ''} ${where}.`);
  } catch (error) {
    show('bp-status', `Error: ${error.message}`);
  } finally {
    running = false;
    $('bp-find').disabled = false;
  }
});
