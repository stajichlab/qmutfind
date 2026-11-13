// The allele confirmation page: map the reads to the genome plus the introduced cassette (in a worker, in WebAssembly) and say where the cassette
// joins the genome and whether it is intact. The files stay in this browser; web addresses are fetched by this page.
const $ = (id) => document.getElementById(id);
const show = (id, text) => { $(id).textContent = text; };
const mb = (bytes) => `${(bytes / 1e6).toFixed(1)} MB`;

// How far a file may unpack when it is loaded whole into memory, and the largest download. The query options lower them (the tests use them).
const LIMIT_BYTES = Math.min(parseFloat(new URLSearchParams(location.search).get('limit_mb')) || 1024, 1024) * 1e6;
const MAX_DOWNLOAD = Math.min(parseFloat(new URLSearchParams(location.search).get('download_limit_mb')) || 2000, 2000) * 1e6;

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

// ---- Choices that show or hide parts of the form --------------------------------------------------------------------------------------
function showSource() {
  const source = $('confirm-genome-source').value;
  for (const group of document.querySelectorAll('[data-source-group]')) {
    const name = group.dataset.sourceGroup;
    group.hidden = name === 'gff-file' ? source === 'ncbi' : name !== source;
  }
}
function showTarget() {
  const kind = $('confirm-target-kind').value;
  for (const group of document.querySelectorAll('[data-target-group]')) group.hidden = group.dataset.targetGroup !== kind;
}
$('confirm-genome-source').addEventListener('change', showSource);
$('confirm-target-kind').addEventListener('change', showTarget);
showSource();
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

// ---- The NCBI list ----------------------------------------------------------------------------------------------------------------------
let genomeList = [];
(async function loadList() {
  const select = $('confirm-ncbi');
  try {
    const response = await fetch('genomes/genomes.json', { cache: 'no-cache' });
    if (!response.ok) throw new Error(`answer ${response.status}`);
    genomeList = (await response.json()).genomes || [];
    select.replaceChildren(new Option('Choose a genome', ''));
    for (const g of genomeList) select.appendChild(new Option(`${g.label} (${g.accession}${g.gff ? '' : ', no gene models'})`, g.accession));
  } catch (error) {
    select.replaceChildren(new Option('The list could not be loaded', ''));
  }
})();

// ---- Downloads --------------------------------------------------------------------------------------------------------------------------
/** An address the page may download from: https (or this site), and no user name or password. */
function checkUrl(text, label) {
  let url;
  try { url = new URL(text, location.href); } catch (error) { throw new Error(`${label}: "${text}" is not an address.`); }
  if (url.protocol !== 'https:' && url.origin !== location.origin) throw new Error(`${label}: only https addresses are accepted, not ${url.protocol}//`);
  if (url.username || url.password) throw new Error(`${label}: addresses with a user name or password are not accepted.`);
  return url;
}

async function download(item, label, signal, budget) {
  const url = checkUrl(item.url, label);
  const response = await fetch(url, { signal });
  if (!response.ok) throw new Error(`${label}: the server answered ${response.status} for ${url.host}`);
  const expected = item.bytes || Number(response.headers.get('Content-Length')) || 0;
  if (expected > budget.left) throw new Error(`${label}: the download is too large (${mb(expected)}; the limit is ${mb(budget.left)}).`);
  const reader = response.body.getReader();
  const parts = [];
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    parts.push(value);
    got += value.length;
    if (got > budget.left) { await reader.cancel(); throw new Error(`${label}: the download is too large (more than ${mb(budget.left)}).`); }
    show('confirm-status', `Downloading ${label}: ${mb(got)}${expected ? ` of ${mb(expected)}` : ''}`);
  }
  if (item.bytes && got !== item.bytes) throw new Error(`${label}: got ${got} bytes, expected ${item.bytes}. The download is incomplete.`);
  budget.left -= got;
  let name = item.name || decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'download');
  if (parts.length && parts[0].length > 1 && parts[0][0] === 0x1f && parts[0][1] === 0x8b && !name.endsWith('.gz')) name += '.gz';
  return new File(parts, name);
}

/** Run one message in a worker of its own. The worker is stopped when the answer comes, or when the run is cancelled. */
function callWorker(message, onProgress) {
  return new Promise((resolve, reject) => {
    const w = new Worker('worker.js?v=__BUILD__', { type: 'module' });
    activeWorker = { w, reject };
    const finish = () => { w.terminate(); if (activeWorker && activeWorker.w === w) activeWorker = null; };
    w.onmessage = (event) => {
      const m = event.data;
      if (m.type === 'status' && message.step === 'confirm') show('confirm-status', m.text);
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
  const source = $('confirm-genome-source').value;
  const refFile = $('confirm-reference').files[0];
  const refUrl = $('confirm-ref-url').value.trim();
  const ncbi = genomeList.find((g) => g.accession === $('confirm-ncbi').value);
  if ((source === 'upload' && !refFile) || (source === 'url' && !refUrl) || (source === 'ncbi' && !ncbi)) {
    show('confirm-status', 'Choose the reference genome first (section 1).');
    return;
  }
  let construct;
  try { construct = await constructText(); } catch (error) { show('confirm-status', `Give the construct: the file cannot be read (${error.message}).`); return; }
  if (!construct) { show('confirm-status', 'Give the construct: a FASTA file or the pasted sequence of the cassette (section 1).'); return; }
  let wanted;
  try { wanted = targetOfForm(); } catch (error) { show('confirm-status', error.message); return; }
  const gffFile = $('confirm-gff').files[0];
  if (wanted.gene && source !== 'ncbi' && !gffFile) { show('confirm-status', 'Choose the GFF3 file that holds the gene: a gene name needs it (section 1).'); return; }
  if (wanted.gene && source === 'ncbi' && !ncbi.gff) { show('confirm-status', 'This genome has no gene models: choose a GFF3 file, or give a region instead.'); return; }
  const localReads = [...$('confirm-reads').files];
  const addresses = $('confirm-reads-urls').value.split('\n').map((x) => x.trim()).filter(Boolean);
  if (!localReads.length && !addresses.length) { show('confirm-status', 'Choose the reads: files from your computer or https addresses (section 3).'); return; }
  let urls;
  try {
    urls = addresses.map((a) => checkUrl(a, 'the reads').href); // all checked before any download
    if (source === 'url') checkUrl(refUrl, 'the reference genome');
  } catch (error) { show('confirm-status', `Error: ${error.message}`); return; }

  running = true;
  cancelled = false;
  $('confirm-run').disabled = true;
  $('confirm-cancel').hidden = false;
  $('confirm-result').hidden = true;
  activeAbort = new AbortController();
  try {
    const budget = { left: MAX_DOWNLOAD };
    let reference = refFile;
    let gff = gffFile || null;
    if (source === 'url') {
      const url = checkUrl(refUrl, 'the reference genome');
      reference = await download({ url: url.href }, 'the reference genome', activeAbort.signal, budget);
    } else if (source === 'ncbi') {
      const base = new URL('genomes/genomes.json', location.href);
      reference = await download({ ...ncbi.fasta, url: new URL(ncbi.fasta.url, base).href }, 'the genome', activeAbort.signal, budget);
      if (ncbi.gff) gff = await download({ ...ncbi.gff, url: new URL(ncbi.gff.url, base).href }, 'the gene models', activeAbort.signal, budget);
    }
    const limitThousands = number('confirm-read-limit', 0);
    const maxReads = Math.round(limitThousands * 1000);
    const reads = [...localReads];
    for (const [i, address] of urls.entries()) {
      const label = `reads file ${i + 1} of ${urls.length}`;
      reads.push(maxReads ? await downloadFirstReads(address, label, maxReads, budget)
        : await download({ url: address }, label, activeAbort.signal, budget));
    }
    show('confirm-status', 'Starting');
    const m = await callWorker({
      step: 'confirm', reference, construct: new TextEncoder().encode(construct), reads, kind: $('confirm-kind').value, gff,
      ...wanted, maxReads,
      options: { minReads: Math.max(1, Math.round(number('confirm-min-reads', 3))), minDepth: number('confirm-min-depth', 10), tolerance: Math.round(number('confirm-tolerance', 20)) },
    }, (p) => show('confirm-status', `Mapping: ${p.reads.toLocaleString()} reads`));
    render(m);
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
  const target = data.target ? `The target is ${data.target.contig}:${data.target.start}-${data.target.end}. ` : 'No target was given. ';
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
