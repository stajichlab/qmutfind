// The page: three steps (call, summarize, draw). The result of each step is kept here and the next step uses it.
const $ = (id) => document.getElementById(id);
const worker = new Worker('worker.js', { type: 'module' });

// The state that carries from step to step.
const state = { vcf: null, annotated: null, depth: null };
let busy = false;
const MAX_ROWS = 500;
// How far a file may unpack when it is loaded whole into memory. `?limit_mb=N` can lower it (the tests use that).
const LIMIT_BYTES = Math.min(parseFloat(new URLSearchParams(location.search).get('limit_mb')) || 1024, 1024) * 1e6;
// A download is limited to 2 GB. `?download_limit_mb=N` can lower the limit (the tests use that).
const MAX_DOWNLOAD = Math.min(parseFloat(new URLSearchParams(location.search).get('download_limit_mb')) || 2000, 2000) * 1e6;
const WARN_DOWNLOAD = 5e8;

const show = (id, text) => { $(id).textContent = text; };
const buttons = ['run', 'summarize', 'annotate', 'draw'];
function setBusy(value) {
  busy = value;
  for (const id of buttons) $(id).disabled = value;
}

function offer(id, name, text, type) {
  const link = $(id);
  if (link.dataset.url) URL.revokeObjectURL(link.dataset.url);
  const url = URL.createObjectURL(new Blob([text], { type }));
  link.dataset.url = url;
  link.href = url;
  link.download = name;
  link.hidden = false;
}

const sampleName = () => $('sample').value || 'sample';

// Files from an example. A file chosen in an input always wins over the example's file.
const example = { reference: null, gff: null, reads: [], expectedReads: null, source: { reference: '', gff: '', reads: '' } };

/** Say, under each input, where its data comes from: a local file, or a download (an example, the menu or the web import). */
function describeSources() {
  const size = (f) => `${(f.size / 1e6).toFixed(f.size < 1e5 ? 3 : 1)} MB`;
  const line = (slot, inputId, loaded) => {
    const chosen = [...$(inputId).files];
    const el = $(`source-${slot}`);
    el.classList.remove('from-web');
    if (chosen.length) {
      el.textContent = `Local file: ${chosen.map((f) => `${f.name} (${size(f)})`).join(', ')}`;
    } else if (loaded.length) {
      el.textContent = `From ${example.source[slot]}: ${loaded.map((f) => `${f.name} (${size(f)})`).join(', ')}`;
      el.classList.add('from-web');
    } else {
      el.textContent = 'Not set. Choose a file, or use an example, the genome menu or the web import.';
    }
  };
  line('reference', 'reference', example.reference ? [example.reference] : []);
  line('gff', 'gff', example.gff ? [example.gff] : []);
  line('reads', 'reads', example.reads);
}
for (const id of ['reference', 'gff', 'reads']) $(id).addEventListener('change', describeSources);
const hostOf = (item, base) => new URL(item.url, base).host;
const refFile = () => $('reference').files[0] || example.reference;
const gffFile = () => $('gff').files[0] || example.gff;
const readFiles = () => ($('reads').files.length ? [...$('reads').files] : example.reads);

/** The VCF the next step should use: a chosen file, else the annotated VCF, else the step 1 result. */
function currentVcf() {
  const file = $('vcf-file').files[0];
  if (file) return file;
  return state.annotated || state.vcf;
}

/** The coverage the next step should use: a chosen file, else the BED from step 1. */
function currentDepth() {
  return $('depth-file').files[0] || state.depth;
}

function describeSource() {
  const file = $('vcf-file').files[0];
  let text;
  if (file) text = `Using the VCF file ${file.name}.`;
  else if (state.annotated) text = 'Using the annotated VCF from step 3.';
  else if (state.vcf) text = 'Using the VCF from step 1.';
  else text = 'No VCF yet. Run step 1 or choose a VCF file.';
  show('vcf-source', text);
}

function send(step, payload, panelStatus) {
  if (busy) return false;
  setBusy(true);
  show(panelStatus, 'Starting');
  worker.postMessage({ step, limit: LIMIT_BYTES, ...payload });
  return true;
}

function needVcf(statusId) {
  const vcf = currentVcf();
  if (!vcf) show(statusId, 'There is no VCF yet. Run step 1, or choose a VCF file.');
  return vcf;
}

// Step 1
$('run').addEventListener('click', () => {
  const reference = refFile();
  const reads = readFiles();
  if (!reference || reads.length === 0) {
    show('status', 'Choose a reference FASTA and at least one reads file.');
    return;
  }
  send('run', { reference, reads, sample: sampleName(), kind: $('kind').value }, 'status');
});

// Step 2
$('summarize').addEventListener('click', () => {
  const vcf = currentVcf();
  const bed = currentDepth();
  if (!vcf && !bed) { show('summarize-status', 'There is no VCF or coverage yet. Run step 1, or choose a VCF or coverage file.'); return; }
  send('summarize', { vcf, bed, reference: refFile() || null, order: $('order').value }, 'summarize-status');
});

function statsTable(s) {
  const rows = [];
  const add = (name, value) => rows.push(`<tr><th scope="row">${name}</th><td>${value}</td></tr>`);
  const n = (x) => Number(x).toLocaleString();
  add('Variants', n(s.total));
  for (const [k, label] of [['SNP', 'SNPs'], ['SMALL_INDEL', 'Small indels (1 to 49 bp)'], ['LARGE_INDEL', 'Large indels (50 bp or more)']]) {
    add(label, n((s.by_class || {})[k] || 0));
  }
  add('Transitions / transversions', s.ts_tv === null ? 'not defined' : `${n(s.transitions)} / ${n(s.transversions)} (Ts/Tv ${s.ts_tv.toFixed(2)})`);
  add('Variants per kb', s.overall_per_kb === null ? 'not defined' : s.overall_per_kb.toFixed(4));
  add('Calls with an allele fraction from 0.2 to 0.8', `${n(s.mid_fraction_calls)}${s.mid_fraction_share === null ? '' : ` (${(100 * s.mid_fraction_share).toFixed(1)}%)`}`);
  add('Calls within 100 bp of another call', n(s.clustered));
  const effects = Object.entries(s.by_effect || {}).filter(([k]) => k !== 'unannotated');
  if (effects.length) {
    add('Effects', effects.sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k}: ${n(v)}`).join(', '));
  }
  return `<table class="stats">${rows.join('')}</table>`;
}

// Step 3
$('annotate').addEventListener('click', () => {
  const vcf = needVcf('annotate-status');
  if (!vcf) return;
  const reference = refFile();
  const gff = gffFile();
  if (!reference || !gff) {
    show('annotate-status', 'Choose the reference genome and the GFF3 file in "Your files".');
    return;
  }
  send('annotate', { vcf, reference, gff }, 'annotate-status');
});

/** A VCF INFO value with its %XX escapes undone (qmutfind writes `;`, `=`, `,`, spaces and `%` that way). */
function decodeInfo(value) {
  if (!value) return '';
  try { return decodeURIComponent(value); } catch (error) { return value; }
}

/** The variants in a VCF as table rows. Uses the INFO tags written by qmutfind annotate. */
function parseVariants(vcfText) {
  const rows = [];
  for (const line of vcfText.split('\n')) {
    if (!line || line[0] === '#') continue;
    const f = line.split('\t');
    if (f.length < 8) continue;
    const info = {};
    for (const kv of f[7].split(';')) {
      const i = kv.indexOf('=');
      if (i > 0) info[kv.slice(0, i)] = kv.slice(i + 1);
    }
    rows.push({
      contig: f[0], pos: Number(f[1]), ref: f[3], alt: f[4].split(',')[0],
      cls: info.CLASS || '', effect: (info.EFFECT || 'unannotated').split(',')[0],
      gene: info.GENE || '', geneName: decodeInfo(info.GENENAME), product: decodeInfo(info.PRODUCT), aa: info.AA && info.AA !== '.' ? info.AA : '', af: info.AF || '', dp: info.DP || '',
    });
  }
  return rows;
}

let variants = [];
const cell = (text) => { const td = document.createElement('td'); td.textContent = text; return td; };

function renderTable() {
  const filter = $('effect-filter').value;
  const shown = variants.filter((v) => !filter || v.effect === filter);
  const table = $('variant-table');
  table.replaceChildren();
  const head = table.createTHead().insertRow();
  for (const name of ['Position', 'Change', 'Class', 'Effect', 'Gene', 'Name', 'Protein change', 'Product', 'AF', 'DP', '']) {
    const th = document.createElement('th');
    th.textContent = name;
    head.appendChild(th);
  }
  const body = table.createTBody();
  for (const v of shown.slice(0, MAX_ROWS)) {
    const row = body.insertRow();
    row.append(cell(`${v.contig}:${v.pos.toLocaleString()}`), cell(`${v.ref.length > 12 ? v.ref.slice(0, 12) + '…' : v.ref} > ${v.alt.length > 12 ? v.alt.slice(0, 12) + '…' : v.alt}`),
      cell(v.cls), cell(v.effect), cell(v.gene), cell(v.geneName), cell(v.aa), cell(v.product), cell(v.af), cell(v.dp));
    const last = row.insertCell();
    const button = document.createElement('button');
    button.type = 'button';
    button.textContent = 'Draw';
    button.className = 'draw-row';
    button.dataset.contig = v.contig;
    button.dataset.pos = String(v.pos);
    last.appendChild(button);
  }
  table.hidden = false;
  const note = shown.length > MAX_ROWS ? ` Showing the first ${MAX_ROWS}.` : '';
  show('annotate-status', `${shown.length.toLocaleString()} of ${variants.length.toLocaleString()} variants.${note}`);
}

$('effect-filter').addEventListener('change', renderTable);

function loadVariants(vcfText) {
  variants = parseVariants(vcfText);
  const counts = {};
  for (const v of variants) counts[v.effect] = (counts[v.effect] || 0) + 1;
  const select = $('effect-filter');
  select.replaceChildren(new Option(`All variants (${variants.length})`, ''));
  for (const [effect, count] of Object.entries(counts).sort((a, b) => b[1] - a[1])) select.add(new Option(`${effect} (${count})`, effect));
  $('variant-controls').hidden = false;
  renderTable();
  const genic = variants.some((v) => !['intergenic', 'unannotated'].includes(v.effect));
  if (variants.length && !genic) {
    show('annotate-status', $('annotate-status').textContent + ' No variant fell in a gene. Check that the GFF3 uses the same contig names as the reference.');
  }
}

function draw(region) {
  const vcf = state.annotated || currentVcf();
  if (!vcf) { show('report-status', 'There is no VCF yet.'); return; }
  send('report', { vcf, gff: gffFile() || null, region }, 'report-status');
}

$('variant-table').addEventListener('click', (event) => {
  const button = event.target.closest('button.draw-row');
  if (!button) return;
  const pos = Number(button.dataset.pos);
  const half = Math.max(50, Math.round((Number($('window-kb').value) || 2.5) * 1000));
  const region = `${button.dataset.contig}:${Math.max(1, pos - half)}-${pos + half}`;
  $('region').value = region;
  draw(region);
});
$('draw').addEventListener('click', () => {
  const region = $('region').value.trim();
  if (!region) { show('report-status', 'Type a region such as Chr1:1,000-6,000, or press Draw on a variant.'); return; }
  draw(region);
});

// A new VCF file replaces the state of the earlier steps.
$('vcf-file').addEventListener('change', () => { state.annotated = null; describeSource(); });

worker.onmessage = (event) => {
  const m = event.data;
  const panel = { run: 'status', summarize: 'summarize-status', annotate: 'annotate-status', report: 'report-status' }[m.step];
  if (m.type === 'status') show(panel, m.text);
  if (m.type === 'progress') show('progress', `${m.reads.toLocaleString()} reads mapped`);
  if (m.type === 'error') {
    show(panel, `Error: ${m.message}`);
    setBusy(false);
    return;
  }
  if (m.type !== 'done') return;
  setBusy(false);
  if (m.step === 'run') {
    state.vcf = m.vcf;
    state.depth = m.depth;
    state.annotated = null;
    let warning = '';
    if (example.expectedReads && !$('reads').files.length && m.reads !== example.expectedReads) {
      warning = ` WARNING: the example should have ${example.expectedReads.toLocaleString()} reads, so a download may be incomplete.`;
    }
    show('status', `Done: ${m.reads.toLocaleString()} reads in ${m.seconds.toFixed(1)} s.${warning}`);
    show('summary', `Calls: ${m.vcf.split('\n').filter((l) => l && l[0] !== '#').length.toLocaleString()}. Go on to step 2 or step 3.`);
    offer('download-vcf', `${sampleName()}.vcf`, m.vcf, 'text/plain');
    offer('download-depth-bed', `${sampleName()}.depth.bed`, m.depth, 'text/plain');
    describeSource();
  }
  if (m.step === 'summarize') {
    const parts = [];
    if (m.stats) {
      const s = JSON.parse(m.stats);
      $('stats-table').innerHTML = statsTable(s);
      parts.push(`${Number(s.total).toLocaleString()} variants`);
      offer('download-stats', `${sampleName()}.stats.json`, m.stats, 'application/json');
    }
    $('genome-plot').innerHTML = m.plot || '';
    if (m.plot) offer('download-plot', `${sampleName()}.genome.svg`, m.plot, 'image/svg+xml');
    $('depth-plot').innerHTML = m.depthPlot || '';
    if (m.depthPlot) { parts.push('coverage plot'); offer('download-depth-plot', `${sampleName()}.depth.svg`, m.depthPlot, 'image/svg+xml'); }
    show('summarize-status', `Done: ${parts.join(', ')}.`);
  }
  if (m.step === 'annotate') {
    state.annotated = m.vcf;
    offer('download-annotated', `${sampleName()}.annotated.vcf`, m.vcf, 'text/plain');
    describeSource();
    loadVariants(m.vcf);
  }
  if (m.step === 'report') {
    $('report').srcdoc = m.page;
    show('report-status', `Drawn: ${m.region}`);
  }
};

describeSource();


// ---- Examples -------------------------------------------------------------------------------------------------
// The list of examples is example/manifest.json. A page address with ?example=PATH uses another manifest, but only
// a path on this site (no http:// or // address), so a link cannot make the page fetch a manifest from elsewhere.
const mb = (bytes) => `${(bytes / 1e6).toFixed(1)} MB`;
let loading = false;

async function download(item, base, label, statusId = 'example-status', signal = undefined) {
  const url = checkUrl(new URL(item.url, base).href, label);
  const response = await fetch(url, { signal });
  if (!response.ok) throw new Error(`${label}: the server answered ${response.status} for ${url.host}`);
  const announced = Number(response.headers.get('Content-Length')) || 0;
  const expected = item.bytes || announced;
  if (expected > MAX_DOWNLOAD) throw new Error(`${label}: the download is too large (${mb(expected)}; the limit is ${mb(MAX_DOWNLOAD)}).`);
  const reader = response.body.getReader();
  const parts = [];
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    parts.push(value);
    got += value.length;
    if (got > MAX_DOWNLOAD) { await reader.cancel(); throw new Error(`${label}: the download is too large (more than ${mb(MAX_DOWNLOAD)}).`); }
    show(statusId, `Downloading ${label}: ${mb(got)}${expected ? ` of ${mb(expected)}` : ''}`);
  }
  if (item.bytes && got !== item.bytes) throw new Error(`${label}: got ${got} bytes, expected ${item.bytes}. The download is incomplete.`);
  // The page reads a file as gzip by its name. A download whose name has no .gz but whose content is gzip gets the name.
  let name = item.name;
  if (parts.length && parts[0].length > 1 && parts[0][0] === 0x1f && parts[0][1] === 0x8b && !name.endsWith('.gz')) name += '.gz';
  return new File(parts, name);
}

/** An address the page may download from: https (or this site), and no user name or password. */
function checkUrl(text, label) {
  let url;
  try { url = new URL(text); } catch (error) { throw new Error(`${label}: "${text}" is not an address.`); }
  if (url.protocol !== 'https:' && url.origin !== location.origin) throw new Error(`${label}: only https addresses are accepted, not ${url.protocol}//`);
  if (url.username || url.password) throw new Error(`${label}: addresses with a user name or password are not accepted.`);
  return url;
}

/**
 * A Dropbox share link (www.dropbox.com/s/... or /scl/fi/...) answers with a redirect that has no CORS header, so a browser
 * page cannot follow it. The same path on dl.dropboxusercontent.com has the header (measured 2026-11-12), so the host is
 * changed and the `dl` and `raw` options are dropped. Any other address is returned as it is.
 */
function directAddress(text) {
  let url;
  try { url = new URL(text); } catch (error) { return text; }
  const shared = /^\/(s|scl\/fi)\//.test(url.pathname);
  if (url.protocol !== 'https:' || !['www.dropbox.com', 'dropbox.com'].includes(url.hostname) || !shared) return text;
  url.hostname = 'dl.dropboxusercontent.com';
  url.searchParams.delete('dl');
  url.searchParams.delete('raw');
  return url.href;
}
window.qmutfindDirectAddress = directAddress;

async function loadExample(manifest, base, chosen) {
  if (loading || busy) return;
  loading = true;
  for (const b of document.querySelectorAll('#example-buttons button')) b.disabled = true;
  try {
    example.reference = await download(manifest.reference, base, 'the reference genome');
    example.gff = await download(manifest.gff, base, 'the gene models');
    example.reads = [];
    for (const [i, item] of chosen.reads.entries()) {
      example.reads.push(await download(item, base, `reads file ${i + 1} of ${chosen.reads.length}`));
    }
    example.expectedReads = chosen.expected_reads || null;
    const origin = (item) => `the example (${item.source ? item.source + ' at ' : ''}${hostOf(item, base)})`;
    example.source = { reference: origin(manifest.reference), gff: origin(manifest.gff), reads: origin(chosen.reads[0]) };
    describeSources();
    if (manifest.sample) $('sample').value = manifest.sample;
    if (manifest.kind) $('kind').value = manifest.kind;
    show('example-status', `Loaded: ${chosen.label}. Press Run in step 1.`);
  } catch (error) {
    example.reference = example.gff = null;
    example.reads = [];
    describeSources();
    show('example-status', `Error: ${error.message} The server may block downloads from web pages, or the network failed.`);
  } finally {
    loading = false;
    for (const b of document.querySelectorAll('#example-buttons button')) b.disabled = false;
  }
}

async function setupExamples() {
  const asked = new URLSearchParams(location.search).get('example');
  const path = asked && !/^([a-z]+:)?\/\//i.test(asked) && !asked.startsWith('/') ? asked : 'example/manifest.json';
  const base = new URL(path, location.href);
  let manifest;
  try {
    const response = await fetch(base, { cache: 'no-cache' });
    if (!response.ok) return;
    manifest = await response.json();
  } catch (error) {
    return;
  }
  $('example-text').textContent = manifest.description || '';
  const holder = $('example-buttons');
  for (const chosen of manifest.examples || []) {
    const button = document.createElement('button');
    button.type = 'button';
    button.textContent = chosen.label;
    button.dataset.example = chosen.id;
    button.addEventListener('click', () => loadExample(manifest, base, chosen));
    holder.appendChild(button);
  }
  $('example-section').hidden = false;
}

setupExamples();


// ---- The genome menu ------------------------------------------------------------------------------------------
// genomes/genomes.json lists NCBI genomes (RefSeq, and some GenBank without gene models) (the page downloads them: NCBI allows that). genomes/fungidb.json lists
// FungiDB organisms. FungiDB blocks downloads from web pages and needs an API key, so for those the page only shows the
// `qmutfind fetch` command. ?genomes=PATH uses another list, but only a path on this site.
const LARGE_GENOME = 50e6;
let genomeList = [];
let genomesBase = location.href; // file addresses in the list are relative to the list
let fungidbList = [];

function fillMenu() {
  const filter = $('genome-filter').value.trim().toLowerCase();
  const menu = $('genome-menu');
  menu.replaceChildren();
  const group = (label) => { const g = document.createElement('optgroup'); g.label = label; return g; };
  const ncbi = group('NCBI (loads in this page)');
  for (const g of genomeList) {
    if (filter && !`${g.label} ${g.accession}`.toLowerCase().includes(filter)) continue;
    ncbi.appendChild(new Option(`${g.label} (${(g.genome_size / 1e6).toFixed(1)} Mb, ${g.gff && g.genes > 0 ? `${g.genes.toLocaleString()} genes` : 'no gene models'}, ${g.level}${g.source === 'GenBank' ? ', GenBank' : ''})`, `ncbi:${g.accession}`));
  }
  const fungidb = group('FungiDB (needs the command line and a key)');
  for (const o of fungidbList) {
    if (filter && !o.name.toLowerCase().includes(filter)) continue;
    fungidb.appendChild(new Option(`${o.name} (${o.genes.toLocaleString()} genes)`, `fungidb:${o.name}`));
  }
  if (ncbi.children.length) menu.appendChild(ncbi);
  if (fungidb.children.length) menu.appendChild(fungidb);
}

async function loadGenome() {
  if (loading || busy) return;
  const value = $('genome-menu').value;
  if (!value) { show('genome-status', 'Choose a genome in the list.'); return; }
  $('genome-fungidb').hidden = true;
  if (value.startsWith('fungidb:')) {
    const name = value.slice('fungidb:'.length);
    show('genome-status', `${name} is in FungiDB.`);
    $('genome-command').textContent = `qmutfind fetch --organism "${name}" --out-dir genome --eupathdb-key-file ~/.eupathdb_keys`;
    $('genome-fungidb').hidden = false;
    return;
  }
  const entry = genomeList.find((g) => `ncbi:${g.accession}` === value);
  if (!entry) return;
  loading = true;
  $('genome-load').disabled = true;
  try {
    const base = genomesBase;
    const reference = await download(entry.fasta, base, 'the genome', 'genome-status');
    const gff = entry.gff ? await download(entry.gff, base, 'the gene models', 'genome-status') : null;
    // A new genome replaces the example, whose reads belong to another genome.
    example.reference = reference;
    example.gff = gff;
    example.reads = [];
    example.expectedReads = null;
    show('example-status', '');
    example.source = { reference: `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.fasta, base)})`, gff: entry.gff ? `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.gff, base)})` : '', reads: '' };
    describeSources();
    let note = entry.gff ? '' : ' NCBI has no gene models for this genome, so step 3 cannot annotate it unless you choose a GFF3 file.';
    note += ' Choose your reads in step 1.';
    if (entry.genome_size > LARGE_GENOME) {
      note += ` This genome is large (${(entry.genome_size / 1e6).toFixed(0)} Mb). Genomes above 50 Mb have not been tested in a browser and may need several GB of memory.`;
    }
    show('genome-status', `Loaded: ${entry.label} (${entry.accession}).${note}`);
  } catch (error) {
    show('genome-status', `Error: ${error.message} The server may block downloads from web pages, or the network failed.`);
  } finally {
    loading = false;
    $('genome-load').disabled = false;
  }
}

async function setupGenomes() {
  const asked = new URLSearchParams(location.search).get('genomes');
  const path = asked && !/^([a-z]+:)?\/\//i.test(asked) && !asked.startsWith('/') ? asked : 'genomes/genomes.json';
  const base = new URL(path, location.href);
  genomesBase = base;
  try {
    const response = await fetch(base, { cache: 'no-cache' });
    if (!response.ok) return;
    genomeList = (await response.json()).genomes || [];
  } catch (error) {
    return;
  }
  try {
    const response = await fetch(new URL('fungidb.json', base), { cache: 'no-cache' });
    if (response.ok) {
      const document_ = await response.json();
      fungidbList = document_.organisms || [];
      $('genome-text').textContent = `${genomeList.length} NCBI genomes (RefSeq, and two GenBank) download in this page. ` +
        `${fungidbList.length} FungiDB organisms (${document_.release || 'FungiDB'}) are listed too; they need the command line.`;
    }
  } catch (error) {
    // The FungiDB list is optional.
  }
  if (!$('genome-text').textContent) $('genome-text').textContent = `${genomeList.length} NCBI genomes download in this page.`;
  fillMenu();
  $('genome-filter').addEventListener('input', fillMenu);
  $('genome-load').addEventListener('click', loadGenome);
  $('genome-menu').addEventListener('dblclick', loadGenome);
  $('genome-section').hidden = false;
}

setupGenomes();


// ---- The locus drawing in its iframe ------------------------------------------------------------------------------
// The drawing is a separate page with its own colours. After it loads, and when the theme changes, the page gives it the
// colours of the panel around it, and a palette with enough contrast for that kind of background (light, dark or mid).
const PALETTES = {
  light: { missense: '#B36B00', nonsense: '#B83A00', frameshift: '#B0457F', inframe: '#007A5A', splice: '#0072B2', utr: '#1F78B4', synonymous: '#707070', intergenic: '#8A7B00', unannotated: '#707070', cds: '#3d5a80', utrfill: '#7a8aa0' },
  dark: { missense: '#F5B942', nonsense: '#FF7F50', frameshift: '#E08FC0', inframe: '#4FD1A5', splice: '#56B4E9', utr: '#6FC3F0', synonymous: '#B0B0B0', intergenic: '#F0E442', unannotated: '#BBBBBB', cds: '#8fb0d8', utrfill: '#7C93B0' },
  mid: { missense: '#6B3F00', nonsense: '#8A2B00', frameshift: '#7A1F55', inframe: '#00573F', splice: '#00467A', utr: '#1A4F7A', synonymous: '#444444', intergenic: '#5F5500', unannotated: '#444444', cds: '#1B2F4A', utrfill: '#3B4E66' },
};
const EFFECT_GROUPS = {
  missense: ['missense'], nonsense: ['nonsense', 'start_lost', 'stop_lost'], frameshift: ['frameshift', 'large_indel_in_cds'],
  inframe: ['inframe_insertion', 'inframe_deletion'], splice: ['splice_site'], utr: ['5UTR_derived', '3UTR_derived', 'intron_derived', 'non_coding'],
  synonymous: ['synonymous'], intergenic: ['intergenic'], unannotated: ['unannotated'],
};

function luminance(rgb) {
  const [r, g, b] = rgb.match(/[\d.]+/g).slice(0, 3).map((v) => { v /= 255; return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4; });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

function themeReport() {
  const doc = $('report').contentDocument;
  if (!doc || !doc.documentElement || !$('report').getAttribute('srcdoc')) return;
  const panel = getComputedStyle($('step3'));
  const bg = panel.backgroundColor;
  const fg = panel.color;
  const l = luminance(bg);
  const p = PALETTES[l > 0.5 ? 'light' : l < 0.2 ? 'dark' : 'mid'];
  let css = `:root{--bg:${bg};--fg:${fg};--line:${fg};--cds:${p.cds};--utr:${p.utrfill}}body{background:${bg};color:${fg}}` +
    `.variant,.variant-legend{stroke:${fg};stroke-width:1}.cds{fill:${p.cds}}.utr{fill:${p.utrfill}}`;
  for (const [key, labels] of Object.entries(EFFECT_GROUPS)) {
    css += labels.map((e) => `[data-effect="${e}"]`).join(',') + `{fill:${p[key]}}`;
  }
  let style = doc.getElementById('qmutfind-theme');
  if (!style) { style = doc.createElement('style'); style.id = 'qmutfind-theme'; doc.head.appendChild(style); }
  style.textContent = css;
}

$('report').addEventListener('load', themeReport);
document.addEventListener('qmutfind-theme', themeReport);


// ---- Load data from the web: an SRA run (looked up in ENA) or https addresses -----------------------------------------
// ENA's lookup service and ENA's file server allow browser downloads (they send Access-Control-Allow-Origin: *).
// ?enaapi=PATH/ uses another lookup service, but only a path on this site (the tests use that).
const ENA_FIELDS = 'run_accession,fastq_ftp,fastq_bytes,read_count,library_layout,instrument_platform,instrument_model,scientific_name,sample_title';
const KIND_OF_PLATFORM = { ILLUMINA: 'short', OXFORD_NANOPORE: 'long' };
let pendingRun = null;
let webAbort = null;

function enaBase() {
  const asked = new URLSearchParams(location.search).get('enaapi');
  if (asked && !/^([a-z]+:)?\/\//i.test(asked) && !asked.startsWith('/')) return new URL(asked, location.href);
  return new URL('https://www.ebi.ac.uk/ena/portal/api/');
}

async function sraLookup() {
  pendingRun = null;
  $('sra-load').hidden = true;
  const accession = $('sra-accession').value.trim().toUpperCase();
  const say = (text) => show('sra-info', text);
  if (/^(PRJ|SAM|SRP|SRS|SRX|ERP|ERS|ERX|DRP|DRS|DRX)/.test(accession)) {
    say(`${accession} is a project, sample or experiment accession. Give a run accession (SRR, ERR or DRR and digits); the ENA page of the project lists its runs.`);
    return;
  }
  if (!/^[SED]RR\d{6,}$/.test(accession)) {
    say(`"${$('sra-accession').value.trim()}" does not look like a run accession (SRR, ERR or DRR followed by at least 6 digits).`);
    return;
  }
  say('Looking the run up in ENA');
  let records;
  try {
    const url = new URL(`filereport?accession=${accession}&result=read_run&fields=${ENA_FIELDS}&format=json`, enaBase());
    const response = await fetch(url, { cache: 'no-cache' });
    if (!response.ok) throw new Error(`ENA answered ${response.status}`);
    records = await response.json();
  } catch (error) {
    say(`The ENA lookup failed: ${error.message}. The network may be down or ENA may be busy.`);
    return;
  }
  const record = records && records[0];
  if (!record) { say(`${accession} was not found in ENA.`); return; }
  const addresses = (record.fastq_ftp || '').split(';').filter(Boolean);
  const sizes = (record.fastq_bytes || '').split(';').filter(Boolean).map(Number);
  if (!addresses.length) {
    say(`ENA has no FASTQ files for ${accession} (it may hold only BAM or SRA-format files). The page cannot use this run.`);
    return;
  }
  const platform = record.instrument_platform;
  if (!KIND_OF_PLATFORM[platform]) {
    say(`${accession} was made on ${platform || 'an unknown platform'}, which is not supported (only ILLUMINA and OXFORD_NANOPORE reads).`);
    return;
  }
  const total = sizes.reduce((a, b) => a + b, 0);
  if (total > MAX_DOWNLOAD) {
    say(`${accession} is too large for this page (${mb(total)}; the limit is ${mb(MAX_DOWNLOAD)}). Use the command-line tool.`);
    return;
  }
  const paired = record.library_layout === 'PAIRED';
  const count = Number(record.read_count) || 0;
  const lines = [
    `${record.scientific_name || 'Unknown organism'}${record.sample_title ? ` (${record.sample_title})` : ''}`,
    `${platform}${record.instrument_model ? `, ${record.instrument_model}` : ''}; library ${record.library_layout}; ${paired ? `${count.toLocaleString()} read pairs` : `${count.toLocaleString()} reads`} in ${addresses.length} file(s), ${mb(total)}.`,
  ];
  if (paired) lines.push('The page maps each mate as a separate read (the command-line tool can map pairs).');
  if (total > WARN_DOWNLOAD) lines.push(`Large download: it can take minutes and the page holds all of it in memory.`);
  say(lines.join('\n'));
  pendingRun = {
    accession, platform, count: paired ? count * 2 : count,
    files: addresses.map((a, i) => ({ name: a.split('/').pop(), url: /^https?:/.test(a) ? a : `https://${a}`, bytes: sizes[i] || 0 })),
  };
  $('sra-load').hidden = false;
}

function webBusy(on) {
  $('sra-load').disabled = $('url-load').disabled = $('sra-lookup').disabled = on;
  $('web-cancel').hidden = !on;
}

async function runWebDownload(task) {
  if (loading || busy) return;
  loading = true;
  webBusy(true);
  webAbort = new AbortController();
  try {
    await task(webAbort.signal);
  } catch (error) {
    show('web-status', error.name === 'AbortError' ? 'Cancelled. Nothing was loaded.' : `Error: ${error.message}`);
  } finally {
    loading = false;
    webAbort = null;
    webBusy(false);
  }
}

$('sra-lookup').addEventListener('click', sraLookup);
$('sra-accession').addEventListener('keydown', (e) => { if (e.key === 'Enter') sraLookup(); });
$('web-cancel').addEventListener('click', () => { if (webAbort) webAbort.abort(); });

$('sra-load').addEventListener('click', () => runWebDownload(async (signal) => {
  const run = pendingRun;
  if (!run) return;
  const files = [];
  for (const [i, item] of run.files.entries()) files.push(await download(item, location.href, `${run.accession} file ${i + 1} of ${run.files.length}`, 'web-status', signal));
  example.reads = files;
  example.expectedReads = run.count;
  example.source.reads = `ENA run ${run.accession} (${new URL(run.files[0].url).host})`;
  $('kind').value = KIND_OF_PLATFORM[run.platform];
  describeSources();
  show('web-status', `Loaded: the reads of ${run.accession} (${mb(files.reduce((a, f) => a + f.size, 0))}). Press Run in step 1.`);
}));

$('url-load').addEventListener('click', () => runWebDownload(async (signal) => {
  const wanted = [
    ['reference', 'the reference genome', $('url-reference').value.trim() ? [$('url-reference').value.trim()] : []],
    ['gff', 'the gene file', $('url-gff').value.trim() ? [$('url-gff').value.trim()] : []],
    ['reads', 'the reads', $('url-reads').value.split('\n').map((x) => x.trim()).filter(Boolean)],
  ];
  if (wanted.every(([, , list]) => !list.length)) throw new Error('Type at least one https address.');
  let changed = false;
  for (const entry of wanted) {
    entry[2] = entry[2].map((u) => { const d = directAddress(u); changed ||= d !== u; return d; });
  }
  const urls = wanted.flatMap(([, label, list]) => list.map((u) => [label, checkUrl(u, label)])); // all checked before any download
  const got = { reference: null, gff: null, reads: [] };
  for (const [slot, label, list] of wanted) {
    for (const text of list) {
      const url = new URL(text);
      const name = decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'download');
      const file = await download({ name, url: url.href, bytes: 0 }, location.href, label, 'web-status', signal);
      if (slot === 'reads') got.reads.push(file); else got[slot] = file;
      example.source[slot] = `the address ${url.host}${url.pathname}`;
    }
  }
  void urls;
  if (got.reference) example.reference = got.reference;
  if (got.gff) example.gff = got.gff;
  if (got.reads.length) { example.reads = got.reads; example.expectedReads = null; }
  describeSources();
  show('web-status', `Loaded: ${wanted.filter(([, , l]) => l.length).map(([, label]) => label).join(', ')}.` +
    (changed ? ' A Dropbox share link was changed to dl.dropboxusercontent.com, the address that allows downloads from web pages.' : ''));
}));

describeSources();
