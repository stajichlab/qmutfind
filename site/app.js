// The page: three steps (call, summarize, draw). The result of each step is kept here and the next step uses it.
import { zipSync, unzipSync, strToU8, strFromU8 } from './vendor/fflate.js?v=__BUILD__';
import { storeAvailable, readStore, writeStore, clearStore, estimateStorage } from './persist.js?v=__BUILD__';

const $ = (id) => document.getElementById(id);
const worker = new Worker('worker.js?v=__BUILD__', { type: 'module' });

// The state that carries from step to step.
const state = { vcf: null, annotated: null, depth: null };
// The last combine of the samples: { table, summary, vcf, parent, minDepth, names }, or null. Declared here because describeSource() runs at load.
let lastCombined = null;
const samples = []; // the sample set, filled by the sample section below; declared here because describeSource() runs at load
const GENIC = '__genic__'; // the value of the "Genic" choice of the Show dropdown: any variant that overlaps a gene
const ALL_SAMPLES = '__all__'; // the value of the "All samples together" choice of the results panel
let busy = false;
const MAX_ROWS = 500;
// How far a file may unpack when it is loaded whole into memory. `?limit_mb=N` can lower it (the tests use that).
const LIMIT_BYTES = Math.min(parseFloat(new URLSearchParams(location.search).get('limit_mb')) || 1024, 1024) * 1e6;
// A download is limited to 2 GB. `?download_limit_mb=N` can lower the limit (the tests use that).
const MAX_DOWNLOAD = Math.min(parseFloat(new URLSearchParams(location.search).get('download_limit_mb')) || 2000, 2000) * 1e6;
const WARN_DOWNLOAD = 5e8;

/**
 * A page address option (?example=, ?genomes=, ?enaapi=) may only name a path on this site. The text is refused when it has a
 * scheme or starts with / or \ (browsers read \ as /, so \\host/x would name another host), and the resolved address must
 * have this page's origin. Returns the text, or `fallback`.
 */
function sitePath(asked, fallback) {
  if (!asked || asked !== asked.trim() || /[\u0000-\u001f]/.test(asked)) return fallback; // the URL parser would drop these characters
  if (/^[\\/]/.test(asked) || /^[a-z][a-z0-9+.-]*:/i.test(asked)) return fallback;
  try {
    return new URL(asked, location.href).origin === location.origin ? asked : fallback;
  } catch (error) {
    return fallback;
  }
}
window.qmutfindSitePath = sitePath;

const show = (id, text) => { $(id).textContent = text; };
const buttons = ['summarize', 'annotate', 'draw'];
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

const sampleName = () => (viewingAll() ? 'combined' : (selectedSample() && selectedSample().name) || 'sample');

// Files from an example. A file chosen in an input always wins over the example's file.
const example = { reference: null, gff: null, source: { reference: '', gff: '' } };

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
      el.textContent = 'Not set. Use the NCBI list, upload files or give addresses above.';
    }
  };
  line('reference', 'reference', example.reference ? [example.reference] : []);
  line('gff', 'gff', example.gff ? [example.gff] : []);
}
for (const id of ['reference', 'gff']) $(id).addEventListener('change', describeSources);
const hostOf = (item, base) => new URL(item.url, base).host;
const refFile = () => $('reference').files[0] || example.reference;
const gffFile = () => $('gff').files[0] || example.gff;

/** The VCF the next step should use: a chosen file, else the annotated VCF, else the step 1 result. */
function currentVcf() {
  const file = $('vcf-file').files[0];
  if (file) return file;
  return state.annotated || state.vcf;
}

/** The coverage the next step should use: a chosen file, else the BED of the chosen sample. */
function currentDepth() {
  return $('depth-file').files[0] || state.depth;
}

function describeSource() {
  const file = $('vcf-file').files[0];
  let text;
  if (file) text = `Using the VCF file ${file.name}.`;
  else if (state.annotated) text = 'Using the annotated VCF from the gene step (section 4).';
  else if (viewingAll()) text = `Using the merged VCF of all ${lastCombined.names.length} samples (one record per distinct variant).`;
  else if (state.vcf) text = `Using the VCF of the sample ${sampleName()}.`;
  else text = 'No VCF yet. Run a sample or choose a VCF file.';
  show('vcf-source', text);
  const showing = showingText();
  show('showing-summarize', showing);
  show('showing-draw', showing);
  updateAfNote();
}

/** What panels 3 and 4 show: the VCF file chosen, all samples together, one sample, or nothing yet. */
function showingText() {
  const file = $('vcf-file').files[0];
  if (file) return `Showing the VCF file ${file.name}. AF is the allele fraction recorded in that file.`;
  const without = subtracting() ? `, without the variants that the parent ${lastCombined.parent} has` : '';
  if (viewingAll()) {
    return `Showing all ${lastCombined.names.length} samples together${without}: ${countRecords(state.vcf).toLocaleString()} distinct variants. AF is the mean over the samples that call each variant.`;
  }
  const s = selectedSample();
  if (s) {
    const c = sharingCounts(state.vcf);
    const shared = c ? `, ${c.one.toLocaleString()} only in this sample and ${c.shared.toLocaleString()} shared with other samples` : '';
    return `Showing the sample ${s.name}${without}: ${countRecords(state.vcf).toLocaleString()} variants${shared}. AF is the allele fraction of this sample.`;
  }
  return 'Nothing to show yet. Run a sample (section 2), combine them, or choose a VCF file below.';
}

/** A short name of what is shown, for titles. */
function viewLabel() {
  const file = $('vcf-file').files[0];
  if (file) return `VCF file ${file.name}`;
  if (viewingAll()) return 'All samples together';
  const s = selectedSample();
  return s ? `Sample ${s.name}` : 'Variants';
}

/** The note under the Show dropdown says whose allele fraction (AF) and depth (DP) the table has. */
function updateAfNote() {
  let text = '';
  if ($('vcf-file').files[0]) text = 'AF is the allele fraction recorded in the VCF file.';
  else if (viewingAll()) text = 'AF is the mean over the samples that call the variant; "Called in" lists them. DP is the mean depth in the same way.';
  else if (selectedSample()) text = 'AF is the allele fraction of the sample shown (alternate reads divided by the depth at the site).';
  show('af-note', text);
}

/** The tooltips of the AF and DP column headers follow the view. */
function columnHint(label) {
  const all = viewingAll() && !$('vcf-file').files[0];
  if (label === 'AF') return all ? 'Mean allele fraction over the samples that call the variant' : 'Allele fraction of this sample (alternate reads divided by depth)';
  if (label === 'DP') return all ? 'Mean depth over the samples that call the variant' : 'Read depth at the site in this sample';
  return '';
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
  if (!vcf) show(statusId, 'There is no VCF yet. Run a sample, or choose a VCF file.');
  return vcf;
}

// Step 2
$('summarize').addEventListener('click', () => {
  const vcf = currentVcf();
  const bed = currentDepth();
  if (!vcf && !bed) { show('summarize-status', 'There is no VCF or coverage yet. Run a sample, or choose a VCF or coverage file.'); return; }
  send('summarize', { vcf, bed, reference: refFile() || null, order: $('order').value }, 'summarize-status');
});

const escapeHtml = (text) => String(text).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

/** The statistics of the merged VCF of all samples: the usual rows, how the variants are shared, and a table per sample. */
function sharingTables(stats) {
  const lines = lastCombined.vcf.split('\n').filter((l) => l && l[0] !== '#');
  const n = lastCombined.names.length;
  let every = 0;
  let one = 0;
  let isNew = 0;
  for (const line of lines) {
    const info = Object.fromEntries(line.split('\t')[7].split(';').map((kv) => { const i = kv.indexOf('='); return i > 0 ? [kv.slice(0, i), kv.slice(i + 1)] : [kv, '']; }));
    const ns = Number(info.NS);
    if (ns === n) every += 1;
    if (ns === 1) one += 1;
    if (info.NEW) isNew += 1;
  }
  const extra = [['Samples combined', String(n)], ['Called in every sample', String(every)], ['Called in one sample only', String(one)]];
  if (lastCombined.parent) extra.push([`New compared with the parent (${lastCombined.parent})`, String(isNew)]);
  const rows = lastCombined.summary.split('\n').filter(Boolean).map((l) => l.split('\t'));
  const head = rows[0].map((h) => `<th>${escapeHtml(h)}</th>`).join('');
  const body = rows.slice(1).map((r) => `<tr>${r.map((c) => `<td>${escapeHtml(c)}</td>`).join('')}</tr>`).join('');
  return `${statsTable(stats, extra)}<table class="stats per-sample"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
}

/** The rows of the summary of one sample that say what is only in it and what is shared (empty with one finished sample). */
function singleSampleSharingRows() {
  const c = sharingCounts(state.vcf);
  if (!c || viewingAll() || $('vcf-file').files[0]) return [];
  return [['Only in this sample', String(c.one)], ['Shared with other samples', String(c.shared)], ['Shared by every sample', String(c.every)]];
}

function statsTable(s, extra = []) {
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
    add('Effects', effects.sort((a, b) => b[1] - a[1]).map(([k, v]) => `${escapeHtml(k)}: ${n(v)}`).join(', '));
  }
  for (const [name, value] of extra) add(escapeHtml(name), escapeHtml(value));
  return `<table class="stats">${rows.join('')}</table>`;
}

// Step 3
$('annotate').addEventListener('click', () => {
  const vcf = needVcf('annotate-status');
  if (!vcf) return;
  const reference = refFile();
  const gff = gffFile();
  if (!reference || !gff) {
    show('annotate-status', 'Choose the reference genome and the GFF3 file in section 1.');
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
      ns: info.NS === undefined ? undefined : Number(info.NS),
      afs: info.AFS ? info.AFS.split(',') : [], // the AF of each sample in `samples`, in the same order
      samples: info.SAMPLES ? info.SAMPLES.split(',').map(decodeInfo) : [],
      newIn: info.NEW ? info.NEW.split(',').map(decodeInfo) : [],
      gene: info.GENE || '', geneName: decodeInfo(info.GENENAME), product: decodeInfo(info.PRODUCT), aa: info.AA && info.AA !== '.' ? info.AA : '',
      // trunc is for a frameshift with a new stop: the new and the reference number of residues
      trunc: info.TRUNC && info.TRUNC !== '.' ? info.TRUNC.split(',')[0] : '', af: info.AF || '', dp: info.DP || '',
    });
  }
  return rows;
}

let variants = [];
const cell = (text) => { const td = document.createElement('td'); td.textContent = text; return td; };

// The columns of the variant table. A header is a button: the first click sorts up, the second sorts down, the third
// restores the order of the file. Empty cells stay last in both directions. All variants are sorted, not only the rows shown.
const short = (text) => (text.length > 12 ? `${text.slice(0, 12)}\u2026` : text);
const BASE_COLUMNS = [
  { label: 'Position', kind: 'position', text: (v) => `${v.contig}:${v.pos.toLocaleString()}` },
  { label: 'Change', kind: 'text', value: (v) => `${v.ref} > ${v.alt}`, text: (v) => `${short(v.ref)} > ${short(v.alt)}` },
  { label: 'Class', kind: 'text', value: (v) => v.cls },
  { label: 'Effect', kind: 'text', value: (v) => v.effect },
  { label: 'Gene', kind: 'text', value: (v) => v.gene },
  { label: 'Name', kind: 'text', value: (v) => v.geneName },
  { label: 'Protein change', kind: 'text', value: (v) => v.aa, text: (v) => (v.trunc && v.aa ? `${v.aa} (protein ${v.trunc.replace('/', ' of ')} aa)` : v.aa) },
  { label: 'Product', kind: 'text', value: (v) => v.product },
  { label: 'AF', kind: 'number', value: (v) => v.af },
  { label: 'DP', kind: 'number', value: (v) => v.dp },
];
// The merged VCF of all samples adds two columns: how many samples carry the variant, and which.
const COMBINED_COLUMNS = [
  { label: 'Samples', kind: 'number', value: (v) => v.ns, text: (v) => (variantMeta.samples ? `${v.ns} of ${variantMeta.samples}` : String(v.ns)) },
  // The names, each with the AF of that sample: the AF column is only their mean.
  { label: 'Called in', kind: 'text', value: (v) => v.samples.join(', '), text: (v) => v.samples.map((n, i) => (Number.isFinite(parseFloat(v.afs[i])) ? `${n} (${parseFloat(v.afs[i]).toFixed(2)})` : n)).join(', ') },
];
const isCombinedTable = () => variants.some((v) => v.ns !== undefined);
const activeColumns = () => (isCombinedTable() ? [...BASE_COLUMNS, ...COMBINED_COLUMNS] : BASE_COLUMNS);
let variantMeta = { samples: 0, parent: '' }; // from the header of the merged VCF
let sortColumn = null; // an index into TABLE_COLUMNS, or null for the order of the file
let sortDirection = 1; // 1 up, -1 down
let contigRank = new Map(); // the order in which the contigs first appear in the file

/** The variants in the order of the chosen column. Array.prototype.sort is stable, so equal cells keep the file order. */
function sortedVariants(list) {
  if (sortColumn === null) return list;
  const column = activeColumns()[sortColumn];
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
  const empty = (v) => (column.kind === 'position' ? false : column.kind === 'number' ? !Number.isFinite(parseFloat(column.value(v))) : !column.value(v));
  const compare = (a, b) => {
    if (column.kind === 'position') return (contigRank.get(a.contig) - contigRank.get(b.contig)) || (a.pos - b.pos);
    if (column.kind === 'number') return parseFloat(column.value(a)) - parseFloat(column.value(b));
    return collator.compare(column.value(a), column.value(b));
  };
  return [...list].sort((a, b) => {
    const ea = empty(a);
    const eb = empty(b);
    if (ea || eb) return ea === eb ? 0 : (ea ? 1 : -1); // empty cells last, whatever the direction
    return sortDirection * compare(a, b);
  });
}

function sortBy(index) {
  if (sortColumn === index) {
    if (sortDirection === 1) sortDirection = -1;
    else sortColumn = null;
  } else {
    sortColumn = index;
    sortDirection = 1;
  }
  renderTable();
}

function renderTable() {
  const filter = $('effect-filter').value;
  const sharing = $('sharing-wrap').hidden ? '' : $('sharing-filter').value;
  const shared = (v) => !sharing || (sharing === 'one' && v.ns === 1) || (sharing === 'shared' && v.ns >= 2) || (sharing === 'all' && v.ns === variantMeta.samples) || (sharing === 'new' && v.newIn.length > 0);
  // The Genic choice is any variant that overlaps a gene, whatever its effect; the other choices are one effect each.
  const wanted = (v) => !filter || (filter === GENIC ? v.gene !== '' : v.effect === filter);
  const shown = sortedVariants(variants.filter((v) => wanted(v) && shared(v)));
  const columns = activeColumns();
  const table = $('variant-table');
  table.replaceChildren();
  const head = table.createTHead().insertRow();
  for (const [index, column] of columns.entries()) {
    const th = document.createElement('th');
    th.setAttribute('aria-sort', sortColumn === index ? (sortDirection === 1 ? 'ascending' : 'descending') : 'none');
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'sort';
    button.textContent = column.label;
    if (columnHint(column.label)) button.title = columnHint(column.label);
    button.addEventListener('click', () => sortBy(index));
    th.appendChild(button);
    head.appendChild(th);
  }
  head.appendChild(document.createElement('th')); // the Draw buttons
  const body = table.createTBody();
  for (const v of shown.slice(0, MAX_ROWS)) {
    const row = body.insertRow();
    row.append(...columns.map((c) => cell(c.text ? c.text(v) : String(c.value(v) ?? ''))));
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
$('sharing-filter').addEventListener('change', renderTable);
window.qmutfindLoadVariants = (vcfText) => loadVariants(vcfText); // for the tests

function loadVariants(vcfText) {
  variants = parseVariants(vcfText);
  sortColumn = null; // a new set of variants starts in the order of the file
  variantMeta = { samples: Number((/^##qmutfind_samples=(\d+)/m.exec(vcfText) || [])[1]) || 0, parent: decodeInfo(((/^##qmutfind_parent=(.*)$/m.exec(vcfText) || [])[1]) || '') };
  const sharingSelect = $('sharing-filter');
  const oneSample = !viewingAll() && !$('vcf-file').files[0] && Boolean(selectedSample()); // one sample among several finished ones
  sharingSelect.replaceChildren(
    new Option('All variants', ''),
    new Option(oneSample ? 'Only in this sample' : 'Called in one sample only', 'one'),
    new Option(oneSample ? 'Shared with at least one other sample' : 'Called in two or more samples', 'shared'),
    new Option(oneSample ? 'Shared by every sample' : 'Called in every sample', 'all'),
  );
  if (variantMeta.parent) sharingSelect.add(new Option(`New compared with the parent (${variantMeta.parent})`, 'new'));
  $('sharing-wrap').hidden = !variants.some((v) => v.ns !== undefined);
  contigRank = new Map();
  for (const v of variants) if (!contigRank.has(v.contig)) contigRank.set(v.contig, contigRank.size);
  const counts = {};
  for (const v of variants) counts[v.effect] = (counts[v.effect] || 0) + 1;
  const select = $('effect-filter');
  select.disabled = false;
  const inGenes = variants.filter((v) => v.gene !== '').length;
  select.replaceChildren(new Option(`All variants (${variants.length})`, ''), new Option(`Genic: in a gene (${inGenes})`, GENIC));
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
  send('report', { vcf, gff: gffFile() || null, region, title: `${viewLabel()}: mutations in ${region}` }, 'report-status');
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
  const panel = { summarize: 'summarize-status', annotate: 'annotate-status', report: 'report-status' }[m.step];
  if (m.type === 'status') show(panel, m.text);
  if (m.type === 'error') {
    show(panel, `Error: ${m.message}`);
    setBusy(false);
    return;
  }
  if (m.type !== 'done') return;
  setBusy(false);
  if (m.step === 'summarize') {
    const parts = [];
    if (m.stats) {
      const s = JSON.parse(m.stats);
      $('stats-table').innerHTML = viewingAll() ? sharingTables(s) : statsTable(s, singleSampleSharingRows());
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

/** `budget` is { left }: the bytes that the downloads of one sample may still take. Its files together stay within the download limit. */
async function download(item, base, label, statusId = 'example-status', signal = undefined, budget = { left: MAX_DOWNLOAD }) {
  const url = checkUrl(new URL(item.url, base).href, label);
  const response = await fetch(url, { signal });
  if (!response.ok) throw new Error(`${label}: the server answered ${response.status} for ${url.host}`);
  const announced = Number(response.headers.get('Content-Length')) || 0;
  const expected = item.bytes || announced;
  const tooLarge = (size, more = '') => new Error(`${label}: the download is too large${budget.left < MAX_DOWNLOAD ? ' together with the other files of this sample' : ''} (${more}${mb(size)}; the limit is ${mb(budget.left < MAX_DOWNLOAD ? budget.left : MAX_DOWNLOAD)}${budget.left < MAX_DOWNLOAD ? ' for what is left' : ''}).`);
  if (expected > budget.left) throw tooLarge(expected);
  const reader = response.body.getReader();
  const parts = [];
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    parts.push(value);
    got += value.length;
    if (got > budget.left) { await reader.cancel(); throw tooLarge(budget.left, 'more than '); }
    show(statusId, `Downloading ${label}: ${mb(got)}${expected ? ` of ${mb(expected)}` : ''}`);
  }
  if (item.bytes && got !== item.bytes) throw new Error(`${label}: got ${got} bytes, expected ${item.bytes}. The download is incomplete.`);
  budget.left -= got;
  // The page reads a file as gzip by its name. A download whose name has no .gz but whose content is gzip gets the name.
  let name = item.name;
  if (parts.length && parts[0].length > 1 && parts[0][0] === 0x1f && parts[0][1] === 0x8b && !name.endsWith('.gz')) name += '.gz';
  return new File(parts, name);
}

/**
 * Download a FASTQ file (plain or gzip) and keep only its first `maxReads` reads, then stop the download. A worker does the
 * unpacking and packing (see `fetchReads` in worker.js), so the page stays responsive. Returns { file, reads, cut }: `cut` is true
 * when the download stopped at the limit. `budget` is { left }, as for `download`.
 */
async function downloadFirstReads(item, base, label, statusId, signal, maxReads, budget = { left: MAX_DOWNLOAD }) {
  const url = checkUrl(new URL(item.url, base).href, label);
  const m = await callWorker({ step: 'fetchreads', url: url.href, maxReads, maxKept: budget.left, maxDownload: budget.left, label },
    (p) => show(statusId, `Reading ${label}: ${mb(p.downloaded)} downloaded, ${p.reads.toLocaleString()} reads kept`));
  budget.left -= m.bytes.length;
  return { file: new File([m.bytes], `${item.name.replace(/\.gz$/, '')}.gz`), reads: m.reads, cut: m.cut };
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
    const reads = [];
    for (const [i, item] of chosen.reads.entries()) {
      reads.push(await download(item, base, `reads file ${i + 1} of ${chosen.reads.length}`));
    }
    const origin = (item) => `the example (${item.source ? item.source + ' at ' : ''}${hostOf(item, base)})`;
    example.source = { reference: origin(manifest.reference), gff: origin(manifest.gff), reads: origin(chosen.reads[0]) };
    describeSources();
    if (manifest.kind) $('kind').value = manifest.kind;
    // The example's reads become a sample of the set. An earlier example sample is replaced.
    dropExampleSamples();
    const blank = samples.length === 1 && isBlank(samples[0]) ? samples[0] : null;
    const added = blank || newSample();
    Object.assign(added, { name: manifest.sample || 'example', source: 'local', files: reads, fromExample: true, expected: chosen.expected_reads || null });
    renderSamples();
    show('example-status', `Loaded: ${chosen.label}. Press "Run the waiting samples" in section 2.`);
  } catch (error) {
    example.reference = example.gff = null;
    describeSources();
    show('example-status', `Error: ${error.message} The server may block downloads from web pages, or the network failed.`);
  } finally {
    loading = false;
    for (const b of document.querySelectorAll('#example-buttons button')) b.disabled = false;
  }
}

async function setupExamples() {
  const asked = new URLSearchParams(location.search).get('example');
  const path = sitePath(asked, 'example/manifest.json');
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
    show('example-status', '');
    dropExampleSamples();
    example.source = { reference: `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.fasta, base)})`, gff: entry.gff ? `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.gff, base)})` : '', reads: '' };
    describeSources();
    let note = entry.gff ? '' : ' NCBI has no gene models for this genome, so the gene step (section 4) cannot annotate it unless you upload a GFF3 file.';
    note += ' Add your samples in section 2.';
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
  const path = sitePath(asked, 'genomes/genomes.json');
  const base = new URL(path, location.href);
  genomesBase = base;
  try {
    const response = await fetch(base, { cache: 'no-cache' });
    if (!response.ok) throw new Error(`the server answered ${response.status}`);
    genomeList = (await response.json()).genomes || [];
  } catch (error) {
    $('genome-text').textContent = `The NCBI genome list could not be loaded (${error.message}). Upload files or give addresses instead.`;
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
}

setupGenomes();

// The reference panel asks one question and shows the controls of the answer: the NCBI list, files, or https addresses.
function applyReferenceMode() {
  const mode = $('reference-source').value;
  $('genome-section').hidden = mode !== 'ncbi';
  $('reference-upload').hidden = mode !== 'upload';
  $('reference-url').hidden = mode !== 'url';
}
$('reference-source').addEventListener('change', applyReferenceMode);
applyReferenceMode();

$('ref-url-load').addEventListener('click', async () => {
  if (loading || busy || sampleRunning) return;
  const fastaText = $('url-reference').value.trim();
  const gffText = $('url-gff').value.trim();
  if (!fastaText) { show('genome-status', 'Type the address of the reference FASTA.'); return; }
  loading = true;
  $('ref-url-load').disabled = true;
  try {
    const wanted = [['reference', 'the reference genome', fastaText], ...(gffText ? [['gff', 'the gene file', gffText]] : [])];
    let changed = false;
    const urls = wanted.map(([, label, text]) => { const d = directAddress(text); changed ||= d !== text; return checkUrl(d, label); }); // all checked before any download
    const got = {};
    for (const [i, [slot, label]] of wanted.entries()) {
      const url = urls[i];
      const name = decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'download');
      got[slot] = await download({ name, url: url.href, bytes: 0 }, location.href, label, 'genome-status');
      example.source[slot] = `the address ${url.host}${url.pathname}`;
    }
    example.reference = got.reference;
    example.gff = got.gff || null;
    if (!got.gff) example.source.gff = '';
    show('example-status', '');
    dropExampleSamples();
    describeSources();
    show('genome-status', `Loaded: the reference genome${got.gff ? ' and the gene file' : ''}.${changed ? ' A Dropbox share link was changed to dl.dropboxusercontent.com, the address that allows downloads from web pages.' : ''}`);
  } catch (error) {
    show('genome-status', `Error: ${error.message}`);
  } finally {
    loading = false;
    $('ref-url-load').disabled = false;
  }
});


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


// ---- The sample set ---------------------------------------------------------------------------------------------
// A sample is a name and one or more reads files, from one source: an SRA run (looked up in ENA), https addresses, or files from
// this computer. The queue maps one sample at a time, each in its own worker, which is stopped afterwards, so the memory of one
// sample is given back before the next starts. A finished sample keeps its VCF, coverage BED and evidence VCF (text).
// `combine_tables` (WebAssembly) joins the finished samples; the catalog is exported as a zip and can be imported again.
// ENA's lookup service and file server allow browser downloads (they send Access-Control-Allow-Origin: *).
// ?enaapi=PATH/ uses another lookup service, but only a path on this site (the tests use that).
const ENA_FIELDS = 'run_accession,fastq_ftp,fastq_bytes,read_count,library_layout,instrument_platform,instrument_model,scientific_name,sample_title';
const KIND_OF_PLATFORM = { ILLUMINA: 'short', OXFORD_NANOPORE: 'long' };
const DEPTH_WINDOW = 1000; // bases per coverage window (the results panel and the combined table use these windows)
const CATALOG_FORMAT = 'qmutfind-catalog';
const CATALOG_VERSION = 1;
// The zip and the total of its files. `?catalog_limit_mb=N` can lower the limit (the tests use that).
const MAX_CATALOG_BYTES = Math.min(parseFloat(new URLSearchParams(location.search).get('catalog_limit_mb')) || 500, 500) * 1e6;

function enaBase() {
  const asked = new URLSearchParams(location.search).get('enaapi');
  if (sitePath(asked, '')) return new URL(asked, location.href);
  return new URL('https://www.ebi.ac.uk/ena/portal/api/');
}

/**
 * Look a run up in ENA. Returns { accession, platform, kind, files: [{ name, url, bytes }], count, perFile, paired, info }
 * (`count` is in reads, each mate counted; `perFile` is the reads or pairs in one file) or throws an Error that says why not.
 * With a read limit (`limit` above 0) a run above the download limit is allowed, because the download stops early.
 */
async function lookupRun(text, limit) {
  const accession = String(text).trim().toUpperCase();
  if (/^(PRJ|SAM|SRP|SRS|SRX|ERP|ERS|ERX|DRP|DRS|DRX)/.test(accession)) {
    throw new Error(`${accession} is a project, sample or experiment accession. Give a run accession (SRR, ERR or DRR and digits); the ENA page of the project lists its runs.`);
  }
  if (!/^[SED]RR\d{6,}$/.test(accession)) {
    throw new Error(`"${String(text).trim()}" does not look like a run accession (SRR, ERR or DRR followed by at least 6 digits).`);
  }
  let records;
  try {
    const url = new URL(`filereport?accession=${accession}&result=read_run&fields=${ENA_FIELDS}&format=json`, enaBase());
    const response = await fetch(url, { cache: 'no-cache' });
    if (!response.ok) throw new Error(`ENA answered ${response.status}`);
    records = await response.json();
  } catch (error) {
    throw new Error(`The ENA lookup failed: ${error.message}. The network may be down or ENA may be busy.`);
  }
  const record = records && records[0];
  if (!record) throw new Error(`${accession} was not found in ENA.`);
  const addresses = (record.fastq_ftp || '').split(';').filter(Boolean);
  const sizes = (record.fastq_bytes || '').split(';').filter(Boolean).map(Number);
  if (!addresses.length) throw new Error(`ENA has no FASTQ files for ${accession} (it may hold only BAM or SRA-format files). The page cannot use this run.`);
  const platform = record.instrument_platform;
  if (!KIND_OF_PLATFORM[platform]) throw new Error(`${accession} was made on ${platform || 'an unknown platform'}, which is not supported (only ILLUMINA and OXFORD_NANOPORE reads).`);
  const total = sizes.reduce((a, b) => a + b, 0);
  if (total > MAX_DOWNLOAD && !limit) {
    throw new Error(`${accession} is too large for this page (${mb(total)}; the limit is ${mb(MAX_DOWNLOAD)}). Use the command-line tool, or set a read limit.`);
  }
  const paired = record.library_layout === 'PAIRED';
  const count = Number(record.read_count) || 0;
  const lines = [
    `${record.scientific_name || 'Unknown organism'}${record.sample_title ? ` (${record.sample_title})` : ''}`,
    `${platform}${record.instrument_model ? `, ${record.instrument_model}` : ''}; library ${record.library_layout}; ${paired ? `${count.toLocaleString()} read pairs` : `${count.toLocaleString()} reads`} in ${addresses.length} file(s), ${mb(total)}.`,
  ];
  if (paired) lines.push('The page maps each mate as a separate read (the command-line tool can map pairs).');
  if (limit) lines.push(`Only the first ${limit.toLocaleString()} reads of each file will be kept; the download stops there.`);
  else if (total > WARN_DOWNLOAD) lines.push('Large download: it can take minutes and the page holds all of it in memory.');
  return {
    accession, platform, kind: KIND_OF_PLATFORM[platform], count: paired ? count * 2 : count, perFile: count, paired,
    files: addresses.map((a, i) => ({ name: a.split('/').pop(), url: /^https?:/.test(a) ? a : `https://${a}`, bytes: sizes[i] || 0 })),
    info: lines.join('\n'),
  };
}

// A sample: { id, name, note, source, accession, urls, files, fromExample, expected, state, info, reads, calls, seconds, vcf, bed, evidence, version }
let nextSampleId = 1;
let sampleRunning = false;
let sampleCancelled = false;
let activeWorker = null;
let activeAbort = null;

function newSample(init = {}) {
  const s = { id: nextSampleId++, name: '', note: '', source: 'local', accession: '', urls: '', files: [], state: 'waiting', ...init };
  samples.push(s);
  return s;
}
const sampleByName = (name) => samples.find((s) => s.name === name);
function isFinished(s) { return Boolean(s.vcf); }
const isBlank = (s) => !s.name && !s.note && !s.accession.trim() && !s.urls.trim() && s.files.length === 0 && !isFinished(s);
const sampleStatus = (text) => show('strain-status', text);

/** The number in the read-limit box is in thousands of reads per file. Returns reads per file; 0 means no limit. */
function readLimit() {
  const text = $('read-limit').value.trim();
  if (text === '') return 0;
  const k = Number(text);
  if (!Number.isFinite(k) || k < 0) throw new Error('The read limit must be a number of thousands of reads per file, 0 or more (0 uses all reads).');
  return Math.round(k * 1000);
}

function describeLimit() {
  let text;
  try {
    const n = readLimit();
    text = n ? `= ${n.toLocaleString()} reads from each file` : 'All reads of each file are used.';
  } catch (error) {
    text = error.message;
  }
  show('read-limit-echo', text);
}
$('read-limit').addEventListener('input', describeLimit);
describeLimit();

/** Run one message in a worker of its own. The worker is stopped when the answer comes, or when the queue is cancelled. */
function callWorker(message, onProgress) {
  return new Promise((resolve, reject) => {
    const w = new Worker('worker.js?v=__BUILD__', { type: 'module' });
    activeWorker = { w, reject };
    const finish = () => { w.terminate(); if (activeWorker && activeWorker.w === w) activeWorker = null; };
    w.onmessage = (event) => {
      const m = event.data;
      if (m.type === 'progress' && onProgress) onProgress(m);
      if (m.type === 'error') { finish(); reject(new Error(m.message)); }
      if (m.type === 'done') { finish(); resolve(m); }
    };
    w.onerror = (event) => { finish(); reject(new Error(event.message || 'the worker stopped')); };
    w.postMessage({ limit: LIMIT_BYTES, ...message });
  });
}

// The rows of the sample table
function refreshChoices() {
  const named = samples.filter((s) => s.name);
  const parent = $('strain-parent');
  const keepParent = parent.value;
  parent.replaceChildren(new Option('No parent', ''));
  for (const s of named) parent.appendChild(new Option(s.name, s.name));
  parent.value = named.some((s) => s.name === keepParent) ? keepParent : '';
  const result = $('result-sample');
  const keepResult = result.value;
  result.replaceChildren();
  const done = samples.filter((s) => s.name && isFinished(s));
  if (!done.length) result.appendChild(new Option('No finished sample yet', ''));
  // After a combine, the merged VCF of all the samples can be summarized and annotated like a sample.
  if (lastCombined && lastCombined.vcf) result.appendChild(new Option('All samples together', ALL_SAMPLES));
  for (const s of done) result.appendChild(new Option(s.name, s.name));
  const stillThere = done.some((s) => s.name === keepResult) || (keepResult === ALL_SAMPLES && lastCombined && lastCombined.vcf);
  result.value = stillThere ? keepResult : (done[0] ? done[0].name : '');
  if (result.value !== keepResult) useSelectedSample();
}

function viewingAll() { // a function declaration: describeSource() calls it at load
  return $('result-sample').value === ALL_SAMPLES && Boolean(lastCombined && lastCombined.vcf);
}
function selectedSample() { return samples.find((s) => s.name && s.name === $('result-sample').value && isFinished(s)) || null; }

/** The results panel follows the chosen sample: its VCF and coverage are what the summary and the gene steps use. */
// ---- One sample among several: which of its variants are only in it and which are shared ----
const infoOf = (field) => Object.fromEntries(field.split(';').map((kv) => { const i = kv.indexOf('='); return i > 0 ? [kv.slice(0, i), kv.slice(i + 1)] : [kv, '']; }));
const variantKey = (f) => `${f[0]}\t${f[1]}\t${f[3]}\t${f[4].split(',')[0]}`;
const escapeSampleName = (name) => name.replace(/[,;=% \x00-\x1f\x7f]/g, (c) => `%${c.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);

/** The finished samples, and for each variant (contig, position, REF, ALT) the samples that call it with their AF. */
function sharingIndex() {
  const finished = samples.filter((s) => s.name && isFinished(s));
  const index = new Map();
  for (const s of finished) {
    for (const line of s.vcf.split('\n')) {
      if (!line || line[0] === '#') continue;
      const f = line.split('\t');
      const af = parseFloat(infoOf(f[7]).AF);
      const key = variantKey(f);
      if (!index.has(key)) index.set(key, []);
      index.get(key).push({ name: s.name, af });
    }
  }
  return { finished, index };
}

/**
 * The VCF of one sample with the sharing tags `NS` (how many finished samples call the variant), `SAMPLES` and `AFS` (the AF in each),
 * like the merged VCF of all samples, so that the table, the summary and the annotated VCF say what is shared. With one finished
 * sample there is nothing to share and the text is returned as it is.
 */
function withSharing(vcfText) {
  const { finished, index } = sharingIndex();
  if (finished.length < 2) return vcfText;
  const header = [
    `##qmutfind_samples=${finished.length}`,
    '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of finished samples with a call">',
    '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Samples with a call (a name has , ; = space and % written as %XX)">',
    '##INFO=<ID=AFS,Number=.,Type=Float,Description="Alternate allele fraction in each sample listed in SAMPLES, in the same order">',
  ];
  return vcfText.split('\n').flatMap((line) => {
    if (line.startsWith('#CHROM')) return [...header, line];
    if (!line || line[0] === '#') return [line];
    const f = line.split('\t');
    const calls = index.get(variantKey(f)) || [];
    const tags = `NS=${calls.length};SAMPLES=${calls.map((c) => escapeSampleName(c.name)).join(',')};AFS=${calls.map((c) => (Number.isFinite(c.af) ? c.af.toFixed(4) : '.')).join(',')}`;
    f[7] = f[7] === '.' || f[7] === '' ? tags : `${f[7]};${tags}`;
    return [f.join('\t')];
  }).join('\n');
}

/** For a VCF with the sharing tags: how many variants are in one sample only, in two or more, and in every sample. */
function sharingCounts(vcfText) {
  const total = Number((/^##qmutfind_samples=(\d+)/m.exec(vcfText || '') || [])[1]) || 0;
  if (!total) return null;
  const c = { total, variants: 0, one: 0, shared: 0, every: 0 };
  for (const line of vcfText.split('\n')) {
    if (!line || line[0] === '#') continue;
    const ns = Number(infoOf(line.split('\t')[7]).NS);
    c.variants += 1;
    if (ns === 1) c.one += 1;
    if (ns >= 2) c.shared += 1;
    if (ns === total) c.every += 1;
  }
  return c;
}

const countRecords = (vcfText) => (vcfText ? vcfText.split('\n').filter((l) => l && l[0] !== '#').length : 0);

/** Subtracting the parent is possible after a combine with a parent, for a sample that is not the parent and for all samples together. */
function subtractAvailable() {
  if (!lastCombined || !lastCombined.parent) return false;
  if (viewingAll()) return true;
  const s = selectedSample();
  return Boolean(s && s.name !== lastCombined.parent);
}
function subtracting() { return $('subtract-parent').checked && subtractAvailable() && !$('vcf-file').files[0]; }

function updateSubtractControl() {
  const box = $('subtract-parent');
  const ok = subtractAvailable();
  box.disabled = !ok;
  if (!ok) box.checked = false;
  show('subtract-note', ok
    ? `Hides the variants that the parent ${lastCombined.parent} has. A variant is also left out when the parent has weak evidence for it or too little coverage (it is "unknown" in the combined table), because then it cannot be called new.`
    : 'Needs a combine with a parent sample (section "Combine the samples"). The parent itself has nothing to subtract.');
}

/** The records of `vcfText` that satisfy `keep(fields)`, with the header lines. */
function filterVcf(vcfText, keep) {
  return vcfText.split('\n').filter((l) => !l || l[0] === '#' || keep(l.split('\t'))).join('\n');
}

/** The records of the merged VCF where the variant is new in the sample `name` (or in any sample when `name` is null). */
function newInParent(vcfText, name) {
  const keyOf = (f) => `${f[0]}\t${f[1]}\t${f[3]}\t${f[4]}`;
  const wanted = new Set();
  for (const line of lastCombined.vcf.split('\n')) {
    if (!line || line[0] === '#') continue;
    const f = line.split('\t');
    const info = Object.fromEntries(f[7].split(';').map((kv) => { const i = kv.indexOf('='); return i > 0 ? [kv.slice(0, i), kv.slice(i + 1)] : [kv, '']; }));
    const news = info.NEW ? info.NEW.split(',').map(decodeInfo) : [];
    if (name === null ? news.length > 0 : news.includes(name)) wanted.add(keyOf(f));
  }
  return filterVcf(vcfText, (f) => wanted.has(keyOf(f)));
}

function useSelectedSample() {
  const all = viewingAll();
  const s = all ? null : selectedSample();
  updateSubtractControl();
  const subtract = subtracting();
  let text = all ? lastCombined.vcf : (s ? s.vcf : null);
  if (text && subtract) text = newInParent(text, all ? null : s.name);
  if (text && s && !all) text = withSharing(text); // what is only in this sample and what is shared with the others
  state.vcf = text;
  state.depth = s ? s.bed : null; // the coverage of separate samples is not merged
  state.annotated = null;
  for (const id of ['download-vcf', 'download-depth-bed']) $(id).hidden = true;
  const suffix = subtract ? '.new_vs_parent' : '';
  if (all) {
    offer('download-vcf', `combined${suffix}.vcf`, text, 'text/plain');
    show('summary', `${countRecords(text).toLocaleString()} distinct variants in ${lastCombined.names.length} samples. Go on to the summary or the gene steps.`);
  } else if (s) {
    offer('download-vcf', `${s.name}${suffix}.vcf`, text, 'text/plain');
    if (s.bed) offer('download-depth-bed', `${s.name}.depth.bed`, s.bed, 'text/plain');
    show('summary', `Calls: ${countRecords(text).toLocaleString()}. Go on to the summary or the gene steps.`);
  } else {
    show('summary', '');
  }
  describeSource();
}
$('result-sample').addEventListener('change', useSelectedSample);
$('subtract-parent').addEventListener('change', useSelectedSample);

function renderSamples() {
  const body = $('strain-table').querySelector('tbody');
  body.replaceChildren();
  for (const s of samples) body.appendChild(sampleRow(s));
  refreshChoices();
}

function sampleRow(s) {
  const row = document.createElement('tr');
  row.dataset.name = s.name;
  row.dataset.id = String(s.id);
  const imported = s.state.startsWith('imported');
  const locked = imported || s.state === 'running' || s.state === 'done' || s.state.startsWith('done');
  const labels = ['Sample name', 'Source', 'Reads', 'State', 'Reads used', 'Calls', 'Seconds', 'Note'];
  const td = (cls) => { const c = document.createElement('td'); if (cls) c.className = cls; c.dataset.label = labels[row.children.length]; row.appendChild(c); return c; };
  const edited = () => { // changing a failed or cancelled sample makes it ready to run again
    if (s.state.startsWith('failed') || s.state === 'cancelled') { s.state = 'waiting'; row.querySelector('.strain-state').textContent = 'waiting'; }
    if (lastCombined) invalidateCombined(); // the combined table is out of date
  };

  const nameCell = td();
  const name = Object.assign(document.createElement('input'), { type: 'text', className: 'sample-name', value: s.name, disabled: locked });
  name.setAttribute('aria-label', 'Sample name');
  name.addEventListener('input', () => { s.name = name.value.trim(); row.dataset.name = s.name; edited(); });
  name.addEventListener('change', () => { refreshChoices(); useSelectedSample(); });
  const nameLine = Object.assign(document.createElement('div'), { className: 'name-line' });
  nameCell.appendChild(nameLine);
  nameLine.appendChild(name);

  const sourceCell = td();
  if (imported) {
    sourceCell.textContent = 'imported';
  } else {
    const source = document.createElement('select');
    source.className = 'sample-source';
    source.setAttribute('aria-label', 'Where the reads come from');
    source.disabled = locked;
    for (const [value, label] of [['sra', 'SRA run'], ['url', 'Web addresses'], ['local', 'Local files']]) source.appendChild(new Option(label, value));
    source.value = s.source;
    sourceCell.appendChild(source);
    source.addEventListener('change', () => { s.source = source.value; edited(); showGroup(); });
  }

  const readsCell = td('sample-reads');
  let showGroup = () => {};
  if (imported) {
    readsCell.textContent = s.files.join(', ');
  } else {
    const sra = document.createElement('div');
    sra.className = 'group-sra';
    const accession = Object.assign(document.createElement('input'), { type: 'text', className: 'sample-sra', placeholder: 'SRR6277229', value: s.accession, disabled: locked, autocomplete: 'off' });
    accession.setAttribute('aria-label', 'SRA run accession');
    const check = Object.assign(document.createElement('button'), { type: 'button', className: 'sample-check', textContent: 'Check', disabled: locked });
    const info = Object.assign(document.createElement('div'), { className: 'sample-info readout' });
    info.textContent = s.info || '';
    accession.addEventListener('input', () => { s.accession = accession.value; s.info = ''; info.textContent = ''; edited(); });
    check.addEventListener('click', async () => {
      info.textContent = 'Looking the run up in ENA';
      try {
        s.info = (await lookupRun(accession.value, readLimit())).info;
      } catch (error) {
        s.info = error.message;
      }
      info.textContent = s.info;
    });
    sra.append(accession, check, info);

    const url = document.createElement('div');
    url.className = 'group-url';
    const urls = Object.assign(document.createElement('textarea'), { className: 'sample-urls', rows: 2, placeholder: 'https://... one address per line (R1, R2, more)', value: s.urls, disabled: locked });
    urls.setAttribute('aria-label', 'Reads file addresses, one per line');
    urls.addEventListener('input', () => { s.urls = urls.value; edited(); });
    url.appendChild(urls);

    const local = document.createElement('div');
    local.className = 'group-local';
    const files = Object.assign(document.createElement('input'), { type: 'file', className: 'sample-files', multiple: true, disabled: locked });
    files.setAttribute('aria-label', 'Reads files from this computer (select several if the sample has more than one)');
    const names = Object.assign(document.createElement('div'), { className: 'sample-files-names hint' });
    const showNames = () => { names.textContent = s.files.length ? `${s.fromExample ? 'From the example: ' : ''}${s.files.map((f) => f.name).join(', ')}` : ''; };
    files.addEventListener('change', () => { s.files = [...files.files]; s.fromExample = false; edited(); showNames(); });
    showNames();
    local.append(files, names);

    showGroup = () => { sra.hidden = s.source !== 'sra'; url.hidden = s.source !== 'url'; local.hidden = s.source !== 'local'; };
    readsCell.append(sra, url, local);
    showGroup();
  }

  const cell = (text, cls) => { const c = td(cls); c.textContent = text; };
  cell(s.state, 'strain-state');
  cell(s.reads === undefined ? '' : String(s.reads), 'strain-reads');
  cell(s.calls === undefined ? '' : String(s.calls), 'strain-calls');
  cell(s.seconds === undefined ? '' : s.seconds.toFixed(1), 'strain-seconds');

  const noteCell = td();
  const note = Object.assign(document.createElement('input'), { type: 'text', className: 'sample-note', value: s.note || '', disabled: imported });
  note.setAttribute('aria-label', 'Note');
  note.addEventListener('input', () => { s.note = note.value.trim(); });
  noteCell.appendChild(note);

  const remove = Object.assign(document.createElement('button'), { type: 'button', className: 'strain-remove', textContent: '−', disabled: s.state === 'running' });
  remove.setAttribute('aria-label', 'Remove this sample');
  remove.addEventListener('click', () => {
    const at = samples.indexOf(s);
    if (at >= 0 && s.state !== 'running') { samples.splice(at, 1); invalidateCombined(); renderSamples(); saveSoon(); }
  });
  nameLine.prepend(remove); // the minus button sits in front of the name
  return row;
}

$('strain-add').addEventListener('click', () => {
  newSample();
  renderSamples();
  const names = $('strain-table').querySelectorAll('.sample-name');
  names[names.length - 1].focus();
});

/** The example's reads become a sample; a new genome from the menu or an address drops them, because they belong to the old one. */
function dropExampleSamples() {
  for (let i = samples.length - 1; i >= 0; i--) if (samples[i].fromExample && !isFinished(samples[i])) samples.splice(i, 1);
  if (!samples.length) newSample();
  renderSamples();
}

/** Get the reads of one sample as files for the worker. Returns { files, expected, kind, anyCut, note }. */
async function prepareReads(s, limit, signal) {
  const status = (text) => sampleStatus(text);
  if (s.source === 'local') return { files: s.files, expected: s.expected || null, kind: null, note: '' };
  if (s.source === 'url') {
    const addresses = s.urls.split('\n').map((x) => x.trim()).filter(Boolean);
    const urls = addresses.map((a) => checkUrl(directAddress(a), `the reads of ${s.name}`)); // all checked before any download
    const files = [];
    const budget = { left: MAX_DOWNLOAD }; // the files of one sample together stay within the download limit
    let anyCut = false;
    for (const [i, url] of urls.entries()) {
      const item = { name: decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'download'), url: url.href, bytes: 0 };
      const label = `${s.name} file ${i + 1} of ${urls.length}`;
      if (limit) {
        const result = await downloadFirstReads(item, location.href, label, 'strain-status', signal, limit, budget);
        files.push(result.file);
        anyCut ||= result.cut;
      } else {
        files.push(await download(item, location.href, label, 'strain-status', signal, budget));
      }
    }
    return { files, expected: null, kind: null, anyCut };
  }
  status(`Looking up ${s.accession.trim()} in ENA`);
  const run = await lookupRun(s.accession, limit);
  s.info = run.info;
  const files = [];
  const budget = { left: MAX_DOWNLOAD };
  let anyCut = false;
  for (const [i, item] of run.files.entries()) {
    const label = `${run.accession} file ${i + 1} of ${run.files.length}`;
    if (limit) {
      const result = await downloadFirstReads(item, location.href, label, 'strain-status', signal, limit, budget);
      files.push(result.file);
      anyCut ||= result.cut;
    } else {
      files.push(await download(item, location.href, label, 'strain-status', signal, budget));
    }
  }
  const perFile = run.paired ? run.perFile : Math.round(run.perFile / run.files.length);
  return { files, expected: limit ? run.files.length * Math.min(limit, perFile) : run.count, kind: run.kind, anyCut };
}

$('strain-cancel').addEventListener('click', () => {
  sampleCancelled = true;
  if (activeAbort) activeAbort.abort();
  if (activeWorker) { activeWorker.w.terminate(); activeWorker.reject(new DOMException('cancelled', 'AbortError')); activeWorker = null; }
});

/** A reason the waiting samples cannot run, or ''. Blank rows are ignored. */
function problemWithSamples() {
  const todo = samples.filter((s) => !isBlank(s));
  const seen = new Set();
  for (const s of todo) {
    if (!s.name) return 'Give every sample a name.';
    if (/[\u0000-\u001f]/.test(s.name)) return `The sample name "${s.name}" must not hold a tab or a line break.`;
    if (seen.has(s.name)) return `The sample name "${s.name}" is used more than once.`;
    seen.add(s.name);
  }
  for (const s of todo.filter((x) => x.state === 'waiting' || x.state === 'cancelled')) {
    const has = s.source === 'sra' ? s.accession.trim() : s.source === 'url' ? s.urls.trim() : s.files.length;
    if (!has) return `The sample "${s.name}" has no reads: ${s.source === 'sra' ? 'give an SRA run accession' : s.source === 'url' ? 'give one or more https addresses' : 'choose one or more files'}.`;
  }
  return '';
}

$('strain-run').addEventListener('click', async () => {
  if (sampleRunning || busy || loading) { sampleStatus('Something else is running. Wait until it has finished.'); return; }
  const reference = refFile();
  if (!reference) { sampleStatus('Choose the reference genome first (section 1).'); return; }
  const problem = problemWithSamples();
  if (problem) { sampleStatus(problem); return; }
  if (!samples.some((s) => !isBlank(s) && (s.state === 'waiting' || s.state === 'cancelled'))) { sampleStatus('No sample is waiting.'); return; }
  let limit;
  try { limit = readLimit(); } catch (error) { sampleStatus(error.message); return; }
  sampleRunning = true;
  sampleCancelled = false;
  setBusy(true);
  $('strain-run').disabled = true;
  $('strain-add').disabled = true;
  $('strain-cancel').hidden = false;
  try {
    for (;;) {
      const next = samples.find((s) => !isBlank(s) && (s.state === 'waiting' || s.state === 'cancelled'));
      if (!next || sampleCancelled) break;
      next.state = 'running';
      renderSamples();
      sampleStatus(`Running ${next.name}`);
      activeAbort = new AbortController();
      try {
        const prepared = await prepareReads(next, limit, activeAbort.signal);
        const m = await callWorker({
          step: 'run', reference, reads: prepared.files, sample: next.name, kind: prepared.kind || $('kind').value,
          depthWindow: DEPTH_WINDOW, evidence: true, maxReads: limit,
        }, (p) => sampleStatus(`Running ${next.name}: ${p.reads.toLocaleString()} reads mapped`));
        const anyCut = prepared.anyCut || (m.files || []).some((f) => f.cut);
        const wrong = prepared.expected && !anyCut && m.reads !== prepared.expected;
        Object.assign(next, {
          state: wrong ? `done: WARNING the page counted ${m.reads.toLocaleString()} reads but ${prepared.expected.toLocaleString()} were expected, so a download may be incomplete` : 'done',
          reads: m.reads, seconds: m.seconds, vcf: m.vcf, bed: m.depth, evidence: m.evidence, version: m.version,
          calls: m.vcf.split('\n').filter((l) => l && l[0] !== '#').length,
        });
        if (!$('result-sample').value) refreshChoices();
        useSelectedSample(); // a sample that finishes changes what is shared
      } catch (error) {
        next.state = error.name === 'AbortError' ? 'cancelled' : `failed: ${error.message}`;
      }
      activeAbort = null;
      invalidateCombined();
      renderSamples();
      saveSoon();
    }
    if (sampleCancelled) for (const s of samples) if (s.state === 'waiting') s.state = 'cancelled';
    sampleStatus(sampleCancelled ? 'Cancelled. Samples that were not finished can run again.' : 'Done. Combine the finished samples, look at one in the results, or export the catalog.');
  } finally {
    sampleRunning = false;
    setBusy(false);
    $('strain-run').disabled = false;
    $('strain-add').disabled = false;
    $('strain-cancel').hidden = true;
    renderSamples();
  }
});

// The combined table
/** The combined table is out of date (a sample changed or was removed): drop it and say so. */
function invalidateCombined() {
  lastCombined = null;
  refreshChoices();
  renderCombined();
  useSelectedSample(); // the subtract choice and the all-samples view depend on the combine
}

// Sorting the combined table: a header button sorts up, down, and back to the order of the file, like the variant table.
// A number column sorts by value; a sample column by the AF of its call (a cell that is not a call counts as empty); empty cells
// (and the dot of "no status") stay last in both directions.
let combinedSort = { column: null, direction: 1 };

function combinedKind(header, index, nSamples) {
  if (header === 'POS' || header === 'strains') return 'number';
  if (index >= 6 && index < 6 + nSamples) return 'af';
  return 'text';
}

function combinedValue(kind, text) {
  if (kind === 'number') { const n = parseFloat(text); return Number.isFinite(n) ? n : null; }
  if (kind === 'af') { const m = /^([0-9.]+):\d+$/.exec(text); return m ? parseFloat(m[1]) : null; }
  return text === '' || text === '.' ? null : text;
}

function sortCombinedRows(rows, header, nSamples) {
  if (combinedSort.column === null) return rows;
  const kind = combinedKind(header[combinedSort.column], combinedSort.column, nSamples);
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
  const c = combinedSort.column;
  return [...rows].sort((a, b) => {
    const va = combinedValue(kind, a[c]);
    const vb = combinedValue(kind, b[c]);
    if (va === null || vb === null) return va === vb ? 0 : (va === null ? 1 : -1);
    return combinedSort.direction * (kind === 'text' ? collator.compare(va, vb) : va - vb);
  });
}

function sortCombinedBy(index) {
  if (combinedSort.column === index) {
    if (combinedSort.direction === 1) combinedSort.direction = -1;
    else combinedSort.column = null;
  } else {
    combinedSort = { column: index, direction: 1 };
  }
  renderCombined();
}

function renderCombined() {
  const head = $('combined-table').querySelector('thead');
  const body = $('combined-table').querySelector('tbody');
  head.replaceChildren();
  body.replaceChildren();
  $('combined-empty').hidden = Boolean(lastCombined);
  if (!lastCombined) { $('combined-title').textContent = ''; return; }
  const total = lastCombined.table.split('\n').filter(Boolean).length - 1;
  $('combined-title').textContent = `Combined table: ${total.toLocaleString()} variants in ${lastCombined.names.length} samples${lastCombined.parent ? ` (parent: ${lastCombined.parent})` : ''}`;
  const lines = lastCombined.table.split('\n').filter(Boolean).map((l) => l.split('\t'));
  const header = lines[0];
  const nSamples = lastCombined.names.length;
  const firstCell = 6;
  const hasStatus = header.length > firstCell + nSamples;
  const filter = $('combined-filter').value;
  const keep = (f) => {
    const cells = f.slice(firstCell, firstCell + nSamples);
    const called = cells.filter((c) => /^[0-9.]+:\d+$/.test(c)).length;
    if (filter === 'shared') return called === nSamples;
    if (filter === 'one') return called === 1;
    if (filter === 'new') return hasStatus && f.slice(firstCell + nSamples).includes('new');
    return true;
  };
  const tr = document.createElement('tr');
  for (const [index, h] of header.entries()) {
    const th = document.createElement('th');
    th.setAttribute('aria-sort', combinedSort.column === index ? (combinedSort.direction === 1 ? 'ascending' : 'descending') : 'none');
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'sort';
    button.textContent = h;
    button.addEventListener('click', () => sortCombinedBy(index));
    th.appendChild(button);
    tr.appendChild(th);
  }
  head.appendChild(tr);
  const rows = sortCombinedRows(lines.slice(1).filter(keep), header, nSamples);
  for (const f of rows.slice(0, MAX_ROWS)) {
    const row = document.createElement('tr');
    for (const text of f) { const c = document.createElement('td'); c.textContent = text; if (text === 'new') c.className = 'new'; row.appendChild(c); }
    body.appendChild(row);
  }
  if (rows.length > MAX_ROWS) {
    const row = document.createElement('tr');
    const c = document.createElement('td');
    c.colSpan = header.length;
    c.textContent = `Showing the first ${MAX_ROWS} of ${rows.length.toLocaleString()} rows. The download has all of them.`;
    row.appendChild(c);
    body.appendChild(row);
  }
}
$('combined-filter').addEventListener('change', renderCombined);

$('strain-combine').addEventListener('click', async () => {
  if (sampleRunning || busy) { show('combine-status', 'Wait until the samples have finished running, then combine.'); return; }
  const done = samples.filter((s) => s.name && isFinished(s));
  if (done.length < 2) { show('combine-status', 'Combine needs at least two finished samples.'); return; }
  const minDepth = Number($('combine-min-depth').value);
  if (!Number.isFinite(minDepth) || minDepth < 0) { show('combine-status', 'The minimum depth must be a number, 0 or more.'); return; }
  const parent = $('strain-parent').value;
  if (parent && !done.some((s) => s.name === parent)) { show('combine-status', `The parent ${parent} has no finished result.`); return; }
  show('combine-status', 'Combining');
  try {
    const m = await callWorker({
      step: 'combine', names: done.map((s) => s.name), vcfs: done.map((s) => s.vcf), beds: done.map((s) => s.bed || ''),
      evidence: done.map((s) => s.evidence || ''), parent: parent || undefined, minDepth,
    });
    lastCombined = { table: m.table, summary: m.summary, vcf: m.vcf, parent, minDepth, names: done.map((s) => s.name) };
    combinedSort = { column: null, direction: 1 }; // a new table starts in the order of the file
    refreshChoices();
    offer('download-combined', 'combined.tsv', m.table, 'text/tab-separated-values');
    renderCombined();
    $('combined-title').scrollIntoView({ block: 'nearest' }); // the table is near the button that made it
    useSelectedSample(); // a parent may now be known
    const n = m.table.split('\n').filter(Boolean).length - 1;
    show('combine-status', `Done: ${n.toLocaleString()} variant(s) in ${done.length} samples.${m.warnings ? ` Warnings: ${m.warnings.split('\n').join(' ')}` : ''}`);
    saveSoon();
  } catch (error) {
    show('combine-status', `Error: ${error.message}`);
  }
});

// Export and import
function catalogFileName(index, name, extension) {
  const safe = name.replace(/[^A-Za-z0-9._-]+/g, '_').slice(0, 60) || 'sample';
  return `strains/${String(index + 1).padStart(2, '0')}_${safe}.${extension}`;
}

/** The files of the catalog zip for the sample list as it is now: { files, count }, or null when no sample has a name. */
function buildCatalog() {
  const listed = samples.filter((s) => s.name);
  if (!listed.length) return null;
  const files = {};
  const entries = listed.map((s, i) => {
    const entry = {
      name: s.name, note: s.note || '',
      files: s.source === 'sra' ? [s.accession.trim()] : s.source === 'url' ? s.urls.split('\n').map((x) => x.trim()).filter(Boolean) : s.files.map((f) => (typeof f === 'string' ? f : f.name)),
      source: s.state.startsWith('imported') ? 'imported' : s.source,
      status: s.state.startsWith('imported') ? s.state.slice('imported: '.length) : s.state.split(':')[0],
    };
    if (s.state.startsWith('failed')) entry.error = s.state.slice('failed: '.length);
    if (isFinished(s)) {
      Object.assign(entry, { reads: s.reads, calls: s.calls, seconds: s.seconds, qmutfind_version: s.version || '' });
      entry.vcf = catalogFileName(i, s.name, 'vcf');
      files[entry.vcf] = strToU8(s.vcf);
      if (s.bed) { entry.bed = catalogFileName(i, s.name, 'depth.bed'); files[entry.bed] = strToU8(s.bed); }
      if (s.evidence) { entry.evidence = catalogFileName(i, s.name, 'evidence.vcf'); files[entry.evidence] = strToU8(s.evidence); }
    }
    return entry;
  });
  const reference = refFile();
  let readLimitReads = null;
  try { readLimitReads = readLimit(); } catch (error) { /* the box holds something else; the catalog then has no limit */ }
  const catalog = {
    format: CATALOG_FORMAT, version: CATALOG_VERSION, created: new Date().toISOString(),
    qmutfind_version: (listed.find((s) => s.version) || {}).version || '',
    reference: reference ? { name: reference.name, size: reference.size } : null,
    settings: { kind: $('kind').value, depth_window: DEPTH_WINDOW, read_limit_per_file: readLimitReads, min_depth: lastCombined ? lastCombined.minDepth : Number($('combine-min-depth').value) },
    parent: $('strain-parent').value || '',
    strains: entries,
  };
  files['catalog.json'] = strToU8(JSON.stringify(catalog, null, 1));
  if (lastCombined) { files['combined.tsv'] = strToU8(lastCombined.table); files['summary.tsv'] = strToU8(lastCombined.summary); }
  return { files, count: listed.length };
}

$('strain-export').addEventListener('click', () => {
  const built = buildCatalog();
  if (!built) { show('catalog-status', 'The list is empty. Add samples first.'); return; }
  offer('download-catalog', 'qmutfind-catalog.zip', zipSync(built.files), 'application/zip');
  show('catalog-status', `Exported ${built.count} sample(s). Use the download link.`);
});

/** Check an imported catalog and return the samples to load. Throws an Error with a reason. */
function readCatalog(archive) {
  const raw = archive['catalog.json'];
  if (!raw) throw new Error('the zip has no catalog.json, so it is not a qmutfind catalog.');
  let catalog;
  try { catalog = JSON.parse(strFromU8(raw)); } catch (error) { throw new Error('catalog.json is not valid JSON (catalog unreadable).'); }
  if (!catalog || catalog.format !== CATALOG_FORMAT) throw new Error(`catalog.json has the format "${catalog && catalog.format}", not ${CATALOG_FORMAT}.`);
  if (catalog.version !== CATALOG_VERSION) throw new Error(`catalog.json has version ${catalog.version}; this page reads version ${CATALOG_VERSION}.`);
  if (!Array.isArray(catalog.strains) || catalog.strains.length > 1000) throw new Error('catalog.json lists no samples, or too many (limit 1000).');
  const text = (path) => {
    if (typeof path !== 'string' || !Object.prototype.hasOwnProperty.call(archive, path)) throw new Error(`catalog.json refers to the file "${path}", which is not in the zip.`);
    return strFromU8(archive[path]);
  };
  const names = new Set();
  const loaded = catalog.strains.map((e) => {
    if (!e || typeof e.name !== 'string' || !e.name.trim() || /[\u0000-\u001f]/.test(e.name)) throw new Error('catalog.json has a sample without a valid name.');
    if (names.has(e.name)) throw new Error(`catalog.json lists the sample name "${e.name}" more than once.`);
    names.add(e.name);
    const s = { name: e.name, note: typeof e.note === 'string' ? e.note : '', files: Array.isArray(e.files) ? e.files.map(String) : [], state: `imported: ${String(e.status || 'unknown')}` };
    if (e.vcf !== undefined) {
      s.vcf = text(e.vcf);
      if (e.bed !== undefined) s.bed = text(e.bed);
      if (e.evidence !== undefined) s.evidence = text(e.evidence);
      for (const k of ['reads', 'calls', 'seconds']) if (typeof e[k] === 'number') s[k] = e[k];
      s.version = typeof e.qmutfind_version === 'string' ? e.qmutfind_version : '';
    }
    return s;
  });
  return { catalog, loaded };
}

/** Unpack a catalog zip (bytes) with the size checks of an import. Throws an Error with a reason. */
function unpackCatalog(bytes) {
  let total = 0;
  try {
    return unzipSync(bytes, {
      filter: (f) => { total += f.originalSize; if (total > MAX_CATALOG_BYTES) throw new Error('too large'); return true; },
    });
  } catch (error) {
    throw new Error(error.message === 'too large' ? `the zip unpacks to more than ${mb(MAX_CATALOG_BYTES)}.` : 'the file is not a zip file that this page can read.');
  }
}

/** Replace the sample list and the settings with a catalog. Returns { catalog, loaded }. */
function applyCatalog(archive) {
  const { catalog, loaded } = readCatalog(archive);
  samples.splice(0, samples.length);
  for (const s of loaded) newSample(s);
  lastCombined = null;
  renderSamples();
  renderCombined();
  const settings = catalog.settings || {};
  if (['short', 'long'].includes(settings.kind)) $('kind').value = settings.kind;
  if (Number.isFinite(settings.min_depth)) $('combine-min-depth').value = String(settings.min_depth);
  if (Number.isFinite(settings.read_limit_per_file)) { $('read-limit').value = String(settings.read_limit_per_file / 1000); describeLimit(); }
  if (typeof catalog.parent === 'string' && sampleByName(catalog.parent)) $('strain-parent').value = catalog.parent;
  return { catalog, loaded };
}

$('catalog-file').addEventListener('change', async () => {
  const file = $('catalog-file').files[0];
  if (!file) return;
  if (sampleRunning || busy || loading) { show('catalog-status', 'Wait until the samples have finished running, then import.'); $('catalog-file').value = ''; return; }
  try {
    if (file.size > MAX_CATALOG_BYTES) throw new Error(`the file is too large (${mb(file.size)}; the limit is ${mb(MAX_CATALOG_BYTES)}).`);
    const archive = unpackCatalog(new Uint8Array(await file.arrayBuffer()));
    const { catalog, loaded } = applyCatalog(archive);
    show('catalog-status', `Loaded ${loaded.length} sample(s) from ${file.name}${catalog.reference ? ` (made with the reference ${catalog.reference.name})` : ''}. Imported samples cannot run again.`);
    saveSoon();
  } catch (error) {
    show('catalog-status', `Error: ${error.message}`);
  } finally {
    $('catalog-file').value = '';
  }
});

// ---- Saved results (IndexedDB) ----------------------------------------------------------------------------------------
// The finished samples, the combine and the reference are saved in the browser after each sample and each combine, and restored
// on request. The reads are never saved. Results that were saved earlier and not restored are not overwritten by new ones.
const KEEP_KEY = 'qmutfind-keep';
// A reference or gene file above this is not saved. `?save_limit_mb=N` can lower the limit (the tests use that).
const MAX_SAVED_FILE = Math.min(parseFloat(new URLSearchParams(location.search).get('save_limit_mb')) || 300, 300) * 1e6;
let savedHeader = null; // { saved, count, bytes, combined, reference, gff } of what is in the store, or null
let saveHeld = false; // saved results exist and were not restored: new results are not saved over them
let saveTimer = null;
let heldFiles = { reference: null, gff: null }; // the File objects whose bytes the store holds, so each is written once
const plural = (n, word) => `${n.toLocaleString()} ${word}${n === 1 ? '' : 's'}`;

function keepResults() {
  return $('keep-results').checked && storeAvailable();
}

async function storageUse() {
  const e = await estimateStorage();
  return e && e.quota ? ` Storage used by this site in this browser: ${mb(e.usage)} of ${mb(e.quota)}.` : '';
}

/** The line under "Keep my results in this browser" and the buttons next to it. */
async function showSaved(note = '') {
  if (!storeAvailable()) {
    show('saved-status', 'Saving in the browser is not available here (the browser blocks or lacks IndexedDB). Use Export to keep the results.');
    for (const id of ['saved-restore', 'saved-delete', 'saved-save-now']) $(id).hidden = true;
    $('keep-results').disabled = true;
    return;
  }
  const off = $('keep-results').checked ? '' : ' Keeping is off: new results are not saved.';
  let text;
  if (!savedHeader) {
    text = `Nothing is saved in this browser.${off}`;
  } else {
    const when = new Date(savedHeader.saved).toLocaleString();
    const what = `${plural(savedHeader.count, 'sample')}${savedHeader.combined ? ', the combined table' : ''}${savedHeader.reference ? `, the reference ${savedHeader.reference}` : ''} (${mb(savedHeader.bytes)}, saved ${when})`;
    text = saveHeld
      ? `Found saved results: ${what}. Restore them to continue, or delete them. New results are not saved until you restore or delete these (or press "Save now").`
      : `Saved in this browser: ${what}.${off}${await storageUse()}`;
  }
  show('saved-status', note ? `${text} ${note}` : text);
  $('saved-restore').hidden = !savedHeader;
  $('saved-delete').hidden = !savedHeader;
  $('saved-save-now').hidden = !(saveHeld || !$('keep-results').checked) || !samples.some(isFinished);
}

function saveSoon() {
  if (!keepResults() || saveHeld) { if (saveHeld && samples.some(isFinished)) showSaved(); return; }
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => saveNow(false), 800);
}

/** Write the finished samples, the combine settings and the reference to the store. */
async function saveNow(force) {
  if (!storeAvailable() || (!force && (!keepResults() || saveHeld))) return;
  const finished = samples.filter(isFinished);
  if (!finished.length) { showSaved('There are no finished samples to save.'); return; } // an empty list never replaces saved results
  try {
    const built = buildCatalog();
    const unpacked = Object.values(built.files).reduce((n, f) => n + f.length, 0);
    if (unpacked > MAX_CATALOG_BYTES) throw new Error(`the results are too large to save (${mb(unpacked)}; the limit is ${mb(MAX_CATALOG_BYTES)}). Use Export`);
    const zip = zipSync(built.files);
    const entries = { results: zip };
    const held = { reference: null, gff: null };
    const names = { reference: '', gff: '' };
    let bytes = zip.length;
    for (const [slot, file] of [['reference', refFile()], ['gff', gffFile()]]) {
      if (!file || file.size > MAX_SAVED_FILE) { entries[slot] = null; continue; } // not saved, and an older copy is removed
      names[slot] = file.name;
      held[slot] = file;
      bytes += file.size;
      if (heldFiles[slot] !== file) entries[slot] = { name: file.name, bytes: new Uint8Array(await file.arrayBuffer()) }; // written only when it changed
    }
    const header = { saved: new Date().toISOString(), count: finished.length, bytes, combined: Boolean(lastCombined), reference: names.reference, gff: names.gff };
    entries.header = header;
    // Without `force` the write is made only if no other tab saved since this page last looked at the store.
    await writeStore(entries, force ? {} : { expect: savedHeader ? savedHeader.saved : null });
    heldFiles = held;
    savedHeader = header;
    saveHeld = false;
    await showSaved();
  } catch (error) {
    if (error && error.code === 'CHANGED') {
      try { const found = await readStore(['header']); savedHeader = found.header || null; } catch (e) { /* the old header stays */ }
      saveHeld = Boolean(savedHeader);
      await showSaved();
      return;
    }
    show('saved-status', `Could not save in this browser: ${error.message || error}. The results are still in the page; use Export to keep them.`);
  }
}

async function restoreSaved() {
  if (sampleRunning || busy || loading) { show('saved-status', 'Wait until the samples have finished running, then restore.'); return; }
  try {
    show('saved-status', 'Restoring');
    const found = await readStore(['results', 'reference', 'gff', 'header']);
    if (!found.results) throw new Error('the saved results are missing');
    const archive = unpackCatalog(found.results);
    const header = found.header || savedHeader || {};
    const { loaded } = applyCatalog(archive);
    const when = header.saved ? new Date(header.saved).toLocaleString() : '';
    for (const slot of ['reference', 'gff']) {
      const entry = found[slot];
      if (!entry) continue;
      $(slot).value = ''; // a file chosen before would win over the restored one (refFile prefers the input)
      example[slot] = new File([entry.bytes], entry.name);
      heldFiles[slot] = example[slot];
      example.source[slot] = `this browser (saved ${when})`;
    }
    if (found.reference) { example.source.reads = ''; }
    describeSources();
    saveHeld = false;
    show('catalog-status', `Restored ${plural(loaded.length, 'sample')} from this browser. Restored samples cannot run again.`);
    if (header.combined && loaded.filter(isFinished).length >= 2) $('strain-combine').click();
    await showSaved();
  } catch (error) {
    show('saved-status', `Error: could not restore: ${error.message || error}`);
  }
}

$('saved-restore').addEventListener('click', restoreSaved);
$('saved-save-now').addEventListener('click', () => saveNow(true));
$('saved-delete').addEventListener('click', async () => {
  try {
    await clearStore();
    savedHeader = null;
    saveHeld = false;
    heldFiles = { reference: null, gff: null };
    await showSaved();
  } catch (error) {
    show('saved-status', `Error: could not delete: ${error.message || error}`);
  }
});
$('keep-results').addEventListener('change', () => {
  try { localStorage.setItem(KEEP_KEY, $('keep-results').checked ? '1' : '0'); } catch (error) { /* the choice is then not remembered */ }
  if ($('keep-results').checked) saveSoon();
  showSaved();
});

(async function startSaved() {
  try { $('keep-results').checked = localStorage.getItem(KEEP_KEY) !== '0'; } catch (error) { /* keeping stays on */ }
  if (!storeAvailable()) { await showSaved(); return; }
  try {
    const found = await readStore(['header']);
    if (found.header) { savedHeader = found.header; saveHeld = true; }
  } catch (error) {
    show('saved-status', `Saving in the browser is not available here (${error.message || error}). Use Export to keep the results.`);
    for (const id of ['saved-restore', 'saved-delete', 'saved-save-now']) $(id).hidden = true;
    $('keep-results').disabled = true;
    return;
  }
  await showSaved();
})();

newSample();
renderSamples();

describeSources();
