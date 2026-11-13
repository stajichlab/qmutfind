// What the main page and the allele page share: the checks on page addresses, the download of a file, and the reference panel (the NCBI list and the
// FungiDB command, files, and web addresses). Each page gives the panel its own state object and a few hooks, and keeps nothing else of it.

/** The most that one download or one set of downloads may take: 2 GB. `?download_limit_mb=N` can lower the limit (the tests use that). */
export const MAX_DOWNLOAD = Math.min(parseFloat(new URLSearchParams(location.search).get('download_limit_mb')) || 2000, 2000) * 1e6;

export const mb = (bytes) => `${(bytes / 1e6).toFixed(1)} MB`;

/**
 * A page address option (?example=, ?genomes=, ?enaapi=) may only name a path on this site. The text is refused when it has a
 * scheme or starts with / or \ (browsers read \ as /, so \\host/x would name another host), and the resolved address must
 * have this page's origin. Returns the text, or `fallback`.
 */
export function sitePath(asked, fallback) {
  if (!asked || asked !== asked.trim() || /[\u0000-\u001f]/.test(asked)) return fallback; // the URL parser would drop these characters
  if (/^[\\/]/.test(asked) || /^[a-z][a-z0-9+.-]*:/i.test(asked)) return fallback;
  try {
    return new URL(asked, location.href).origin === location.origin ? asked : fallback;
  } catch (error) {
    return fallback;
  }
}

/** An address the page may download from: https (or this site), and no user name or password. Returns the URL, or throws an Error that names `label`. */
export function checkUrl(text, label) {
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
export function directAddress(text) {
  let url;
  try { url = new URL(text); } catch (error) { return text; }
  const shared = /^\/(s|scl\/fi)\//.test(url.pathname);
  if (url.protocol !== 'https:' || !['www.dropbox.com', 'dropbox.com'].includes(url.hostname) || !shared) return text;
  url.hostname = 'dl.dropboxusercontent.com';
  url.searchParams.delete('dl');
  url.searchParams.delete('raw');
  return url.href;
}

export const hostOf = (item, base) => new URL(item.url, base).host;

/**
 * Download a file into memory and return it as a File. `item` is { url, name?, bytes? } (`url` can be relative to `base`; `bytes` is the size that the file must
 * have). `budget` is { left }: the bytes that the downloads of one task may still take; the file's bytes are taken off it. `onProgress(text)` gets a line for the
 * status display. Throws an Error that names `label`.
 */
export async function download(item, { label, base = location.href, signal, budget = { left: MAX_DOWNLOAD }, onProgress = () => {} }) {
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
    onProgress(`Downloading ${label}: ${mb(got)}${expected ? ` of ${mb(expected)}` : ''}`);
  }
  if (item.bytes && got !== item.bytes) throw new Error(`${label}: got ${got} bytes, expected ${item.bytes}. The download is incomplete.`);
  budget.left -= got;
  // The page reads a file as gzip by its name. A download whose name has no .gz but whose content is gzip gets the name.
  let name = item.name || decodeURIComponent(url.pathname.split('/').filter(Boolean).pop() || 'download');
  if (parts.length && parts[0].length > 1 && parts[0][0] === 0x1f && parts[0][1] === 0x8b && !name.endsWith('.gz')) name += '.gz';
  return new File(parts, name);
}

const LARGE_GENOME = 50e6;

/**
 * The reference panel: the NCBI list with its filter and the FungiDB command, the upload of a FASTA and a GFF3, and the download from https addresses. The page needs
 * the elements of the panel (`reference-source`, `genome-text`, `genome-filter`, `genome-menu`, `genome-load`, `genome-fungidb`, `genome-command`, `reference`, `gff`,
 * `url-reference`, `url-gff`, `ref-url-load`, `genome-status`, `source-reference`, `source-gff`, and the three groups `genome-section`, `reference-upload` and `reference-url`).
 *
 * - `state`: { reference, gff, source: { reference, gff } }. The panel puts the files that it loads here; the page reads them (a file chosen in an input wins).
 * - `canLoad()`: false while the page is busy with something else.
 * - `setLoading(flag)`: the panel is loading (true) or has finished (false).
 * - `onLoaded()`: called after a genome has loaded, before the source lines are drawn.
 * - `texts`: { noGenes, next }: what the page adds to the "Loaded" line for a genome without gene models, and the sentence that closes it.
 * Returns { describe }, which draws the source lines (call it after the page changes `state` itself).
 */
export function createReferencePanel({ $, show, state, canLoad, setLoading, onLoaded, texts }) {
  let genomeList = [];
  let genomesBase = location.href; // file addresses in the list are relative to the list
  let fungidbList = [];

  /** Say, under each input, where its data comes from: a local file, or a download (an example, the menu, the web addresses or a saved copy). */
  function describe() {
    const size = (f) => `${(f.size / 1e6).toFixed(f.size < 1e5 ? 3 : 1)} MB`;
    const line = (slot, inputId, file) => {
      const chosen = [...$(inputId).files];
      const element = $(`source-${slot}`);
      element.classList.remove('from-web');
      if (chosen.length) {
        element.textContent = `Local file: ${chosen.map((f) => `${f.name} (${size(f)})`).join(', ')}`;
      } else if (file) {
        element.textContent = `From ${state.source[slot]}: ${file.name} (${size(file)})`;
        element.classList.add('from-web');
      } else {
        element.textContent = 'Not set. Use the NCBI list, upload files or give addresses above.';
      }
    };
    line('reference', 'reference', state.reference);
    line('gff', 'gff', state.gff);
  }
  for (const id of ['reference', 'gff']) $(id).addEventListener('change', describe);

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
    if (!canLoad()) return;
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
    setLoading(true);
    $('genome-load').disabled = true;
    try {
      const base = genomesBase;
      const progress = (text) => show('genome-status', text);
      const reference = await download(entry.fasta, { label: 'the genome', base, onProgress: progress });
      const gff = entry.gff ? await download(entry.gff, { label: 'the gene models', base, onProgress: progress }) : null;
      state.reference = reference;
      state.gff = gff;
      state.source = {
        reference: `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.fasta, base)})`,
        gff: entry.gff ? `the genome menu (${entry.accession}, ${entry.label}, at ${hostOf(entry.gff, base)})` : '',
        reads: '',
      };
      onLoaded();
      describe();
      let note = entry.gff ? '' : texts.noGenes;
      note += texts.next;
      if (entry.genome_size > LARGE_GENOME) {
        note += ` This genome is large (${(entry.genome_size / 1e6).toFixed(0)} Mb). Genomes above 50 Mb have not been tested in a browser and may need several GB of memory.`;
      }
      show('genome-status', `Loaded: ${entry.label} (${entry.accession}).${note}`);
    } catch (error) {
      show('genome-status', `Error: ${error.message} The server may block downloads from web pages, or the network failed.`);
    } finally {
      setLoading(false);
      $('genome-load').disabled = false;
    }
  }

  // genomes/genomes.json lists NCBI genomes (RefSeq, and some GenBank without gene models) (the page downloads them: NCBI allows that). genomes/fungidb.json lists
  // FungiDB organisms. FungiDB blocks downloads from web pages and needs an API key, so for those the page only shows the `qmutfind fetch` command.
  // ?genomes=PATH uses another list, but only a path on this site.
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
    if (!canLoad()) return;
    const fastaText = $('url-reference').value.trim();
    const gffText = $('url-gff').value.trim();
    if (!fastaText) { show('genome-status', 'Type the address of the reference FASTA.'); return; }
    setLoading(true);
    $('ref-url-load').disabled = true;
    try {
      const wanted = [['reference', 'the reference genome', fastaText], ...(gffText ? [['gff', 'the gene file', gffText]] : [])];
      let changed = false;
      const urls = wanted.map(([, label, text]) => { const d = directAddress(text); changed ||= d !== text; return checkUrl(d, label); }); // all checked before any download
      const got = {};
      const source = { reads: '' };
      for (const [i, [slot, label]] of wanted.entries()) {
        const url = urls[i];
        got[slot] = await download({ url: url.href }, { label, onProgress: (text) => show('genome-status', text) });
        source[slot] = `the address ${url.host}${url.pathname}`;
      }
      state.reference = got.reference;
      state.gff = got.gff || null;
      state.source = { reference: source.reference, gff: got.gff ? source.gff : '', reads: '' };
      onLoaded();
      describe();
      show('genome-status', `Loaded: the reference genome${got.gff ? ' and the gene file' : ''}.${changed ? ' A Dropbox share link was changed to dl.dropboxusercontent.com, the address that allows downloads from web pages.' : ''}`);
    } catch (error) {
      show('genome-status', `Error: ${error.message}`);
    } finally {
      setLoading(false);
      $('ref-url-load').disabled = false;
    }
  });

  describe();
  return { describe };
}
