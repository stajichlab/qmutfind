// A module worker: it runs qmutfind on the files the page gives it, in WebAssembly, off the main thread.
// It answers four kinds of message: run, summarize, annotate and report. Every answer has the same `step`.
// Build the package first with `pixi run build-wasm-web`.
import init, { Session, annotate, stats, report, genome_plot, depth_plot, combine_tables, version } from './pkg-web/qmutfind_wasm.js?v=__BUILD__';

// fflate (MIT, vendor/fflate.LICENSE) reads gzip files with several members. The built-in
// DecompressionStream stops with an error on them in Firefox, and tools such as parallel-fastq-dump write them.
// gz.js adds the check of the CRC-32 and the length of every member, which fflate does not make.
import { VerifyingGunzip } from './gz.js?v=__BUILD__';
import { keepFirstReads } from './reads.js?v=__BUILD__';

// The wasm file is named with the build stamp too, so a cached old file cannot meet a new glue file.
const start = () => init({ module_or_path: new URL('./pkg-web/qmutfind_wasm_bg.wasm?v=__BUILD__', import.meta.url) });

/** Yield the bytes of a file in chunks. A .gz file is decompressed on the way. */
async function* chunksOf(file) {
  const reader = file.stream().getReader();
  if (!file.name.endsWith('.gz')) {
    for (;;) {
      const { done, value } = await reader.read();
      if (done) return;
      yield value;
    }
  }
  const pending = [];
  const gunzip = new VerifyingGunzip((chunk) => pending.push(chunk));
  for (;;) {
    const { done, value } = await reader.read();
    gunzip.push(value || new Uint8Array(0), done);
    while (pending.length) yield pending.shift();
    if (done) return;
  }
}

// A file that is loaded whole into memory (reference, genes, VCF, coverage) may unpack to at most this many bytes. A
// small gzip file can unpack to many gigabytes, which would stop the tab. The page sends the limit with each message.
let maxBytes = 1e9;

async function readAll(file) {
  const parts = [];
  let size = 0;
  for await (const chunk of chunksOf(file)) {
    parts.push(chunk);
    size += chunk.length;
    if (size > maxBytes) {
      throw new Error(`${file.name} is too large: it unpacks to more than ${Math.round(maxBytes / 1e6)} MB, which the page does not load into memory.`);
    }
  }
  const out = new Uint8Array(size);
  let offset = 0;
  for (const part of parts) {
    out.set(part, offset);
    offset += part.length;
  }
  return out;
}

// The last reference and gene files are kept, because steps 2 and 3 use them again.
const kept = new Map();
async function bytesOf(name, file) {
  const have = kept.get(name);
  if (have && have.file === file) return have.bytes;
  const bytes = await readAll(file);
  kept.set(name, { file, bytes });
  return bytes;
}

async function vcfText(source) {
  return typeof source === 'string' ? source : new TextDecoder().decode(await readAll(source));
}

async function run({ reference, reads, sample, kind, depthWindow, evidence, maxReads }) {
  const started = performance.now();
  postMessage({ step: 'run', type: 'status', text: 'Reading the reference and building the index' });
  const fasta = await bytesOf('reference', reference);
  const session = new Session(fasta, kind, sample);
  let reads_mapped, vcf, depth, evidenceVcf;
  const fileReports = []; // per file: the reads taken and whether the read limit cut it
  try {
    let lastReport = 0;
    for (const file of reads) {
      postMessage({ step: 'run', type: 'status', text: `Mapping ${file.name}` });
      // With a read limit, only the first `maxReads` reads of each file are mapped: four lines are one read.
      const lineLimit = maxReads ? maxReads * 4 : 0;
      let lines = 0;
      let cut = false;
      for await (const chunk of chunksOf(file)) {
        let bytes = chunk;
        if (lineLimit) {
          let seen = 0;
          for (let i = 0; i < bytes.length; i++) {
            if (bytes[i] === 10 && ++seen + lines === lineLimit) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
          }
          lines += seen;
        }
        session.add_reads(bytes);
        if (performance.now() - lastReport > 1000) {
          lastReport = performance.now();
          postMessage({ step: 'run', type: 'progress', reads: session.read_count() });
        }
        if (cut) break;
      }
      session.end_file();
      fileReports.push({ name: file.name, cut });
    }
    reads_mapped = session.read_count();
    vcf = session.vcf();
    depth = session.depth_bed(depthWindow || 10000);
    if (evidence) evidenceVcf = session.evidence_vcf();
  } finally {
    session.free(); // also after an error: a failed run must not keep its index and counts (about 460 MB for Af293)
  }
  postMessage({ step: 'run', type: 'done', vcf, depth, evidence: evidenceVcf, version: version(), files: fileReports, reads: reads_mapped, seconds: (performance.now() - started) / 1000 });
}

/**
 * Download a FASTQ file (plain or gzip) and keep only its first `maxReads` reads, then stop the download (see reads.js). This runs here,
 * not on the main thread, because unpacking and packing take seconds for a million reads. Answers with { bytes, reads, cut }.
 * `maxKept` is the most packed bytes allowed and `maxDownload` the most bytes taken from the server.
 */
async function fetchReads({ url, maxReads, maxKept, maxDownload, label }) {
  let address;
  try { address = new URL(url); } catch (error) { throw new Error(`${label}: "${url}" is not an address.`); }
  if (address.protocol !== 'https:' && address.origin !== self.location.origin) throw new Error(`${label}: only https addresses are accepted.`);
  if (address.username || address.password) throw new Error(`${label}: addresses with a user name or password are not accepted.`);
  const response = await fetch(address);
  if (!response.ok) throw new Error(`${label}: the server answered ${response.status} for ${address.host}`);
  const result = await keepFirstReads(response.body.getReader(), {
    maxReads, maxKept, maxDownload, label,
    onProgress: (p) => postMessage({ step: 'fetchreads', type: 'progress', ...p }),
  });
  postMessage({ step: 'fetchreads', type: 'done', ...result }, [result.bytes.buffer]);
}

async function summarize({ vcf, bed, reference, order }) {
  postMessage({ step: 'summarize', type: 'status', text: 'Summarizing' });
  const text = vcf ? await vcfText(vcf) : null;
  const fasta = reference ? await bytesOf('reference', reference) : undefined;
  const answer = { step: 'summarize', type: 'done', stats: null, plot: null, depthPlot: null };
  if (text !== null) {
    answer.stats = stats(text, 100);
    answer.plot = genome_plot(text, fasta, order);
  }
  if (bed) {
    const bedText = typeof bed === 'string' ? bed : new TextDecoder().decode(await readAll(bed));
    answer.depthPlot = depth_plot(bedText, text === null ? undefined : text, fasta, order);
  }
  postMessage(answer);
}

async function annotateStep({ vcf, reference, gff }) {
  postMessage({ step: 'annotate', type: 'status', text: 'Adding the effect of each variant' });
  const text = await vcfText(vcf);
  const fasta = await bytesOf('reference', reference);
  const gffBytes = await bytesOf('gff', gff);
  postMessage({ step: 'annotate', type: 'done', vcf: annotate(text, fasta, gffBytes) });
}

async function reportStep({ vcf, gff, region, title }) {
  const text = await vcfText(vcf);
  const gffBytes = gff ? await bytesOf('gff', gff) : undefined;
  postMessage({ step: 'report', type: 'done', page: report(text, gffBytes, region, title || `Mutations in ${region}`), region });
}

// The size of the WebAssembly memory, for the tests (it can grow but never shrinks).
async function memory() {
  const exports = await start();
  postMessage({ step: 'memory', type: 'done', bytes: exports.memory.buffer.byteLength });
}

// The combined table of several strains (see `combine_tables`).
async function combineStep({ names, vcfs, beds, evidence, parent, minDepth }) {
  const [table, summary, warnings, vcf] = combine_tables(names, vcfs, beds, evidence, parent, minDepth);
  postMessage({ step: 'combine', type: 'done', table, summary, warnings, vcf });
}

const handlers = { run, fetchreads: fetchReads, summarize, annotate: annotateStep, report: reportStep, memory, combine: combineStep };

self.onmessage = async (event) => {
  const { step } = event.data;
  if (event.data.limit) maxBytes = event.data.limit;
  try {
    if (step !== 'fetchreads') await start(); // the download step needs no WebAssembly
    await handlers[step](event.data);
  } catch (error) {
    postMessage({ step, type: 'error', message: String((error && error.message) || error) });
  }
};
