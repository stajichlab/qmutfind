// A module worker: it runs qmutfind on the files the page gives it, in WebAssembly, off the main thread.
// It answers four kinds of message: run, summarize, annotate and report. Every answer has the same `step`.
// Build the package first with `pixi run build-wasm-web`.
import init, { Session, annotate, stats, report, genome_plot, depth_plot } from './pkg-web/qmutfind_wasm.js';

// fflate (MIT, vendor/fflate.LICENSE) reads gzip files with several members. The built-in
// DecompressionStream stops with an error on them in Firefox, and tools such as parallel-fastq-dump write them.
import { Gunzip } from './vendor/fflate.js';

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
  const gunzip = new Gunzip((chunk) => pending.push(chunk));
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

async function run({ reference, reads, sample, kind }) {
  const started = performance.now();
  postMessage({ step: 'run', type: 'status', text: 'Reading the reference and building the index' });
  const fasta = await bytesOf('reference', reference);
  const session = new Session(fasta, kind, sample);
  let lastReport = 0;
  for (const file of reads) {
    postMessage({ step: 'run', type: 'status', text: `Mapping ${file.name}` });
    for await (const chunk of chunksOf(file)) {
      session.add_reads(chunk);
      if (performance.now() - lastReport > 1000) {
        lastReport = performance.now();
        postMessage({ step: 'run', type: 'progress', reads: session.read_count() });
      }
    }
    session.end_file();
  }
  const reads_mapped = session.read_count();
  const vcf = session.vcf();
  const depth = session.depth_bed(10000);
  session.free();
  postMessage({ step: 'run', type: 'done', vcf, depth, reads: reads_mapped, seconds: (performance.now() - started) / 1000 });
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

async function reportStep({ vcf, gff, region }) {
  const text = await vcfText(vcf);
  const gffBytes = gff ? await bytesOf('gff', gff) : undefined;
  postMessage({ step: 'report', type: 'done', page: report(text, gffBytes, region, `Mutations in ${region}`), region });
}

const handlers = { run, summarize, annotate: annotateStep, report: reportStep };

self.onmessage = async (event) => {
  const { step } = event.data;
  if (event.data.limit) maxBytes = event.data.limit;
  try {
    await init();
    await handlers[step](event.data);
  } catch (error) {
    postMessage({ step, type: 'error', message: String((error && error.message) || error) });
  }
};
