// A module worker: it runs qmutfind on the files the page gives it, in WebAssembly, off the main thread.
// Build the package first with `pixi run build-wasm-web`.
import init, { Session, annotate, stats, report } from '../../pkg-web/qmutfind_wasm.js';

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

async function readAll(file) {
  const parts = [];
  let size = 0;
  for await (const chunk of chunksOf(file)) {
    parts.push(chunk);
    size += chunk.length;
  }
  const out = new Uint8Array(size);
  let offset = 0;
  for (const part of parts) {
    out.set(part, offset);
    offset += part.length;
  }
  return out;
}

self.onmessage = async (event) => {
  try {
    const { reference, reads, gff, sample, kind, region } = event.data;
    await init();
    const started = performance.now();

    postMessage({ type: 'status', text: 'Reading the reference and building the index' });
    const fasta = await readAll(reference);
    const session = new Session(fasta, kind, sample);

    let lastReport = 0;
    for (const file of reads) {
      postMessage({ type: 'status', text: `Mapping ${file.name}` });
      for await (const chunk of chunksOf(file)) {
        session.add_reads(chunk);
        if (performance.now() - lastReport > 1000) {
          lastReport = performance.now();
          postMessage({ type: 'progress', reads: session.read_count() });
        }
      }
      session.end_file();
    }
    const reads_mapped = session.read_count();
    let vcf = session.vcf();
    session.free();

    let gffBytes;
    if (gff) {
      postMessage({ type: 'status', text: 'Adding the effect of each variant' });
      gffBytes = await readAll(gff);
      vcf = annotate(vcf, fasta, gffBytes);
    }
    const statistics = stats(vcf, 100);
    const page = region ? report(vcf, gffBytes, region, `Mutations in ${region}`) : null;
    postMessage({
      type: 'done',
      vcf,
      stats: statistics,
      report: page,
      reads: reads_mapped,
      seconds: (performance.now() - started) / 1000,
    });
  } catch (error) {
    postMessage({ type: 'error', message: String((error && error.message) || error) });
  }
};
