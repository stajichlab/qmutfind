// Run the WebAssembly build on gzip or plain files from Node and report time and memory.
//   node scripts/wasm-run.mjs --reference ref.fa[.gz] --reads r1.fq[.gz] [--reads r2.fq[.gz]]
//        --sample NAME --kind short|long --out out.vcf [--pkg pkg-node]
//        [--min-depth N --min-snp-fraction X --min-indel-fraction X]
// The last line printed is a JSON object with the read count, the times and the memory in use.
import { createRequire } from 'node:module';
import { createReadStream, readFileSync, writeFileSync } from 'node:fs';
import { createGunzip, gunzipSync } from 'node:zlib';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

function parseArgs(argv) {
  const opts = { reads: [] };
  for (let i = 0; i < argv.length; i += 2) {
    const key = argv[i];
    if (!key.startsWith('--') || i + 1 >= argv.length) throw new Error(`bad argument: ${key}`);
    const name = key.slice(2);
    if (name === 'reads') opts.reads.push(argv[i + 1]);
    else opts[name] = argv[i + 1];
  }
  for (const required of ['reference', 'sample', 'kind', 'out']) {
    if (!opts[required]) throw new Error(`--${required} is required`);
  }
  if (opts.reads.length === 0) throw new Error('--reads is required');
  return opts;
}

function openStream(path) {
  const stream = createReadStream(path, { highWaterMark: 4 << 20 });
  return path.endsWith('.gz') ? stream.pipe(createGunzip()) : stream;
}

try {
  const opts = parseArgs(process.argv.slice(2));
  const require = createRequire(import.meta.url);
  // The default package sits next to the scripts folder, wherever the script is run from.
  const pkgDir = opts.pkg ? resolve(opts.pkg) : fileURLToPath(new URL('../pkg-node', import.meta.url));
  const wasm = require(resolve(pkgDir, 'qmutfind_wasm.js'));
  const mb = (bytes) => Math.round(bytes / 1048576);

  const t0 = process.hrtime.bigint();
  const raw = readFileSync(opts.reference);
  const fasta = opts.reference.endsWith('.gz') ? gunzipSync(raw) : raw;
  const session = new wasm.Session(new Uint8Array(fasta.buffer, fasta.byteOffset, fasta.length), opts.kind, opts.sample);
  if (opts['min-depth'] || opts['min-snp-fraction'] || opts['min-indel-fraction']) {
    session.set_thresholds(Number(opts['min-depth'] ?? 5), Number(opts['min-snp-fraction'] ?? 0.8),
      Number(opts['min-indel-fraction'] ?? 0.5));
  }
  const t1 = process.hrtime.bigint();

  for (const path of opts.reads) {
    for await (const chunk of openStream(path)) {
      session.add_reads(new Uint8Array(chunk.buffer, chunk.byteOffset, chunk.length));
    }
    session.end_file();
  }
  const t2 = process.hrtime.bigint();

  writeFileSync(opts.out, session.vcf());
  const reads = session.read_count();
  const wasmBytes = wasm.__wasm?.memory?.buffer?.byteLength ?? null;
  session.free();
  console.log(JSON.stringify({
    reads,
    seconds_index: Number(t1 - t0) / 1e9,
    seconds_map: Number(t2 - t1) / 1e9,
    seconds_total: Number(process.hrtime.bigint() - t0) / 1e9,
    rss_mb: mb(process.memoryUsage().rss),
    wasm_memory_mb: wasmBytes === null ? null : mb(wasmBytes),
  }));
} catch (e) {
  console.error(`wasm-run: ${e.message}`);
  process.exit(1);
}
