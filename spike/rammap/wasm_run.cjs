// Stream a reference and a FASTQ through rammap's AlignSession (WASM) in Node.
// Usage: node wasm_run.cjs <pkg_dir> <ref.fa|.gz> <reads.fq|.gz> <preset> <out.sam>
const fs = require('fs');
const zlib = require('zlib');
const path = require('path');

const [pkgDir, refPath, readsPath, preset, outPath] = process.argv.slice(2);
const mod = require(path.join(pkgDir, 'rammap.js'));
if (mod._set_panic_hook) mod._set_panic_hook();

function openStream(p) {
  const s = fs.createReadStream(p, { highWaterMark: 4 << 20 });
  return p.endsWith('.gz') ? s.pipe(zlib.createGunzip()) : s;
}

(async () => {
  const t0 = process.hrtime.bigint();
  const session = new mod.AlignSession(preset, true, true);
  let refBytes = 0;
  for await (const chunk of openStream(refPath)) {
    refBytes += chunk.length;
    session.append_ref(new Uint8Array(chunk.buffer, chunk.byteOffset, chunk.length));
  }
  session.finalize_ref();
  const t1 = process.hrtime.bigint();
  const out = fs.createWriteStream(outPath);
  let readBytes = 0;
  let tail = '';
  for await (const chunk of openStream(readsPath)) {
    readBytes += chunk.length;
    const s = session.append_query(new Uint8Array(chunk.buffer, chunk.byteOffset, chunk.length));
    if (s.length) out.write(s);
  }
  const fin = session.finalize();
  const t2 = process.hrtime.bigint();
  const idx = fin.indexOf('---LOG---');
  if (idx > 0) out.write(fin.slice(0, idx));
  tail = idx >= 0 ? fin.slice(idx) : '';
  await new Promise((r) => out.end(r));
  const ms = (a, b) => Number(b - a) / 1e6;
  console.log(JSON.stringify({
    index_ms: Math.round(ms(t0, t1)),
    align_ms: Math.round(ms(t1, t2)),
    ref_bytes: refBytes,
    read_bytes: readBytes,
    rss_mb: Math.round(process.memoryUsage().rss / 1048576),
  }));
  console.log(tail.split('\n').slice(0, 8).join('\n'));
})().catch((e) => { console.error('FAILED', e); process.exit(1); });
