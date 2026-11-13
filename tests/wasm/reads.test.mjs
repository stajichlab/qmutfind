// site/reads.js: keep the first reads of a FASTQ stream (plain or gzip), packed as gzip.
import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync, gunzipSync } from 'node:zlib';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const { keepFirstReads } = await import(pathToFileURL(join(process.cwd(), 'site', 'reads.js')).href);

const fastq = (n, prefix = 'r') => Array.from({ length: n }, (_, i) => `@${prefix}${i}\nACGTACGTAC\n+\nIIIIIIIIII\n`).join('');

/** A reader over `bytes` in pieces of the given sizes (the last size repeats). It records whether it was cancelled. */
function readerOf(bytes, sizes) {
  const buf = Buffer.from(bytes);
  let at = 0, k = 0;
  const reader = {
    cancelled: false,
    async read() {
      if (at >= buf.length) return { done: true, value: undefined };
      const size = sizes[Math.min(k++, sizes.length - 1)];
      const value = new Uint8Array(buf.subarray(at, at + size));
      at += size;
      return { done: false, value };
    },
    async cancel() { reader.cancelled = true; },
  };
  return reader;
}
const options = (extra = {}) => ({ maxReads: 1000, maxKept: 1e9, maxDownload: 1e9, label: 'test file', onProgress: () => {}, ...extra });
const unpack = (bytes) => gunzipSync(Buffer.from(bytes)).toString();

test('a plain file is cut at the limit and the download is cancelled', async () => {
  const reader = readerOf(fastq(100), [50]);
  const r = await keepFirstReads(reader, options({ maxReads: 7 }));
  assert.equal(r.reads, 7); assert.equal(r.cut, true); assert.equal(reader.cancelled, true);
  assert.equal(unpack(r.bytes), fastq(7));
});

test('a gzip file is cut at the limit', async () => {
  const r = await keepFirstReads(readerOf(gzipSync(fastq(100)), [64]), options({ maxReads: 7 }));
  assert.equal(r.reads, 7); assert.equal(r.cut, true);
  assert.equal(unpack(r.bytes), fastq(7));
});

test('a limit above the number of reads keeps them all and does not cut', async () => {
  const reader = readerOf(gzipSync(fastq(10)), [64]);
  const r = await keepFirstReads(reader, options({ maxReads: 100 }));
  assert.equal(r.reads, 10); assert.equal(r.cut, false); assert.equal(reader.cancelled, false);
  assert.equal(unpack(r.bytes), fastq(10));
});

test('a gzip file whose first chunk is one byte is still read as gzip', async () => {
  const r = await keepFirstReads(readerOf(gzipSync(fastq(20)), [1, 1, 3, 50]), options({ maxReads: 5 }));
  assert.equal(r.reads, 5);
  assert.equal(unpack(r.bytes), fastq(5));
});

test('a damaged gzip file read to its end is an error that names the checksum', async () => {
  const gz = Buffer.from(gzipSync(fastq(200), { level: 0 }));
  gz[gz.length >> 1] ^= 1;
  await assert.rejects(keepFirstReads(readerOf(gz, [100]), options({ maxReads: 100000 })), /checksum/i);
});

test('damage after the cut is never read', async () => {
  const gz = Buffer.concat([gzipSync(fastq(50)), Buffer.alloc(5000, 7)]);
  const r = await keepFirstReads(readerOf(gz, [100]), options({ maxReads: 3 }));
  assert.equal(r.reads, 3); assert.equal(r.cut, true);
});

test('kept reads above the limit are an error, also when only the last flush crosses it', async () => {
  let seed = 12345; const rnd = () => ((seed = (Math.imul(seed, 1103515245) + 12345) >>> 0) >>> 16) % 4;
  const body = Array.from({ length: 200 }, (_, i) => `@q${i}\n${Array.from({ length: 100 }, () => 'ACGT'[rnd()]).join('')}\n+\n${'I'.repeat(100)}\n`).join('');
  const whole = await keepFirstReads(readerOf(body, [500]), options({ maxReads: 1000 }));
  const size = whole.bytes.length;
  assert.ok(size > 2000, `the test data must not pack small (${size} bytes)`);
  // The last bytes come only with the final flush, so one byte below the full size can only be caught there.
  for (const maxKept of [50, size >> 1, size - 1]) {
    await assert.rejects(keepFirstReads(readerOf(body, [500]), options({ maxReads: 1000, maxKept })), /too large/i, `maxKept ${maxKept}`);
  }
  const ok = await keepFirstReads(readerOf(body, [500]), options({ maxReads: 1000, maxKept: size }));
  assert.equal(ok.bytes.length, size);
});

test('the download itself is bounded, also for input that packs very small', async () => {
  const body = 'A'.repeat(5_000_000); // no line breaks: nothing is ever kept
  await assert.rejects(keepFirstReads(readerOf(body, [64 * 1024]), options({ maxReads: 1000, maxDownload: 1_000_000 })), /too large/i);
});

test('progress is reported with the bytes downloaded and the reads kept', async () => {
  const seen = [];
  const reader = readerOf(fastq(5000), [4096]);
  await keepFirstReads(reader, options({ maxReads: 100000, onProgress: (p) => seen.push(p), progressEvery: 0 }));
  assert.ok(seen.length > 0 && seen.at(-1).reads > 0 && seen.at(-1).downloaded > 0, JSON.stringify(seen.at(-1)));
});
