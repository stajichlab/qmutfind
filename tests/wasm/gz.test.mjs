// site/gz.js: a streaming gzip reader that checks the CRC-32 and the length of every member (fflate does not).
import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const { VerifyingGunzip, crc32 } = await import(pathToFileURL(join(process.cwd(), 'site', 'gz.js')).href);

function text(n, seed) {
  let s = '';
  for (let i = 0; i < n; i++) s += `line ${seed} ${i} ${(i * 2654435761 % 97)}\n`;
  return Buffer.from(s);
}

/** Feed `bytes` in pieces of `size` and return what comes out. A pushed final chunk is empty, as in the page. */
function run(bytes, size) {
  const out = [];
  const g = new VerifyingGunzip((chunk) => out.push(Buffer.from(chunk)));
  for (let i = 0; i < bytes.length; i += size) g.push(bytes.subarray(i, i + size), false);
  g.push(new Uint8Array(0), true);
  return Buffer.concat(out);
}

test('crc32 of the standard check string, in one piece and in two', () => {
  const data = Buffer.from('123456789');
  assert.equal(crc32(data), 0xcbf43926);
  assert.equal(crc32(data.subarray(4), crc32(data.subarray(0, 4))), 0xcbf43926);
});

test('a good file with one member comes out whole for every chunk size', () => {
  const data = text(5000, 1);
  const gz = gzipSync(data);
  for (const size of [1, 7, 100, 4096, gz.length]) assert.ok(run(gz, size).equals(data), `chunk size ${size}`);
});

test('a good file with several members comes out whole for every chunk size', () => {
  const parts = [text(3000, 1), text(10, 2), text(4000, 3)];
  const gz = Buffer.concat(parts.map((p) => gzipSync(p)));
  const data = Buffer.concat(parts);
  for (const size of [1, 5, 63, 1000, gz.length]) assert.ok(run(gz, size).equals(data), `chunk size ${size}`);
});

test('a changed data byte in a stored block is an error that names the checksum', () => {
  const data = text(2000, 1);
  const gz = Buffer.from(gzipSync(data, { level: 0 })); // stored: the inflater cannot notice a changed byte
  gz[gz.length - 8 - 100] ^= 0x01;
  assert.throws(() => run(gz, 513), /checksum/i);
});

test('a wrong CRC in the trailer of the last member is an error', () => {
  const gz = Buffer.from(gzipSync(text(2000, 1)));
  gz[gz.length - 8] ^= 0xff;
  assert.throws(() => run(gz, 100), /checksum/i);
});

test('a wrong CRC in a member before the last is an error, even with a good last member', () => {
  const a = Buffer.from(gzipSync(text(2000, 1)));
  a[a.length - 6] ^= 0x10;
  const gz = Buffer.concat([a, gzipSync(text(500, 2))]);
  for (const size of [1, 64, gz.length]) assert.throws(() => run(gz, size), /checksum/i, `chunk size ${size}`);
});

test('a wrong length in the trailer is an error', () => {
  const gz = Buffer.from(gzipSync(text(2000, 1)));
  gz[gz.length - 1] ^= 0x01;
  assert.throws(() => run(gz, 100), /length|size|checksum/i);
});

test('a file cut inside the data is an error', () => {
  const gz = gzipSync(text(20000, 1));
  assert.throws(() => run(gz.subarray(0, gz.length - 200), 1000));
});

test('a file with its trailer cut off is an error', () => {
  const gz = gzipSync(text(2000, 1));
  assert.throws(() => run(gz.subarray(0, gz.length - 5), 1000));
});

test('output of a large file with a good checksum is not delayed to the end', () => {
  const data = text(200000, 1); // about 5 MB
  const gz = gzipSync(data);
  const seen = [];
  const g = new VerifyingGunzip((chunk) => seen.push(chunk.length));
  g.push(gz.subarray(0, Math.floor(gz.length / 2)), false);
  assert.ok(seen.length > 0, 'data comes out while the file is still arriving');
  g.push(gz.subarray(Math.floor(gz.length / 2)), false);
  g.push(new Uint8Array(0), true);
  assert.equal(seen.reduce((a, b) => a + b, 0), data.length);
});
