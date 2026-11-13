// site/gz.js: a streaming gzip reader that checks the CRC-32 and the length of every member (fflate does not).
import test from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync, deflateRawSync } from 'node:zlib';
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

// A chunk can be larger than the tail that the reader keeps, and a member can start anywhere in a chunk.
function bigMembers() {
  const parts = [text(20000, 1), text(20000, 2), text(5, 3)];
  const members = parts.map((p) => gzipSync(p, { level: 0 })); // stored: members of about 500 KB
  return { gz: Buffer.concat(members), data: Buffer.concat(parts), firstEnd: members[0].length, secondEnd: members[0].length + members[1].length };
}

test('chunks larger than the tail window are fine, up to the whole file in one push', () => {
  const { gz, data } = bigMembers();
  for (const size of [65536, 70000, 262144, 1 << 20, gz.length]) assert.ok(run(gz, size).equals(data), `chunk size ${size}`);
});

test('a member that starts within 12 bytes of a chunk edge, with chunks above 64 KiB, is fine', () => {
  const { gz, data, firstEnd, secondEnd } = bigMembers();
  for (const end of [firstEnd, secondEnd]) {
    for (let d = -12; d <= 12; d++) assert.ok(run(gz, end + d).equals(data), `chunk size ${end + d} (member boundary at ${end})`);
  }
});

/** A gzip member of exactly `length` bytes that holds `data` (the header name field pads it). */
function memberOfLength(data, length) {
  const raw = deflateRawSync(data, { level: 0 });
  const trailer = gzipSync(data).subarray(-8); // CRC-32 and length of `data`
  const name = length - (10 + 1 + raw.length + 8);
  assert.ok(name >= 0, 'the data is too long for this length');
  const header = Buffer.concat([Buffer.from([0x1f, 0x8b, 8, 8, 0, 0, 0, 0, 0, 3]), Buffer.alloc(name, 0x61), Buffer.from([0])]);
  return Buffer.concat([header, raw, trailer]);
}

test('chunks of exactly 64 KiB, with the end of a member 0 to 24 bytes before or after a chunk edge, are fine', () => {
  const first = Buffer.alloc(100000, 0x41);
  const second = text(300, 2);
  const tail = gzipSync(second);
  const data = Buffer.concat([first, second]);
  for (const edge of [65536 * 2, 65536 * 3]) {
    for (let d = -24; d <= 24; d++) {
      const gz = Buffer.concat([memberOfLength(first, edge + d), tail]);
      assert.ok(run(gz, 65536).equals(data), `first member ends ${d} bytes from the edge at ${edge}`);
    }
  }
});

test('an empty file is an error that says the file is empty, not a header error', () => {
  assert.throws(() => run(Buffer.alloc(0), 10), /empty/i);
});

test('bytes after the last member give an error that mentions them', () => {
  const gz = Buffer.concat([gzipSync(text(500, 1)), Buffer.alloc(2)]);
  assert.throws(() => run(gz, 100), /after the last member|extra|padding/i);
});
