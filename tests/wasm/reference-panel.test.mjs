// site/reference-panel.js: what the main page and the allele page share: address checks, the download, and the reference panel. The panel itself
// (the menu, the files, the addresses) needs a page, so the browser tests cover it on both pages; this file tests the helpers.
import test from 'node:test';
import assert from 'node:assert/strict';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

globalThis.location = new URL('https://qmutfind.example.org/dir/page.html?x=1');
const MODULE = pathToFileURL(join(process.cwd(), 'site', 'reference-panel.js')).href;
const m = await import(MODULE);

test('sitePath accepts a path on this site and refuses everything else', () => {
  assert.equal(m.sitePath('genomes/list.json', 'FALLBACK'), 'genomes/list.json');
  assert.equal(m.sitePath('../up/list.json', 'FALLBACK'), '../up/list.json');
  for (const bad of [null, '', ' list.json', 'list.json ', 'https://evil.example/x', '//evil.example/x', '/abs/path', '\\\\evil.example/x', 'javascript:alert(1)', 'a\u0000b', 'ftp://x/y']) {
    assert.equal(m.sitePath(bad, 'FALLBACK'), 'FALLBACK', JSON.stringify(bad));
  }
});

test('directAddress turns a Dropbox share link into the host that allows downloads, and leaves other addresses alone', () => {
  assert.equal(m.directAddress('https://www.dropbox.com/s/abc123/file.fa.gz?dl=0'), 'https://dl.dropboxusercontent.com/s/abc123/file.fa.gz');
  assert.equal(m.directAddress('https://www.dropbox.com/scl/fi/xyz/file.fa?rlkey=k&dl=1&raw=1'), 'https://dl.dropboxusercontent.com/scl/fi/xyz/file.fa?rlkey=k');
  for (const same of ['https://ftp.ncbi.nlm.nih.gov/genomes/x.fna.gz', 'https://www.dropbox.com/home', 'http://www.dropbox.com/s/a/b', 'not an address']) {
    assert.equal(m.directAddress(same), same);
  }
});

test('checkUrl accepts https and this site, and names the reason for everything else', () => {
  assert.equal(m.checkUrl('https://ftp.ncbi.nlm.nih.gov/x', 'the genome').host, 'ftp.ncbi.nlm.nih.gov');
  assert.equal(m.checkUrl('https://qmutfind.example.org/a', 'the genome').host, 'qmutfind.example.org');
  const keep = globalThis.location;
  try {
    globalThis.location = new URL('http://127.0.0.1:8000/page.html'); // a page served over plain http (the tests do this): its own address is fine
    assert.equal(m.checkUrl('http://127.0.0.1:8000/a', 'the genome').host, '127.0.0.1:8000');
    assert.throws(() => m.checkUrl('http://127.0.0.1:9000/a', 'the genome'), /only https addresses are accepted/, 'another port is another site');
  } finally {
    globalThis.location = keep;
  }
  assert.throws(() => m.checkUrl('http://other.example/a', 'the genome'), /the genome: only https addresses are accepted, not http:\/\//);
  assert.throws(() => m.checkUrl('https://user:pw@host.example/a', 'the genome'), /user name or password/);
  assert.throws(() => m.checkUrl('nonsense', 'the genome'), /"nonsense" is not an address/);
});

function fakeFetch(bytes, { status = 200, length = bytes.length } = {}) {
  return async () => new Response(new ReadableStream({
    start(c) {
      for (let i = 0; i < bytes.length; i += 4) c.enqueue(bytes.subarray(i, i + 4));
      c.close();
    },
  }), { status, headers: length ? { 'Content-Length': String(length) } : {} });
}
const bytesOf = async (file) => new Uint8Array(await file.arrayBuffer());

test('download gives a File with the bytes and the name, reports progress and takes the bytes off the budget', async () => {
  globalThis.fetch = fakeFetch(new TextEncoder().encode('>chr\nACGTACGT\n'));
  const budget = { left: 1000 };
  const seen = [];
  const file = await m.download({ name: 'ref.fa', url: 'https://host.example/ref.fa', bytes: 14 }, { label: 'the genome', budget, onProgress: (t) => seen.push(t) });
  assert.equal(file.name, 'ref.fa');
  assert.equal(new TextDecoder().decode(await bytesOf(file)), '>chr\nACGTACGT\n');
  assert.equal(budget.left, 1000 - 14);
  assert.ok(seen.length > 0 && seen.at(-1).startsWith('Downloading the genome: '), seen);
});

test('download names a gzip file that has no .gz in its address, and takes the name from the address when none is given', async () => {
  globalThis.fetch = fakeFetch(new Uint8Array([0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3]));
  const gz = await m.download({ url: 'https://host.example/dir/reads.fq' }, { label: 'reads' });
  assert.equal(gz.name, 'reads.fq.gz');
  globalThis.fetch = fakeFetch(new TextEncoder().encode('x'));
  const plain = await m.download({ url: 'https://host.example/dir/my%20file.txt' }, { label: 'a file' });
  assert.equal(plain.name, 'my file.txt');
});

test('download refuses a file above the budget from its announced size, and when it grows past the budget, and when the size is wrong', async () => {
  const data = new TextEncoder().encode('0123456789');
  globalThis.fetch = fakeFetch(data);
  await assert.rejects(m.download({ name: 'a', url: 'https://h.example/a' }, { label: 'the genome', budget: { left: 5 } }), /the genome: the download is too large .*the limit is/);
  globalThis.fetch = fakeFetch(data, { length: 0 }); // no Content-Length: the bytes are counted as they come
  await assert.rejects(m.download({ name: 'a', url: 'https://h.example/a' }, { label: 'the genome', budget: { left: 5 } }), /too large.*\(more than/);
  globalThis.fetch = fakeFetch(data);
  await assert.rejects(m.download({ name: 'a', url: 'https://h.example/a', bytes: 99 }, { label: 'the genome' }), /got 10 bytes, expected 99/);
  globalThis.fetch = fakeFetch(data, { status: 404 });
  await assert.rejects(m.download({ name: 'a', url: 'https://h.example/a' }, { label: 'the genome' }), /the server answered 404 for h\.example/);
});

test('a budget below the limit says that the other files of the sample take their share', async () => {
  globalThis.fetch = fakeFetch(new TextEncoder().encode('0123456789'));
  await assert.rejects(m.download({ name: 'a', url: 'https://h.example/a' }, { label: 'reads file 2 of 3', budget: { left: 5 } }), /together with the other files of this sample/);
});

test('the download limit is 2 GB, and ?download_limit_mb can only lower it', async () => {
  assert.equal(m.MAX_DOWNLOAD, 2000e6);
  const keep = globalThis.location;
  try {
    globalThis.location = new URL('https://qmutfind.example.org/?download_limit_mb=0.001');
    assert.equal((await import(`${MODULE}?low`)).MAX_DOWNLOAD, 1000);
    globalThis.location = new URL('https://qmutfind.example.org/?download_limit_mb=99999');
    assert.equal((await import(`${MODULE}?high`)).MAX_DOWNLOAD, 2000e6);
  } finally {
    globalThis.location = keep;
  }
});

test('mb writes megabytes with one decimal', () => {
  assert.equal(m.mb(1_500_000), '1.5 MB');
  assert.equal(m.mb(0), '0.0 MB');
});
