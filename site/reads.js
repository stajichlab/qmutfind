// Keep the first reads of a FASTQ stream. The stream may be plain or gzip. The kept reads are packed again as gzip, so the page holds the
// compressed reads only. A record is four lines. This runs in the worker (see worker.js), because unpacking and packing take seconds for a
// million reads. It takes a reader (`stream.getReader()`) so that tests can feed it by hand.
import { Gzip } from './vendor/fflate.js?v=__BUILD__';
import { VerifyingGunzip } from './gz.js?v=__BUILD__';

const mb = (bytes) => `${Math.round(bytes / 1e6)} MB`;

/**
 * Read from `reader` until `maxReads` reads are kept or the stream ends. Gives { bytes, reads, cut }: the packed reads, how many, and
 * whether the stream was cut at the limit (the reader is then cancelled). A gzip stream that is read to its end is checked (CRC-32 and
 * length of every member); one that is cut is not read to its end, so it is not checked.
 * `maxKept` limits the packed bytes, `maxDownload` the bytes taken from the reader. Both give an Error that says "too large".
 * `onProgress({ downloaded, reads })` is called about every half second.
 */
export async function keepFirstReads(reader, { maxReads, maxKept, maxDownload, label, onProgress = () => {}, progressEvery = 500 }) {
  const target = maxReads * 4; // lines
  const packed = [];
  let packedBytes = 0;
  const gzip = new Gzip({ level: 3 }, (chunk) => { packed.push(chunk); packedBytes += chunk.length; });
  const tooManyKept = () => new Error(`${label}: the kept reads are too large (more than ${mb(maxKept)}).`);
  let lines = 0;
  let cut = false;
  let downloaded = 0;
  let lastReport = 0;
  let gunzip = null;
  let head = new Uint8Array(0); // the first bytes, until two are there to tell gzip from plain text
  let decided = false;

  const feed = (bytes) => {
    if (cut) return;
    let seen = 0;
    for (let i = 0; i < bytes.length; i++) {
      if (bytes[i] === 10 && ++seen + lines === target) { bytes = bytes.subarray(0, i + 1); cut = true; break; }
    }
    lines += seen;
    gzip.push(bytes, false);
    if (packedBytes > maxKept) throw tooManyKept();
  };
  const take = (value) => {
    if (!decided) {
      const joined = new Uint8Array(head.length + value.length);
      joined.set(head, 0);
      joined.set(value, head.length);
      if (joined.length < 2) { head = joined; return; }
      decided = true;
      head = new Uint8Array(0);
      if (joined[0] === 0x1f && joined[1] === 0x8b) gunzip = new VerifyingGunzip((chunk) => feed(chunk));
      value = joined;
    }
    if (gunzip) gunzip.push(value, false); else feed(value);
  };

  while (!cut) {
    const { done, value } = await reader.read();
    if (done) {
      if (!decided && head.length) take(new Uint8Array(0)); // a file of one byte is plain text
      if (!decided && head.length) { decided = true; feed(head); }
      if (gunzip) gunzip.push(new Uint8Array(0), true); // the end of the file: the last member is checked, and a file that ends early is an error
      break;
    }
    downloaded += value.length;
    if (downloaded > maxDownload) {
      await reader.cancel();
      throw new Error(`${label}: the download is too large (more than ${mb(maxDownload)}).`);
    }
    try {
      take(value);
    } catch (error) {
      if (!cut) throw error; // data after the limit may be anything, so an error there is ignored
    }
    if (progressEvery === 0 || performance.now() - lastReport > progressEvery) {
      lastReport = performance.now();
      onProgress({ downloaded, reads: Math.floor(lines / 4) });
    }
  }
  if (cut) await reader.cancel();
  gzip.push(new Uint8Array(0), true);
  if (packedBytes > maxKept) throw tooManyKept(); // the last flush can cross the limit
  const bytes = new Uint8Array(packedBytes);
  let offset = 0;
  for (const part of packed) { bytes.set(part, offset); offset += part.length; }
  return { bytes, reads: Math.floor(lines / 4), cut };
}
