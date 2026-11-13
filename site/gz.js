// A streaming gzip reader that checks every member. fflate (vendor/fflate.js) reads files with several members, which the built-in
// DecompressionStream does not do in Firefox, but it never compares the CRC-32 and the length that end each member, so a damaged
// download that still decodes would give wrong reads without a word. This wrapper makes that check.
import { Gunzip } from './vendor/fflate.js?v=__BUILD__';

const TABLE = (() => {
  const t = new Int32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c;
  }
  return t;
})();

/** The CRC-32 of `bytes`, continuing from `previous` (a value this function returned). */
export function crc32(bytes, previous = 0) {
  let c = ~previous;
  for (let i = 0; i < bytes.length; i++) c = TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return ~c >>> 0;
}

const KEEP = 1 << 16; // the last bytes of the compressed input that are kept, to read a trailer after the fact

/**
 * Like fflate's Gunzip: `push(chunk, final)`, and `ondata(chunk)` gets the unpacked bytes. When a member ends, its CRC-32 and
 * length are compared with the 8 bytes that end it, and an Error is thrown on a difference. The last member is checked when
 * `final` is true. A file that is not read to its end (the page stops at a read limit) is not checked.
 */
export class VerifyingGunzip {
  constructor(ondata) {
    this.ondata = ondata;
    this.crc = 0;
    this.length = 0; // unpacked bytes of this member
    this.members = 0;
    this.total = 0; // compressed bytes pushed so far
    this.tail = new Uint8Array(0); // the last KEEP compressed bytes
    this.gunzip = new Gunzip((chunk) => {
      this.crc = crc32(chunk, this.crc);
      this.length += chunk.length;
      this.ondata(chunk);
    });
    // fflate calls this when a new member starts; `offset` is where its header begins, so the trailer of the member before is
    // the 8 bytes just before it.
    this.gunzip.onmember = (offset) => this.check(offset);
  }

  /** Compare the member that just ended with the 8 bytes before byte `end` of the compressed input. */
  check(end) {
    const from = end - 8 - (this.total - this.tail.length);
    if (from < 0 || end > this.total) throw new Error('the gzip file could not be checked: a member header is too long to read back');
    const t = this.tail.subarray(from, from + 8);
    if (t.length < 8) throw new Error('the gzip file is cut inside the end of a member');
    const crc = (t[0] | (t[1] << 8) | (t[2] << 16) | (t[3] << 24)) >>> 0;
    const size = (t[4] | (t[5] << 8) | (t[6] << 16) | (t[7] << 24)) >>> 0;
    this.members++;
    if (crc !== this.crc) {
      throw new Error(`the gzip checksum does not match in member ${this.members} (the file is damaged; CRC-32 is ${this.crc.toString(16).padStart(8, '0')}, the file says ${crc.toString(16).padStart(8, '0')})`);
    }
    if (size !== this.length >>> 0) throw new Error(`the gzip length does not match in member ${this.members} (the file is damaged; ${this.length >>> 0} bytes unpacked, the file says ${size})`);
    this.crc = 0;
    this.length = 0;
  }

  push(chunk, final = false) {
    // Keep the compressed bytes before they are given on, so that a trailer can be read when fflate reports the next member.
    const joined = new Uint8Array(this.tail.length + chunk.length);
    joined.set(this.tail, 0);
    joined.set(chunk, this.tail.length);
    this.total += chunk.length;
    this.tail = joined.length > KEEP ? joined.subarray(joined.length - KEEP) : joined;
    this.gunzip.push(chunk, final);
    if (final) this.check(this.total); // the last member ends the file
  }
}
