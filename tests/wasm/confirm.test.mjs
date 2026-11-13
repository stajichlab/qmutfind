// The WASM `Confirmation` must give the same report, JSON and tables as `qmutfind confirm`.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const require = createRequire(import.meta.url);
const wasm = require(join(process.cwd(), 'pkg-node', 'qmutfind_wasm.js'));
const CLI = process.env.QMUTFIND_BIN;

function lcg(seed) {
  let state = BigInt(seed);
  const mask = (1n << 64n) - 1n;
  return () => { state = (state * 6364136223846793005n + 1442695040888963407n) & mask; return state; };
}
const bases = (next, n) => Array.from({ length: n }, () => 'ACGT'[Number((next() >> 33n) % 4n)]).join('');
const revcomp = (s) => [...s].reverse().map((c) => ({ A: 'T', C: 'G', G: 'C', T: 'A' })[c]).join('');

const next = lcg(21);
const GENOME = bases(next, 120000);
const CASSETTE = bases(next, 2000);
const MODIFIED = GENOME.slice(0, 50000) + CASSETTE + GENOME.slice(51500);

function reads(seq, depth, seed) {
  const n = Math.floor((seq.length * depth) / 150);
  const r = lcg(seed);
  let out = '';
  for (let i = 0; i < n; i++) {
    const v = r() >> 20n;
    const start = Number(v % BigInt(seq.length - 150));
    let read = seq.slice(start, start + 150);
    if (((r() >> 40n) & 1n) === 1n) read = revcomp(read);
    out += `@r${i}\n${read}\n+\n${'I'.repeat(150)}\n`;
  }
  return out;
}
const FASTQ = reads(MODIFIED, 30, 7);
const enc = (s) => new TextEncoder().encode(s);

function cli(extra) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-conf-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), `>chr\n${GENOME}\n`);
    writeFileSync(join(dir, 'hph.fa'), `>hph\n${CASSETTE}\n`);
    writeFileSync(join(dir, 'reads.fq'), FASTQ);
    const r = spawnSync(CLI, ['confirm', '--reference', 'ref.fa', '--construct', 'hph.fa', '--reads', 'reads.fq', '--preset', 'short', '--out-prefix', 'o', '--always-exit-zero', ...extra], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    const read = (n) => readFileSync(join(dir, `o.${n}`), 'utf8');
    return { text: read('confirm.txt'), json: read('confirm.json'), junctions: read('junctions.tsv'), depth: read('cassette_depth.tsv') };
  } finally { rmSync(dir, { recursive: true, force: true }); }
}

function wasmRun(target, chunk) {
  const c = new wasm.Confirmation(enc(`>chr\n${GENOME}\n`), enc(`>hph\n${CASSETTE}\n`), 'short', target ? target.contig : undefined, target ? target.start : 0, target ? target.end : 0, 3, 10, 20);
  try {
    const bytes = enc(FASTQ);
    for (let i = 0; i < bytes.length; i += chunk) c.add_reads(bytes.subarray(i, i + chunk));
    c.end_file();
    return { text: c.report_text(), json: c.report_json(), junctions: c.junctions_tsv(), depth: c.depth_tsv(), reads: c.read_count() };
  } finally { c.free(); }
}

test('confirm: WASM text, JSON and tables equal the CLI, with a target', () => {
  const expected = cli(['--target', 'chr:50001-51500']);
  const got = wasmRun({ contig: 'chr', start: 50001, end: 51500 }, 1 << 18);
  assert.equal(got.text, expected.text);
  assert.equal(got.json, expected.json);
  assert.equal(got.junctions, expected.junctions);
  assert.equal(got.depth, expected.depth);
  assert.ok(got.text.includes('CONFIRMED') && got.text.includes('left flank 50000'), got.text);
});

test('confirm: the chunk size does not change the result, and no target is allowed', () => {
  const a = wasmRun(null, 4096);
  const b = wasmRun(null, 1 << 20);
  assert.equal(a.text, b.text);
  assert.equal(a.text, cli([]).text);
});

test('confirm: errors for a construct with the name of a genome contig, an empty construct, a bad kind and a missing target contig', () => {
  const genome = enc(`>chr\n${GENOME}\n`);
  assert.throws(() => new wasm.Confirmation(genome, enc(`>chr\n${CASSETTE}\n`), 'short', undefined, 0, 0, 3, 10, 20), /same name/);
  assert.throws(() => new wasm.Confirmation(genome, enc(''), 'short', undefined, 0, 0, 3, 10, 20), /no sequence/);
  assert.throws(() => new wasm.Confirmation(genome, enc(`>hph\n${CASSETTE}\n`), 'medium', undefined, 0, 0, 3, 10, 20), /short or long/);
  assert.throws(() => new wasm.Confirmation(genome, enc(`>hph\n${CASSETTE}\n`), 'short', 'nope', 1, 100, 3, 10, 20), /nope/);
});

test('confirm: find_gene in WASM returns the region of a gene by name', () => {
  const gff = '##gff-version 3\nchr\ts\tgene\t50001\t51500\t.\t+\t.\tID=gene-G1;Name=G1;locus_tag=G1\n';
  assert.deepEqual(JSON.parse(wasm.find_gene(enc(gff), 'G1')), { contig: 'chr', start: 50001, end: 51500 });
  assert.throws(() => wasm.find_gene(enc(gff), 'NOPE'), /NOPE/);
});
