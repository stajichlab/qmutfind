// The WebAssembly build must give the same VCF as the native CLI for the same input.
// Build first with `pixi run build-wasm-node`, then run `pixi run test-wasm`.
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

// The same deterministic reference as the Rust tests (a 64-bit LCG).
function lcgBases(seed, n) {
  let state = BigInt(seed);
  const mask = (1n << 64n) - 1n;
  let out = '';
  for (let i = 0; i < n; i++) {
    state = (state * 6364136223846793005n + 1442695040888963407n) & mask;
    out += 'ACGT'[Number((state >> 33n) % 4n)];
  }
  return out;
}
const chrA = lcgBases(42, 120000);
const chrB = lcgBases(42, 200000).slice(120000); // continues the same stream as the Rust helper
const FASTA = `>chrA test\n${chrA}\n>chrB\n${chrB}\n`;

function reads(altCount, total) {
  const alt = chrA[50000] === 'G' ? 'T' : 'G';
  let fq = '';
  for (let i = 0; i < total; i++) {
    const start = 49900 + i * 10;
    const seq = chrA.slice(start, start + 150).split('');
    if (i < altCount) seq[50000 - start] = alt;
    fq += `@r${i}\n${seq.join('')}\n+\n${'I'.repeat(150)}\n`;
  }
  return fq;
}

function runWasm(fastq, chunkSize, thresholds, minMapq, minIndelAltReads) {
  const enc = new TextEncoder();
  const session = new wasm.Session(enc.encode(FASTA), 'short', 's1');
  if (thresholds) session.set_thresholds(...thresholds);
  if (minMapq !== undefined) session.set_min_mapq(minMapq);
  if (minIndelAltReads !== undefined) session.set_min_indel_alt_reads(minIndelAltReads);
  const bytes = enc.encode(fastq);
  for (let i = 0; i < bytes.length; i += chunkSize) session.add_reads(bytes.subarray(i, i + chunkSize));
  session.end_file();
  const vcf = session.vcf();
  session.free();
  return vcf;
}

function runCli(fastq, extraArgs = []) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-wasm-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), FASTA);
    writeFileSync(join(dir, 'reads.fq'), fastq);
    const r = spawnSync(CLI, ['call', '--reference', 'ref.fa', '--reads', 'reads.fq', '--sample', 's1',
      '--preset', 'short', '--out', 'out.vcf', ...extraArgs], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    return readFileSync(join(dir, 'out.vcf'), 'utf8');
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

test('the QMUTFIND_BIN environment variable points at the native CLI', () => {
  assert.ok(CLI, 'set QMUTFIND_BIN to the qmutfind binary (pixi run test-wasm does this)');
});

test('the WASM VCF is byte-identical to the native CLI VCF', () => {
  const fastq = reads(10, 10);
  const fromWasm = runWasm(fastq, 1 << 20);
  const fromCli = runCli(fastq);
  assert.equal(fromWasm, fromCli);
  assert.match(fromWasm, /SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1\.0000/);
});

test('the chunk size does not change the result', () => {
  const fastq = reads(10, 10);
  assert.equal(runWasm(fastq, 13), runWasm(fastq, 1 << 20));
});

test('set_thresholds matches the CLI flags', () => {
  const fastq = reads(7, 10); // fraction 0.7, below the short-read default of 0.8
  assert.equal(runWasm(fastq, 4096), runCli(fastq));
  assert.equal(runWasm(fastq, 4096, [5, 0.5, 0.5]), runCli(fastq, ['--min-snp-fraction', '0.5', '--min-indel-fraction', '0.5']));
  assert.notEqual(runWasm(fastq, 4096), runWasm(fastq, 4096, [5, 0.5, 0.5]));
});

test('an unknown read kind is an error that names the allowed values', () => {
  assert.throws(() => new wasm.Session(new TextEncoder().encode(FASTA), 'medium', 's1'), /short or long/);
});

test('set_min_mapq matches the --min-mapq flag', () => {
  const fastq = reads(10, 10);
  assert.equal(runWasm(fastq, 4096, undefined, 61), runCli(fastq, ['--min-mapq', '61']));
  assert.notEqual(runWasm(fastq, 4096, undefined, 61), runWasm(fastq, 4096));
});

function deletionReads(total) {
  let fq = '';
  for (let i = 0; i < total; i++) {
    const start = 49900 + i * 10;
    const seq = chrA.slice(start, start + 152).split('');
    const at = 50001 - start;
    seq.splice(at, 2);
    fq += `@d${i}\n${seq.join('')}\n+\n${'I'.repeat(150)}\n`;
  }
  return fq;
}

test('set_min_indel_alt_reads matches the --min-indel-alt-reads flag', () => {
  const fastq = deletionReads(10);
  assert.equal(runWasm(fastq, 4096, undefined, undefined, 10), runCli(fastq, ['--min-indel-alt-reads', '10']));
  assert.equal(runWasm(fastq, 4096, undefined, undefined, 11), runCli(fastq, ['--min-indel-alt-reads', '11']));
  assert.notEqual(runWasm(fastq, 4096, undefined, undefined, 10), runWasm(fastq, 4096, undefined, undefined, 11));
});
