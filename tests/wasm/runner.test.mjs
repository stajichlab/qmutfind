// scripts/wasm-run.mjs streams gzip files through the WASM build. Its VCF must equal the CLI VCF.
import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, readFileSync, rmSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const CLI = process.env.QMUTFIND_BIN;

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
const FASTA = `>chrA test\n${chrA}\n`;

function reads(prefix, first, count) {
  const alt = chrA[50000] === 'G' ? 'T' : 'G';
  let fq = '';
  for (let i = first; i < first + count; i++) {
    const start = 49900 + i * 10;
    const seq = chrA.slice(start, start + 150).split('');
    seq[50000 - start] = alt;
    fq += `@${prefix}${i}\n${seq.join('')}\n+\n${'I'.repeat(150)}\n`;
  }
  return fq;
}

test('the runner streams gzip reference and two gzip read files and matches the CLI', () => {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-run-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), FASTA);
    writeFileSync(join(dir, 'ref.fa.gz'), gzipSync(FASTA));
    const r1 = reads('a', 0, 5), r2 = reads('b', 5, 5);
    writeFileSync(join(dir, 'r1.fq'), r1);
    writeFileSync(join(dir, 'r2.fq'), r2);
    writeFileSync(join(dir, 'r1.fq.gz'), gzipSync(r1));
    writeFileSync(join(dir, 'r2.fq.gz'), gzipSync(r2));

    const cli = spawnSync(CLI, ['call', '--reference', 'ref.fa', '--reads', 'r1.fq', '--reads2', 'r2.fq',
      '--sample', 's1', '--preset', 'short', '--out', 'cli.vcf'], { cwd: dir, encoding: 'utf8' });
    assert.equal(cli.status, 0, cli.stderr);

    const run = spawnSync('node', [join(process.cwd(), 'scripts', 'wasm-run.mjs'),
      '--reference', 'ref.fa.gz', '--reads', 'r1.fq.gz', '--reads', 'r2.fq.gz',
      '--sample', 's1', '--kind', 'short', '--out', 'wasm.vcf'], { cwd: dir, encoding: 'utf8' });
    assert.equal(run.status, 0, run.stderr);

    assert.equal(readFileSync(join(dir, 'wasm.vcf'), 'utf8'), readFileSync(join(dir, 'cli.vcf'), 'utf8'));
    const stats = JSON.parse(run.stdout.trim().split('\n').pop());
    assert.equal(stats.reads, 10);
    assert.ok(stats.seconds_total >= 0 && stats.rss_mb > 0);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test('a missing required option is reported and exits non-zero', () => {
  const run = spawnSync('node', [join(process.cwd(), 'scripts', 'wasm-run.mjs'), '--kind', 'short'], { encoding: 'utf8' });
  assert.notEqual(run.status, 0);
  assert.match(run.stderr, /--reference is required/);
});
