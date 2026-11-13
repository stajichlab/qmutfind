// The WASM combine function must give the same tables as `qmutfind combine`.
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

const vcf = (records) => '##fileformat=VCFv4.2\n##contig=<ID=chrA,length=3000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' +
  records.map(([pos, ref, alt, dp, af]) => `chrA\t${pos}\t.\t${ref}\t${alt}\t50.0\tPASS\tSAMPLE=s;CLASS=SNP;DP=${dp};AO=1;AF=${af}\n`).join('');
const bed = (depth) => [0, 1, 2].map((i) => `chrA\t${i * 1000}\t${(i + 1) * 1000}\t${depth}\n`).join('');

const PARENT = vcf([[500, 'A', 'G', 30, '1.0000']]);
const EVOLVED = vcf([[500, 'A', 'G', 28, '1.0000'], [1500, 'C', 'T', 25, '0.9600']]);

function cli(parent) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-combine-'));
  try {
    writeFileSync(join(dir, 'parent.vcf'), PARENT);
    writeFileSync(join(dir, 'evolved.vcf'), EVOLVED);
    writeFileSync(join(dir, 'parent.bed'), bed(30));
    const args = ['combine', '--strain', 'parent=parent.vcf,parent.bed', '--strain', 'evolved=evolved.vcf', '--out-prefix', 'out'];
    if (parent) args.push('--parent', parent);
    const result = spawnSync(CLI, args, { cwd: dir, encoding: 'utf8' });
    assert.equal(result.status, 0, result.stderr);
    return [readFileSync(join(dir, 'out.combined.tsv'), 'utf8'), readFileSync(join(dir, 'out.summary.tsv'), 'utf8')];
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

test('combine in WASM equals the CLI tables, with and without a parent (a missing coverage file is an empty string)', () => {
  for (const parent of [undefined, 'parent']) {
    const [table, summary] = wasm.combine_tables(['parent', 'evolved'], [PARENT, EVOLVED], [bed(30), ''], ['', ''], parent, 10);
    const [cliTable, cliSummary] = cli(parent);
    assert.equal(table, cliTable);
    assert.equal(summary, cliSummary);
  }
});

test('combine in WASM reports a bad parent as an error', () => {
  assert.throws(() => wasm.combine_tables(['a'], [PARENT], [''], [''], 'zzz', 10), /parent/);
});

test('combine in WASM needs one VCF and one coverage text per name', () => {
  assert.throws(() => wasm.combine_tables(['a', 'b'], [PARENT], ['', ''], ['', ''], undefined, 10), /same number/);
});

test('combine in WASM uses the evidence text of a strain (an empty string means none)', () => {
  const parent = vcf([]);
  const evidence = vcf([[1500, 'C', 'T', 30, '0.2000']]);
  const child = vcf([[1500, 'C', 'T', 25, '0.9600']]);
  const b = bed(30);
  const [table] = wasm.combine_tables(['parent', 'child'], [parent, child], [b, b], [evidence, ''], 'parent', 5);
  assert.match(table.split('\n')[1], /weak:0\.20:30\t0\.96:25\tparent\tunknown$/);
});

test('the session evidence VCF equals the CLI --evidence-out file', () => {
  let state = 42n;
  const next = () => { state = (state * 6364136223846793005n + 1442695040888963407n) & 0xFFFFFFFFFFFFFFFFn; return 'ACGT'[Number((state >> 33n) % 4n)]; };
  const seq = Array.from({ length: 20000 }, next).join('');
  const alt = seq[10000] === 'G' ? 'T' : 'G';
  let fastq = '';
  for (let i = 0; i < 10; i++) {
    const start = 9900 + i * 10;
    const read = seq.slice(start, start + 150).split('');
    if (i < 3) read[10000 - start] = alt; // 3 of 10 reads: below the 0.8 SNP threshold, above the 0.1 evidence threshold
    fastq += `@r${i}\n${read.join('')}\n+\n${'I'.repeat(150)}\n`;
  }
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-evidence-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), `>chrA\n${seq}\n`);
    writeFileSync(join(dir, 'reads.fq'), fastq);
    const result = spawnSync(CLI, ['call', '--reference', 'ref.fa', '--reads', 'reads.fq', '--sample', 's1', '--preset', 'short', '--out', 'out.vcf', '--evidence-out', 'evidence.vcf'], { cwd: dir, encoding: 'utf8' });
    assert.equal(result.status, 0, result.stderr);
    const session = new wasm.Session(new TextEncoder().encode(`>chrA\n${seq}\n`), 'short', 's1');
    session.add_reads(new TextEncoder().encode(fastq));
    session.end_file();
    const evidence = session.evidence_vcf();
    session.free();
    assert.equal(evidence, readFileSync(join(dir, 'evidence.vcf'), 'utf8'));
    assert.match(evidence, /\t10001\t/);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test('combine_tables returns the warnings as a third text, and version() gives the package version', () => {
  const bad = '##fileformat=VCFv4.2\n##contig=<ID=chrA,length=3000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n';
  const [, , warnings] = wasm.combine_tables(['a'], [bad], ['chrZ\t0\t1000\t30\n'], [''], undefined, 5);
  assert.match(warnings, /chrZ/);
  const [, , none] = wasm.combine_tables(['a'], [bad], [''], [''], undefined, 5);
  assert.equal(none, '');
  assert.match(wasm.version(), /^\d+\.\d+\.\d+/);
});
