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
    const [table, summary] = wasm.combine_tables(['parent', 'evolved'], [PARENT, EVOLVED], [bed(30), ''], parent, 10);
    const [cliTable, cliSummary] = cli(parent);
    assert.equal(table, cliTable);
    assert.equal(summary, cliSummary);
  }
});

test('combine in WASM reports a bad parent as an error', () => {
  assert.throws(() => wasm.combine_tables(['a'], [PARENT], [''], 'zzz', 10), /parent/);
});

test('combine in WASM needs one VCF and one coverage text per name', () => {
  assert.throws(() => wasm.combine_tables(['a', 'b'], [PARENT], ['', ''], undefined, 10), /same number/);
});
