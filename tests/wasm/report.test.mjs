// The WASM report function must give the same page as `qmutfind report`.
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
const GFF = `##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
`;
const VCF = `##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t28\t.\tC\tA\t.\tPASS\tEFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D
c1\t60\t.\tA\tC\t.\tPASS\tEFFECT=intron_derived;GENE=G1;TRANSCRIPT=tx1
`;
const enc = new TextEncoder();

test('the WASM report equals the CLI report', () => {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-report-'));
  try {
    writeFileSync(join(dir, 'in.vcf'), VCF);
    writeFileSync(join(dir, 'genes.gff3'), GFF);
    const r = spawnSync(CLI, ['report', '--vcf', 'in.vcf', '--gff', 'genes.gff3', '--region', 'c1:1-200',
      '--title', 'Locus', '--out', 'r.html'], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    const fromCli = readFileSync(join(dir, 'r.html'), 'utf8');
    const fromWasm = wasm.report(VCF, enc.encode(GFF), 'c1:1-200', 'Locus');
    assert.equal(fromWasm, fromCli);
    assert.match(fromWasm, /<svg/);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test('a bad region is an error that shows the expected form', () => {
  assert.throws(() => wasm.report(VCF, undefined, 'c1:9-3', 'x'), /contig:start-end/);
});
