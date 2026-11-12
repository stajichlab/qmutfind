// The WASM annotate function must give the same VCF as `qmutfind annotate`.
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

// The hand-checked gene model of qmutfind-core/tests/effects.rs (contig c1).
const C1 = 'T'.repeat(20) + 'CAC' + 'ATGGCTGAACGTAAATTAGATCCCTGG' + 'GT' + 'A'.repeat(26) + 'AG' +
  'GGCAAACTGGACCATGAATGCTAA' + 'CCCCCC' + 'T'.repeat(90);
const FASTA = `>c1\n${C1}\n`;
const GFF = `##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
`;
const VCF = `##fileformat=VCFv4.2
##contig=<ID=c1,length=200>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t28\t.\tC\tA\t.\tPASS\tDP=20
c1\t30\t.\tG\tGT\t.\tPASS\tDP=20
c1\t60\t.\tA\tC\t.\tPASS\t.
c1\t150\t.\tT\tA\t.\tPASS\tDP=20
`;

function runCli(args) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-annot-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), FASTA);
    writeFileSync(join(dir, 'genes.gff3'), GFF);
    writeFileSync(join(dir, 'in.vcf'), VCF);
    const r = spawnSync(CLI, ['annotate', '--vcf', 'in.vcf', '--reference', 'ref.fa', ...args, '--out', 'out.vcf'],
      { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    return readFileSync(join(dir, 'out.vcf'), 'utf8');
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

const enc = new TextEncoder();

test('annotate with a GFF equals the CLI output', () => {
  const fromWasm = wasm.annotate(VCF, enc.encode(FASTA), enc.encode(GFF));
  assert.equal(fromWasm, runCli(['--gff', 'genes.gff3']));
  assert.match(fromWasm, /EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D;NTX=1/);
  assert.match(fromWasm, /EFFECT=frameshift/);
});

test('annotate without a GFF equals the CLI output and marks every record unannotated', () => {
  const fromWasm = wasm.annotate(VCF, enc.encode(FASTA), undefined);
  assert.equal(fromWasm, runCli([]));
  assert.equal((fromWasm.match(/EFFECT=unannotated/g) || []).length, 4);
});

test('a bad GFF is an error that names the GFF', () => {
  assert.throws(() => wasm.annotate(VCF, enc.encode(FASTA), enc.encode('c1\tsrc\tgene\tx\t5\t.\t+\t.\tID=g\n')), /GFF/);
});
