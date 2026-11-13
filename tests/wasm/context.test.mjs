// The WASM context function must give the same text and HTML as `qmutfind context`.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const require = createRequire(import.meta.url);
const wasm = require(join(process.cwd(), 'pkg-node', 'qmutfind_wasm.js'));
const CLI = process.env.QMUTFIND_BIN;

// The hand-checked gene model of qmutfind-core/tests/effects.rs (contig c1): exon 21-50 (CDS from 24), intron, exon 81-110.
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

function cli(args) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-ctx-'));
  try {
    writeFileSync(join(dir, 'ref.fa'), FASTA);
    writeFileSync(join(dir, 'genes.gff3'), GFF);
    const r = spawnSync(CLI, ['context', '--reference', 'ref.fa', '--gff', 'genes.gff3', ...args], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    return r.stdout;
  } finally { rmSync(dir, { recursive: true, force: true }); }
}
const enc = (s) => new TextEncoder().encode(s);

for (const [name, pos, ref, alt] of [
  ['an SNV in the first exon', 30, C1[29], C1[29] === 'G' ? 'T' : 'G'],
  ['an insertion in the CDS', 30, C1[29], C1[29] + 'T'],
  ['a deletion over the intron edge is refused or shown the same by both', 50, C1.slice(49, 52), C1[49]],
  ['an intron SNV', 60, C1[59], C1[59] === 'A' ? 'C' : 'A'],
  ['a variant outside the gene', 150, C1[149], 'A'],
]) {
  test(`context: text and HTML from WASM equal the CLI for ${name}`, () => {
    const variant = `c1:${pos}:${ref}:${alt}`;
    const text = cli(['--variant', variant, '--flank', '15']);
    assert.equal(wasm.context(enc(FASTA), enc(GFF), 'c1', pos, ref, alt, 15, false), text);
    const html = cli(['--variant', variant, '--flank', '15', '--html']);
    assert.equal(wasm.context(enc(FASTA), enc(GFF), 'c1', pos, ref, alt, 15, true), html);
  });
}

test('context: a wrong reference allele is an error in WASM, as in the CLI', () => {
  assert.throws(() => wasm.context(enc(FASTA), enc(GFF), 'c1', 30, C1[29] === 'A' ? 'C' : 'A', 'T', 15, false), /reference allele/);
});
