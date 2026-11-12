// The WASM genome_plot function must give the same SVG as `qmutfind plot`.
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
const VCF = `##fileformat=VCFv4.2
##contig=<ID=c1,length=1000>
##contig=<ID=c2,length=2000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000
c2\t500\t.\tCT\tC\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.6000
`;
const BARE = `##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t5\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000
`;
const enc = new TextEncoder();

function cli(args, files) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-plot-'));
  try {
    for (const [name, text] of Object.entries(files)) writeFileSync(join(dir, name), text);
    const r = spawnSync(CLI, ['plot', ...args, '--out', 'p.svg'], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    return readFileSync(join(dir, 'p.svg'), 'utf8');
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

test('the WASM plot equals the CLI plot (contig lines in the VCF)', () => {
  const fromCli = cli(['--vcf', 'in.vcf'], { 'in.vcf': VCF });
  const fromWasm = wasm.genome_plot(VCF, undefined);
  assert.equal(fromWasm, fromCli);
  assert.equal(fromWasm.match(/class="call"/g).length, 2);
});

test('the WASM plot equals the CLI plot (contig lengths from the reference)', () => {
  const fasta = '>c1 first\nACGTACGTAC\n>c2\nACGTACGT\n';
  const fromCli = cli(['--vcf', 'in.vcf', '--reference', 'ref.fa'], { 'in.vcf': BARE, 'ref.fa': fasta });
  const fromWasm = wasm.genome_plot(BARE, enc.encode(fasta));
  assert.equal(fromWasm, fromCli);
  assert.match(fromWasm, /data-contig-label="c2"/);
});

test('a VCF with no contig lengths and no reference is an error that names the fix', () => {
  assert.throws(() => wasm.genome_plot(BARE, undefined), /reference/);
});
