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

test('the WASM variant plot honours the contig order like the CLI', () => {
  const fromCli = cli(['--vcf', 'in.vcf', '--order', 'length'], { 'in.vcf': VCF });
  assert.equal(wasm.genome_plot(VCF, undefined, 'length'), fromCli);
  assert.throws(() => wasm.genome_plot(VCF, undefined, 'sideways'), /input, length or name/);
});

const BED = 'c1\t0\t500\t12.50\nc1\t500\t1000\t13.00\nc2\t0\t2000\t11.00\n';

function cliDepthPlot(args, files) {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-depth-'));
  try {
    for (const [name, text] of Object.entries(files)) writeFileSync(join(dir, name), text);
    const r = spawnSync(CLI, ['depth-plot', ...args, '--out', 'd.svg'], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    return readFileSync(join(dir, 'd.svg'), 'utf8');
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

test('the WASM depth plot equals the CLI depth plot, with variant ticks and a contig order', () => {
  const fromCli = cliDepthPlot(['--bed', 'd.bed', '--vcf', 'in.vcf', '--order', 'length'], { 'd.bed': BED, 'in.vcf': VCF });
  const fromWasm = wasm.depth_plot(BED, VCF, undefined, 'length');
  assert.equal(fromWasm, fromCli);
  assert.equal(fromWasm.match(/class="vtick"/g).length, 2);
});

test('the Session depth BED equals the CLI --depth-out file', () => {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-depthbed-'));
  try {
    const seq = Array.from({ length: 2000 }, (_, i) => 'ACGT'[(i * 7 + Math.floor(i / 3)) % 4]).join('');
    const reads = Array.from({ length: 8 }, (_, i) => `@r${i}\n${seq.slice(100 + i * 20, 250 + i * 20)}\n+\n${'I'.repeat(150)}\n`).join('');
    writeFileSync(join(dir, 'ref.fa'), `>c1\n${seq}\n`);
    writeFileSync(join(dir, 'reads.fq'), reads);
    const r = spawnSync(CLI, ['call', '--reference', 'ref.fa', '--reads', 'reads.fq', '--sample', 's', '--preset', 'short',
      '--out', 'o.vcf', '--depth-out', 'd.bed', '--depth-window', '500'], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    const session = new wasm.Session(enc.encode(`>c1\n${seq}\n`), 'short', 's');
    session.add_reads(enc.encode(reads));
    session.end_file();
    assert.equal(session.depth_bed(500), readFileSync(join(dir, 'd.bed'), 'utf8'));
    session.free();
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
