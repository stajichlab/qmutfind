// The WASM stats function must give the same JSON as `qmutfind stats`.
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
##contig=<ID=c1,length=100000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000;EFFECT=missense
c1\t180\t.\tC\tT\t.\tPASS\tCLASS=SNP;AF=0.9000
c1\t5000\t.\tCAA\tC\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.9500
`;

test('stats JSON from WASM equals the CLI JSON and parses', () => {
  const dir = mkdtempSync(join(tmpdir(), 'qmutfind-stats-'));
  try {
    writeFileSync(join(dir, 'in.vcf'), VCF);
    const r = spawnSync(CLI, ['stats', '--vcf', 'in.vcf', '--out-prefix', 'o', '--window', '100'], { cwd: dir, encoding: 'utf8' });
    assert.equal(r.status, 0, r.stderr);
    const fromCli = readFileSync(join(dir, 'o.stats.json'), 'utf8');
    const fromWasm = wasm.stats(VCF, 100);
    assert.equal(fromWasm, fromCli);
    const parsed = JSON.parse(fromWasm);
    assert.equal(parsed.total, 3);
    assert.equal(parsed.indel_lengths['-2'], 1);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});

test('a malformed VCF record is an error', () => {
  assert.throws(() => wasm.stats('c1\t5\t.\tA\n', 100), /fields/);
});
