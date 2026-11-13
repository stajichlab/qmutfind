// The browser page cannot run in a real browser from Node. This checks what Node can check:
// the JavaScript parses as ES modules, the page refers to files that exist, and the worker
// imports the package path that `pixi run build-wasm-web` writes. It is not a browser test.
import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { copyFileSync, existsSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const DIR = join(process.cwd(), 'site');

for (const file of ['worker.js', 'app.js']) {
  test(`${file} parses as an ES module`, () => {
    const dir = mkdtempSync(join(tmpdir(), 'qmutfind-example-'));
    try {
      const copy = join(dir, file.replace(/\.js$/, '.mjs'));
      copyFileSync(join(DIR, file), copy);
      const r = spawnSync(process.execPath, ['--check', copy], { encoding: 'utf8' });
      assert.equal(r.status, 0, r.stderr);
    } finally {
      rmSync(dir, { recursive: true, force: true });
    }
  });
}

test('the page loads the app script and the app starts the worker, and both files exist', () => {
  const html = readFileSync(join(DIR, 'index.html'), 'utf8');
  assert.match(html, /<script type="module" src="app\.js\?v=__BUILD__">/);
  assert.ok(existsSync(join(DIR, 'app.js')));
  assert.match(readFileSync(join(DIR, 'app.js'), 'utf8'), /new Worker\('worker\.js\?v=__BUILD__', \{ type: 'module' \}\)/);
});

test('the worker imports the web package from the path the build writes and uses the exports of the three steps', () => {
  const worker = readFileSync(join(DIR, 'worker.js'), 'utf8');
  assert.match(worker, /from '\.\/pkg-web\/qmutfind_wasm\.js\?v=__BUILD__'/);
  for (const name of ['Session', 'annotate', 'stats', 'report', 'genome_plot']) {
    assert.match(worker, new RegExp(`\\b${name}\\b`), `worker does not use ${name}`);
  }
});
