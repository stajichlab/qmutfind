// Run a WASI module under Node. Usage: node run_wasi.cjs <module.wasm>
const fs = require('fs');
const { WASI } = require('node:wasi');
const wasi = new WASI({ version: 'preview1', args: ['probe'], env: {} });
const bytes = fs.readFileSync(process.argv[2]);
WebAssembly.compile(bytes).then(async (mod) => {
  const imports = wasi.getImportObject();
  const inst = await WebAssembly.instantiate(mod, imports);
  const t0 = process.hrtime.bigint();
  wasi.start(inst);
  console.error('wasm_run_ms', Number(process.hrtime.bigint() - t0) / 1e6);
}).catch((e) => { console.error('FAILED', e.message); process.exit(1); });
