#!/usr/bin/env bash
# Build the Node package and the native CLI, then run the WASM parity tests.
# Run from the repository root (the pixi task does). No BASH_SOURCE is used.
set -euo pipefail
target="${CARGO_TARGET_DIR:-$PWD/target}"
bash scripts/build-wasm-pkg.sh nodejs pkg-node
cargo build -p qmutfind-cli
QMUTFIND_BIN="$target/debug/qmutfind" node --test tests/wasm/engine.test.mjs tests/wasm/runner.test.mjs tests/wasm/annotate.test.mjs tests/wasm/stats.test.mjs tests/wasm/report.test.mjs
