#!/usr/bin/env bash
# Build the web package and the native CLI, then run the browser test in the `browser` environment.
# Run from the repository root (the pixi task does). No BASH_SOURCE is used.
set -euo pipefail
target="${CARGO_TARGET_DIR:-$PWD/target}"
bash scripts/build-wasm-pkg.sh web pkg-web
bash scripts/build-site.sh > /dev/null
cargo build -p qmutfind-cli
QMUTFIND_BIN="$target/debug/qmutfind" pixi run -e browser python -m pytest tests/browser -q
