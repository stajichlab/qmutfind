#!/usr/bin/env bash
# Build the WebAssembly package: scripts/build-wasm-pkg.sh <nodejs|web> <output dir>
# Run from the repository root (the pixi tasks do). No BASH_SOURCE is used.
# The wasm-bindgen CLI version must equal the library version pinned in qmutfind-wasm/Cargo.toml.
set -euo pipefail
target_kind="${1:?usage: build-wasm-pkg.sh <nodejs|web> <output dir>}"
out_dir="${2:?usage: build-wasm-pkg.sh <nodejs|web> <output dir>}"
target="${CARGO_TARGET_DIR:-$PWD/target}"

# The SIMD128 flag is set here because Cargo reads config files from the current directory only.
RUSTFLAGS="-C target-feature=+simd128" cargo build --release --target wasm32-unknown-unknown -p qmutfind-wasm
rm -rf "$out_dir"
wasm-bindgen --target "$target_kind" --out-dir "$out_dir" "$target/wasm32-unknown-unknown/release/qmutfind_wasm.wasm"
ls -l "$out_dir"
