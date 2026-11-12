#!/usr/bin/env bash
# Assemble the static site in dist-site/: the page, its script and worker, the vendored gzip reader
# and the web package. Run from the repository root.
# Build the web package first: pixi run build-wasm-web. No BASH_SOURCE is used.
set -euo pipefail
test -f pkg-web/qmutfind_wasm.js || { echo "pkg-web is missing: run pixi run build-wasm-web first" >&2; exit 1; }
rm -rf dist-site
mkdir -p dist-site
cp site/index.html site/style.css site/theme.js dist-site/
cp site/app.js site/worker.js dist-site/
cp -r site/vendor dist-site/vendor
mkdir -p dist-site/pkg-web
cp pkg-web/qmutfind_wasm.js pkg-web/qmutfind_wasm_bg.wasm dist-site/pkg-web/
grep -q "./pkg-web/qmutfind_wasm.js" dist-site/worker.js || { echo "worker.js does not import ./pkg-web/qmutfind_wasm.js" >&2; exit 1; }
if [ -n "${SITE_CNAME:-}" ]; then echo "$SITE_CNAME" > dist-site/CNAME; fi
ls -l dist-site
