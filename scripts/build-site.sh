#!/usr/bin/env bash
# Assemble the static site in dist-site/: the landing page, the example page script and worker,
# the vendored gzip reader and the web package. Run from the repository root.
# Build the web package first: pixi run build-wasm-web. No BASH_SOURCE is used.
set -euo pipefail
test -f pkg-web/qmutfind_wasm.js || { echo "pkg-web is missing: run pixi run build-wasm-web first" >&2; exit 1; }
rm -rf dist-site
mkdir -p dist-site
cp site/index.html site/style.css site/theme.js dist-site/
cp examples/browser/app.js examples/browser/worker.js dist-site/
cp -r examples/browser/vendor dist-site/vendor
mkdir -p dist-site/pkg-web
cp pkg-web/qmutfind_wasm.js pkg-web/qmutfind_wasm_bg.wasm dist-site/pkg-web/
# The worker imports the package from ../../pkg-web in the repository. In the site it is next to the page.
sed -i "s#'../../pkg-web/qmutfind_wasm.js'#'./pkg-web/qmutfind_wasm.js'#" dist-site/worker.js
grep -q "./pkg-web/qmutfind_wasm.js" dist-site/worker.js || { echo "worker.js import path was not rewritten" >&2; exit 1; }
if [ -n "${SITE_CNAME:-}" ]; then echo "$SITE_CNAME" > dist-site/CNAME; fi
ls -l dist-site
