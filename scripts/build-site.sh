#!/usr/bin/env bash
# Assemble the static site in dist-site/: the page, its script and worker, the vendored gzip reader
# and the web package. Run from the repository root.
# Build the web package first: pixi run build-wasm-web. No BASH_SOURCE is used.
set -euo pipefail
test -f pkg-web/qmutfind_wasm.js || { echo "pkg-web is missing: run pixi run build-wasm-web first" >&2; exit 1; }
rm -rf dist-site
mkdir -p dist-site
cp site/index.html site/confirm_allele.html site/confirm.html site/style.css site/theme.js dist-site/
cp site/app.js site/confirm.js site/worker.js site/gz.js site/reads.js site/persist.js site/reference-panel.js dist-site/
cp -r site/vendor dist-site/vendor
cp -r site/example dist-site/example
cp -r site/genomes dist-site/genomes
mkdir -p dist-site/pkg-web
cp pkg-web/qmutfind_wasm.js pkg-web/qmutfind_wasm_bg.wasm dist-site/pkg-web/
grep -q "./pkg-web/qmutfind_wasm.js" dist-site/worker.js || { echo "worker.js does not import ./pkg-web/qmutfind_wasm.js" >&2; exit 1; }
# Name every file the page loads with one stamp (a hash of their content). GitHub Pages lets a browser keep a file for
# 10 minutes, so an old script could otherwise meet a new page. The source files hold the token __BUILD__ for the stamp.
STAMP=$(cat site/index.html site/confirm_allele.html site/confirm.html site/style.css site/theme.js site/app.js site/confirm.js site/worker.js site/gz.js site/reads.js site/persist.js site/reference-panel.js site/vendor/fflate.js pkg-web/qmutfind_wasm.js pkg-web/qmutfind_wasm_bg.wasm | sha1sum | cut -c1-10)
sed -i "s/__BUILD__/$STAMP/g" dist-site/index.html dist-site/confirm_allele.html dist-site/confirm.js dist-site/app.js dist-site/worker.js dist-site/gz.js dist-site/reads.js dist-site/persist.js dist-site/reference-panel.js
if grep -rl "__BUILD__" dist-site --include=*.html --include=*.js >/dev/null; then echo "a __BUILD__ token is left in dist-site" >&2; exit 1; fi
echo "build stamp: $STAMP"
if [ -n "${SITE_CNAME:-}" ]; then echo "$SITE_CNAME" > dist-site/CNAME; fi
ls -l dist-site
