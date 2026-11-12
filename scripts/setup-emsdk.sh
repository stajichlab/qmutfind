#!/usr/bin/env bash
# Install an exact Emscripten version from the official emsdk repository.
# Run from the repository root (the pixi task does). Safe to run again: it skips a matching install.
#   EMSDK_VERSION  default 3.1.64 (the genomics-wasm ports need this version)
#   EMSDK_DIR      default $PWD/.emsdk
# No BASH_SOURCE is used, so this works in SLURM jobs.
set -euo pipefail

version="${EMSDK_VERSION:-3.1.64}"
dir="${EMSDK_DIR:-$PWD/.emsdk}"

matching_install() {
    [ -f "$dir/emsdk_env.sh" ] || return 1
    bash -c "source '$dir/emsdk_env.sh' >/dev/null 2>&1 && emcc --version" 2>/dev/null | head -1 | grep -q " ${version}\b"
}

if matching_install; then
    echo "Emscripten ${version} is already installed in ${dir}"
    exit 0
fi

if [ ! -d "$dir/.git" ]; then
    # The emsdk tag has the same name as the Emscripten version.
    git clone --depth 1 --branch "$version" https://github.com/emscripten-core/emsdk.git "$dir"
fi
"$dir/emsdk" install "$version"
"$dir/emsdk" activate "$version"

matching_install || { echo "Emscripten ${version} did not install correctly" >&2; exit 1; }
echo "Installed Emscripten ${version} in ${dir}"
