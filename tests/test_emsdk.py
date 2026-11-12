"""Emscripten toolchain checks. The genomics-wasm ports need Emscripten 3.1.64 exactly.

Install it first with `pixi run setup-emsdk`. These tests fail, and do not skip, when it
is missing, so a CI run cannot pass without the toolchain.
"""
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EMSDK_DIR = Path(os.environ.get("EMSDK_DIR", ROOT / ".emsdk"))
EXPECTED_VERSION = "3.1.64"


def _in_emsdk(command, cwd=None):
    """Run a shell command with the Emscripten environment activated."""
    env_script = EMSDK_DIR / "emsdk_env.sh"
    assert env_script.exists(), f"Emscripten is not installed: {env_script} is missing. Run: pixi run setup-emsdk"
    return subprocess.run(
        ["bash", "-c", f'source "{env_script}" >/dev/null 2>&1 && {command}'],
        capture_output=True, text=True, cwd=cwd,
    )


def test_emcc_is_exactly_the_pinned_version():
    result = _in_emsdk("emcc --version")

    assert result.returncode == 0, result.stderr
    assert f"emcc (Emscripten gcc/clang-like replacement + linker emulating GNU ld) {EXPECTED_VERSION}" in result.stdout, result.stdout


def test_a_c_program_compiles_to_webassembly_and_runs_under_node(tmp_path):
    (tmp_path / "hello.c").write_text('#include <stdio.h>\nint main(void){ puts("hello from wasm"); return 0; }\n')

    build = _in_emsdk("emcc hello.c -O2 -o hello.js", cwd=tmp_path)
    assert build.returncode == 0, build.stderr
    assert (tmp_path / "hello.wasm").exists()

    run = subprocess.run(["node", str(tmp_path / "hello.js")], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert run.stdout.strip() == "hello from wasm"
