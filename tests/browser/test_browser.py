"""The example page run in a real browser (headless Firefox) against the native CLI.

Run with `pixi run test-browser`. It needs the web package (`pixi run build-wasm-web`) and the CLI
binary (QMUTFIND_BIN). Firefox uses a different JavaScript engine from Node (SpiderMonkey, not V8).
"""
import gzip
import os
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parents[2]
CLI = os.environ.get("QMUTFIND_BIN")


def lcg_bases(seed, n):
    """The deterministic 64-bit generator the Rust and Node tests use."""
    state, mask, out = seed, (1 << 64) - 1, []
    for _ in range(n):
        state = (state * 6364136223846793005 + 1442695040888963407) & mask
        out.append("ACGT"[(state >> 33) % 4])
    return "".join(out)


CHR_A = lcg_bases(42, 120_000)
FASTA = f">chrA test\n{CHR_A}\n"
GFF = (
    "##gff-version 3\n"
    "chrA\tsrc\tgene\t49000\t51000\t.\t+\t.\tID=gene1;locus_tag=G1\n"
    "chrA\tsrc\tmRNA\t49000\t51000\t.\t+\t.\tID=tx1;Parent=gene1\n"
    "chrA\tsrc\texon\t49000\t51000\t.\t+\t.\tParent=tx1\n"
    "chrA\tsrc\tCDS\t49002\t50999\t.\t+\t0\tParent=tx1\n"
)
REGION = "chrA:48,500-51,500"


def snp_reads(count=10):
    alt = "T" if CHR_A[50_000] == "G" else "G"
    fastq = ""
    for i in range(count):
        start = 49_900 + i * 10
        read = list(CHR_A[start : start + 150])
        read[50_000 - start] = alt
        fastq += f"@r{i}\n{''.join(read)}\n+\n{'I' * 150}\n"
    return fastq


def run_cli(*args, cwd):
    result = subprocess.run([CLI, *args], cwd=cwd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result


@pytest.fixture(scope="session")
def server():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    proc = subprocess.Popen(
        [sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1", "--directory", str(ROOT)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    time.sleep(1.0)
    yield f"http://127.0.0.1:{port}"
    proc.terminate()


@pytest.fixture(scope="session")
def driver():
    assert CLI, "set QMUTFIND_BIN to the native qmutfind binary (pixi run test-browser does this)"
    assert (ROOT / "pkg-web" / "qmutfind_wasm.js").exists(), "build the web package first: pixi run build-wasm-web"
    options = Options()
    options.add_argument("-headless")
    options.binary_location = shutil.which("firefox")
    d = webdriver.Firefox(options=options, service=Service(executable_path=shutil.which("geckodriver")))
    d.set_page_load_timeout(60)
    yield d
    d.quit()


def run_page(driver, server, files, tmp_path, kind="short", region=REGION, wait=180):
    """Upload the files, click Run, wait for the result. Returns the texts the page produced."""
    driver.get(f"{server}/examples/browser/index.html")
    driver.find_element(By.ID, "reference").send_keys(str(files["reference"]))
    driver.find_element(By.ID, "reads").send_keys("\n".join(str(p) for p in files["reads"]))
    if files.get("gff"):
        driver.find_element(By.ID, "gff").send_keys(str(files["gff"]))
    driver.find_element(By.ID, "sample").clear()
    driver.find_element(By.ID, "sample").send_keys("s1")
    driver.find_element(By.ID, "region").send_keys(region)
    driver.find_element(By.ID, "run").click()
    WebDriverWait(driver, wait).until(
        lambda d: d.find_element(By.ID, "status").text.startswith(("Done", "Error"))
    )
    status = driver.find_element(By.ID, "status").text
    out = {"status": status}
    if status.startswith("Done"):
        fetch = "const cb = arguments[arguments.length - 1]; fetch(arguments[0]).then(r => r.text()).then(cb);"
        for name in ("vcf", "stats"):
            href = driver.find_element(By.ID, f"download-{name}").get_attribute("href")
            out[name] = driver.execute_async_script(fetch, href)
        out["report"] = driver.find_element(By.ID, "report").get_attribute("srcdoc")
        out["summary"] = driver.find_element(By.ID, "summary").text
    return out


def expected_from_cli(tmp_path):
    """The VCF, statistics and report the native CLI gives for the same inputs."""
    run_cli("call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1", "--preset", "short",
            "--out", "called.vcf", cwd=tmp_path)
    run_cli("annotate", "--vcf", "called.vcf", "--reference", "ref.fa", "--gff", "genes.gff3", "--out", "annot.vcf", cwd=tmp_path)
    run_cli("stats", "--vcf", "annot.vcf", "--out-prefix", "st", cwd=tmp_path)
    run_cli("report", "--vcf", "annot.vcf", "--gff", "genes.gff3", "--region", REGION, "--out", "report.html", cwd=tmp_path)
    return {
        "vcf": (tmp_path / "annot.vcf").read_text(),
        "stats": (tmp_path / "st.stats.json").read_text(),
        "report": (tmp_path / "report.html").read_text(),
    }


@pytest.fixture
def inputs(tmp_path):
    (tmp_path / "ref.fa").write_text(FASTA)
    (tmp_path / "reads.fq").write_text(snp_reads())
    (tmp_path / "genes.gff3").write_text(GFF)
    (tmp_path / "ref.fa.gz").write_bytes(gzip.compress(FASTA.encode()))
    (tmp_path / "reads.fq.gz").write_bytes(gzip.compress(snp_reads().encode()))
    return tmp_path


def test_the_page_gives_the_same_vcf_statistics_and_report_as_the_cli(driver, server, inputs):
    expected = expected_from_cli(inputs)

    got = run_page(driver, server, {"reference": inputs / "ref.fa", "reads": [inputs / "reads.fq"], "gff": inputs / "genes.gff3"}, inputs)

    assert got["status"].startswith("Done"), got["status"]
    assert got["vcf"] == expected["vcf"]
    assert got["stats"] == expected["stats"]
    assert got["report"] == expected["report"]
    assert "<svg" in got["report"] and "Calls: 1." in got["summary"]


def test_gzip_files_are_decompressed_in_the_browser_and_give_the_same_vcf(driver, server, inputs):
    expected = expected_from_cli(inputs)

    got = run_page(driver, server, {"reference": inputs / "ref.fa.gz", "reads": [inputs / "reads.fq.gz"], "gff": inputs / "genes.gff3"}, inputs)

    assert got["status"].startswith("Done"), got["status"]
    assert got["vcf"] == expected["vcf"]


def test_a_truncated_reads_file_is_reported_as_an_error_not_a_silent_success(driver, server, inputs):
    (inputs / "bad.fq").write_text("@r1\nACGT\n+\nIIII\n@r2\nGGGG\n")

    got = run_page(driver, server, {"reference": inputs / "ref.fa", "reads": [inputs / "bad.fq"]}, inputs, region="")

    assert got["status"].startswith("Error"), got["status"]
    assert "truncated" in got["status"]


def test_a_gzip_file_made_of_several_members_is_read_completely(driver, server, inputs):
    # parallel-fastq-dump and `cat a.gz b.gz` make multi-member gzip files. The browser's own
    # DecompressionStream rejects them in Firefox ("Unexpected input after the end of stream").
    expected = expected_from_cli(inputs)
    text = snp_reads()
    cut = text.index("@r5")
    (inputs / "multi.fq.gz").write_bytes(gzip.compress(text[:cut].encode()) + gzip.compress(text[cut:].encode()))

    got = run_page(driver, server, {"reference": inputs / "ref.fa", "reads": [inputs / "multi.fq.gz"], "gff": inputs / "genes.gff3"}, inputs)

    assert got["status"].startswith("Done"), got["status"]
    assert got["vcf"] == expected["vcf"]
