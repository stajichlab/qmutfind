"""The site run in a real browser (headless Firefox) against the native CLI.

Run with `pixi run test-browser`. It needs the built site (`bash scripts/build-site.sh`) and the CLI binary
(QMUTFIND_BIN). Firefox uses a different JavaScript engine from Node (SpiderMonkey, not V8).
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


SITE = "dist-site/index.html"
THEMES = ["lab", "grey", "mono", "slate", "quest"]
VARIANT_POS = 50_001
HALF = 2500  # the page default window: 2.5 kb on each side
LOCUS = f"chrA:{VARIANT_POS - HALF}-{VARIANT_POS + HALF}"


def open_site(driver, server, theme=None):
    driver.get(f"{server}/{SITE}" + (f"?theme={theme}" if theme else ""))


def choose(driver, element_id, *paths):
    driver.find_element(By.ID, element_id).send_keys("\n".join(str(p) for p in paths))


def press(driver, button_id, status_id, wait=180):
    """Click a button and wait until the status line starts with Done, Drawn or Error. Returns the status text."""
    driver.find_element(By.ID, button_id).click()
    WebDriverWait(driver, wait).until(
        lambda d: d.find_element(By.ID, status_id).text.startswith(("Done", "Drawn", "Error", "There is no", "Choose", "Type a"))
        or "variants." in d.find_element(By.ID, status_id).text
    )
    return driver.find_element(By.ID, status_id).text


def download_text(driver, link_id):
    href = driver.find_element(By.ID, link_id).get_attribute("href")
    return driver.execute_async_script("const cb = arguments[arguments.length - 1]; fetch(arguments[0]).then(r => r.text()).then(cb);", href)


def step1(driver, server, inputs, reference="ref.fa", reads="reads.fq", theme=None):
    open_site(driver, server, theme)
    choose(driver, "reference", inputs / reference)
    choose(driver, "reads", inputs / reads)
    driver.find_element(By.ID, "sample").clear()
    driver.find_element(By.ID, "sample").send_keys("s1")
    return press(driver, "run", "status")


def cli_results(tmp_path):
    """What the native CLI gives for the same inputs: the VCF, the annotated VCF, statistics, plot and locus page."""
    run_cli("call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1", "--preset", "short", "--out", "called.vcf", cwd=tmp_path)
    run_cli("annotate", "--vcf", "called.vcf", "--reference", "ref.fa", "--gff", "genes.gff3", "--out", "annot.vcf", cwd=tmp_path)
    run_cli("stats", "--vcf", "called.vcf", "--out-prefix", "st", cwd=tmp_path)
    run_cli("stats", "--vcf", "annot.vcf", "--out-prefix", "sta", cwd=tmp_path)
    run_cli("plot", "--vcf", "called.vcf", "--out", "plot.svg", cwd=tmp_path)
    run_cli("report", "--vcf", "annot.vcf", "--gff", "genes.gff3", "--region", LOCUS, "--out", "report.html", cwd=tmp_path)
    read = lambda name: (tmp_path / name).read_text()
    return {"vcf": read("called.vcf"), "annot": read("annot.vcf"), "stats": read("st.stats.json"), "stats_annot": read("sta.stats.json"),
            "plot": read("plot.svg"), "report": read("report.html")}


@pytest.fixture
def inputs(tmp_path):
    (tmp_path / "ref.fa").write_text(FASTA)
    (tmp_path / "reads.fq").write_text(snp_reads())
    (tmp_path / "genes.gff3").write_text(GFF)
    (tmp_path / "ref.fa.gz").write_bytes(gzip.compress(FASTA.encode()))
    (tmp_path / "reads.fq.gz").write_bytes(gzip.compress(snp_reads().encode()))
    return tmp_path


# Step 1: call

def test_step_1_gives_the_same_vcf_as_the_cli(driver, server, inputs):
    expected = cli_results(inputs)

    status = step1(driver, server, inputs)

    assert status.startswith("Done"), status
    assert download_text(driver, "download-vcf") == expected["vcf"]
    assert "Calls: 1." in driver.find_element(By.ID, "summary").text


def test_gzip_files_are_decompressed_in_the_browser_and_give_the_same_vcf(driver, server, inputs):
    expected = cli_results(inputs)

    status = step1(driver, server, inputs, reference="ref.fa.gz", reads="reads.fq.gz")

    assert status.startswith("Done"), status
    assert download_text(driver, "download-vcf") == expected["vcf"]


def test_a_gzip_file_made_of_several_members_is_read_completely(driver, server, inputs):
    # parallel-fastq-dump and `cat a.gz b.gz` make multi-member gzip files. The browser's own
    # DecompressionStream rejects them in Firefox ("Unexpected input after the end of stream").
    expected = cli_results(inputs)
    text = snp_reads()
    cut = text.index("@r5")
    (inputs / "multi.fq.gz").write_bytes(gzip.compress(text[:cut].encode()) + gzip.compress(text[cut:].encode()))

    status = step1(driver, server, inputs, reads="multi.fq.gz")

    assert status.startswith("Done"), status
    assert download_text(driver, "download-vcf") == expected["vcf"]


def test_a_truncated_reads_file_is_reported_as_an_error_not_a_silent_success(driver, server, inputs):
    (inputs / "bad.fq").write_text("@r1\nACGT\n+\nIIII\n@r2\nGGGG\n")

    status = step1(driver, server, inputs, reads="bad.fq")

    assert status.startswith("Error") and "truncated" in status, status


# Step 2: summarize, using the VCF that step 1 left in the page

def test_step_2_uses_the_vcf_from_step_1_and_gives_the_cli_statistics_and_plot(driver, server, inputs):
    expected = cli_results(inputs)
    assert step1(driver, server, inputs).startswith("Done")

    status = press(driver, "summarize", "summarize-status")

    assert status.startswith("Done"), status
    assert download_text(driver, "download-stats") == expected["stats"]
    assert download_text(driver, "download-plot") == expected["plot"]
    assert len(driver.find_elements(By.CSS_SELECTOR, "#genome-plot svg .call")) == 1
    assert "Variants" in driver.find_element(By.ID, "stats-table").text


def test_step_2_can_start_from_a_vcf_file_with_no_step_1(driver, server, inputs):
    expected = cli_results(inputs)
    open_site(driver, server)
    choose(driver, "vcf-file", inputs / "called.vcf")

    status = press(driver, "summarize", "summarize-status")

    assert status.startswith("Done"), status
    assert download_text(driver, "download-plot") == expected["plot"]


def test_step_2_asks_for_the_reference_when_the_vcf_has_no_contig_lines(driver, server, inputs):
    bare = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrA\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1\n"
    (inputs / "bare.vcf").write_text(bare)
    open_site(driver, server)
    choose(driver, "vcf-file", inputs / "bare.vcf")

    status = press(driver, "summarize", "summarize-status")
    assert status.startswith("Error") and "reference" in status, status

    choose(driver, "reference", inputs / "ref.fa")
    status = press(driver, "summarize", "summarize-status")
    assert status.startswith("Done"), status


def test_step_2_without_any_vcf_says_so(driver, server):
    open_site(driver, server)

    status = press(driver, "summarize", "summarize-status")

    assert "no VCF" in status


# Step 3: annotate, pick a variant, draw its locus

def annotate_after_step_1(driver, server, inputs):
    assert step1(driver, server, inputs).startswith("Done")
    choose(driver, "gff", inputs / "genes.gff3")
    return press(driver, "annotate", "annotate-status")


def test_step_3_lists_each_variant_with_its_effect_as_the_cli_annotates_it(driver, server, inputs):
    expected = cli_results(inputs)

    status = annotate_after_step_1(driver, server, inputs)

    assert "1 of 1 variants" in status, status
    assert download_text(driver, "download-annotated") == expected["annot"]
    effect = next(f for f in expected["annot"].splitlines()[-1].split("\t")[7].split(";") if f.startswith("EFFECT=")).split("=")[1]
    row = driver.find_element(By.CSS_SELECTOR, "#variant-table tbody tr").text
    assert "chrA:50,001" in row and effect in row and "G1" in row, row
    assert effect in [o.get_attribute("value") for o in driver.find_elements(By.CSS_SELECTOR, "#effect-filter option")]


def test_pressing_draw_on_a_variant_shows_the_cli_locus_page_for_the_window_around_it(driver, server, inputs):
    expected = cli_results(inputs)
    assert annotate_after_step_1(driver, server, inputs)

    driver.find_element(By.CSS_SELECTOR, "#variant-table button.draw-row").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "report-status").text.startswith(("Drawn", "Error")))

    assert driver.find_element(By.ID, "report-status").text == f"Drawn: {LOCUS}"
    assert driver.find_element(By.ID, "region").get_attribute("value") == LOCUS
    assert driver.find_element(By.ID, "report").get_attribute("srcdoc") == expected["report"]


def test_a_region_can_also_be_typed(driver, server, inputs):
    expected = cli_results(inputs)
    assert annotate_after_step_1(driver, server, inputs)
    driver.find_element(By.ID, "region").send_keys(LOCUS)

    status = press(driver, "draw", "report-status")

    assert status == f"Drawn: {LOCUS}"
    assert driver.find_element(By.ID, "report").get_attribute("srcdoc") == expected["report"]


def test_step_2_after_step_3_carries_the_annotation_into_the_statistics(driver, server, inputs):
    expected = cli_results(inputs)
    assert annotate_after_step_1(driver, server, inputs)

    assert press(driver, "summarize", "summarize-status").startswith("Done")

    assert download_text(driver, "download-stats") == expected["stats_annot"]
    assert "Effects" in driver.find_element(By.ID, "stats-table").text


def test_a_new_run_of_step_1_replaces_the_annotated_vcf_so_step_2_never_uses_stale_results(driver, server, inputs):
    assert annotate_after_step_1(driver, server, inputs)
    # Reads copied from the reference carry no variant, so the new VCF has no calls.
    clean = "".join(f"@c{i}\n{CHR_A[49_900 + i * 10 : 50_050 + i * 10]}\n+\n{'I' * 150}\n" for i in range(10))
    (inputs / "clean.fq").write_text(clean)
    choose(driver, "reads", inputs / "clean.fq")
    assert press(driver, "run", "status").startswith("Done")

    assert press(driver, "summarize", "summarize-status").startswith("Done")

    import json
    assert json.loads(download_text(driver, "download-stats"))["total"] == 0
    assert "step 1" in driver.find_element(By.ID, "vcf-source").text


def test_step_3_without_the_gff_file_says_what_is_missing(driver, server, inputs):
    assert step1(driver, server, inputs).startswith("Done")

    status = press(driver, "annotate", "annotate-status")

    assert status.startswith("Choose") and "GFF3" in status, status


def test_step_3_warns_when_no_variant_falls_in_a_gene_because_the_contig_names_differ(driver, server, inputs):
    (inputs / "other.gff3").write_text(GFF.replace("chrA", "NC_999.1"))
    assert step1(driver, server, inputs).startswith("Done")
    choose(driver, "gff", inputs / "other.gff3")

    status = press(driver, "annotate", "annotate-status")

    assert "contig names" in status, status


# The page itself

def test_each_theme_can_be_chosen_by_the_address_and_has_its_own_colours(driver, server):
    seen = {}
    for theme in THEMES:
        open_site(driver, server, theme)
        assert driver.find_element(By.TAG_NAME, "html").get_attribute("data-theme") == theme
        seen[theme] = (
            driver.execute_script("return getComputedStyle(document.body).backgroundColor"),
            driver.execute_script("return getComputedStyle(document.body).color"),
        )
    assert len(set(seen.values())) == len(THEMES), seen


def test_the_theme_buttons_switch_the_theme_and_the_choice_is_kept_after_a_reload(driver, server):
    open_site(driver, server, "slate")
    driver.execute_script("localStorage.clear()")
    open_site(driver, server)
    driver.find_element(By.CSS_SELECTOR, "[data-set-theme=mono]").click()
    assert driver.find_element(By.TAG_NAME, "html").get_attribute("data-theme") == "mono"

    open_site(driver, server)

    assert driver.find_element(By.TAG_NAME, "html").get_attribute("data-theme") == "mono"


def test_an_unknown_theme_in_the_address_falls_back_to_the_default(driver, server):
    open_site(driver, server, "not-a-theme")
    driver.execute_script("localStorage.clear()")
    open_site(driver, server, "not-a-theme")

    assert driver.find_element(By.TAG_NAME, "html").get_attribute("data-theme") == "slate"


def test_the_page_states_its_limits_and_links_to_the_source(driver, server):
    open_site(driver, server)
    text = driver.find_element(By.TAG_NAME, "body").text
    assert "Limits" in text and "simulated" in text
    links = [a.get_attribute("href") for a in driver.find_elements(By.TAG_NAME, "a")]
    assert "https://github.com/stajichlab/qmutfind" in links
