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

    assert "no VCF" in status or "no VCF or coverage" in status


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
    assert 'class="variant-label"' in expected["report"] and '<table class="variants">' in expected["report"]


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


# Coverage: step 1 keeps a depth profile, step 2 draws it, and a mosdepth file can be used alone

def test_step_1_keeps_the_coverage_and_step_2_draws_it_like_the_cli(driver, server, inputs):
    run_cli("call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1", "--preset", "short", "--out", "called.vcf",
            "--depth-out", "depth.bed", "--depth-window", "10000", cwd=inputs)
    run_cli("depth-plot", "--bed", "depth.bed", "--vcf", "called.vcf", "--order", "name", "--out", "depth.svg", cwd=inputs)
    assert step1(driver, server, inputs).startswith("Done")
    assert download_text(driver, "download-depth-bed") == (inputs / "depth.bed").read_text()

    status = press(driver, "summarize", "summarize-status")

    assert status.startswith("Done") and "coverage plot" in status, status
    assert download_text(driver, "download-depth-plot") == (inputs / "depth.svg").read_text()
    assert len(driver.find_elements(By.CSS_SELECTOR, "#depth-plot svg .depth")) > 0
    assert len(driver.find_elements(By.CSS_SELECTOR, "#depth-plot svg .vtick")) == 1


def test_a_mosdepth_file_alone_gives_a_coverage_plot_and_the_contig_order_can_be_chosen(driver, server, inputs):
    (inputs / "m.regions.bed.gz").write_bytes(gzip.compress(b"chr10\t0\t5000\t30\nchr2\t0\t9000\t28\nchr1\t0\t7000\t31\n"))
    open_site(driver, server)
    choose(driver, "depth-file", inputs / "m.regions.bed.gz")
    labels = lambda: [t.get_attribute("data-contig-label") for t in driver.find_elements(By.CSS_SELECTOR, "#depth-plot svg text[data-contig-label]")]

    assert press(driver, "summarize", "summarize-status").startswith("Done")
    assert labels() == ["chr1", "chr2", "chr10"]  # the default is by name

    driver.find_element(By.CSS_SELECTOR, "#order option[value=length]").click()
    assert press(driver, "summarize", "summarize-status").startswith("Done")
    assert labels() == ["chr2", "chr1", "chr10"]

    driver.find_element(By.CSS_SELECTOR, "#order option[value=input]").click()
    assert press(driver, "summarize", "summarize-status").startswith("Done")
    assert labels() == ["chr10", "chr2", "chr1"]


# Examples: files that the page downloads by itself, listed in a manifest

import json
import shutil as _shutil

EXAMPLE_DIR = ROOT / "dist-site" / "_test_example"


@pytest.fixture
def example_site(inputs):
    """A manifest and files on the test server, so the tests never use NCBI or ENA."""
    EXAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    for name in ("ref.fa.gz", "genes.gff3", "reads.fq.gz"):
        _shutil.copy(inputs / name, EXAMPLE_DIR / name)

    def write(expected_reads=10, reads_bytes=None, reads_name="reads.fq.gz"):
        size = lambda n: (EXAMPLE_DIR / n).stat().st_size
        manifest = {
            "description": "A test example.", "sample": "ex", "kind": "short",
            "reference": {"name": "ref.fa.gz", "url": "ref.fa.gz", "bytes": size("ref.fa.gz")},
            "gff": {"name": "genes.gff3", "url": "genes.gff3", "bytes": size("genes.gff3")},
            "examples": [{"id": "small", "label": "Small test example", "expected_reads": expected_reads,
                          "reads": [{"name": "reads.fq.gz", "url": reads_name, "bytes": reads_bytes or size("reads.fq.gz")}]}],
        }
        (EXAMPLE_DIR / "manifest.json").write_text(json.dumps(manifest))

    write()
    yield write
    _shutil.rmtree(EXAMPLE_DIR, ignore_errors=True)


def load_example(driver, server):
    # A new query string each time: the test server answers 304 for a file rewritten within the same second.
    driver.get(f"{server}/{SITE}?example=_test_example/manifest.json%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "example-section").is_displayed())
    driver.find_element(By.CSS_SELECTOR, "[data-example=small]").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "example-status").text.startswith(("Loaded", "Error")))
    return driver.find_element(By.ID, "example-status").text


def test_an_example_fills_in_the_files_and_step_1_runs_without_choosing_any(driver, server, inputs, example_site):
    expected = cli_results(inputs)

    status = load_example(driver, server)
    assert status.startswith("Loaded"), status
    assert driver.find_element(By.ID, "sample").get_attribute("value") == "ex"
    result = press(driver, "run", "status")

    assert result.startswith("Done") and "WARNING" not in result, result
    assert download_text(driver, "download-vcf").replace("SAMPLE=ex", "SAMPLE=s1") == expected["vcf"].replace("SAMPLE=s1", "SAMPLE=s1")


def test_the_gene_file_of_the_example_is_used_by_step_3(driver, server, inputs, example_site):
    assert load_example(driver, server).startswith("Loaded")
    assert press(driver, "run", "status").startswith("Done")

    status = press(driver, "annotate", "annotate-status")

    assert "1 of 1 variants" in status, status


def test_a_file_chosen_by_the_user_replaces_the_examples_file(driver, server, inputs, example_site):
    clean = "".join(f"@c{i}\n{CHR_A[49_900 + i * 10 : 50_050 + i * 10]}\n+\n{'I' * 150}\n" for i in range(10))
    (inputs / "clean.fq").write_text(clean)
    assert load_example(driver, server).startswith("Loaded")
    choose(driver, "reads", inputs / "clean.fq")

    result = press(driver, "run", "status")

    assert result.startswith("Done"), result
    assert "Calls: 0." in driver.find_element(By.ID, "summary").text


def test_a_read_count_that_differs_from_the_manifest_gives_a_warning(driver, server, inputs, example_site):
    example_site(expected_reads=999)
    assert load_example(driver, server).startswith("Loaded")

    result = press(driver, "run", "status")

    assert "WARNING" in result and "999" in result, result


def test_a_download_with_the_wrong_size_is_an_error_and_loads_nothing(driver, server, inputs, example_site):
    example_site(reads_bytes=12345)

    status = load_example(driver, server)

    assert status.startswith("Error") and "incomplete" in status, status
    assert "Choose" in press(driver, "run", "status")


def test_a_missing_example_file_is_an_error_that_names_the_status_code(driver, server, inputs, example_site):
    example_site(reads_name="missing.fq.gz")

    status = load_example(driver, server)

    assert status.startswith("Error") and "404" in status, status


def test_the_manifest_address_must_be_a_path_on_the_site(driver, server):
    driver.get(f"{server}/{SITE}?example=https://example.org/manifest.json")
    # The default manifest is used instead, so the real example buttons appear and the other address is never fetched.
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "example-section").is_displayed())
    assert len(driver.find_elements(By.CSS_SELECTOR, "[data-example]")) == 2


def test_the_real_manifest_lists_https_files_on_ncbi_and_ena_and_a_site_file_of_the_stated_size():
    manifest = json.loads((ROOT / "site" / "example" / "manifest.json").read_text())
    hosts = {"ftp.ncbi.nlm.nih.gov", "ftp.sra.ebi.ac.uk"}
    from urllib.parse import urlparse
    for item in [manifest["reference"], manifest["gff"]] + [r for e in manifest["examples"] for r in e["reads"] if r["url"].startswith("http")]:
        parsed = urlparse(item["url"])
        assert parsed.scheme == "https" and parsed.netloc in hosts, item["url"]
        assert item["bytes"] > 0
    local = [r for e in manifest["examples"] for r in e["reads"] if not r["url"].startswith("http")]
    assert local and all((ROOT / "site" / "example" / r["url"]).stat().st_size == r["bytes"] for r in local)
    assert [e["id"] for e in manifest["examples"]] == ["small", "full"]


# The genome menu: NCBI genomes load in the page; FungiDB genomes only show the command that downloads them

GENOME_DIR = ROOT / "dist-site" / "_test_genomes"


@pytest.fixture
def genome_site(inputs):
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    for name in ("ref.fa.gz", "genes.gff3"):
        _shutil.copy(inputs / name, GENOME_DIR / name)
    size = lambda n: (GENOME_DIR / n).stat().st_size

    def write(fasta_bytes=None):
        entry = lambda acc, label, genome_size, fa_url: {
            "accession": acc, "label": label, "organism": label.split(" (")[0], "level": "Chromosome", "genome_size": genome_size,
            "genes": 100, "fasta": {"name": "ref.fa.gz", "url": fa_url, "bytes": fasta_bytes or size("ref.fa.gz")},
            "gff": {"name": "genes.gff3", "url": "genes.gff3", "bytes": size("genes.gff3")}}
        (GENOME_DIR / "genomes.json").write_text(json.dumps({"genomes": [
            entry("GCF_1.1", "Aspergillus testii (strain A)", 120_000, "ref.fa.gz"),
            entry("GCF_2.1", "Candida testii", 80_000_000, "ref.fa.gz"),
        ]}))
        (GENOME_DIR / "fungidb.json").write_text(json.dumps({"release": "FungiDB-71", "organisms": [
            {"name": "Aspergillus fakeii X1", "species": "Aspergillus fakeii", "genes": 5000},
            {"name": "Zymo fakeii", "species": "Zymo fakeii", "genes": 4000},
        ]}))

    write()
    yield write
    _shutil.rmtree(GENOME_DIR, ignore_errors=True)


def open_menu(driver, server):
    driver.get(f"{server}/{SITE}?genomes=_test_genomes/genomes.json%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "genome-section").is_displayed())


def options(driver):
    return [o.text for o in driver.find_elements(By.CSS_SELECTOR, "#genome-menu option")]


def choose_genome(driver, value_prefix):
    driver.find_element(By.CSS_SELECTOR, f"#genome-menu option[value^='{value_prefix}']").click()
    driver.find_element(By.ID, "genome-load").click()


def test_the_menu_lists_ncbi_genomes_and_fungidb_organisms_in_two_groups(driver, server, genome_site):
    open_menu(driver, server)

    groups = [g.get_attribute("label") for g in driver.find_elements(By.CSS_SELECTOR, "#genome-menu optgroup")]
    assert len(groups) == 2 and groups[0].startswith("NCBI") and groups[1].startswith("FungiDB"), groups
    texts = options(driver)
    assert any("Aspergillus testii" in t and "0.1 Mb" in t for t in texts), texts
    assert any("Zymo fakeii" in t for t in texts)


def test_the_filter_narrows_the_menu_to_matching_names(driver, server, genome_site):
    open_menu(driver, server)

    driver.find_element(By.ID, "genome-filter").send_keys("aspergillus")

    texts = options(driver)
    assert len(texts) == 2 and all("Aspergillus" in t for t in texts), texts


def test_loading_an_ncbi_genome_fills_in_the_reference_and_genes_and_step_1_then_runs_with_chosen_reads(driver, server, inputs, genome_site):
    expected = cli_results(inputs)
    open_menu(driver, server)

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))
    assert driver.find_element(By.ID, "genome-status").text.startswith("Loaded"), driver.find_element(By.ID, "genome-status").text
    choose(driver, "reads", inputs / "reads.fq")
    driver.find_element(By.ID, "sample").clear()
    driver.find_element(By.ID, "sample").send_keys("s1")
    assert press(driver, "run", "status").startswith("Done")
    assert download_text(driver, "download-vcf") == expected["vcf"]

    assert press(driver, "annotate", "annotate-status").startswith("1 of 1")  # the genome's gene file was loaded too


def test_a_large_genome_gets_a_warning_that_it_was_not_tested(driver, server, genome_site):
    open_menu(driver, server)

    choose_genome(driver, "ncbi:GCF_2.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))

    assert "not been tested" in driver.find_element(By.ID, "genome-status").text


def test_a_fungidb_entry_shows_the_fetch_command_and_downloads_nothing(driver, server, inputs, genome_site):
    open_menu(driver, server)

    choose_genome(driver, "fungidb:")

    command = driver.find_element(By.ID, "genome-command").text
    assert command.startswith('qmutfind fetch --organism "Aspergillus fakeii X1"') and "--out-dir" in command, command
    assert driver.find_element(By.ID, "genome-fungidb").is_displayed()
    assert "Choose a reference" in press(driver, "run", "status")  # nothing was loaded


def test_a_wrong_download_size_is_an_error_and_loads_nothing(driver, server, inputs, genome_site):
    genome_site(fasta_bytes=777)
    open_menu(driver, server)

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))

    assert driver.find_element(By.ID, "genome-status").text.startswith("Error") and "incomplete" in driver.find_element(By.ID, "genome-status").text


def test_the_real_genome_lists_are_complete_and_point_at_https_ncbi_files():
    from urllib.parse import urlparse
    document = json.loads((ROOT / "site" / "genomes" / "genomes.json").read_text())
    genomes = document["genomes"]
    assert 45 <= len(genomes) <= 60
    assert len({g["accession"] for g in genomes}) == len(genomes)
    for g in genomes:
        assert g["accession"].startswith("GCF_") and g["genes"] > 0 and g["genome_size"] > 0, g["accession"]
        for key in ("fasta", "gff"):
            assert urlparse(g[key]["url"]).scheme == "https" and urlparse(g[key]["url"]).netloc == "ftp.ncbi.nlm.nih.gov"
            assert g[key]["bytes"] > 0 and "//_" not in g[key]["url"]
    fungidb = json.loads((ROOT / "site" / "genomes" / "fungidb.json").read_text())
    assert fungidb["release"].startswith("FungiDB-") and len(fungidb["organisms"]) > 100


def test_choosing_a_genome_after_loading_an_example_drops_the_examples_reads(driver, server, inputs, example_site, genome_site):
    driver.get(f"{server}/{SITE}?example=_test_example/manifest.json%3Fv%3D{time.time_ns()}&genomes=_test_genomes/genomes.json%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "example-section").is_displayed() and d.find_element(By.ID, "genome-section").is_displayed())
    driver.find_element(By.CSS_SELECTOR, "[data-example=small]").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "example-status").text.startswith("Loaded"))

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith("Loaded"))

    # The example's reads belong to another genome, so step 1 must ask for reads again.
    assert "Choose a reference FASTA and at least one reads file" in press(driver, "run", "status")


# Where the data came from: every input says so, and the web import (an SRA run or https addresses) checks what it gets

ENA_DIR = ROOT / "dist-site" / "_test_ena"


def source_text(driver, slot):
    return driver.find_element(By.ID, f"source-{slot}").text


def test_every_input_says_where_its_data_comes_from(driver, server, inputs, example_site):
    open_site(driver, server)
    assert all(source_text(driver, s).startswith("Not set") for s in ("reference", "gff", "reads"))

    choose(driver, "reference", inputs / "ref.fa")
    assert source_text(driver, "reference").startswith("Local file: ref.fa")

    assert load_example(driver, server).startswith("Loaded")
    assert source_text(driver, "reference").startswith("From the example")
    assert "ref.fa.gz" in source_text(driver, "reference") and "127.0.0.1" in source_text(driver, "reference")
    assert source_text(driver, "reads").startswith("From the example") and "reads.fq.gz" in source_text(driver, "reads")

    choose(driver, "reads", inputs / "reads.fq")
    assert source_text(driver, "reads").startswith("Local file: reads.fq")  # a chosen file replaces the example's


def test_a_genome_from_the_menu_names_its_accession_and_the_ncbi_host_in_the_source_line(driver, server, genome_site):
    open_menu(driver, server)

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith("Loaded"))

    assert source_text(driver, "reference").startswith("From the genome menu") and "GCF_1.1" in source_text(driver, "reference")
    assert source_text(driver, "gff").startswith("From the genome menu")


@pytest.fixture
def ena_site(inputs):
    ENA_DIR.mkdir(parents=True, exist_ok=True)
    _shutil.copy(inputs / "reads.fq.gz", ENA_DIR / "reads.fq.gz")
    size = (ENA_DIR / "reads.fq.gz").stat().st_size

    def write(**changes):
        record = {"run_accession": "SRR9999991", "fastq_ftp": f"{HOST_PLACEHOLDER}/_test_ena/reads.fq.gz", "fastq_bytes": str(size),
                  "read_count": "10", "library_layout": "SINGLE", "instrument_platform": "ILLUMINA", "instrument_model": "Illumina Test",
                  "scientific_name": "Testus fungus", "sample_title": "a test sample"}
        record.update(changes)
        (ENA_DIR / "filereport").write_text(json.dumps([record]))

    write()
    yield write
    _shutil.rmtree(ENA_DIR, ignore_errors=True)


HOST_PLACEHOLDER = "{HOST}"


def open_web(driver, server, ena_dir="_test_ena"):
    driver.get(f"{server}/{SITE}?enaapi={ena_dir}/%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "web-section").is_displayed())


def fix_host(server):
    """The mock's file address must point at the test server."""
    path = ENA_DIR / "filereport"
    path.write_text(path.read_text().replace(HOST_PLACEHOLDER, server))


def sra_lookup(driver, accession):
    box = driver.find_element(By.ID, "sra-accession")
    box.clear()
    box.send_keys(accession)
    driver.find_element(By.ID, "sra-lookup").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "sra-info").text.strip() != "")
    return driver.find_element(By.ID, "sra-info").text


def test_an_sra_run_is_looked_up_shown_and_downloaded_with_its_checks(driver, server, inputs, ena_site):
    expected = cli_results(inputs)
    fix_host(server)
    open_web(driver, server)

    info = sra_lookup(driver, "srr9999991")

    assert "Testus fungus" in info and "ILLUMINA" in info and "SINGLE" in info and "10 reads" in info, info
    driver.find_element(By.ID, "sra-load").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "web-status").text.startswith(("Loaded", "Error")))
    assert driver.find_element(By.ID, "web-status").text.startswith("Loaded"), driver.find_element(By.ID, "web-status").text
    assert "SRR9999991" in source_text(driver, "reads") and "ENA" in source_text(driver, "reads")
    assert driver.find_element(By.ID, "kind").get_attribute("value") == "short"

    choose(driver, "reference", inputs / "ref.fa")
    driver.find_element(By.ID, "sample").clear()
    driver.find_element(By.ID, "sample").send_keys("s1")
    result = press(driver, "run", "status")

    assert result.startswith("Done") and "WARNING" not in result, result
    assert download_text(driver, "download-vcf") == expected["vcf"]


def test_a_read_count_that_differs_from_the_enas_gives_a_warning_after_mapping(driver, server, inputs, ena_site):
    ena_site(read_count="999")
    fix_host(server)
    open_web(driver, server)
    sra_lookup(driver, "SRR9999991")
    driver.find_element(By.ID, "sra-load").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "web-status").text.startswith("Loaded"))
    choose(driver, "reference", inputs / "ref.fa")

    result = press(driver, "run", "status")

    assert "WARNING" in result and "999" in result, result


@pytest.mark.parametrize("accession,fragment", [("hello", "does not look like"), ("SRR1", "does not look like"), ("PRJNA417720", "run accession")])
def test_something_that_is_not_a_run_accession_is_refused_before_any_request(driver, server, ena_site, accession, fragment):
    open_web(driver, server)

    assert fragment in sra_lookup(driver, accession)


def test_an_unknown_run_a_run_without_fastq_an_unsupported_platform_and_a_huge_run_are_each_refused_with_a_reason(driver, server, ena_site):
    fix_host(server)
    open_web(driver, server)
    cases = [
        ({"fastq_ftp": "", "fastq_bytes": ""}, "no FASTQ"),
        ({"instrument_platform": "PACBIO_SMRT"}, "not supported"),
        ({"fastq_bytes": str(5 * 10**9)}, "too large"),
    ]
    for changes, fragment in cases:
        ena_site(**changes)
        fix_host(server)
        open_web(driver, server)
        info = sra_lookup(driver, "SRR9999991")
        assert fragment in info, (changes, info)
        assert not driver.find_element(By.ID, "sra-load").is_displayed(), changes

    (ENA_DIR / "filereport").write_text("[]")
    open_web(driver, server)
    assert "not found" in sra_lookup(driver, "SRR9999992")


def test_a_paired_run_is_flagged_because_the_page_maps_each_mate_as_a_separate_read(driver, server, ena_site):
    ena_site(library_layout="PAIRED")
    fix_host(server)
    open_web(driver, server)

    info = sra_lookup(driver, "SRR9999991")

    assert "PAIRED" in info and "separate" in info, info


def test_https_addresses_load_each_file_and_say_where_they_came_from(driver, server, inputs, example_site):
    expected = cli_results(inputs)
    _shutil.copy(inputs / "ref.fa.gz", EXAMPLE_DIR / "refdata")  # no .gz in the name, but gzip inside
    open_web(driver, server)
    driver.find_element(By.ID, "url-reference").send_keys(f"{server}/dist-site/_test_example/refdata")
    driver.find_element(By.ID, "url-gff").send_keys(f"{server}/dist-site/_test_example/genes.gff3")
    driver.find_element(By.ID, "url-reads").send_keys(f"{server}/dist-site/_test_example/reads.fq.gz")

    driver.find_element(By.ID, "url-load").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "web-status").text.startswith(("Loaded", "Error")))
    assert driver.find_element(By.ID, "web-status").text.startswith("Loaded"), driver.find_element(By.ID, "web-status").text
    assert source_text(driver, "reference").startswith("From the address") and "refdata" in source_text(driver, "reference")
    driver.find_element(By.ID, "sample").clear()
    driver.find_element(By.ID, "sample").send_keys("s1")

    assert press(driver, "run", "status").startswith("Done")
    assert download_text(driver, "download-vcf") == expected["vcf"]


@pytest.mark.parametrize("address,fragment", [
    ("http://example.org/ref.fa.gz", "https"),
    ("ftp://example.org/ref.fa.gz", "https"),
    ("https://user:secret@example.org/ref.fa.gz", "user name"),
    ("not a url", "address"),
])
def test_an_unsafe_or_malformed_address_is_refused_before_any_request(driver, server, address, fragment):
    open_web(driver, server)
    driver.find_element(By.ID, "url-reference").send_keys(address)

    driver.find_element(By.ID, "url-load").click()
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "web-status").text.startswith("Error"))

    assert fragment in driver.find_element(By.ID, "web-status").text


def test_a_file_that_unpacks_to_more_than_the_limit_is_refused_instead_of_filling_the_memory(driver, server, inputs):
    big = "".join(f">c{i}\n{'ACGT' * 25_000}\n" for i in range(20))  # about 2 MB of bases
    (inputs / "big.fa").write_text(big)
    driver.get(f"{server}/{SITE}?limit_mb=1")
    choose(driver, "reference", inputs / "big.fa")
    choose(driver, "reads", inputs / "reads.fq")

    result = press(driver, "run", "status")

    assert result.startswith("Error") and "too large" in result, result


# Gene names and products in the variant table

def test_the_variant_table_shows_the_gene_name_and_the_product_from_the_gff(driver, server, inputs):
    named = GFF.replace("ID=gene1;locus_tag=G1", "ID=gene1;locus_tag=G1;Name=hrmA;gene=hrmA").replace(
        "ID=tx1;Parent=gene1", "ID=tx1;Parent=gene1;product=protein%20hrmA%3B%20test")
    (inputs / "named.gff3").write_text(named)
    assert step1(driver, server, inputs).startswith("Done")
    choose(driver, "gff", inputs / "named.gff3")

    status = press(driver, "annotate", "annotate-status")

    assert "1 of 1 variants" in status, status
    headers = [h.text for h in driver.find_elements(By.CSS_SELECTOR, "#variant-table thead th")]
    assert "Name" in headers and "Product" in headers, headers
    cells = {h: c.text for h, c in zip(headers, driver.find_elements(By.CSS_SELECTOR, "#variant-table tbody tr td"))}
    assert cells["Gene"] == "G1" and cells["Name"] == "hrmA" and cells["Product"] == "protein hrmA; test", cells
    driver.find_element(By.CSS_SELECTOR, "#variant-table button.draw-row").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "report-status").text.startswith("Drawn"))
    assert "hrmA (G1)" in driver.find_element(By.ID, "report").get_attribute("srcdoc")


# Plots must be readable in every theme

CONTRAST_JS = """
const parse = (c) => { const m = c.match(/rgba?\\(([^)]+)\\)/); const v = m[1].split(',').map(Number); return {r: v[0], g: v[1], b: v[2], a: v.length > 3 ? v[3] : 1}; };
const lin = (x) => { x /= 255; return x <= 0.03928 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4); };
const lum = (c) => 0.2126 * lin(c.r) + 0.7152 * lin(c.g) + 0.0722 * lin(c.b);
const ratio = (a, b) => { const [hi, lo] = [lum(a), lum(b)].sort((p, q) => q - p); return (hi + 0.05) / (lo + 0.05); };
const backdrop = (el) => { for (let e = el; e; e = e.parentElement) { const c = parse(e.ownerDocument.defaultView.getComputedStyle(e).backgroundColor); if (c.a > 0.5) return c; } return {r: 255, g: 255, b: 255, a: 1}; };
const blend = (f, o, bg) => ({r: f.r * o + bg.r * (1 - o), g: f.g * o + bg.g * (1 - o), b: f.b * o + bg.b * (1 - o)});
const worst = (els, kind) => {
  let min = Infinity, what = '';
  for (const el of els) {
    const cs = el.ownerDocument.defaultView.getComputedStyle(el);
    const bg = backdrop(el);
    let fg = parse(cs.fill);
    if (kind === 'mark') fg = blend(fg, parseFloat(cs.fillOpacity || '1'), bg);
    const r = ratio(fg, bg);
    if (r < min) { min = r; what = (el.getAttribute('class') || el.tagName) + ' ' + (el.getAttribute('data-class') || el.getAttribute('data-effect') || ''); }
  }
  return [min, what, els.length];
};
const q = (root, sel) => [...root.querySelectorAll(sel)];
const frame = document.getElementById('report').contentDocument;
const out = {};
out.genome_text = worst(q(document, '#genome-plot svg text'), 'text');
out.depth_text = worst(q(document, '#depth-plot svg text'), 'text');
out.genome_marks = worst(q(document, '#genome-plot svg .call, #genome-plot svg .bin'), 'mark');
out.depth_marks = worst(q(document, '#depth-plot svg .depth'), 'mark');
out.locus_text = worst(q(frame, 'svg text'), 'text');
out.locus_marks = worst(q(frame, 'svg .variant, svg .cds'), 'mark');
const body = frame.body, cs = frame.defaultView.getComputedStyle(body);
out.locus_body = [ratio(parse(cs.color), parse(cs.backgroundColor)), 'body text', 1];
out.locus_vs_panel = [ratio(parse(cs.backgroundColor), backdrop(document.getElementById('report'))), 'iframe background against its panel (1 = same)', 1];
return out;
"""


def test_plot_text_and_marks_are_readable_against_their_background_in_every_theme(driver, server, inputs):
    refs = CHR_A
    snp_ref = refs[50_000]
    snp_alt = "A" if snp_ref != "A" else "C"
    vcf = ("##fileformat=VCFv4.2\n##contig=<ID=chrA,length=120000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
           f"chrA\t50001\t.\t{snp_ref}\t{snp_alt}\t.\tPASS\tCLASS=SNP;AF=1.0000;DP=10\n"
           f"chrA\t50101\t.\t{refs[50100:50102]}\t{refs[50100]}\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.6000;DP=10\n"
           f"chrA\t10000\t.\t{refs[9999]}\t{'A' if refs[9999] != 'A' else 'C'}\t.\tPASS\tCLASS=SNP;AF=0.9000;DP=10\n"
           f"chrA\t60000\t.\t{refs[59999]}\t<INS>\t.\tPASS\tCLASS=LARGE_INDEL;AF=0.8000;DP=10;SVLEN=120\n")
    (inputs / "marks.vcf").write_text(vcf)
    (inputs / "cov.bed").write_text("".join(f"chrA\t{s}\t{s + 10000}\t{8 + (s // 10000) % 5}\n" for s in range(0, 120_000, 10_000)))
    open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")
    choose(driver, "gff", inputs / "genes.gff3")
    choose(driver, "vcf-file", inputs / "marks.vcf")
    choose(driver, "depth-file", inputs / "cov.bed")
    assert press(driver, "summarize", "summarize-status").startswith("Done")
    assert press(driver, "annotate", "annotate-status").startswith("4 of 4")
    driver.find_element(By.CSS_SELECTOR, "#variant-table button.draw-row[data-pos='50001']").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "report-status").text.startswith("Drawn"))
    time.sleep(0.5)

    problems = []
    for theme in THEMES:
        driver.find_element(By.CSS_SELECTOR, f"[data-set-theme={theme}]").click()
        time.sleep(0.4)
        result = driver.execute_script(CONTRAST_JS)
        for name, (ratio, what, count) in result.items():
            assert count > 0 or name == "locus_vs_panel", (theme, name, "nothing measured")
            needed = 1.0 if name == "locus_vs_panel" else (3.0 if name.endswith("marks") else 4.5)
            if name == "locus_vs_panel":
                continue
            if ratio < needed:
                problems.append(f"{theme}: {name} {ratio:.2f} < {needed} ({what})")
    assert not problems, "\n".join(problems)
