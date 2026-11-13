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


REFERENCE_MODE_OF = {"reference": "upload", "gff": "upload", "url-reference": "url", "url-gff": "url"}


def set_reference_mode(driver, mode):
    driver.find_element(By.CSS_SELECTOR, f"#reference-source option[value='{mode}']").click()


def choose(driver, element_id, *paths):
    if element_id in REFERENCE_MODE_OF:  # the reference panel shows one way to give the genome at a time
        set_reference_mode(driver, REFERENCE_MODE_OF[element_id])
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


# The sample set: rows with a name, a source (SRA, web addresses or local files) and the reads. These helpers fill and run it.

def _blank_last_row(driver):
    return driver.execute_script("""
      const rows = document.querySelectorAll('#strain-table tbody tr'); const r = rows[rows.length - 1];
      if (!r) return false;
      return r.querySelector('.sample-name').value === '' && !r.querySelector('.sample-sra').value &&
        !r.querySelector('.sample-urls').value && !r.querySelector('.sample-files-names').textContent;
    """)


def add_sample(driver, name, source="local", value=None, note=""):
    """Fill the blank last row of the sample table (or add a row with +). `value` is a path or list of paths (local), a list of
    addresses (url) or an accession (sra)."""
    if not _blank_last_row(driver):
        driver.find_element(By.ID, "strain-add").click()
    row = driver.find_elements(By.CSS_SELECTOR, "#strain-table tbody tr")[-1]
    if name:
        row.find_element(By.CSS_SELECTOR, ".sample-name").send_keys(name)
    row.find_element(By.CSS_SELECTOR, f".sample-source option[value='{source}']").click()
    if value is not None:
        values = value if isinstance(value, (list, tuple)) else [value]
        if source == "local":
            row.find_element(By.CSS_SELECTOR, ".sample-files").send_keys("\n".join(str(v) for v in values))
        elif source == "url":
            row.find_element(By.CSS_SELECTOR, ".sample-urls").send_keys("\n".join(str(v) for v in values))
        else:
            row.find_element(By.CSS_SELECTOR, ".sample-sra").send_keys(str(values[0]))
    if note:
        row.find_element(By.CSS_SELECTOR, ".sample-note").send_keys(note)
    return row


def strain_states(driver):
    """{name: state text} of the named rows, read in one script call (the page rebuilds the rows, so held elements go stale)."""
    return {n: t for n, t in driver.execute_script(
        "return Array.from(document.querySelectorAll('#strain-table tbody tr')).map(r => [r.dataset.name, r.querySelector('.strain-state').textContent]);") if n}


def strain_state(driver, name):
    return strain_states(driver)[name]


def strain_rows(driver):
    return {r.get_attribute("data-name"): r for r in driver.find_elements(By.CSS_SELECTOR, "#strain-table tbody tr") if r.get_attribute("data-name")}


def wait_all_settled(driver, wait=180):
    WebDriverWait(driver, wait).until(lambda d: all(
        state.split(":")[0] in ("done", "failed", "cancelled", "imported") for state in strain_states(d).values())
        and d.find_element(By.ID, "strain-run").is_enabled())


def run_samples(driver, wait=180):
    """Press Run and wait until every named sample has settled. Returns {name: 'Done: N reads' | 'Error: ...' | state}."""
    driver.find_element(By.ID, "strain-run").click()
    wait_all_settled(driver, wait)
    reads = {n: t for n, t in driver.execute_script(
        "return Array.from(document.querySelectorAll('#strain-table tbody tr')).map(r => [r.dataset.name, r.querySelector('.strain-reads').textContent]);") if n}
    out = {}
    for name, state in strain_states(driver).items():
        if state.startswith("done"):
            out[name] = f"Done: {reads[name]} reads" + (" " + state.split(":", 1)[1].strip() if ":" in state else "")
        elif state.startswith("failed"):
            out[name] = "Error: " + state[len("failed: "):]
        else:
            out[name] = state
    return out


def sample_status(driver):
    return driver.find_element(By.ID, "strain-status").text


def set_read_limit_thousands(driver, value):
    box = driver.find_element(By.ID, "read-limit")
    box.clear()
    box.send_keys(str(value))


def step1(driver, server, inputs, reference="ref.fa", reads="reads.fq", theme=None):
    """Open the page, give the reference file and one local sample 's1', run it. Returns 'Done: N reads' or 'Error: ...'."""
    open_site(driver, server, theme)
    choose(driver, "reference", inputs / reference)
    add_sample(driver, "s1", "local", inputs / reads)
    return run_samples(driver)["s1"]


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


def test_choosing_another_sample_replaces_the_annotated_vcf_so_the_summary_never_uses_stale_results(driver, server, inputs):
    assert annotate_after_step_1(driver, server, inputs)
    # Reads copied from the reference carry no variant, so the second sample has no calls.
    clean = "".join(f"@c{i}\n{CHR_A[49_900 + i * 10 : 50_050 + i * 10]}\n+\n{'I' * 150}\n" for i in range(10))
    (inputs / "clean.fq").write_text(clean)
    add_sample(driver, "clean", "local", inputs / "clean.fq")
    assert run_samples(driver)["clean"].startswith("Done")
    driver.find_element(By.CSS_SELECTOR, "#result-sample option[value='clean']").click()

    assert press(driver, "summarize", "summarize-status").startswith("Done")

    import json
    assert json.loads(download_text(driver, "download-stats"))["total"] == 0
    assert "clean" in driver.find_element(By.ID, "vcf-source").text


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
            "--depth-out", "depth.bed", "--depth-window", "1000", cwd=inputs)  # the sample set uses 1 kb windows
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
    assert list(strain_states(driver)) == ["ex"] and strain_state(driver, "ex") == "waiting"
    assert "reads.fq.gz" in driver.find_element(By.CSS_SELECTOR, "#strain-table .sample-files-names").text
    result = run_samples(driver)["ex"]

    assert result.startswith("Done") and "WARNING" not in result, result
    assert download_text(driver, "download-vcf").replace("SAMPLE=ex", "SAMPLE=s1") == expected["vcf"].replace("SAMPLE=s1", "SAMPLE=s1")


def test_the_gene_file_of_the_example_is_used_by_step_3(driver, server, inputs, example_site):
    assert load_example(driver, server).startswith("Loaded")
    assert run_samples(driver)["ex"].startswith("Done")

    status = press(driver, "annotate", "annotate-status")

    assert "1 of 1 variants" in status, status


def test_the_user_can_add_a_sample_of_their_own_next_to_the_example(driver, server, inputs, example_site):
    clean = "".join(f"@c{i}\n{CHR_A[49_900 + i * 10 : 50_050 + i * 10]}\n+\n{'I' * 150}\n" for i in range(10))
    (inputs / "clean.fq").write_text(clean)
    assert load_example(driver, server).startswith("Loaded")
    add_sample(driver, "mine", "local", inputs / "clean.fq")

    results = run_samples(driver)

    assert results["ex"].startswith("Done") and results["mine"].startswith("Done"), results
    assert strain_rows(driver)["ex"].find_element(By.CSS_SELECTOR, ".strain-calls").text == "1"
    assert strain_rows(driver)["mine"].find_element(By.CSS_SELECTOR, ".strain-calls").text == "0"


def test_a_read_count_that_differs_from_the_manifest_gives_a_warning(driver, server, inputs, example_site):
    example_site(expected_reads=999)
    assert load_example(driver, server).startswith("Loaded")

    result = run_samples(driver)["ex"]

    assert "WARNING" in result and "999" in result, result


def test_a_download_with_the_wrong_size_is_an_error_and_loads_nothing(driver, server, inputs, example_site):
    example_site(reads_bytes=12345)

    status = load_example(driver, server)

    assert status.startswith("Error") and "incomplete" in status, status
    driver.find_element(By.ID, "strain-run").click()
    assert "reference" in sample_status(driver).lower()  # nothing was loaded


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
            {**entry("GCA_3.1", "Penicillium nogffii", 130_000, "ref.fa.gz"), "genes": 0, "gff": None, "source": "GenBank"},
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
    WebDriverWait(driver, 20).until(lambda d: len(d.find_elements(By.CSS_SELECTOR, "#genome-menu option")) > 0)


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
    add_sample(driver, "s1", "local", inputs / "reads.fq")
    assert run_samples(driver)["s1"].startswith("Done")
    assert download_text(driver, "download-vcf") == expected["vcf"]

    assert press(driver, "annotate", "annotate-status").startswith("1 of 1")  # the genome's gene file was loaded too


def test_a_genome_without_gene_models_loads_the_reference_only_and_says_so(driver, server, genome_site):
    open_menu(driver, server)

    label = [o.text for o in driver.find_elements(By.CSS_SELECTOR, "#genome-menu option[value='ncbi:GCA_3.1']")][0]
    assert "0 genes" not in label and "undefined" not in label and "no gene models" in label
    choose_genome(driver, "ncbi:GCA_3.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))

    status = driver.find_element(By.ID, "genome-status").text
    assert status.startswith("Loaded") and "no gene models" in status, status
    assert source_text(driver, "reference").startswith("From the genome menu (GCA_3.1")
    assert source_text(driver, "gff").startswith("Not set")


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
    add_sample(driver, "x", "local", inputs / "reads.fq")
    driver.find_element(By.ID, "strain-run").click()
    assert "reference" in sample_status(driver).lower()  # nothing was loaded


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
    assert 130 <= len(genomes) <= 250
    assert len({g["accession"] for g in genomes}) == len(genomes)
    by_accession = {g["accession"]: g for g in genomes}
    assert {"GCF_000002655.1", "GCA_000150145.1", "GCA_051225625.1"} <= set(by_accession)  # Af293, A1163, CEA10
    for g in genomes:
        refseq = g["accession"].startswith("GCF_")
        assert (refseq and g["source"] == "RefSeq") or (g["accession"].startswith("GCA_") and g["source"] == "GenBank"), g["accession"]
        assert g["genome_size"] > 0 and (g["genes"] > 0 or not refseq), g["accession"]
        for key in ("fasta", "gff"):
            if g[key] is None:  # only a GenBank genome may lack gene models
                assert key == "gff" and not refseq, g["accession"]
                continue
            assert urlparse(g[key]["url"]).scheme == "https" and urlparse(g[key]["url"]).netloc == "ftp.ncbi.nlm.nih.gov"
            assert g[key]["bytes"] > 0 and "//_" not in g[key]["url"]
    assert by_accession["GCA_051225625.1"]["gff"] is None and by_accession["GCA_000150145.1"]["gff"]["bytes"] > 0
    fungidb = json.loads((ROOT / "site" / "genomes" / "fungidb.json").read_text())
    assert fungidb["release"].startswith("FungiDB-") and len(fungidb["organisms"]) > 100


def test_choosing_a_genome_after_loading_an_example_drops_the_example_sample(driver, server, inputs, example_site, genome_site):
    driver.get(f"{server}/{SITE}?example=_test_example/manifest.json%3Fv%3D{time.time_ns()}&genomes=_test_genomes/genomes.json%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "example-section").is_displayed() and d.find_element(By.ID, "genome-section").is_displayed())
    driver.find_element(By.CSS_SELECTOR, "[data-example=small]").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "example-status").text.startswith("Loaded"))

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith("Loaded"))

    # The example's reads belong to another genome, so the example sample is gone from the list.
    assert strain_states(driver) == {}


# Where the data came from: every input says so, and the web import (an SRA run or https addresses) checks what it gets

ENA_DIR = ROOT / "dist-site" / "_test_ena"


def source_text(driver, slot):
    return driver.find_element(By.ID, f"source-{slot}").text


def test_every_input_says_where_its_data_comes_from(driver, server, inputs, example_site):
    open_site(driver, server)
    assert all(source_text(driver, s).startswith("Not set") for s in ("reference", "gff"))

    choose(driver, "reference", inputs / "ref.fa")
    assert source_text(driver, "reference").startswith("Local file: ref.fa")

    assert load_example(driver, server).startswith("Loaded")
    assert source_text(driver, "reference").startswith("From the example")
    assert "ref.fa.gz" in source_text(driver, "reference") and "127.0.0.1" in source_text(driver, "reference")
    assert "reads.fq.gz" in driver.find_element(By.CSS_SELECTOR, "#strain-table .sample-files-names").text

    choose(driver, "reference", inputs / "ref.fa")
    assert source_text(driver, "reference").startswith("Local file: ref.fa")  # a chosen file replaces the example's


def test_a_genome_from_the_menu_names_its_accession_and_the_ncbi_host_in_the_source_line(driver, server, genome_site):
    open_menu(driver, server)

    choose_genome(driver, "ncbi:GCF_1.1")
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith("Loaded"))

    assert source_text(driver, "reference").startswith("From the genome menu") and "GCF_1.1" in source_text(driver, "reference")
    assert source_text(driver, "gff").startswith("From the genome menu")


_clock = [time.time()]


def _write_fresh(path, text):
    """Write a file and give it a new modification time. The test server answers 304 for a file that was rewritten within
    the same second, and the browser would then keep the old copy."""
    path.write_text(text)
    _clock[0] += 10
    os.utime(path, (_clock[0], _clock[0]))


@pytest.fixture
def ena_site(inputs):
    ENA_DIR.mkdir(parents=True, exist_ok=True)
    _shutil.copy(inputs / "reads.fq.gz", ENA_DIR / "reads.fq.gz")
    size = (ENA_DIR / "reads.fq.gz").stat().st_size

    def write(**changes):
        record = {"run_accession": "SRR9999991", "fastq_ftp": f"{HOST_PLACEHOLDER}/dist-site/_test_ena/reads.fq.gz", "fastq_bytes": str(size),
                  "read_count": "10", "library_layout": "SINGLE", "instrument_platform": "ILLUMINA", "instrument_model": "Illumina Test",
                  "scientific_name": "Testus fungus", "sample_title": "a test sample"}
        record.update(changes)
        _write_fresh(ENA_DIR / "filereport", json.dumps([record]))

    write()
    yield write
    _shutil.rmtree(ENA_DIR, ignore_errors=True)


HOST_PLACEHOLDER = "{HOST}"


def open_web(driver, server, ena_dir="_test_ena"):
    """Open the page with a mock ENA service (a path on the test site)."""
    driver.get(f"{server}/{SITE}?enaapi={ena_dir}/%3Fv%3D{time.time_ns()}")
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "strain-table").is_displayed())


def fix_host(server):
    """The mock's file address must point at the test server."""
    path = ENA_DIR / "filereport"
    _write_fresh(path, path.read_text().replace(HOST_PLACEHOLDER, server))


def sra_check(driver, accession, row=None):
    """Type an accession in the (blank) row and press Check. Returns the text under it."""
    row = row or add_sample(driver, "", "sra", accession)
    row.find_element(By.CSS_SELECTOR, ".sample-check").click()
    WebDriverWait(driver, 30).until(lambda d: row.find_element(By.CSS_SELECTOR, ".sample-info").text.strip() not in ("", "Looking the run up in ENA"))
    return row.find_element(By.CSS_SELECTOR, ".sample-info").text


def open_with_reference(driver, server, inputs, ena=False):
    if ena:
        fix_host(server)
        open_web(driver, server)
    else:
        open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")


# SRA samples: looked up in ENA, downloaded, checked

def test_an_sra_sample_is_checked_in_ena_and_run_with_its_checks(driver, server, inputs, ena_site):
    expected = cli_results(inputs)
    open_with_reference(driver, server, inputs, ena=True)
    row = add_sample(driver, "s1", "sra", "srr9999991")

    info = sra_check(driver, "", row)

    assert "Testus fungus" in info and "ILLUMINA" in info and "SINGLE" in info and "10 reads" in info, info
    result = run_samples(driver)["s1"]
    assert result.startswith("Done: 10 reads") and "WARNING" not in result, result
    assert download_text(driver, "download-vcf") == expected["vcf"]


def test_a_read_count_that_differs_from_the_enas_gives_a_warning_after_mapping(driver, server, inputs, ena_site):
    ena_site(read_count="999")
    open_with_reference(driver, server, inputs, ena=True)
    add_sample(driver, "s1", "sra", "SRR9999991")

    result = run_samples(driver)["s1"]

    assert "WARNING" in result and "999" in result, result


@pytest.mark.parametrize("accession,fragment", [("hello", "does not look like"), ("SRR1", "does not look like"), ("PRJNA417720", "run accession")])
def test_something_that_is_not_a_run_accession_is_refused_before_any_request(driver, server, ena_site, accession, fragment):
    open_web(driver, server)

    assert fragment in sra_check(driver, accession)


def test_an_unknown_run_a_run_without_fastq_an_unsupported_platform_and_a_huge_run_are_each_refused_with_a_reason(driver, server, ena_site):
    cases = [
        ({"fastq_ftp": "", "fastq_bytes": ""}, "no FASTQ"),
        ({"instrument_platform": "PACBIO_SMRT"}, "not supported"),
        ({"fastq_bytes": str(5 * 10**9)}, "too large"),
    ]
    for changes, fragment in cases:
        ena_site(**changes)
        fix_host(server)
        open_web(driver, server)
        set_read_limit_thousands(driver, 0)  # with a read limit a big run is allowed (tested below)
        info = sra_check(driver, "SRR9999991")
        assert fragment in info, (changes, info)

    _write_fresh(ENA_DIR / "filereport", "[]")
    open_web(driver, server)
    assert "not found" in sra_check(driver, "SRR9999992")


def test_a_paired_run_is_flagged_because_the_page_maps_each_mate_as_a_separate_read(driver, server, ena_site):
    ena_site(library_layout="PAIRED")
    fix_host(server)
    open_web(driver, server)

    info = sra_check(driver, "SRR9999991")

    assert "PAIRED" in info and "separate" in info, info


def test_an_sra_run_is_cut_to_the_limit_and_the_read_count_check_uses_the_cut_count(driver, server, inputs, ena_site):
    open_with_reference(driver, server, inputs, ena=True)
    set_read_limit_thousands(driver, "0.004")  # 4 reads
    add_sample(driver, "s1", "sra", "SRR9999991")

    result = run_samples(driver)["s1"]

    assert result.startswith("Done: 4 reads") and "WARNING" not in result, result


def test_a_run_above_the_size_limit_is_allowed_when_a_read_limit_is_set(driver, server, ena_site):
    ena_site(fastq_bytes=str(5 * 10**9))
    fix_host(server)
    open_web(driver, server)
    set_read_limit_thousands(driver, 1)

    info = sra_check(driver, "SRR9999991")

    assert "too large" not in info and "first 1,000 reads" in info, info


# The read limit: in thousands of reads per file; the page keeps the first reads of each file and downloads stop there

def test_the_read_limit_box_has_a_default_in_1000s_and_says_what_it_means(driver, server):
    open_site(driver, server)

    box = driver.find_element(By.ID, "read-limit")
    assert box.get_attribute("value") == "2000"
    assert "1000s" in driver.find_element(By.CSS_SELECTOR, "label[for=read-limit]").text
    assert driver.find_element(By.ID, "read-limit-echo").text == "= 2,000,000 reads from each file"
    set_read_limit_thousands(driver, 500)
    assert driver.find_element(By.ID, "read-limit-echo").text == "= 500,000 reads from each file"
    set_read_limit_thousands(driver, 0)
    assert "All reads" in driver.find_element(By.ID, "read-limit-echo").text


def run_one(driver, server, inputs, source, value, limit=None, name="u"):
    open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")
    if limit is not None:
        set_read_limit_thousands(driver, limit)
    add_sample(driver, name, source, value)
    return run_samples(driver)[name]


def test_a_local_file_is_cut_to_the_first_n_reads(driver, server, inputs):
    assert run_one(driver, server, inputs, "local", inputs / "reads.fq", "0.004").startswith("Done: 4 reads")


def test_a_gzip_file_from_a_web_address_is_cut_to_the_first_n_reads(driver, server, inputs, example_site):
    result = run_one(driver, server, inputs, "url", f"{server}/dist-site/_test_example/reads.fq.gz", "0.004")

    assert result.startswith("Done: 4 reads"), result


def test_a_limit_above_the_number_of_reads_keeps_them_all(driver, server, inputs, example_site):
    result = run_one(driver, server, inputs, "url", f"{server}/dist-site/_test_example/reads.fq.gz", 100)

    assert result.startswith("Done: 10 reads"), result


def test_a_plain_fastq_from_a_web_address_is_cut_the_same_way(driver, server, inputs, example_site):
    _shutil.copy(inputs / "reads.fq", EXAMPLE_DIR / "reads.fq")

    result = run_one(driver, server, inputs, "url", f"{server}/dist-site/_test_example/reads.fq", "0.003")

    assert result.startswith("Done: 3 reads"), result


def test_the_download_stops_at_the_limit_so_data_after_it_is_never_read(driver, server, inputs, example_site):
    # A valid gzip member of 10 reads followed by 3 MB of noise. Reading past the reads would fail.
    (EXAMPLE_DIR / "tail.fq.gz").write_bytes(gzip.compress(snp_reads().encode()) + os.urandom(3_000_000))

    result = run_one(driver, server, inputs, "url", f"{server}/dist-site/_test_example/tail.fq.gz", "0.004")

    assert result.startswith("Done: 4 reads"), result


def test_zero_or_an_empty_box_means_no_limit_and_a_bad_value_is_refused(driver, server, inputs):
    assert run_one(driver, server, inputs, "local", inputs / "reads.fq", 0).startswith("Done: 10 reads")
    assert run_one(driver, server, inputs, "local", inputs / "reads.fq", "").startswith("Done: 10 reads")
    open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")
    set_read_limit_thousands(driver, "-5")
    add_sample(driver, "u", "local", inputs / "reads.fq")
    driver.find_element(By.ID, "strain-run").click()
    assert "read limit" in sample_status(driver), sample_status(driver)
    assert strain_state(driver, "u") == "waiting"


# Web addresses: several files make one sample

def test_several_web_addresses_are_downloaded_and_make_one_sample(driver, server, inputs, example_site):
    expected = cli_results(inputs)
    _shutil.copy(inputs / "reads.fq.gz", EXAMPLE_DIR / "reads_b.fq.gz")
    base = f"{server}/dist-site/_test_example"

    result = run_one(driver, server, inputs, "url", [f"{base}/reads.fq.gz", f"{base}/reads_b.fq.gz"], name="s1")

    assert result.startswith("Done: 20 reads"), result  # both files, mapped as separate reads
    assert strain_rows(driver)["s1"].find_element(By.CSS_SELECTOR, ".strain-calls").text == "1"
    assert "SAMPLE=s1" in download_text(driver, "download-vcf") and expected["vcf"].count("\n") == download_text(driver, "download-vcf").count("\n")


@pytest.mark.parametrize("address,fragment", [
    ("http://example.org/reads.fq.gz", "https"),
    ("ftp://example.org/reads.fq.gz", "https"),
    ("https://user:secret@example.org/reads.fq.gz", "user name"),
    ("not a url", "address"),
])
def test_an_unsafe_or_malformed_reads_address_fails_the_sample_before_any_request(driver, server, inputs, address, fragment):
    result = run_one(driver, server, inputs, "url", address)

    assert result.startswith("Error") and fragment in result, result


def test_a_reads_address_without_a_file_has_no_reads_and_the_page_says_so(driver, server, inputs):
    open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")
    add_sample(driver, "u", "url", None)

    driver.find_element(By.ID, "strain-run").click()

    assert "no reads" in sample_status(driver) and "address" in sample_status(driver), sample_status(driver)


@pytest.mark.parametrize("given,expected", [
    ("https://www.dropbox.com/scl/fi/abc/r2.fastq.gz?rlkey=K&st=S&dl=0", "https://dl.dropboxusercontent.com/scl/fi/abc/r2.fastq.gz?rlkey=K&st=S"),
    ("https://dropbox.com/scl/fi/abc/r2.fastq.gz?rlkey=K&dl=1", "https://dl.dropboxusercontent.com/scl/fi/abc/r2.fastq.gz?rlkey=K"),
    ("https://www.dropbox.com/s/abc/r2.fastq.gz?dl=0", "https://dl.dropboxusercontent.com/s/abc/r2.fastq.gz"),
    ("https://dl.dropboxusercontent.com/scl/fi/abc/r2.fastq.gz?rlkey=K", "https://dl.dropboxusercontent.com/scl/fi/abc/r2.fastq.gz?rlkey=K"),
    ("https://ftp.ncbi.nlm.nih.gov/x/y.fna.gz", "https://ftp.ncbi.nlm.nih.gov/x/y.fna.gz"),
    ("https://www.dropbox.com.example.org/s/abc/r2.fq.gz?dl=0", "https://www.dropbox.com.example.org/s/abc/r2.fq.gz?dl=0"),  # not Dropbox
    ("https://www.dropbox.com/home/folder", "https://www.dropbox.com/home/folder"),  # not a shared file
    ("not a url", "not a url"),
])
def test_a_dropbox_share_link_is_changed_to_the_direct_address_that_allows_browser_downloads(driver, server, given, expected):
    driver.get(f"{server}/{SITE}")

    assert driver.execute_script("return window.qmutfindDirectAddress(arguments[0]);", given) == expected


def test_a_file_that_unpacks_to_more_than_the_limit_is_refused_instead_of_filling_the_memory(driver, server, inputs):
    big = "".join(f">c{i}\n{'ACGT' * 25_000}\n" for i in range(20))  # about 2 MB of bases
    (inputs / "big.fa").write_text(big)
    driver.get(f"{server}/{SITE}?limit_mb=1")
    choose(driver, "reference", inputs / "big.fa")
    add_sample(driver, "s1", "local", inputs / "reads.fq")

    result = run_samples(driver)["s1"]

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


def test_a_download_larger_than_the_limit_is_refused_from_its_announced_size_before_it_is_read(driver, server, inputs, example_site):
    driver.get(f"{server}/{SITE}?download_limit_mb=0.001")  # a limit of about 1 KB
    choose(driver, "url-reference", f"{server}/dist-site/_test_example/ref.fa.gz")

    driver.find_element(By.ID, "ref-url-load").click()
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Error", "Loaded")))

    text = driver.find_element(By.ID, "genome-status").text
    # "the limit is" is the message of the check on the announced size, which runs before any byte is read.
    assert text.startswith("Error") and "too large" in text and "the limit is" in text, text
    assert source_text(driver, "reference").startswith("Not set")


def test_the_sample_section_lists_the_limits_and_risks_of_loading_reads_from_the_web(driver, server):
    open_site(driver, server)

    text = driver.find_element(By.CSS_SELECTOR, "#sample-section details").get_attribute("textContent")

    for needle in ["https", "2 GB", "MD5", "access-controlled", "Nothing is sent to this site"]:
        assert needle in text, needle


# Review fixes (2026-11-12): a failed run frees its session, same-site paths are checked by origin, VCF text is never HTML

def test_failed_runs_do_not_leak_memory_in_the_worker(driver, server, inputs):
    # Reads that end in a cut record fail at the end of the file. Before the fix each failed run kept its whole session.
    open_site(driver, server)
    driver.set_script_timeout(180)
    reference = (inputs / "ref.fa").read_text()
    bad_reads = (inputs / "reads.fq").read_text() + "@cut\nACGTACGT\n"
    sizes = driver.execute_async_script("""
        const [reference, reads, done] = arguments;
        const w = new Worker('worker.js', { type: 'module' });
        const next = (wanted) => new Promise((resolve) => {
          const on = (e) => { if (e.data.type === wanted || e.data.type === 'error') { w.removeEventListener('message', on); resolve(e.data); } };
          w.addEventListener('message', on);
        });
        (async () => {
          const out = [];
          for (let i = 0; i < 6; i++) {
            const answer = next('done');
            w.postMessage({ step: 'run', reference: new File([reference], 'ref.fa'), reads: [new File([reads], 'bad.fq')], sample: 's', kind: 'short' });
            const r = await answer;
            if (r.type !== 'error') { done({ failed: 'the run did not fail: ' + r.type }); return; }
            const m = next('done');
            w.postMessage({ step: 'memory' });
            out.push((await m).bytes);
          }
          done({ sizes: out });
        })();
    """, reference, bad_reads)

    assert "sizes" in sizes, sizes
    sizes = sizes["sizes"]
    assert sizes[-1] - sizes[0] <= 1_000_000, sizes  # six failed runs must not add more than 1 MB to the first one's size


@pytest.mark.parametrize("asked,expected", [
    ("_test_example/manifest.json", "_test_example/manifest.json"),
    ("sub/dir/list.json?v=1", "sub/dir/list.json?v=1"),
    ("\\\\evil.example/m.json", "FALLBACK"),
    ("\\evil.example/m.json", "FALLBACK"),
    ("//evil.example/m.json", "FALLBACK"),
    ("https://evil.example/m.json", "FALLBACK"),
    ("/etc/passwd", "FALLBACK"),
    ("javascript:alert(1)", "FALLBACK"),
    ("data:text/plain,x", "FALLBACK"),
    ("\t//evil.example/m.json", "FALLBACK"),  # the URL parser removes a leading tab, leaving //evil.example
    (" //evil.example/m.json", "FALLBACK"),
    ("\n\\evil.example/m.json", "FALLBACK"),
    ("", "FALLBACK"),
])
def test_a_page_address_option_may_only_name_a_path_on_this_site(driver, server, asked, expected):
    driver.get(f"{server}/{SITE}")

    got = driver.execute_script("return window.qmutfindSitePath(arguments[0], 'FALLBACK');", asked)

    assert got == expected


def test_effect_names_in_a_vcf_are_shown_as_text_not_as_html(driver, server, inputs):
    hostile = "<img src=x onerror=window.__pwned=1>"
    (inputs / "hostile.vcf").write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chrA,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chrA\t500\t.\tA\tG\t50\tPASS\tCLASS=SNP;DP=30;AO=30;AF=1.0;EFFECT={hostile}\n")
    open_site(driver, server)
    choose(driver, "vcf-file", inputs / "hostile.vcf")

    status = press(driver, "summarize", "summarize-status")

    assert status.startswith("Done"), status
    table = driver.find_element(By.ID, "stats-table")
    assert hostile in table.text
    assert len(table.find_elements(By.CSS_SELECTOR, "img")) == 0
    assert driver.execute_script("return window.__pwned === undefined;")


# The sample set: rows with + and -, three sources, a queue, a catalog, a combined table, an export and an import

def plain_reads(count=10):
    """Reads from the reference with no change: a sample without the SNP of `snp_reads`."""
    fastq = ""
    for i in range(count):
        start = 49_900 + i * 10
        fastq += f"@p{i}\n{CHR_A[start : start + 150]}\n+\n{'I' * 150}\n"
    return fastq


@pytest.fixture
def strain_files(inputs):
    (inputs / "parent.fq").write_text(plain_reads())
    (inputs / "garbage.fq").write_text("this is not a fastq file\n")
    return inputs


def open_strains(driver, server, inputs):
    open_site(driver, server)
    choose(driver, "reference", inputs / "ref.fa")


def row_count(driver):
    return len(driver.find_elements(By.CSS_SELECTOR, "#strain-table tbody tr"))


def test_the_sample_table_starts_with_one_blank_row_and_plus_and_minus_add_and_remove_rows(driver, server):
    open_site(driver, server)
    assert row_count(driver) == 1

    driver.find_element(By.ID, "strain-add").click()
    driver.find_element(By.ID, "strain-add").click()
    assert row_count(driver) == 3

    driver.find_elements(By.CSS_SELECTOR, "button.strain-remove")[0].click()
    assert row_count(driver) == 2
    assert driver.find_element(By.ID, "strain-add").text.startswith("+")
    assert driver.find_element(By.CSS_SELECTOR, "button.strain-remove").text == "\u2212"


def test_each_source_shows_only_its_own_reads_fields(driver, server):
    open_site(driver, server)
    row = driver.find_elements(By.CSS_SELECTOR, "#strain-table tbody tr")[0]
    shown = lambda: [row.find_element(By.CSS_SELECTOR, sel).is_displayed() for sel in (".sample-sra", ".sample-urls", ".sample-files")]
    options = [o.get_attribute("value") for o in row.find_elements(By.CSS_SELECTOR, ".sample-source option")]
    assert options == ["sra", "url", "local"]

    assert shown() == [False, False, True]  # the default source is local files
    row.find_element(By.CSS_SELECTOR, ".sample-source option[value='sra']").click()
    assert shown() == [True, False, False]
    row.find_element(By.CSS_SELECTOR, ".sample-source option[value='url']").click()
    assert shown() == [False, True, False]
    assert row.find_element(By.CSS_SELECTOR, ".sample-files").get_attribute("multiple") is not None


def test_a_sample_is_filled_with_a_name_a_source_and_reads_and_listed_as_waiting(driver, server, strain_files):
    open_strains(driver, server, strain_files)

    add_sample(driver, "parent", "local", strain_files / "parent.fq", note="the ancestor")

    assert list(strain_rows(driver)) == ["parent"]
    row = strain_rows(driver)["parent"]
    assert "parent.fq" in row.find_element(By.CSS_SELECTOR, ".sample-files-names").text
    assert row.find_element(By.CSS_SELECTOR, ".sample-note").get_attribute("value") == "the ancestor"
    assert strain_state(driver, "parent") == "waiting"


def test_a_sample_can_group_several_files(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "two", "local", [strain_files / "reads.fq", strain_files / "parent.fq"])

    result = run_samples(driver)["two"]

    assert result.startswith("Done: 20 reads"), result  # each file is mapped as separate reads
    names = strain_rows(driver)["two"].find_element(By.CSS_SELECTOR, ".sample-files-names").text
    assert "reads.fq" in names and "parent.fq" in names


def test_run_needs_a_name_for_every_sample_reads_for_each_and_names_that_are_not_used_twice(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    run = lambda: (driver.find_element(By.ID, "strain-run").click(), sample_status(driver))[1]

    add_sample(driver, "", "local", strain_files / "reads.fq")  # a reads file but no name
    assert "name" in run().lower()
    driver.find_element(By.CSS_SELECTOR, ".sample-name").send_keys("one")
    add_sample(driver, "one", "local", strain_files / "reads.fq")  # the same name twice
    assert "more than once" in run()
    driver.find_elements(By.CSS_SELECTOR, "button.strain-remove")[1].click()
    add_sample(driver, "empty", "local", None)  # a name but no reads
    assert "no reads" in run() and "empty" in run()
    assert set(strain_states(driver).values()) == {"waiting"}


def test_a_blank_row_is_ignored_when_the_samples_run(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "only", "local", strain_files / "reads.fq")
    driver.find_element(By.ID, "strain-add").click()  # a blank row stays

    results = run_samples(driver)

    assert results["only"].startswith("Done: 10 reads")


def test_a_sample_name_is_shown_as_text_and_never_as_html(driver, server, strain_files):
    hostile = "<img src=x onerror=window.__pwned=1>"
    open_strains(driver, server, strain_files)

    add_sample(driver, hostile, "local", strain_files / "reads.fq", note="<b>bold</b>")
    run_samples(driver)

    assert len(driver.find_elements(By.CSS_SELECTOR, "#strain-table img, #strain-table b, #result-sample img")) == 0
    assert driver.find_element(By.CSS_SELECTOR, "#strain-table .sample-name").get_attribute("value") == hostile
    assert hostile in [o.text for o in driver.find_elements(By.CSS_SELECTOR, "#result-sample option")]
    assert driver.execute_script("return window.__pwned === undefined;")


def test_run_maps_each_sample_in_turn_and_a_failed_sample_does_not_stop_the_others(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "bad", "local", strain_files / "garbage.fq")
    add_sample(driver, "child", "local", strain_files / "reads.fq")
    add_sample(driver, "parent", "local", strain_files / "parent.fq")

    run_samples(driver)

    assert strain_state(driver, "bad").startswith("failed")
    assert strain_state(driver, "child") == "done"
    assert strain_state(driver, "parent") == "done"
    assert strain_rows(driver)["child"].find_element(By.CSS_SELECTOR, ".strain-calls").text == "1"
    assert strain_rows(driver)["parent"].find_element(By.CSS_SELECTOR, ".strain-calls").text == "0"
    assert strain_rows(driver)["child"].find_element(By.CSS_SELECTOR, ".strain-reads").text == "10"


def test_a_failed_sample_that_is_edited_is_ready_to_run_again(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "x", "local", strain_files / "garbage.fq")
    run_samples(driver)
    assert strain_state(driver, "x").startswith("failed")

    strain_rows(driver)["x"].find_element(By.CSS_SELECTOR, ".sample-files").send_keys(str(strain_files / "reads.fq"))
    assert strain_state(driver, "x") == "waiting"
    assert run_samples(driver)["x"].startswith("Done: 10 reads")


def test_run_without_a_reference_says_so(driver, server, strain_files):
    open_site(driver, server)
    add_sample(driver, "x", "local", strain_files / "reads.fq")

    driver.find_element(By.ID, "strain-run").click()

    assert "reference" in sample_status(driver).lower()
    assert strain_state(driver, "x") == "waiting"


def test_cancel_stops_the_queue_and_leaves_no_sample_running(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    for name in ("a", "b", "c", "d"):
        add_sample(driver, name, "local", strain_files / "reads.fq")

    driver.find_element(By.ID, "strain-run").click()
    driver.find_element(By.ID, "strain-cancel").click()
    wait_all_settled(driver)

    states = {n: t.split(":")[0] for n, t in strain_states(driver).items()}
    assert set(states.values()) <= {"done", "cancelled"}, states
    assert "cancelled" in states.values(), states


def run_pair(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "parent", "local", strain_files / "parent.fq")
    add_sample(driver, "child", "local", strain_files / "reads.fq")
    run_samples(driver)


def combine(driver, parent="parent", min_depth="0"):  # the test reads are few, so a window mean is below 1
    driver.find_element(By.CSS_SELECTOR, f"#strain-parent option[value='{parent}']").click()
    box = driver.find_element(By.ID, "combine-min-depth")
    box.clear()
    box.send_keys(min_depth)
    driver.find_element(By.ID, "strain-combine").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "combine-status").text.startswith(("Done", "Error")))
    return driver.find_element(By.ID, "combine-status").text


def test_the_results_panel_follows_the_chosen_sample(driver, server, strain_files):
    run_pair(driver, server, strain_files)
    options = [o.get_attribute("value") for o in driver.find_elements(By.CSS_SELECTOR, "#result-sample option")]
    assert options == ["parent", "child"]  # in the order of the list

    driver.find_element(By.CSS_SELECTOR, "#result-sample option[value='child']").click()
    assert "Calls: 1." in driver.find_element(By.ID, "summary").text
    assert "child" in driver.find_element(By.ID, "vcf-source").text
    assert "SAMPLE=child" in download_text(driver, "download-vcf")
    driver.find_element(By.CSS_SELECTOR, "#result-sample option[value='parent']").click()
    assert "Calls: 0." in driver.find_element(By.ID, "summary").text
    vcf = download_text(driver, "download-vcf")
    assert "SAMPLE=child" not in vcf and not [l for l in vcf.splitlines() if l and not l.startswith("#")]  # the parent has no call


def test_removing_a_sample_removes_it_from_the_choices_too(driver, server, strain_files):
    run_pair(driver, server, strain_files)

    strain_rows(driver)["child"].find_element(By.CSS_SELECTOR, "button.strain-remove").click()

    assert list(strain_rows(driver)) == ["parent"]
    assert [o.get_attribute("value") for o in driver.find_elements(By.CSS_SELECTOR, "#result-sample option")] == ["parent"]
    assert [o.get_attribute("value") for o in driver.find_elements(By.CSS_SELECTOR, "#strain-parent option")] == ["", "parent"]
    driver.find_element(By.ID, "strain-combine").click()
    assert "two" in driver.find_element(By.ID, "combine-status").text.lower()


def test_combine_shows_the_new_variant_of_the_child_against_the_parent(driver, server, strain_files):
    run_pair(driver, server, strain_files)

    status = combine(driver, "parent")

    assert status.startswith("Done"), status
    rows = driver.find_elements(By.CSS_SELECTOR, "#combined-table tbody tr")
    assert len(rows) == 1
    text = rows[0].text
    assert "50001" in text and "new" in text, text
    tsv = download_text(driver, "download-combined")
    assert tsv.splitlines()[0].startswith("CHROM\tPOS\tREF\tALT\tCLASS\tstrains\tparent\tchild\tparent.status\tchild.status")
    assert tsv.splitlines()[1].endswith("\tparent\tnew")


def test_combine_without_a_parent_gives_a_plain_table_and_the_filter_narrows_it(driver, server, strain_files):
    run_pair(driver, server, strain_files)

    assert combine(driver, "").startswith("Done")

    assert "status" not in download_text(driver, "download-combined").splitlines()[0]
    driver.find_element(By.CSS_SELECTOR, "#combined-filter option[value='shared']").click()
    assert len(driver.find_elements(By.CSS_SELECTOR, "#combined-table tbody tr")) == 0  # the SNP is in one sample only
    driver.find_element(By.CSS_SELECTOR, "#combined-filter option[value='one']").click()
    assert len(driver.find_elements(By.CSS_SELECTOR, "#combined-table tbody tr")) == 1


def test_combine_with_fewer_than_two_finished_samples_says_so(driver, server, strain_files):
    open_strains(driver, server, strain_files)
    add_sample(driver, "only", "local", strain_files / "reads.fq")

    driver.find_element(By.ID, "strain-combine").click()

    assert "two" in driver.find_element(By.ID, "combine-status").text.lower()


def fetch_bytes(driver, link_id):
    href = driver.find_element(By.ID, link_id).get_attribute("href")
    return bytes(driver.execute_async_script(
        "const cb = arguments[arguments.length - 1]; fetch(arguments[0]).then(r => r.arrayBuffer()).then(b => cb(Array.from(new Uint8Array(b))));", href))


def test_the_catalog_export_is_a_zip_with_the_catalog_the_files_and_the_combined_table(driver, server, strain_files):
    import io
    import zipfile

    run_pair(driver, server, strain_files)
    combine(driver, "parent")
    driver.find_element(By.ID, "strain-export").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "download-catalog").get_attribute("href"))

    archive = zipfile.ZipFile(io.BytesIO(fetch_bytes(driver, "download-catalog")))

    names = archive.namelist()
    assert "catalog.json" in names and "combined.tsv" in names and "summary.tsv" in names
    catalog = json.loads(archive.read("catalog.json"))
    assert catalog["format"] == "qmutfind-catalog" and catalog["version"] == 1
    assert catalog["parent"] == "parent" and [s["name"] for s in catalog["strains"]] == ["parent", "child"]
    child = catalog["strains"][1]
    assert child["status"] == "done" and child["calls"] == 1 and child["reads"] == 10 and child["files"] == ["reads.fq"]
    assert child["source"] == "local"
    assert archive.read(child["vcf"]).decode().count("\n") > 5
    assert archive.read(child["bed"]) and archive.read(child["evidence"])
    assert catalog["settings"]["kind"] == "short" and catalog["settings"]["min_depth"] == 0 and catalog["settings"]["read_limit_per_file"] == 2_000_000
    assert catalog["reference"]["name"] == "ref.fa" and catalog["qmutfind_version"]
    assert archive.read("combined.tsv").decode() == download_text(driver, "download-combined")


def test_an_exported_catalog_can_be_imported_again_and_combined_to_the_same_table(driver, server, strain_files, tmp_path):
    run_pair(driver, server, strain_files)
    first = combine(driver, "parent")
    expected = download_text(driver, "download-combined")
    driver.find_element(By.ID, "strain-export").click()
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "download-catalog").get_attribute("href"))
    (tmp_path / "catalog.zip").write_bytes(fetch_bytes(driver, "download-catalog"))

    open_site(driver, server)
    choose(driver, "catalog-file", tmp_path / "catalog.zip")
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "catalog-status").text.startswith(("Loaded", "Error")))

    assert driver.find_element(By.ID, "catalog-status").text.startswith("Loaded"), driver.find_element(By.ID, "catalog-status").text
    assert set(strain_rows(driver)) == {"parent", "child"}
    assert strain_state(driver, "child") == "imported: done"
    assert driver.find_element(By.ID, "strain-parent").get_attribute("value") == "parent"
    assert combine(driver, "parent").startswith("Done") and first.startswith("Done")
    assert download_text(driver, "download-combined") == expected


@pytest.mark.parametrize("make,fragment", [
    (lambda z: None, "zip"),  # not a zip at all (the file is written as text below)
    (lambda z: z.writestr("readme.txt", "nothing"), "catalog.json"),
    (lambda z: z.writestr("catalog.json", "{not json"), "catalog"),
    (lambda z: z.writestr("catalog.json", json.dumps({"format": "other", "version": 1, "strains": []})), "format"),
    (lambda z: z.writestr("catalog.json", json.dumps({"format": "qmutfind-catalog", "version": 99, "strains": []})), "version"),
    (lambda z: z.writestr("catalog.json", json.dumps({"format": "qmutfind-catalog", "version": 1, "strains": [
        {"name": "x", "status": "done", "vcf": "../../etc/passwd"}]})), "file"),
])
def test_a_bad_catalog_file_is_refused_with_a_reason_and_changes_nothing(driver, server, strain_files, tmp_path, make, fragment):
    import zipfile

    path = tmp_path / "bad.zip"
    if fragment == "zip":
        path.write_text("this is not a zip file")
    else:
        with zipfile.ZipFile(path, "w") as z:
            make(z)
    open_site(driver, server)

    choose(driver, "catalog-file", path)
    WebDriverWait(driver, 30).until(lambda d: d.find_element(By.ID, "catalog-status").text.startswith(("Loaded", "Error")))

    status = driver.find_element(By.ID, "catalog-status").text
    assert status.startswith("Error") and fragment in status, status
    assert strain_rows(driver) == {}


# The reference panel: one question, three ways to answer it

def test_the_reference_panel_asks_where_the_genome_is_and_shows_only_the_matching_controls(driver, server):
    open_site(driver, server)
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "genome-menu").is_displayed())

    options = [o.get_attribute("value") for o in driver.find_elements(By.CSS_SELECTOR, "#reference-source option")]
    assert options == ["ncbi", "upload", "url"]
    assert driver.find_element(By.ID, "reference-source").get_attribute("value") == "ncbi"
    shown = lambda *ids: [driver.find_element(By.ID, i).is_displayed() for i in ids]
    # NCBI list: the menu; nothing of the other two.
    assert shown("genome-menu", "reference", "gff", "url-reference", "url-gff") == [True, False, False, False, False]
    set_reference_mode(driver, "upload")
    assert shown("genome-menu", "reference", "gff", "url-reference", "url-gff") == [False, True, True, False, False]
    set_reference_mode(driver, "url")
    assert shown("genome-menu", "reference", "gff", "url-reference", "url-gff") == [False, False, False, True, True]
    assert driver.find_element(By.ID, "ref-url-load").is_displayed()


def test_a_reference_address_without_a_gff_loads_the_genome_alone(driver, server, inputs, example_site):
    _shutil.copy(inputs / "ref.fa.gz", EXAMPLE_DIR / "only.fa.gz")
    open_web(driver, server)
    choose(driver, "url-reference", f"{server}/dist-site/_test_example/only.fa.gz")

    driver.find_element(By.ID, "ref-url-load").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))

    assert driver.find_element(By.ID, "genome-status").text.startswith("Loaded"), driver.find_element(By.ID, "genome-status").text
    assert source_text(driver, "reference").startswith("From the address") and "only.fa.gz" in source_text(driver, "reference")
    assert source_text(driver, "gff").startswith("Not set")


def test_a_reference_address_needs_the_fasta_address(driver, server):
    open_site(driver, server)
    set_reference_mode(driver, "url")

    driver.find_element(By.ID, "ref-url-load").click()

    assert "address" in driver.find_element(By.ID, "genome-status").text.lower()


def test_a_reference_and_gene_file_from_https_addresses_load_and_say_where_they_came_from(driver, server, inputs, example_site):
    expected = cli_results(inputs)
    _shutil.copy(inputs / "ref.fa.gz", EXAMPLE_DIR / "refdata")  # no .gz in the name, but gzip inside
    open_site(driver, server)
    choose(driver, "url-reference", f"{server}/dist-site/_test_example/refdata")
    driver.find_element(By.ID, "url-gff").send_keys(f"{server}/dist-site/_test_example/genes.gff3")

    driver.find_element(By.ID, "ref-url-load").click()
    WebDriverWait(driver, 60).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith(("Loaded", "Error")))

    assert driver.find_element(By.ID, "genome-status").text.startswith("Loaded"), driver.find_element(By.ID, "genome-status").text
    assert source_text(driver, "reference").startswith("From the address") and "refdata" in source_text(driver, "reference")
    assert source_text(driver, "gff").startswith("From the address") and "genes.gff3" in source_text(driver, "gff")
    add_sample(driver, "s1", "local", inputs / "reads.fq")
    assert run_samples(driver)["s1"].startswith("Done")
    assert download_text(driver, "download-vcf") == expected["vcf"]


@pytest.mark.parametrize("address,fragment", [
    ("http://example.org/ref.fa.gz", "https"),
    ("ftp://example.org/ref.fa.gz", "https"),
    ("https://user:secret@example.org/ref.fa.gz", "user name"),
    ("not a url", "address"),
])
def test_an_unsafe_or_malformed_reference_address_is_refused_before_any_request(driver, server, address, fragment):
    open_site(driver, server)
    choose(driver, "url-reference", address)

    driver.find_element(By.ID, "ref-url-load").click()
    WebDriverWait(driver, 20).until(lambda d: d.find_element(By.ID, "genome-status").text.startswith("Error"))

    assert fragment in driver.find_element(By.ID, "genome-status").text


# A cached old file must never mix with a new page: every file the page loads is named with the build stamp

def test_the_files_of_the_page_carry_one_build_stamp(driver, server):
    import re

    open_site(driver, server)
    html = driver.page_source
    stamps = set(re.findall(r'(?:style\.css|theme\.js|app\.js)\?v=([0-9a-f]{10})\b', html))
    assert len(stamps) == 1, stamps
    stamp = stamps.pop()
    get = lambda path: driver.execute_async_script("const cb = arguments[arguments.length - 1]; fetch(arguments[0]).then(r => r.text()).then(cb);", f"{server}/{SITE.rsplit('/', 1)[0]}/{path}".replace("//index.html", "/index.html"))
    app, worker = get("app.js"), get("worker.js")
    assert "__BUILD__" not in html and "__BUILD__" not in app and "__BUILD__" not in worker
    assert app.count(f"worker.js?v={stamp}") == 2  # both places that start a worker
    assert f"./pkg-web/qmutfind_wasm.js?v={stamp}" in worker and f"./vendor/fflate.js?v={stamp}" in worker
    assert f"qmutfind_wasm_bg.wasm?v={stamp}" in worker


def test_the_scripts_the_page_imports_are_loaded_with_the_stamp(driver, server):
    open_site(driver, server)
    names = driver.execute_script("return performance.getEntriesByType('resource').map(e => e.name);")

    assert any("app.js?v=" in n for n in names) and any("vendor/fflate.js?v=" in n for n in names), names
    assert not any(n.endswith("vendor/fflate.js") for n in names), names
