"""Make site/genomes/genomes.json: the menu of reference genomes of the site.

The genomes are chosen by hand (ACCESSIONS below). Everything else comes from the NCBI RefSeq fungi assembly table
(https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt) and from HEAD requests to the two files of each
genome (the FASTA and the GFF3), which give the file sizes the page checks after each download.

    python3 make_genome_list.py --out ../site/genomes/genomes.json [--summary assembly_summary.txt]

A genome is kept only if it is the latest version, has an annotation (an annotation provider and a gene count) and is not
larger than --max-size. Anything else is reported and the script stops with an error, so the list is never silently
shorter than planned.
"""
import argparse
import concurrent.futures
import datetime
import json
import sys
import urllib.request

SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt"
ORIGIN = "https://qmutfind.fungalgenomes.org"

# Chosen by hand: model fungi, pathogens and plant pathogens that people run mutation analysis on. Each is a RefSeq
# genome with an annotation, at most 70 Mb.
ACCESSIONS = [
    "GCF_000002655.1",  # Aspergillus fumigatus Af293
    "GCF_000011425.1",  # Aspergillus nidulans FGSC A4
    "GCF_000002855.4",  # Aspergillus niger CBS 513.88
    "GCF_000184455.2",  # Aspergillus oryzae RIB40
    "GCF_009017415.1",  # Aspergillus flavus NRRL 3357
    "GCF_000149615.1",  # Aspergillus terreus NIH2624
    "GCF_000203795.2",  # Batrachochytrium dendrobatidis JAM81
    "GCF_000143535.2",  # Botrytis cinerea B05.10
    "GCF_000182965.3",  # Candida albicans SC5314
    "GCF_000026945.1",  # Candida dubliniensis CD36
    "GCF_000006335.3",  # Candida tropicalis MYA-3404
    "GCF_010111755.1",  # Nakaseomyces glabratus (Candida glabrata)
    "GCF_002759435.1",  # Candidozyma auris B8441
    "GCF_000182765.1",  # Lodderomyces parapsilosis CDC317 (Candida parapsilosis)
    "GCF_000143365.1",  # Chaetomium globosum CBS 148.51
    "GCF_000149335.2",  # Coccidioides immitis RS
    "GCF_018416015.2",  # Coccidioides posadasii Silveira
    "GCF_001672515.1",  # Colletotrichum higginsianum IMI 349063
    "GCF_000182895.1",  # Coprinopsis cinerea okayama7#130
    "GCF_000149245.1",  # Cryptococcus neoformans H99
    "GCF_000091045.1",  # Cryptococcus deneoformans JEC21
    "GCF_000185945.1",  # Cryptococcus gattii WM276
    "GCF_002954075.1",  # Cryptococcus deuterogattii R265
    "GCF_000006445.2",  # Debaryomyces hansenii CBS767
    "GCF_000091025.4",  # Eremothecium gossypii (Ashbya gossypii) ATCC 10895
    "GCF_000240135.3",  # Fusarium graminearum PH-1
    "GCF_000149955.1",  # Fusarium oxysporum f. sp. lycopersici 4287
    "GCF_000149555.1",  # Fusarium verticillioides 7600
    "GCF_000150115.1",  # Histoplasma capsulatum G186AR
    "GCF_000002515.2",  # Kluyveromyces lactis NRRL Y-1140
    "GCF_000027005.1",  # Komagataella phaffii GS115 (Pichia pastoris)
    "GCF_000143565.1",  # Laccaria bicolor S238N-H82
    "GCF_000181695.2",  # Malassezia globosa CBS 7966
    "GCF_000328475.2",  # Mycosarcoma maydis 521 (Ustilago maydis)
    "GCF_000182925.2",  # Neurospora crassa OR74A
    "GCF_000146915.1",  # Parastagonospora nodorum SN15
    "GCF_028828025.1",  # Penicillium rubens IBT 27055
    "GCF_000226545.1",  # Podospora anserina S mat+
    "GCF_000002495.2",  # Pyricularia oryzae 70-15 (Magnaporthe oryzae)
    "GCF_000149305.1",  # Rhizopus delemar RA 99-880
    "GCF_000146045.2",  # Saccharomyces cerevisiae S288C
    "GCF_002079055.1",  # Saccharomyces paradoxus CBS432
    "GCF_000143185.2",  # Schizophyllum commune H4-8
    "GCF_000002945.2",  # Schizosaccharomyces pombe 972h-
    "GCF_000146945.2",  # Sclerotinia sclerotiorum 1980
    "GCF_009556855.1",  # Talaromyces marneffei
    "GCF_000167675.1",  # Trichoderma reesei QM6a
    "GCF_000150675.1",  # Verticillium dahliae VdLs.17
    "GCF_001761485.1",  # Yarrowia lipolytica CLIB89(W29)
    "GCF_000219625.1",  # Zymoseptoria tritici IPO323
]


def parse_summary(text):
    """The assembly table as {accession: {column: value}}. Comment lines start with ##; the header starts with one #."""
    columns, rows = None, {}
    for line in text.splitlines():
        if line.startswith("##") or not line.strip():
            continue
        if line.startswith("#"):
            columns = line[1:].split("\t")
            continue
        values = line.split("\t")
        row = dict(zip(columns, values))
        rows[row["assembly_accession"]] = row
    return rows


def _int(text):
    try:
        return int(text)
    except ValueError:
        return 0


def make_entry(row):
    """One menu entry from a table row. File sizes are filled in later by `add_file_sizes`."""
    organism = row["organism_name"]
    strain = row["infraspecific_name"].split("=", 1)[-1] if "=" in row["infraspecific_name"] else ""
    label = organism if (not strain or strain in organism) else f"{organism} {strain}"
    base = row["ftp_path"].replace("ftp://", "https://", 1).rstrip("/")
    stem = f"{base}/{base.rsplit('/', 1)[-1]}"
    name = base.rsplit("/", 1)[-1]
    return {
        "accession": row["assembly_accession"],
        "label": label,
        "organism": organism,
        "strain": strain,
        "level": row["assembly_level"],
        "version_status": row["version_status"],
        "assembly": row["asm_name"],
        "released": row["seq_rel_date"],
        "genome_size": _int(row["genome_size"]),
        "genes": _int(row["total_gene_count"]),
        "protein_coding_genes": _int(row["protein_coding_gene_count"]),
        "annotation": " ".join(x for x in (row["annotation_provider"], row["annotation_name"], row["annotation_date"]) if x and x != "na"),
        "fasta": {"name": f"{name}_genomic.fna.gz", "url": f"{stem}_genomic.fna.gz", "bytes": 0},
        "gff": {"name": f"{name}_genomic.gff.gz", "url": f"{stem}_genomic.gff.gz", "bytes": 0},
    }


def check_entry(entry, max_size):
    """The reasons an entry cannot be used. An empty list means it can."""
    problems = []
    if entry["version_status"] != "latest":
        problems.append(f"not the latest version ({entry['version_status']})")
    if not entry["annotation"] or entry["genes"] <= 0:
        problems.append("no annotation (no annotation provider or no genes)")
    if entry["genome_size"] <= 0 or entry["genome_size"] > max_size:
        problems.append(f"genome size {entry['genome_size']} is outside 1 to {max_size}")
    return problems


def head(url):
    """(status, content length, Access-Control-Allow-Origin) of a HEAD request."""
    request = urllib.request.Request(url, method="HEAD", headers={"Origin": ORIGIN})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.status, int(response.headers.get("Content-Length", 0)), response.headers.get("Access-Control-Allow-Origin")


def add_file_sizes(entry):
    """Fill in the byte sizes of the two files. Returns the problems found (a missing file, no CORS header)."""
    problems = []
    for key in ("fasta", "gff"):
        try:
            status, size, cors = head(entry[key]["url"])
        except Exception as error:  # a network or HTTP error
            problems.append(f"{key}: {error}")
            continue
        entry[key]["bytes"] = size
        if status != 200 or size <= 0:
            problems.append(f"{key}: HTTP {status}, {size} bytes")
        if cors not in ("*", ORIGIN):
            problems.append(f"{key}: no Access-Control-Allow-Origin header, so a browser could not download it")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description="Make the genome menu of the site.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", help="a saved assembly_summary.txt (default: download it)")
    ap.add_argument("--max-size", type=int, default=70_000_000, help="largest genome in bases (default 70 Mb)")
    args = ap.parse_args(argv)
    if args.summary:
        text = open(args.summary).read()
    else:
        with urllib.request.urlopen(SUMMARY_URL, timeout=120) as response:
            text = response.read().decode()
    rows = parse_summary(text)
    entries, failures = [], []
    for accession in ACCESSIONS:
        if accession not in rows:
            failures.append(f"{accession}: not in the RefSeq fungi table")
            continue
        entry = make_entry(rows[accession])
        problems = check_entry(entry, args.max_size)
        if problems:
            failures.append(f"{accession} {entry['label']}: " + "; ".join(problems))
        else:
            entries.append(entry)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for entry, problems in zip(entries, pool.map(add_file_sizes, entries)):
            if problems:
                failures.append(f"{entry['accession']} {entry['label']}: " + "; ".join(problems))
    if failures:
        print("PROBLEMS:\n  " + "\n  ".join(failures), file=sys.stderr)
        return 1
    entries.sort(key=lambda e: e["label"].lower())
    document = {
        "source": "NCBI RefSeq fungi assembly table (" + SUMMARY_URL + ") and HEAD requests to the files",
        "generated": datetime.date.today().isoformat(),
        "max_genome_size": args.max_size,
        "genomes": entries,
    }
    with open(args.out, "w") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")
    print(f"wrote {len(entries)} genomes to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
