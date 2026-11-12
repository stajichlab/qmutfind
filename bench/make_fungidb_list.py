"""Make site/genomes/fungidb.json: the FungiDB organisms that the site menu lists.

FungiDB does not let a web page download its files (no Access-Control-Allow-Origin header, and the download needs an API
key), so the site only lists these organisms and shows the `qmutfind fetch` command that downloads one.

    python3 make_fungidb_list.py --out ../site/genomes/fungidb.json

The list comes from the public organism search of FungiDB (no key needed). Organisms without genes are left out.
"""
import argparse
import datetime
import json
import re
import sys
import urllib.parse
import urllib.request

URL = "https://fungidb.org/fungidb/service/record-types/organism/searches/GenomeDataTypes/reports/standard"
CONFIG = {"attributes": ["primary_key", "species", "genecount", "URLGenomeFasta"], "tables": [], "attributeFormat": "text"}


def parse_organisms(text):
    """(organisms, release) from the JSON of the organism search. The release (such as FungiDB-71) is in the file names."""
    records = json.loads(text)["records"]
    organisms, release = [], None
    for record in records:
        a = record["attributes"]
        match = re.search(r"/(FungiDB-\d+)_", a.get("URLGenomeFasta", ""))
        if match and release is None:
            release = match.group(1)
        try:
            genes = int(a.get("genecount") or 0)
        except ValueError:
            genes = 0
        if genes > 0:
            organisms.append({"name": a["primary_key"], "species": a.get("species", ""), "genes": genes})
    organisms.sort(key=lambda o: o["name"].lower())
    return organisms, release


def main(argv=None):
    ap = argparse.ArgumentParser(description="List the FungiDB organisms for the site menu.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--json", help="a saved organism search (default: download it)")
    args = ap.parse_args(argv)
    if args.json:
        text = open(args.json).read()
    else:
        url = URL + "?reportConfig=" + urllib.parse.quote(json.dumps(CONFIG))
        with urllib.request.urlopen(url, timeout=120) as response:
            text = response.read().decode()
    organisms, release = parse_organisms(text)
    document = {
        "source": "FungiDB organism search (" + URL + ")",
        "release": release,
        "generated": datetime.date.today().isoformat(),
        "organisms": organisms,
    }
    with open(args.out, "w") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")
    print(f"wrote {len(organisms)} organisms ({release}) to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
