import json

from make_fungidb_list import parse_organisms

SAMPLE = json.dumps({"records": [
    {"attributes": {"primary_key": "Aspergillus fumigatus Af293", "species": "Aspergillus fumigatus", "genecount": "10130",
                    "URLGenomeFasta": "http://FungiDB.org/a/service/raw-files/Current_Release/AfumigatusAf293/fasta/data/FungiDB-71_AfumigatusAf293_Genome.fasta"}},
    {"attributes": {"primary_key": "Candida albicans SC5314", "species": "Candida albicans", "genecount": "6455",
                    "URLGenomeFasta": "http://FungiDB.org/a/service/raw-files/Current_Release/CalbicansSC5314/fasta/data/FungiDB-71_CalbicansSC5314_Genome.fasta"}},
    {"attributes": {"primary_key": "No genes", "species": "X", "genecount": "0", "URLGenomeFasta": "http://x/FungiDB-71_X_Genome.fasta"}},
    {"attributes": {"primary_key": "Odd <name> & \"quotes\"", "species": "Y y", "genecount": "50", "URLGenomeFasta": "http://x/FungiDB-71_Y_Genome.fasta"}},
]})


def test_organisms_are_listed_by_name_with_gene_counts_and_the_release_is_taken_from_the_file_names():
    organisms, release = parse_organisms(SAMPLE)

    assert release == "FungiDB-71"
    assert [o["name"] for o in organisms][:2] == ["Aspergillus fumigatus Af293", "Candida albicans SC5314"]
    assert organisms[0]["genes"] == 10130 and organisms[0]["species"] == "Aspergillus fumigatus"


def test_an_organism_without_genes_is_left_out_because_it_has_no_annotation():
    organisms, _ = parse_organisms(SAMPLE)

    assert "No genes" not in [o["name"] for o in organisms]


def test_the_organisms_are_sorted_by_name_ignoring_case():
    organisms, _ = parse_organisms(SAMPLE)

    names = [o["name"] for o in organisms]
    assert names == sorted(names, key=str.lower)
