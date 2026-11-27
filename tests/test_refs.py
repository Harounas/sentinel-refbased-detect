from sentinel_refbased_detect import DEFAULT_REF_TSV
from sentinel_refbased_detect.candidates import load_refs
from sentinel_refbased_detect.cli import build_parser


def test_bundled_table_loads():
    assert DEFAULT_REF_TSV.exists()
    refs = load_refs(DEFAULT_REF_TSV)
    assert len(refs) > 250                       # ~288 taxids in V15
    assert refs["54290"][0][:3] == ("GBV-C", "GB virus C", ["NC_001710.1"])


def test_ref_tsv_defaults_to_bundled():
    args = build_parser().parse_args(["run", "-k", "db", "-i", "fq", "-o", "out"])
    assert args.ref_tsv == str(DEFAULT_REF_TSV)
