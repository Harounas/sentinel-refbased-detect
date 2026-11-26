import csv

from sentinel_refbased_detect.align import longest_called_run, parse_tag
from sentinel_refbased_detect.calls import final_calls
from sentinel_refbased_detect.candidates import rank_depth, select_candidates
from sentinel_refbased_detect.kraken import library_size
from sentinel_refbased_detect.params import Params
from sentinel_refbased_detect.summarize import SUMMARY_COLS, detection_rule

P = Params()


def region(longest, n_blocks, covered=None):
    return {"longest_block_bp": longest, "n_blocks": n_blocks,
            "covered_bp_min_depth": covered if covered is not None else longest}


def test_detection_rule_contiguous():
    det, why, _ = detection_rule([region(1000, 1)], P)
    assert det and why.startswith("contiguous 1000 bp")


def test_detection_rule_two_regions():
    det, why, n = detection_rule([region(250, 1), region(300, 1), region(0, 0)], P)
    assert det and n == 2 and "2 regions" in why


def test_detection_rule_unsegmented_short():
    det, why, _ = detection_rule([region(999, 3, covered=5000)], P)
    assert not det
    assert why == ("longest contiguous 999 bp (5000 bp covered in all), "
                   "1 region(s) with a block >= 200 bp")


def test_parse_tag_with_underscored_code():
    assert parse_tag("3052310__LASV__0") == ("3052310", "LASV", "0")
    assert parse_tag("10239__Some__virus__12") == ("10239", "Some__virus", "12")


def test_longest_called_run(tmp_path):
    fa = tmp_path / "c.fa"
    fa.write_text(">x\nNNACGTNNACGTACGTN\nAA\n")
    assert longest_called_run(fa) == (9, 8, "ACGTACGT")


def test_rank_depth():
    assert rank_depth("G") < rank_depth("S") < rank_depth("S1")
    assert rank_depth("F") < rank_depth("G")
    assert rank_depth("X") == -1.0


REPORT = """\
 10.00\t1000\t1000\tU\t0\tunclassified
 90.00\t9000\t10\tR\t1\troot
  5.00\t500\t0\tG\t11618\t    Mammarenavirus
  4.00\t400\t400\tS\t11620\t      Lassa mammarenavirus
  1.00\t100\t0\tF\t11617\t    Arenaviridae
"""


def _write(path, rows, header=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        if header:
            w.writerow(header)
        w.writerows(rows)
    return path


def test_select_candidates_legacy_alias(tmp_path):
    rep = tmp_path / "r.txt"
    rep.write_text(REPORT)
    refs = _write(tmp_path / "refs.tsv", [["3052310", "LASV", "Lassa virus Josiah", "NC_004296.1:NC_004297.1"]])
    aliases = _write(tmp_path / "a.tsv", [["3052310", "3052310", "self"],
                                          ["11620", "3052310", "legacy"],
                                          ["11617", "3052310", "descendant:1"]],
                     ["report_taxid", "table_taxid", "relation"])
    names = _write(tmp_path / "n.tsv", [], ["normalized_name", "table_taxid"])
    levels = _write(tmp_path / "l.tsv", [["3052310", "3052310", "Mammarenavirus lassaense",
                                          "3052310,11620", "11618", "Mammarenavirus", "11618", "3052310"]],
                    ["table_taxid", "species_taxid", "species_name", "species_lookup",
                     "genus_taxid", "genus_name", "genus_lookup", "fallback_taxid"])
    out = tmp_path / "cand.tsv"
    select_candidates("S1", rep, refs, aliases, names, levels, out, P)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert len(rows) == 1
    r = rows[0]
    # The family line (F) is above min_match_rank G and must not win.
    assert r["match_basis"] == "taxid_alias" and r["report_taxid"] == "11620"
    assert r["kraken_reads"] == "400" and r["n_segments"] == "2"
    assert r["species_direct_reads"] == "400" and r["genus_clade_reads"] == "500"


def test_library_size(tmp_path):
    rep = tmp_path / "r.txt"
    rep.write_text(REPORT)
    assert library_size(rep, paired=True) == 20000
    assert library_size(rep, paired=False) == 10000


def test_final_calls_control_background(tmp_path):
    def row(sample, reads, det="detected"):
        r = {c: "" for c in SUMMARY_COLS}
        r.update(sample_id=sample, taxid="11250", virus_code="RSV", strain_name="RSV",
                 total_unique_reads_mapped=str(reads), detection=det,
                 detection_basis="x", longest_covered_block_bp="1500")
        return r
    summary = tmp_path / "s.tsv"
    with open(summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t")
        w.writeheader()
        w.writerows([row("NTC", 100), row("A", 500), row("B", 5000), row("C", 3, "not detected")])
    final, bg = tmp_path / "f.tsv", tmp_path / "bg.tsv"
    lib = {"NTC": 1_000_000, "A": 1_000_000, "B": 1_000_000, "C": 1_000_000}
    final_calls(summary, final, bg, P, library_sizes=lib)
    calls = {r["sample_id"]: r["final_call"] for r in csv.DictReader(open(final), delimiter="\t")}
    assert calls == {"NTC": "control sample", "A": "suspect: control background",
                     "B": "detected, unconfirmed", "C": "not detected"}
