import csv

from sentinel_refbased_detect.align import longest_called_run, parse_tag
from sentinel_refbased_detect.calls import final_calls
from sentinel_refbased_detect.candidates import rank_depth, select_candidates
from sentinel_refbased_detect.kraken import library_size
from sentinel_refbased_detect.params import Params
from sentinel_refbased_detect.summarize import SUMMARY_COLS, detection_rule

P = Params()


def region(longest, n_blocks, covered=None, unique_reads=0, n_loci=0):
    return {"longest_block_bp": longest, "n_blocks": n_blocks,
            "covered_bp_min_depth": covered if covered is not None else longest,
            "unique_reads": unique_reads, "n_loci": n_loci}


def test_detection_rule_contiguous():
    det, why, _, rule = detection_rule([region(1000, 1)], P)
    assert det and rule == "contiguous" and why.startswith("contiguous 1000 bp")


def test_detection_rule_two_regions():
    det, why, n, rule = detection_rule(
        [region(250, 1), region(300, 1), region(0, 0)], P)
    assert det and n == 2 and rule == "multi-region" and "2 regions" in why


def test_two_reads_at_separate_loci_detect():
    """The weakest rule: evidence that is not one piece of genome read twice."""
    det, _, _, rule = detection_rule([region(120, 0, unique_reads=2, n_loci=2)], P)
    assert det and rule == "two-position"


def test_reads_piled_on_one_locus_do_not_detect():
    """Ten reads stacked on one spot are one observation, not ten."""
    det, why, _, rule = detection_rule(
        [region(140, 0, covered=140, unique_reads=10, n_loci=1)], P)
    assert not det and rule == "none"
    assert "10 unique read(s) at 1 position(s)" in why


def test_detection_rule_unsegmented_short():
    det, why, _, rule = detection_rule([region(999, 3, covered=5000)], P)
    assert not det and rule == "none"
    assert why.startswith("longest contiguous 999 bp (5000 bp covered in all), "
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


def test_ref_span_counts_reference_bases_only():
    from sentinel_refbased_detect.align import ref_span
    assert ref_span("150M") == 150
    assert ref_span("10S140M") == 140          # soft-clipped bases are not on the reference
    assert ref_span("100M5I45M") == 145        # an insertion consumes no reference
    assert ref_span("100M5D45M") == 150        # a deletion does


def test_overlapping_reads_count_as_one_locus(tmp_path):
    """What separates evidence from repetition, on a real BAM."""
    import shutil
    import subprocess
    if not (shutil.which("samtools")):
        import pytest
        pytest.skip("samtools not installed")
    from sentinel_refbased_detect.align import reads_by_ref

    sam = tmp_path / "x.sam"
    reads = [
        # Three reads piled on one spot: one locus.
        (100, "r1"), (120, "r2"), (140, "r3"),
        # Two more, far away and overlapping each other: a second locus.
        (5000, "r4"), (5080, "r5"),
        # One alone: a third.
        (9000, "r6"),
    ]
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:20000\n")
        for pos, name in reads:
            fh.write("%s\t0\tref\t%d\t60\t150M\t*\t0\t0\t%s\t%s\n"
                     % (name, pos, "A" * 150, "I" * 150))
    bam = tmp_path / "x.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)

    n, mapq, starts, loci = reads_by_ref(bam)["ref"]
    assert n == 6 and starts == 6               # six reads, six distinct starts
    assert loci == 3                            # but only three separate places
    assert mapq == 60.0
