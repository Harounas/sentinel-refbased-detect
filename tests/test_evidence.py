import csv

import pytest

from sentinel_refbased_detect import evidence
from sentinel_refbased_detect.params import Params
from sentinel_refbased_detect.summarize import SUMMARY_COLS

P = Params()


def test_aligned_length_counts_ref_and_insert_bases():
    assert evidence.aligned_length("10S140M") == 140
    assert evidence.aligned_length("100M2I48M") == 150
    assert evidence.aligned_length("100M5D45M") == 150


def _summary(tmp_path, rows):
    path = tmp_path / "summary.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            full = {c: "" for c in SUMMARY_COLS}
            full.update(r)
            w.writerow(full)
    return path


def _row(sample, taxid="11620", **kw):
    r = {"sample_id": sample, "taxid": taxid, "virus_code": "LASV",
         "virus_importance": "priority", "detection": "not detected",
         "kraken_reads": "40", "total_unique_reads_mapped": "0",
         "breadth_of_coverage": "0.0", "best_reference": "", "detection_basis": "x"}
    r.update(kw)
    return r


def test_row_with_no_mapped_reads(tmp_path):
    s = _summary(tmp_path, [_row("S1")])
    v = evidence.review(s, tmp_path / "rev.tsv", tmp_path, P)
    assert v[("S1", "11620")]["evidence_verdict"] == "no mapped reads"


def test_detected_and_discovered_rows_are_not_reviewed(tmp_path):
    s = _summary(tmp_path, [
        _row("S1", detection="detected"),
        _row("S2", taxid="99", virus_importance="discovered"),
        _row("NTC"),                      # a control is never a finding
        _row("S3", kraken_reads="0"),     # under evidence_min_kraken_reads
    ])
    v = evidence.review(s, tmp_path / "rev.tsv", tmp_path, P)
    assert v == {}


def test_control_and_background_criteria(tmp_path):
    """A virus sitting at control level, and a run-mate 1000x louder, both fail."""
    rows = [_row("S1", total_unique_reads_mapped="12"),
            _row("NTC", total_unique_reads_mapped="500"),
            _row("S9", total_unique_reads_mapped="20000")]
    s = _summary(tmp_path, rows)
    lib = {"S1": 1_000_000, "NTC": 1_000_000, "S9": 1_000_000}
    evidence.review(s, tmp_path / "rev.tsv", tmp_path, P, library_sizes=lib)
    out = {r["sample_id"]: r for r in csv.DictReader(open(tmp_path / "rev.tsv"), delimiter="\t")}
    r = out["S1"]
    assert r["control_verdict"] == "fail" and r["control_sample"] == "NTC"
    assert r["background_verdict"] == "fail" and r["max_other_sample"] == "S9"
    assert r["hopping_ratio"] == "1666.7"


def test_annotate_marks_unreviewed_rows(tmp_path):
    s = _summary(tmp_path, [_row("S1"), _row("S2", detection="detected")])
    v = evidence.review(s, tmp_path / "rev.tsv", tmp_path, P)
    out = evidence.annotate(s, v, tmp_path / "final.tsv")
    rows = {r["sample_id"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert rows["S1"]["evidence_verdict"] == "no mapped reads"
    assert rows["S2"]["evidence_verdict"] == "not reviewed"
    # The coverage call itself is never rewritten.
    assert rows["S1"]["detection"] == "not detected"
    assert rows["S2"]["detection"] == "detected"


def test_specificity_is_not_scored_when_nothing_competes(tmp_path, monkeypatch):
    """A pathogen with no relatives, as ASFV has none, must not pass by default.

    Asfarviridae holds one genus holding one species, so no sibling genome
    exists to offer. Attribution is then trivially 100% and would read as a
    pass while having tested nothing.
    """
    from sentinel_refbased_detect import evidence as ev

    monkeypatch.setattr(ev, "alignment_metrics", lambda *a, **k: {
        "n_reads": 10, "pct_identity_mean": 99.0, "pct_identity_median": 99.0,
        "mean_mapq": 60.0, "pct_multimapping": 0.0, "n_distinct_starts": 10,
        "n_bins_covered": 5, "n_bins_total": 20, "pct_busiest_bin": 30.0,
        "n_regions_with_reads": 1, "reference_bp": 10000})
    monkeypatch.setattr(ev, "reference_specificity", lambda *a, **k: {
        "pct_target": 100.0, "pct_identity_best": 99.0, "n_unique_target": 10,
        "n_unique_other": 0, "n_ambiguous": 0, "reads_by_taxon": "this taxon:10.0",
        "competing_best": "", "n_competitors": 0, "competitors": "",
        "reads_by_reference": ""})

    s = _summary(tmp_path, [_row("S1", taxid="10497", total_unique_reads_mapped="10",
                                 best_reference="10497__ASFV__0",
                                 breadth_of_coverage="0.02")])
    ev.review(s, tmp_path / "rev.tsv", tmp_path, P)
    r = list(csv.DictReader(open(tmp_path / "rev.tsv"), delimiter="\t"))[0]
    assert r["pct_reads_attributed_to_target"] == "100.0"
    assert r["specificity_verdict"] == "not assessable"
