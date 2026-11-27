"""The taxonomic decision: what the reads can and cannot say about species."""
import csv

import pytest

from sentinel_refbased_detect import rescue
from sentinel_refbased_detect.params import Params

P = Params()
TARGET = "__target__"
# The target is resolved by its species name, not by the placeholder label,
# which is what lets a tie with a sibling read as evidence for the genus.
GENUS = {"Lassa mammarenavirus": "Mammarenavirus", "Mopeia virus": "Mammarenavirus",
         "Human adenovirus": "Mastadenovirus", "Homo sapiens": None}


def genus_of(label):
    return GENUS.get(label)


def _read(hits):
    """hits: {label: (score, reference, position)}"""
    return hits


def _rows(per_read, control_verdict="pass"):
    return [{"_spec": {"per_read": per_read, "target_label": TARGET},
             "_row": {"species_name": "Lassa mammarenavirus"},
             "_report": {"sample_id": "S1", "taxid": "11620"},
             "control_verdict": control_verdict,
             "strain_name": "Lassa virus"}]


def test_reads_that_only_one_species_explains_give_a_species_call():
    per_read = {"r%d" % i: _read({TARGET: (150, "L", i * 500),
                                  "Mopeia virus": (90, "L", i * 500)})
                for i in range(6)}
    out = rescue.run(_rows(per_read), P, genus_lookup=genus_of)[0]
    assert out["taxonomic_outcome"] == rescue.OUTCOME_SPECIES
    assert out["reads_species_informative"] == 6
    assert out["best_species_positions"] == 6
    assert out["confidence"] == "high"


def test_reads_that_fit_both_relatives_give_a_genus_call():
    """The case the spec is built around: no read tells the species apart."""
    per_read = {"r%d" % i: _read({TARGET: (150, "L", i * 500),
                                  "Mopeia virus": (150, "L", i * 500)})
                for i in range(8)}
    out = rescue.run(_rows(per_read), P, genus_lookup=genus_of)[0]
    assert out["taxonomic_outcome"] == rescue.OUTCOME_GENUS
    assert out["reads_genus_informative"] == 8
    assert out["reads_species_informative"] == 0
    assert out["taxonomic_call"] == "Mammarenavirus sp."


def test_a_marginal_score_lead_does_not_buy_a_species():
    """One point of alignment score is not discrimination.

    With the default margin a one-point lead is a win, so the guard that has
    to hold is the advantage rule: 3 reads against 2 is not enough to name a
    species over its relative.
    """
    per_read = {}
    for i in range(3):
        per_read["a%d" % i] = _read({TARGET: (150, "L", i * 500),
                                     "Mopeia virus": (149, "L", i * 500)})
    for i in range(2):
        per_read["b%d" % i] = _read({TARGET: (149, "L", i * 700),
                                     "Mopeia virus": (150, "L", i * 700)})
    out = rescue.run(_rows(per_read), P, genus_lookup=genus_of)[0]
    assert out["taxonomic_outcome"] != rescue.OUTCOME_SPECIES


def test_both_species_with_their_own_spread_reads_is_a_mixed_infection():
    per_read = {}
    for i in range(4):
        per_read["a%d" % i] = _read({TARGET: (150, "L", i * 900),
                                     "Mopeia virus": (80, "L", i * 900)})
        per_read["b%d" % i] = _read({"Mopeia virus": (150, "L", 5000 + i * 900),
                                     TARGET: (80, "L", 5000 + i * 900)})
    out = rescue.run(_rows(per_read), P, genus_lookup=genus_of)[0]
    assert out["taxonomic_outcome"] == rescue.OUTCOME_MIXED
    assert "+" in out["taxonomic_call"]


def test_reads_piled_on_one_position_are_one_observation():
    """Twenty reads at one spot are a repeat, not twenty pieces of evidence."""
    per_read = {"r%d" % i: _read({TARGET: (150, "L", 4000),
                                  "Mopeia virus": (90, "L", 4000)})
                for i in range(20)}
    out = rescue.run(_rows(per_read), P, genus_lookup=genus_of)[0]
    assert out["best_species_positions"] == 1
    assert out["taxonomic_outcome"] != rescue.OUTCOME_SPECIES


def test_control_background_overrides_everything():
    per_read = {"r%d" % i: _read({TARGET: (150, "L", i * 500)}) for i in range(10)}
    out = rescue.run(_rows(per_read, control_verdict="fail"), P, genus_lookup=genus_of)[0]
    assert out["taxonomic_outcome"] == rescue.OUTCOME_INSUFFICIENT
    assert "control" in out["decision_reason"]


def test_report_has_every_section(tmp_path):
    per_read = {"r%d" % i: _read({TARGET: (150, "L", i * 500)}) for i in range(5)}
    rescue.run(_rows(per_read), P, genus_lookup=genus_of, out_path=tmp_path / "r.tsv")
    cols = list(csv.DictReader(open(tmp_path / "r.tsv"), delimiter="\t").fieldnames)
    for need in ("primary_detection", "reads_species_informative", "best_species",
                 "mean_mapq", "breadth_of_coverage", "competing_species_detail",
                 "control_sample", "taxonomic_call", "decision_reason"):
        assert need in cols
