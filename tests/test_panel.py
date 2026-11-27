"""Ways a virus can silently fall out of the panel.

Every step after the genomes are fetched makes the panel smaller: clustering,
the representative cap, dedup across shards. Each of these is a check that the
shrinking stopped short of removing a taxon altogether.
"""
import csv

from sentinel_refbased_detect import panel


def _row(acc, role="sibling", taxid="100", related="1", name="Some virus"):
    return {"taxid": taxid, "code": "X", "name": name, "accessions": acc,
            "role": role, "related_to": related, "note": ""}


def test_dedup_merges_guards_instead_of_dropping_a_row(tmp_path):
    """The same genome reached from two genera guards both of them.

    Keeping only the first copy would take that decoy away from every priority
    taxon but one, which is a gap nothing else would report.
    """
    out = tmp_path / "panel.tsv"
    panel.write_panel([_row("NC_1.1", related="11620"),
                       _row("NC_1.1", related="11621"),
                       _row("NC_2.1", related="11620")], out)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert len(rows) == 2
    merged = next(r for r in rows if r["accessions"] == "NC_1.1")
    assert sorted(merged["related_to"].split(",")) == ["11620", "11621"]


def test_same_genome_in_both_roles_is_kept_once_each(tmp_path):
    out = tmp_path / "panel.tsv"
    panel.write_panel([_row("NC_1.1", role="sibling"),
                       _row("NC_1.1", role="strain")], out)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert {r["role"] for r in rows} == {"sibling", "strain"}


class _FakeEntrez:
    """Returns the one genome asked for, so no network is needed."""

    def __init__(self, seqs):
        self.seqs = seqs

    def efetch_fasta(self, accs):
        return "".join(">%s\n%s\n" % (a, self.seqs[a]) for a in accs)


def test_a_species_with_one_genome_is_still_included(tmp_path):
    """A lone deposit is not a clustering problem, but it is the only genome
    that species has; dropping it leaves the species out of the panel."""
    recs = [{"acc": "NC_9.1", "taxid": "100", "slen": 9000,
             "title": "Some virus, complete genome"}]
    rows = panel.strains_for_species(
        "100", {"code": "SV", "name": "Some virus", "guards": ["100"]},
        recs, _FakeEntrez({"NC_9.1": "ACGT" * 2250}), tmp_path)
    assert len(rows) == 1
    assert rows[0]["accessions"] == "NC_9.1"
    assert "only genome available" in rows[0]["note"]


def test_a_segment_with_one_sequence_is_not_dropped(tmp_path):
    """Lassa deposits L and S separately. A species with several S segments and
    a single L must still get an L representative."""
    recs = [{"acc": "L1.1", "taxid": "100", "slen": 7000,
             "title": "Virus isolate A segment L, complete sequence"},
            {"acc": "S1.1", "taxid": "100", "slen": 3400,
             "title": "Virus isolate A segment S, complete sequence"},
            {"acc": "S2.1", "taxid": "100", "slen": 3400,
             "title": "Virus isolate B segment S, complete sequence"}]
    seqs = {"L1.1": "ACGTTGCA" * 875, "S1.1": "TTGACCAG" * 425, "S2.1": "TTGACCAG" * 425}
    rows = panel.strains_for_species(
        "100", {"code": "V", "name": "Virus", "guards": ["100"]},
        recs, _FakeEntrez(seqs), tmp_path)
    parts = {r["note"].split(";")[0] for r in rows}
    assert "segment L" in parts and "segment S" in parts


def test_verify_names_the_pathogens_left_without_decoys(tmp_path):
    ref = tmp_path / "refs.tsv"
    ref.write_text("11620\tLASV\tLassa\tNC_1.1\n3052505\tMARV\tMarburg\tNC_2.1\n")
    pan = tmp_path / "panel.tsv"
    panel.write_panel([_row("NC_8.1", related="11620", name="Mopeia virus")], pan)
    rows, missing = panel.verify(pan, ref, out_path=tmp_path / "coverage.tsv")
    assert [m[1] for m in missing] == ["MARV"]
    by_code = {r["code"]: r for r in rows}
    assert by_code["LASV"]["status"] == "ok" and by_code["MARV"]["status"] == "NO DECOYS"


def test_efetch_batches_large_accession_lists(monkeypatch):
    """A few hundred accessions in a URL make NCBI answer 414.

    The real run asked for 1000 West Nile genomes, which is an 11,000-character
    URL; every such shard failed and, because a failing shard writes an empty
    file to keep the scatter whole, the workflow reported success with 45 of
    83 species missing.
    """
    from sentinel_refbased_detect.ncbi import Entrez

    seen = []

    class _Resp:
        def __init__(self, text):
            self._t = text

        def read(self):
            return self._t.encode()

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(req, timeout=None):
        # The id list must travel in the body, not the query string.
        assert req.data, "accessions were sent in the URL, not by POST"
        body = req.data.decode()
        ids = [p for p in body.split("&") if p.startswith("id=")][0]
        n = ids.count("%2C") + 1
        seen.append(n)
        return _Resp("".join(">acc%d\nACGT\n" % i for i in range(n)))

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    e = Entrez("test@example.org")
    e.gap = 0
    text = e.efetch_fasta(["NC_%06d.1" % i for i in range(650)], chunk=200)
    assert text.count(">") == 650
    assert seen == [200, 200, 200, 50]
