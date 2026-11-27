"""The weighted specificity check, on reads whose origin is known.

Built from synthetic genomes so the right answer is not in doubt: a target, a
near-identical strain of the same taxon, and a sibling species that shares one
conserved stretch with the target.
"""
import gzip
import random
import shutil
import subprocess

import pytest

from sentinel_refbased_detect.evidence import reference_specificity
from sentinel_refbased_detect.utils import run

pytestmark = pytest.mark.skipif(
    not (shutil.which("bwa") and shutil.which("samtools")),
    reason="bwa/samtools not installed")

COMP = str.maketrans("ACGT", "TGCA")


def _seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng, seq, rate):
    return "".join(rng.choice([x for x in "ACGT" if x != c]) if rng.random() < rate else c
                   for c in seq)


def _fasta(path, name, seq):
    path.write_text(">%s\n%s\n" % (name, "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70))))


def _bam(tmp_path, ref_fa, reads, out_dir):
    """Align reads to ref and leave a unique.bam, as align_and_consensus does."""
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ref_fa, out_dir / "ref.fasta")
    fq = tmp_path / "reads.fastq"
    fq.write_text("".join("@r%d\n%s\n+\n%s\n" % (i, s, "I" * len(s))
                          for i, s in enumerate(reads)))
    run(["bwa", "index", "ref.fasta"], cwd=out_dir, stderr=out_dir / "log")
    run(["samtools", "faidx", "ref.fasta"], cwd=out_dir, stderr=out_dir / "log")
    sam = out_dir / "aln.sam"
    run(["bwa", "mem", str(out_dir / "ref.fasta"), fq], stdout=sam, stderr=out_dir / "log")
    run(["samtools", "sort", "-o", str(out_dir / "unique.bam"), str(sam)], stderr=out_dir / "log")
    run(["samtools", "index", str(out_dir / "unique.bam")], stderr=out_dir / "log")


@pytest.fixture
def world(tmp_path):
    rng = random.Random(7)
    conserved = _seq(rng, 1200)
    target = _seq(rng, 6000) + conserved + _seq(rng, 6000)
    strain = _mutate(rng, target, 0.004)          # same taxon, near identical
    sibling = _seq(rng, 6000) + conserved + _seq(rng, 6000)   # shares only the conserved part

    sample = tmp_path / "sample"
    (sample / "refs").mkdir(parents=True)
    _fasta(sample / "refs" / "900__TGT__0.fasta", "target", target)
    # The panel lives outside the sample, as it does in a real run: these
    # genomes are offered as decoys, not as the sample's own candidates.
    panel = tmp_path / "panel"
    panel.mkdir()
    _fasta(panel / "panel_strain_X1.fasta", "strain", strain)
    _fasta(panel / "panel_sibling_S1.fasta", "sibling", sibling)
    return tmp_path, sample, rng, target, conserved


def _reads(rng, source, n, start=0, span=None, length=150):
    span = span or (len(source) - length)
    return [source[p:p + length] for p in
            (rng.randrange(start, start + span) for _ in range(n))]


def test_reads_unique_to_target_are_attributed_to_it(world, tmp_path):
    root, sample, rng, target, _ = world
    reads = _reads(rng, target, 30, start=0, span=5000)       # the unique part
    _bam(root, sample / "refs" / "900__TGT__0.fasta", reads, sample / "align" / "900__TGT__0")
    r = reference_specificity(sample / "align" / "900__TGT__0", sample, "900__TGT__0",
                              root / "work",
                              decoys=[root / "panel" / "panel_sibling_S1.fasta"],
                              tag_taxon={"panel_sibling_S1": "Sibling species"})
    assert r["pct_target"] == 100.0
    assert r["n_unique_target"] == 30 and r["n_unique_other"] == 0


def test_conserved_reads_are_shared_out_not_lost(world, tmp_path):
    """Reads from a stretch the sibling shares are ambiguous, not a failure.

    With unique reads backing the target, the ambiguous ones are apportioned
    to it, so a true finding keeps a high share instead of being halved.
    """
    root, sample, rng, target, conserved = world
    unique_reads = _reads(rng, target, 20, start=0, span=5000)
    shared_reads = _reads(rng, target, 10, start=6000, span=1000)
    _bam(root, sample / "refs" / "900__TGT__0.fasta", unique_reads + shared_reads,
         sample / "align" / "900__TGT__0")
    r = reference_specificity(sample / "align" / "900__TGT__0", sample, "900__TGT__0",
                              root / "work",
                              decoys=[root / "panel" / "panel_sibling_S1.fasta"],
                              tag_taxon={"panel_sibling_S1": "Sibling species"})
    assert r["n_ambiguous"] > 0                      # the conserved stretch is ambiguous
    assert r["n_unique_target"] >= 15
    assert r["pct_target"] > 90                      # ambiguity does not sink a true finding
    assert "this taxon" in r["reads_by_taxon"]


def test_strain_of_same_taxon_counts_as_on_target(world, tmp_path):
    """A read fitting another strain of the same taxon is not evidence against."""
    root, sample, rng, target, _ = world
    reads = _reads(rng, target, 25, start=0, span=5000)
    _bam(root, sample / "refs" / "900__TGT__0.fasta", reads, sample / "align" / "900__TGT__0")
    r = reference_specificity(sample / "align" / "900__TGT__0", sample, "900__TGT__0",
                              root / "work",
                              same_taxon_tags=["panel_strain_X1"],
                              decoys=[root / "panel" / "panel_strain_X1.fasta",
                                      root / "panel" / "panel_sibling_S1.fasta"],
                              tag_taxon={"panel_sibling_S1": "Sibling species"})
    assert r["pct_target"] == 100.0
    assert r["n_competitors"] == 1          # only the sibling truly competes


def test_divergent_lineage_is_judged_against_its_own_clade(tmp_path):
    """A Lassa-like case: reads from a clade far from the chosen reference.

    Mapped against a distant lineage they look poor; the panel holds their own
    clade, so identity is measured there and the read is not failed merely for
    belonging to an under-represented lineage.
    """
    rng = random.Random(11)
    clade_a = _seq(rng, 12000)                 # the reference the pipeline picked
    clade_b = _mutate(rng, clade_a, 0.18)      # a lineage ~82% identical to it

    sample = tmp_path / "sample"
    (sample / "refs").mkdir(parents=True)
    _fasta(sample / "refs" / "900__LASV__0.fasta", "clade_a", clade_a)
    panel = tmp_path / "panel"
    panel.mkdir()
    _fasta(panel / "panel_strain_B.fasta", "clade_b", clade_b)

    reads = _reads(rng, clade_b, 25, start=0, span=10000)   # truly from clade B
    _bam(tmp_path, sample / "refs" / "900__LASV__0.fasta", reads,
         sample / "align" / "900__LASV__0")
    r = reference_specificity(sample / "align" / "900__LASV__0", sample, "900__LASV__0",
                              tmp_path / "work",
                              same_taxon_tags=["panel_strain_B"],
                              decoys=[panel / "panel_strain_B.fasta"])
    # Reads are near-identical to their own clade, which the panel supplies...
    assert r["pct_identity_best"] > 95
    # ...and the taxon still owns them.
    assert r["pct_target"] == 100.0


def test_background_gives_an_unlisted_virus_somewhere_true_to_go(tmp_path):
    """The case the priority table cannot cover.

    Reads really come from a virus nobody listed. Against the table alone they
    are attributed to whatever they were aligned against; with the background
    layer they go to the virus they actually came from.
    """
    rng = random.Random(31)
    listed = _seq(rng, 12000)          # the priority pathogen the reads hit
    unlisted = _seq(rng, 12000)        # what the reads are really from

    sample = tmp_path / "sample"
    (sample / "refs").mkdir(parents=True)
    _fasta(sample / "refs" / "900__LISTED__0.fasta", "listed", listed)

    bg_dir = tmp_path / "bg"
    bg_dir.mkdir()
    bg_fa = bg_dir / "background.fasta"
    _fasta(bg_fa, "NC_999.1", unlisted)
    (bg_dir / "background.taxa.tsv").write_text(
        "sequence_id\ttaxid\tname\nNC_999.1\t555\tSome unlisted virus\n")
    run(["bwa", "index", str(bg_fa)], stderr=bg_dir / "log")

    # A handful of reads from the unlisted virus, plus a few that genuinely
    # match the listed one, so the comparison is not trivially one-sided.
    reads = _reads(rng, unlisted, 20, span=10000) + _reads(rng, listed, 5, span=10000)
    _bam(tmp_path, sample / "refs" / "900__LISTED__0.fasta", reads,
         sample / "align" / "900__LISTED__0")

    from sentinel_refbased_detect.evidence import load_background
    r = reference_specificity(sample / "align" / "900__LISTED__0", sample,
                              "900__LISTED__0", tmp_path / "work",
                              background=load_background(bg_fa))
    assert "Some unlisted virus" in r["reads_by_taxon"]
    assert r["n_competitors"] >= 1
    # Most reads belong elsewhere, so the listed pathogen must not keep them.
    assert r["pct_target"] < 50
