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
