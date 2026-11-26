import gzip

import pytest

from sentinel_refbased_detect.samples import discover, sample_id_from_name
from sentinel_refbased_detect.utils import PipelineError


@pytest.mark.parametrize("name,expected", [
    ("1_04_24_0586_S67_R1_001.fastq.gz", ("1_04_24_0586", 1)),
    ("1-02-26-0109_S5_L001_R2_001.fastq.gz", ("1-02-26-0109", 2)),
    ("NEC_S96_R1_001.fastq.gz", ("NEC", 1)),
    ("sampleA_1.fq.gz", ("sampleA", 1)),
    ("sampleA_2.fq", ("sampleA", 2)),
    ("sampleB.R1.fastq", ("sampleB", 1)),
    ("1-02-26-0109.l1.AAF5NT5HV.1.cleaned_R1.fastq.gz", ("1-02-26-0109", 1)),
    ("single_end.fastq.gz", ("single_end", None)),
])
def test_sample_id_from_name(name, expected):
    assert sample_id_from_name(name) == expected


def _fq(path, n=1):
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write("@r%d\nACGT\n+\nIIII\n" % i)


def test_discover_pairs_and_lanes(tmp_path):
    for lane in ("L001", "L002"):
        _fq(tmp_path / f"S1_S1_{lane}_R1_001.fastq.gz")
        _fq(tmp_path / f"S1_S1_{lane}_R2_001.fastq.gz")
    _fq(tmp_path / "NTC_S2_R1_001.fastq.gz")
    _fq(tmp_path / "NTC_S2_R2_001.fastq.gz")
    _fq(tmp_path / "SE.fastq.gz")
    samples = {s.sample_id: s for s in discover(tmp_path)}
    assert set(samples) == {"S1", "NTC", "SE"}
    assert len(samples["S1"].r1) == 2 and samples["S1"].paired
    assert not samples["SE"].paired


def test_discover_unbalanced(tmp_path):
    _fq(tmp_path / "X_R1.fastq.gz")
    _fq(tmp_path / "X_L002_R1.fastq.gz")
    _fq(tmp_path / "X_R2.fastq.gz")
    with pytest.raises(PipelineError):
        discover(tmp_path)
