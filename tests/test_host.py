import gzip
import random
import shutil

import pytest

from sentinel_refbased_detect.host import deplete, host_index

pytestmark = pytest.mark.skipif(
    not (shutil.which("minimap2") and shutil.which("samtools")),
    reason="minimap2/samtools not installed")

COMP = str.maketrans("ACGT", "TGCA")


def _genome(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _write_pairs(fh1, fh2, genome, rng, n, tag, start=0):
    for i in range(n):
        p = rng.randrange(0, len(genome) - 400)
        frag = genome[p:p + 400]
        name = "%s%d" % (tag, start + i)
        fh1.write("@%s/1\n%s\n+\n%s\n" % (name, frag[:150], "I" * 150))
        fh2.write("@%s/2\n%s\n+\n%s\n" % (name, frag[-150:].translate(COMP)[::-1], "I" * 150))


def test_host_pairs_removed_others_kept(tmp_path):
    rng = random.Random(1)
    host, virus = _genome(rng, 20000), _genome(rng, 5000)
    (tmp_path / "host.fa").write_text(">host\n%s\n" % host)
    r1, r2 = tmp_path / "s_R1.fastq.gz", tmp_path / "s_R2.fastq.gz"
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        _write_pairs(f1, f2, host, rng, 200, "h")
        _write_pairs(f1, f2, virus, rng, 50, "v")

    idx = host_index(tmp_path / "host.fa", tmp_path / "cache", threads=1)
    assert idx.exists()
    o1, o2, st = deplete("s", r1, r2, idx, tmp_path / "out", threads=1)
    assert st["reads_in"] == 500 and st["reads_kept"] == 100 and st["pct_host"] == "80.00"
    with gzip.open(o1, "rt") as fh:
        names = {l[1:].strip() for i, l in enumerate(fh) if i % 4 == 0}
    assert all(n.startswith("v") for n in names) and len(names) == 50
    # Resume reuses the result without recomputing.
    assert deplete("s", r1, r2, idx, tmp_path / "out", threads=1)[2] == st
