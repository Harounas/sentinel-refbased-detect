"""Optional host depletion (the step that produced the cleaned BAMs on Terra).

Reads are aligned to a host genome with minimap2 (-x sr). A pair is kept only
when neither mate aligns to the host; a single-end read when it does not
align. The kept reads are written to OUTDIR/host_depleted/ and everything
downstream (Kraken2, alignment, library sizes) runs on them.
"""

import csv
import gzip
from pathlib import Path

from .utils import PipelineError, is_done, log, pipeline, run, write_done

MMI_SUFFIX = ".mmi"


def host_index(host_ref, cache_dir, threads=4):
    """Return a minimap2 index for host_ref, building it once if given a FASTA."""
    host_ref = Path(host_ref).expanduser()
    if not host_ref.exists():
        raise PipelineError("host reference not found: %s" % host_ref)
    if host_ref.name.endswith(MMI_SUFFIX):
        return host_ref.resolve()
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    stem = host_ref.name
    for ext in (".gz", ".fasta", ".fa", ".fna"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    mmi = cache_dir / (stem + ".sr" + MMI_SUFFIX)
    if not mmi.exists():
        log.info("building minimap2 index of %s (once; a human genome takes ~10 min "
                 "and ~16 GB RAM)", host_ref)
        tmp = mmi.with_name(mmi.name + ".part")
        run(["minimap2", "-x", "sr", "-t", threads, "-d", tmp, host_ref],
            stderr=cache_dir / "minimap2_index.log")
        tmp.rename(mmi)
    return mmi


def count_reads(path):
    if not path:
        return 0
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for n, _ in enumerate(fh, 1):
            pass
    return n // 4


def deplete(sample_id, r1, r2, index, out_dir, threads=4, resume=True):
    """Return (r1, r2, stats) for the host-depleted reads of one sample."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sdir = out_dir / sample_id
    o1 = out_dir / ("%s_R1.fastq.gz" % sample_id)
    o2 = out_dir / ("%s_R2.fastq.gz" % sample_id) if r2 else None
    stats_path = sdir / "stats.tsv"

    if not (resume and is_done(sdir) and o1.exists()):
        sdir.mkdir(parents=True, exist_ok=True)
        logf = sdir / "host_depletion.log"
        t = str(threads)
        mm = ["minimap2", "-a", "-x", "sr", "-t", t, index, r1] + ([r2] if r2 else [])
        if r2:
            # -f 12: read AND mate unmapped, so a pair with either mate on the
            # host is removed. -F 0x900: no secondary/supplementary records.
            fq = ["samtools", "fastq", "-f", "12", "-F", "0x900", "-c", "6",
                  "-1", o1, "-2", o2, "-0", "/dev/null", "-s", "/dev/null", "-n", "-"]
        else:
            fq = ["samtools", "fastq", "-f", "4", "-F", "0x900", "-c", "6", "-0", o1, "-"]
        log.info("host depletion: %s", sample_id)
        pipeline([mm, fq], stdout=logf, stderr=logf)

        before = count_reads(r1) + count_reads(r2)
        after = count_reads(o1) + count_reads(o2)
        with open(stats_path, "w") as fh:
            fh.write("%d\t%d\n" % (before, after))
        write_done(sdir)

    before, after = (int(x) for x in stats_path.read_text().split())
    host = before - after
    stats = {"sample_id": sample_id, "reads_in": before, "reads_kept": after,
             "reads_host": host,
             "pct_host": "%.2f" % (100.0 * host / before) if before else "0.00"}
    log.info("host depletion: %s kept %d of %d reads (%s%% host)",
             sample_id, after, before, stats["pct_host"])
    return o1, o2, stats


def write_summary(rows, path):
    cols = ["sample_id", "reads_in", "reads_kept", "reads_host", "pct_host"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
