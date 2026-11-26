"""Task kraken2_classify: one Kraken2 report per sample."""

import os
import shutil
import tarfile
from pathlib import Path

from .samples import has_reads
from .utils import PipelineError, log, run

K2D = ("hash.k2d", "opts.k2d", "taxo.k2d")


def prepare_db(db, work_dir):
    """Accept a directory holding hash/opts/taxo.k2d, or a .tar.gz/.tar.zst
    archive of one (unpacked once into work_dir)."""
    db = Path(db).expanduser()
    if db.is_dir():
        if all((db / f).exists() for f in K2D):
            return db.resolve()
        found = next(db.rglob("hash.k2d"), None)
        if found and all((found.parent / f).exists() for f in K2D):
            return found.parent.resolve()
        raise PipelineError("%s does not contain %s" % (db, ", ".join(K2D)))
    if not db.exists():
        raise PipelineError("Kraken2 database not found: %s" % db)

    out = Path(work_dir) / "kraken2_db"
    if not all((out / f).exists() for f in K2D):
        out.mkdir(parents=True, exist_ok=True)
        log.info("unpacking Kraken2 database %s -> %s", db, out)
        if db.name.endswith(".zst"):
            if not shutil.which("zstd"):
                raise PipelineError("zstd is needed to unpack %s" % db)
            run(["sh", "-c", 'zstd -d --stdout "$0" | tar -x -C "$1"', db, out])
        else:
            with tarfile.open(db) as tf:
                tf.extractall(out)
        if not (out / "hash.k2d").exists():
            inner = next(out.rglob("hash.k2d"), None)
            if inner is None:
                raise PipelineError("no hash.k2d inside %s" % db)
            for f in K2D:
                shutil.move(str(inner.parent / f), str(out / f))
    for f in K2D:
        if not (out / f).stat().st_size:
            raise PipelineError("kraken2 db file %s is empty" % f)
    return out


def classify(sample, r1, r2, db_dir, out_dir, threads=8, memory_mapping=False,
             keep_output=False):
    """Run Kraken2 and return the report path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / ("%s.kraken2.report.txt" % sample.sample_id)
    output = (out_dir / ("%s.kraken2.output.txt" % sample.sample_id)) if keep_output else Path(os.devnull)
    logf = out_dir / ("%s.kraken2.log" % sample.sample_id)
    tmp_report = report.with_name(report.name + ".part")

    reads = [f for f in (r1, r2) if f]
    if not any(has_reads(f) for f in reads):
        # An empty barcode gives an empty report rather than failing the run.
        log.warning("%s has no reads; emitting an empty report", sample.sample_id)
        report.write_text("")
        return report

    cmd = ["kraken2", "--db", db_dir, "--threads", threads,
           "--report", tmp_report, "--output", output]
    if memory_mapping:
        cmd.append("--memory-mapping")
    if any(str(f).endswith(".gz") for f in reads):
        cmd.append("--gzip-compressed")
    if r1 and r2 and has_reads(r1) and has_reads(r2):
        cmd += ["--paired", r1, r2]
    else:
        cmd += [f for f in reads if has_reads(f)]
    log.info("kraken2: %s", sample.sample_id)
    run(cmd, stderr=logf)
    if not tmp_report.exists():
        tmp_report.write_text("")
    tmp_report.rename(report)
    if keep_output and output.exists():
        run(["gzip", "-f", output])
    return report


def library_size(report_path, paired):
    """Reads in the library, from a Kraken2 report: unclassified + root clade
    counts give fragments; x2 for paired input so the number matches the
    Terra read_counts_depleted column (reads in the cleaned BAM)."""
    total = 0
    try:
        with open(report_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                taxid = f[6].strip() if len(f) >= 8 and f[6].strip().isdigit() else f[4].strip()
                if taxid in ("0", "1"):
                    try:
                        total += int(f[1])
                    except ValueError:
                        pass
    except FileNotFoundError:
        return 0
    return total * (2 if paired else 1)
