"""Find samples in a FASTQ directory (or read them from a sample sheet).

Sample ids follow the Terra convention: the leading field of the file name up
to the first '.', after Illumina suffixes (_S1, _L001, _R1, _001) are removed.
E.g. 1-02-26-0109_S5_L001_R1_001.fastq.gz -> 1-02-26-0109. This matters
because negative controls are recognised by their id (NEC..., NTC..., ...).
"""

import csv
import gzip
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from .utils import PipelineError

FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
MATE = re.compile(
    r"^(?P<stem>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?[_.](?P<mate>R?[12])(?:_\d{3})?$",
    re.I)
SINGLE = re.compile(r"^(?P<stem>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?(?:_\d{3})?$", re.I)


@dataclass
class Sample:
    sample_id: str
    r1: List[Path] = field(default_factory=list)
    r2: List[Path] = field(default_factory=list)

    @property
    def paired(self):
        return bool(self.r1) and bool(self.r2)

    @property
    def files(self):
        return list(self.r1) + list(self.r2)


def sample_id_from_name(filename):
    """Return (sample_id, mate) for a FASTQ filename; mate is 1, 2 or None."""
    base = FASTQ_EXT.sub("", Path(filename).name)
    m = MATE.match(base)
    if m:
        mate = 2 if m.group("mate")[-1] == "2" else 1
        stem = m.group("stem")
    else:
        mate = None
        stem = SINGLE.match(base).group("stem")
    return stem.split(".")[0], mate


def discover(fastq_dir, recursive=False):
    d = Path(fastq_dir).expanduser()
    if not d.is_dir():
        raise PipelineError("FASTQ directory not found: %s" % fastq_dir)
    it = d.rglob("*") if recursive else d.iterdir()
    files = sorted(p for p in it if p.is_file() and FASTQ_EXT.search(p.name)
                   and not p.name.startswith("."))
    if not files:
        raise PipelineError("no *.fastq[.gz] / *.fq[.gz] files in %s" % d)

    samples = {}
    singles = {}
    for f in files:
        sid, mate = sample_id_from_name(f.name)
        s = samples.setdefault(sid, Sample(sid))
        if mate == 2:
            s.r2.append(f.resolve())
        elif mate == 1:
            s.r1.append(f.resolve())
        else:
            singles.setdefault(sid, []).append(f.resolve())

    for sid, fs in singles.items():
        s = samples[sid]
        if s.r1 or s.r2:
            raise PipelineError(
                "sample %s mixes paired files with unpaired file(s) %s; "
                "use --samplesheet to say which is which" % (sid, [x.name for x in fs]))
        s.r1.extend(fs)

    for s in samples.values():
        if s.r2 and len(s.r1) != len(s.r2):
            raise PipelineError(
                "sample %s: %d R1 file(s) but %d R2 file(s): %s"
                % (s.sample_id, len(s.r1), len(s.r2), [p.name for p in s.files]))
        if s.r2 and not s.r1:
            raise PipelineError("sample %s has R2 but no R1" % s.sample_id)
    return [samples[k] for k in sorted(samples)]


def from_samplesheet(path):
    """TSV/CSV with columns sample_id, fastq_1[, fastq_2]. Repeated sample ids
    (e.g. one row per lane) are merged."""
    path = Path(path)
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    samples = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter=delim):
            sid = (row.get("sample_id") or row.get("sample") or "").strip()
            r1 = (row.get("fastq_1") or row.get("r1") or "").strip()
            r2 = (row.get("fastq_2") or row.get("r2") or "").strip()
            if not sid or not r1:
                continue
            s = samples.setdefault(sid, Sample(sid))
            p1 = (path.parent / r1).resolve() if not Path(r1).is_absolute() else Path(r1)
            s.r1.append(p1)
            if r2:
                p2 = (path.parent / r2).resolve() if not Path(r2).is_absolute() else Path(r2)
                s.r2.append(p2)
    for s in samples.values():
        for f in s.files:
            if not f.exists():
                raise PipelineError("sample %s: file not found: %s" % (s.sample_id, f))
        if s.r2 and len(s.r1) != len(s.r2):
            raise PipelineError("sample %s: unequal R1/R2 file counts" % s.sample_id)
    if not samples:
        raise PipelineError("no samples in %s (need columns sample_id, fastq_1[, fastq_2])" % path)
    return [samples[k] for k in sorted(samples)]


def has_reads(path):
    path = str(path)
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rb") as fh:
            return bool(fh.read(1))
    except (OSError, EOFError):
        return False


def merge_lanes(sample, work_dir):
    """bwa mem takes one file per mate, so concatenate multi-lane input once.
    Concatenated gzip members are themselves a valid gzip file."""
    if len(sample.r1) <= 1:
        return sample.r1[0] if sample.r1 else None, (sample.r2[0] if sample.r2 else None)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    out = []
    for mate, files in (("R1", sample.r1), ("R2", sample.r2)):
        if not files:
            out.append(None)
            continue
        gz = all(str(f).endswith(".gz") for f in files)
        dest = work_dir / ("%s_%s.fastq%s" % (sample.sample_id, mate, ".gz" if gz else ""))
        if not dest.exists():
            tmp = dest.with_name(dest.name + ".part")
            with open(tmp, "wb") as o:
                for f in files:
                    opener = gzip.open if (not gz and str(f).endswith(".gz")) else open
                    with opener(f, "rb") as i:
                        shutil.copyfileobj(i, o, 1 << 20)
            tmp.rename(dest)
        out.append(dest)
    return out[0], out[1]
