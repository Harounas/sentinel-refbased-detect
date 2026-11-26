"""Task align_and_consensus: align one sample to one candidate reference.

Writes into its own directory:
  region_stats.tsv   tag, region, length, unique_reads, mean_mapq,
                     covered_bp(>=min_depth), longest_block_bp, n_blocks(>=min_block_bp)
  aln_stats.tsv      tag, taxid, code, ref_index, n_mapped, n_proper, n_dup, consensus_status
  <tag>.consensus.fasta, <tag>.block.fasta, unique.bam, flagstat.txt
"""

import os
import shutil
import subprocess
from pathlib import Path

from .utils import capture, file_nonempty, is_done, log, pipeline, run, write_done

INTERMEDIATES = ("namesort.bam", "fixmate.bam", "positionsort.bam", "marked.bam",
                 "marked.bam.bai", "region.bam", "region.bam.bai")


def parse_tag(tag):
    """<taxid>__<code>__<index>; the code itself may contain '__'."""
    taxid, rest = tag.split("__", 1)
    code, idx = rest.rsplit("__", 1)
    return taxid, code, idx


def coverage_blocks(bam, min_depth, min_block_bp):
    """Collapse `samtools depth -a` into per-region (covered, longest, n_blocks).

    Only regions with at least one base at min_depth are returned, as in the
    WDL's awk."""
    covered, longest, nblocks = {}, {}, {}
    cur, run_len = None, 0

    def close_run():
        nonlocal run_len
        if run_len > 0:
            if run_len > longest.get(cur, 0):
                longest[cur] = run_len
            if run_len >= min_block_bp:
                nblocks[cur] = nblocks.get(cur, 0) + 1
        run_len = 0

    proc = subprocess.Popen(["samtools", "depth", "-a", "-Q", "0", str(bam)],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t")
        if len(f) < 3:
            continue
        if f[0] != cur:
            close_run()
            cur = f[0]
        if int(f[2]) >= min_depth:
            run_len += 1
            covered[cur] = covered.get(cur, 0) + 1
        else:
            close_run()
    close_run()
    if proc.wait() != 0:
        raise RuntimeError("samtools depth failed on %s" % bam)
    return {r: (covered[r], longest.get(r, 0), nblocks.get(r, 0)) for r in covered}


def mapq_by_ref(bam):
    n, s = {}, {}
    proc = subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.split("\t", 5)
        n[f[2]] = n.get(f[2], 0) + 1
        s[f[2]] = s.get(f[2], 0) + int(f[4])
    if proc.wait() != 0:
        raise RuntimeError("samtools view failed on %s" % bam)
    return {r: (n[r], float("%.2f" % (s[r] / n[r]))) for r in n}


def longest_called_run(fasta_path):
    """Longest unbroken run of called bases (not N / -) in an ivar consensus.
    Returns (start_1based, length, sequence)."""
    seq = []
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    seq = "".join(seq)
    best = bstart = start = run_len = 0
    for i, c in enumerate(seq, 1):
        if c in "N-":
            if run_len > best:
                best, bstart = run_len, start
            run_len = 0
        else:
            if run_len == 0:
                start = i
            run_len += 1
    if run_len > best:
        best, bstart = run_len, start
    return bstart, best, seq[bstart - 1:bstart - 1 + best] if best else ""


def align_and_consensus(sample_id, r1, r2, reference_fasta, out_dir, params,
                        threads=4, keep_bam=True, keep_intermediate=False, resume=True):
    tag = Path(reference_fasta).name[:-len(".fasta")]
    wd = Path(out_dir)
    if resume and is_done(wd):
        return wd
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)
    logf = wd / "align.log"
    t = str(threads)

    shutil.copyfile(reference_fasta, wd / "ref.fasta")
    run(["bwa", "index", "ref.fasta"], cwd=wd, stderr=logf)
    run(["samtools", "faidx", "ref.fasta"], cwd=wd, stderr=logf)

    reads = [str(Path(r1).resolve())]
    if r2:
        reads.append(str(Path(r2).resolve()))
    # bwa streams straight into the name sort, so no SAM file hits the disk.
    pipeline([["bwa", "mem", "-t", t, "ref.fasta", *reads],
              ["samtools", "sort", "-n", "-@", t, "-o", "namesort.bam", "-"]],
             cwd=wd, stderr=logf)
    run(["samtools", "fixmate", "-m", "-@", t, "namesort.bam", "fixmate.bam"], cwd=wd, stderr=logf)
    run(["samtools", "sort", "-@", t, "-o", "positionsort.bam", "fixmate.bam"], cwd=wd, stderr=logf)
    run(["samtools", "markdup", "-@", t, "positionsort.bam", "marked.bam"], cwd=wd, stderr=logf)
    run(["samtools", "index", "marked.bam"], cwd=wd, stderr=logf)

    # Unique = mapped, primary, non-supplementary, non-duplicate, MAPQ >= floor.
    run(["samtools", "view", "-b", "-F", "0xD04", "-q", str(params.min_mapq),
         "-o", "unique.bam", "marked.bam"], cwd=wd, stderr=logf)
    run(["samtools", "index", "unique.bam"], cwd=wd, stderr=logf)
    run(["samtools", "flagstat", "marked.bam"], cwd=wd, stdout=wd / "flagstat.txt", stderr=logf)

    ref_lengths = []
    with open(wd / "ref.fasta.fai") as fh:
        for line in fh:
            f = line.split("\t")
            ref_lengths.append((f[0], int(f[1])))

    mq = mapq_by_ref(wd / "unique.bam")
    blocks = coverage_blocks(wd / "unique.bam", params.min_depth, params.min_block_bp)

    def count(*flags):
        return int(capture(["samtools", "view", "-c", *flags, "marked.bam"], cwd=wd).strip())

    n_mapped = count("-F", "0xD04")
    n_proper = count("-F", "0xD04", "-f", "0x2")
    n_dup = count("-F", "0x904", "-f", "0x400")

    # ivar consensus can segfault on some pileups; the consensus is secondary,
    # so try unlimited depth, then a capped depth, then record the failure.
    def run_consensus(depth, prefix, bam="unique.bam", region=None):
        mp = ["samtools", "mpileup", "-aa", "-A", "-d", str(depth), "-Q", "0"]
        if region:
            mp += ["-r", region]
        mp += ["--reference", "ref.fasta", bam]
        iv = ["ivar", "consensus", "-p", prefix, "-q", str(params.consensus_min_q),
              "-t", str(params.consensus_min_freq), "-m", str(params.consensus_min_depth),
              "-n", params.consensus_ambig_char]
        with open(wd / "ivar.log", "a") as lf:
            rc = pipeline([mp, iv], cwd=wd, stdout=lf, stderr=lf, check=False)
        # Report a signal the way bash does (segfault -11 -> 139), as on Terra.
        return 128 - rc if rc < 0 else rc

    status = "ok"
    rc = run_consensus(0, "ivar_out")
    if rc != 0 or not file_nonempty(wd / "ivar_out.fa"):
        log.warning("%s %s: ivar failed at unlimited depth (rc=%s); retrying at -d %d",
                    sample_id, tag, rc, params.consensus_retry_depth)
        if (wd / "ivar_out.fa").exists():
            os.remove(wd / "ivar_out.fa")
        rc = run_consensus(params.consensus_retry_depth, "ivar_out")
        status = ("failed_rc%s" % rc) if (rc != 0 or not file_nonempty(wd / "ivar_out.fa")) \
            else "ok_depth_capped"

    cons = wd / (tag + ".consensus.fasta")
    if file_nonempty(wd / "ivar_out.fa"):
        shutil.move(str(wd / "ivar_out.fa"), str(cons))
    else:
        cons.write_text("")

    # The longest covered block, as the sample's own sequence, for DIAMOND.
    block_fa = wd / (tag + ".block.fasta")
    block_fa.write_text("")
    if blocks:
        best_region = max(blocks, key=lambda r: (blocks[r][1], r))
        best_len = blocks[best_region][1]
        if best_len >= params.min_block_bp:
            ok = run(["samtools", "view", "-b", "-o", "region.bam", "unique.bam", best_region],
                     cwd=wd, stderr=logf, check=False) == 0
            ok = ok and run(["samtools", "index", "region.bam"], cwd=wd, stderr=logf, check=False) == 0
            brc = run_consensus(0, "block_out", bam="region.bam", region=best_region) if ok else 1
            if brc == 0 and file_nonempty(wd / "block_out.fa"):
                start, length, seq = longest_called_run(wd / "block_out.fa")
                if length >= params.min_block_bp:
                    block_fa.write_text(">%s|%s|%d-%d|%d\n%s\n" % (
                        tag, best_region, start, start + length - 1, length, seq))
            else:
                log.warning("%s %s: block consensus failed for %s", sample_id, tag, best_region)

    with open(wd / "region_stats.tsv", "w") as out:
        for r, length in ref_lengths:
            ur, mean_mq = mq.get(r, (0, 0.0))
            cov, lng, nb = blocks.get(r, (0, 0, 0))
            out.write("%s\t%s\t%d\t%d\t%.2f\t%d\t%d\t%d\n" % (tag, r, length, ur, mean_mq, cov, lng, nb))

    taxid, code, idx = parse_tag(tag)
    with open(wd / "aln_stats.tsv", "w") as out:
        out.write("\t".join([tag, taxid, code, idx, str(n_mapped), str(n_proper),
                             str(n_dup), status]) + "\n")

    if not keep_intermediate:
        for f in INTERMEDIATES:
            p = wd / f
            if p.exists():
                os.remove(p)
        for p in wd.glob("ref.fasta.*"):
            if not p.name.endswith(".fai"):
                os.remove(p)
        for p in wd.glob("*_out.qual.txt"):
            os.remove(p)
    if not keep_bam:
        for f in ("unique.bam", "unique.bam.bai"):
            if (wd / f).exists():
                os.remove(wd / f)
    write_done(wd)
    return wd
