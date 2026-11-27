"""Evidence review for priority pathogens that fail the coverage rule.

The coverage rule answers one question: is there enough contiguous sequence to
call this virus present? For a priority pathogen that Kraken2 named but that
falls short, "not detected" is the right call and stays the call. It does not,
however, say whether the few reads there are look like real pathogen sequence
or like an artefact, and that is what a reviewer has to decide by hand.

This module makes that assessment explicit and reproducible. For every priority
row that failed, it scores seven criteria off the alignment that already
exists, and returns an advisory verdict. It never changes `detection`, and the
verdict is evidence for a human decision, not a positive result:

  sequence_identity     percent identity of the mapped reads to the reference
  mapping_quality       whether the reads map uniquely and reliably
  genome_distribution   spread across the genome vs piled in one spot
  genome_coverage       breadth of the reference covered
  negative_controls     whether the same virus sits in this run's controls
  background            whether a much louder sample could have bled into this one
  reference_specificity whether the reads stay on this virus when related
                        genomes from the same run compete for them
"""

import csv
import os
import re
import shutil
import subprocess
from pathlib import Path

from .utils import PipelineError, capture, log, pipeline, run

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

VERDICTS = ("consistent with pathogen sequence", "equivocal",
            "not supported", "no mapped reads")

REVIEW_COLS = [
    "sample_id", "taxid", "virus_code", "strain_name", "best_reference",
    "kraken_reads", "detection_basis",
    # 1 sequence identity
    "n_reads_reviewed", "pct_identity_mean", "pct_identity_median", "identity_verdict",
    # 2 mapping quality
    "mean_mapq", "pct_reads_passing_mapq", "pct_reads_multimapping", "mapping_verdict",
    # 3 genome distribution
    "n_distinct_start_positions", "n_bins_covered", "n_bins_total",
    "pct_reads_in_busiest_bin", "n_regions_with_reads", "distribution_verdict",
    # 4 genome coverage
    "breadth_of_coverage", "covered_bp", "reference_bp", "longest_block_bp", "coverage_verdict",
    # 5 negative controls
    "control_reads", "control_sample", "control_signal_rpm", "control_verdict",
    # 6 background / index hopping
    "sample_signal_rpm", "max_other_sample_rpm", "max_other_sample", "hopping_ratio",
    "background_verdict",
    # 7 reference specificity
    "n_reads_retained_on_target", "pct_reads_retained_on_target",
    "competing_best_reference", "specificity_verdict",
    # overall
    "criteria_passed", "criteria_assessable", "evidence_verdict", "evidence_reason",
]


def aligned_length(cigar):
    """Reference-consuming + insertion bases, the denominator for identity."""
    return sum(int(n) for n, op in CIGAR.findall(cigar) if op in "MID=X")


def read_alignments(bam):
    """Yield (rname, pos, mapq, cigar, nm, has_alt) for each record in a BAM."""
    proc = subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE, text=True)
    try:
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            nm, alt = None, False
            for tag in f[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                elif tag.startswith(("XA:Z:", "SA:Z:")):
                    alt = True
            yield f[2], int(f[3]), int(f[4]), f[5], nm, alt
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise PipelineError("samtools view failed on %s" % bam)


def reference_lengths(fai):
    out = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    if not n:
        return 0.0
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def alignment_metrics(align_dir, bin_bp):
    """Criteria 1-4, computed from one reference's unique.bam."""
    align_dir = Path(align_dir)
    bam, fai = align_dir / "unique.bam", align_dir / "ref.fasta.fai"
    if not bam.exists() or not fai.exists():
        return None
    lengths = reference_lengths(fai)

    idents, mapqs, starts, bins, per_region = [], [], set(), {}, {}
    n_multi = n = 0
    for rname, pos, mapq, cigar, nm, alt in read_alignments(bam):
        n += 1
        mapqs.append(mapq)
        starts.add((rname, pos))
        per_region[rname] = per_region.get(rname, 0) + 1
        key = (rname, pos // bin_bp)
        bins[key] = bins.get(key, 0) + 1
        if alt:
            n_multi += 1
        alen = aligned_length(cigar)
        if nm is not None and alen:
            idents.append(100.0 * (alen - nm) / alen)

    total_bp = sum(lengths.values())
    n_bins_total = sum(max(1, (L + bin_bp - 1) // bin_bp) for L in lengths.values())
    return {
        "n_reads": n,
        "pct_identity_mean": round(sum(idents) / len(idents), 2) if idents else "",
        "pct_identity_median": round(_median(idents), 2) if idents else "",
        "mean_mapq": round(sum(mapqs) / len(mapqs), 2) if mapqs else "",
        "pct_multimapping": round(100.0 * n_multi / n, 2) if n else "",
        "n_distinct_starts": len(starts),
        "n_bins_covered": len(bins),
        "n_bins_total": n_bins_total,
        "pct_busiest_bin": round(100.0 * max(bins.values()) / n, 2) if bins else "",
        "n_regions_with_reads": len(per_region),
        "reference_bp": total_bp,
    }


def build_competitor_index(sample_dir, target_tag, work_dir, same_taxon_tags=()):
    """One bwa index of every OTHER candidate reference of this sample.

    Related genomes from the same run are the competitors that matter: if the
    reads came from a neighbouring species, or from a conserved stretch that
    several genomes share, they will not stay on this one once the others are
    offered. Built once per (sample, target) and cached.

    References of the SAME species count as on-target. Three adenovirus C
    references trading reads among themselves says nothing about whether the
    sequence is adenovirus; the question is whether it is distinguishable from
    other taxa.
    """
    refs = sorted(p for p in Path(sample_dir, "refs").glob("*.fasta")
                  if p.name[:-len(".fasta")] != target_tag)
    if not refs:
        return None, {}
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    combined = work_dir / "competitors.fasta"
    target_fa = Path(sample_dir, "refs", target_tag + ".fasta")
    on_target = set()
    if not (work_dir / "competitors.fasta.bwt").exists():
        with open(combined, "w") as out:
            for fa in [target_fa] + refs:
                tag = fa.name[:-len(".fasta")]
                with open(fa) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            name = line[1:].split()[0]
                            out.write(">%s|%s\n" % (tag, name))
                        else:
                            out.write(line)
        run(["bwa", "index", combined], stderr=work_dir / "bwa.log")
    for tag in (target_tag,) + tuple(same_taxon_tags):
        fa = Path(sample_dir, "refs", tag + ".fasta")
        if not fa.exists():
            continue
        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    on_target.add("%s|%s" % (tag, line[1:].split()[0]))
    return combined, on_target


def reference_specificity(align_dir, sample_dir, target_tag, work_dir, threads=2,
                          same_taxon_tags=(), tag_labels=None):
    """Criterion 7: do the mapped reads stay on this virus when its neighbours
    from the same run compete for them?"""
    align_dir, work_dir = Path(align_dir), Path(work_dir)
    bam = align_dir / "unique.bam"
    if not bam.exists():
        return None
    index, on_target = build_competitor_index(sample_dir, target_tag, work_dir,
                                              same_taxon_tags=same_taxon_tags)
    if index is None:
        return {"n_retained": "", "pct_retained": "", "competing_best": "no competing reference"}

    work_dir.mkdir(parents=True, exist_ok=True)
    fq = work_dir / ("%s.reads.fastq" % target_tag)
    logf = work_dir / "specificity.log"
    # Everything on stdout, mates included: these reads are realigned as
    # single-end, so pairing does not matter and must not reach the console.
    run(["samtools", "fastq", "-n", bam], stdout=fq, stderr=logf)
    if not os.path.getsize(fq):
        return {"n_retained": 0, "pct_retained": "", "competing_best": ""}

    sam = work_dir / ("%s.competitors.sam" % target_tag)
    run(["bwa", "mem", "-t", str(threads), index, fq], stdout=sam, stderr=logf)
    kept = moved = 0
    elsewhere = {}
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t")
            flag = int(f[1])
            if flag & 0x904:            # secondary, supplementary or unmapped
                continue
            if f[2] in on_target:
                kept += 1
            else:
                moved += 1
                tag = f[2].split("|")[0]
                elsewhere[tag] = elsewhere.get(tag, 0) + 1
    os.remove(sam)
    os.remove(fq)
    total = kept + moved
    best = max(elsewhere, key=elsewhere.get) if elsewhere else ""
    label = (tag_labels or {}).get(best, best)
    return {"n_retained": kept,
            "pct_retained": round(100.0 * kept / total, 2) if total else "",
            "competing_best": ("%s (%d reads)" % (label, elsewhere[best])) if best else ""}


def _f(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def review(summary_tsv, out_path, outdir, params, library_sizes=None, threads=2,
           control_regex=None):
    """Score every failing priority row; write the review table.

    Returns {(sample_id, taxid): {"evidence_verdict": ..., "evidence_reason": ...}}
    so the caller can annotate the main summary.
    """
    control_re = re.compile(control_regex or params.control_sample_regex, re.I)
    lib = {s: int(n) for s, n in (library_sizes or {}).items() if n and int(n) > 0}

    with open(summary_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    def rpm(row):
        n = _f(row.get("total_unique_reads_mapped"))
        total = lib.get(row["sample_id"], 0)
        return (1e6 * n / total) if total else n

    controls = {r["sample_id"] for r in rows if control_re.search(r["sample_id"])}

    # Which candidate references of a sample belong to the same taxon as a
    # given row, so the specificity check does not treat two references of one
    # species as rivals. A candidate table lists every reference considered for
    # that sample, not only the winners in the summary.
    def taxon_key(row):
        return (row.get("species_taxid") or "").strip() or row["taxid"]

    same_taxon, tag_labels = {}, {}
    for sample_id in {r["sample_id"] for r in rows}:
        cand = Path(outdir) / "samples" / sample_id / "candidates.tsv"
        if not cand.exists():
            continue
        groups = {}
        with open(cand) as fh:
            for c in csv.DictReader(fh, delimiter="\t"):
                tag = "%s__%s__%s" % (c["taxid"], c["virus_code"], c["ref_index"])
                tag = "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in tag)
                groups.setdefault(taxon_key(c), []).append(tag)
                tag_labels[tag] = (c.get("species_name") or c.get("strain_name")
                                   or c.get("virus_code") or tag)
        same_taxon[sample_id] = groups

    # Per taxon: the loudest control, and the loudest other sample on the run
    # (a plausible index-hopping source).
    by_taxid = {}
    for r in rows:
        by_taxid.setdefault(r["taxid"], []).append(r)

    todo = [r for r in rows
            if r.get("virus_importance", "priority") == "priority"
            and r.get("detection") != "detected"
            and r["sample_id"] not in controls
            and _f(r.get("kraken_reads")) >= params.evidence_min_kraken_reads]

    log.info("evidence review: %d priority row(s) failed the coverage rule with "
             ">= %d Kraken2 reads", len(todo), params.evidence_min_kraken_reads)

    out_rows, verdicts = [], {}
    work_root = Path(outdir) / "work" / "evidence"
    for i, r in enumerate(todo, 1):
        sample_id, taxid, tag = r["sample_id"], r["taxid"], r.get("best_reference", "")
        sample_dir = Path(outdir) / "samples" / sample_id
        align_dir = sample_dir / "align" / tag
        log.info("evidence review [%d/%d] %s %s", i, len(todo), sample_id,
                 r.get("virus_code") or taxid)

        m = alignment_metrics(align_dir, params.evidence_bin_bp) if tag else None
        spec = None
        if m and m["n_reads"] >= params.evidence_min_reads_for_specificity:
            kin = [t for t in same_taxon.get(sample_id, {}).get(taxon_key(r), []) if t != tag]
            try:
                spec = reference_specificity(align_dir, sample_dir, tag,
                                             work_root / sample_id / tag, threads=threads,
                                             same_taxon_tags=kin, tag_labels=tag_labels)
            except PipelineError as exc:
                log.warning("  specificity check failed: %s", exc)

        # --- criterion verdicts ------------------------------------------
        n_reads = m["n_reads"] if m else 0
        passed, assessable, notes = 0, 0, []

        def score(ok, label, note):
            nonlocal passed, assessable
            if ok is None:
                return "not assessable"
            assessable += 1
            if ok:
                passed += 1
                return "pass"
            notes.append(note)
            return "fail"

        ident = m["pct_identity_mean"] if m else ""
        identity_v = score(None if ident == "" else ident >= params.evidence_min_identity,
                           "identity",
                           "reads sit at %s%% identity, under %g%%"
                           % (ident, params.evidence_min_identity))

        # n_mapped (no MAPQ floor) vs the unique reads that survived it.
        n_mapped = _f(r.get("total_unique_reads_mapped")) + _f(r.get("n_duplicate_reads"))
        pct_pass_mapq = ""
        aln = align_dir / "aln_stats.tsv"
        if aln.exists():
            f = aln.read_text().rstrip("\n").split("\t")
            if len(f) >= 5 and _f(f[4]):
                pct_pass_mapq = round(100.0 * n_reads / _f(f[4]), 2)
        multi = m["pct_multimapping"] if m else ""
        mapping_ok = None
        if m and n_reads:
            mapping_ok = (_f(m["mean_mapq"]) >= params.evidence_min_mean_mapq
                          and _f(multi) <= params.evidence_max_multimapping)
        mapping_v = score(mapping_ok, "mapping",
                          "mean MAPQ %s with %s%% multi-mapping reads"
                          % (m["mean_mapq"] if m else "", multi))

        dist_ok = None
        if m and n_reads:
            dist_ok = (m["n_bins_covered"] >= params.evidence_min_bins
                       and _f(m["pct_busiest_bin"]) <= params.evidence_max_busiest_bin)
        distribution_v = score(dist_ok, "distribution",
                               "reads sit in %s %d-bp window(s), %s%% of them in one"
                               % (m["n_bins_covered"] if m else 0, params.evidence_bin_bp,
                                  m["pct_busiest_bin"] if m else ""))

        breadth = _f(r.get("breadth_of_coverage"))
        coverage_v = score(breadth >= params.evidence_min_breadth, "coverage",
                           "%.3f%% of the genome covered, under %g%%"
                           % (100 * breadth, 100 * params.evidence_min_breadth))

        # Controls and background, over the same taxon elsewhere on the run.
        ctrl_rows = [x for x in by_taxid.get(taxid, []) if x["sample_id"] in controls]
        ctrl = max(ctrl_rows, key=lambda x: rpm(x)) if ctrl_rows else None
        ctrl_reads = _f(ctrl.get("total_unique_reads_mapped")) if ctrl else 0
        this_rpm = rpm(r)
        control_ok = None
        if controls:
            control_ok = (ctrl is None or ctrl_reads < params.min_control_support_reads
                          or rpm(ctrl) <= 0
                          or this_rpm >= params.min_control_fold * rpm(ctrl))
        control_v = score(control_ok, "controls",
                          "control %s carries %d reads of this virus"
                          % (ctrl["sample_id"] if ctrl else "", ctrl_reads))

        others = [x for x in by_taxid.get(taxid, [])
                  if x["sample_id"] != sample_id and x["sample_id"] not in controls]
        loudest = max(others, key=lambda x: rpm(x)) if others else None
        max_other = rpm(loudest) if loudest else 0.0
        ratio = round(max_other / this_rpm, 1) if this_rpm > 0 and max_other else ""
        hop_ok = None
        if others:
            hop_ok = not (max_other > 0 and this_rpm > 0
                          and max_other / this_rpm >= params.evidence_hopping_ratio)
        background_v = score(hop_ok, "background",
                             "%s on this run carries %sx more of it, so cross-talk "
                             "could explain this signal"
                             % (loudest["sample_id"] if loudest else "", ratio))

        spec_ok = None
        if spec and spec["pct_retained"] != "":
            spec_ok = _f(spec["pct_retained"]) >= params.evidence_min_retained
        specificity_v = score(spec_ok, "specificity",
                              "only %s%% of the reads stay on this virus when related "
                              "references compete; %s takes the rest"
                              % (spec["pct_retained"] if spec else "",
                                 spec["competing_best"] if spec else ""))

        # --- overall ------------------------------------------------------
        if n_reads == 0:
            verdict = "no mapped reads"
            reason = "Kraken2 assigned %s read(s) but none aligned to the reference" % r.get("kraken_reads")
        elif assessable and passed == assessable:
            verdict = "consistent with pathogen sequence"
            reason = ("all %d assessable criteria pass on %d read(s); too little sequence "
                      "for a call, worth confirming by a targeted assay" % (assessable, n_reads))
        elif passed >= assessable - 1 and passed >= params.evidence_min_criteria:
            verdict = "equivocal"
            reason = ("%d of %d criteria pass, %s" % (passed, assessable, "; ".join(notes)))
        else:
            verdict = "not supported"
            reason = "%d of %d criteria pass; %s" % (passed, assessable, "; ".join(notes))

        verdicts[(sample_id, taxid)] = {"evidence_verdict": verdict, "evidence_reason": reason}
        out_rows.append({
            "sample_id": sample_id, "taxid": taxid, "virus_code": r.get("virus_code", ""),
            "strain_name": r.get("strain_name", ""), "best_reference": tag,
            "kraken_reads": r.get("kraken_reads", ""),
            "detection_basis": r.get("detection_basis", ""),
            "n_reads_reviewed": n_reads,
            "pct_identity_mean": ident,
            "pct_identity_median": m["pct_identity_median"] if m else "",
            "identity_verdict": identity_v,
            "mean_mapq": m["mean_mapq"] if m else "",
            "pct_reads_passing_mapq": pct_pass_mapq,
            "pct_reads_multimapping": multi,
            "mapping_verdict": mapping_v,
            "n_distinct_start_positions": m["n_distinct_starts"] if m else 0,
            "n_bins_covered": m["n_bins_covered"] if m else 0,
            "n_bins_total": m["n_bins_total"] if m else "",
            "pct_reads_in_busiest_bin": m["pct_busiest_bin"] if m else "",
            "n_regions_with_reads": m["n_regions_with_reads"] if m else 0,
            "distribution_verdict": distribution_v,
            "breadth_of_coverage": r.get("breadth_of_coverage", ""),
            "covered_bp": round(breadth * (m["reference_bp"] if m else 0)) if m else "",
            "reference_bp": m["reference_bp"] if m else "",
            "longest_block_bp": r.get("longest_covered_block_bp", ""),
            "coverage_verdict": coverage_v,
            "control_reads": int(ctrl_reads) if ctrl else "",
            "control_sample": ctrl["sample_id"] if ctrl else "",
            "control_signal_rpm": round(rpm(ctrl), 4) if ctrl else "",
            "control_verdict": control_v,
            "sample_signal_rpm": round(this_rpm, 4),
            "max_other_sample_rpm": round(max_other, 4) if loudest else "",
            "max_other_sample": loudest["sample_id"] if loudest else "",
            "hopping_ratio": ratio,
            "background_verdict": background_v,
            "n_reads_retained_on_target": spec["n_retained"] if spec else "",
            "pct_reads_retained_on_target": spec["pct_retained"] if spec else "",
            "competing_best_reference": spec["competing_best"] if spec else "",
            "specificity_verdict": specificity_v,
            "criteria_passed": passed, "criteria_assessable": assessable,
            "evidence_verdict": verdict, "evidence_reason": reason,
        })

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=REVIEW_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)

    shutil.rmtree(work_root, ignore_errors=True)
    by_verdict = {}
    for r in out_rows:
        by_verdict[r["evidence_verdict"]] = by_verdict.get(r["evidence_verdict"], 0) + 1
    for v in VERDICTS:
        if by_verdict.get(v):
            log.info("  %-35s %d", v, by_verdict[v])
    return verdicts


def annotate(summary_tsv, verdicts, out_path):
    """Add evidence_verdict / evidence_reason to the final summary."""
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        cols = list(reader.fieldnames or [])
    new = ["evidence_verdict", "evidence_reason"]
    for r in rows:
        v = verdicts.get((r["sample_id"], r["taxid"]))
        r["evidence_verdict"] = v["evidence_verdict"] if v else "not reviewed"
        r["evidence_reason"] = v["evidence_reason"] if v else ""
    cols = [c for c in cols if c not in new] + new
    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return out_path
