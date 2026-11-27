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
    "n_reads_reviewed", "pct_identity_mean", "pct_identity_median",
    "pct_identity_closest_genome", "identity_verdict",
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
    "pct_reads_attributed_to_target", "n_reads_unique_to_target",
    "n_reads_unique_elsewhere", "n_reads_ambiguous",
    "n_competing_references", "competing_references_tested",
    "reads_by_taxon", "reads_by_reference", "competing_best_reference",
    "specificity_verdict",
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


def build_competitor_index(sample_dir, target_tag, work_dir, same_taxon_tags=(), decoys=()):
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
    decoys = [Path(d) for d in decoys]
    # The index always holds the target, even when nothing else is on offer:
    # the background layer supplies the competitors in that case, and the
    # reads still have to be scored against their own reference.
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    combined = work_dir / "competitors.fasta"
    target_fa = Path(sample_dir, "refs", target_tag + ".fasta")
    on_target = set()
    if not (work_dir / "competitors.fasta.bwt").exists():
        with open(combined, "w") as out:
            for fa in [target_fa] + refs + decoys:
                tag = fa.name[:-len(".fasta")]
                with open(fa) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            name = line[1:].split()[0]
                            out.write(">%s|%s\n" % (tag, name))
                        else:
                            out.write(line)
        run(["bwa", "index", combined], stderr=work_dir / "bwa.log")
    decoy_by_tag = {p.name[:-len(".fasta")]: p for p in decoys}
    for tag in (target_tag,) + tuple(same_taxon_tags):
        fa = decoy_by_tag.get(tag) or Path(sample_dir, "refs", tag + ".fasta")
        if not fa.exists():
            continue
        with open(fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    on_target.add("%s|%s" % (tag, line[1:].split()[0]))
    competitors = [p.name[:-len(".fasta")] for p in refs + decoys
                   if p.name[:-len(".fasta")] not in set(same_taxon_tags)]
    return combined, on_target, competitors


def align_to_background(fastq, background, work_dir, threads=2):
    """Best background hit per read: (taxon label, alignment score).

    One shared index of every RefSeq viral species, the host and the usual
    reagent organisms, built once rather than per pathogen. It is kept separate
    from the small per-pathogen index because 3.6 Gb cannot be rebuilt for
    every row under review.
    """
    index = Path(background["index"])
    if not Path(str(index) + ".bwt").exists():
        log.warning("background index %s is not built; skipping the background layer", index)
        return {}
    labels = background.get("labels") or {}
    work_dir = Path(work_dir)
    sam = work_dir / "background.sam"
    run(["bwa", "mem", "-t", str(threads), str(index), str(fastq)],
        stdout=sam, stderr=work_dir / "background.log")
    best = {}
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or int(f[1]) & 0x804:
                continue
            score = None
            for t in f[11:]:
                if t.startswith("AS:i:"):
                    score = int(t[5:])
                    break
            if score is None:
                continue
            label = labels.get(f[2], f[2])
            cur = best.get(f[0])
            if cur is None or score > cur[1]:
                best[f[0]] = (label, score)
    os.remove(sam)
    return best


def load_background(path):
    """A prebuilt background index plus its sequence-to-taxon labels.

    `path` is the background FASTA that bwa indexed; background.taxa.tsv beside
    it names each sequence. Without the labels every background hit would be
    reported as a bare accession.
    """
    index = Path(path)
    labels = {}
    taxa = index.with_name("background.taxa.tsv")
    if not taxa.exists():
        taxa = index.with_suffix(".taxa.tsv")
    if taxa.exists():
        with open(taxa) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    labels[f[0]] = f[2] or f[1] or f[0]
    log.info("background layer: %s (%d labelled sequences)", index, len(labels))
    return {"index": index, "labels": labels}


def reference_specificity(align_dir, sample_dir, target_tag, work_dir, threads=2,
                          same_taxon_tags=(), tag_labels=None, decoys=(),
                          tag_taxon=None, margin=0, em_rounds=20, background=None):
    """Criterion 7: can these reads be attributed to this taxon?

    Not "did they stay put": with sibling species and strain genomes in the
    pool a true read often fits several references equally well, and counting
    that as a loss would fail real findings. Reads that fit one taxon alone are
    counted there; reads that tie are shared out in proportion to that
    unambiguous evidence, which is what the user's weighting asks for.
    """
    align_dir, work_dir = Path(align_dir), Path(work_dir)
    bam = align_dir / "unique.bam"
    if not bam.exists():
        return None
    labels = tag_labels or {}
    blank = {"pct_target": "", "n_unique_target": "", "n_unique_other": "",
             "n_ambiguous": "", "competing_best": "", "n_competitors": 0,
             "competitors": "", "reads_by_taxon": "", "reads_by_reference": ""}
    built = build_competitor_index(sample_dir, target_tag, work_dir,
                                   same_taxon_tags=same_taxon_tags, decoys=decoys)
    index, on_target, competitors = built

    def described(tags, cap=6):
        """Name the references, by taxon rather than by file tag."""
        names, seen = [], set()
        for t in tags:
            nm = labels.get(t, t)
            if nm not in seen:
                seen.add(nm)
                names.append(nm)
        if len(names) > cap:
            return "; ".join(names[:cap]) + "; +%d more" % (len(names) - cap)
        return "; ".join(names)

    comp_desc = described(competitors)
    work_dir.mkdir(parents=True, exist_ok=True)
    fq = work_dir / ("%s.reads.fastq" % target_tag)
    logf = work_dir / "specificity.log"
    # Everything on stdout, mates included: these reads are realigned as
    # single-end, so pairing does not matter and must not reach the console.
    # -N, not -n: both mates otherwise carry the same name. bwa then strips the
    # /1 and /2 it adds, so rewrite them to a suffix it keeps; without this the
    # two mates merge into one observation and every count below is halved.
    raw = work_dir / ("%s.raw.fastq" % target_tag)
    run(["samtools", "fastq", "-N", bam], stdout=raw, stderr=logf)
    with open(raw) as src, open(fq, "w") as dst:
        for i, line in enumerate(src):
            if i % 4 == 0 and line.startswith("@"):
                name = line[1:].rstrip("\n")
                dst.write("@%s\n" % re.sub(r"/([12])$", r"_mate\1", name))
            else:
                dst.write(line)
    os.remove(raw)
    if not os.path.getsize(fq):
        return dict(blank, n_unique_target=0, n_competitors=len(competitors),
                    competitors=comp_desc)

    # -a reports every alignment, not just the best one, which is what makes
    # the ambiguity visible: a read that fits two strains equally well has to
    # be seen fitting both before it can be apportioned between them.
    sam = work_dir / ("%s.competitors.sam" % target_tag)
    run(["bwa", "mem", "-a", "-t", str(threads), index, fq], stdout=sam, stderr=logf)

    taxon_of = dict(tag_taxon or {})
    TARGET = "__target__"
    for t in (target_tag,) + tuple(same_taxon_tags):
        taxon_of[t] = TARGET

    per_read, landed, best_ident = {}, {}, {}
    background_hits = (align_to_background(fq, background, work_dir, threads=threads)
                       if background else {})
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11 or int(f[1]) & 0x804:     # unmapped or supplementary
                continue
            score, nm = None, None
            for t in f[11:]:
                if t.startswith("AS:i:"):
                    score = int(t[5:])
                elif t.startswith("NM:i:"):
                    nm = int(t[5:])
            if score is None:
                continue
            tag = f[2].split("|")[0]
            hits = per_read.setdefault(f[0], {})
            # Best score per reference, then per taxon below.
            if tag not in hits or score > hits[tag]:
                hits[tag] = score
            # Identity to the closest genome of this taxon. A divergent
            # lineage read against a distant reference looks like a poor match
            # on that reference while being a perfect match to its own clade,
            # so the question worth asking is how well it fits the nearest
            # genome the taxon has, not the one the pipeline happened to pick.
            if nm is not None and taxon_of.get(tag, tag) == TARGET:
                alen = aligned_length(f[5])
                if alen:
                    ident = 100.0 * (alen - nm) / alen
                    best_ident[f[0]] = max(best_ident.get(f[0], 0.0), ident)

    # Whatever the background offered, as one more competitor per read. A read
    # that fits an unlisted virus, the host or a reagent organism better than
    # anything on the priority table now has somewhere true to go.
    for name, (label, score) in background_hits.items():
        hits = per_read.setdefault(name, {})
        key = "bg:" + label
        if key not in hits or score > hits[key]:
            hits[key] = score
        taxon_of[key] = label

    for name, hits in per_read.items():
        best_tag = max(hits, key=hits.get)
        landed[best_tag] = landed.get(best_tag, 0) + 1

    # A read counts for a taxon when that taxon holds an alignment within
    # `margin` of the read's best score anywhere. Ties stay ties instead of
    # being handed to whichever reference bwa happened to report first.
    unique, ambiguous = {}, []
    for name, hits in per_read.items():
        best = max(hits.values())
        taxa = {taxon_of.get(t, t) for t, s in hits.items() if s >= best - margin}
        if len(taxa) == 1:
            k = next(iter(taxa))
            unique[k] = unique.get(k, 0) + 1
        else:
            ambiguous.append(taxa)

    # Reads that fit several taxa equally are split between them in proportion
    # to the reads that fit only one, so a strain carrying real unique support
    # claims the share of the ambiguous reads that its own evidence earns.
    weight = {k: float(v) for k, v in unique.items()}
    for taxa in ambiguous:
        for k in taxa:
            weight.setdefault(k, 0.0)
    if ambiguous and weight:
        for _ in range(em_rounds):
            new = {k: float(unique.get(k, 0)) for k in weight}
            for taxa in ambiguous:
                total_w = sum(weight[k] for k in taxa)
                if total_w > 0:
                    for k in taxa:
                        new[k] += weight[k] / total_w
                else:                       # no unique evidence anywhere: split evenly
                    for k in taxa:
                        new[k] += 1.0 / len(taxa)
            weight = new

    total_w = sum(weight.values())
    pct_target = round(100.0 * weight.get(TARGET, 0.0) / total_w, 2) if total_w else ""
    others = {k: v for k, v in weight.items() if k != TARGET and v > 0}
    best_taxon = max(others, key=others.get) if others else ""
    by_taxon = ";".join(
        "%s:%.1f" % ("this taxon" if k == TARGET else k, v)
        for k, v in sorted(weight.items(), key=lambda kv: (kv[0] != TARGET, -kv[1], kv[0]))
        if v > 0)
    elsewhere = {t: n for t, n in landed.items() if taxon_of.get(t, t) != TARGET}
    best = max(elsewhere, key=elsewhere.get) if elsewhere else ""
    # The background is a real competitor: every RefSeq viral species, the host
    # and the reagent organisms were offered, so a taxon with no siblings of
    # its own is no longer tested against an empty field.
    bg_taxa = sorted({lab for lab, _ in background_hits.values()})
    n_competitors = len(competitors) + len(bg_taxa)
    if bg_taxa:
        comp_desc = ((comp_desc + "; ") if comp_desc else "") + \
                    "background: " + "; ".join(bg_taxa[:4]) + \
                    ("; +%d more" % (len(bg_taxa) - 4) if len(bg_taxa) > 4 else "")
    # Where every read ended up, target first, so a reviewer can see what the
    # reads were weighed against and not only the headline percentage.
    # Keyed by reference tag (<taxid>__<code>__<index>), not by taxon name:
    # several references of one species share a name, and the point here is
    # which reference each read preferred.
    by_ref = ";".join(
        "%s%s:%d" % (t, " [target]" if t == target_tag else
                     (" [same taxon]" if t in set(same_taxon_tags) else " [other taxon]"), n)
        for t, n in sorted(landed.items(),
                           key=lambda kv: (kv[0] != target_tag, -kv[1], kv[0])))
    return {"pct_target": pct_target,
            "pct_identity_best": (round(sum(best_ident.values()) / len(best_ident), 2)
                                  if best_ident else ""),
            "n_closest_genome": len(best_ident),
            "n_unique_target": unique.get(TARGET, 0),
            "n_unique_other": sum(v for k, v in unique.items() if k != TARGET),
            "n_ambiguous": len(ambiguous),
            "reads_by_taxon": by_taxon,
            "competing_best": ("%s = %s (%d reads)"
                               % (best, labels.get(best, best), elsewhere[best]))
                              if best else "",
            "n_competitors": n_competitors,
            "competitors": comp_desc,
            "reads_by_reference": by_ref}


def _f(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def decoy_references(taxid, panel_by_taxid, cache_dir, entrez):
    """Fetch this taxon's panel genomes, cached across runs.

    Returns (paths, tag -> taxon label). Panel rows are the sibling species of
    the same genus, plus clustered strain representatives of the taxon itself,
    so a true read has somewhere right to go and a read from a neighbour has
    somewhere better to go.
    """
    rows = panel_by_taxid.get(str(taxid), [])
    if not rows or entrez is None:
        return [], {}, {}
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    paths, labels, taxon_of = [], {}, {}
    for row in rows:
        accs = [a for a in (row.get("accessions") or "").split(":") if a]
        if not accs:
            continue
        tag = "panel_%s_%s" % (row.get("role", "sibling"), accs[0].replace(".", "_"))
        fa = cache_dir / (tag + ".fasta")
        if not fa.exists():
            try:
                text = entrez.efetch_fasta(accs)
            except Exception as exc:  # noqa: BLE001
                log.warning("  panel genome %s unavailable: %s", accs[0], exc)
                continue
            tmp = fa.with_name(fa.name + ".part")
            tmp.write_text(text if text.endswith("\n") else text + "\n")
            tmp.rename(fa)
        paths.append(fa)
        labels[tag] = row.get("name") or row.get("taxid", "")
        # Strain genomes of the taxon itself are on-target; siblings are not.
        taxon_of[tag] = ("__same__" if row.get("role") == "strain"
                         else (row.get("name") or row.get("taxid", "")))
    return paths, labels, taxon_of


def review(summary_tsv, out_path, outdir, params, library_sizes=None, threads=2,
           control_regex=None, panel_by_taxid=None, entrez=None, cache_dir=None,
           background=None):
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
            decoys, dec_labels, dec_taxon = decoy_references(
                taxid, panel_by_taxid or {}, Path(cache_dir or outdir) / "panel", entrez)
            # A strain genome of this same taxon counts as on-target.
            kin = kin + [t for t, v in dec_taxon.items() if v == "__same__"]
            labels = dict(tag_labels, **dec_labels)
            taxon_of = {t: v for t, v in dec_taxon.items() if v != "__same__"}
            try:
                spec = reference_specificity(align_dir, sample_dir, tag,
                                             work_root / sample_id / tag, threads=threads,
                                             same_taxon_tags=kin, tag_labels=labels,
                                             decoys=decoys, tag_taxon=taxon_of,
                                             margin=params.evidence_score_margin,
                                             em_rounds=params.evidence_em_rounds,
                                             background=background)
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

        # Identity to the reference that was aligned, and to the closest genome
        # of the same taxon in the panel. For a divergent species the two differ
        # sharply: a Lassa clade II read is ~80% identical to a Josiah reference
        # and ~99% to its own clade. Judge on the closest genome, so a true read
        # of an under-represented lineage is not failed for being divergent;
        # fall back to the aligned reference when no panel genome exists.
        ident_ref = m["pct_identity_mean"] if m else ""
        ident_best = (spec or {}).get("pct_identity_best", "")
        ident = ident_best if ident_best != "" else ident_ref
        identity_v = score(None if ident == "" else ident >= params.evidence_min_identity,
                           "identity",
                           "reads sit at %s%% identity to the closest genome of this taxon, "
                           "under %g%%" % (ident, params.evidence_min_identity))

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

        # With nothing to compete, every read is trivially "attributable" and
        # the criterion would pass while testing nothing. ASFV is the real
        # case: Asfarviridae holds one genus holding one species, so no sibling
        # exists anywhere to offer. Say so instead of scoring it.
        spec_ok = None
        if spec and spec["pct_target"] != "" and spec.get("n_competitors", 0) > 0:
            spec_ok = _f(spec["pct_target"]) >= params.evidence_min_retained
        specificity_v = score(spec_ok, "specificity",
                              "only %s%% of the reads are attributable to this taxon once "
                              "related genomes compete; %s takes the rest"
                              % (spec["pct_target"] if spec else "",
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
            "pct_identity_closest_genome": ident_best,
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
            "pct_reads_attributed_to_target": spec["pct_target"] if spec else "",
            "n_reads_unique_to_target": spec["n_unique_target"] if spec else "",
            "n_reads_unique_elsewhere": spec["n_unique_other"] if spec else "",
            "n_reads_ambiguous": spec["n_ambiguous"] if spec else "",
            "reads_by_taxon": spec["reads_by_taxon"] if spec else "",
            "n_competing_references": spec["n_competitors"] if spec else "",
            "competing_references_tested": spec["competitors"] if spec else "",
            "reads_by_reference": spec["reads_by_reference"] if spec else "",
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
