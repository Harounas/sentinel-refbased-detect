"""Run the whole workflow locally, mirroring the WDL's call graph."""

import csv
import json
import logging
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from . import __version__, WDL_METHOD
from .align import align_and_consensus
from .calls import final_calls
from .candidates import select_candidates
from .discover import discover_extra_taxa
from .kraken import classify, library_size, prepare_db
from .ncbi import Entrez
from .references import fetch_references
from .samples import discover, from_samplesheet, merge_lanes
from .summarize import merge_tables, summarize_sample
from .taxonomy import extract_dmp, get_taxdump, resolve_taxid_aliases
from .utils import CORE_TOOLS, PipelineError, localize, log, require_tools


def setup_logging(outdir, verbose=False):
    outdir.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S")
    log.setLevel(logging.DEBUG)
    log.handlers.clear()
    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG if verbose else logging.INFO)
    sh.setFormatter(fmt)
    fh = logging.FileHandler(outdir / "pipeline.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    log.addHandler(sh)
    log.addHandler(fh)


def find_report(report_dir, sample_id):
    """Pick <sample_id>*report*.txt from a directory of existing Kraken2 reports."""
    cands = sorted(p for p in Path(report_dir).iterdir()
                   if p.is_file() and "report" in p.name
                   and p.name.split(".")[0] == sample_id)
    return cands[0] if cands else None


def read_counts_file(path):
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and row[1].strip().isdigit():
                out[row[0].strip()] = int(row[1])
    return out


def run_pipeline(args, params):
    t0 = time.time()
    outdir = Path(args.outdir).expanduser().resolve()
    setup_logging(outdir, args.verbose)
    log.info("sentinel-refbased-detect %s (port of %s)", __version__, WDL_METHOD)

    work = outdir / "work"
    cache = Path(args.cache_dir).expanduser().resolve() if args.cache_dir else outdir / "cache"
    for d in (work, cache):
        d.mkdir(parents=True, exist_ok=True)

    tools = list(CORE_TOOLS)
    if args.kraken_reports:
        tools.remove("kraken2")
    if args.diamond_db:
        tools.append("diamond")
    if args.host_ref:
        tools.append("minimap2")
    require_tools(tools)

    entrez = Entrez(args.ncbi_email, args.ncbi_api_key)

    # ---- inputs ------------------------------------------------------------
    samples = from_samplesheet(args.samplesheet) if args.samplesheet else \
        discover(args.fastq_dir, recursive=args.recursive)
    log.info("%d samples: %s", len(samples), ", ".join(s.sample_id for s in samples))
    ref_tsv = localize(args.ref_tsv, cache / "inputs")
    taxdump_src = localize(args.taxdump, cache / "inputs") if args.taxdump else None

    with open(outdir / "run_params.json", "w") as fh:
        json.dump({"version": __version__, "wdl": WDL_METHOD,
                   "inputs": {"kraken_db": str(args.kraken_db or ""),
                              "kraken_reports": str(args.kraken_reports or ""),
                              "ref_tsv": str(args.ref_tsv), "fastq_dir": str(args.fastq_dir or ""),
                              "samplesheet": str(args.samplesheet or ""),
                              "host_ref": str(args.host_ref or ""),
                              "diamond_db": str(args.diamond_db or "")},
                   "samples": {s.sample_id: [str(f) for f in s.files] for s in samples},
                   "params": asdict(params)}, fh, indent=2)

    # ---- resolve_priority --------------------------------------------------
    tax_dir = outdir / "taxonomy"
    taxdump = get_taxdump(taxdump_src, cache / "taxonomy")
    dmp = extract_dmp(taxdump, cache / "taxonomy" / "dmp")
    prio = resolve_taxid_aliases(ref_tsv, tax_dir / "priority", dmp)

    # ---- kraken2_classify --------------------------------------------------
    reads = {}
    for s in samples:
        reads[s.sample_id] = merge_lanes(s, work / "reads" / s.sample_id)

    # ---- optional host depletion -------------------------------------------
    if args.host_ref:
        from .host import deplete, host_index, write_summary
        index = host_index(args.host_ref, cache / "host", threads=args.threads)
        stats = []
        for s in samples:
            r1, r2 = reads[s.sample_id]
            d1, d2, st = deplete(s.sample_id, r1, r2, index, outdir / "host_depleted",
                                 threads=args.threads, resume=args.resume)
            reads[s.sample_id] = (d1, d2)
            stats.append(st)
        write_summary(stats, outdir / "host_depletion.tsv")

    kdir = outdir / "kraken2"
    reports = {}
    if args.kraken_reports:
        for s in samples:
            rep = find_report(args.kraken_reports, s.sample_id)
            if rep is None:
                raise PipelineError("no Kraken2 report for %s in %s" % (s.sample_id, args.kraken_reports))
            kdir.mkdir(parents=True, exist_ok=True)
            dest = kdir / ("%s.kraken2.report.txt" % s.sample_id)
            shutil.copyfile(rep, dest)
            reports[s.sample_id] = dest
    else:
        db = prepare_db(args.kraken_db, cache)
        for i, s in enumerate(samples, 1):
            dest = kdir / ("%s.kraken2.report.txt" % s.sample_id)
            if args.resume and dest.exists():
                log.info("kraken2 [%d/%d] %s: report exists, skipping", i, len(samples), s.sample_id)
                reports[s.sample_id] = dest
                continue
            log.info("kraken2 [%d/%d] %s", i, len(samples), s.sample_id)
            r1, r2 = reads[s.sample_id]
            reports[s.sample_id] = classify(s, r1, r2, db, kdir, threads=args.threads,
                                            memory_mapping=args.kraken_memory_mapping,
                                            keep_output=args.keep_kraken_output)

    # ---- discover_extra_taxa + resolve_all -----------------------------------
    disc_dir = outdir / "discovery"
    disc_dir.mkdir(parents=True, exist_ok=True)
    extra, discovered = discover_extra_taxa([str(reports[s.sample_id]) for s in samples],
                                            ref_tsv, prio["taxid_aliases"], disc_dir, entrez, params)
    allres = resolve_taxid_aliases(ref_tsv, tax_dir / "all", dmp, extra_refs_tsv=extra)

    # ---- select_candidates + fetch_references ---------------------------------
    jobs = []
    per_sample = {}
    for s in samples:
        sdir = outdir / "samples" / s.sample_id
        sdir.mkdir(parents=True, exist_ok=True)
        cand = select_candidates(s.sample_id, reports[s.sample_id], allres["combined_refs"],
                                 allres["taxid_aliases"], allres["name_aliases"],
                                 allres["taxon_levels"], sdir / "candidates.tsv", params,
                                 name_regex_tsv=args.name_regex_tsv)
        refs_dir = sdir / "refs"
        if refs_dir.exists():
            shutil.rmtree(refs_dir)
        fastas = fetch_references(cand, refs_dir, cache / "references", entrez)
        per_sample[s.sample_id] = {"dir": sdir, "candidates": cand, "align_dirs": []}
        r1, r2 = reads[s.sample_id]
        for fa in fastas:
            adir = sdir / "align" / fa.name[:-len(".fasta")]
            per_sample[s.sample_id]["align_dirs"].append(adir)
            jobs.append((s.sample_id, r1, r2, fa, adir))

    # ---- align_and_consensus (parallel) ------------------------------------
    log.info("%d alignment jobs, %d at a time, %d threads each", len(jobs), args.jobs, args.threads)
    failures = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(align_and_consensus, sid, r1, r2, fa, adir, params,
                          threads=args.threads, keep_bam=not args.no_bam,
                          keep_intermediate=args.keep_intermediate,
                          resume=args.resume): (sid, adir)
                for sid, r1, r2, fa, adir in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            sid, adir = futs[fut]
            try:
                fut.result()
                log.info("align [%d/%d] %s %s done", n, len(jobs), sid, adir.name)
            except Exception as exc:  # noqa: BLE001
                log.error("align [%d/%d] %s %s FAILED: %s (see %s/align.log)",
                          n, len(jobs), sid, adir.name, exc, adir)
                failures.append((sid, adir))
    if failures and not args.keep_going:
        raise PipelineError("%d alignment job(s) failed; fix and rerun (finished jobs are kept), "
                            "or pass --keep-going to summarise without them" % len(failures))
    failed = {a for _, a in failures}

    # ---- summarize_sample + merge_results ----------------------------------
    summaries = {}
    for s in samples:
        ps = per_sample[s.sample_id]
        dirs = [d for d in ps["align_dirs"] if d not in failed]
        summaries[s.sample_id] = summarize_sample(s.sample_id, dirs, ps["candidates"], ps["dir"], params)

    precall = outdir / "detection_summary.tsv"
    by_region = outdir / "detection_by_region.tsv"
    n_rows = merge_tables([summaries[s.sample_id]["summary"] for s in samples], precall)
    merge_tables([summaries[s.sample_id]["by_region"] for s in samples], by_region)
    log.info("%d summary rows", n_rows)

    # ---- DIAMOND (optional) -------------------------------------------------
    summary_before_calls = precall
    if args.diamond_db:
        from . import diamond as dm
        ddir = outdir / "diamond"
        ddir.mkdir(parents=True, exist_ok=True)
        fasta, pmap = dm.pool_blocks([summaries[s.sample_id]["blocks_fasta"] for s in samples],
                                     [summaries[s.sample_id]["block_map"] for s in samples], ddir)
        hits = dm.blastx(fasta, Path(args.diamond_db).expanduser(), ddir, params, threads=args.threads)
        if args.diamond_taxdump:
            dtd = extract_dmp(localize(args.diamond_taxdump, cache / "inputs"),
                              cache / "taxonomy" / "diamond_dmp")
        else:
            dtd = dmp
        acc = dm.accession_taxids(hits, args.accession2taxid, ddir)
        ctax, _ = dm.assign_taxonomy(hits, pmap, acc, dtd, params.diamond_top_pct, ddir)
        summary_before_calls = dm.annotate_summary(precall, ctax,
                                                   outdir / "detection_summary_annotated.tsv")

    # ---- final_calls --------------------------------------------------------
    if args.no_normalize:
        lib = {}
    elif args.read_counts:
        lib = read_counts_file(args.read_counts)
    else:
        lib = {s.sample_id: library_size(reports[s.sample_id], reads[s.sample_id][1] is not None)
               for s in samples}
    with open(outdir / "library_sizes.tsv", "w") as fh:
        for s in samples:
            fh.write("%s\t%s\n" % (s.sample_id, lib.get(s.sample_id, "")))

    final = outdir / "detection_summary_final.tsv"
    if n_rows == 0:
        log.warning("no virus reached alignment in any sample; final table is empty")
        shutil.copyfile(summary_before_calls, final)
        (outdir / "control_background.tsv").write_text("")
        counts = {}
    else:
        counts = final_calls(summary_before_calls, final, outdir / "control_background.tsv",
                             params, library_sizes=lib)

    # ---- evidence review of failing priority rows ---------------------------
    # Runs after the calls, on the alignments that already exist. It reads the
    # final table and writes its verdict back into it, leaving `detection` and
    # `final_call` untouched.
    if args.evidence_review and n_rows and not args.no_bam:
        from . import evidence, panel
        panel_rows = panel.load_panel(args.decoy_panel) if args.decoy_panel else {}
        if args.decoy_panel and not panel_rows:
            log.info("no decoy panel at %s; the specificity check uses this sample's own "
                     "references only. Build one with: sentinel-detect panel build",
                     args.decoy_panel)
        verdicts = evidence.review(final, outdir / "evidence_review.tsv", outdir, params,
                                   library_sizes=lib, threads=args.threads,
                                   panel_by_taxid=panel_rows, entrez=entrez, cache_dir=cache)
        # Always annotate, so the columns exist whether or not anything was
        # reviewed; rows that were not reviewed say so.
        evidence.annotate(final, verdicts, outdir / "detection_summary_final.annotated.tsv")
        shutil.move(str(outdir / "detection_summary_final.annotated.tsv"), str(final))
    elif args.evidence_review and args.no_bam:
        log.warning("evidence review skipped: it needs the per-reference BAMs that --no-bam removes")

    if not args.keep_intermediate:
        shutil.rmtree(work / "reads", ignore_errors=True)

    log.info("done in %.1f min. Final table: %s", (time.time() - t0) / 60, final)
    return final, counts
