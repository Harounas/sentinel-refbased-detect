"""Command-line interface: `sentinel-detect run ...` and `sentinel-detect check`."""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path
from dataclasses import MISSING

from . import DEFAULT_REF_TSV, __version__, WDL_METHOD
from .params import Params, param_fields
from .utils import PipelineError

GROUP_TITLES = {
    "selection": "Kraken2 candidate selection",
    "discovery": "discovery of non-priority viruses",
    "detection": "alignment and detection rule",
    "consensus": "ivar consensus",
    "controls": "negative controls and final call",
    "diamond": "DIAMOND (only with --diamond-db)",
}


def _add_param_args(parser):
    groups = {}
    for f in param_fields():
        g = f.metadata.get("group", "detection")
        if g not in groups:
            groups[g] = parser.add_argument_group(GROUP_TITLES.get(g, g))
        flag = "--" + f.name.replace("_", "-")
        default = f.default if f.default is not MISSING else None
        help_ = f.metadata.get("help", "").replace("%", "%%")
        if f.type is bool:
            groups[g].add_argument(flag, dest=f.name, default=default,
                                   action=argparse.BooleanOptionalAction, help=help_)
        else:
            groups[g].add_argument(flag, dest=f.name, default=default, type=f.type, help=help_)


def build_parser():
    p = argparse.ArgumentParser(
        prog="sentinel-detect",
        description="Reference-based detection of priority pathogens from FASTQ "
                    "(local port of Terra %s, config V15)." % WDL_METHOD)
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="run the pipeline",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    req = r.add_argument_group("required inputs")
    req.add_argument("-k", "--kraken-db",
                     help="Kraken2 database: a directory with hash.k2d/opts.k2d/taxo.k2d, "
                          "or a .tar.gz/.tar.zst archive of one")
    req.add_argument("-r", "--ref-tsv", default=str(DEFAULT_REF_TSV),
                     help="priority pathogen table: taxid, code, name, colon-joined accessions; "
                          "no header. Local path, gs:// or https://. Default: the bundled "
                          "Refgenomes_V15.tsv")
    req.add_argument("-i", "--fastq-dir", help="directory of FASTQ files (*.fastq.gz, *_R1_001.fastq.gz, ...)")
    req.add_argument("-o", "--outdir", required=True, help="output directory")
    req.add_argument("--ncbi-email", default=os.environ.get("NCBI_EMAIL"),
                     help="e-mail sent to NCBI E-utilities (or set $NCBI_EMAIL)")

    opt = r.add_argument_group("optional inputs")
    opt.add_argument("--samplesheet", help="TSV/CSV with sample_id, fastq_1[, fastq_2]; use instead of --fastq-dir")
    opt.add_argument("--recursive", action="store_true", help="search --fastq-dir recursively")
    opt.add_argument("--ncbi-api-key", default=os.environ.get("NCBI_API_KEY"),
                     help="NCBI API key (or $NCBI_API_KEY); raises the rate limit")
    opt.add_argument("--host-ref",
                     help="host genome (FASTA[.gz], or a prebuilt minimap2 .mmi) to deplete "
                          "before Kraken2: read pairs where either mate maps are removed. "
                          "Leave unset if your FASTQs are already host-depleted")
    opt.add_argument("--taxdump",help="NCBI taxdump.tar.gz; default downloads the current one")
    opt.add_argument("--name-regex-tsv", help="optional regex<TAB>taxid overrides for report names")
    opt.add_argument("--kraken-reports",
                     help="directory of existing Kraken2 reports (<sample_id>.*report*.txt); "
                          "skips Kraken2, so --kraken-db is not needed")
    opt.add_argument("--read-counts", help="TSV sample_id<TAB>library reads, for reads-per-million "
                                           "control comparison; default derives them from the Kraken2 reports")
    opt.add_argument("--no-normalize", action="store_true",
                     help="compare samples to controls on raw read counts")
    opt.add_argument("--diamond-db", help="DIAMOND nr .dmnd; turns on blastx identification of detected blocks")
    opt.add_argument("--accession2taxid", nargs="*", default=[],
                     help="prot.accession2taxid(.gz) files; downloaded from NCBI if omitted")
    opt.add_argument("--diamond-taxdump", help="taxdump used for the DIAMOND LCA (default: --taxdump)")

    rt = r.add_argument_group("resources and behaviour")
    rt.add_argument("-t", "--threads", type=int, default=min(8, os.cpu_count() or 1),
                    help="threads per tool call")
    rt.add_argument("-j", "--jobs", type=int, default=2, help="alignment jobs run at once")
    rt.add_argument("--kraken-memory-mapping", action="store_true",
                    help="kraken2 --memory-mapping (less RAM, slower)")
    rt.add_argument("--cache-dir", help="where references/taxdump are cached (default: OUTDIR/cache)")
    rt.add_argument("--resume", default=True, action=argparse.BooleanOptionalAction,
                    help="reuse finished Kraken2 reports and alignments")
    rt.add_argument("--keep-going", action="store_true",
                    help="summarise even if some alignment jobs failed")
    rt.add_argument("--keep-intermediate", action="store_true", help="keep intermediate BAMs/indexes")
    rt.add_argument("--no-bam", action="store_true", help="do not keep unique.bam per alignment")
    rt.add_argument("--keep-kraken-output", action="store_true", help="keep per-read Kraken2 output (large)")
    rt.add_argument("-v", "--verbose", action="store_true")

    _add_param_args(r)

    sub.add_parser("check", help="check that external tools are installed")

    rf = sub.add_parser("refs", help="show or export the bundled reference table")
    rf.add_argument("--export", metavar="PATH",
                    help="copy the bundled Refgenomes_V15.tsv to PATH (to edit it)")
    return p


def refs_command(args):
    if args.export:
        shutil.copyfile(DEFAULT_REF_TSV, args.export)
        print("copied to %s" % args.export)
        return
    n, taxids = 0, set()
    with open(DEFAULT_REF_TSV) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f and f[0].strip().isdigit():
                n += 1
                taxids.add(f[0].strip())
    print(DEFAULT_REF_TSV)
    print("%d reference rows, %d taxids" % (n, len(taxids)))


def check_tools():
    ok = True
    for tool, args in [("kraken2", ["--version"]), ("bwa", []), ("samtools", ["--version"]),
                       ("ivar", ["version"]), ("diamond", ["version"]),
                       ("minimap2", ["--version"])]:
        path = shutil.which(tool)
        if not path:
            need = {"diamond": "--diamond-db", "minimap2": "--host-ref"}.get(tool)
            optional = need is not None
            print("%-9s MISSING%s" % (tool, " (optional, needed for %s)" % need if optional else ""))
            ok = ok and optional
            continue
        try:
            res = subprocess.run([tool, *args], capture_output=True, text=True, timeout=30)
            text = (res.stdout or res.stderr).strip().splitlines()
            ver = next((l for l in text if "ersion" in l or l[:1].isdigit()), text[0] if text else "")
        except Exception:  # noqa: BLE001
            ver = ""
        print("%-9s %s  %s" % (tool, path, ver.strip()))
    return ok


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "check":
        sys.exit(0 if check_tools() else 1)
    if args.command == "refs":
        refs_command(args)
        return
    if not args.fastq_dir and not args.samplesheet:
        sys.exit("error: give --fastq-dir or --samplesheet")
    if not args.kraken_db and not args.kraken_reports:
        sys.exit("error: give --kraken-db (or --kraken-reports to reuse existing reports)")
    if not args.ncbi_email:
        sys.exit("error: --ncbi-email (or $NCBI_EMAIL) is required by NCBI E-utilities")

    params = Params(**{f.name: getattr(args, f.name) for f in param_fields()})
    from .pipeline import run_pipeline
    try:
        run_pipeline(args, params)
    except PipelineError as exc:
        sys.exit("ERROR: %s" % exc)
    except KeyboardInterrupt:
        sys.exit("interrupted; rerun the same command to resume")


if __name__ == "__main__":
    main()
