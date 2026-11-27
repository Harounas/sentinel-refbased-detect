"""The background layer: what a read can be, independent of the priority table.

The sibling and strain layers are both derived from Refgenomes_V15.tsv, and
that is their limitation. They can say "these reads fit Lassa better than its
neighbours", but a read from a virus nobody put on the list, or from the host,
or from a reagent organism, has nowhere to go and is attributed to whatever
reference it happened to be aligned against. At ten or twenty reads that is the
commonest way a false signal survives.

This layer is built from NCBI directly and never consults the priority table:

  viral        every RefSeq viral genome, 14,844 species in about 580 Mb,
               taken from the aggregated release files rather than 15,865
               separate assembly downloads
  host         the human genome, where most unexplained reads actually come
               from
  contaminant  phiX and the usual reagent organisms

It is one index, built once and cached, not one per pathogen: 3.6 Gb is far too
much to rebuild for every row under review.
"""

import gzip
import os
import re
import shutil
import urllib.request
from pathlib import Path

from .utils import PipelineError, download, log, run

REFSEQ_VIRAL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
ACC2TAXID = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"

# Sequencing controls and the organisms that turn up in extraction kits. phiX
# is spiked into most Illumina runs; the rest are the reagent contaminants
# repeatedly reported from low-biomass samples.
CONTAMINANTS = {
    "NC_001422.1": "Escherichia phage phiX174 (sequencing spike-in)",
    "NC_002695.2": "Escherichia coli O157:H7",
    "NC_000913.3": "Escherichia coli K-12",
    "NC_003098.1": "Streptococcus pneumoniae",
    "NZ_CP009617.1": "Cutibacterium acnes",
    "NC_004129.6": "Pseudomonas fluorescens",
    "NC_007204.1": "Psychrobacter arcticus",
    "NC_008044.1": "Ruegeria sp. TM1040",
    "NC_010162.1": "Sorangium cellulosum",
}


def viral_release_files(index_url=REFSEQ_VIRAL):
    """The aggregated RefSeq viral FASTA files, discovered from the directory.

    There are a handful of viral.N.1.genomic.fna.gz holding every RefSeq viral
    sequence between them. Downloading those beats walking 15,865 assemblies.
    """
    with urllib.request.urlopen(index_url, timeout=120) as fh:
        html = fh.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r'href="(viral\.\d+\.1\.genomic\.fna\.gz)"', html)))
    if not names:
        raise PipelineError("no viral.N.1.genomic.fna.gz found at %s" % index_url)
    return [index_url + n for n in names]


def fetch_viral(out_dir, index_url=REFSEQ_VIRAL):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for url in viral_release_files(index_url):
        dest = out_dir / url.rsplit("/", 1)[-1]
        if not dest.exists():
            tmp = dest.with_name(dest.name + ".part")
            download(url, tmp)
            tmp.rename(dest)
        paths.append(dest)
    return paths


def fetch_contaminants(out_dir, entrez, accessions=None):
    """Contaminant genomes, by accession, through E-utilities."""
    accessions = list(accessions or CONTAMINANTS)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "contaminants.fasta"
    if dest.exists() and dest.stat().st_size:
        return dest
    chunks = []
    for acc in accessions:
        try:
            chunks.append(entrez.efetch_fasta([acc]))
        except Exception as exc:  # noqa: BLE001
            log.warning("contaminant %s unavailable: %s", acc, exc)
    if not chunks:
        raise PipelineError("no contaminant genomes could be fetched")
    dest.write_text("".join(chunks))
    return dest


def _open_any(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def concatenate(sources, out_fasta, label_of=None):
    """Write one FASTA, recording every sequence id as it goes.

    Sequence ids are left exactly as they are: the accession is what the
    taxid map is keyed by, and rewriting headers here would break it.
    """
    out_fasta = Path(out_fasta)
    seen = set()
    n = 0
    with open(out_fasta, "w") as out:
        for src in sources:
            with _open_any(src) as fh:
                keep = True
                for line in fh:
                    if line.startswith(">"):
                        sid = line[1:].split()[0]
                        # A genome present in two sources would otherwise be
                        # indexed twice and split its own reads.
                        keep = sid not in seen
                        if keep:
                            seen.add(sid)
                            n += 1
                    if keep:
                        out.write(line)
            log.info("  %s: %d sequences so far", Path(src).name, n)
    return out_fasta, seen


def accession_taxids(accessions, out_path, acc2taxid_files=None, cache_dir=None,
                     url=ACC2TAXID):
    """accession -> taxid, for the accessions in the background only.

    nucl_gb.accession2taxid has hundreds of millions of lines, so it is
    streamed and filtered rather than loaded.
    """
    out_path = Path(out_path)
    files = list(acc2taxid_files or [])
    if not files:
        cache_dir = Path(cache_dir or out_path.parent)
        cache_dir.mkdir(parents=True, exist_ok=True)
        dest = cache_dir / "nucl_gb.accession2taxid.gz"
        if not dest.exists():
            tmp = dest.with_name(dest.name + ".part")
            download(url, tmp)
            tmp.rename(dest)
        files = [dest]

    wanted = set(accessions)
    found = {}
    for path in files:
        with _open_any(path) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) >= 3 and f[1] in wanted:
                    found[f[1]] = f[2].strip()
        log.info("  %s: %d of %d accessions mapped", Path(path).name, len(found), len(wanted))
    with open(out_path, "w") as out:
        for acc in sorted(wanted):
            out.write("%s\t%s\n" % (acc, found.get(acc, "")))
    missing = len(wanted) - len(found)
    if missing:
        log.warning("%d background sequences have no taxid; they still compete "
                    "for reads but are reported by accession", missing)
    return out_path, found


def name_map(acc_taxid, dmp_dir, out_path):
    """accession -> taxid -> scientific name, so a competitor can be named."""
    want = {t for t in acc_taxid.values() if t}
    names = {}
    with open(Path(dmp_dir) / "names.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = [x.strip() for x in line.split("|")]
            if len(p) >= 4 and p[3] == "scientific name" and p[0] in want:
                names[p[0]] = p[1]
    with open(out_path, "w") as out:
        out.write("sequence_id\ttaxid\tname\n")
        for acc, taxid in sorted(acc_taxid.items()):
            out.write("%s\t%s\t%s\n" % (acc, taxid, names.get(taxid, "")))
    return out_path


def build(out_dir, entrez=None, host_fasta=None, include_viral=True,
          include_contaminants=True, acc2taxid_files=None, dmp_dir=None,
          cache_dir=None, index=True, threads=4):
    """Assemble the background FASTA, its taxid map, and (optionally) the index."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(cache_dir or out_dir / "downloads")

    sources = []
    if include_viral:
        log.info("fetching RefSeq viral release files")
        sources += fetch_viral(cache_dir / "viral")
    if include_contaminants:
        if entrez is None:
            log.warning("no NCBI contact given; skipping contaminant genomes")
        else:
            log.info("fetching %d contaminant genomes", len(CONTAMINANTS))
            sources.append(fetch_contaminants(cache_dir / "contaminants", entrez))
    if host_fasta:
        log.info("including host genome %s", host_fasta)
        sources.append(Path(host_fasta))
    if not sources:
        raise PipelineError("nothing to build: enable viral, contaminants or a host genome")

    fasta = out_dir / "background.fasta"
    log.info("writing %s", fasta)
    fasta, ids = concatenate(sources, fasta)
    log.info("%d sequences, %.1f Gb", len(ids), fasta.stat().st_size / 1e9)

    acc_path, acc_taxid = accession_taxids(ids, out_dir / "background.acc2taxid.tsv",
                                           acc2taxid_files=acc2taxid_files,
                                           cache_dir=cache_dir)
    if dmp_dir:
        name_map(acc_taxid, dmp_dir, out_dir / "background.taxa.tsv")

    if index:
        log.info("building the bwa index; this is the slow step and happens once")
        run(["bwa", "index", str(fasta)], stderr=out_dir / "bwa_index.log")
    return fasta
