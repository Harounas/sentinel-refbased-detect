"""Task discover_extra_taxa.

Viruses outside the priority table that the run keeps showing: seen in many
samples, or seen once with real depth. A RefSeq reference is pulled for each,
and from then on they are treated like priority taxa, with
virus_importance = "discovered" telling them apart.
"""

import csv
import os
import re
from pathlib import Path

from .utils import log

VIRUSES_TAXID = "10239"


def parse_viral_species(path):
    """Species nodes inside the Viruses subtree, with their clade reads,
    read off the report's own indentation."""
    out = {}
    viral_indent = None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            if len(f) >= 8 and f[6].strip().isdigit():
                taxid, rank, name_field = f[6].strip(), f[5].strip(), f[7]
            elif f[4].strip().isdigit():
                taxid, rank, name_field = f[4].strip(), f[3].strip(), f[5]
            else:
                continue
            indent = len(name_field) - len(name_field.lstrip(" "))
            if taxid == VIRUSES_TAXID:
                viral_indent = indent
                continue
            if viral_indent is None:
                continue
            if indent <= viral_indent:
                viral_indent = None
                continue
            if rank != "S":
                continue
            try:
                out[taxid] = (int(f[1]), name_field.strip())
            except ValueError:
                pass
    return out


def pick_reference(entrez, taxid):
    """One RefSeq genome for a taxon, or every segment when it is segmented."""
    term = "txid%s[Organism:exp] AND srcdb_refseq[PROP] AND 1000:1000000[SLEN]" % taxid
    res = entrez.json("esearch.fcgi", db="nuccore", term=term, retmax=50)
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return [], "no refseq genome"
    summ = entrez.json("esummary.fcgi", db="nuccore", id=",".join(ids))
    recs = []
    for uid in ids:
        d = summ.get("result", {}).get(uid)
        if not d:
            continue
        recs.append((d.get("accessionversion") or d.get("caption", ""),
                     int(d.get("slen") or 0), d.get("title", "")))
    recs = [r for r in recs if r[0]]
    if not recs:
        return [], "no usable esummary record"

    seg_of = {}
    for acc, slen, title in recs:
        m = re.search(r"segment\s+([A-Za-z0-9]+)", title, re.I)
        if m:
            key = m.group(1).upper()
            if key not in seg_of or slen > seg_of[key][1]:
                seg_of[key] = (acc, slen)
    if len(seg_of) >= 2:
        accs = [seg_of[k][0] for k in sorted(seg_of)]
        return accs, "refseq, %d segments" % len(accs)
    best = max(recs, key=lambda r: r[1])
    return [best[0]], "refseq %d bp" % best[1]


def code_for(taxid, name):
    tag = re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")[:24]
    return tag or ("taxid" + taxid)


def discover_extra_taxa(kraken_reports, ref_genomes_tsv, taxid_aliases_tsv, out_dir,
                        entrez, params):
    out_dir = Path(out_dir)
    known = set()
    with open(taxid_aliases_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            known.add(row["report_taxid"])
            known.add(row["table_taxid"])
    with open(ref_genomes_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].strip().isdigit():
                known.add(row[0].strip())

    n_samples_with, max_reads, names, total_reads = {}, {}, {}, {}
    for path in kraken_reports:
        if not os.path.exists(path):
            continue
        for taxid, (reads, name) in parse_viral_species(path).items():
            if taxid in known or reads < params.discover_present_min_reads:
                continue
            n_samples_with[taxid] = n_samples_with.get(taxid, 0) + 1
            max_reads[taxid] = max(max_reads.get(taxid, 0), reads)
            total_reads[taxid] = total_reads.get(taxid, 0) + reads
            names[taxid] = name

    qualified = []
    for taxid, n in n_samples_with.items():
        if n >= params.discover_min_samples:
            why = "in %d samples" % n
        elif max_reads[taxid] >= params.discover_min_reads_one_sample:
            why = "%d reads in one sample" % max_reads[taxid]
        else:
            continue
        qualified.append((taxid, n, max_reads[taxid], total_reads[taxid], names[taxid], why))

    qualified.sort(key=lambda r: (-r[1], -r[2]))
    dropped = qualified[params.discover_max_novel_taxa:]
    qualified = qualified[:params.discover_max_novel_taxa]

    extra_rows, audit = [], []
    for taxid, n, mx, tot, name, why in qualified:
        try:
            accs, how = pick_reference(entrez, taxid)
        except Exception as exc:  # noqa: BLE001
            accs, how = [], "lookup failed: %s" % exc
        audit.append({"taxid": taxid, "name": name, "n_samples": n,
                      "max_reads_one_sample": mx, "total_reads": tot,
                      "qualified_by": why, "reference": ":".join(accs) or "",
                      "reference_source": how,
                      "status": "included" if accs else "skipped"})
        if accs:
            extra_rows.append([taxid, code_for(taxid, name), name, ":".join(accs), "discovered"])

    for taxid, n, mx, tot, name, why in dropped:
        audit.append({"taxid": taxid, "name": name, "n_samples": n,
                      "max_reads_one_sample": mx, "total_reads": tot,
                      "qualified_by": why, "reference": "", "reference_source": "",
                      "status": "over max_novel_taxa cap"})

    extra_path = out_dir / "extra_refs.tsv"
    with open(extra_path, "w", newline="") as out:
        csv.writer(out, delimiter="\t").writerows(extra_rows)

    cols = ["taxid", "name", "n_samples", "max_reads_one_sample", "total_reads",
            "qualified_by", "reference", "reference_source", "status"]
    audit_path = out_dir / "discovered_taxa.tsv"
    with open(audit_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(audit)

    log.info("%d reports, %d non-priority viral species seen, %d qualified, %d with a reference",
             len(kraken_reports), len(n_samples_with), len(qualified) + len(dropped), len(extra_rows))
    for r in audit[:40]:
        log.info("  %-9s %-44s samples=%-3d max=%-8d %-24s %s", r["taxid"], r["name"][:44],
                 r["n_samples"], r["max_reads_one_sample"], r["qualified_by"], r["status"])
    return extra_path, audit_path
