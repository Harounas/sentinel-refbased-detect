"""Optional DIAMOND stage: pool_blocks, diamond_blastx_nr, assign_taxonomy,
annotate_summary. Switched on by --diamond-db."""

import csv
import gzip
import os
import shutil
from collections import defaultdict
from pathlib import Path

from .utils import download, log, run

NCBI_TAXONOMY = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"


# --------------------------------------------------------------------------- pool
def pool_blocks(block_fastas, block_maps, out_dir):
    out_dir = Path(out_dir)
    fasta, mapf = out_dir / "pooled_blocks.fasta", out_dir / "pooled_map.tsv"
    n = 0
    with open(fasta, "w") as out:
        for path in block_fastas:
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    out.write(line)
    rows, seen = [], set()
    for path in block_maps:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["query_id"] in seen:
                    continue
                seen.add(row["query_id"])
                rows.append(row)
    with open(mapf, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["query_id", "sample_id", "contig_id", "contig_length"])
        for r in rows:
            w.writerow([r["query_id"], r["sample_id"], r["contig_id"], r["contig_length"]])
    log.info("%d blocks pooled from %d samples", n, len(block_fastas))
    return fasta, mapf


# ------------------------------------------------------------------------ blastx
def blastx(contigs, diamond_db, out_dir, params, threads=8, basename="detected_blocks"):
    out_dir = Path(out_dir)
    hits = out_dir / ("%s.nr_hits.tsv" % basename)
    if not os.path.getsize(contigs):
        hits.write_text("")
    else:
        tmp = out_dir / "tmp"
        tmp.mkdir(exist_ok=True)
        cmd = ["diamond", "blastx", "--db", diamond_db, "--query", contigs, "--out", hits,
               "--outfmt", "6", "qseqid", "qlen", "sseqid", "pident", "length",
               "qstart", "qend", "evalue", "bitscore"]
        if params.diamond_mode:
            cmd.append(params.diamond_mode)
        cmd += ["--top", params.diamond_top_pct, "--max-hsps", "1",
                "--evalue", params.diamond_evalue,
                "--block-size", params.diamond_block_size,
                "--index-chunks", params.diamond_index_chunks,
                "--threads", threads, "--tmpdir", tmp]
        log.info("diamond blastx on %s", contigs)
        run(cmd, stderr=out_dir / "diamond.log")
        shutil.rmtree(tmp, ignore_errors=True)
    run(["gzip", "-f", hits])
    return Path(str(hits) + ".gz")


# ----------------------------------------------------------------- taxonomy/LCA
LINEAGE = ["domain", "phylum", "class", "order", "family", "genus", "species"]
DOMAIN_RANKS = {"superkingdom", "domain", "acellular root"}
PROTOZOA = {"Apicomplexa", "Euglenozoa", "Metamonada", "Amoebozoa", "Sar", "Discoba",
            "Parabasalia", "Cryptophyceae", "Haptista", "Ciliophora"}
INVERT = {"Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca", "Annelida"}


def load_taxdump(d):
    parent, rank, name = {}, {}, {}
    with open(os.path.join(d, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            t = int(f[0])
            parent[t], rank[t] = int(f[1]), f[2]
    with open(os.path.join(d, "names.dmp"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if f[3] == "scientific name":
                name[int(f[0])] = f[1]
    merged = {}
    p = os.path.join(d, "merged.dmp")
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\t|\n").split("\t|\t")
                if len(f) >= 2:
                    merged[int(f[0])] = int(f[1])
    return parent, rank, name, merged


class Taxonomy:
    def __init__(self, d):
        self.parent, self.rank, self.name, self.merged = load_taxdump(d)
        self._path = {}

    def resolve(self, t):
        t = self.merged.get(t, t)
        return t if t in self.parent else None

    def path(self, t):
        if t not in self._path:
            p, x = [], t
            while True:
                p.append(x)
                if x == 1 or self.parent[x] == x:
                    break
                x = self.parent[x]
            self._path[t] = p[::-1]
        return self._path[t]

    def lca(self, taxids):
        paths = [self.path(t) for t in taxids]
        out = 1
        for level in zip(*paths):
            if all(x == level[0] for x in level):
                out = level[0]
            else:
                break
        return out

    def lineage(self, t):
        cols = {c: "" for c in LINEAGE}
        if not t:
            return cols
        for x in self.path(t):
            r = self.rank[x]
            if r in DOMAIN_RANKS and not cols["domain"]:
                cols["domain"] = self.name.get(x, "")
            elif r in cols and r != "domain":
                cols[r] = self.name.get(x, "")
        return cols

    def category(self, t):
        if not t:
            return "Unclassified"
        names = {self.name.get(x, "") for x in self.path(t)}
        nm = self.name.get(t, "")
        if nm in ("root", "cellular organisms"):
            return "Root"
        for key, cat in (("Viruses", "Virus"), ("Bacteria", "Bacteria"), ("Archaea", "Archaea")):
            if key in names:
                return cat
        if "Eukaryota" in names:
            if "Homo sapiens" in names:
                return "Human"
            if "Chordata" in names:
                return "Vertebrate"
            if "Fungi" in names:
                return "Fungi"
            if names & PROTOZOA:
                return "Protozoa"
            if names & INVERT:
                return "Invertebrate"
            if "Viridiplantae" in names:
                return "Plant"
            return "Other eukaryote"
        return "Other"


def _open_any(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def accession_taxids(hits_gz, accession2taxid, out_dir):
    """accession -> taxid for the accessions that were hit. Uses awk, because
    prot.accession2taxid has over a billion lines."""
    out_dir = Path(out_dir)
    accs = set()
    with _open_any(hits_gz) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                accs.add(f[2])
    hit_accs = out_dir / "hit_accs.txt"
    hit_accs.write_text("".join(a + "\n" for a in sorted(accs)))
    log.info("%d distinct hit accessions", len(accs))

    files = list(accession2taxid or [])
    if not files and accs:
        cache = out_dir / "accession2taxid"
        cache.mkdir(exist_ok=True)
        for f in ("prot", "pdb", "dead_prot"):
            dest = cache / ("%s.accession2taxid.gz" % f)
            if not dest.exists():
                tmp = cache / (dest.name + ".part")
                download("%s/accession2taxid/%s" % (NCBI_TAXONOMY, dest.name), tmp)
                tmp.rename(dest)
            files.append(dest)

    acc_taxid = out_dir / "acc_taxid.tsv"
    acc_taxid.write_text("")
    if accs:
        for f in files:
            with open(acc_taxid, "a") as out:
                run(["sh", "-c",
                     "gzip -dcf \"$1\" | awk -F'\\t' 'NR==FNR{h[$1]; next} ($2 in h){print $2\"\\t\"$3}' \"$0\" -",
                     hit_accs, f], stdout=out)
    return acc_taxid


def assign_taxonomy(hits_gz, contig_map, acc_taxid_tsv, taxdump_dir, top_pct, out_dir,
                    prefix="detected_blocks"):
    out_dir = Path(out_dir)
    top = float(top_pct) / 100.0
    tax = Taxonomy(taxdump_dir)

    acc2tax = {}
    with open(acc_taxid_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                acc2tax[f[0]] = int(f[1])

    hits = defaultdict(list)
    with _open_any(hits_gz) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            hits[f[0]].append((float(f[8]), f[2], float(f[3])))

    with open(contig_map) as fh:
        contigs = list(csv.DictReader(fh, delimiter="\t"))

    ctab_path = out_dir / ("%s.contig_taxonomy.tsv" % prefix)
    agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    sample_bp = defaultdict(int)
    no_taxid_acc = set()
    with open(ctab_path, "w", newline="") as ctab:
        cw = csv.writer(ctab, delimiter="\t", lineterminator="\n")
        cw.writerow(["sample_id", "contig_id", "contig_length", "n_hits", "n_hits_lca",
                     "best_hit", "best_hit_taxid", "best_hit_taxon", "best_pident", "best_bitscore",
                     "best_hit_category", "taxid", "taxon_name", "rank", "category"] + LINEAGE)
        for c in contigs:
            q, sid, clen = c["query_id"], c["sample_id"], int(c["contig_length"])
            sample_bp[sid] += clen
            h = sorted(hits.get(q, []), reverse=True)
            best = h[0] if h else None
            tids, used = [], 0
            if h:
                cutoff = best[0] * (1 - top)
                for bs, s, _ in h:
                    if bs < cutoff:
                        break
                    t = acc2tax.get(s)
                    t = tax.resolve(t) if t else None
                    if t is None:
                        no_taxid_acc.add(s)
                        continue
                    tids.append(t)
                    used += 1
            lca = tax.lca(tids) if tids else 0
            bt = tax.resolve(acc2tax.get(best[1], 0)) if best else None
            agg[sid][lca][0] += 1
            agg[sid][lca][1] += clen
            lin = tax.lineage(lca)
            cw.writerow([sid, c["contig_id"], clen, len(h), used,
                         best[1] if best else "", bt or "", tax.name.get(bt, "") if bt else "",
                         best[2] if best else "", best[0] if best else "",
                         tax.category(bt) if bt else "",
                         lca, tax.name.get(lca, "unclassified") if lca else "unclassified",
                         tax.rank.get(lca, "unclassified") if lca else "unclassified",
                         tax.category(lca)] + [lin[k] for k in LINEAGE])

    long_path = out_dir / ("%s.diamond_long.tsv" % prefix)
    with open(long_path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "percentage_bp", "contigs_clade", "contigs_taxon", "bp_clade",
                    "bp_taxon", "taxid", "taxon_name", "category", "rank"] + LINEAGE)
        for sid in sorted(agg):
            clade = defaultdict(lambda: [0, 0])
            for t, (n, bp) in agg[sid].items():
                for x in (tax.path(t) if t else [0]):
                    clade[x][0] += n
                    clade[x][1] += bp
            total = sample_bp[sid] or 1
            children = defaultdict(list)
            for t in clade:
                if t not in (0, 1):
                    children[tax.parent[t]].append(t)
            order = [0] if 0 in clade else []
            stack = [1] if 1 in clade else []
            while stack:
                t = stack.pop()
                order.append(t)
                stack.extend(sorted(children[t], key=lambda x: clade[x][1]))
            for t in order:
                n, bp = clade[t]
                own = agg[sid].get(t, [0, 0])
                lin = tax.lineage(t)
                w.writerow([sid, round(100.0 * bp / total, 4), n, own[0], bp, own[1], t,
                            tax.name.get(t, "unclassified") if t else "unclassified",
                            tax.category(t), tax.rank.get(t, "unclassified") if t else "unclassified"]
                           + [lin[k] for k in LINEAGE])

    log.info("%d blocks, %d with hits, %d hit accessions without taxid",
             len(contigs), sum(1 for c in contigs if c["query_id"] in hits), len(no_taxid_acc))
    run(["gzip", "-f", ctab_path])
    return Path(str(ctab_path) + ".gz"), long_path


# ---------------------------------------------------------------- annotate
ANNOT_COLS = ["blast_block_bp", "blast_region", "blast_n_hits", "blast_n_hits_lca",
              "blast_best_subject", "blast_best_subject_taxon",
              "blast_best_subject_category", "blast_pct_identity",
              "blast_bitscore", "diamond_taxid", "diamond_taxon", "diamond_rank",
              "diamond_category", "diamond_species", "diamond_genus",
              "diamond_family", "blast_agreement"]
VAGUE = {"", "root", "cellular organisms", "Viruses", "unclassified"}


def _norm(s):
    return (s or "").strip().lower()


def agreement(row, t):
    if t is None:
        return "no nr hit"
    name = (t.get("taxon_name") or "").strip()
    if name in VAGUE:
        return "unresolved (%s)" % (name or "no LCA")
    d_sp, d_gn = _norm(t.get("species")), _norm(t.get("genus"))
    k_sp, k_gn = _norm(row.get("species_name")), _norm(row.get("genus_name"))
    if d_sp and k_sp and d_sp == k_sp:
        return "species match"
    if d_gn and k_gn and d_gn == k_gn:
        return "genus match"
    if _norm(t.get("category")) != "virus":
        return "non-viral hit (%s)" % (t.get("category") or "?")
    return "different taxon"


def annotate_summary(detection_summary_tsv, contig_taxonomy_gz, out_path):
    tax = {}
    with gzip.open(contig_taxonomy_gz, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            parts = r["contig_id"].split("|")
            if len(parts) >= 2:
                tax[(r["sample_id"], parts[1])] = r
    with open(detection_summary_tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = (list(rows[0].keys()) if rows else []) + ANNOT_COLS
    n_hit = 0
    for row in rows:
        t = tax.get((row["sample_id"], row["taxid"]))
        if t is not None:
            n_hit += 1
        g = t or {}
        parts = g.get("contig_id", "").split("|")
        row.update({
            "blast_block_bp": g.get("contig_length", ""),
            "blast_region": parts[3] if len(parts) > 3 else "",
            "blast_n_hits": g.get("n_hits", ""),
            "blast_n_hits_lca": g.get("n_hits_lca", ""),
            "blast_best_subject": g.get("best_hit", ""),
            "blast_best_subject_taxon": g.get("best_hit_taxon", ""),
            "blast_best_subject_category": g.get("best_hit_category", ""),
            "blast_pct_identity": g.get("best_pident", ""),
            "blast_bitscore": g.get("best_bitscore", ""),
            "diamond_taxid": g.get("taxid", ""),
            "diamond_taxon": g.get("taxon_name", ""),
            "diamond_rank": g.get("rank", ""),
            "diamond_category": g.get("category", ""),
            "diamond_species": g.get("species", ""),
            "diamond_genus": g.get("genus", ""),
            "diamond_family": g.get("family", ""),
            "blast_agreement": agreement(row, t),
        })
    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    log.info("%d of %d summary rows carry a DIAMOND result", n_hit, len(rows))
    return out_path
