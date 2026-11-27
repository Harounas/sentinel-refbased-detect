"""Build the decoy panel the specificity criterion maps against.

Two kinds of row, both generated rather than hand-curated, so that 288 priority
taxa do not become 288 lists for somebody to maintain:

  sibling  one RefSeq genome per OTHER species of the same genus. These are the
           genomes a read from a neighbouring virus would really prefer, and
           they are what makes "can this read be attributed to this taxon?" a
           question with an answer.
  strain   genomes of the priority species itself, clustered so that near
           identical deposits collapse to one representative. Without them a
           true read looks ambiguous merely because the panel holds one
           arbitrary lineage; with hundreds of them the evidence for one taxon
           gets split across copies of the same sequence.

The result is a frozen, versioned TSV that ships with the package. Regenerate
it deliberately, read the diff, bump the version.
"""

import csv
import re
from pathlib import Path

from .cluster import cluster_sequences
from .taxonomy import extract_dmp, get_taxdump
from .utils import log

PANEL_COLS = ["taxid", "code", "name", "accessions", "role", "related_to", "note"]
SEGMENT = re.compile(r"segment\s+([A-Za-z0-9]+)", re.I)


def load_panel(path):
    """Panel rows grouped by the priority taxid each one guards."""
    by_taxid = {}
    if not path or not Path(path).exists():
        return by_taxid
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for t in (row.get("related_to") or "").split(","):
                t = t.strip()
                if t:
                    by_taxid.setdefault(t, []).append(row)
    return by_taxid


def _taxonomy(dmp_dir):
    parent, rank, name = {}, {}, {}
    with open(Path(dmp_dir) / "nodes.dmp") as fh:
        for line in fh:
            p = [x.strip() for x in line.split("|")]
            if len(p) >= 3 and p[0].isdigit():
                parent[p[0]], rank[p[0]] = p[1], p[2]
    with open(Path(dmp_dir) / "names.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = [x.strip() for x in line.split("|")]
            if len(p) >= 4 and p[3] == "scientific name":
                name[p[0]] = p[1]
    return parent, rank, name


def _lineage(taxid, parent):
    out, seen = [], set()
    while taxid in parent and taxid not in seen:
        seen.add(taxid)
        out.append(taxid)
        if parent[taxid] == taxid:
            break
        taxid = parent[taxid]
    return out


def _ancestor_at(taxid, parent, rank, want):
    for a in _lineage(taxid, parent):
        if rank.get(a) == want:
            return a
    return None


def _records(entrez, term, retmax):
    res = entrez.json("esearch.fcgi", db="nuccore", term=term, retmax=retmax)
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    out = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        summ = entrez.json("esummary.fcgi", db="nuccore", id=",".join(chunk))
        for uid in chunk:
            d = summ.get("result", {}).get(uid)
            if not d:
                continue
            acc = d.get("accessionversion") or d.get("caption", "")
            if acc:
                out.append({"acc": acc, "taxid": str(d.get("taxid", "")),
                            "slen": int(d.get("slen") or 0), "title": d.get("title", "")})
    return out


def _pick_genome(recs):
    """One genome per taxon: every segment if segmented, else the longest."""
    seg = {}
    for r in recs:
        m = SEGMENT.search(r["title"])
        if m:
            key = m.group(1).upper()
            if key not in seg or r["slen"] > seg[key]["slen"]:
                seg[key] = r
    if len(seg) >= 2:
        return [seg[k]["acc"] for k in sorted(seg)], "%d segments" % len(seg)
    best = max(recs, key=lambda r: r["slen"])
    return [best["acc"]], "%d bp" % best["slen"]


def _code(name, fallback):
    tag = re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")[:24]
    return tag or fallback


def build_siblings(ref_tsv, entrez, dmp_dir, max_per_genus=200):
    """One RefSeq genome per other species of each priority taxon's genus."""
    parent, rank, name = _taxonomy(dmp_dir)

    priority = []
    with open(ref_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].strip().isdigit():
                if not (len(row) > 4 and row[4].strip() == "discovered"):
                    priority.append(row[0].strip())

    # Genus -> the priority taxa that sit under it. One query serves them all.
    genera, no_genus = {}, []
    for taxid in sorted(set(priority), key=int):
        genus = _ancestor_at(taxid, parent, rank, "genus")
        if genus:
            genera.setdefault(genus, []).append(taxid)
        else:
            # Fall back to family, where a genus is missing or the taxon sits
            # above one. Rare, and better than guarding nothing.
            fam = _ancestor_at(taxid, parent, rank, "family")
            (genera.setdefault(fam, []).append(taxid) if fam else no_genus.append(taxid))
    if no_genus:
        log.info("%d priority taxa have neither genus nor family: %s",
                 len(no_genus), ", ".join(no_genus[:10]))

    rows = []
    for n, (group, guarded) in enumerate(sorted(genera.items(), key=lambda kv: int(kv[0])), 1):
        gname = name.get(group, group)
        log.info("[%d/%d] %s (%s), guarding %d priority taxa",
                 n, len(genera), gname, group, len(guarded))
        term = ("txid%s[Organism:exp] AND srcdb_refseq[PROP] AND 500:3000000[SLEN]" % group)
        try:
            recs = _records(entrez, term, max_per_genus)
        except Exception as exc:  # noqa: BLE001
            log.warning("  lookup failed: %s", exc)
            continue

        # A record belongs to a priority taxon when that taxon is on its
        # lineage; everything else in the genus is a sibling worth offering.
        by_species = {}
        for r in recs:
            if not r["taxid"]:
                continue
            lin = set(_lineage(r["taxid"], parent))
            if lin & set(guarded):
                continue
            sp = _ancestor_at(r["taxid"], parent, rank, "species") or r["taxid"]
            by_species.setdefault(sp, []).append(r)

        for sp, sp_recs in sorted(by_species.items(), key=lambda kv: int(kv[0])):
            accs, how = _pick_genome(sp_recs)
            sp_name = name.get(sp, sp_recs[0]["title"][:60])
            rows.append({"taxid": sp, "code": _code(sp_name, "taxid" + sp), "name": sp_name,
                         "accessions": ":".join(accs), "role": "sibling",
                         "related_to": ",".join(sorted(guarded, key=int)),
                         "note": "%s of %s; %s" % (rank.get(sp, "species"), gname, how)})
        log.info("  %d sibling species", len(by_species))
    return rows


def build_strains(ref_tsv, entrez, dmp_dir, cache_dir, min_identity=0.95,
                  max_candidates=300, max_representatives=12):
    """Strain genomes of each priority species, clustered to representatives."""
    parent, rank, name = _taxonomy(dmp_dir)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    priority = []
    with open(ref_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].strip().isdigit():
                if not (len(row) > 4 and row[4].strip() == "discovered"):
                    priority.append((row[0].strip(), (row[1] or "").strip()))

    # One job per species: strains of a strain-level row belong to its species.
    jobs = {}
    for taxid, code in priority:
        sp = _ancestor_at(taxid, parent, rank, "species") or taxid
        jobs.setdefault(sp, {"code": code or "NA", "guards": []})["guards"].append(taxid)

    rows = []
    for n, (sp, job) in enumerate(sorted(jobs.items(), key=lambda kv: int(kv[0])), 1):
        sp_name = name.get(sp, sp)
        log.info("[%d/%d] %s (%s)", n, len(jobs), sp_name, sp)
        term = ("txid%s[Organism:exp] AND complete genome[Title] "
                "AND 500:3000000[SLEN]" % sp)
        try:
            recs = _records(entrez, term, max_candidates)
        except Exception as exc:  # noqa: BLE001
            log.warning("  lookup failed: %s", exc)
            continue
        # Segmented genomes are not clustered: a segment is not a genome, and
        # comparing one segment against another says nothing useful.
        if any(SEGMENT.search(r["title"]) for r in recs):
            log.info("  %d records look segmented; keeping the reference genome only", len(recs))
            continue
        recs = [r for r in recs if r["slen"] > 0][:max_candidates]
        if len(recs) < 2:
            log.info("  %d genome(s); nothing to cluster", len(recs))
            continue

        fasta = cache_dir / ("%s.fasta" % sp)
        if not fasta.exists():
            text = entrez.efetch_fasta([r["acc"] for r in recs])
            fasta.write_text(text if text.endswith("\n") else text + "\n")
        seqs = {}
        label, buf = None, []
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                if label:
                    seqs[label] = "".join(buf)
                label, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
        if label:
            seqs[label] = "".join(buf)

        clusters = cluster_sequences(seqs.items(), min_identity=min_identity)
        log.info("  %d genomes -> %d cluster(s) at %.0f%% identity",
                 len(seqs), len(clusters), 100 * min_identity)
        for c in clusters[:max_representatives]:
            rows.append({
                "taxid": sp, "code": _code(sp_name, job["code"]), "name": sp_name,
                "accessions": c["representative"], "role": "strain",
                "related_to": ",".join(sorted(job["guards"], key=int)),
                "note": "%d bp; represents %d genome(s) at %.0f%% identity"
                        % (c["length"], len(c["members"]), 100 * min_identity)})
    return rows


def write_panel(rows, out_path):
    seen, unique = set(), []
    for r in rows:
        key = (r["accessions"], r["role"])
        if key not in seen:
            seen.add(key)
            unique.append(r)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PANEL_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(unique)
    log.info("wrote %d panel rows to %s", len(unique), out_path)
    return out_path


def build(ref_tsv, out_path, entrez, cache_dir, kinds=("sibling",), taxdump=None,
          min_identity=0.95, max_candidates=300):
    dmp = extract_dmp(get_taxdump(taxdump, Path(cache_dir) / "taxonomy"),
                      Path(cache_dir) / "taxonomy" / "dmp")
    rows = []
    if "sibling" in kinds:
        rows += build_siblings(ref_tsv, entrez, dmp)
    if "strain" in kinds:
        rows += build_strains(ref_tsv, entrez, dmp, Path(cache_dir) / "strains",
                              min_identity=min_identity, max_candidates=max_candidates)
    return write_panel(rows, out_path)
