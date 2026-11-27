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

from .cluster import cluster_sequences, pairwise_identity, select_diverse
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


def _priority_rows(ref_tsv):
    """(taxid, code) of every non-discovered row of the priority table."""
    out = []
    with open(ref_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].strip().isdigit():
                if not (len(row) > 4 and row[4].strip() == "discovered"):
                    out.append((row[0].strip(), (row[1] or "").strip()))
    return out


def _priority_taxids(ref_tsv):
    return sorted({t for t, _ in _priority_rows(ref_tsv)}, key=int)


def siblings_for_genus(group, guarded, entrez, parent, rank, name, max_per_genus=200):
    """One RefSeq genome per species of this genus that is not a guarded taxon."""
    gname = name.get(group, group)
    term = "txid%s[Organism:exp] AND srcdb_refseq[PROP] AND 500:3000000[SLEN]" % group
    try:
        recs = _records(entrez, term, max_per_genus)
    except Exception as exc:  # noqa: BLE001
        log.warning("  lookup failed for %s: %s", gname, exc)
        return []

    by_species = {}
    for r in recs:
        if not r["taxid"]:
            continue
        if set(_lineage(r["taxid"], parent)) & set(guarded):
            continue
        sp = _ancestor_at(r["taxid"], parent, rank, "species") or r["taxid"]
        by_species.setdefault(sp, []).append(r)

    rows = []
    for sp, sp_recs in sorted(by_species.items(), key=lambda kv: int(kv[0])):
        accs, how = _pick_genome(sp_recs)
        sp_name = name.get(sp, sp_recs[0]["title"][:60])
        rows.append({"taxid": sp, "code": _code(sp_name, "taxid" + sp), "name": sp_name,
                     "accessions": ":".join(accs), "role": "sibling",
                     "related_to": ",".join(sorted(guarded, key=int)),
                     "note": "%s of %s; %s" % (rank.get(sp, "species"), gname, how)})
    log.info("  %d sibling species", len(by_species))
    return rows


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
        log.info("[%d/%d] %s (%s), guarding %d priority taxa",
                 n, len(genera), name.get(group, group), group, len(guarded))
        rows += siblings_for_genus(group, guarded, entrez, parent, rank, name,
                                   max_per_genus=max_per_genus)
    return rows


def build_strains(ref_tsv, entrez, dmp_dir, cache_dir, min_identity=0.99,
                  max_candidates=300, max_representatives=12,
                  coverage_identity=0.90, min_members=2):
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
        jobs.setdefault(sp, {"code": code or "NA", "name": name.get(sp, sp),
                             "guards": []})["guards"].append(taxid)

    rows = []
    for n, (sp, job) in enumerate(sorted(jobs.items(), key=lambda kv: int(kv[0])), 1):
        sp_name = name.get(sp, sp)
        log.info("[%d/%d] %s (%s)", n, len(jobs), sp_name, sp)
        # "complete sequence" as well as "complete genome": a segmented virus
        # deposits one record per segment and titles it "segment L, complete
        # sequence", so the narrower phrase returns nothing at all for Lassa,
        # CCHF and influenza -- the divergent, segmented pathogens that most
        # need strain representatives.
        term = ("txid%s[Organism:exp] AND (complete genome[Title] OR "
                "complete sequence[Title]) AND 500:3000000[SLEN]" % sp)
        try:
            recs = _records(entrez, term, max_candidates)
        except Exception as exc:  # noqa: BLE001
            log.warning("  lookup failed: %s", exc)
            continue
        rows += strains_for_species(sp, job, recs, entrez, cache_dir,
                                    min_identity=min_identity,
                                    max_representatives=max_representatives,
                                    coverage_identity=coverage_identity,
                                    min_members=min_members)
    return rows


def strains_for_species(sp, job, recs, entrez, cache_dir, min_identity=0.99,
                        max_representatives=12, coverage_identity=0.90,
                        min_members=2):
    """Cluster one species' genomes into lineage representatives.

    Two things this has to get right, both of which a flat clustering gets
    wrong:

    Segments. Lassa, CCHF and influenza deposit one record per segment. An L
    segment and an S segment of the same virus share almost nothing, so
    clustering them together would call every genome unique and keep nothing
    useful. Each segment is therefore clustered on its own.

    Divergence. Lassa's clades sit far enough apart that they survive
    clustering as separate representatives, which is what should happen; the
    risk is the cap afterwards. Representatives are chosen to span the
    diversity, so a cap trims the crowded middle rather than deleting a clade.
    """
    sp_name = job.get("name") or sp
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    recs = [r for r in recs if r["slen"] > 0]
    if len(recs) < 2:
        log.info("  %d genome(s); nothing to cluster", len(recs))
        return []

    by_acc = {r["acc"]: r for r in recs}
    fasta = cache_dir / ("%s.fasta" % sp)
    if not fasta.exists():
        text = entrez.efetch_fasta([r["acc"] for r in recs])
        tmp = fasta.with_name(fasta.name + ".part")
        tmp.write_text(text if text.endswith("\n") else text + "\n")
        tmp.rename(fasta)
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

    # Group by segment where the records name one; everything else is a whole
    # genome and is clustered as a single group.
    groups = {}
    for acc, seq in seqs.items():
        title = by_acc.get(acc, {}).get("title", "")
        m = SEGMENT.search(title)
        groups.setdefault(m.group(1).upper() if m else "", []).append((acc, seq))

    out = []
    for seg, members in sorted(groups.items()):
        if len(members) < 2:
            continue
        clusters = cluster_sequences(members, min_identity=min_identity)
        lo, mid = pairwise_identity(clusters)
        kept, covered = select_diverse(clusters, max_representatives,
                                       min_members=min_members,
                                       coverage_identity=coverage_identity)
        label = ("segment %s" % seg) if seg else "genome"
        capped = len(kept) >= max_representatives and covered < coverage_identity
        log.info("  %s: %d sequences -> %d cluster(s) at %.0f%% dedup, keeping %d; "
                 "every genome within %.0f%% of a representative%s "
                 "(lineages %.0f%%-%.0f%% apart)",
                 label, len(members), len(clusters), 100 * min_identity, len(kept),
                 100 * covered,
                 " [capped: raise max_representatives]" if capped else "",
                 100 * lo, 100 * mid)
        for c in kept:
            out.append({
                "taxid": sp, "code": _code(sp_name, job.get("code") or "NA"),
                "name": sp_name, "accessions": c["representative"], "role": "strain",
                "related_to": ",".join(sorted(job["guards"], key=int)),
                "note": "%s; %d bp; represents %d of %d sequence(s) deduplicated at "
                        "%.0f%%; covers every genome within %.0f%%%s; "
                        "lineage spread %.0f%%-%.0f%%"
                        % (label, c["length"], len(c["members"]), len(members),
                           100 * min_identity, 100 * covered,
                           "; CAPPED" if capped else "", 100 * lo, 100 * mid)})
    return out


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


def list_jobs(ref_tsv, dmp_dir, kinds=("sibling", "strain")):
    """The shards a workflow scatters over: one genus, or one species, each.

    Grouping first matters for cost as well as correctness. 288 priority taxa
    collapse to about 50 genera, because many of them are siblings already, and
    a genus is one query however many of its species are on the list.
    """
    parent, rank, name = _taxonomy(dmp_dir)
    jobs = []
    if "sibling" in kinds:
        genera = {}
        for taxid in _priority_taxids(ref_tsv):
            group = (_ancestor_at(taxid, parent, rank, "genus")
                     or _ancestor_at(taxid, parent, rank, "family"))
            if group:
                genera.setdefault(group, []).append(taxid)
        for group, guarded in sorted(genera.items(), key=lambda kv: int(kv[0])):
            jobs.append({"kind": "sibling", "id": group, "name": name.get(group, group),
                         "guards": ",".join(sorted(set(guarded), key=int))})
    if "strain" in kinds:
        species = {}
        for taxid, code in _priority_rows(ref_tsv):
            sp = _ancestor_at(taxid, parent, rank, "species") or taxid
            species.setdefault(sp, {"code": code or "NA", "guards": []})["guards"].append(taxid)
        for sp, job in sorted(species.items(), key=lambda kv: int(kv[0])):
            jobs.append({"kind": "strain", "id": sp, "name": name.get(sp, sp),
                         "guards": ",".join(sorted(set(job["guards"]), key=int))})
    return jobs


def build_shard(kind, job_id, guards, ref_tsv, entrez, dmp_dir, cache_dir,
                min_identity=0.99, max_candidates=300, max_representatives=12,
                coverage_identity=0.90, min_members=2):
    """One shard of the panel: the siblings of a genus, or a species' strains."""
    parent, rank, name = _taxonomy(dmp_dir)
    guarded = [g for g in (guards or "").split(",") if g]
    if kind == "sibling":
        return siblings_for_genus(job_id, guarded, entrez, parent, rank, name)
    term = ("txid%s[Organism:exp] AND (complete genome[Title] OR "
            "complete sequence[Title]) AND 500:3000000[SLEN]" % job_id)
    recs = _records(entrez, term, max_candidates)
    job = {"code": "NA", "name": name.get(job_id, job_id), "guards": guarded or [job_id]}
    return strains_for_species(job_id, job, recs, entrez, cache_dir,
                               min_identity=min_identity,
                               max_representatives=max_representatives,
                               coverage_identity=coverage_identity,
                               min_members=min_members)


def build(ref_tsv, out_path, entrez, cache_dir, kinds=("sibling",), taxdump=None,
          min_identity=0.99, max_candidates=300, coverage_identity=0.90):
    dmp = extract_dmp(get_taxdump(taxdump, Path(cache_dir) / "taxonomy"),
                      Path(cache_dir) / "taxonomy" / "dmp")
    rows = []
    if "sibling" in kinds:
        rows += build_siblings(ref_tsv, entrez, dmp)
    if "strain" in kinds:
        rows += build_strains(ref_tsv, entrez, dmp, Path(cache_dir) / "strains",
                              min_identity=min_identity, max_candidates=max_candidates,
                              coverage_identity=coverage_identity)
    return write_panel(rows, out_path)
