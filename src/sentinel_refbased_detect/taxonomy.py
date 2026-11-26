"""Task resolve_taxid_aliases.

The reference table carries current NCBI taxids, but a Kraken2 database is
pinned to the taxonomy of its build date and reports the legacy taxid for
anything renumbered since (Lassa 11620 -> 3052310). Resolve those aliases, the
ancestor-to-reference links, and the species/genus levels of every reference.
"""

import csv
import re
import shutil
import tarfile
from pathlib import Path

from .utils import PipelineError, download, log

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def get_taxdump(taxdump_tar_gz, out_dir, url=TAXDUMP_URL):
    """Return a local taxdump.tar.gz, downloading the current one if needed."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / "taxdump.tar.gz"
    if taxdump_tar_gz:
        src = Path(taxdump_tar_gz)
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
    elif not dest.exists():
        tmp = out_dir / "taxdump.tar.gz.part"
        download(url, tmp)
        tmp.rename(dest)
    return dest


def extract_dmp(taxdump_tar_gz, out_dir, members=("merged.dmp", "names.dmp", "nodes.dmp")):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if all((out_dir / m).exists() for m in members):
        return out_dir
    with tarfile.open(taxdump_tar_gz) as tf:
        names = set(tf.getnames())
        for m in members:
            cand = m if m in names else next((n for n in names if n.endswith("/" + m)), None)
            if cand is None:
                if m == "merged.dmp":
                    (out_dir / m).write_text("")
                    continue
                raise PipelineError("%s not found in %s" % (m, taxdump_tar_gz))
            with tf.extractfile(cand) as src, open(out_dir / m, "wb") as dst:
                shutil.copyfileobj(src, dst)
    return out_dir


def build_combined_refs(ref_genomes_tsv, extra_refs_tsv, out_path):
    shutil.copyfile(ref_genomes_tsv, out_path)
    with open(out_path, "rb") as fh:
        data = fh.read()
    with open(out_path, "ab") as fh:
        # Refgenomes tables can end without a newline; appending blind would
        # weld the discovered row onto the last row.
        if data and not data.endswith(b"\n"):
            fh.write(b"\n")
        if extra_refs_tsv and Path(extra_refs_tsv).exists() and Path(extra_refs_tsv).stat().st_size:
            fh.write(Path(extra_refs_tsv).read_bytes())
    return out_path


def resolve_taxid_aliases(ref_genomes_tsv, out_dir, dmp_dir, extra_refs_tsv=None):
    """Write taxid_aliases.tsv, name_aliases.tsv, taxon_levels.tsv and
    combined_refs.tsv into out_dir and return their paths as a dict."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dmp_dir = Path(dmp_dir)
    combined = build_combined_refs(ref_genomes_tsv, extra_refs_tsv, out_dir / "combined_refs.tsv")

    old_to_new = {}
    with open(dmp_dir / "merged.dmp") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                old_to_new[parts[0]] = parts[1]

    new_to_old = {}
    for old, new in old_to_new.items():
        new_to_old.setdefault(new, []).append(old)

    table_taxids = set()
    priority_taxids = set()
    table_rows_per_taxid = {}
    with open(combined) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and row[0].strip().isdigit():
                t = row[0].strip()
                table_taxids.add(t)
                table_rows_per_taxid[t] = table_rows_per_taxid.get(t, 0) + 1
                if not (len(row) > 4 and row[4].strip() == "discovered"):
                    priority_taxids.add(t)

    # A discovered virus is reachable only by the taxid Kraken2 already named;
    # it takes no part in the alias, synonym or ancestor machinery.
    rows = []
    for taxid in sorted(table_taxids, key=int):
        rows.append((taxid, taxid, "self"))
        if taxid not in priority_taxids:
            continue
        for old in sorted(new_to_old.get(taxid, []), key=int):
            rows.append((old, taxid, "legacy"))
        if taxid in old_to_new:
            rows.append((old_to_new[taxid], taxid, "successor"))

    parent, rank_of = {}, {}
    with open(dmp_dir / "nodes.dmp") as fh:
        for line in fh:
            p = [x.strip() for x in line.split("|")]
            if len(p) >= 2 and p[0].isdigit():
                parent[p[0]] = p[1]
                if len(p) >= 3:
                    rank_of[p[0]] = p[2]

    def lineage(taxid):
        out, seen = [], set()
        while taxid in parent and taxid not in seen:
            seen.add(taxid)
            out.append(taxid)
            if parent[taxid] == taxid:
                break
            taxid = parent[taxid]
        return out

    # How many reference rows each ancestor would pull in.
    rows_under = {}
    for taxid in priority_taxids:
        n_rows = table_rows_per_taxid.get(taxid, 1)
        for anc in lineage(taxid)[1:]:
            rows_under.setdefault(anc, {})[taxid] = n_rows

    for anc in sorted(rows_under, key=int):
        total = sum(rows_under[anc].values())
        for taxid in sorted(rows_under[anc], key=int):
            rows.append((anc, taxid, "descendant:%d" % total))

    def legacy_ids(node):
        return [node] + sorted(new_to_old.get(node, []), key=int)

    levels = {}
    need_sci = set()
    for taxid in sorted(table_taxids, key=int):
        path = lineage(taxid)
        i_sp = next((i for i, a in enumerate(path) if rank_of.get(a) == "species"), None)
        i_gn = next((i for i, a in enumerate(path) if rank_of.get(a) == "genus"), None)

        def chain(idx, lo, path=path):
            if idx is None:
                return "", []
            node = path[idx]
            out = legacy_ids(node)
            for n in reversed(path[lo:idx]):
                if rank_of.get(n, "no rank") == "no rank":
                    out.extend(legacy_ids(n))
            return node, out

        sp, sp_chain = chain(i_sp, 1)
        gn, gn_chain = chain(i_gn, 0 if i_sp is None else i_sp + 1)
        levels[taxid] = (sp, sp_chain, gn, gn_chain, taxid)
        need_sci.update(x for x in (sp, gn) if x)

    aliases_path = out_dir / "taxid_aliases.tsv"
    with open(aliases_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["report_taxid", "table_taxid", "relation"])
        w.writerows(rows)

    n_alias = sum(1 for r in rows if r[2] in ("legacy", "successor"))
    n_desc = sum(1 for r in rows if r[2].startswith("descendant"))
    log.info("%d table taxids, %d alias taxids, %d ancestor-to-reference links",
             len(table_taxids), n_alias, n_desc)

    # Only the taxon's own identities, never an ancestor's.
    want = {r[0]: r[1] for r in rows
            if r[2] in ("self", "legacy", "successor") and r[1] in priority_taxids}

    def norm(s):
        return re.sub(r"[^0-9a-z]+", " ", s.lower()).strip()

    names = {}
    sci_name = {}
    with open(dmp_dir / "names.dmp", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = [x.strip() for x in line.split("|")]
            if len(p) < 4:
                continue
            if p[0] in need_sci and p[3] == "scientific name":
                sci_name[p[0]] = p[1]
            if p[0] not in want:
                continue
            nm = norm(p[1])
            if len(nm) >= 4:
                names.setdefault(nm, set()).add(want[p[0]])

    levels_path = out_dir / "taxon_levels.tsv"
    with open(levels_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["table_taxid", "species_taxid", "species_name",
                    "species_lookup", "genus_taxid", "genus_name",
                    "genus_lookup", "fallback_taxid"])
        for taxid in sorted(levels, key=int):
            sp, sp_chain, gn, gn_chain, fb = levels[taxid]
            w.writerow([taxid, sp, sci_name.get(sp, ""), ",".join(sp_chain),
                        gn, sci_name.get(gn, ""), ",".join(gn_chain), fb])

    ambiguous = sorted(n for n, t in names.items() if len(t) > 1)
    name_path = out_dir / "name_aliases.tsv"
    with open(name_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["normalized_name", "table_taxid"])
        for nm in sorted(names):
            if len(names[nm]) == 1:
                w.writerow([nm, next(iter(names[nm]))])
    log.info("%d unambiguous names, %d dropped as ambiguous",
             len(names) - len(ambiguous), len(ambiguous))
    for nm in ambiguous:
        log.debug("  ambiguous, not used: %s -> %s", nm, sorted(names[nm]))

    return {"taxid_aliases": aliases_path, "name_aliases": name_path,
            "taxon_levels": levels_path, "combined_refs": combined}
