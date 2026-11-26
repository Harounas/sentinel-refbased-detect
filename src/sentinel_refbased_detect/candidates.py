"""Task select_candidates: which reference rows to align for one sample."""

import csv
import os
import re

from .utils import PipelineError, log

ROUTE_RANK = {"taxid": 4, "taxid_alias": 3, "name_synonym": 2, "regex": 1}
RANK_ORDER = "URDKPCOFGS"

CANDIDATE_COLS = [
    "sample_id", "taxid", "report_taxid", "report_name", "report_rank", "match_basis",
    "virus_code", "strain_name", "virus_importance", "ref_index",
    "accessions", "n_segments", "kraken_reads", "kraken_pct",
    "species_taxid", "species_name", "species_direct_reads",
    "species_clade_reads", "species_report_taxid",
    "genus_taxid", "genus_name", "genus_direct_reads",
    "genus_clade_reads", "genus_report_taxid"]


def norm(s):
    return re.sub(r"[^0-9a-z]+", " ", s.lower()).strip()


def specificity(basis):
    return ROUTE_RANK.get(basis.split(":")[0], 0)


def rank_depth(code):
    """Kraken2 rank codes, shallowest first; S1/D2 suffixes count as tenths."""
    if not code:
        return -1.0
    base = RANK_ORDER.find(code[0])
    if base < 0:
        return -1.0
    try:
        return base + int(code[1:] or 0) / 10.0
    except ValueError:
        return float(base)


def load_refs(ref_genomes_tsv):
    """Priority table: taxid, code, name, colon-joined accessions[, importance]."""
    refs = {}
    with open(ref_genomes_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 4 or not row[0].strip().isdigit():
                continue
            taxid = row[0].strip()
            code = (row[1] or "").strip() or "NA"
            name = (row[2] or "").strip()
            accs = [a.strip() for a in row[3].split(":") if a.strip()]
            importance = row[4].strip() if len(row) > 4 and row[4].strip() else "priority"
            if accs:
                refs.setdefault(taxid, []).append((code, name, accs, importance))
    return refs


def load_regex_rules(path):
    rules = []
    if not path or not os.path.exists(path):
        return rules
    with open(path) as fh:
        for raw in fh:
            raw = raw.rstrip("\n")
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            f = raw.split("\t")
            if len(f) < 2 or not f[1].strip().isdigit():
                continue
            if f[0].strip().lower() in ("pattern", "regex"):
                continue
            try:
                rules.append((re.compile(f[0].strip(), re.I), f[1].strip()))
            except re.error as exc:
                raise PipelineError("bad regex in name_regex_tsv: %r: %s" % (f[0], exc))
    return rules


def select_candidates(sample_id, kraken_report, ref_genomes_tsv, taxid_aliases_tsv,
                      name_aliases_tsv, taxon_levels_tsv, out_path, params,
                      name_regex_tsv=None):
    refs = load_refs(ref_genomes_tsv)

    alias, descendants = {}, {}
    with open(taxid_aliases_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rel = row["relation"]
            if rel.startswith("descendant"):
                descendants.setdefault(row["report_taxid"], []).append(
                    (row["table_taxid"], int(rel.split(":")[1])))
            else:
                alias[row["report_taxid"]] = row["table_taxid"]

    name_alias = {}
    with open(name_aliases_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name_alias[row["normalized_name"]] = row["table_taxid"]

    levels = {}
    with open(taxon_levels_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            levels[row["table_taxid"]] = row

    regex_rules = load_regex_rules(name_regex_tsv)

    min_rank_depth = rank_depth(params.min_match_rank)
    if min_rank_depth < 0:
        raise PipelineError("min_match_rank %r is not a Kraken2 rank code" % params.min_match_rank)

    hits, report_index = {}, {}
    with open(kraken_report) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            if len(f) >= 8 and f[6].strip().isdigit():
                taxid, rank = f[6].strip(), f[5].strip()
            elif f[4].strip().isdigit():
                taxid, rank = f[4].strip(), f[3].strip()
            else:
                continue
            try:
                pct = float(f[0])
                clade_reads = int(f[1])
                direct_reads = int(f[2])
            except ValueError:
                continue
            report_name = f[-1].strip()

            report_index[taxid] = {"rank": rank, "direct": direct_reads,
                                   "clade": clade_reads, "pct": pct, "name": report_name}

            this_depth = rank_depth(rank)
            if 0 <= this_depth < min_rank_depth:
                continue

            found = []
            if taxid in refs:
                found = [(taxid, "taxid")]
            elif taxid in alias and alias[taxid] in refs:
                found = [(alias[taxid], "taxid_alias")]
            else:
                hit = name_alias.get(norm(report_name))
                if hit in refs:
                    found = [(hit, "name_synonym")]
                else:
                    for pattern, target in regex_rules:
                        if target in refs and pattern.search(report_name):
                            found = [(target, "regex:" + pattern.pattern)]
                            break

            if params.expand_descendants and taxid in descendants:
                below_species = rank.startswith("S")
                for table_taxid, n_rows in descendants[taxid]:
                    if table_taxid not in refs:
                        continue
                    if not below_species and n_rows > params.max_rows_above_species:
                        continue
                    if any(t == table_taxid for t, _ in found):
                        continue
                    found.append((table_taxid, "descendant_of:%s" % taxid))

            # Most specific route first, then deepest rank, then read count.
            for key, basis in found:
                prev = hits.get(key)
                cur_key = (specificity(basis), rank_depth(rank), clade_reads)
                if prev is None or cur_key > (specificity(prev[3]), rank_depth(prev[5]), prev[1]):
                    hits[key] = (pct, clade_reads, taxid, basis, report_name, rank)

    def level_counts(node, lookup_csv, fallback):
        ids = [x for x in lookup_csv.split(",") if x] or ([node] if node else [])
        present = [(i, report_index[i]) for i in ids if i in report_index]
        if present:
            direct = sum(e["direct"] for _, e in present)
            return str(direct), str(present[0][1]["clade"]), ",".join(i for i, _ in present)
        e = report_index.get(fallback)
        if e is not None:
            return str(e["direct"]), str(e["clade"]), fallback
        return "", "", ""

    rows = []
    for taxid, (pct, clade_reads, report_taxid, basis, report_name, report_rank) in sorted(
            hits.items(), key=lambda kv: -kv[1][1]):
        discovered = refs[taxid][0][3] == "discovered"
        floor = params.min_kraken_reads_discovered if discovered else params.min_kraken_reads
        if clade_reads < floor or pct < params.min_kraken_pct:
            continue
        lv = levels.get(taxid, {})
        sp_taxid, gn_taxid = lv.get("species_taxid", ""), lv.get("genus_taxid", "")
        fb = lv.get("fallback_taxid", "")
        sp_direct, sp_clade, sp_seen = level_counts(sp_taxid, lv.get("species_lookup", ""), fb)
        gn_direct, gn_clade, gn_seen = level_counts(gn_taxid, lv.get("genus_lookup", ""), "")
        for idx, (code, name, accs, importance) in enumerate(refs[taxid][:params.max_refs_per_virus]):
            rows.append({
                "sample_id": sample_id, "taxid": taxid, "report_taxid": report_taxid,
                "report_name": report_name, "report_rank": report_rank,
                "match_basis": basis, "virus_code": code, "strain_name": name,
                "virus_importance": importance, "ref_index": str(idx),
                "accessions": ":".join(accs), "n_segments": str(len(accs)),
                "kraken_reads": str(clade_reads), "kraken_pct": f"{pct:.4f}",
                "species_taxid": sp_taxid, "species_name": lv.get("species_name", ""),
                "species_direct_reads": sp_direct, "species_clade_reads": sp_clade,
                "species_report_taxid": sp_seen,
                "genus_taxid": gn_taxid, "genus_name": lv.get("genus_name", ""),
                "genus_direct_reads": gn_direct, "genus_clade_reads": gn_clade,
                "genus_report_taxid": gn_seen,
            })

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=CANDIDATE_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    basis_counts = {}
    for r in rows:
        b = r["match_basis"].split(":")[0]
        basis_counts[b] = basis_counts.get(b, 0) + 1
    log.info("%s: %d candidate references across %d viruses; matched by %s",
             sample_id, len(rows), len({r["taxid"] for r in rows}), basis_counts or "nothing")
    return out_path
