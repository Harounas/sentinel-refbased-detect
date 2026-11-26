"""Tasks summarize_sample and merge_results."""

import csv
import os
import shutil
from pathlib import Path

from .utils import log

LEVEL_COLS = ["species_taxid", "species_name", "species_direct_reads",
              "species_clade_reads", "species_report_taxid",
              "genus_taxid", "genus_name", "genus_direct_reads",
              "genus_clade_reads", "genus_report_taxid"]

SUMMARY_COLS = ["sample_id", "taxid", "report_taxid", "report_name", "report_rank", "match_basis",
                "virus_code", "strain_name", "virus_importance", "best_reference",
                "accessions", "n_regions", "kraken_reads", "kraken_pct"] + LEVEL_COLS + \
               ["mapq_by_region", "unique_reads_by_region", "covered_bp_by_region",
                "total_unique_reads_mapped", "n_duplicate_reads",
                "duplicate_fraction", "breadth_of_coverage",
                "paired_read_concordance", "longest_covered_block_bp",
                "n_covered_blocks", "n_qualifying_regions",
                "detection", "detection_basis",
                "consensus_status", "n_references_compared"]

REGION_COLS = ["sample_id", "taxid", "virus_code", "virus_importance",
               "best_reference", "region",
               "region_length", "mean_mapq", "unique_reads",
               "covered_bp_min_depth", "longest_block_bp", "n_blocks",
               "breadth_of_coverage", "total_unique_reads_mapped",
               "paired_read_concordance", "detection"]


def detection_rule(region_rows, params):
    """Return (detected, reason, n_qualifying_regions).

    detected if one region holds an unbroken run >= min_region_bp at min_depth,
    or if >= min_qualifying_regions regions each hold a block >= min_block_bp.
    """
    n_qual = sum(1 for r in region_rows if r["n_blocks"] >= 1)
    best_block = max((r["longest_block_bp"] for r in region_rows), default=0)
    best_region = max((r["covered_bp_min_depth"] for r in region_rows), default=0)
    if best_block >= params.min_region_bp:
        return True, (f"contiguous {best_block} bp >= {params.min_region_bp} bp "
                      f"at {params.min_depth}x"), n_qual
    if n_qual >= params.min_qualifying_regions:
        return True, (f"{n_qual} regions each with a block >= "
                      f"{params.min_block_bp} bp at {params.min_depth}x"), n_qual
    return False, (f"longest contiguous {best_block} bp ({best_region} bp covered "
                   f"in all), {n_qual} region(s) with a block >= "
                   f"{params.min_block_bp} bp"), n_qual


def summarize_sample(sample_id, align_dirs, candidates_tsv, out_dir, params):
    out_dir = Path(out_dir)
    cons_dir = out_dir / "consensus"
    cons_dir.mkdir(parents=True, exist_ok=True)

    blocks, aln, by_ref, consensus_paths = {}, {}, {}, []
    for d in align_dirs:
        d = Path(d)
        for path in d.glob("*.block.fasta"):
            if path.stat().st_size == 0:
                continue
            head, seq = "", []
            with open(path) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        head = line[1:]
                    elif line:
                        seq.append(line)
            parts = head.split("|")
            if head and len(parts) >= 4 and seq:
                blocks[parts[0]] = {"region": parts[1], "coords": parts[2],
                                    "length": parts[3], "seq": "".join(seq)}
        consensus_paths.extend(d.glob("*.consensus.fasta"))
        p = d / "aln_stats.tsv"
        if p.exists():
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 7:
                        aln[f[0]] = {"taxid": f[1], "virus_code": f[2], "ref_index": f[3],
                                     "n_mapped": int(f[4]), "n_proper": int(f[5]),
                                     "n_duplicate": int(f[6]),
                                     "consensus_status": f[7] if len(f) > 7 else ""}
        p = d / "region_stats.tsv"
        if p.exists():
            with open(p) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) != 8:
                        continue
                    by_ref.setdefault(f[0], []).append({
                        "region": f[1], "region_length": int(f[2]),
                        "unique_reads": int(f[3]), "mean_mapq": float(f[4]),
                        "covered_bp_min_depth": int(f[5]),
                        "longest_block_bp": int(f[6]), "n_blocks": int(f[7])})

    meta = {}
    with open(candidates_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[(row["taxid"], row["ref_index"])] = row

    records = []
    for tag, region_rows in by_ref.items():
        st = aln.get(tag)
        if st is None:
            continue
        total_len = sum(r["region_length"] for r in region_rows)
        total_cov = sum(r["covered_bp_min_depth"] for r in region_rows)
        detected, reason, n_qual = detection_rule(region_rows, params)
        cand = meta.get((st["taxid"], st["ref_index"]), {})
        records.append({
            "ref_tag": tag, "taxid": st["taxid"], "virus_code": st["virus_code"],
            "region_rows": region_rows,
            "total_unique": sum(r["unique_reads"] for r in region_rows),
            "n_mapped": st["n_mapped"], "n_duplicate": st["n_duplicate"],
            "breadth": (total_cov / total_len) if total_len else 0.0,
            "concordance": (st["n_proper"] / st["n_mapped"]) if st["n_mapped"] else 0.0,
            "detected": detected, "reason": reason, "total_cov": total_cov,
            "n_blocks": sum(r["n_blocks"] for r in region_rows),
            "n_qual_regions": n_qual,
            "consensus_status": st.get("consensus_status", ""),
            "strain_name": cand.get("strain_name", ""),
            "virus_importance": cand.get("virus_importance", "priority"),
            "accessions": cand.get("accessions", ""),
            "report_taxid": cand.get("report_taxid", ""),
            "report_name": cand.get("report_name", ""),
            "report_rank": cand.get("report_rank", ""),
            "match_basis": cand.get("match_basis", ""),
            "kraken_reads": cand.get("kraken_reads", ""),
            "kraken_pct": cand.get("kraken_pct", ""),
            "levels": {c: cand.get(c, "") for c in LEVEL_COLS},
        })

    # Best reference per virus: most bases covered at min_depth, then most reads.
    best = {}
    for rec in records:
        key = (rec["total_cov"], rec["total_unique"])
        if rec["taxid"] not in best or key > best[rec["taxid"]][0]:
            best[rec["taxid"]] = (key, rec)

    summary_rows, region_out, block_records = [], [], []
    for taxid, (_, rec) in sorted(best.items()):
        rows = rec["region_rows"]
        detection = "detected" if rec["detected"] else "not detected"
        denom = rec["n_mapped"] + rec["n_duplicate"]
        summary_rows.append({
            "sample_id": sample_id, "taxid": taxid,
            "report_taxid": rec["report_taxid"], "report_name": rec["report_name"],
            "report_rank": rec["report_rank"], "match_basis": rec["match_basis"],
            "virus_code": rec["virus_code"], "strain_name": rec["strain_name"],
            "virus_importance": rec["virus_importance"], "best_reference": rec["ref_tag"],
            "accessions": rec["accessions"], "n_regions": len(rows),
            "kraken_reads": rec["kraken_reads"], "kraken_pct": rec["kraken_pct"],
            **rec["levels"],
            "mapq_by_region": ";".join(f"{r['region']}:{r['mean_mapq']}" for r in rows),
            "unique_reads_by_region": ";".join(f"{r['region']}:{r['unique_reads']}" for r in rows),
            "covered_bp_by_region": ";".join(f"{r['region']}:{r['covered_bp_min_depth']}" for r in rows),
            "total_unique_reads_mapped": rec["total_unique"],
            "n_duplicate_reads": rec["n_duplicate"],
            "duplicate_fraction": round(rec["n_duplicate"] / denom, 6) if denom else "",
            "breadth_of_coverage": round(rec["breadth"], 6),
            "paired_read_concordance": round(rec["concordance"], 6),
            "longest_covered_block_bp": max((r["longest_block_bp"] for r in rows), default=0),
            "n_covered_blocks": sum(r["n_blocks"] for r in rows),
            "n_qualifying_regions": rec["n_qual_regions"],
            "detection": detection, "detection_basis": rec["reason"],
            "consensus_status": rec["consensus_status"],
            "n_references_compared": sum(1 for r in records if r["taxid"] == taxid),
        })
        for r in rows:
            region_out.append({
                "sample_id": sample_id, "taxid": taxid, "virus_code": rec["virus_code"],
                "virus_importance": rec["virus_importance"], "best_reference": rec["ref_tag"],
                "region": r["region"], "region_length": r["region_length"],
                "mean_mapq": r["mean_mapq"], "unique_reads": r["unique_reads"],
                "covered_bp_min_depth": r["covered_bp_min_depth"],
                "longest_block_bp": r["longest_block_bp"], "n_blocks": r["n_blocks"],
                "breadth_of_coverage": round(rec["breadth"], 6),
                "total_unique_reads_mapped": rec["total_unique"],
                "paired_read_concordance": round(rec["concordance"], 6),
                "detection": detection,
            })

        # Only detected viruses carry their longest block forward to DIAMOND.
        blk = blocks.get(rec["ref_tag"])
        if rec["detected"] and blk:
            qid = "|".join([sample_id, taxid, rec["ref_tag"], blk["region"], blk["coords"]])
            block_records.append((qid, blk))

        wanted = rec["ref_tag"] + ".consensus.fasta"
        src = next((p for p in consensus_paths if os.path.basename(p) == wanted), None)
        if src and os.path.getsize(src) > 0:
            shutil.copyfile(src, cons_dir / f"{sample_id}.{rec['ref_tag']}.consensus.fasta")

    paths = {"summary": out_dir / "summary.tsv", "by_region": out_dir / "by_region.tsv",
             "blocks_fasta": out_dir / "best_blocks.fasta", "block_map": out_dir / "block_map.tsv"}
    with open(paths["summary"], "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=SUMMARY_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    with open(paths["by_region"], "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=REGION_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(region_out)
    with open(paths["blocks_fasta"], "w") as out:
        for qid, blk in block_records:
            out.write(">%s\n" % qid)
            seq = blk["seq"]
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    with open(paths["block_map"], "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["query_id", "sample_id", "contig_id", "contig_length"])
        for qid, blk in block_records:
            w.writerow([qid, sample_id, qid, blk["length"]])

    log.info("%s: %d viruses, %d detected, %d blocks for DIAMOND", sample_id, len(summary_rows),
             sum(1 for r in summary_rows if r["detection"] == "detected"), len(block_records))
    return paths


def merge_tables(paths, out_path):
    """Concatenate TSVs that share a header; returns the number of data rows."""
    header_written, n = False, 0
    with open(out_path, "w") as out:
        for path in paths:
            if not path or not os.path.exists(path):
                continue
            with open(path) as fh:
                lines = fh.read().splitlines()
            if not lines:
                continue
            if not header_written:
                out.write(lines[0] + "\n")
                header_written = True
            for line in lines[1:]:
                if line.strip():
                    out.write(line + "\n")
                    n += 1
    return n
