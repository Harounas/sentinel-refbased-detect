"""Genus-level rescue: what the reads can and cannot say about species.

The coverage rule answers "is there enough sequence to call this virus
present?". When a priority pathogen fails it on a handful of reads, a second
question remains, and it is a taxonomic one: those reads came from something,
so what can be said about what?

The trap is to answer it with arithmetic. Two species of one genus share most
of their genome, so a read that fits both lands on whichever reference scores a
point higher, and summing those points names a species on no evidence at all.
This module refuses that. A read counts for a species only where it could not
equally have come from a relative, and when no read can tell the two apart the
honest answer is the genus.

Four outcomes, and three of them are not a species:

  species level        one species holds discriminating reads the others do not
  genus level          the reads fit the genus; nothing distinguishes a species
  possible mixed       two species each hold their own discriminating reads,
                       spread across the genome rather than piled in one spot
  insufficient         too little, or indistinguishable from the controls
"""

import csv
from pathlib import Path

from .utils import log

SPECIES_INFORMATIVE = "species-informative"
GENUS_INFORMATIVE = "genus-informative"
AMBIGUOUS = "ambiguous"

OUTCOME_SPECIES = "species-level identification"
OUTCOME_GENUS = "genus-level identification"
OUTCOME_MIXED = "possible mixed infection"
OUTCOME_INSUFFICIENT = "insufficient evidence"

REPORT_COLS = [
    # 1. primary detection
    "sample_id", "taxid", "virus_code", "strain_name", "primary_detection",
    "primary_detection_basis", "kraken_reads",
    # 2. rescue candidate and panel
    "rescue_triggered", "reference_panel_size", "panel_species",
    # 3. competitive mapping
    "reads_mapped", "reads_species_informative", "reads_genus_informative",
    "reads_ambiguous",
    # 4. species-specific counts
    "best_species", "best_species_informative_reads", "best_species_positions",
    "second_species", "second_species_informative_reads", "second_species_positions",
    "score_margin_median",
    # 5. alignment quality
    "mean_mapq", "pct_identity_closest_genome", "pct_reads_multimapping",
    # 6. coverage and distribution
    "breadth_of_coverage", "longest_block_bp", "n_bins_covered",
    "pct_reads_in_busiest_bin",
    # 7. competing species
    "competing_species_detail",
    # 8. controls
    "control_sample", "control_reads", "control_verdict",
    # 9-10. decision and why
    "taxonomic_outcome", "confidence", "taxonomic_call", "decision_reason",
]


def classify_reads(per_read, target_label, genus_of, margin=0):
    """Sort every read into species-informative, genus-informative or ambiguous.

    A read is species-informative when one species beats the rest by more than
    `margin`: no relative explains it as well. It is genus-informative when the
    species it cannot choose between all sit in one genus, which is real
    evidence about the genus and none about the species. Anything that spans
    genera, or that the background explains as well, is ambiguous.
    """
    out = {}
    for read, hits in per_read.items():
        if not hits:
            continue
        best = max(s for s, _, _ in hits.values())
        tied = {taxon for taxon, (s, _, _) in hits.items() if s >= best - margin}
        winner = max(hits, key=lambda t: hits[t][0])
        if len(tied) == 1:
            kind = SPECIES_INFORMATIVE
        else:
            genera = {genus_of(t) for t in tied}
            # An unknown genus cannot be vouched for, so it is not folded in.
            kind = (GENUS_INFORMATIVE
                    if len(genera) == 1 and None not in genera else AMBIGUOUS)
        out[read] = {"kind": kind, "winner": winner, "tied": tied,
                     "best": best, "hits": hits,
                     "second": sorted((s for s, _, _ in hits.values()), reverse=True)[1]
                              if len(hits) > 1 else None}
    return out


def species_evidence(classified, label):
    """Discriminating reads for one species, and how spread out they are.

    Positions matter as much as counts. Twenty reads stacked on one spot are
    one observation copied twenty times, most likely a conserved or repetitive
    stretch; five reads at five positions are five independent observations.
    """
    reads = [c for c in classified.values()
             if c["kind"] == SPECIES_INFORMATIVE and c["winner"] == label]
    positions = set()
    margins = []
    for c in reads:
        score, rname, pos = c["hits"][label]
        positions.add((rname, pos // 100))     # 100 bp, so near-duplicates merge
        if c["second"] is not None:
            margins.append(c["best"] - c["second"])
    margins.sort()
    return {"reads": len(reads), "positions": len(positions),
            "median_margin": margins[len(margins) // 2] if margins else 0}


def decide(target, best, second, counts, params, control_ok=True,
           target_label=None):
    """The decision table, in the order that protects against the usual error.

    Never the species with the most reads or a point more score: that is how a
    relative's reads get a name. Species is claimed only on discriminating
    reads that the runner-up does not have.
    """
    n_species = counts["species_informative"]
    n_genus = counts["genus_informative"]
    n_total = counts["mapped"]

    if not control_ok:
        return (OUTCOME_INSUFFICIENT, "low",
                "signal is not distinguishable from the negative controls")

    if n_total < params.rescue_min_reads:
        return (OUTCOME_INSUFFICIENT, "low",
                "%d mapped read(s), under rescue_min_reads=%d"
                % (n_total, params.rescue_min_reads))

    b_reads, b_pos = best["reads"], best["positions"]
    s_reads, s_pos = second["reads"], second["positions"]

    # Both species carry their own discriminating reads, each spread over the
    # genome. Neither explains the other's reads, so reporting one would mean
    # discarding evidence for the second.
    if (b_reads >= params.rescue_min_species_reads
            and s_reads >= params.rescue_min_species_reads
            and b_pos >= params.rescue_min_positions
            and s_pos >= params.rescue_min_positions):
        return (OUTCOME_MIXED, "medium",
                "%s has %d discriminating read(s) at %d position(s) and %s has "
                "%d at %d; each stands without the other"
                % (best["label"], b_reads, b_pos, second["label"], s_reads, s_pos))

    # One species is told apart from its relatives by reads they cannot claim.
    if (b_reads >= params.rescue_min_species_reads
            and b_pos >= params.rescue_min_positions
            and b_reads >= params.rescue_species_advantage * max(s_reads, 1)):
        conf = "high" if b_pos >= 2 * params.rescue_min_positions else "medium"
        return (OUTCOME_SPECIES, conf,
                "%d discriminating read(s) at %d independent position(s), against "
                "%d for %s" % (b_reads, b_pos, s_reads, second["label"] or "the runner-up"))

    # Reads belong to the genus and no further. Reporting a species here would
    # be a guess dressed as a result.
    if n_genus >= params.rescue_min_genus_reads:
        return (OUTCOME_GENUS, "medium",
                "%d read(s) fit the genus but no species is distinguished "
                "(%d discriminating read(s) for %s, %d for %s)"
                % (n_genus, b_reads, best["label"], s_reads, second["label"] or "none"))

    return (OUTCOME_INSUFFICIENT, "low",
            "%d mapped read(s): %d species-informative, %d genus-informative; "
            "not enough to place the signal" % (n_total, n_species, n_genus))


def run(review_rows, params, genus_lookup=None, out_path=None):
    """Score every reviewed row that carried a competitive mapping.

    `review_rows` are the rows the evidence review produced, each with the
    per-read detail attached under "_spec".
    """
    genus_of = genus_lookup or (lambda t: None)
    out = []
    for row in review_rows:
        spec = row.get("_spec") or {}
        per_read = spec.get("per_read") or {}
        target_label = spec.get("target_label", "__target__")

        if not per_read:
            out.append(dict(row.get("_report", {}),
                            primary_detection="not detected",
                            primary_detection_basis=(row.get("_report") or {}).get("detection_basis", ""),
                            rescue_triggered="no", taxonomic_outcome=OUTCOME_INSUFFICIENT,
                            confidence="low", taxonomic_call="",
                            decision_reason="no competitive mapping was possible"))
            continue

        # The target is carried as a placeholder label, so its genus has to be
        # looked up under its real name: without this, a tie between the target
        # and a sibling of its own genus reads as ambiguous rather than as
        # evidence for the genus, which is exactly the case this step exists for.
        species_name = ((row.get("_row") or {}).get("species_name")
                        or row.get("strain_name") or "")

        def genus_here(label, _sp=species_name):
            return genus_of(_sp) if label == target_label else genus_of(label)

        classified = classify_reads(per_read, target_label, genus_here,
                                    margin=params.evidence_score_margin)
        counts = {
            "mapped": len(classified),
            "species_informative": sum(1 for c in classified.values()
                                       if c["kind"] == SPECIES_INFORMATIVE),
            "genus_informative": sum(1 for c in classified.values()
                                     if c["kind"] == GENUS_INFORMATIVE),
            "ambiguous": sum(1 for c in classified.values() if c["kind"] == AMBIGUOUS),
        }

        labels = {c["winner"] for c in classified.values()}
        scored = []
        for label in labels:
            e = species_evidence(classified, label)
            e["label"] = label
            scored.append(e)
        scored.sort(key=lambda e: (-e["reads"], -e["positions"], str(e["label"])))
        best = scored[0] if scored else {"label": None, "reads": 0, "positions": 0,
                                         "median_margin": 0}
        second = scored[1] if len(scored) > 1 else {"label": None, "reads": 0,
                                                    "positions": 0, "median_margin": 0}

        control_ok = row.get("control_verdict") != "fail"
        outcome, confidence, why = decide(target_label, best, second, counts, params,
                                          control_ok=control_ok,
                                          target_label=target_label)

        def naming(label):
            return (row.get("strain_name") or row.get("virus_code") or "this taxon") \
                if label == target_label else (label or "")

        if outcome == OUTCOME_SPECIES:
            call = naming(best["label"])
        elif outcome == OUTCOME_GENUS:
            g = genus_here(best["label"]) or genus_here(target_label)
            call = (g + " sp.") if g else "genus level, genus unresolved"
        elif outcome == OUTCOME_MIXED:
            call = "%s + %s" % (naming(best["label"]), naming(second["label"]))
        else:
            call = ""

        detail = "; ".join("%s: %d discriminating read(s), %d position(s)"
                           % (naming(e["label"]), e["reads"], e["positions"])
                           for e in scored[:4])

        # Sections 1 and 2 of the report come from the primary path, under the
        # names the evidence review gave them.
        rep = dict(row.get("_report", {}))
        rep.setdefault("primary_detection", "not detected")
        rep["primary_detection_basis"] = rep.get("detection_basis", "")
        rep["reference_panel_size"] = rep.get("n_competing_references", "")
        rep["panel_species"] = rep.get("competing_references_tested", "")

        out.append(dict(
            rep,
            rescue_triggered="yes",
            reads_mapped=counts["mapped"],
            reads_species_informative=counts["species_informative"],
            reads_genus_informative=counts["genus_informative"],
            reads_ambiguous=counts["ambiguous"],
            best_species=naming(best["label"]),
            best_species_informative_reads=best["reads"],
            best_species_positions=best["positions"],
            second_species=naming(second["label"]),
            second_species_informative_reads=second["reads"],
            second_species_positions=second["positions"],
            score_margin_median=best["median_margin"],
            competing_species_detail=detail,
            taxonomic_outcome=outcome,
            confidence=confidence,
            taxonomic_call=call,
            decision_reason=why,
        ))

    if out_path:
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=REPORT_COLS, delimiter="\t",
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for r in out:
                w.writerow({k: r.get(k, "") for k in REPORT_COLS})
    by_outcome = {}
    for r in out:
        by_outcome[r["taxonomic_outcome"]] = by_outcome.get(r["taxonomic_outcome"], 0) + 1
    for k in (OUTCOME_SPECIES, OUTCOME_GENUS, OUTCOME_MIXED, OUTCOME_INSUFFICIENT):
        if by_outcome.get(k):
            log.info("  %-32s %d", k, by_outcome[k])
    return out
