"""Task final_calls: negative-control comparison and the combined call."""

import csv
import re

from .utils import PipelineError, log

NONVIRAL = {"Bacteria", "Archaea", "Human", "Vertebrate", "Fungi", "Protozoa",
            "Invertebrate", "Plant", "Other eukaryote", "Other"}
NEW_COLS = ["is_control", "control_signal", "control_reads", "control_sample",
            "control_fold", "control_status", "final_call", "final_call_reason"]
BACKGROUND_COLS = ["taxid", "virus_code", "strain_name", "control_sample",
                   "control_reads", "control_signal", "control_detection",
                   "control_longest_block_bp", "n_detected_at_background"]


def final_calls(detection_summary_tsv, final_out, background_out, params, library_sizes=None):
    """library_sizes: {sample_id: reads}; empty/None -> compare raw read counts."""
    try:
        control_re = re.compile(params.control_sample_regex, re.I)
    except re.error as exc:
        raise PipelineError("control_sample_regex is not a valid regex: %s" % exc)

    with open(detection_summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        header = list(reader.fieldnames or [])

    lib = {s: int(n) for s, n in (library_sizes or {}).items() if n and int(n) > 0}
    normalized = bool(lib)

    def reads(row):
        try:
            return int(row.get("total_unique_reads_mapped") or 0)
        except ValueError:
            return 0

    def signal(row):
        n = reads(row)
        total = lib.get(row["sample_id"], 0)
        return (1e6 * n / total) if (normalized and total) else float(n)

    samples = sorted({r["sample_id"] for r in rows})
    controls = {s for s in samples if control_re.search(s)}
    if rows and not controls:
        log.warning("no sample matched %r; every row is reported as 'no control on flowcell'",
                    control_re.pattern)

    shallow = set()
    if normalized:
        shallow = {s for s in controls if lib.get(s, 0) < params.min_control_library_reads}
        for s in sorted(shallow):
            log.warning("control %s has %d reads, under min_control_library_reads=%d; "
                        "reported but not used as background",
                        s, lib.get(s, 0), params.min_control_library_reads)
    background_controls = controls - shallow
    log.info("%d samples, %d controls: %s", len(samples), len(controls),
             ", ".join(sorted(controls)) or "none")
    log.info("%d usable as background: %s", len(background_controls),
             ", ".join(sorted(background_controls)) or "none")

    bg = {}
    for r in rows:
        if r["sample_id"] not in background_controls:
            continue
        if reads(r) < params.min_control_support_reads:
            continue
        cur = bg.get(r["taxid"])
        s = signal(r)
        if cur is None or s > cur["signal"]:
            bg[r["taxid"]] = {"signal": s, "reads": reads(r), "sample": r["sample_id"],
                              "detection": r.get("detection", ""),
                              "block": r.get("longest_covered_block_bp", "")}

    def pident(row):
        try:
            return float(row.get("blast_pct_identity") or 0)
        except ValueError:
            return 0.0

    n_by_call = {}
    for r in rows:
        is_control = r["sample_id"] in controls
        b = bg.get(r["taxid"])
        s = signal(r)
        fold = ""
        if not is_control and b and b["signal"] > 0:
            fold = s / b["signal"]

        if is_control:
            control_status = "control sample"
        elif not background_controls:
            control_status = ("no control on flowcell" if not controls
                              else "no control deep enough to use")
        elif b is None or b["signal"] <= 0:
            control_status = "absent from controls"
        elif fold != "" and fold < params.min_control_fold:
            control_status = "at or near control background"
        else:
            control_status = "above control background"

        detected = r.get("detection") == "detected"
        cat = (r.get("blast_best_subject_category") or "").strip()
        taxon = (r.get("blast_best_subject_taxon") or "").strip()
        agree = (r.get("blast_agreement") or "").strip()

        if is_control:
            call, why = "control sample", "negative or process control"
        elif not detected:
            call, why = "not detected", r.get("detection_basis", "")
        elif control_status == "at or near control background":
            call = "suspect: control background"
            why = ("%.3g vs %.3g in %s (%.1fx, under %gx)"
                   % (s, b["signal"], b["sample"], fold, params.min_control_fold))
        elif params.min_unique_reads and reads(r) < params.min_unique_reads:
            call = "suspect: too few reads"
            why = "%d unique reads, under min_unique_reads=%d" % (reads(r), params.min_unique_reads)
        elif cat in NONVIRAL and pident(r) >= params.min_nonviral_pident:
            call = "suspect: non-viral block"
            why = ("best nr hit is %s (%s) at %.1f%% identity: the covered block "
                   "looks like host sequence, not virus" % (taxon or "?", cat, pident(r)))
        elif agree == "different taxon":
            call = "suspect: conflicting taxon"
            why = "nr names %s where the reference claims %s" % (
                r.get("diamond_taxon") or "?",
                r.get("species_name") or r.get("strain_name") or "?")
        elif agree in ("species match", "genus match"):
            call = "confirmed"
            why = "coverage rule plus nr %s (%s)" % (agree, r.get("diamond_taxon") or "?")
        else:
            call = "detected, unconfirmed"
            why = "coverage rule only; nr gave %s" % (agree or "no result")

        if control_status.startswith("no control") and call not in ("not detected", "control sample"):
            why += " [%s]" % control_status

        r.update({
            "is_control": "yes" if is_control else "no",
            "control_signal": ("%.4f" % s) if normalized else str(reads(r)),
            "control_reads": str(b["reads"]) if b else "",
            "control_sample": b["sample"] if b else "",
            "control_fold": ("%.2f" % fold) if fold != "" else "",
            "control_status": control_status,
            "final_call": call,
            "final_call_reason": why,
        })
        n_by_call[call] = n_by_call.get(call, 0) + 1

    cols = [c for c in header if c not in NEW_COLS] + NEW_COLS
    with open(final_out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    with open(background_out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(BACKGROUND_COLS)
        name = {}
        for r in rows:
            name.setdefault(r["taxid"], (r.get("virus_code", ""), r.get("strain_name", "")))
        for taxid, b in sorted(bg.items(), key=lambda kv: -kv[1]["signal"]):
            near = sum(1 for r in rows if r["taxid"] == taxid
                       and r["final_call"] == "suspect: control background")
            code, strain = name.get(taxid, ("", ""))
            w.writerow([taxid, code, strain, b["sample"], b["reads"],
                        "%.4f" % b["signal"], b["detection"], b["block"], near])

    log.info("normalization: %s", "reads per million" if normalized
             else "raw read counts (no library sizes)")
    for call in sorted(n_by_call, key=lambda c: -n_by_call[c]):
        log.info("  %-30s %d", call, n_by_call[call])
    return n_by_call
