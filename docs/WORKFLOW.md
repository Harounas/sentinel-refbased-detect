# The workflow, end to end

This repository is the first of two. It answers one question — *is there
enough sequence to call this virus present?* — the way the Terra method does.

The rows it cannot settle, a priority pathogen that carried reads but failed
the coverage rule, are picked up by
[sentinel-genus-rescue](https://github.com/Harounas/sentinel-genus-rescue),
which runs over a finished output directory and writes its own files. Nothing
there changes a call made here.

## Primary detection

```
FASTQ
  │
  ├─ host depletion (optional, --host-ref)
  │
Kraken2
  │
  ├─ candidate selection ──── priority table + legacy/alias taxids
  │                           (Refgenomes_V15.tsv, bundled)
  └─ discovery of recurrent non-priority viruses
  │
NCBI references (cached locally)
  │
bwa mem ─► fixmate ─► markdup ─► MAPQ >= 20, primary, non-duplicate
  │
  ├─ samtools depth ──► covered bp, longest block, n blocks
  ├─ samtools view  ──► unique reads, mean MAPQ, distinct starts, separate loci
  └─ ivar consensus ──► consensus + longest covered block as FASTA
  │
best reference per virus
  │   covered bp ▸ separate loci ▸ unique reads ▸ mean MAPQ ▸ name
  │   runner-up within 5% on coverage and equal on loci
  │     ──► reported in n_references_tied / tied_references, not silently broken
  │
DETECTION RULE          fires           column: detection_rule
  contiguous            one unbroken block >= 1000 bp at depth >= 2   ] the Terra
  multi-region          >= 2 regions each holding a block >= 200 bp   ] rules
  two-position          >= 2 unique reads at >= 2 separate loci
                        (overlapping reads merge into ONE locus)
                        OFF by default; --min-specific-reads 2 adds it
  │
DIAMOND blastx of the block vs nr (optional, --diamond-db)
  │
negative-control comparison (reads per million vs the NTC/NEC)
  │
final_call   detected, confirmed / detected, unconfirmed /
             suspect: control background / not detected / control sample
```

**The two Terra rules are the default.** The package reproduces the V15 run as
it stands; the two-position rule is an addition you ask for, not a change to
what the first workflow calls.

**Why loci and not reads.** Ten reads stacked on one spot are one piece of
genome seen ten times — a repeat, a conserved stretch, a PCR duplicate that
survived marking. Two reads at unconnected loci are two independent
observations. The two-position rule counts places, so a pile-up cannot clear
it.

## What happens next, elsewhere

In [sentinel-genus-rescue](https://github.com/Harounas/sentinel-genus-rescue),
run as `sentinel-rescue review -d OUTDIR`:

1. **Evidence review** — seven criteria scored off the alignments this run
   already wrote: identity, mapping quality, genome distribution, coverage,
   negative controls, background and index hopping, reference specificity.
   Verdict: consistent / equivocal / not supported.
2. **Genus-level rescue** — competitive mapping against sibling species,
   clustered strains and a background database, with every read classified as
   species-informative, genus-informative or ambiguous before anything is
   named. Outcome: species level / genus level (`Genus sp.`) / possible mixed
   infection / insufficient evidence.

The panels those steps map against, and the two Terra WDLs that build them,
live in that repository too.

## Where each column comes from

| Column | Written by |
|---|---|
| `detection`, `detection_basis`, `detection_rule` | `summarize.py` |
| `n_distinct_start_positions`, `n_separate_loci` | `align.py` |
| `n_references_tied`, `tied_references` | best-reference selection |
| `final_call` | `calls.py`, after the control comparison |
| `evidence_verdict` and the seven criteria | sentinel-genus-rescue, `evidence.py` |
| `taxonomic_outcome`, `taxonomic_call`, `confidence` | sentinel-genus-rescue, `rescue.py` |
