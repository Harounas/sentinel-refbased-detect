# The workflow, end to end

Three paths, run in order. The first asks *is this virus present?*, the second
asks *is the evidence any good?*, the third asks *what is it evidence of?*
Nothing in path 2 or 3 overturns a detection; they explain the cases the
coverage rule cannot settle on its own.

## Path 1 — primary detection

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
  contiguous            one unbroken block >= 1000 bp at depth >= 2
  multi-region          >= 2 regions each holding a block >= 200 bp
  two-position          >= 2 unique reads at >= 2 separate loci
                        (overlapping reads merge into ONE locus)
  │
DIAMOND blastx of the block vs nr (optional, --diamond-db)
  │
negative-control comparison (reads per million vs the NTC/NEC)
  │
final_call   detected, confirmed / detected, unconfirmed /
             suspect: control background / not detected / control sample
```

**Why loci and not reads.** Ten reads stacked on one spot are one piece of
genome seen ten times — a repeat, a conserved stretch, a PCR duplicate that
survived marking. Two reads at unconnected loci are two independent
observations. The two-position rule counts places, so a pile-up cannot clear
it.

## Path 2 — evidence review

Runs for a priority pathogen that carried reads but failed the rule. Seven
criteria, scored independently:

```
sequence identity to the closest genome
mapping quality and multi-mapping rate
distribution across the genome (bins covered, busiest bin)
genomic coverage (breadth, longest block)
negative controls
background and index hopping
reference specificity ──── competitive mapping, EM-weighted attribution
          │
  consistent (>= 4 of the applicable criteria pass)
  equivocal  (exactly one failure)
  not supported
```

The verdict is advisory. It never flips `detection`; it tells you whether the
reads deserve a second look.

## Path 3 — genus-level rescue

```
competitive mapping against three layers
  │
  ├─ siblings    genus panel, 50 genera, 858 genomes (bundled)
  ├─ strains     per-species clustering, 83 species (built on Terra)
  └─ background  all RefSeq viral + host + reagent contaminants
  │
every read classified by the margin between its best and second-best hit
  │
  ├─ species-informative   no relative explains it as well
  ├─ genus-informative     the species it cannot choose between share a genus
  └─ ambiguous             spans genera, or the background explains it too
  │
discriminating reads per species, and how many distinct positions they sit at
  │
OUTCOME
  species-level identification   one species holds reads the others cannot claim
  genus-level identification     reads fit the genus, no species distinguished
                                 ──► reported as "Genus sp."
  possible mixed infection       two species each hold their own, spread out
  insufficient evidence          too little, or matches the controls
```

A species is never claimed on read counts or a point of score. It is claimed
only on reads a relative could not have produced, at positions that are not the
same position twice.

## The databases

Built once, on Terra, then reused. Both are WDLs in [`wdl/`](../wdl).

| Workflow | What it builds |
|---|---|
| `build_pathogen_panel.wdl` | 133 shards → genus sibling panel + clustered strain panel. Redundancy removed by Mash distance (k=16) with a measured coverage radius per species, not one fixed threshold; singletons and two-member clusters retained so no species drops out. |
| `build_background_db.wdl` | 14,844 viral species from the RefSeq viral release, plus the host genome and reagent contaminants (~3.6 Gb). |

The panel is not drawn only from the priority table: every genus in the table
is expanded to its siblings, so a read from a relative has somewhere to go.

## Where each column comes from

| Column | Written by |
|---|---|
| `detection`, `detection_basis`, `detection_rule` | path 1, `summarize.py` |
| `n_distinct_start_positions`, `n_separate_loci` | path 1, `align.py` |
| `n_references_tied`, `tied_references` | path 1, best-reference selection |
| `final_call` | path 1, `calls.py`, after the control comparison |
| `evidence_verdict` and the seven criteria | path 2, `evidence.py` |
| `taxonomic_outcome`, `taxonomic_call`, `confidence` | path 3, `rescue.py` |
