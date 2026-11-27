# sentinel-refbased-detect

Run the Terra workflow **`warn-id/sentinel_refbased_detect`** (method v19, config
`sentinel_refbased_detect_0109_V15`) on your own machine. It does
reference-based detection of priority pathogens from FASTQ files.

The package ports every WDL task to Python (`wdl/sentinel_refbased_detect.wdl` is
the original, kept for reference). It keeps the same detection rule, the same
defaults as the V15 configuration, and the same output columns and file names.
Results can be compared directly with the Terra outputs.

```
FASTQ ─► Kraken2 ─► candidate selection ─► NCBI references ─► bwa + samtools + ivar
            │        (priority table +         (cached)          per reference
            │         alias/legacy taxids)                               │
            └─► discovery of recurrent non-priority viruses              ▼
                                         best reference per virus + detection rule
                                                         │
                       (optional) DIAMOND blastx of the longest covered block vs nr
                                                         │
                              negative-control comparison ─► final call
                                                         │
                     priority pathogen that failed the rule on a few reads?
                                                         ▼
                        evidence review: identity, mapping quality, genome
                        distribution, coverage, controls, background,
                        reference specificity ─► advisory verdict
```

## Install

The Python code uses only the standard library. The external tools come from bioconda:

```bash
git clone git@github.com:Harounas/sentinel-refbased-detect.git
cd sentinel-refbased-detect
mamba env create -f environment.yml      # kraken2, bwa, samtools, ivar, diamond + this package
conda activate sentinel-detect
sentinel-detect check                    # confirms the tools are on PATH
```

To use an existing environment instead: `pip install git+https://github.com/Harounas/sentinel-refbased-detect.git`
and install `kraken2 bwa samtools ivar` (and `diamond` if you want it) yourself.

On Apple Silicon, if a bioconda package has no arm64 build, create the env with
`CONDA_SUBDIR=osx-64 mamba env create -f environment.yml`.

## Inputs

| Input | Flag | Notes |
|---|---|---|
| Kraken2 database | `-k/--kraken-db` | A directory containing `hash.k2d`, `opts.k2d` and `taxo.k2d` (the same as `gs://krakendb/`), or a `.tar.gz`/`.tar.zst` archive of one |
| Reference table | `-r/--ref-tsv` | **Optional.** `Refgenomes_V15.tsv` (the same as `gs://lasvdata/Refgenomes_V15.tsv`) is bundled in [`src/sentinel_refbased_detect/data/`](src/sentinel_refbased_detect/data/Refgenomes_V15.tsv) and used by default. Pass your own table to override it: taxid, code, name, colon-joined accessions, no header; a local path, `gs://…` or `https://…` |
| FASTQ directory | `-i/--fastq-dir` | `*_R1_001.fastq.gz`/`*_R2_001.fastq.gz`, `*_1.fq.gz`, single-end files, multiple lanes |
| NCBI e-mail | `--ncbi-email` or `$NCBI_EMAIL` | Required by NCBI E-utilities. `--ncbi-api-key`/`$NCBI_API_KEY` raises the rate limit |

**Sample IDs** follow the Terra convention: the file name up to the first `.`, with
Illumina suffixes removed. For example `1-02-26-0109_S5_L001_R1_001.fastq.gz` becomes
`1-02-26-0109`. Negative controls are recognised by ID (`--control-sample-regex`,
default `^(nec|ntc|neg|negative|blank|water|h2o)([-_.]|$)`), so name control
files `NEC…`, `NTC…` and so on. If your file names don't fit this pattern, use
`--samplesheet samples.tsv` with columns `sample_id  fastq_1  fastq_2`.

## Run

```bash
export NCBI_EMAIL=you@example.org

sentinel-detect run \
    --kraken-db /data/krakendb \
    --fastq-dir /data/run_261029/fastq \
    --outdir    results_261029 \
    --threads 8 --jobs 4
```

This uses the bundled `Refgenomes_V15.tsv`. To change the pathogen panel:

```bash
sentinel-detect refs                          # where the bundled table is, and how many rows
sentinel-detect refs --export my_panel.tsv    # copy it, edit it, then run with:
sentinel-detect run --ref-tsv my_panel.tsv ...
```

Optional DIAMOND confirmation of detected blocks (the same nr stage as on Terra):

```bash
    --diamond-db /data/nr.dmnd \
    --accession2taxid /data/prot.accession2taxid.gz   # downloaded from NCBI if omitted
```

### Host depletion (`--host-ref`)

If your FASTQs are raw (not yet host-depleted), add a host genome. Reads are aligned
to it with minimap2 (`-x sr`). A read pair is kept only when **neither** mate maps
to the host, and everything downstream runs on the kept reads.

```bash
    --host-ref GRCh38.fa.gz          # or T2T CHM13; FASTA(.gz) or a prebuilt .mmi
```

* With a FASTA, the minimap2 index is built once and cached in `OUTDIR/cache/host/`
  (or `--cache-dir`). For a human genome this takes about 10 minutes and about 16 GB of RAM.
  To build it yourself ahead of time: `minimap2 -x sr -d GRCh38.sr.mmi GRCh38.fa.gz`.
* The depleted reads are kept in `OUTDIR/host_depleted/`, and `host_depletion.tsv`
  lists reads in, kept and removed (% host) for each sample.
* Library sizes for the control comparison are computed **after** depletion, as on
  Terra (`read_counts_depleted`).
* Depletion is deliberately conservative. Only reads that are near-identical to the host
  are removed, so divergent viral reads are never lost to it.
* Leave `--host-ref` unset if your reads are already cleaned, for example FASTQs made
  from Terra's `*.cleaned.bam`.

Other useful options:

* `--kraken-reports DIR`: reuse existing Kraken2 reports (`<sample_id>.*report*.txt`, for
  example ones downloaded from Terra) and skip Kraken2. `--kraken-db` is then not needed.
* `--taxdump taxdump.tar.gz`: pin the NCBI taxonomy. By default the current one is downloaded once and cached.
* Resuming is on by default. Rerun the same command after an interruption and finished
  Kraken2 reports and alignments are reused. Use `--no-resume` to recompute everything.
* `--cache-dir DIR`: share the reference and taxonomy cache between runs.
* Every V15 parameter is a flag, for example `--min-depth`, `--min-region-bp`,
  `--min-block-bp`, `--min-qualifying-regions`, `--min-mapq`, `--min-control-fold`.
  See `sentinel-detect run --help`.

### Resources

Kraken2 loads the whole database for each sample, so it needs as much RAM as the
database is large (`--kraken-memory-mapping` trades speed for memory). Samples are
classified one at a time. Alignments run `--jobs` at a time with `--threads` each.

## Detection rule (unchanged from the WDL)

A virus is **detected** when the best of its candidate references meets either condition:

* one region (segment or genome) has an unbroken run of at least `min_region_bp` (1000) bases at depth `min_depth` (2) or more, using unique reads (MAPQ ≥ 20, duplicates removed); or
* at least `min_qualifying_regions` (2) separate regions each have a block of at least `min_block_bp` (200) such bases.

The final call (`final_call`) then applies these checks in order:

| Order | Check | `final_call` |
|---|---|---|
| 1 | Sample is a negative control | `control sample` |
| 2 | Detection rule not met | `not detected` |
| 3 | Signal is less than `min_control_fold` (10) times the loudest control | `suspect: control background` |
| 4 | Fewer than `min_unique_reads` unique reads (off by default) | `suspect: too few reads` |
| 5 | Best nr hit is non-viral at ≥ `min_nonviral_pident` (80%) identity | `suspect: non-viral block` |
| 6 | nr names a different virus | `suspect: conflicting taxon` |
| 7 | nr agrees at species or genus level | `confirmed` |
| 8 | Otherwise (for example, no DIAMOND run) | `detected, unconfirmed` |

## Evidence review of priority pathogens that fail the coverage rule

A priority pathogen that Kraken2 names but that falls short of the coverage rule is
reported `not detected`, and that call stands. It does not say whether the few reads
that exist look like real pathogen sequence or like an artefact, which is the judgement
a reviewer would otherwise make by hand.

The review runs automatically after the calls (turn it off with `--no-evidence-review`),
on alignments that already exist, so it costs almost nothing. Rows are reviewed when the
virus is **priority**, `detection` is **not detected**, the sample is **not a control**,
and Kraken2 assigned at least `--evidence-min-kraken-reads` reads. Each row is scored
against seven criteria:

| Criterion | What is measured | Passes when (default) |
|---|---|---|
| Sequence identity | Percent identity of the mapped reads to the reference, from the NM tag | mean ≥ 95% |
| Mapping quality | Mean MAPQ, and reads carrying an alternative alignment (XA/SA) | MAPQ ≥ 30 and ≤ 20% multi-mapping |
| Genome distribution | Spread over 500 bp windows, vs piled in one spot | ≥ 3 windows and ≤ 80% of reads in one |
| Genome coverage | Breadth of the reference covered | ≥ 0.1% |
| Negative controls | The same virus in this run's controls | sample ≥ 10× the loudest control |
| Background | A much louder sample on the run that could have bled in | no run-mate ≥ 100× louder |
| Reference specificity | Whether the reads can be attributed to this taxon when sibling species and strain genomes compete for them | ≥ 80% of reads attributed |

Specificity treats references of the **same species** as on-target: three adenovirus C
references trading reads among themselves says nothing about whether the sequence is
adenovirus. What matters is whether it is distinguishable from other taxa.

### What the reads are weighed against: the decoy panel

Comparing only against the references Kraken2 happened to flag in that same sample has a
blind spot: if the reads really come from a relative the classifier never called, nothing
competes and specificity passes by default. So the check also maps against a **decoy
panel**, bundled with the package and generated rather than hand-curated:

| Role | What it is | Why |
|---|---|---|
| `sibling` | One RefSeq genome per **other species of the same genus** | These are the genomes a read from a neighbouring virus would really prefer |
| `strain` | Genomes of the **priority species itself**, clustered to representatives | A true read must have somewhere right to go; without strains it looks ambiguous merely because the panel holds one arbitrary lineage |

Genus, not family: genus holds the genuinely confusable genomes, while family mostly adds
distant ones that would never win a read, at a large cost in index size. Where a taxon has
no genus, the builder falls back to family.

Strain genomes are deduplicated by **MinHash clustering** (Mash distances, no external
tool), keeping the longest genome of each group above the identity threshold. Hundreds of
near-identical deposits would otherwise split one taxon's evidence across copies of the
same sequence.

```bash
sentinel-detect panel show                    # what the bundled panel covers
sentinel-detect panel build                   # regenerate the sibling panel
sentinel-detect panel build --kinds sibling strain --identity 0.95
sentinel-detect run --decoy-panel ''          # or switch it off
```

#### Divergent species: Lassa, the coronaviruses

Some priority pathogens are far more diverse inside one species than others, and a
clustering that ignores that would quietly do the wrong thing. Three places where it is
handled:

**Segments are clustered separately.** Lassa, CCHF and influenza deposit one record per
segment. An L segment and an S segment share almost nothing, so clustering them together
calls every genome unique and keeps nothing useful. (NCBI titles these "segment L,
complete **sequence**", not "complete genome", so the query has to ask for both — with the
narrower phrase Lassa returns zero records and is skipped entirely.)

**Representatives span the diversity.** Lassa's clades sit around 65% identical to one
another, so they survive clustering as separate representatives, which is what should
happen. The risk is the cap afterwards: keeping the first N would drop a whole clade. The
builder takes the longest genome first, then repeatedly the representative furthest from
everything already chosen, so a cap trims the crowded middle of a clade instead.

**One identity threshold does not fit all viruses**, so the panel does not use one for
this. `--identity` (99%) only removes redundancy: two genomes that similar are the same
thing. How many representatives a species keeps is decided instead by
`--coverage-identity` (90%), a promise about the result: *every known genome of this
species is at least this similar to some representative*. Selection stops as soon as that
holds, so the number of representatives follows the species rather than a fixed rule:

| Species | Genomes | After 99% dedup | Kept | Coverage reached |
|---|---|---|---|---|
| SARS-CoV-2 | 120 | 2 clusters | **1** | every genome within 99% |
| Lassa, segment L | 40 | 31 clusters | **12** (cap) | within 83%, flagged `CAPPED` |
| Lassa, segment S | 53 | 42 clusters | **12** (cap) | within 82%, flagged `CAPPED` |

A single 95% rule would have been wrong in both directions here: twelve near-identical
SARS-CoV-2 genomes carrying no extra information, and Lassa silently under-covered with no
sign that anything was missing. `panel_diversity.tsv` sorts species by the coverage
actually reached, worst first, and marks the ones that hit the cap — those are the species
to give a higher `max_representatives`.

**Every pathogen keeps at least one decoy.** Each step after the genomes are fetched makes
the panel smaller — clustering, the cap, dedup across shards — and a pathogen can lose its
decoys at any of them without anything failing. Three places where that was happening, and
what each does now:

| Where a virus could vanish | What happens instead |
|---|---|
| A genus holding **only** priority species, such as Orthomarburgvirus (Marburg alone) or Lentivirus (HIV-1 and HIV-2): excluding the guarded taxa left nothing | The search widens to the family, so Marburg is offered the other filoviruses |
| A species or a **segment with a single deposit** was skipped as "nothing to cluster" | A lone genome is kept: it is not a clustering problem, it is the only genome that exists |
| Dedup across shards kept the **first copy** of a genome reached from two genera | Duplicate rows are merged and their guard lists combined, so neither taxon loses the decoy |

```bash
sentinel-detect panel verify     # one row per priority pathogen: decoys found, or NO DECOYS
```

`panel_coverage.tsv` lists every priority row with how many sibling species and strain
representatives survived for it. Run it after any rebuild: it is the check that the
shrinking stopped short of removing something. The bundled panel stands at **858 genomes
covering 354 of the 355 priority rows**.

The one exception is real and not a bug: **African swine fever virus** has no siblings
anywhere, because *Asfarviridae* holds a single genus holding a single species. Where
nothing can compete, the specificity criterion reports `not assessable` rather than
passing — a pathogen with no relatives must not score a free pass on a test that had
nothing to test it against. For ASFV the strain layer is the only decoy available.

**Lone outliers cannot displace a real lineage.** Furthest-first selection reaches for
whatever is most unlike everything else, which is exactly what a misannotated, chimeric or
heavily-N'd deposit looks like. A cluster therefore needs `--min-cluster-members`
(2 by default) before it can be picked for diversity; singletons are still taken when
nothing corroborated remains.

**The measurement itself has to reach that far.** MinHash uses k=16 rather than Mash's
usual 21: a 21-mer survives in two sequences 76% identical only 0.4% of the time, so the
sketch saturates and every divergent pair reads as 0% identity — exactly the species the
number was meant to describe. At k=16 Lassa's clades measure 65%, which is about right.

Coronaviruses need none of this: SARS-CoV-2, MERS and the seasonal ones are separate
species with their own taxids and their own rows, so each gets its own references and its
own siblings.

The same divergence affects the **identity criterion**. A Lassa clade II read mapped
against a Josiah reference is ~80% identical and would fail a 95% floor, though it is a
perfect match to its own clade. So identity is measured against the **closest genome of
the taxon** in the panel (`pct_identity_closest_genome`), falling back to the aligned
reference when the panel has nothing. The detection side already handles this a different
way: `Refgenomes_V15.tsv` carries five Lassa clade references under one taxid, and the
pipeline aligns to all of them and keeps the best.

#### Building the panel on Terra

For the full database, `wdl/build_pathogen_panel.wdl` does the same work as a WDL
workflow: one shard per genus for siblings and one per species for strains, which on the
V15 table is 133 shards running in parallel instead of hours in series on one machine.

```
list_jobs ──► scatter over 50 genera + 83 species ──► merge_panel
  (one          sibling: 1 RefSeq genome per             decoy_panel.tsv
  taxdump       other species of the genus               panel_diversity.tsv
  for all)      strain:  cluster, per segment,
                         keep diverse representatives
```

Register it in Terra like the detection workflow, give it `ref_genomes_tsv` and
`ncbi_email`, and point `--decoy-panel` at the `panel_tsv` it writes. `panel_diversity.tsv`
lists each species' lineage spread, lowest first, so the species that need more
representatives say so themselves. A shard that fails leaves a gap in the decoys rather
than taking the whole panel down.

Rebuild it deliberately, read the diff, bump the version. Set `NCBI_EMAIL` to your own
address first; the bundled panel was generated with a generic contact.

### How the criterion is scored

Not "did the reads stay put". With sibling species and strain genomes in the pool, a
genuine read often fits several genomes equally well, and counting that as a loss would
fail real findings. The question is **can these reads be attributed to this taxon**:

1. Every read is realigned against the whole pool with `bwa mem -a`, keeping all alignments.
2. A read that fits **one taxon alone** is counted there (`n_reads_unique_to_target`,
   `n_reads_unique_elsewhere`).
3. A read that fits **several taxa equally** is ambiguous (`n_reads_ambiguous`) and is
   shared out between them in proportion to the reads that fit only one — so a strain with
   real unique support claims the share its own evidence earns, which is the weighting by
   reads mapped that this is built around.
4. `pct_reads_attributed_to_target` is the resulting share, and the criterion passes at
   `--evidence-min-retained` (80% by default).

Strain genomes of the same taxon count as on-target, as do several references of one
species: three adenovirus C references trading reads among themselves says nothing about
whether the sequence is adenovirus. `--evidence-score-margin` widens what counts as an
equal fit, if you want near-ties treated as ambiguous too.

The review records what the reads were weighed against, not just the headline percentage:

| Column | Example |
|---|---|
| `pct_reads_attributed_to_target` | `92.5` |
| `n_reads_unique_to_target` / `n_reads_unique_elsewhere` / `n_reads_ambiguous` | `20` / `0` / `10` |
| `reads_by_taxon` | `this taxon:28.5;Mopeia mammarenavirus:1.5` |
| `n_competing_references` | `2` |
| `competing_references_tested` | `Human mastadenovirus C; Human adenovirus 5` |
| `reads_by_reference` | `10533__HAdV-1__0 [target]:3;129951__HAdV-C__0 [other taxon]:7` |
| `competing_best_reference` | `129951__HAdV-C__0 = Human mastadenovirus C (7 reads)` |

`reads_by_reference` is keyed by reference tag (`<taxid>__<code>__<index>`), because
several references of one species share a name and the question is which reference each
read preferred. `n_competing_references` of `0` means every other reference offered was
of the same taxon, so nothing genuinely competed.

A criterion that cannot be judged, such as controls on a run that has none, is recorded
as `not assessable` and left out of the tally. The overall `evidence_verdict` is:

| Verdict | Meaning |
|---|---|
| `consistent with pathogen sequence` | Every assessable criterion passes. Too little sequence to call it present, but worth confirming with a targeted assay such as PCR |
| `equivocal` | One criterion fails |
| `not supported` | Two or more fail |
| `no mapped reads` | Kraken2 assigned reads but none aligned to the reference |
| `not reviewed` | The row did not qualify, for example it was already detected |

**The verdict is advisory and never changes `detection` or `final_call`.** It is evidence
for a human decision, not a positive result.

Full per-criterion numbers go to `evidence_review.tsv`; `evidence_verdict` and
`evidence_reason` are added to the final summary. Every threshold is a flag, for example
`--evidence-min-identity`, `--evidence-min-bins`, `--evidence-min-retained`. See
`sentinel-detect run --help`.

## Outputs (`--outdir`)

| File | Terra output |
|---|---|
| `detection_summary_final.tsv` | `detection_summary_tsv` (the main result) |
| `detection_summary.tsv` / `detection_summary_annotated.tsv` | `detection_summary_precall_tsv` |
| `detection_by_region.tsv` | `detection_by_region_tsv` |
| `control_background.tsv` | `control_background_tsv` |
| `evidence_review.tsv` | *(new)* per-criterion review of failing priority rows |
| `kraken2/<sample>.kraken2.report.txt` | `kraken2_reports` |
| `samples/<sample>/candidates.tsv` | `candidate_tables` |
| `samples/<sample>/consensus/*.fasta` | `consensus_tars` (unpacked) |
| `samples/<sample>/align/<ref>/unique.bam` | per-reference alignment |
| `discovery/discovered_taxa.tsv` | `discovered_taxa_tsv` |
| `taxonomy/all/combined_refs.tsv` | `combined_refs_tsv` |
| `diamond/` | `block_blast_hits`, `block_taxonomy`, `blocks_fasta` |
| `host_depleted/`, `host_depletion.tsv` | cleaned reads and per-sample host fraction (with `--host-ref`) |
| `library_sizes.tsv`, `run_params.json`, `pipeline.log` | run provenance |

## Differences from the Terra run

* **Input is FASTQ instead of cleaned unaligned BAMs.** Give it host-depleted reads (the
  `cleaned` reads), or raw reads with `--host-ref`. Terra's depletion workflow used other
  tools, so the exact reads removed can differ slightly. To turn a Terra BAM into FASTQ:
  `samtools fastq -1 X_R1.fastq.gz -2 X_R2.fastq.gz -0 /dev/null -s /dev/null -n X.cleaned.bam`.
* **Library sizes** for the reads-per-million control comparison (`sample_read_counts`
  on Terra) are taken from each Kraken2 report (classified + unclassified fragments,
  ×2 for paired reads). You can supply your own with `--read-counts` or switch
  normalization off with `--no-normalize`.
* NCBI references are cached and downloaded once per run, not once per sample.
* bwa output is streamed straight into `samtools sort`, so no SAM file is written. The
  coverage-block counting that the WDL did in `awk` is done in Python with the same logic.
* If no virus reaches alignment in any sample, the run writes empty final tables
  instead of failing.

## Tests

```bash
pip install -e ".[test]" && pytest
```
