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
| Reference specificity | Whether the reads stay on this virus when the run's other references compete for them, realigned with bwa | ≥ 80% retained |

Specificity treats references of the **same species** as on-target: three adenovirus C
references trading reads among themselves says nothing about whether the sequence is
adenovirus. What matters is whether it is distinguishable from other taxa.

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
