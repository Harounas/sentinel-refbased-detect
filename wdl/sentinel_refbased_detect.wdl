version 1.0

# Reference-based detection of priority pathogens in sentinel NGS data.
#
# For each cleaned, unaligned BAM on a flowcell:
#   1. Kraken2 taxonomic classification.
#   2. Keep taxids Kraken2 reports that appear in the priority-pathogen table
#      (gs://lasvdata/Refgenomes_V1.tsv), plus viruses outside it that recur
#      across the flowcell or carry real depth in one sample.
#   3. Fetch every candidate reference for those taxids from NCBI.
#   4. Align the reads to each candidate reference, call a consensus with ivar.
#   5. Pick the best reference per virus and apply the detection rule.
#
# Detection rule (min_depth = 2 by default):
#   detected  if any region holds one unbroken run of >= min_region_bp
#             bases at >= min_depth
#   detected  if, when no run reaches min_region_bp, at least
#             min_qualifying_regions separate regions each hold a covered
#             block of >= min_block_bp bases at >= min_depth
#   not detected otherwise
#
# "Region" means one sequence of the reference: a segment for Lassa, CCHF or
# influenza, and the single genome for an unsegmented virus such as MPXV.

workflow sentinel_refbased_detect {
    input {
        Array[File]   cleaned_reads_unaligned_bams
        File?         kraken2_db_tgz
        File?         kraken2_hash_k2d
        File?         kraken2_opts_k2d
        File?         kraken2_taxo_k2d
        Boolean       kraken2_memory_mapping = false
        File          ref_genomes_tsv
        File?         taxdump_tar_gz
        File?         name_regex_tsv
        String        ncbi_email
        String?       ncbi_api_key

        Int           min_kraken_reads   = 1
        Int           min_kraken_reads_discovered = 100
        Float         min_kraken_pct     = 0.0
        Int           min_depth          = 2
        Int           min_region_bp      = 1000
        Int           min_block_bp       = 200
        Int           min_qualifying_regions = 2
        Int           min_mapq           = 20
        Int           max_refs_per_virus = 20
        Boolean       expand_descendants = true
        Int           max_rows_above_species = 8
        String        min_match_rank     = "G"

        # DIAMOND identification of the longest covered block. Supply
        # diamond_db to switch it on; leave it unset and the workflow behaves
        # exactly as before.
        File?         diamond_db
        Array[File]   accession2taxid    = []
        File?         diamond_taxdump_tar_gz
        Float         diamond_top_pct    = 5.0

        Int           discover_min_samples          = 10
        Int           discover_min_reads_one_sample = 1000
        Int           discover_present_min_reads    = 1
        Int           discover_max_novel_taxa       = 50

        # Negative-control comparison and the combined call.
        String        control_sample_regex = "^(nec|ntc|neg|negative|blank|water|h2o)([-_.]|$)"
        Float         min_control_fold     = 10.0
        Int           min_control_library_reads = 50000
        Int           min_control_support_reads = 10
        Int           min_unique_reads     = 0
        Float         min_nonviral_pident  = 80.0
        Array[Int]    sample_read_counts   = []

        Int           consensus_min_q     = 20
        Float         consensus_min_freq  = 0.75
        Int           consensus_min_depth = 2
        String        consensus_ambig_char = "N"
        Int           consensus_retry_depth = 10000
    }

    parameter_meta {
        cleaned_reads_unaligned_bams: "Cleaned, host-depleted unaligned BAMs. On the flowcell table this is cleaned_reads_unaligned_bams."
        kraken2_db_tgz:               "Kraken2 database as a .tar.zst/.tar.gz archive. Supply this OR the three .k2d files, not both."
        kraken2_hash_k2d:             "hash.k2d of an already unpacked database, e.g. gs://krakendb/hash.k2d. Supply all three .k2d files together."
        ref_genomes_tsv:              "Priority pathogen table: taxid, code, name, colon-joined accessions. No header."
        taxdump_tar_gz:               "NCBI taxdump holding merged.dmp and names.dmp. Leave unset to download the current one. Must be at least as new as the taxonomy behind the reference table."
        name_regex_tsv:               "Optional overrides, one per line: regex <tab> taxid <tab> comment. Tried against the Kraken2 report name only after taxid, merged-taxid and NCBI-synonym matching have all failed."
        min_kraken_reads:             "Reads Kraken2 must assign to a priority taxid (clade total) before its references are aligned. One read by default: whether a priority virus is really there is decided by coverage, not by the classifier."
        min_kraken_reads_discovered:  "The same floor for a virus that is not on the priority list. Higher on purpose: a trace hit on something nobody asked about is not worth a VM."
        min_region_bp:                "Length of one unbroken run of bases at >= min_depth that qualifies a virus on its own. Contiguous, not a total: scattered coverage adding to this does not count."
        min_block_bp:                 "Smallest covered block that makes its region count toward the fallback rule. A block shorter than this is noise, not evidence."
        min_qualifying_regions:       "Regions that must each hold a block of min_block_bp or more for the fallback rule to fire. An unsegmented virus has one region, so the fallback never applies to it."
        max_refs_per_virus:           "Cap on candidate references aligned per virus. Controls cost."
        diamond_db:                   "DIAMOND nr database (.dmnd). Supplying it turns on blastx identification of each detected virus's longest covered block. Unset, the stage is skipped."
        accession2taxid:              "NCBI prot.accession2taxid files, used to turn DIAMOND subject accessions into taxids."
        diamond_top_pct:              "DIAMOND --top: hits within this percent of the best bitscore are kept, and the LCA is taken over them."
        min_match_rank:               "Shallowest Kraken2 rank code a report line may carry and still be matched to a priority taxon. G keeps genus and everything below it and refuses family and above, where a line covers too many viruses to be evidence for one. Set to F to restore the old behaviour."

        discover_min_samples:         "A virus outside the priority table is assembled when it is present in at least this many samples on the flowcell."
        discover_min_reads_one_sample: "...or when a single sample carries at least this many reads for it."
        discover_present_min_reads:   "Reads that make a virus count as present in a sample, for discover_min_samples. One read is enough by default; the read-count branch is what asks for depth."
        discover_max_novel_taxa:      "Cap on discovered viruses, most widespread first. Bounds the cost of a dirty flowcell."

        control_sample_regex:         "Which sample ids are negative controls, matched case-insensitively against the sample id. These are compared against, never reported as findings. Matching nothing is allowed and leaves every row marked 'no control on flowcell'. Spike-in and positive controls are deliberately not in the default: an ERCC well is expected to carry material, so treating it as a background floor would veto real findings."
        min_control_library_reads:    "Smallest control library whose composition is taken as the background. A control sequenced far shallower than the samples cannot say what the background is: on one flowcell a 2,516-read ERCC well put a phage at 5% of its reads, which as reads per million outranks almost every real sample. Such a control is reported but not used to veto calls."
        min_control_support_reads:    "Reads a control must carry for a virus before that virus has a background at all. A background of one or two reads is not a measurement, and letting it set a floor would reject findings on noise."
        min_control_fold:             "How far above the loudest control a sample must sit before its detection is trusted. Below this the call becomes 'suspect: control background'. Measured in reads per million when sample_read_counts_tsv is supplied, in raw reads otherwise."
        min_unique_reads:             "Optional floor on unique mapped reads for any positive call. 0 leaves the coverage rule as the only gate, which is the previous behaviour."
        min_nonviral_pident:          "Identity at or above which a non-viral nr best hit is taken as proof the covered block is host sequence rather than virus. A weak non-viral hit says nothing and is ignored."
        sample_read_counts:           "Optional library sizes, in the same order as cleaned_reads_unaligned_bams, for comparing samples and controls as reads per million. On the flowcell table this is read_counts_depleted. Left empty the comparison uses raw read counts, which flatters a sample whose control library was small. A length that does not match the BAM list is refused rather than silently misaligned."
    }

    # The reference table carries current NCBI taxids, but a Kraken2 database
    # is pinned to the taxonomy of its build date and reports the legacy taxid
    # for anything renumbered since. Lassa (11620 -> 3052310) is one of 44 such
    # entries in Refgenomes_V1.tsv, so matching on the table taxid alone would
    # silently miss them. Resolve the aliases once, up front.
    call resolve_taxid_aliases as resolve_priority {
        input:
            ref_genomes_tsv = ref_genomes_tsv,
            taxdump_tar_gz  = taxdump_tar_gz
    }

    scatter (reads_bam in cleaned_reads_unaligned_bams) {
        call kraken2_classify {
            input:
                reads_bam      = reads_bam,
                kraken2_db_tgz   = kraken2_db_tgz,
                kraken2_hash_k2d = kraken2_hash_k2d,
                kraken2_opts_k2d = kraken2_opts_k2d,
                kraken2_taxo_k2d = kraken2_taxo_k2d,
                memory_mapping   = kraken2_memory_mapping
        }
    }

    # Which non-priority viruses are worth assembling is a question about the
    # flowcell, not about one sample, so it can only be asked once every
    # report exists. The answer is appended to the reference table and the
    # aliases are resolved again over the combined set.
    call discover_extra_taxa {
        input:
            kraken2_reports      = kraken2_classify.kraken2_report,
            ref_genomes_tsv      = ref_genomes_tsv,
            taxid_aliases_tsv    = resolve_priority.taxid_aliases_tsv,
            min_samples          = discover_min_samples,
            min_reads_one_sample = discover_min_reads_one_sample,
            present_min_reads    = discover_present_min_reads,
            max_novel_taxa       = discover_max_novel_taxa,
            ncbi_email           = ncbi_email,
            ncbi_api_key         = ncbi_api_key
    }

    call resolve_taxid_aliases as resolve_all {
        input:
            ref_genomes_tsv = ref_genomes_tsv,
            extra_refs_tsv  = discover_extra_taxa.extra_refs_tsv,
            taxdump_tar_gz  = resolve_priority.taxdump_out
    }

    scatter (pair in zip(cleaned_reads_unaligned_bams,
                         zip(kraken2_classify.sample_id,
                             kraken2_classify.kraken2_report))) {
        File   reads_bam      = pair.left
        String sample_id      = pair.right.left
        File   kraken2_report = pair.right.right

        call select_candidates {
            input:
                sample_id          = sample_id,
                kraken2_report     = kraken2_report,
                ref_genomes_tsv    = resolve_all.combined_refs_tsv,
                taxid_aliases_tsv  = resolve_all.taxid_aliases_tsv,
                name_aliases_tsv   = resolve_all.name_aliases_tsv,
                taxon_levels_tsv   = resolve_all.taxon_levels_tsv,
                name_regex_tsv     = name_regex_tsv,
                min_kraken_reads   = min_kraken_reads,
                min_kraken_reads_discovered = min_kraken_reads_discovered,
                min_kraken_pct     = min_kraken_pct,
                max_refs_per_virus = max_refs_per_virus,
                expand_descendants = expand_descendants,
                max_rows_above_species = max_rows_above_species,
                min_match_rank     = min_match_rank
        }

        call fetch_references {
            input:
                sample_id      = sample_id,
                candidates_tsv = select_candidates.candidates_tsv,
                ncbi_email     = ncbi_email,
                ncbi_api_key   = ncbi_api_key
        }

        scatter (ref_fasta in fetch_references.reference_fastas) {
            call align_and_consensus {
                input:
                    sample_id            = sample_id,
                    reads_bam            = reads_bam,
                    reference_fasta      = ref_fasta,
                    min_depth            = min_depth,
                    min_block_bp         = min_block_bp,
                    min_mapq             = min_mapq,
                    consensus_min_q      = consensus_min_q,
                    consensus_min_freq   = consensus_min_freq,
                    consensus_min_depth  = consensus_min_depth,
                    consensus_ambig_char = consensus_ambig_char,
                    consensus_retry_depth = consensus_retry_depth
            }
        }

        call summarize_sample {
            input:
                sample_id        = sample_id,
                region_stats     = align_and_consensus.region_stats,
                aln_stats        = align_and_consensus.aln_stats,
                consensus_fastas = align_and_consensus.consensus_fasta,
                block_fastas   = align_and_consensus.block_fasta,
                candidates_tsv = select_candidates.candidates_tsv,
                min_depth      = min_depth,
                min_region_bp  = min_region_bp,
                min_block_bp   = min_block_bp,
                min_qualifying_regions = min_qualifying_regions
        }
    }

    call merge_results {
        input:
            summary_tsvs   = summarize_sample.summary_tsv,
            by_region_tsvs = summarize_sample.by_region_tsv
    }

    # One DIAMOND run for the whole flowcell, not one per sample. nr.dmnd is
    # around 370 GB and localizing it dominates the cost of this stage, so the
    # blocks from every sample are pooled and queried together.
    if (defined(diamond_db)) {
        call pool_blocks {
            input:
                block_fastas = summarize_sample.blocks_fasta,
                block_maps   = summarize_sample.block_map_tsv
        }

        call diamond_blastx_nr {
            input:
                contigs      = pool_blocks.pooled_fasta,
                diamond_db   = select_first([diamond_db]),
                top_pct      = diamond_top_pct,
                out_basename = "detected_blocks"
        }

        call assign_taxonomy {
            input:
                hits            = diamond_blastx_nr.hits,
                contig_map      = pool_blocks.pooled_map,
                accession2taxid = accession2taxid,
                taxdump_tar_gz  = diamond_taxdump_tar_gz,
                top_pct         = diamond_top_pct,
                out_basename    = "detected_blocks"
        }

        call annotate_summary {
            input:
                detection_summary_tsv = merge_results.detection_summary_tsv,
                contig_taxonomy       = assign_taxonomy.contig_taxonomy
        }
    }

    # The DIAMOND annotation is optional, so the call step takes whichever
    # summary exists. Without it the blast columns are simply absent and the
    # combined call rests on coverage and the controls alone.
    File summary_before_calls = select_first([annotate_summary.annotated_tsv,
                                              merge_results.detection_summary_tsv])

    call final_calls {
        input:
            detection_summary_tsv = summary_before_calls,
            control_sample_regex  = control_sample_regex,
            min_control_fold      = min_control_fold,
            min_control_library_reads = min_control_library_reads,
            min_control_support_reads = min_control_support_reads,
            min_unique_reads      = min_unique_reads,
            min_nonviral_pident   = min_nonviral_pident,
            sample_ids            = kraken2_classify.sample_id,
            sample_read_counts    = sample_read_counts
    }

    output {
        File         detection_summary_tsv   = final_calls.final_summary_tsv
        File         control_background_tsv  = final_calls.control_background_tsv
        File         detection_summary_precall_tsv = summary_before_calls
        File?        block_blast_hits        = diamond_blastx_nr.hits
        File?        block_taxonomy          = assign_taxonomy.contig_taxonomy
        File?        blocks_fasta            = pool_blocks.pooled_fasta
        File         detection_by_region_tsv = merge_results.detection_by_region_tsv
        Array[File]  kraken2_reports         = kraken2_classify.kraken2_report
        Array[File]  candidate_tables        = select_candidates.candidates_tsv
        Array[File]  consensus_tars          = summarize_sample.consensus_tar
        File         discovered_taxa_tsv     = discover_extra_taxa.discovered_taxa_tsv
        File         combined_refs_tsv       = resolve_all.combined_refs_tsv
    }
}


task kraken2_classify {
    input {
        File    reads_bam

        # Either an archived database, or the three .k2d files of an already
        # unpacked one (as in gs://krakendb/). Supply exactly one form.
        File?   kraken2_db_tgz
        File?   kraken2_hash_k2d
        File?   kraken2_opts_k2d
        File?   kraken2_taxo_k2d

        Boolean memory_mapping = false
        Int     cpu       = 8
        Int     mem_gb    = 120
        Int     disk_gb   = 750
        String  docker    = "quay.io/broadinstitute/viral-classify:2.2.4.0"
    }

    String sample_basename = basename(reads_bam, ".bam")

    command <<<
        set -e -o pipefail

        # Sample id is the leading field of the cleaned BAM name,
        # e.g. 1_01_23_0002.AACJWCVHV.1.cleaned.bam -> 1_01_23_0002
        echo "~{sample_basename}" | cut -d. -f1 > SAMPLE_ID

        mkdir -p kdb
        TGZ="~{default='' kraken2_db_tgz}"
        HASH="~{default='' kraken2_hash_k2d}"

        if [ -n "${HASH}" ]; then
            # Unpacked database: link the three files Kraken2 needs into one
            # directory. Links, not copies, so a 90 GB hash is not duplicated.
            ln -s "${HASH}" kdb/hash.k2d
            ln -s "~{default='' kraken2_opts_k2d}" kdb/opts.k2d
            ln -s "~{default='' kraken2_taxo_k2d}" kdb/taxo.k2d
        elif [ -n "${TGZ}" ]; then
            if [[ "${TGZ}" == *.zst ]]; then
                zstd -d --stdout "${TGZ}" | tar -x -C kdb
            else
                tar -xz -f "${TGZ}" -C kdb
            fi
            # The archive may or may not contain a top-level directory.
            if [ ! -f kdb/hash.k2d ]; then
                inner=$(find kdb -name hash.k2d -print -quit)
                mv "$(dirname "$inner")"/*.k2d kdb/
            fi
        else
            echo "ERROR: supply kraken2_db_tgz, or all three .k2d files" >&2
            exit 1
        fi

        for f in hash.k2d opts.k2d taxo.k2d; do
            if [ ! -s "kdb/$f" ]; then
                echo "ERROR: kdb/$f missing or empty" >&2
                exit 1
            fi
        done

        MMAP=""
        if [ "~{true='yes' false='no' memory_mapping}" = "yes" ]; then
            MMAP="--memory-mapping"
        fi

        samtools fastq -1 r1.fastq -2 r2.fastq -0 /dev/null -s singles.fastq -n "~{reads_bam}"

        if [ ! -s r1.fastq ] && [ ! -s r2.fastq ] && [ ! -s singles.fastq ]; then
            # A barcode that caught no reads leaves a header-only BAM. Kraken2
            # writes no --output file for an empty input, so the gzip below would
            # fail the whole shard. Emit empty files instead: the sample goes
            # through the workflow with zero reads at every taxon, which is the
            # truth about it, rather than taking the whole flowcell down with it.
            echo "WARNING: ~{sample_basename} has no reads; emitting an empty report" >&2
            : > kraken2.report.txt
            : > kraken2.output.txt
        elif [ -s r1.fastq ] && [ -s r2.fastq ]; then
            kraken2 --db kdb --threads ~{cpu} --paired r1.fastq r2.fastq \
                --report kraken2.report.txt --output kraken2.output.txt
        else
            cat r1.fastq singles.fastq > single.fastq
            kraken2 --db kdb --threads ~{cpu} single.fastq \
                --report kraken2.report.txt --output kraken2.output.txt
        fi

        # Kraken2 can also return without writing --output on a near-empty input.
        touch kraken2.report.txt kraken2.output.txt
        gzip -f kraken2.output.txt
    >>>

    output {
        String sample_id       = read_string("SAMPLE_ID")
        File   kraken2_report  = "kraken2.report.txt"
        File   kraken2_output  = "kraken2.output.txt.gz"
    }

    runtime {
        docker:   docker
        cpu:      cpu
        memory:   mem_gb + " GB"
        disks:    "local-disk " + disk_gb + " LOCAL"
        preemptible: 1
        maxRetries: 1
    }
}


task resolve_taxid_aliases {
    input {
        File    ref_genomes_tsv
        File?   extra_refs_tsv
        File?   taxdump_tar_gz
        String  taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        String  docker      = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        # The discovered taxa are appended to the priority table and from here
        # on are treated identically: same alias resolution, same levels, same
        # selection. Only the importance column in column 5 tells them apart.
        cp "~{ref_genomes_tsv}" combined_refs.tsv
        # Refgenomes_V1.tsv ends without a newline. Appending to it blind
        # welds the discovered row onto the HIV row, which silently corrupts
        # an accession rather than failing, so close the line first.
        if [ -s combined_refs.tsv ] && [ -n "$(tail -c1 combined_refs.tsv)" ]; then
            printf '\n' >> combined_refs.tsv
        fi
        EXTRA="~{default='' extra_refs_tsv}"
        if [ -n "${EXTRA}" ] && [ -s "${EXTRA}" ]; then
            cat "${EXTRA}" >> combined_refs.tsv
        fi

        if [ -n "~{default='' taxdump_tar_gz}" ]; then
            cp "~{default='' taxdump_tar_gz}" taxdump.tar.gz
        else
            python3 - <<'PYDL'
import time, urllib.request
url = "~{taxdump_url}"
last = None
for attempt in range(4):
    try:
        urllib.request.urlretrieve(url, "taxdump.tar.gz")
        break
    except Exception as exc:
        last = exc
        time.sleep(2 ** attempt)
else:
    raise last
PYDL
        fi

        tar -xzf taxdump.tar.gz merged.dmp names.dmp nodes.dmp

        python3 <<'PYCODE'
import csv, re

# merged.dmp: "old_taxid | new_taxid |"
old_to_new = {}
with open("merged.dmp") as fh:
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
with open("combined_refs.tsv") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if row and row[0].strip().isdigit():
            t = row[0].strip()
            table_taxids.add(t)
            table_rows_per_taxid[t] = table_rows_per_taxid.get(t, 0) + 1
            if not (len(row) > 4 and row[4].strip() == "discovered"):
                priority_taxids.add(t)

# A discovered virus is reported, not promoted. It is reachable only by the
# taxid Kraken2 already named, which is how it was discovered in the first
# place, and it takes no part in the alias, synonym or ancestor machinery.
# That machinery decides which references a priority virus reaches, and
# whether an ancestor stays under max_rows_above_species, so letting this
# flowcell's incidental finds into it would make the priority results depend
# on what else happened to be sequenced alongside them.
rows = []
for taxid in sorted(table_taxids, key=int):
    # A report may name the table taxid itself,
    rows.append((taxid, taxid, "self"))
    if taxid not in priority_taxids:
        continue
    # a legacy taxid merged into it,
    for old in sorted(new_to_old.get(taxid, []), key=int):
        rows.append((old, taxid, "legacy"))
    # or, if the table ever falls behind, its own successor.
    if taxid in old_to_new:
        rows.append((old_to_new[taxid], taxid, "successor"))

# A Kraken2 database usually carries the species node but not the strain
# nodes the reference table names. Against gs://krakendb, 200 of 332
# reference rows have a taxid that is not a node in the database at all,
# while every one of them has an ancestor that is: Kraken2 can only ever say
# "Zaire ebolavirus", never "Zaire ebolavirus Makona/2014". So a report line
# must also reach the references sitting below it.
parent = {}
rank_of = {}
with open("nodes.dmp") as fh:
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

# How many reference rows each ancestor would pull in. Expanding from a high
# node is ruinous: the Enterovirus genus covers 200 reference rows and
# "Viruses" covers all of them, so select_candidates gates on this count.
rows_under = {}
for taxid in priority_taxids:
    n_rows = table_rows_per_taxid.get(taxid, 1)
    for anc in lineage(taxid)[1:]:
        rows_under.setdefault(anc, {})[taxid] = n_rows

for anc in sorted(rows_under, key=int):
    total = sum(rows_under[anc].values())
    for taxid in sorted(rows_under[anc], key=int):
        rows.append((anc, taxid, "descendant:%d" % total))

# Where each reference taxid sits in the taxonomy. Kraken2 assigns a read to
# the lowest node it can resolve, so a strain-level reference row says nothing
# about how many reads the classifier actually placed on that virus. The
# species node, and the genus above it, are the two levels a Kraken2 report
# reliably carries, so record both and let select_candidates read their counts
# straight out of the report.
# The 2023 ICTV renaming did not always merge the old species id away. For
# EBV, NCBI minted a new species 3050299 "Lymphocryptovirus humangamma4" and
# kept the old 10376 as a no-rank node underneath it, so a database built
# before then reports 10376 and merged.dmp says nothing about the pair. Every
# no-rank node between the reference taxid and its species is a former name
# for that species, so each becomes a fallback the report can be looked up by.
def legacy_ids(node):
    return [node] + sorted(new_to_old.get(node, []), key=int)

levels = {}
need_sci = set()
for taxid in sorted(table_taxids, key=int):
    path = lineage(taxid)
    i_sp = next((i for i, a in enumerate(path)
                 if rank_of.get(a) == "species"), None)
    i_gn = next((i for i, a in enumerate(path)
                 if rank_of.get(a) == "genus"), None)

    def chain(idx, lo):
        """The node at idx, plus the unranked former names below it.

        All of these are the same taxon under different ids, so a count may
        be summed across them: NCBI kept 10376 under the new EBV species
        3050299, and Kraken2 puts most of the reads on the old id.
        """
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
    # The reference taxid itself is a last resort, never part of the sum.
    # Often it is the pre-2023 species node (CMV 10359 is now a no-rank child
    # of species 3050295) and an older database reports nothing else, but
    # where the row really is a strain its reads are not the species' reads.
    levels[taxid] = (sp, sp_chain, gn, gn_chain, taxid)
    need_sci.update(x for x in (sp, gn) if x)

with open("taxid_aliases.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["report_taxid", "table_taxid", "relation"])
    w.writerows(rows)

n_alias = sum(1 for r in rows if r[2] in ("legacy", "successor"))
n_desc = sum(1 for r in rows if r[2].startswith("descendant"))
print(f"{len(table_taxids)} table taxids, {n_alias} alias taxids, "
      f"{n_desc} ancestor-to-reference links")

# names.dmp carries every scientific name, synonym, equivalent name and
# acronym NCBI knows for a taxon. A Kraken2 database may print any of them,
# e.g. "Lassa mammarenavirus" where the table says "Lassa virus", so resolve
# report names through this list rather than by matching text against the
# table's strain names. Those strain names contain place names and strain
# codes ("Yellow fever Angola CNYF01/2016") that collide with unrelated
# organisms, which makes fuzzy matching on them actively dangerous.
# Only the taxon's own identities, never an ancestor's. Including ancestors
# here would make "root" and "Viruses" names for a priority taxon and let the
# name route bypass the rank guard that descendant expansion depends on.
want = {r[0]: r[1] for r in rows
        if r[2] in ("self", "legacy", "successor") and r[1] in priority_taxids}

def norm(s):
    return re.sub(r"[^0-9a-z]+", " ", s.lower()).strip()

names = {}
sci_name = {}
with open("names.dmp", encoding="utf-8", errors="replace") as fh:
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

with open("taxon_levels.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["table_taxid", "species_taxid", "species_name",
                "species_lookup", "genus_taxid", "genus_name",
                "genus_lookup", "fallback_taxid"])
    for taxid in sorted(levels, key=int):
        sp, sp_chain, gn, gn_chain, fb = levels[taxid]
        w.writerow([taxid, sp, sci_name.get(sp, ""), ",".join(sp_chain),
                    gn, sci_name.get(gn, ""), ",".join(gn_chain), fb])

print("%d of %d reference taxids have a species node, %d a genus node"
      % (sum(1 for v in levels.values() if v[0]),
         len(levels),
         sum(1 for v in levels.values() if v[1])))

# A name that points at more than one priority taxon cannot support a call on
# its own ("human parainfluenza virus" covers three). Drop it; an explicit
# name_regex_tsv row can override if a site wants one of them.
ambiguous = sorted(n for n, t in names.items() if len(t) > 1)
with open("name_aliases.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["normalized_name", "table_taxid"])
    for nm in sorted(names):
        if len(names[nm]) == 1:
            w.writerow([nm, next(iter(names[nm]))])

print(f"{len(names) - len(ambiguous)} unambiguous names, "
      f"{len(ambiguous)} dropped as ambiguous")
for nm in ambiguous:
    print(f"  ambiguous, not used: {nm} -> {sorted(names[nm])}")
PYCODE
    >>>

    output {
        File taxid_aliases_tsv = "taxid_aliases.tsv"
        File name_aliases_tsv  = "name_aliases.tsv"
        File taxon_levels_tsv  = "taxon_levels.tsv"
        File combined_refs_tsv = "combined_refs.tsv"
        File taxdump_out       = "taxdump.tar.gz"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "4 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


# Viruses outside the priority table that the flowcell keeps showing. Two ways
# in: seen in many samples, or seen once but with real depth behind it. The
# first catches something circulating, the second a single heavy infection.
# A reference is pulled from RefSeq for each, and from that point on they are
# treated exactly like priority taxa, with virus_importance telling them apart.
task discover_extra_taxa {
    input {
        Array[File] kraken2_reports
        File        ref_genomes_tsv
        File        taxid_aliases_tsv
        Int         min_samples         = 10
        Int         min_reads_one_sample = 1000
        Int         present_min_reads   = 1
        Int         max_novel_taxa      = 50
        String      ncbi_email
        String?     ncbi_api_key
        String      docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        python3 <<'PYCODE'
import csv, json, os, re, sys, time, urllib.parse, urllib.request

MIN_SAMPLES   = ~{min_samples}
MIN_READS_ONE = ~{min_reads_one_sample}
PRESENT_READS = ~{present_min_reads}
MAX_NOVEL     = ~{max_novel_taxa}
EMAIL         = "~{ncbi_email}"
API_KEY       = "~{default='' ncbi_api_key}"
VIRUSES_TAXID = "10239"

reports = [p for p in
           open("~{write_lines(kraken2_reports)}").read().split("\n") if p]

# Anything the priority table already reaches, directly or through an alias or
# an ancestor link, is not a discovery. taxid_aliases.tsv holds exactly that
# set of report taxids, which is why this task runs after the first pass.
known = set()
with open("~{taxid_aliases_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        known.add(row["report_taxid"])
        known.add(row["table_taxid"])
with open("~{ref_genomes_tsv}") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if row and row[0].strip().isdigit():
            known.add(row[0].strip())

def parse_report(path):
    """Species nodes inside the Viruses subtree, with their clade reads.

    The subtree is read off the report's own indentation rather than from the
    taxonomy, so it needs no taxdump and always agrees with the database that
    produced the file.
    """
    out = {}
    viral_indent = None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            if len(f) >= 8 and f[6].strip().isdigit():
                taxid, rank, name_field = f[6].strip(), f[5].strip(), f[7]
            elif f[4].strip().isdigit():
                taxid, rank, name_field = f[4].strip(), f[3].strip(), f[5]
            else:
                continue
            indent = len(name_field) - len(name_field.lstrip(" "))
            if taxid == VIRUSES_TAXID:
                viral_indent = indent
                continue
            if viral_indent is None:
                continue
            if indent <= viral_indent:
                viral_indent = None          # left the Viruses subtree
                continue
            if rank != "S":
                continue
            try:
                out[taxid] = (int(f[1]), name_field.strip())
            except ValueError:
                pass
    return out

n_samples_with, max_reads, names, total_reads = {}, {}, {}, {}
for path in reports:
    if not os.path.exists(path):
        continue
    for taxid, (reads, name) in parse_report(path).items():
        if taxid in known or reads < PRESENT_READS:
            continue
        n_samples_with[taxid] = n_samples_with.get(taxid, 0) + 1
        max_reads[taxid] = max(max_reads.get(taxid, 0), reads)
        total_reads[taxid] = total_reads.get(taxid, 0) + reads
        names[taxid] = name

qualified = []
for taxid, n in n_samples_with.items():
    if n >= MIN_SAMPLES:
        why = "in %d samples" % n
    elif max_reads[taxid] >= MIN_READS_ONE:
        why = "%d reads in one sample" % max_reads[taxid]
    else:
        continue
    qualified.append((taxid, n, max_reads[taxid], total_reads[taxid],
                      names[taxid], why))

# Most widespread first, then deepest, so the cap keeps what matters most.
qualified.sort(key=lambda r: (-r[1], -r[2]))
dropped = qualified[MAX_NOVEL:]
qualified = qualified[:MAX_NOVEL]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def eutils(endpoint, **params):
    params.update(tool="sentinel_refbased_detect", email=EMAIL,
                  retmode="json")
    if API_KEY:
        params["api_key"] = API_KEY
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh)
        except Exception as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise last

def pick_reference(taxid):
    """One RefSeq genome for a taxon, or every segment when it is segmented.

    Restricted to RefSeq so the choice is reproducible: for a virus with
    hundreds of GenBank submissions, "the longest record" is whichever strain
    someone happened to deposit, while RefSeq is one curated genome per
    species. Taxa without a RefSeq genome are reported and skipped rather
    than guessed at.
    """
    term = "txid%s[Organism:exp] AND srcdb_refseq[PROP] AND 1000:1000000[SLEN]" % taxid
    res = eutils("esearch.fcgi", db="nuccore", term=term, retmax=50)
    ids = res.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return [], "no refseq genome"
    time.sleep(0.4)
    summ = eutils("esummary.fcgi", db="nuccore", id=",".join(ids))
    recs = []
    for uid in ids:
        d = summ.get("result", {}).get(uid)
        if not d:
            continue
        recs.append((d.get("accessionversion") or d.get("caption", ""),
                     int(d.get("slen") or 0), d.get("title", "")))
    recs = [r for r in recs if r[0]]
    if not recs:
        return [], "no usable esummary record"

    # A segmented genome is several RefSeq records that differ by segment.
    # Keep the longest record per segment; otherwise keep the longest record.
    seg_of = {}
    for acc, slen, title in recs:
        m = re.search(r"segment\s+([A-Za-z0-9]+)", title, re.I)
        if m:
            key = m.group(1).upper()
            if key not in seg_of or slen > seg_of[key][1]:
                seg_of[key] = (acc, slen)
    if len(seg_of) >= 2:
        accs = [seg_of[k][0] for k in sorted(seg_of)]
        return accs, "refseq, %d segments" % len(accs)
    best = max(recs, key=lambda r: r[1])
    return [best[0]], "refseq %d bp" % best[1]

def code_for(taxid, name):
    tag = re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")[:24]
    return tag or ("taxid" + taxid)

extra_rows, audit = [], []
for taxid, n, mx, tot, name, why in qualified:
    try:
        accs, how = pick_reference(taxid)
    except Exception as exc:
        accs, how = [], "lookup failed: %s" % exc
    time.sleep(0.4)
    audit.append({"taxid": taxid, "name": name, "n_samples": n,
                  "max_reads_one_sample": mx, "total_reads": tot,
                  "qualified_by": why, "reference": ":".join(accs) or "",
                  "reference_source": how,
                  "status": "included" if accs else "skipped"})
    if accs:
        extra_rows.append([taxid, code_for(taxid, name), name,
                           ":".join(accs), "discovered"])

for taxid, n, mx, tot, name, why in dropped:
    audit.append({"taxid": taxid, "name": name, "n_samples": n,
                  "max_reads_one_sample": mx, "total_reads": tot,
                  "qualified_by": why, "reference": "",
                  "reference_source": "",
                  "status": "over max_novel_taxa cap"})

with open("extra_refs.tsv", "w", newline="") as out:
    csv.writer(out, delimiter="\t").writerows(extra_rows)

cols = ["taxid", "name", "n_samples", "max_reads_one_sample", "total_reads",
        "qualified_by", "reference", "reference_source", "status"]
with open("discovered_taxa.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(audit)

print("%d reports, %d non-priority viral species seen, %d qualified, "
      "%d with a reference" % (len(reports), len(n_samples_with),
                               len(qualified) + len(dropped), len(extra_rows)))
for r in audit[:40]:
    print("  %-9s %-44s samples=%-3d max=%-8d %-24s %s"
          % (r["taxid"], r["name"][:44], r["n_samples"],
             r["max_reads_one_sample"], r["qualified_by"], r["status"]))
PYCODE
    >>>

    output {
        File extra_refs_tsv      = "extra_refs.tsv"
        File discovered_taxa_tsv = "discovered_taxa.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "4 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 1
        maxRetries: 2
    }
}

task select_candidates {
    input {
        String  sample_id
        File    kraken2_report
        File    ref_genomes_tsv
        File    taxid_aliases_tsv
        File    name_aliases_tsv
        File    taxon_levels_tsv
        File?   name_regex_tsv
        Int     min_kraken_reads
        Int     min_kraken_reads_discovered
        Float   min_kraken_pct
        Int     max_refs_per_virus
        Boolean expand_descendants
        Int     max_rows_above_species
        String  min_match_rank
        String  docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        python3 <<'PYCODE'
import csv, os, re

MIN_READS = ~{min_kraken_reads}
MIN_READS_DISCOVERED = ~{min_kraken_reads_discovered}
EXPAND_DESCENDANTS      = ~{true='True' false='False' expand_descendants}
MAX_ROWS_ABOVE_SPECIES  = ~{max_rows_above_species}
MIN_MATCH_RANK          = "~{min_match_rank}"
MIN_PCT   = ~{min_kraken_pct}
MAX_REFS  = ~{max_refs_per_virus}

# Priority pathogen table: taxid, code, name, colon-joined accessions.
refs = {}
with open("~{ref_genomes_tsv}") as fh:
    for row in csv.reader(fh, delimiter="\t"):
        if len(row) < 4 or not row[0].strip().isdigit():
            continue
        taxid = row[0].strip()
        code  = (row[1] or "").strip() or "NA"
        name  = (row[2] or "").strip()
        accs  = [a.strip() for a in row[3].split(":") if a.strip()]
        # Column 5 is absent in the priority table as supplied and present
        # only on rows discover_extra_taxa appended.
        importance = (row[4].strip() if len(row) > 4 and row[4].strip()
                      else "priority")
        if accs:
            refs.setdefault(taxid, []).append((code, name, accs, importance))

# report taxid -> table taxid, covering legacy ids the Kraken2 database
# may still use (e.g. 11620 for Lassa, now 3052310).
alias = {}
descendants = {}
with open("~{taxid_aliases_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rel = row["relation"]
        if rel.startswith("descendant"):
            n_rows = int(rel.split(":")[1])
            descendants.setdefault(row["report_taxid"], []).append(
                (row["table_taxid"], n_rows))
        else:
            alias[row["report_taxid"]] = row["table_taxid"]

# normalized NCBI name (scientific name, synonym, equivalent name, acronym)
# -> table taxid, for databases whose taxids do not line up at all.
def norm(s):
    return re.sub(r"[^0-9a-z]+", " ", s.lower()).strip()

name_alias = {}
with open("~{name_aliases_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        name_alias[row["normalized_name"]] = row["table_taxid"]

# table taxid -> the species and genus node above it, so the report can be
# read at those two levels regardless of which node Kraken2 happened to name.
levels = {}
with open("~{taxon_levels_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        levels[row["table_taxid"]] = row

# Optional site-maintained overrides: regex <tab> taxid [<tab> comment].
# Applied to the report name, case-insensitively, when nothing else matched.
# These are explicit on purpose: text matched against the table's own strain
# names produces nonsense ("Yersinia pestis Angola" hitting a yellow fever
# reference named "Yellow fever Angola CNYF01/2016").
regex_rules = []
regex_path = "~{default='' name_regex_tsv}"
if regex_path and os.path.exists(regex_path):
    with open(regex_path) as fh:
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
                regex_rules.append((re.compile(f[0].strip(), re.I), f[1].strip()))
            except re.error as exc:
                raise SystemExit(f"bad regex in name_regex_tsv: {f[0]!r}: {exc}")

ROUTE_RANK = {"taxid": 4, "taxid_alias": 3, "name_synonym": 2, "regex": 1}

def specificity(basis):
    return ROUTE_RANK.get(basis.split(":")[0], 0)

# Kraken2 rank codes, shallowest first, with R1/S1/D2 style suffixes counting
# as fractions of a step. Used to prefer the deepest report line that reaches
# a reference, so a species node beats the realm sitting above it.
RANK_ORDER = "URDKPCOFGS"

def rank_depth(code):
    if not code:
        return -1.0
    base = RANK_ORDER.find(code[0])
    if base < 0:
        return -1.0
    try:
        return base + int(code[1:] or 0) / 10.0
    except ValueError:
        return float(base)

# A report line above genus is not evidence for any single virus. The
# subfamily Alphaherpesvirinae covers HSV-1, HSV-2 and VZV alike, and
# Orthoretrovirinae covers every retrovirus in the database, so routing one
# of them to a panel entry hands that entry a count earned by its siblings.
# On flowcell 261029_VH00635_18_AAF5NT5HV, 309 of 1144 rows matched above
# genus; the best breadth any of them reached was 0.13% and the longest block
# was 49 bp, so none could ever have been detected. They are refused here
# before any reference is fetched for them. Set min_match_rank to "F" to
# restore the previous behaviour.
MIN_MATCH_RANK_DEPTH = rank_depth(MIN_MATCH_RANK)
if MIN_MATCH_RANK_DEPTH < 0:
    raise SystemExit("min_match_rank %r is not a Kraken2 rank code" % MIN_MATCH_RANK)

# Kraken2 report: pct, clade reads, direct reads, rank, taxid, name.
# With --report-minimizer-data two extra columns sit before rank, which
# pushes taxid to index 6.
hits = {}
report_index = {}
with open("~{kraken2_report}") as fh:
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
            pct          = float(f[0])
            clade_reads  = int(f[1])
            direct_reads = int(f[2])
        except ValueError:
            continue
        report_name = f[-1].strip()

        # Every line, keyed by exactly the node Kraken2 named. The species and
        # genus lookups carry their own legacy ids, so nothing is aliased here
        # and the taxid reported back is the one that was really in the file.
        report_index[taxid] = {"rank": rank, "direct": direct_reads,
                               "clade": clade_reads, "pct": pct,
                               "name": report_name}

        # Indexed above whatever its rank, so the species and genus lookups
        # can still read its counts, but only matched when it sits at
        # min_match_rank or below. A code outside RANK_ORDER scores -1 and is
        # kept rather than dropped: an unplaceable rank is not proof the line
        # is too shallow, and silently discarding it would lose real hits.
        this_depth = rank_depth(rank)
        if 0 <= this_depth < MIN_MATCH_RANK_DEPTH:
            continue

        # Resolve the report line onto table taxids, most trustworthy route
        # first: the taxid itself, then a taxid NCBI has merged, then an NCBI
        # name or synonym, then a site-supplied regex.
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

        # Reach references that sit below this node, since Kraken2 names the
        # species while the table often names a strain. Expanding from a
        # species node is safe; expanding from a genus or higher is not, so
        # above species we only allow it when few references hang below.
        if EXPAND_DESCENDANTS and taxid in descendants:
            below_species = rank.startswith("S")
            for table_taxid, n_rows in descendants[taxid]:
                if table_taxid not in refs:
                    continue
                if not below_species and n_rows > MAX_ROWS_ABOVE_SPECIES:
                    continue
                if any(t == table_taxid for t, _ in found):
                    continue
                found.append((table_taxid, "descendant_of:%s" % taxid))

        # Several report lines can fold onto one table taxid. Keep the most
        # specific one rather than the largest clade: the realm Duplodnaviria
        # carries more reads than the EBV species node below it, and crediting
        # EBV with the realm's count is how a number stops meaning anything.
        # Route first, then depth of rank, and only then read count.
        for key, basis in found:
            prev = hits.get(key)
            cur_key = (specificity(basis), rank_depth(rank), clade_reads)
            if prev is None or cur_key > (specificity(prev[3]),
                                          rank_depth(prev[5]), prev[1]):
                hits[key] = (pct, clade_reads, taxid, basis, report_name, rank)

def level_counts(node, lookup_csv, fallback):
    """Reads the report puts on one taxon: direct, and clade.

    The taxon can appear under several ids at once. Kraken2 splits EBV across
    the species 3050299, which it credits with no reads of its own, and the
    former name 10376 below it, which holds 8816. They are one taxon, so the
    direct counts are summed; the clade comes from the highest id present,
    which already contains the rest.

    Empty, not zero, when the report has no line for any of them. Kraken2
    omits nodes it placed no reads on, and a database need not carry the
    species at all.
    """
    ids = [x for x in lookup_csv.split(",") if x] or ([node] if node else [])
    present = [(i, report_index[i]) for i in ids if i in report_index]
    if present:
        direct = sum(e["direct"] for _, e in present)
        return str(direct), str(present[0][1]["clade"]), \
               ",".join(i for i, _ in present)
    e = report_index.get(fallback)
    if e is not None:
        return str(e["direct"]), str(e["clade"]), fallback
    return "", "", ""

rows = []
for taxid, (pct, clade_reads, report_taxid, basis, report_name,
            report_rank) in sorted(hits.items(), key=lambda kv: -kv[1][1]):
    # A priority virus is aligned on a single read, because coverage is what
    # decides whether it is there. A virus nobody put on the list has to earn
    # the VM, so it meets a higher floor in the sample at hand, on top of
    # having already qualified across the flowcell.
    discovered = refs[taxid][0][3] == "discovered"
    floor = MIN_READS_DISCOVERED if discovered else MIN_READS
    if clade_reads < floor or pct < MIN_PCT:
        continue
    lv = levels.get(taxid, {})
    sp_taxid, gn_taxid = lv.get("species_taxid", ""), lv.get("genus_taxid", "")
    fb = lv.get("fallback_taxid", "")
    sp_direct, sp_clade, sp_seen = level_counts(
        sp_taxid, lv.get("species_lookup", ""), fb)
    gn_direct, gn_clade, gn_seen = level_counts(
        gn_taxid, lv.get("genus_lookup", ""), "")
    for idx, (code, name, accs, importance) in enumerate(refs[taxid][:MAX_REFS]):
        rows.append({
            "sample_id":     "~{sample_id}",
            "taxid":         taxid,
            "report_taxid":  report_taxid,
            "report_name":   report_name,
            "report_rank":   report_rank,
            "match_basis":   basis,
            "virus_code":    code,
            "strain_name":   name,
            "virus_importance": importance,
            "ref_index":     str(idx),
            "accessions":    ":".join(accs),
            "n_segments":    str(len(accs)),
            "kraken_reads":  str(clade_reads),
            "kraken_pct":    f"{pct:.4f}",
            "species_taxid":        sp_taxid,
            "species_name":         lv.get("species_name", ""),
            "species_direct_reads": sp_direct,
            "species_clade_reads":  sp_clade,
            "species_report_taxid": sp_seen,
            "genus_taxid":          gn_taxid,
            "genus_name":           lv.get("genus_name", ""),
            "genus_direct_reads":   gn_direct,
            "genus_clade_reads":    gn_clade,
            "genus_report_taxid":   gn_seen,
        })

cols = ["sample_id","taxid","report_taxid","report_name","report_rank","match_basis",
        "virus_code","strain_name","virus_importance","ref_index",
        "accessions","n_segments","kraken_reads","kraken_pct",
        "species_taxid","species_name","species_direct_reads",
        "species_clade_reads","species_report_taxid",
        "genus_taxid","genus_name","genus_direct_reads",
        "genus_clade_reads","genus_report_taxid"]
with open("candidates.tsv","w",newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

basis_counts = {}
for r in rows:
    b = r["match_basis"].split(":")[0]
    basis_counts[b] = basis_counts.get(b, 0) + 1
print(f"{len(rows)} candidate references across "
      f"{len({r['taxid'] for r in rows})} viruses; "
      f"matched by {basis_counts or 'nothing'}")
PYCODE
    >>>

    output {
        File candidates_tsv = "candidates.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "2 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task fetch_references {
    input {
        String  sample_id
        File    candidates_tsv
        String  ncbi_email
        String? ncbi_api_key
        String  docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail
        mkdir -p refs

        python3 <<'PYCODE'
import csv, os, time, urllib.parse, urllib.request

EMAIL   = "~{ncbi_email}"
API_KEY = "~{default='' ncbi_api_key}"
EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def efetch(accessions):
    params = {
        "db":      "nuccore",
        "id":      ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool":    "sentinel_refbased_detect",
        "email":   EMAIL,
    }
    if API_KEY:
        params["api_key"] = API_KEY
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                text = resp.read().decode()
            if text.lstrip().startswith(">"):
                return text
            last = RuntimeError("no FASTA returned: " + text[:200])
        except Exception as exc:
            last = exc
        time.sleep(2 ** attempt)
    raise last

written = 0
with open("~{candidates_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        accs = [a for a in row["accessions"].split(":") if a]
        # One FASTA per candidate reference; segmented viruses keep all
        # their segments in the same file so they are aligned together.
        tag  = "{}__{}__{}".format(row["taxid"], row["virus_code"], row["ref_index"])
        tag  = "".join(c if (c.isalnum() or c in "._-") else "_" for c in tag)
        path = os.path.join("refs", tag + ".fasta")
        text = efetch(accs)
        n_seq = text.count(">")
        if n_seq != len(accs):
            raise RuntimeError(
                f"{tag}: asked for {len(accs)} sequences, NCBI returned {n_seq}")
        with open(path, "w") as out:
            out.write(text if text.endswith("\n") else text + "\n")
        written += 1
        # Stay under the NCBI rate limit: 10/s with a key, 3/s without.
        time.sleep(0.12 if API_KEY else 0.4)

print(f"wrote {written} reference FASTA files")
PYCODE

        # Cromwell needs at least one file to glob even when nothing matched.
        if [ -z "$(ls -A refs)" ]; then
            echo "no candidate references" > refs/.placeholder
        fi
    >>>

    output {
        Array[File] reference_fastas = glob("refs/*.fasta")
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "2 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task align_and_consensus {
    input {
        String  sample_id
        File    reads_bam
        File    reference_fasta
        Int     min_depth
        Int     min_block_bp
        Int     min_mapq
        Int     consensus_min_q
        Float   consensus_min_freq
        Int     consensus_min_depth
        String  consensus_ambig_char
        Int     consensus_retry_depth

        Int     cpu     = 4
        Int     mem_gb  = 16
        Int     disk_gb = 100
        String  docker  = "andersenlabapps/ivar:1.3.1"
    }

    String ref_tag = basename(reference_fasta, ".fasta")

    command <<<
        set -e -o pipefail

        REF_TAG="~{ref_tag}"
        cp "~{reference_fasta}" ref.fasta
        bwa index ref.fasta
        samtools faidx ref.fasta

        samtools fastq -1 r1.fastq -2 r2.fastq -0 /dev/null -s singles.fastq -n "~{reads_bam}"

        if [ -s r1.fastq ] && [ -s r2.fastq ]; then
            bwa mem -t ~{cpu} ref.fasta r1.fastq r2.fastq > aln.sam
        else
            cat r1.fastq singles.fastq > single.fastq
            bwa mem -t ~{cpu} ref.fasta single.fastq > aln.sam
        fi

        # Mark duplicates so "unique reads" excludes PCR and optical copies.
        samtools sort -n -@ ~{cpu} -o namesort.bam aln.sam
        samtools fixmate -m -@ ~{cpu} namesort.bam fixmate.bam
        samtools sort -@ ~{cpu} -o positionsort.bam fixmate.bam
        samtools markdup -@ ~{cpu} positionsort.bam marked.bam
        samtools index marked.bam

        # Unique = mapped, primary, non-supplementary, non-duplicate, above the
        # mapping-quality floor. 0xD04 = unmapped + secondary + duplicate + supplementary.
        samtools view -b -F 0xD04 -q ~{min_mapq} marked.bam > unique.bam
        samtools index unique.bam

        samtools flagstat marked.bam > flagstat.txt
        samtools idxstats unique.bam > idxstats.txt
        cut -f1,2 ref.fasta.fai > ref_lengths.txt

        # Mean MAPQ and unique read count per reference sequence.
        samtools view unique.bam \
            | awk -F'\t' '{ s[$3]+=$5; n[$3]++ }
                          END { for (r in n) printf "%s\t%d\t%.2f\n", r, n[r], s[r]/n[r] }' \
            > mapq.txt

        # Collapse per-base depth into coverage blocks here, so the large
        # depth stream never leaves this task. samtools depth -a emits every
        # position of every reference, in order.
        samtools depth -a -Q 0 unique.bam \
            | awk -F'\t' -v md=~{min_depth} -v mb=~{min_block_bp} '
                function close_run() {
                    if (run > 0) {
                        if (run > longest[cur]) longest[cur] = run
                        if (run >= mb) nblocks[cur]++
                    }
                    run = 0
                }
                { if ($1 != cur) { close_run(); cur = $1 }
                  if ($3 >= md) { run++; covered[cur]++ } else { close_run() } }
                END { close_run()
                      for (r in covered)
                          printf "%s\t%d\t%d\t%d\n", r, covered[r], longest[r], nblocks[r] }' \
            > blocks.txt
        touch blocks.txt

        # Properly paired fraction, over the same unique reads counted above.
        samtools view -c -F 0xD04 marked.bam          > n_mapped.txt
        samtools view -c -F 0xD04 -f 0x2 marked.bam   > n_proper.txt
        samtools view -c -F 0x904 -f 0x400 marked.bam > n_dup.txt

        # ivar consensus is the one step here that is known to die on some
        # inputs: 1.3.x segfaults on particular pileups, and it did so on this
        # workspace's data for EBV and GBV-C while succeeding on the same
        # references in other samples. The consensus is a secondary output, so
        # a crash must not cost us the detection call for the whole flowcell.
        # Try unlimited depth first, then a bounded retry, then give up and
        # record that we did.
        run_consensus() {
            samtools mpileup -aa -A -d "$1" -Q 0 --reference ref.fasta unique.bam \
                | ivar consensus \
                    -p ivar_out \
                    -q ~{consensus_min_q} \
                    -t ~{consensus_min_freq} \
                    -m ~{consensus_min_depth} \
                    -n ~{consensus_ambig_char}
        }

        set +e +o pipefail
        run_consensus 0 > ivar.log 2>&1
        IVAR_RC=$?
        CONSENSUS_STATUS="ok"
        if [ "${IVAR_RC}" -ne 0 ] || [ ! -s ivar_out.fa ]; then
            echo "ivar consensus failed at unlimited depth (rc=${IVAR_RC}); retrying at -d ~{consensus_retry_depth}" >&2
            rm -f ivar_out.fa
            run_consensus ~{consensus_retry_depth} >> ivar.log 2>&1
            IVAR_RC=$?
            if [ "${IVAR_RC}" -ne 0 ] || [ ! -s ivar_out.fa ]; then
                CONSENSUS_STATUS="failed_rc${IVAR_RC}"
            else
                CONSENSUS_STATUS="ok_depth_capped"
            fi
        fi
        set -e -o pipefail

        tail -5 ivar.log >&2 || true

        # ivar writes <prefix>.fa; give it a name unique to this reference so
        # the summarize step can tell the candidates apart. An empty file is a
        # valid outcome, and keeps the output wiring intact.
        if [ -s ivar_out.fa ]; then
            mv ivar_out.fa "${REF_TAG}.consensus.fasta"
        else
            : > "${REF_TAG}.consensus.fasta"
        fi
        echo "${CONSENSUS_STATUS}" > CONSENSUS_STATUS

        # The longest covered block, as the sample's own sequence, for the
        # DIAMOND stage. Blasting the reference here would be circular: it
        # would only confirm what the reference already claims to be. The
        # consensus is called again restricted to the one region that carries
        # the longest block, because ivar names its output after the -p prefix
        # rather than after the reference sequence, so a multi-segment
        # consensus cannot be attributed back to a segment.
        BEST_LINE="$(sort -k3,3nr blocks.txt | head -1 || true)"
        BEST_REGION="$(printf '%s' "${BEST_LINE}" | cut -f1)"
        BEST_LEN="$(printf '%s' "${BEST_LINE}" | cut -f3)"
        : > "${REF_TAG}.block.fasta"

        if [ -n "${BEST_REGION}" ] && [ "${BEST_LEN:-0}" -ge ~{min_block_bp} ]; then
            set +e +o pipefail
            samtools view -b unique.bam "${BEST_REGION}" > region.bam 2>/dev/null
            samtools index region.bam 2>/dev/null
            samtools mpileup -aa -A -d 0 -Q 0 -r "${BEST_REGION}" \
                    --reference ref.fasta region.bam \
                | ivar consensus -p block_out \
                    -q ~{consensus_min_q} -t ~{consensus_min_freq} \
                    -m ~{consensus_min_depth} -n ~{consensus_ambig_char} \
                    > block_ivar.log 2>&1
            BLOCK_RC=$?
            set -e -o pipefail

            if [ "${BLOCK_RC}" -eq 0 ] && [ -s block_out.fa ]; then
                # Longest unbroken run of called bases. ivar writes N wherever
                # depth fell below the floor, so this run is the block itself,
                # with its coordinates taken on the consensus rather than on
                # the reference.
                awk -v tag="${REF_TAG}" -v reg="${BEST_REGION}" -v minb=~{min_block_bp} '
                    /^>/ { next }
                    { seq = seq toupper($0) }
                    END {
                        L = length(seq); run = 0; best = 0; bstart = 0; start = 0
                        for (i = 1; i <= L; i++) {
                            c = substr(seq, i, 1)
                            if (c == "N" || c == "-") {
                                if (run > best) { best = run; bstart = start }
                                run = 0
                            } else {
                                if (run == 0) start = i
                                run++
                            }
                        }
                        if (run > best) { best = run; bstart = start }
                        if (best >= minb)
                            printf(">%s|%s|%d-%d|%d\n%s\n", tag, reg,
                                   bstart, bstart + best - 1, best,
                                   substr(seq, bstart, best))
                    }' block_out.fa > "${REF_TAG}.block.fasta"
            else
                echo "WARNING: block consensus failed for ${REF_TAG} ${BEST_REGION} (rc=${BLOCK_RC})" >&2
                tail -3 block_ivar.log >&2 || true
            fi
        fi

        # One small row per region: tag, region, length, unique reads, mean
        # MAPQ, bases at min_depth, longest block, qualifying blocks. Regions
        # with no coverage still get a row, taken from ref_lengths.txt.
        awk -F'\t' -v tag="${REF_TAG}" '
            FILENAME == "ref_lengths.txt" { len[$1] = $2; order[++k] = $1; next }
            FILENAME == "mapq.txt"        { ur[$1] = $2; mq[$1] = $3; next }
            FILENAME == "blocks.txt"      { cov[$1] = $2; lng[$1] = $3; nb[$1] = $4; next }
            END { for (i = 1; i <= k; i++) { r = order[i]
                      printf "%s\t%s\t%d\t%d\t%.2f\t%d\t%d\t%d\n",
                             tag, r, len[r], ur[r] + 0, mq[r] + 0,
                             cov[r] + 0, lng[r] + 0, nb[r] + 0 } }' \
            ref_lengths.txt mapq.txt blocks.txt > region_stats.tsv

        # ref_tag is <taxid>__<code>__<index>.
        TAXID="${REF_TAG%%__*}"
        REST="${REF_TAG#*__}"
        REF_INDEX="${REST##*__}"
        VIRUS_CODE="${REST%__*}"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${REF_TAG}" "${TAXID}" "${VIRUS_CODE}" "${REF_INDEX}" \
            "$(cat n_mapped.txt)" "$(cat n_proper.txt)" "$(cat n_dup.txt)" \
            "$(cat CONSENSUS_STATUS)" \
            > aln_stats.tsv
    >>>

    output {
        File   region_stats    = "region_stats.tsv"
        File   aln_stats       = "aln_stats.tsv"
        File   consensus_fasta = ref_tag + ".consensus.fasta"
        File   block_fasta     = ref_tag + ".block.fasta"
        File   alignment_bam   = "unique.bam"
        File   flagstat        = "flagstat.txt"
    }

    runtime {
        docker:   docker
        cpu:      cpu
        memory:   mem_gb + " GB"
        # SSD, not HDD. A 100 GB standard persistent disk sustains about
        # 12 MB/s, so a multi-GB read set costs minutes just to be read once,
        # and bwa plus the samtools sort pass over it more than once. On
        # flowcell 261029_VH00635_18_AAF5NT5HV the largest sample's shards took
        # a median of 142 minutes each, against 24 minutes for a whole median
        # sample.
        disks:    "local-disk " + disk_gb + " SSD"
        # Was preemptible 1. Every one of that sample's 88 shards was killed
        # twice before it finished, so three times the compute was spent and
        # thrown away. A shard that runs for hours is worth more than the
        # preemptible discount.
        preemptible: 0
        maxRetries: 2
    }
}


task summarize_sample {
    input {
        String       sample_id
        Array[File]  region_stats
        Array[File]  aln_stats
        Array[File]  consensus_fastas
        Array[File]  block_fastas
        File         candidates_tsv
        Int          min_depth
        Int          min_region_bp
        Int          min_block_bp
        Int          min_qualifying_regions
        String       docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail
        mkdir -p consensus

        python3 <<'PYCODE'
import csv, os, shutil

MIN_DEPTH     = ~{min_depth}
MIN_REGION_BP = ~{min_region_bp}
MIN_BLOCK_BP  = ~{min_block_bp}
MIN_QUALIFYING_REGIONS = ~{min_qualifying_regions}
SAMPLE_ID     = "~{sample_id}"

# Kraken2 support at the two levels a report reliably carries, passed through
# from the candidate table unchanged.
LEVEL_COLS = ["species_taxid","species_name","species_direct_reads",
              "species_clade_reads","species_report_taxid",
              "genus_taxid","genus_name","genus_direct_reads",
              "genus_clade_reads","genus_report_taxid"]

# Paths arrive as files, not as an interpolated string: a joined list would
# break this source the moment a path contained a quote or newline.
def read_paths(list_file):
    with open(list_file) as fh:
        return [l.strip() for l in fh if l.strip()]

region_paths    = read_paths("~{write_lines(region_stats)}")
aln_paths       = read_paths("~{write_lines(aln_stats)}")
consensus_paths = read_paths("~{write_lines(consensus_fastas)}")
block_paths     = read_paths("~{write_lines(block_fastas)}")

# One record per reference that produced a block, keyed by ref_tag. The header
# written upstream is  >ref_tag|region|start-end|length.
blocks = {}
for path in block_paths:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        continue
    with open(path) as fh:
        head, seq = "", []
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

meta = {}
with open("~{candidates_tsv}") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        meta[(row["taxid"], row["ref_index"])] = row

# aln_stats.tsv: ref_tag, taxid, virus_code, ref_index, n_mapped, n_proper, n_dup
aln = {}
for path in aln_paths:
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7:
                aln[f[0]] = {
                    "taxid": f[1], "virus_code": f[2], "ref_index": f[3],
                    "n_mapped": int(f[4]), "n_proper": int(f[5]),
                    "n_duplicate": int(f[6]),
                    "consensus_status": f[7] if len(f) > 7 else "",
                }

# region_stats.tsv: ref_tag, region, length, unique_reads, mean_mapq,
#                   covered_bp, longest_block_bp, n_blocks
by_ref = {}
for path in region_paths:
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 8:
                continue
            by_ref.setdefault(f[0], []).append({
                "region":               f[1],
                "region_length":        int(f[2]),
                "unique_reads":         int(f[3]),
                "mean_mapq":            float(f[4]),
                "covered_bp_min_depth": int(f[5]),
                "longest_block_bp":     int(f[6]),
                "n_blocks":             int(f[7]),
            })

records = []
for ref_tag, region_rows in by_ref.items():
    st = aln.get(ref_tag)
    if st is None:
        continue

    total_len = sum(r["region_length"] for r in region_rows)
    total_cov = sum(r["covered_bp_min_depth"] for r in region_rows)
    n_qualifying_blocks = sum(r["n_blocks"] for r in region_rows)
    # A region qualifies when it holds at least one block of MIN_BLOCK_BP or
    # more at MIN_DEPTH, which is exactly n_blocks >= 1 since the awk that
    # counted them already applied the floor.
    n_qualifying_regions = sum(1 for r in region_rows if r["n_blocks"] >= 1)

    # Detection rule. The primary rule asks for one unbroken stretch of
    # MIN_REGION_BP bases at MIN_DEPTH, not that many bases scattered over
    # the region: a run with no gap in it is a sequence you could call,
    # while the same total spread over hundreds of fragments is not.
    # The fallback counts regions, not blocks, because two blocks a few
    # hundred bases apart on one segment are the same stretch of sequence
    # read twice, while two separate regions are independent evidence. For
    # an unsegmented virus there is only one region, so the fallback cannot
    # fire and the contiguous rule is the only way in.
    best_block_bp  = max((r["longest_block_bp"] for r in region_rows), default=0)
    best_region_bp = max((r["covered_bp_min_depth"] for r in region_rows), default=0)
    if best_block_bp >= MIN_REGION_BP:
        detected, reason = True, (
            f"contiguous {best_block_bp} bp >= {MIN_REGION_BP} bp at {MIN_DEPTH}x")
    elif n_qualifying_regions >= MIN_QUALIFYING_REGIONS:
        detected, reason = True, (
            f"{n_qualifying_regions} regions each with a block >= "
            f"{MIN_BLOCK_BP} bp at {MIN_DEPTH}x")
    else:
        detected, reason = False, (
            f"longest contiguous {best_block_bp} bp ({best_region_bp} bp covered "
            f"in all), {n_qualifying_regions} region(s) with a block >= "
            f"{MIN_BLOCK_BP} bp")

    total_unique = sum(r["unique_reads"] for r in region_rows)
    concordance  = (st["n_proper"] / st["n_mapped"]) if st["n_mapped"] else 0.0
    breadth      = (total_cov / total_len) if total_len else 0.0
    cand         = meta.get((st["taxid"], st["ref_index"]), {})

    records.append({
        "ref_tag":      ref_tag,
        "taxid":        st["taxid"],
        "virus_code":   st["virus_code"],
        "region_rows":  region_rows,
        "total_unique": total_unique,
        # n_mapped, not total_unique: both of these come from marked.bam under
        # -F 0xD04 with no mapping-quality floor, while total_unique is counted
        # off the MAPQ-filtered BAM. Dividing by that would mix two populations
        # and overstate duplication.
        "n_mapped":     st["n_mapped"],
        "n_duplicate":  st["n_duplicate"],
        "breadth":      breadth,
        "concordance":  concordance,
        "detected":     detected,
        "reason":       reason,
        "total_cov":    total_cov,
        "n_blocks":     n_qualifying_blocks,
        "n_qual_regions": n_qualifying_regions,
        "consensus_status": st.get("consensus_status", ""),
        "strain_name":  cand.get("strain_name", ""),
        "virus_importance": cand.get("virus_importance", "priority"),
        "accessions":   cand.get("accessions", ""),
        "report_taxid": cand.get("report_taxid", ""),
        "report_name":  cand.get("report_name", ""),
        "report_rank":  cand.get("report_rank", ""),
        "match_basis":  cand.get("match_basis", ""),
        "kraken_reads": cand.get("kraken_reads", ""),
        "kraken_pct":   cand.get("kraken_pct", ""),
        "levels":       {c: cand.get(c, "") for c in LEVEL_COLS},
    })

# Best reference per virus: most bases covered at MIN_DEPTH, then most reads.
best = {}
for rec in records:
    taxid = rec["taxid"]
    key   = (rec["total_cov"], rec["total_unique"])
    if taxid not in best or key > best[taxid][0]:
        best[taxid] = (key, rec)

summary_cols = ["sample_id","taxid","report_taxid","report_name","report_rank","match_basis",
                "virus_code","strain_name","virus_importance","best_reference",
                "accessions","n_regions","kraken_reads","kraken_pct"] + \
               LEVEL_COLS + \
               ["mapq_by_region","unique_reads_by_region","covered_bp_by_region",
                "total_unique_reads_mapped","n_duplicate_reads",
                "duplicate_fraction","breadth_of_coverage",
                "paired_read_concordance","longest_covered_block_bp",
                "n_covered_blocks","n_qualifying_regions",
                "detection","detection_basis",
                "consensus_status","n_references_compared"]
region_cols  = ["sample_id","taxid","virus_code","virus_importance",
                "best_reference","region",
                "region_length","mean_mapq","unique_reads",
                "covered_bp_min_depth","longest_block_bp","n_blocks",
                "breadth_of_coverage","total_unique_reads_mapped",
                "paired_read_concordance","detection"]

summary_rows, region_out, block_records = [], [], []
for taxid, (_, rec) in sorted(best.items()):
    rows       = rec["region_rows"]
    n_compared = sum(1 for r in records if r["taxid"] == taxid)
    detection  = "detected" if rec["detected"] else "not detected"

    summary_rows.append({
        "sample_id":                 SAMPLE_ID,
        "taxid":                     taxid,
        "report_taxid":              rec["report_taxid"],
        "report_name":               rec["report_name"],
        "report_rank":               rec["report_rank"],
        "match_basis":               rec["match_basis"],
        "virus_code":                rec["virus_code"],
        "strain_name":               rec["strain_name"],
        "virus_importance":          rec["virus_importance"],
        "best_reference":            rec["ref_tag"],
        "accessions":                rec["accessions"],
        "n_regions":                 len(rows),
        "kraken_reads":              rec["kraken_reads"],
        "kraken_pct":                rec["kraken_pct"],
        **rec["levels"],
        "mapq_by_region":            ";".join(f"{r['region']}:{r['mean_mapq']}" for r in rows),
        "unique_reads_by_region":    ";".join(f"{r['region']}:{r['unique_reads']}" for r in rows),
        "covered_bp_by_region":      ";".join(f"{r['region']}:{r['covered_bp_min_depth']}" for r in rows),
        "total_unique_reads_mapped": rec["total_unique"],
        # n_duplicate counts only reads flagged as duplicates, while n_mapped
        # already excludes them, so the two are disjoint and the denominator is
        # their sum, never the kept count on its own.
        "n_duplicate_reads":         rec["n_duplicate"],
        "duplicate_fraction":        (round(rec["n_duplicate"] /
                                            (rec["n_mapped"] + rec["n_duplicate"]), 6)
                                      if (rec["n_mapped"] + rec["n_duplicate"]) else ""),
        "breadth_of_coverage":       round(rec["breadth"], 6),
        "paired_read_concordance":   round(rec["concordance"], 6),
        "longest_covered_block_bp":  max((r["longest_block_bp"] for r in rows), default=0),
        "n_covered_blocks":          sum(r["n_blocks"] for r in rows),
        "n_qualifying_regions":      rec["n_qual_regions"],
        "detection":                 detection,
        "detection_basis":           rec["reason"],
        "consensus_status":          rec["consensus_status"],
        "n_references_compared":     n_compared,
    })

    for r in rows:
        region_out.append({
            "sample_id":                 SAMPLE_ID,
            "taxid":                     taxid,
            "virus_code":                rec["virus_code"],
            "virus_importance":          rec["virus_importance"],
            "best_reference":            rec["ref_tag"],
            "region":                    r["region"],
            "region_length":             r["region_length"],
            "mean_mapq":                 r["mean_mapq"],
            "unique_reads":              r["unique_reads"],
            "covered_bp_min_depth":      r["covered_bp_min_depth"],
            "longest_block_bp":          r["longest_block_bp"],
            "n_blocks":                  r["n_blocks"],
            "breadth_of_coverage":       round(rec["breadth"], 6),
            "total_unique_reads_mapped": rec["total_unique"],
            "paired_read_concordance":   round(rec["concordance"], 6),
            "detection":                 detection,
        })

    # The longest block of the winning reference, for DIAMOND. Only detected
    # viruses are carried forward: a block that failed the detection rule is
    # too short to align meaningfully, and nr is far too expensive to spend on
    # sequence the coverage rule has already rejected.
    blk = blocks.get(rec["ref_tag"])
    if rec["detected"] and blk:
        qid = "|".join([SAMPLE_ID, taxid, rec["ref_tag"],
                        blk["region"], blk["coords"]])
        block_records.append((qid, blk, taxid, rec))

    # Keep the consensus of the winning reference only. Each align task names
    # its output <ref_tag>.consensus.fasta, so match on basename.
    wanted = rec["ref_tag"] + ".consensus.fasta"
    src = next((p for p in consensus_paths
                if p and os.path.basename(p) == wanted), "")
    if src and os.path.exists(src) and os.path.getsize(src) > 0:
        shutil.copyfile(src, os.path.join(
            "consensus", f"{SAMPLE_ID}.{rec['ref_tag']}.consensus.fasta"))

with open("summary.tsv","w",newline="") as out:
    w = csv.DictWriter(out, fieldnames=summary_cols, delimiter="\t")
    w.writeheader(); w.writerows(summary_rows)

with open("by_region.tsv","w",newline="") as out:
    w = csv.DictWriter(out, fieldnames=region_cols, delimiter="\t")
    w.writeheader(); w.writerows(region_out)

# query_id must be whitespace-free: DIAMOND takes qseqid as the first token.
with open("best_blocks.fasta", "w") as out:
    for qid, blk, taxid, rec in block_records:
        out.write(">%s\n" % qid)
        seq = blk["seq"]
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")

with open("block_map.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["query_id", "sample_id", "contig_id", "contig_length"])
    for qid, blk, taxid, rec in block_records:
        w.writerow([qid, SAMPLE_ID, qid, blk["length"]])

print(f"{SAMPLE_ID}: {len(summary_rows)} viruses, "
      f"{sum(1 for r in summary_rows if r['detection']=='detected')} detected, "
      f"{len(block_records)} blocks for DIAMOND")
PYCODE

        tar -czf consensus.tar.gz consensus
    >>>

    output {
        File summary_tsv   = "summary.tsv"
        File by_region_tsv = "by_region.tsv"
        File consensus_tar = "consensus.tar.gz"
        File blocks_fasta  = "best_blocks.fasta"
        File block_map_tsv = "block_map.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "8 GB"
        disks:    "local-disk 50 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task merge_results {
    input {
        Array[File] summary_tsvs
        Array[File] by_region_tsvs
        String      docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        python3 <<'PYCODE'
def merge(paths, out_path):
    header_written = False
    n = 0
    with open(out_path, "w") as out:
        for path in paths:
            if not path:
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

def paths(list_file):
    with open(list_file) as fh:
        return [l.strip() for l in fh if l.strip()]

s = merge(paths("~{write_lines(summary_tsvs)}"),   "detection_summary.tsv")
r = merge(paths("~{write_lines(by_region_tsvs)}"), "detection_by_region.tsv")
print(f"{s} summary rows, {r} region rows")
PYCODE
    >>>

    output {
        File detection_summary_tsv   = "detection_summary.tsv"
        File detection_by_region_tsv = "detection_by_region.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "4 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task pool_blocks {
    input {
        Array[File] block_fastas
        Array[File] block_maps
        String      docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        python3 <<'PYCODE'
import csv, os

def read_paths(list_file):
    with open(list_file) as fh:
        return [l.strip() for l in fh if l.strip()]

fastas = read_paths("~{write_lines(block_fastas)}")
maps   = read_paths("~{write_lines(block_maps)}")

n = 0
with open("pooled_blocks.fasta", "w") as out:
    for path in fastas:
        if not path or not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
                out.write(line)

rows, seen = [], set()
for path in maps:
    if not path or not os.path.exists(path):
        continue
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["query_id"] in seen:
                continue
            seen.add(row["query_id"])
            rows.append(row)

with open("pooled_map.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["query_id", "sample_id", "contig_id", "contig_length"])
    for r in rows:
        w.writerow([r["query_id"], r["sample_id"], r["contig_id"],
                    r["contig_length"]])

print(f"{n} blocks pooled from {len(fastas)} samples, {len(rows)} map rows")
PYCODE
    >>>

    output {
        File pooled_fasta = "pooled_blocks.fasta"
        File pooled_map   = "pooled_map.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "4 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task annotate_summary {
    input {
        File   detection_summary_tsv
        File   contig_taxonomy
        String docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        python3 <<'PYCODE'
import csv, gzip

# contig_id is the query id written by summarize_sample:
#   sample_id|taxid|ref_tag|region|start-end
tax = {}
with gzip.open("~{contig_taxonomy}", "rt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        parts = r["contig_id"].split("|")
        if len(parts) >= 2:
            tax[(r["sample_id"], parts[1])] = r

NEW = ["blast_block_bp", "blast_region", "blast_n_hits", "blast_n_hits_lca",
       "blast_best_subject", "blast_best_subject_taxon",
       "blast_best_subject_category", "blast_pct_identity",
       "blast_bitscore", "diamond_taxid", "diamond_taxon", "diamond_rank",
       "diamond_category", "diamond_species", "diamond_genus",
       "diamond_family", "blast_agreement"]

# An LCA of root or cellular organisms means the hits spread too widely for a
# call. That is an absence of information, not a contradiction, so it is
# reported as unresolved and never as a disagreement.
VAGUE = {"", "root", "cellular organisms", "Viruses", "unclassified"}

def norm(s):
    return (s or "").strip().lower()

def agreement(row, t):
    if t is None:
        return "no nr hit"
    name = (t.get("taxon_name") or "").strip()
    if name in VAGUE:
        return "unresolved (%s)" % (name or "no LCA")
    d_sp, d_gn = norm(t.get("species")), norm(t.get("genus"))
    k_sp, k_gn = norm(row.get("species_name")), norm(row.get("genus_name"))
    if d_sp and k_sp and d_sp == k_sp:
        return "species match"
    if d_gn and k_gn and d_gn == k_gn:
        return "genus match"
    if norm(t.get("category")) != "virus":
        return "non-viral hit (%s)" % (t.get("category") or "?")
    return "different taxon"

with open("~{detection_summary_tsv}") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
cols = (list(rows[0].keys()) if rows else []) + NEW

n_hit = 0
for row in rows:
    t = tax.get((row["sample_id"], row["taxid"]))
    if t is not None:
        n_hit += 1
    parts = (t or {}).get("contig_id", "").split("|")
    row.update({
        "blast_block_bp":           (t or {}).get("contig_length", ""),
        "blast_region":             parts[3] if len(parts) > 3 else "",
        "blast_n_hits":             (t or {}).get("n_hits", ""),
        "blast_n_hits_lca":         (t or {}).get("n_hits_lca", ""),
        "blast_best_subject":       (t or {}).get("best_hit", ""),
        "blast_best_subject_taxon": (t or {}).get("best_hit_taxon", ""),
        "blast_best_subject_category": (t or {}).get("best_hit_category", ""),
        "blast_pct_identity":       (t or {}).get("best_pident", ""),
        "blast_bitscore":           (t or {}).get("best_bitscore", ""),
        "diamond_taxid":            (t or {}).get("taxid", ""),
        "diamond_taxon":            (t or {}).get("taxon_name", ""),
        "diamond_rank":             (t or {}).get("rank", ""),
        "diamond_category":         (t or {}).get("category", ""),
        "diamond_species":          (t or {}).get("species", ""),
        "diamond_genus":            (t or {}).get("genus", ""),
        "diamond_family":           (t or {}).get("family", ""),
        "blast_agreement":          agreement(row, t),
    })

with open("detection_summary_annotated.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t",
                       lineterminator="\n")
    w.writeheader(); w.writerows(rows)

print("%d of %d summary rows carry a DIAMOND result" % (n_hit, len(rows)))
PYCODE
    >>>

    output {
        File annotated_tsv = "detection_summary_annotated.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "8 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


# Negative-control comparison and the combined call.
#
# Everything this task needs is already in the merged summary: the controls are
# ordinary samples on the flowcell and went through Kraken2, alignment and the
# coverage rule exactly like the rest. So the comparison is pure table work on
# one small file, with no extra alignment and no second pass over any BAM.
#
# Two things are decided here that the coverage rule alone cannot decide:
#
#   1. Whether a virus is distinguishable from the flowcell's own background.
#      Index hopping and reagent contamination put a few hundred reads of a
#      real virus into the water control, and a sample carrying no more than
#      the control does is not evidence of infection.
#   2. Whether the covered block is actually viral. A retroviral reference
#      carries a host oncogene, so reads from the host genome pile onto that
#      part of the reference and clear the contiguity rule. The nr best hit
#      names the host outright, while the LCA over the top hits collapses to
#      root and reports nothing.
task final_calls {
    input {
        File    detection_summary_tsv
        String  control_sample_regex
        Float   min_control_fold
        Int          min_control_library_reads
        Int          min_control_support_reads
        Int          min_unique_reads
        Float        min_nonviral_pident
        Array[String] sample_ids
        Array[Int]   sample_read_counts
        String       docker = "python:3.11-slim"
    }

    command <<<
        set -e -o pipefail

        # Through the environment, not interpolated into the source: a regex is
        # full of quotes and backslashes, and one of them inside the heredoc
        # would rewrite this program rather than fail.
        export CONTROL_REGEX='~{control_sample_regex}'

        python3 <<'PYCODE'
import csv, os, re, sys

try:
    CONTROL_RE = re.compile(os.environ["CONTROL_REGEX"], re.I)
except re.error as exc:
    sys.exit("control_sample_regex is not a valid regex: %s" % exc)
MIN_FOLD     = ~{min_control_fold}
MIN_CONTROL_LIBRARY = ~{min_control_library_reads}
MIN_CONTROL_SUPPORT = ~{min_control_support_reads}
MIN_UNIQUE   = ~{min_unique_reads}
MIN_NONVIRAL_PIDENT = ~{min_nonviral_pident}

with open("~{detection_summary_tsv}") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

if not rows:
    sys.exit("detection summary is empty")

# Optional library sizes, so a 1M-read control is not compared head to head
# with a 20M-read sample. Without them the comparison falls back to raw reads,
# which is conservative in the wrong direction: a small control makes a sample
# look cleaner than it is.
def read_lines(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]

ids    = read_lines("~{write_lines(sample_ids)}")
counts = read_lines("~{write_lines(sample_read_counts)}")

lib = {}
if counts:
    if len(counts) != len(ids):
        sys.exit("sample_read_counts has %d entries but there are %d samples; "
                 "they are paired by position, so a mismatch would attach the "
                 "wrong library size to every sample after the first gap"
                 % (len(counts), len(ids)))
    for sid, n in zip(ids, counts):
        try:
            if int(n) > 0:
                lib[sid] = int(n)
        except ValueError:
            pass
NORMALIZED = bool(lib)

def reads(row):
    try:
        return int(row.get("total_unique_reads_mapped") or 0)
    except ValueError:
        return 0

def signal(row):
    """Reads per million of library, or raw reads when sizes are unknown."""
    n = reads(row)
    total = lib.get(row["sample_id"], 0)
    return (1e6 * n / total) if (NORMALIZED and total) else float(n)

samples  = sorted({r["sample_id"] for r in rows})
controls = {s for s in samples if CONTROL_RE.search(s)}
if not controls:
    print("WARNING: no sample matched %r; every row is reported as "
          "'no control on flowcell'" % CONTROL_RE.pattern, file=sys.stderr)

# A control sequenced far shallower than the samples is still a control: its
# rows are never reported as findings. It just cannot say what the background
# is, because one stray read in a tiny library is a large share of it.
shallow = set()
if NORMALIZED:
    shallow = {s for s in controls if lib.get(s, 0) < MIN_CONTROL_LIBRARY}
    for s in sorted(shallow):
        print("WARNING: control %s has %d reads, under min_control_library_reads"
              "=%d; reported but not used as background"
              % (s, lib.get(s, 0), MIN_CONTROL_LIBRARY), file=sys.stderr)
background_controls = controls - shallow

print("%d samples, %d controls: %s"
      % (len(samples), len(controls), ", ".join(sorted(controls)) or "none"))
print("%d usable as background: %s"
      % (len(background_controls), ", ".join(sorted(background_controls)) or "none"))

# Background per virus: the loudest control on the flowcell. Taken over every
# control row, detected or not, because a control that falls short of the
# coverage rule still shows how much of that virus is moving between wells.
bg = {}
for r in rows:
    if r["sample_id"] not in background_controls:
        continue
    # A background resting on a couple of reads is noise, and letting it set a
    # floor would reject real findings on the strength of it.
    if reads(r) < MIN_CONTROL_SUPPORT:
        continue
    key = r["taxid"]
    cur = bg.get(key)
    s = signal(r)
    if cur is None or s > cur["signal"]:
        bg[key] = {"signal": s, "reads": reads(r), "sample": r["sample_id"],
                   "detection": r.get("detection", ""),
                   "block": r.get("longest_covered_block_bp", "")}

NONVIRAL = {"Bacteria", "Archaea", "Human", "Vertebrate", "Fungi", "Protozoa",
            "Invertebrate", "Plant", "Other eukaryote", "Other"}

def pident(row):
    try:
        return float(row.get("blast_pct_identity") or 0)
    except ValueError:
        return 0.0

NEW = ["is_control", "control_signal", "control_reads", "control_sample",
       "control_fold", "control_status", "final_call", "final_call_reason"]

n_by_call = {}
for r in rows:
    is_control = r["sample_id"] in controls
    b = bg.get(r["taxid"])
    s = signal(r)

    # Fold over the control. Blank for a control row, and for a virus no
    # control carried: a background of zero gives no fold, and calling it
    # infinite would be a number that only looks like a measurement.
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
    elif fold != "" and fold < MIN_FOLD:
        control_status = "at or near control background"
    else:
        control_status = "above control background"

    detected = r.get("detection") == "detected"
    cat  = (r.get("blast_best_subject_category") or "").strip()
    taxon = (r.get("blast_best_subject_taxon") or "").strip()
    agree = (r.get("blast_agreement") or "").strip()

    # Ordered by which error does the most damage to a surveillance call.
    # Contamination first: a virus indistinguishable from the water control is
    # not reportable however well nr confirms its identity, because what nr
    # confirms is the identity of the contaminant.
    if is_control:
        call, why = "control sample", "negative or process control"
    elif not detected:
        call, why = "not detected", r.get("detection_basis", "")
    elif control_status == "at or near control background":
        call = "suspect: control background"
        why = ("%.3g vs %.3g in %s (%.1fx, under %gx)"
               % (s, b["signal"], b["sample"], fold, MIN_FOLD))
    elif MIN_UNIQUE and reads(r) < MIN_UNIQUE:
        call = "suspect: too few reads"
        why = "%d unique reads, under min_unique_reads=%d" % (reads(r), MIN_UNIQUE)
    elif cat in NONVIRAL and pident(r) >= MIN_NONVIRAL_PIDENT:
        call = "suspect: non-viral block"
        why = ("best nr hit is %s (%s) at %.1f%% identity: the covered block "
               "looks like host sequence, not virus" % (taxon or "?", cat, pident(r)))
    elif agree == "different taxon":
        call = "suspect: conflicting taxon"
        why = "nr names %s where the reference claims %s" % (
            r.get("diamond_taxon") or "?", r.get("species_name") or r.get("strain_name") or "?")
    elif agree in ("species match", "genus match"):
        call = "confirmed"
        why = "coverage rule plus nr %s (%s)" % (agree, r.get("diamond_taxon") or "?")
    else:
        call = "detected, unconfirmed"
        why = "coverage rule only; nr gave %s" % (agree or "no result")

    if control_status.startswith("no control") and call not in (
            "not detected", "control sample"):
        why += " [%s]" % control_status

    r.update({
        "is_control":     "yes" if is_control else "no",
        "control_signal": ("%.4f" % s) if NORMALIZED else str(reads(r)),
        "control_reads":  str(b["reads"]) if b else "",
        "control_sample": b["sample"] if b else "",
        "control_fold":   ("%.2f" % fold) if fold != "" else "",
        "control_status": control_status,
        "final_call":     call,
        "final_call_reason": why,
    })
    n_by_call[call] = n_by_call.get(call, 0) + 1

cols = [c for c in rows[0].keys() if c not in NEW] + NEW
with open("detection_summary_final.tsv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)

# The background itself, so a reviewer can see what the controls carried
# without reading the whole summary.
with open("control_background.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["taxid", "virus_code", "strain_name", "control_sample",
                "control_reads", "control_signal", "control_detection",
                "control_longest_block_bp", "n_detected_at_background"])
    name = {}
    for r in rows:
        name.setdefault(r["taxid"], (r.get("virus_code", ""), r.get("strain_name", "")))
    for taxid, b in sorted(bg.items(), key=lambda kv: -kv[1]["signal"]):
        # Only calls this background actually threatens: a virus absent from a
        # sample is trivially at background, and counting those would bury the
        # number that matters.
        near = sum(1 for r in rows
                   if r["taxid"] == taxid
                   and r["final_call"] == "suspect: control background")
        code, strain = name.get(taxid, ("", ""))
        w.writerow([taxid, code, strain, b["sample"], b["reads"],
                    "%.4f" % b["signal"], b["detection"], b["block"], near])

print("normalization: %s" % ("reads per million" if NORMALIZED
                             else "raw read counts (no library sizes given)"))
for call in sorted(n_by_call, key=lambda c: -n_by_call[c]):
    print("  %-30s %d" % (call, n_by_call[call]))
PYCODE
    >>>

    output {
        File final_summary_tsv      = "detection_summary_final.tsv"
        File control_background_tsv = "control_background.tsv"
    }

    runtime {
        docker:   docker
        cpu:      1
        memory:   "4 GB"
        disks:    "local-disk 20 HDD"
        preemptible: 2
        maxRetries: 2
    }
}


task diamond_blastx_nr {
    input {
        File   contigs
        File   diamond_db
        Float  top_pct
        String out_basename
        String diamond_mode = "--mid-sensitive"
        Float  evalue       = 0.00001
        Int    block_size   = 6
        Int    index_chunks = 1
        Int    cpu          = 32
        Int    mem_gb       = 128
        String docker       = "quay.io/biocontainers/diamond:2.1.11--h5ca1c30_2"
    }
    Int disk_gb = ceil(size(diamond_db, "GB") + 20 * size(contigs, "GB")) + 100

    command <<<
        set -euo pipefail
        diamond version
        mkdir -p tmp
        if [ ! -s ~{contigs} ]; then
          : > "~{out_basename}.nr_hits.tsv"
        else
          diamond blastx \
            --db ~{diamond_db} \
            --query ~{contigs} \
            --out "~{out_basename}.nr_hits.tsv" \
            --outfmt 6 qseqid qlen sseqid pident length qstart qend evalue bitscore \
            ~{diamond_mode} \
            --top ~{top_pct} --max-hsps 1 \
            --evalue ~{evalue} \
            --block-size ~{block_size} --index-chunks ~{index_chunks} \
            --threads ~{cpu} --tmpdir tmp
        fi
        gzip "~{out_basename}.nr_hits.tsv"
    >>>

    output {
        File hits = "~{out_basename}.nr_hits.tsv.gz"
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: 0
        maxRetries: 0
    }
}

task assign_taxonomy {
    input {
        File        hits
        File        contig_map
        Array[File] accession2taxid
        File?       taxdump_tar_gz
        Float       top_pct
        String      out_basename
        String      ncbi = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
        Int         mem_gb = 32
        String      docker = "python:3.11"
    }
    Int disk_gb = ceil(size(accession2taxid, "GB") + 5 * size(hits, "GB")) + 60

    command <<<
        set -euo pipefail
        cat > assign_taxonomy.py <<'EOF_PY'
#!/usr/bin/env python3
"""Assign a taxon to each contig from DIAMOND blastx hits (LCA of top hits) and
build a per-sample long table with contig counts / bp at every rank.

usage: assign_taxonomy.py HITS.tsv[.gz] CONTIG_MAP.tsv ACC_TAXID.tsv TAXDUMP_DIR TOP_PCT OUT_PREFIX

HITS columns (diamond --outfmt 6): qseqid qlen sseqid pident length qstart qend evalue bitscore
ACC_TAXID: accession.version <tab> taxid (only accessions that were hit)
"""
import csv
import gzip
import os
import sys
from collections import defaultdict

LINEAGE = ["domain", "phylum", "class", "order", "family", "genus", "species"]
DOMAIN_RANKS = {"superkingdom", "domain", "acellular root"}
PROTOZOA = {"Apicomplexa", "Euglenozoa", "Metamonada", "Amoebozoa", "Sar", "Discoba",
            "Parabasalia", "Cryptophyceae", "Haptista", "Ciliophora"}
INVERT = {"Nematoda", "Platyhelminthes", "Arthropoda", "Mollusca", "Annelida"}


def load_taxdump(d):
    parent, rank, name = {}, {}, {}
    with open(os.path.join(d, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            t = int(f[0])
            parent[t], rank[t] = int(f[1]), f[2]
    with open(os.path.join(d, "names.dmp")) as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if f[3] == "scientific name":
                name[int(f[0])] = f[1]
    merged = {}
    p = os.path.join(d, "merged.dmp")
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip("\t|\n").split("\t|\t")
                merged[int(f[0])] = int(f[1])
    return parent, rank, name, merged


class Taxonomy:
    def __init__(self, d):
        self.parent, self.rank, self.name, self.merged = load_taxdump(d)
        self._path = {}

    def resolve(self, t):
        t = self.merged.get(t, t)
        return t if t in self.parent else None

    def path(self, t):
        """root -> t"""
        if t not in self._path:
            p, x = [], t
            while True:
                p.append(x)
                if x == 1 or self.parent[x] == x:
                    break
                x = self.parent[x]
            self._path[t] = p[::-1]
        return self._path[t]

    def lca(self, taxids):
        paths = [self.path(t) for t in taxids]
        out = 1
        for level in zip(*paths):
            if all(x == level[0] for x in level):
                out = level[0]
            else:
                break
        return out

    def lineage(self, t):
        cols = {c: "" for c in LINEAGE}
        if not t:
            return cols
        for x in self.path(t):
            r = self.rank[x]
            if r in DOMAIN_RANKS and not cols["domain"]:
                cols["domain"] = self.name.get(x, "")
            elif r in cols and r != "domain":
                cols[r] = self.name.get(x, "")
        return cols

    def category(self, t):
        if not t:
            return "Unclassified"
        names = {self.name.get(x, "") for x in self.path(t)}
        nm = self.name.get(t, "")
        if nm in ("root", "cellular organisms"):
            return "Root"
        for key, cat in (("Viruses", "Virus"), ("Bacteria", "Bacteria"), ("Archaea", "Archaea")):
            if key in names:
                return cat
        if "Eukaryota" in names:
            if "Homo sapiens" in names:
                return "Human"
            if "Chordata" in names:
                return "Vertebrate"
            if "Fungi" in names:
                return "Fungi"
            if names & PROTOZOA:
                return "Protozoa"
            if names & INVERT:
                return "Invertebrate"
            if "Viridiplantae" in names:
                return "Plant"
            return "Other eukaryote"
        return "Other"


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    hits_f, map_f, acc_f, tax_d, top_pct, prefix = sys.argv[1:7]
    top = float(top_pct) / 100.0
    tax = Taxonomy(tax_d)

    acc2tax = {}
    with open(acc_f) as fh:
        for line in fh:
            a, t = line.rstrip("\n").split("\t")[:2]
            acc2tax[a] = int(t)

    # best-scoring hits per query
    hits = defaultdict(list)  # q -> [(bitscore, sseqid, pident)]
    with open_any(hits_f) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            hits[f[0]].append((float(f[8]), f[2], float(f[3])))

    contigs = []
    with open(map_f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            contigs.append(r)

    ctab = open(f"{prefix}.contig_taxonomy.tsv", "w")
    cw = csv.writer(ctab, delimiter="\t", lineterminator="\n")
    # best_hit_category is what the single best-scoring subject is, independent
    # of the LCA. A retroviral block that carries a host oncogene hits host
    # protein at 99% identity while the LCA over the top hits collapses to
    # root, so the LCA alone reports "unresolved" where the best hit is saying
    # plainly that this sequence is not viral.
    cw.writerow(["sample_id", "contig_id", "contig_length", "n_hits", "n_hits_lca",
                 "best_hit", "best_hit_taxid", "best_hit_taxon", "best_pident", "best_bitscore",
                 "best_hit_category",
                 "taxid", "taxon_name", "rank", "category"] + LINEAGE)

    # per sample aggregation: taxid -> [contigs_taxon, bp_taxon]
    agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    sample_bp = defaultdict(int)
    no_taxid_acc = set()

    for c in contigs:
        q, sid, clen = c["query_id"], c["sample_id"], int(c["contig_length"])
        sample_bp[sid] += clen
        h = sorted(hits.get(q, []), reverse=True)
        best = h[0] if h else None
        tids, used = [], 0
        if h:
            cutoff = best[0] * (1 - top)
            for bs, s, _ in h:
                if bs < cutoff:
                    break
                t = acc2tax.get(s)
                t = tax.resolve(t) if t else None
                if t is None:
                    no_taxid_acc.add(s)
                    continue
                tids.append(t)
                used += 1
        lca = tax.lca(tids) if tids else 0
        bt = tax.resolve(acc2tax.get(best[1], 0)) if best else None
        agg[sid][lca][0] += 1
        agg[sid][lca][1] += clen
        lin = tax.lineage(lca)
        cw.writerow([sid, c["contig_id"], clen, len(h), used,
                     best[1] if best else "", bt or "", tax.name.get(bt, "") if bt else "",
                     best[2] if best else "", best[0] if best else "",
                     tax.category(bt) if bt else "",
                     lca, tax.name.get(lca, "unclassified") if lca else "unclassified",
                     tax.rank.get(lca, "unclassified") if lca else "unclassified",
                     tax.category(lca)] + [lin[k] for k in LINEAGE])
    ctab.close()

    with open(f"{prefix}.diamond_long.tsv", "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "percentage_bp", "contigs_clade", "contigs_taxon", "bp_clade",
                    "bp_taxon", "taxid", "taxon_name", "category", "rank"] + LINEAGE)
        for sid in sorted(agg):
            clade = defaultdict(lambda: [0, 0])
            for t, (n, bp) in agg[sid].items():
                for x in (tax.path(t) if t else [0]):
                    clade[x][0] += n
                    clade[x][1] += bp
            total = sample_bp[sid] or 1
            # depth-first order like a kraken report
            children = defaultdict(list)
            for t in clade:
                if t not in (0, 1):
                    children[tax.parent[t]].append(t)
            order = [0] if 0 in clade else []
            stack = [1] if 1 in clade else []
            while stack:
                t = stack.pop()
                order.append(t)
                stack.extend(sorted(children[t], key=lambda x: clade[x][1]))
            for t in order:
                n, bp = clade[t]
                own = agg[sid].get(t, [0, 0])
                lin = tax.lineage(t)
                w.writerow([sid, round(100.0 * bp / total, 4), n, own[0], bp, own[1], t,
                            tax.name.get(t, "unclassified") if t else "unclassified",
                            tax.category(t), tax.rank.get(t, "unclassified") if t else "unclassified"]
                           + [lin[k] for k in LINEAGE])

    print(f"{len(contigs)} contigs, {sum(1 for c in contigs if c['query_id'] in hits)} with hits, "
          f"{len(no_taxid_acc)} hit accessions without taxid", file=sys.stderr)


if __name__ == "__main__":
    main()
EOF_PY

        # taxonomy dump
        mkdir -p taxdump
        if [ -n "~{taxdump_tar_gz}" ]; then
          tar -xzf "~{taxdump_tar_gz}" -C taxdump
        else
          curl -fsSL --retry 5 -o taxdump.tar.gz "~{ncbi}/taxdump.tar.gz"
          tar -xzf taxdump.tar.gz -C taxdump nodes.dmp names.dmp merged.dmp
        fi

        # accession -> taxid, only for accessions that were hit
        gzip -dc ~{hits} | cut -f3 | sort -u > hit_accs.txt
        echo "$(wc -l < hit_accs.txt) distinct hit accessions"
        A2T=( ~{sep=" " accession2taxid} )
        if [ ${#A2T[@]} -eq 0 ]; then
          for f in prot pdb dead_prot; do
            curl -fsSL --retry 5 -C - -o "$f.accession2taxid.gz" "~{ncbi}/accession2taxid/$f.accession2taxid.gz"
            A2T+=( "$f.accession2taxid.gz" )
          done
        fi
        : > acc_taxid.tsv
        for f in "${A2T[@]}"; do
          awk -F'\t' 'NR==FNR{h[$1]; next} ($2 in h){print $2"\t"$3}' hit_accs.txt <(gzip -dcf "$f") >> acc_taxid.tsv
        done
        echo "$(cut -f1 acc_taxid.tsv | sort -u | wc -l) accessions mapped to taxids"

        python3 assign_taxonomy.py ~{hits} ~{contig_map} acc_taxid.tsv taxdump ~{top_pct} "~{out_basename}"
        gzip -f "~{out_basename}.contig_taxonomy.tsv"
        # Field 12 is taxid: best_hit_category was inserted at field 11.
        awk -F'\t' 'NR>1 && $12 != "0"' <(gzip -dc "~{out_basename}.contig_taxonomy.tsv.gz") | wc -l > n_classified.txt
    >>>

    output {
        File contig_taxonomy = "~{out_basename}.contig_taxonomy.tsv.gz"
        File long_table      = "~{out_basename}.diamond_long.tsv"
        Int  n_classified    = read_int("n_classified.txt")
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu: 2
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: 1
    }
}

