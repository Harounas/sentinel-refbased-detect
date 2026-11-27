version 1.0

# Build the decoy panel that sentinel_refbased_detect's evidence review maps
# against, on Terra rather than on a laptop.
#
# The panel answers one question for a priority pathogen that Kraken2 named but
# that failed the coverage rule: can the few reads be attributed to this taxon?
# That needs two kinds of genome, and the work divides naturally:
#
#   sibling  one RefSeq genome per OTHER species of the pathogen's genus. These
#            are the genomes a read from a neighbouring virus would prefer, so
#            they are what makes the question answerable. One shard per genus,
#            because a genus is a single NCBI query however many of its species
#            are on the priority list: 288 taxa collapse to about 50 shards.
#   strain   genomes of the pathogen itself, clustered to representatives. One
#            shard per species, about 83 of them, and the slow ones: a species
#            with hundreds of deposits is downloaded and sketched here.
#
# Divergence is the thing this has to get right. Lassa's clades sit around 65%
# identical to each other, and SARS-CoV-2, MERS and the seasonal coronaviruses
# are separate species entirely. So:
#   * each segment is clustered on its own, because an L segment and an S
#     segment share almost nothing and clustering them together keeps nothing;
#   * representatives are chosen to span the diversity, so the cap trims the
#     crowded middle of a clade rather than deleting a clade outright;
#   * the panel records the lineage spread it measured, so a species that needs
#     more representatives says so in its own rows.

workflow build_pathogen_panel {
    input {
        File          ref_genomes_tsv
        String        ncbi_email
        String?       ncbi_api_key
        File?         taxdump_tar_gz

        Array[String] kinds               = ["sibling", "strain"]
        Float         min_identity        = 0.99
        Float         coverage_identity   = 0.90
        Int           max_candidates      = 300
        Int           max_representatives = 12
        Int           min_cluster_members = 2

        # Install the same code the detection workflow runs, by tag or branch.
        String        package_url =
            "https://github.com/Harounas/sentinel-refbased-detect/archive/refs/heads/main.tar.gz"
        String        docker      = "python:3.11-slim"
    }

    parameter_meta {
        ref_genomes_tsv:     "Priority pathogen table, the same Refgenomes_V15.tsv the detection workflow uses."
        ncbi_email:          "Contact address sent to NCBI E-utilities, as they ask of automated clients."
        taxdump_tar_gz:      "NCBI taxdump. Left unset it is downloaded once and passed to every shard, rather than each shard fetching its own copy."
        kinds:               "sibling, strain, or both. Siblings alone take minutes; strains dominate the cost."
        min_identity:        "Identity above which two genomes of one species are the same thing and collapse to one representative. Redundancy removal only."
        coverage_identity:   "Every genome of a species must end up at least this similar to some representative. This, not a fixed cluster threshold, is what decides how many a species keeps: a near-clonal species such as SARS-CoV-2 stops at one, while Lassa keeps many."
        min_cluster_members: "Genomes a cluster needs before it may be chosen as a diversity representative. Guards against a lone misannotated deposit, which is exactly what furthest-first selection reaches for."
        max_representatives: "Representatives kept per species and per segment, chosen to span the lineage spread rather than taken in order."
        max_candidates:      "Genomes considered per species before clustering."
    }

    call list_jobs {
        input:
            ref_genomes_tsv = ref_genomes_tsv,
            kinds           = kinds,
            taxdump_tar_gz  = taxdump_tar_gz,
            ncbi_email      = ncbi_email,
            ncbi_api_key    = ncbi_api_key,
            package_url     = package_url,
            docker          = docker
    }

    scatter (job in read_tsv(list_jobs.job_rows)) {
        call build_shard {
            input:
                kind                = job[0],
                job_id              = job[1],
                job_name            = job[2],
                guards              = job[3],
                ref_genomes_tsv     = ref_genomes_tsv,
                taxdump_tar_gz      = list_jobs.taxdump_out,
                ncbi_email          = ncbi_email,
                ncbi_api_key        = ncbi_api_key,
                min_identity        = min_identity,
                coverage_identity   = coverage_identity,
                max_candidates      = max_candidates,
                max_representatives = max_representatives,
                min_cluster_members = min_cluster_members,
                package_url         = package_url,
                docker              = docker
        }
    }

    call merge_panel {
        input:
            shards      = build_shard.panel_rows,
            package_url = package_url,
            docker      = docker
    }

    output {
        File panel_tsv     = merge_panel.panel_tsv
        File diversity_tsv = merge_panel.diversity_tsv
        File empty_shards  = merge_panel.empty_shards
        File jobs_tsv      = list_jobs.jobs_tsv
    }
}


task list_jobs {
    input {
        File          ref_genomes_tsv
        Array[String] kinds
        File?         taxdump_tar_gz
        String        ncbi_email
        String?       ncbi_api_key
        String        package_url
        String        docker
    }

    command <<<
        set -e -o pipefail
        pip install --quiet "~{package_url}"

        sentinel-detect panel jobs \
            --ref-tsv "~{ref_genomes_tsv}" \
            --kinds ~{sep=' ' kinds} \
            ~{if defined(taxdump_tar_gz) then "--taxdump " + taxdump_tar_gz else ""} \
            --cache-dir cache \
            --ncbi-email "~{ncbi_email}" \
            ~{if defined(ncbi_api_key) then "--ncbi-api-key " + ncbi_api_key else ""} \
            -o jobs.tsv

        # Headerless, because a WDL scatter over read_tsv has no way to skip
        # the header row and would otherwise run a shard called "kind".
        tail -n +2 jobs.tsv > job_rows.tsv

        # The taxdump every shard needs, fetched once here instead of ~130
        # times. It is the largest download in the workflow.
        cp cache/taxonomy/taxdump.tar.gz taxdump.tar.gz
    >>>

    output {
        File jobs_tsv    = "jobs.tsv"
        File job_rows    = "job_rows.tsv"
        File taxdump_out = "taxdump.tar.gz"
    }

    runtime {
        docker:      docker
        cpu:         1
        memory:      "4 GB"
        disks:       "local-disk 20 HDD"
        preemptible: 2
        maxRetries:  2
    }
}


task build_shard {
    input {
        String  kind
        String  job_id
        String  job_name
        String  guards
        File    ref_genomes_tsv
        File    taxdump_tar_gz
        String  ncbi_email
        String? ncbi_api_key
        Float   min_identity
        Float   coverage_identity
        Int     max_candidates
        Int     max_representatives
        Int     min_cluster_members
        String  package_url
        String  docker
    }

    command <<<
        set -e -o pipefail
        pip install --quiet "~{package_url}"
        echo "shard: ~{kind} ~{job_id} ~{job_name}"

        # A shard that cannot be built must not take the panel down with it:
        # one genus whose NCBI query fails is a gap in the decoys, while a
        # failed workflow is no panel at all. The empty file keeps the scatter
        # whole and merge_panel simply finds nothing in it.
        if ! sentinel-detect panel shard \
                --kind "~{kind}" \
                --id "~{job_id}" \
                --guards "~{guards}" \
                --ref-tsv "~{ref_genomes_tsv}" \
                --taxdump "~{taxdump_tar_gz}" \
                --cache-dir cache \
                --ncbi-email "~{ncbi_email}" \
                ~{if defined(ncbi_api_key) then "--ncbi-api-key " + ncbi_api_key else ""} \
                --identity ~{min_identity} \
                --coverage-identity ~{coverage_identity} \
                --max-candidates ~{max_candidates} \
                --max-representatives ~{max_representatives} \
                --min-cluster-members ~{min_cluster_members} \
                -o shard.tsv 2> shard.log; then
            echo "WARNING: shard ~{kind} ~{job_id} (~{job_name}) failed" >&2
            tail -20 shard.log >&2 || true
            printf 'taxid\tcode\tname\taccessions\trole\trelated_to\tnote\n' > shard.tsv
        fi

        tail -5 shard.log >&2 || true
        wc -l < shard.tsv
    >>>

    output {
        File panel_rows = "shard.tsv"
        File shard_log  = "shard.log"
    }

    runtime {
        docker:      docker
        cpu:         1
        # A species with hundreds of genomes is sketched in memory; viral
        # genomes are small, so this is generous rather than tight.
        memory:      "8 GB"
        disks:       "local-disk 30 HDD"
        preemptible: 2
        maxRetries:  2
    }
}


task merge_panel {
    input {
        Array[File] shards
        String      package_url
        String      docker
    }

    command <<<
        set -e -o pipefail
        pip install --quiet "~{package_url}"

        sentinel-detect panel merge \
            --shards ~{sep=' ' shards} \
            -o decoy_panel.tsv

        # What the clustering found, species by species: the lineage spread is
        # the number that says whether one genome can stand for a species or
        # whether it needs several. A species sitting near the bottom of this
        # table is one to look at before trusting a single representative.
        python3 - ~{sep=' ' shards} <<'PYCODE'
import csv, re

rows = list(csv.DictReader(open("decoy_panel.tsv"), delimiter="\t"))
spread = {}
for r in rows:
    if r["role"] != "strain":
        continue
    note = r["note"] or ""
    seg = re.match(r"(segment \w+|genome)", note)
    key = (r["taxid"], r["name"], seg.group(1) if seg else "")
    e = spread.setdefault(key, {"representatives": 0, "covered": "", "capped": "no",
                                "lowest": "", "median": ""})
    e["representatives"] += 1
    m = re.search(r"covers every genome within (\d+)%", note)
    if m:
        e["covered"] = m.group(1)
    if "CAPPED" in note:
        e["capped"] = "yes"
    m = re.search(r"lineage spread (\d+)%-(\d+)%", note)
    if m:
        e["lowest"], e["median"] = m.group(1), m.group(2)

# Sorted by the coverage actually achieved, worst first. This is the column to
# read: it is what the panel can promise for that species, measured rather
# than assumed. A row marked capped ran out of representatives before it
# reached the target, and wants a higher max_representatives.
with open("panel_diversity.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["taxid", "name", "part", "representatives",
                "genomes_covered_within_pct", "hit_cap",
                "lowest_identity_pct", "median_identity_pct"])
    for (taxid, name, part), e in sorted(
            spread.items(), key=lambda kv: int(kv[1]["covered"] or 100)):
        w.writerow([taxid, name, part, e["representatives"], e["covered"],
                    e["capped"], e["lowest"], e["median"]])

capped = [k for k, e in spread.items() if e["capped"] == "yes"]
if capped:
    print("%d species/segment(s) hit the representative cap: %s"
          % (len(capped), ", ".join("%s %s" % (k[1], k[2]) for k in capped[:8])))

n_sib = sum(1 for r in rows if r["role"] == "sibling")
print("%d panel rows: %d sibling, %d strain, over %d species"
      % (len(rows), n_sib, len(rows) - n_sib, len({r["taxid"] for r in rows})))

# A shard that failed wrote a header and nothing else, so the scatter stays
# whole. That is right for keeping the run alive and wrong for believing the
# result: count them here, where the number is impossible to miss.
import os, sys
empty = [p for p in sys.argv[1:] if os.path.getsize(p) < 80]
if empty:
    print("WARNING: %d of %d shards produced no rows. Their logs say why; a run "
          "with many empty shards is not a complete panel."
          % (len(empty), len(sys.argv) - 1))
with open("empty_shards.txt", "w") as fh:
    fh.write("%d of %d shards empty\n" % (len(empty), len(sys.argv) - 1))
    for p in empty:
        fh.write(p + "\n")
PYCODE
    >>>

    output {
        File panel_tsv     = "decoy_panel.tsv"
        File diversity_tsv = "panel_diversity.tsv"
        File empty_shards  = "empty_shards.txt"
    }

    runtime {
        docker:      docker
        cpu:         1
        memory:      "4 GB"
        disks:       "local-disk 20 HDD"
        preemptible: 2
        maxRetries:  2
    }
}
