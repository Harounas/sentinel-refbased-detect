version 1.0

# Build the background layer and leave it in a bucket, indexed and ready.
#
# The sibling and strain layers are both derived from the priority table, and
# that is their limit: they can say a read fits Lassa better than its
# neighbours, but a read from a virus nobody listed, from the host, or from a
# reagent organism has nowhere to go and is credited to whatever reference it
# was aligned against. At ten or twenty reads that is the commonest way a false
# signal survives.
#
# So this layer is assembled from NCBI directly and never consults the table:
# every RefSeq viral genome (14,844 species, about 580 Mb, taken from the
# aggregated release files rather than 15,865 separate assembly downloads),
# the host genome, and the usual reagent contaminants.
#
# It is built once. The index is 3.6 Gb and is far too large to rebuild per
# pathogen, so the output belongs in a bucket beside the Kraken2 and nr
# databases, and every later run points at it.

workflow build_background_db {
    input {
        String  ncbi_email
        String? ncbi_api_key

        # GRCh38 or T2T. Optional only in the sense that the workflow runs
        # without it; leaving it out means host reads stay unexplained.
        File?   host_fasta
        File?   taxdump_tar_gz
        Array[File] accession2taxid = []

        Boolean include_viral        = true
        Boolean include_contaminants = true

        Int     cpu     = 8
        Int     mem_gb  = 32
        Int     disk_gb = 200
        String  package_url =
            "https://github.com/Harounas/sentinel-refbased-detect/archive/refs/heads/main.tar.gz"
        String  docker = "quay.io/biocontainers/bwa:0.7.18--he4a0461_1"
    }

    parameter_meta {
        host_fasta:      "Host genome, e.g. gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta. Most reads that cannot be explained by any virus are host, so leaving this out blunts the layer considerably."
        accession2taxid: "nucl_gb.accession2taxid.gz, used to name each background sequence. Downloaded if not supplied; it is large, so passing a copy already in your bucket is faster."
        include_viral:   "Every RefSeq viral genome. This is what gives a read from an unlisted virus somewhere true to go."
    }

    call build {
        input:
            ncbi_email           = ncbi_email,
            ncbi_api_key         = ncbi_api_key,
            host_fasta           = host_fasta,
            taxdump_tar_gz       = taxdump_tar_gz,
            accession2taxid      = accession2taxid,
            include_viral        = include_viral,
            include_contaminants = include_contaminants,
            cpu                  = cpu,
            mem_gb               = mem_gb,
            disk_gb              = disk_gb,
            package_url          = package_url,
            docker               = docker
    }

    output {
        # The whole directory as one archive: bwa needs the .amb/.ann/.bwt/.pac
        # /.sa files beside the FASTA, and shipping them separately invites a
        # half-copied index that fails only when a run tries to use it.
        File background_tar  = build.background_tar
        File background_taxa = build.background_taxa
        File build_log       = build.build_log
    }
}


task build {
    input {
        String      ncbi_email
        String?     ncbi_api_key
        File?       host_fasta
        File?       taxdump_tar_gz
        Array[File] accession2taxid
        Boolean     include_viral
        Boolean     include_contaminants
        Int         cpu
        Int         mem_gb
        Int         disk_gb
        String      package_url
        String      docker
    }

    command <<<
        set -e -o pipefail

        # The bwa image has no pip; install into a virtualenv on top of it so
        # bwa itself stays the version this index will be read back with.
        apt-get update -qq && apt-get install -y -qq python3 python3-venv curl
        python3 -m venv /opt/venv
        /opt/venv/bin/pip install --quiet "~{package_url}"

        /opt/venv/bin/sentinel-detect background build \
            --out-dir background \
            --ncbi-email "~{ncbi_email}" \
            ~{if defined(ncbi_api_key) then "--ncbi-api-key " + ncbi_api_key else ""} \
            ~{if defined(host_fasta) then "--host " + host_fasta else ""} \
            ~{if defined(taxdump_tar_gz) then "--taxdump " + taxdump_tar_gz else ""} \
            ~{if length(accession2taxid) > 0 then "--accession2taxid " else ""}~{sep=' ' accession2taxid} \
            ~{if include_viral then "" else "--no-viral"} \
            ~{if include_contaminants then "" else "--no-contaminants"} \
            --threads ~{cpu} \
            --verbose 2>&1 | tee build.log

        ls -la background/
        cp background/background.taxa.tsv .
        cp background/pipeline.log build.log || true

        # The FASTA is 3.6 Gb and compresses poorly; the index files are what
        # matter on the far side, so archive the lot without compression.
        tar -cf background.tar -C background .
        ls -la background.tar
    >>>

    output {
        File background_tar  = "background.tar"
        File background_taxa = "background.taxa.tsv"
        File build_log       = "build.log"
    }

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      mem_gb + " GB"
        disks:       "local-disk " + disk_gb + " SSD"
        # Hours of indexing thrown away by a preemption costs more than the
        # discount saves, and this runs once.
        preemptible: 0
        maxRetries:  1
    }
}
