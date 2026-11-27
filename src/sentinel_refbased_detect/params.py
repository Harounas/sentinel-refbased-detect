"""Pipeline parameters. Defaults are those of the Terra method configuration
warn-id/sentinel_refbased_detect_0109_V15 (method version 19)."""

from dataclasses import dataclass, field, fields


def _p(default, help, group="detection"):
    return field(default=default, metadata={"help": help, "group": group})


@dataclass
class Params:
    # --- Kraken2 / candidate selection -------------------------------------
    min_kraken_reads: int = _p(1, "Clade reads Kraken2 must assign to a priority taxid before its references are aligned.", "selection")
    min_kraken_reads_discovered: int = _p(100, "Same floor for a discovered (non-priority) virus.", "selection")
    min_kraken_pct: float = _p(0.0, "Minimum Kraken2 percentage for a candidate.", "selection")
    max_refs_per_virus: int = _p(20, "Cap on candidate references aligned per virus.", "selection")
    expand_descendants: bool = _p(True, "Let a report line reach reference rows below it in the taxonomy.", "selection")
    max_rows_above_species: int = _p(8, "Expansion from above species only when at most this many reference rows hang below.", "selection")
    min_match_rank: str = _p("G", "Shallowest Kraken2 rank code that may be matched to a priority taxon (G = genus).", "selection")

    # --- discovery of non-priority viruses ---------------------------------
    discover_min_samples: int = _p(10, "Assemble a non-priority virus present in at least this many samples.", "discovery")
    discover_min_reads_one_sample: int = _p(1000, "...or carrying at least this many reads in one sample.", "discovery")
    discover_present_min_reads: int = _p(1, "Reads that make a virus count as present in a sample.", "discovery")
    discover_max_novel_taxa: int = _p(50, "Cap on discovered viruses, most widespread first.", "discovery")

    # --- alignment / detection rule ----------------------------------------
    min_depth: int = _p(2, "Depth a base needs to count as covered.")
    min_region_bp: int = _p(1000, "One unbroken covered run of this length detects a virus on its own.")
    min_block_bp: int = _p(200, "Smallest covered block that makes its region count toward the fallback rule.")
    min_qualifying_regions: int = _p(2, "Regions that must each hold a block >= min_block_bp for the fallback rule.")
    min_mapq: int = _p(20, "Mapping-quality floor for unique reads.")

    # --- consensus ----------------------------------------------------------
    consensus_min_q: int = _p(20, "ivar consensus -q.", "consensus")
    consensus_min_freq: float = _p(0.75, "ivar consensus -t.", "consensus")
    consensus_min_depth: int = _p(2, "ivar consensus -m.", "consensus")
    consensus_ambig_char: str = _p("N", "ivar consensus -n.", "consensus")
    consensus_retry_depth: int = _p(10000, "mpileup -d used when ivar fails at unlimited depth.", "consensus")

    # --- negative controls / final call -------------------------------------
    control_sample_regex: str = _p(r"^(nec|ntc|neg|negative|blank|water|h2o)([-_.]|$)", "Case-insensitive regex identifying negative-control sample ids.", "controls")
    min_control_fold: float = _p(10.0, "How far above the loudest control a sample must sit.", "controls")
    min_control_library_reads: int = _p(50000, "Smallest control library used as background.", "controls")
    min_control_support_reads: int = _p(10, "Reads a control must carry for a virus before it sets a background.", "controls")
    min_unique_reads: int = _p(0, "Optional floor on unique mapped reads for any positive call.", "controls")
    min_nonviral_pident: float = _p(80.0, "Identity at/above which a non-viral nr best hit marks the block as host.", "controls")

    # --- evidence review of failing priority rows ----------------------------
    evidence_min_kraken_reads: int = _p(1, "Kraken2 clade reads a failing priority row needs before it is reviewed.", "evidence")
    evidence_min_identity: float = _p(95.0, "Mean percent identity of mapped reads to the reference for the identity criterion to pass.", "evidence")
    evidence_min_mean_mapq: float = _p(30.0, "Mean MAPQ of the mapped reads for the mapping-quality criterion to pass.", "evidence")
    evidence_max_multimapping: float = _p(20.0, "Percent of reads carrying an alternative alignment above which mapping quality fails.", "evidence")
    evidence_bin_bp: int = _p(500, "Window size used to judge how the reads are spread over the genome.", "evidence")
    evidence_min_bins: int = _p(3, "Windows that must carry reads for the distribution criterion to pass.", "evidence")
    evidence_max_busiest_bin: float = _p(80.0, "Percent of reads allowed in a single window before distribution fails.", "evidence")
    evidence_min_breadth: float = _p(0.001, "Fraction of the reference covered for the coverage criterion to pass.", "evidence")
    evidence_hopping_ratio: float = _p(100.0, "How many times louder another sample must be before cross-talk is a plausible explanation.", "evidence")
    evidence_min_reads_for_specificity: int = _p(2, "Mapped reads needed before the competing-reference check is run.", "evidence")
    evidence_min_retained: float = _p(80.0, "Percent of reads that must be attributable to the target taxon once related genomes compete. Reads that fit several taxa equally are shared out in proportion to the reads that fit only one.", "evidence")
    evidence_score_margin: int = _p(0, "Alignment-score margin within which a second taxon counts as fitting a read equally well. 0 means an exact tie; raise it to treat near-ties as ambiguous too.", "evidence")
    evidence_em_rounds: int = _p(20, "Rounds used to apportion ambiguous reads between taxa.", "evidence")
    evidence_min_criteria: int = _p(4, "Criteria that must pass for a verdict of 'equivocal' rather than 'not supported'. One failure out of enough passes is equivocal; two or more is not supported.", "evidence")

    # --- DIAMOND -------------------------------------------------------------
    diamond_top_pct: float = _p(10.0, "DIAMOND --top; LCA is taken over hits within this percent of the best.", "diamond")
    diamond_mode: str = _p("--mid-sensitive", "DIAMOND sensitivity flag.", "diamond")
    diamond_evalue: float = _p(0.00001, "DIAMOND --evalue.", "diamond")
    diamond_block_size: float = _p(6, "DIAMOND --block-size.", "diamond")
    diamond_index_chunks: int = _p(1, "DIAMOND --index-chunks.", "diamond")


def param_fields():
    return fields(Params)
