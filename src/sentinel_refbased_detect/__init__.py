"""Local port of the Terra workflow warn-id/sentinel_refbased_detect (config V15).

Reference-based detection of priority pathogens in sentinel NGS data:

1. Kraken2 taxonomic classification of each sample.
2. Keep taxids from the priority-pathogen table that Kraken2 reports, plus
   viruses outside it that recur across the run or carry real depth in one sample.
3. Fetch every candidate reference for those taxids from NCBI.
4. Align reads to each candidate reference (bwa), call a consensus (ivar).
5. Pick the best reference per virus and apply the coverage detection rule.
6. Optionally identify each detected block with DIAMOND blastx against nr.
7. Compare against negative controls and make the combined final call.
"""

__version__ = "0.15.0"
WDL_METHOD = "agora://warn-id/sentinel_refbased_detect/19"
