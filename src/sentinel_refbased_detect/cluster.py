"""MinHash sketching and greedy clustering, to drop redundant strain genomes.

A priority pathogen can have hundreds of near-identical GenBank genomes. Putting
all of them in a panel costs index size and, worse, splits a read's evidence
across copies of the same thing. Clustering keeps one representative per group
of genomes that are alike above a chosen identity.

The sketch is a bottom-k MinHash over canonical k-mers and the distance is
Mash's, so the identity threshold means what it does in Mash: no external tool,
no extra dependency, and the arithmetic is the published one.
"""

import hashlib
import math

# k=16, not the 21 Mash uses for bacteria. A k-mer survives a substitution
# only if none of its k bases is hit, so with k=21 two sequences 76% identical
# share 0.76**21, under half a percent, of their k-mers: the sketch saturates
# and every divergent pair reads as 0% identity. Lassa's clades sit near that
# range, so the one number that should tell us a species needs several
# representatives would be the one number we could not measure. k=16 keeps
# resolution down to roughly 70% identity, which covers viral lineages.
K = 16
SKETCH_SIZE = 1000
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def canonical_kmers(seq, k=K):
    seq = seq.upper()
    rc = seq.translate(COMPLEMENT)[::-1]
    n = len(seq)
    for i in range(n - k + 1):
        fwd = seq[i:i + k]
        if "N" in fwd:
            continue
        rev = rc[n - k - i:n - i]
        yield fwd if fwd <= rev else rev


def sketch(seq, k=K, size=SKETCH_SIZE):
    """Bottom-`size` hashes of the canonical k-mers, as a sorted tuple."""
    seen = set()
    for kmer in canonical_kmers(seq, k):
        seen.add(int.from_bytes(hashlib.blake2b(kmer.encode(), digest_size=8).digest(), "big"))
    return tuple(sorted(seen)[:size])


def jaccard(a, b, size=SKETCH_SIZE):
    """Jaccard of two bottom-k sketches, estimated over their merged bottom-k."""
    if not a or not b:
        return 0.0
    merged = sorted(set(a) | set(b))[:size]
    if not merged:
        return 0.0
    sa, sb = set(a), set(b)
    shared = sum(1 for h in merged if h in sa and h in sb)
    return shared / len(merged)


def mash_identity(j, k=K):
    """Mash distance turned into an identity, as 1 - D."""
    if j <= 0:
        return 0.0
    if j >= 1:
        return 1.0
    return 1.0 + math.log(2 * j / (1 + j)) / k


def cluster_sequences(records, min_identity=0.95, k=K, size=SKETCH_SIZE):
    """Greedy clustering, longest genome first.

    `records` is an iterable of (name, sequence). Returns a list of
    {"representative", "length", "members"}: the longest genome of each group
    represents it, which keeps the most complete sequence rather than whichever
    happened to be deposited first.
    """
    items = []
    for name, seq in records:
        seq = "".join(seq.split())
        items.append({"name": name, "length": len(seq), "sketch": sketch(seq, k, size)})
    items.sort(key=lambda r: (-r["length"], r["name"]))

    clusters = []
    for item in items:
        for c in clusters:
            j = jaccard(c["sketch"], item["sketch"], size)
            if mash_identity(j, k) >= min_identity:
                c["members"].append(item["name"])
                break
        else:
            clusters.append({"representative": item["name"], "length": item["length"],
                             "sketch": item["sketch"], "members": [item["name"]]})
    return clusters


def pairwise_identity(clusters, k=K, size=SKETCH_SIZE):
    """Min and median identity between cluster representatives.

    A species whose lineages sit far apart, as Lassa's clades do, shows a low
    minimum here. That is the number that says whether one genome can stand
    for the species or whether several must.
    """
    sk = [c["sketch"] for c in clusters if c.get("sketch")]
    if len(sk) < 2:
        return 1.0, 1.0
    ids = [mash_identity(jaccard(sk[i], sk[j], size), k)
           for i in range(len(sk)) for j in range(i + 1, len(sk))]
    ids.sort()
    return ids[0], ids[len(ids) // 2]


def select_diverse(clusters, max_n, k=K, size=SKETCH_SIZE):
    """Pick at most max_n representatives that span the diversity.

    Taking the first max_n clusters would keep whichever genomes happen to be
    longest and could drop a whole lineage: for Lassa that could mean losing
    clade III entirely. This takes the longest genome first, then repeatedly
    the representative furthest from everything already chosen, so the kept set
    covers the spread rather than one corner of it.
    """
    if len(clusters) <= max_n:
        return list(clusters)
    remaining = list(clusters)
    chosen = [remaining.pop(0)]            # already sorted longest first
    while remaining and len(chosen) < max_n:
        best, best_dist = None, None
        for c in remaining:
            d = min(1.0 - mash_identity(jaccard(c["sketch"], x["sketch"], size), k)
                    for x in chosen)
            if best_dist is None or d > best_dist:
                best, best_dist = c, d
        chosen.append(best)
        remaining.remove(best)
    return chosen


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        yield name, "".join(seq)
