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


def select_diverse(clusters, max_n, k=K, size=SKETCH_SIZE, min_members=2,
                   coverage_identity=None):
    """Representatives that cover the species, and the coverage they achieve.

    Returns (chosen, achieved_identity): every genome of the species is at
    least `achieved_identity` similar to one of the chosen representatives.
    That number is the guarantee the panel can actually make, and it is
    measured rather than assumed.

    Taking the first max_n clusters would keep whichever genomes happen to be
    longest and could drop a whole lineage: for Lassa that could mean losing
    clade III entirely. This takes the longest genome first, then repeatedly
    the representative furthest from everything already chosen, which is the
    standard greedy k-centre selection and bounds exactly that worst case.

    Two guards on top of it:

    `coverage_identity` stops early once every genome sits within that identity
    of some representative. A species whose genomes are all alike needs one;
    SARS-CoV-2 deposits differ by well under a percent, and taking twelve of
    them would be twelve copies of the same thing.

    `min_members` keeps a cluster of one from being picked while clusters with
    corroboration are still available. Furthest-first reaches for whatever is
    most unlike everything else, which is precisely what a misannotated,
    chimeric or heavily-N'd deposit looks like. Requiring a second genome to
    agree costs little: a real lineage rarely has exactly one deposit, and
    singletons are still taken when nothing better remains.
    """
    if not clusters:
        return [], 1.0

    def dist(a, b):
        return 1.0 - mash_identity(jaccard(a["sketch"], b["sketch"], size), k)

    pool = list(clusters)
    chosen = [pool.pop(0)]                  # already sorted longest first
    while pool and len(chosen) < max_n:
        scored = [(min(dist(c, x) for x in chosen), c) for c in pool]
        if coverage_identity is not None and max(d for d, _ in scored) <= 1.0 - coverage_identity:
            break                            # everything is already covered
        corroborated = [t for t in scored if len(t[1]["members"]) >= min_members]
        _, best = max(corroborated or scored, key=lambda t: t[0])
        chosen.append(best)
        pool.remove(best)

    achieved = max((min(dist(c, x) for x in chosen) for c in pool), default=0.0)
    return chosen, 1.0 - achieved


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
