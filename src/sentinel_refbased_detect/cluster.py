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

K = 21
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
    for c in clusters:
        c.pop("sketch", None)
    return clusters


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
