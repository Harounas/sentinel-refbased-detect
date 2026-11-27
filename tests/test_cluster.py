import random

from sentinel_refbased_detect.cluster import (cluster_sequences, jaccard,
                                              mash_identity, sketch)


def _genome(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng, seq, rate):
    out = []
    for c in seq:
        if rng.random() < rate:
            out.append(rng.choice([x for x in "ACGT" if x != c]))
        else:
            out.append(c)
    return "".join(out)


def test_sketch_is_strand_agnostic():
    rng = random.Random(0)
    seq = _genome(rng, 5000)
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    assert sketch(seq) == sketch(rc)


def test_identity_tracks_divergence():
    rng = random.Random(1)
    seq = _genome(rng, 20000)
    near = mash_identity(jaccard(sketch(seq), sketch(_mutate(rng, seq, 0.01))))
    far = mash_identity(jaccard(sketch(seq), sketch(_mutate(rng, seq, 0.15))))
    unrelated = mash_identity(jaccard(sketch(seq), sketch(_genome(rng, 20000))))
    assert near > 0.95 > far > unrelated
    assert unrelated < 0.75


def test_clustering_collapses_near_identical_keeps_distinct():
    rng = random.Random(2)
    a = _genome(rng, 15000)
    b = _genome(rng, 15000)                      # an unrelated genome
    records = [("a_ref", a),
               ("a_copy1", _mutate(rng, a, 0.005)),
               ("a_copy2", _mutate(rng, a, 0.01)),
               ("b_ref", b)]
    clusters = cluster_sequences(records, min_identity=0.95)
    assert len(clusters) == 2
    sizes = sorted(len(c["members"]) for c in clusters)
    assert sizes == [1, 3]
    # The longest genome of a group represents it.
    big = next(c for c in clusters if len(c["members"]) == 3)
    assert big["representative"] in {"a_ref", "a_copy1", "a_copy2"}
    assert big["length"] == max(len(x) for _, x in records[:3])


def test_strict_threshold_keeps_strains_apart():
    rng = random.Random(3)
    a = _genome(rng, 15000)
    records = [("a", a), ("a_3pct", _mutate(rng, a, 0.03))]
    assert len(cluster_sequences(records, min_identity=0.99)) == 2
    assert len(cluster_sequences(records, min_identity=0.90)) == 1
