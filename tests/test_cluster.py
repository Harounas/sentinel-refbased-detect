import random

from sentinel_refbased_detect.cluster import (cluster_sequences, jaccard,
                                              select_diverse,
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


def _clusters(rng, specs):
    """specs: list of (name, base_seq, n_members)."""
    from sentinel_refbased_detect.cluster import sketch as _sk
    out = []
    for name, seq, n in specs:
        out.append({"representative": name, "length": len(seq), "sketch": _sk(seq),
                    "members": [name] * n})
    return out


def test_selection_stops_once_the_species_is_covered():
    """A near-clonal species keeps one representative, not the whole cap."""
    rng = random.Random(21)
    base = _genome(rng, 12000)
    clusters = _clusters(rng, [("a", base, 5),
                               ("b", _mutate(rng, base, 0.002), 4),
                               ("c", _mutate(rng, base, 0.003), 3)])
    kept, covered = select_diverse(clusters, 12, coverage_identity=0.90)
    assert len(kept) == 1 and covered > 0.95


def test_diverse_species_keeps_several_and_reports_coverage():
    rng = random.Random(22)
    a = _genome(rng, 12000)
    clusters = _clusters(rng, [("clade_a", a, 5),
                               ("clade_b", _mutate(rng, a, 0.20), 4),
                               ("clade_c", _mutate(rng, a, 0.25), 3)])
    kept, covered = select_diverse(clusters, 12, coverage_identity=0.90)
    assert len(kept) == 3 and covered > 0.90


def test_cap_reports_the_coverage_it_could_not_reach():
    rng = random.Random(23)
    a = _genome(rng, 12000)
    clusters = _clusters(rng, [("c%d" % i, _mutate(rng, a, 0.05 * i), 3) for i in range(5)])
    kept, covered = select_diverse(clusters, 2, coverage_identity=0.99)
    assert len(kept) == 2 and covered < 0.99      # honest about falling short


def test_singleton_outlier_does_not_displace_a_real_lineage():
    """Furthest-first reaches for the most unlike thing, which is what a
    misannotated deposit looks like. The two-member guard holds it back."""
    rng = random.Random(24)
    a = _genome(rng, 12000)
    clusters = _clusters(rng, [("ref", a, 6),
                               ("real_lineage", _mutate(rng, a, 0.12), 4),
                               ("junk_singleton", _genome(rng, 12000), 1)])
    kept, _ = select_diverse(clusters, 2, min_members=2, coverage_identity=0.50)
    assert [c["representative"] for c in kept] == ["ref", "real_lineage"]
    # With no corroborated cluster left, a singleton is still allowed.
    kept3, _ = select_diverse(clusters, 3, min_members=2, coverage_identity=0.50)
    assert "junk_singleton" in [c["representative"] for c in kept3]
