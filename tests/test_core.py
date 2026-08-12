from __future__ import annotations

import pytest

from umi_collapse import (
    Cluster,
    collapse,
    hamming_distance,
    indexed_adjacency_edges,
    naive_adjacency_edges,
    one_edit_neighbors,
)


def test_hamming_distance() -> None:
    assert hamming_distance("AAAA", "AATA") == 1
    assert hamming_distance("ACGT", "TGCA") == 4
    with pytest.raises(ValueError, match=r"equal-length"):
        hamming_distance("A", "AA")
    with pytest.raises(
        ValueError, match=r"^Hamming distance requires equal-length strings$"
    ):
        hamming_distance("A", "AA")
    # The explicit length check makes zip truncation impossible by construction.


def test_one_edit_neighbors_are_complete_and_stable() -> None:
    neighbors = list(one_edit_neighbors("AC"))
    assert neighbors == ["CC", "GC", "TC", "AA", "AG", "AT"]


def test_candidate_generators_match() -> None:
    umis = {"AAAA", "AAAT", "AATT", "TTTT"}
    expected = {("AAAA", "AAAT"), ("AAAT", "AATT")}
    assert naive_adjacency_edges(umis) == expected
    assert indexed_adjacency_edges(umis) == expected


def test_candidate_generators_reject_mixed_lengths_consistently() -> None:
    for generator in (naive_adjacency_edges, indexed_adjacency_edges):
        with pytest.raises(ValueError, match=r"equal-length"):
            generator({"AAA", "AAAT"})


def test_adjacency_connects_through_lexically_later_center() -> None:
    counts = {"ACC": 2, "CAC": 3, "CCC": 5, "TTT": 7}
    assert collapse(counts, mode="adjacency") == [
        Cluster("CCC", 10, ("CCC", "CAC", "ACC")),
        Cluster("TTT", 7, ("TTT",)),
    ]
    assert indexed_adjacency_edges({"AA"}) == set()


def test_adjacency_transitively_connects_and_orders() -> None:
    counts = {"AAAA": 10, "AAAT": 3, "AATT": 4, "TTTT": 8}
    expected = [
        Cluster("AAAA", 17, ("AAAA", "AATT", "AAAT")),
        Cluster("TTTT", 8, ("TTTT",)),
    ]
    assert collapse(counts, mode="adjacency") == expected
    assert collapse(counts, mode="adjacency", indexed=False) == expected


def test_directional_threshold_chains_but_not_equal_high_counts() -> None:
    counts = {"AAAA": 10, "AAAT": 5, "AATT": 3, "TTTT": 10, "TTTA": 10}
    assert collapse(counts) == [
        Cluster("AAAA", 18, ("AAAA", "AAAT", "AATT")),
        Cluster("TTTA", 10, ("TTTA",)),
        Cluster("TTTT", 10, ("TTTT",)),
    ]


def test_directional_uses_count_order_bidirectional_edges_and_exact_threshold() -> None:
    assert collapse({"AA": 1, "AC": 2}) == [Cluster("AC", 3, ("AC", "AA"))]
    assert collapse({"AA": 3, "AC": 4}) == [
        Cluster("AC", 4, ("AC",)),
        Cluster("AA", 3, ("AA",)),
    ]
    assert collapse({"AA": 2, "AC": 2}) == [
        Cluster("AA", 2, ("AA",)),
        Cluster("AC", 2, ("AC",)),
    ]
    assert collapse({"AA": 5, "AC": 3, "AG": 2}) == [
        Cluster("AA", 10, ("AA", "AC", "AG"))
    ]


def test_deterministic_ties_and_case_normalization() -> None:
    assert collapse({"aaaa": 2, "AAAT": 2}, mode="adjacency") == [
        Cluster("AAAA", 4, ("AAAA", "AAAT"))
    ]
    with pytest.raises(ValueError, match=r"duplicate UMI after normalization"):
        collapse({"aaaa": 1, "AAAA": 2})


def test_empty_input_and_bad_mode() -> None:
    assert collapse({}, mode="adjacency") == []
    assert collapse({}, mode="directional") == []
    with pytest.raises(ValueError, match=r"unknown collapse mode"):
        collapse({}, mode="other")


@pytest.mark.parametrize(
    ("counts", "message"),
    [
        ({"": 1}, "nonempty"),
        ({1: 1}, "nonempty"),
        ({"AAN": 1}, "invalid DNA"),
        ({"AAA": 1, "AAAA": 2}, "same length"),
        ({"AAA": 0}, "positive integer"),
        ({"AAA": -1}, "positive integer"),
        ({"AAA": True}, "positive integer"),
        ({"AAA": 1.5}, "positive integer"),
    ],
)
def test_invalid_counts(counts: dict[object, object], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        collapse(counts)  # type: ignore[arg-type]


def test_validation_error_messages_are_stable() -> None:
    cases = [
        ({"": 1}, "UMIs must be nonempty DNA strings"),
        ({"AAA": 1, "AAAA": 2}, "all UMIs must have the same length"),
    ]
    for counts, message in cases:
        with pytest.raises(ValueError) as caught:
            collapse(counts)
        assert str(caught.value) == message
    with pytest.raises(
        ValueError, match=r"^adjacency edges require equal-length UMIs$"
    ):
        indexed_adjacency_edges({"A", "AA"})
