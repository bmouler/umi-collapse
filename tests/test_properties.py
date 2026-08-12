from __future__ import annotations

from hypothesis import given
from hypothesis import strategies as st

from umi_collapse import (
    collapse,
    hamming_distance,
    indexed_adjacency_edges,
    naive_adjacency_edges,
)


@st.composite
def umi_lists(draw: st.DrawFn) -> list[str]:
    length = draw(st.integers(min_value=1, max_value=8))
    return draw(
        st.lists(
            st.text(alphabet="ACGT", min_size=length, max_size=length),
            unique=True,
            max_size=40,
        )
    )


@st.composite
def umi_counts(draw: st.DrawFn) -> dict[str, int]:
    umis = draw(umi_lists())
    counts = draw(
        st.lists(
            st.integers(min_value=1, max_value=10_000),
            min_size=len(umis),
            max_size=len(umis),
        )
    )
    return dict(zip(umis, counts, strict=True))


@st.composite
def umi_pairs(draw: st.DrawFn) -> tuple[str, str]:
    length = draw(st.integers(min_value=1, max_value=8))
    umi = st.text(alphabet="ACGT", min_size=length, max_size=length)
    return draw(umi), draw(umi)


@given(umis=umi_lists())
def test_indexed_adjacency_matches_naive(umis: list[str]) -> None:
    assert indexed_adjacency_edges(umis) == naive_adjacency_edges(umis)


@given(counts=umi_counts())
def test_collapse_conserves_counts_and_partitions_members(
    counts: dict[str, int],
) -> None:
    for mode in ("adjacency", "directional"):
        clusters = collapse(counts, mode=mode)
        members = [member for cluster in clusters for member in cluster.members]

        assert sum(cluster.total for cluster in clusters) == sum(counts.values())
        assert len(members) == len(set(members))
        assert set(members) == set(counts)


@given(pair=umi_pairs())
def test_hamming_distance_matches_definition_and_is_symmetric(
    pair: tuple[str, str],
) -> None:
    left, right = pair
    expected = sum(
        left_base != right_base
        for left_base, right_base in zip(left, right, strict=True)
    )

    assert hamming_distance(left, right) == expected
    assert hamming_distance(left, right) == hamming_distance(right, left)
