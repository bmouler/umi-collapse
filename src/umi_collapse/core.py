"""Deterministic UMI clustering algorithms."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Literal, TypeAlias

DNA = "ACGT"

CollapseMode: TypeAlias = Literal["adjacency", "directional"]
Edge: TypeAlias = tuple[str, str]
EdgeSet: TypeAlias = set[Edge]


@dataclass(frozen=True)
class Cluster:
    """A collapsed UMI cluster."""

    representative: str
    total: int
    members: tuple[str, ...]


def hamming_distance(left: str, right: str) -> int:
    """Return the number of differing positions in equal-length strings."""
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(a != b for a, b in zip(left, right, strict=False))


def _validate_counts(counts: Mapping[str, int]) -> dict[str, int]:
    if not counts:
        return {}
    normalized: dict[str, int] = {}
    length: int | None = None
    for umi, count in counts.items():
        if not isinstance(umi, str) or not umi:
            raise ValueError("UMIs must be nonempty DNA strings")
        sequence = umi.upper()
        if any(base not in DNA for base in sequence):
            raise ValueError(f"invalid DNA UMI: {umi!r}")
        if length is None:
            length = len(sequence)
        elif len(sequence) != length:
            raise ValueError("all UMIs must have the same length")
        if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise ValueError(f"count for {umi!r} must be a positive integer")
        if sequence in normalized:
            raise ValueError(f"duplicate UMI after normalization: {sequence}")
        normalized[sequence] = count
    return normalized


def one_edit_neighbors(umi: str) -> Iterator[str]:
    """Yield every DNA string at Hamming distance one in stable order."""
    for index, original in enumerate(umi):
        for base in DNA:
            if base != original:
                yield f"{umi[:index]}{base}{umi[index + 1 :]}"


def _ordered_equal_length_umis(umis: Iterable[str]) -> list[str]:
    ordered = sorted(umis)
    if ordered and any(len(umi) != len(ordered[0]) for umi in ordered[1:]):
        raise ValueError("adjacency edges require equal-length UMIs")
    return ordered


def naive_adjacency_edges(umis: Iterable[str]) -> EdgeSet:
    """Generate radius-one edges by all-pairs comparison."""
    ordered = _ordered_equal_length_umis(umis)
    return {
        (left, right)
        for index, left in enumerate(ordered)
        for right in ordered[index + 1 :]
        if hamming_distance(left, right) <= 1
    }


def indexed_adjacency_edges(umis: Iterable[str]) -> EdgeSet:
    """Generate radius-one edges by enumerating substitutions."""
    present = set(_ordered_equal_length_umis(umis))
    edges: EdgeSet = set()
    for umi in sorted(present):
        for neighbor in one_edit_neighbors(umi):
            if neighbor in present and umi < neighbor:
                edges.add((umi, neighbor))
    return edges


def _components(nodes: Iterable[str], edges: Iterable[Edge]) -> list[set[str]]:
    graph: dict[str, set[str]] = {node: set() for node in nodes}
    for left, right in edges:
        graph[left].add(right)
        graph[right].add(left)
    components: list[set[str]] = []
    unseen = set(graph)
    while unseen:
        start = min(unseen)
        stack = [start]
        component: set[str] = set()
        unseen.remove(start)
        while stack:
            current = stack.pop()
            component.add(current)
            new = unseen.intersection(graph[current])
            unseen.difference_update(new)
            stack.extend(new)
        components.append(component)
    return components


def _make_cluster(members: set[str], counts: Mapping[str, int]) -> Cluster:
    ordered = tuple(sorted(members, key=lambda umi: (-counts[umi], umi)))
    return Cluster(ordered[0], sum(counts[umi] for umi in ordered), ordered)


def _adjacency(counts: Mapping[str, int], indexed: bool) -> list[Cluster]:
    edge_function = indexed_adjacency_edges if indexed else naive_adjacency_edges
    return [
        _make_cluster(group, counts)
        for group in _components(counts, edge_function(counts))
    ]


def _directional(counts: Mapping[str, int]) -> list[Cluster]:
    neighbors: dict[str, set[str]] = {umi: set() for umi in counts}
    for high, low in indexed_adjacency_edges(counts):
        neighbors[high].add(low)
        neighbors[low].add(high)

    remaining = set(counts)
    clusters: list[Cluster] = []
    for root in sorted(counts, key=lambda umi: (-counts[umi], umi)):
        if root not in remaining:
            continue
        remaining.remove(root)
        members = {root}
        frontier = [root]
        while frontier:
            current = frontier.pop()
            absorbable = [
                candidate
                for candidate in remaining.intersection(neighbors[current])
                if counts[current] >= 2 * counts[candidate] - 1
            ]
            remaining.difference_update(absorbable)
            members.update(absorbable)
            frontier.extend(absorbable)
        clusters.append(_make_cluster(members, counts))
    return clusters


def collapse(
    counts: Mapping[str, int],
    *,
    mode: CollapseMode = "directional",
    indexed: bool = True,
) -> list[Cluster]:
    """Collapse UMI counts using deterministic radius-one clustering.

    ``mode`` is ``"adjacency"`` or ``"directional"``. The ``indexed`` flag
    selects indexed or all-pairs candidate generation for adjacency mode.
    """
    validated = _validate_counts(counts)
    if mode == "adjacency":
        clusters = _adjacency(validated, indexed)
    elif mode == "directional":
        clusters = _directional(validated)
    else:
        raise ValueError(f"unknown collapse mode: {mode!r}")
    return sorted(
        clusters, key=lambda cluster: (-cluster.total, cluster.representative)
    )
