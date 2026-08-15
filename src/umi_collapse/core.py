"""Deterministic UMI clustering algorithms."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Literal, TypeAlias

DNA = "ACGT"
_DNA_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}

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


def _encode_umi(umi: str) -> int:
    encoded = 0
    for base in umi:
        encoded = (encoded << 2) | _DNA_CODE[base]
    return encoded


def _indexed_neighbor_indices(ordered: list[str]) -> list[list[int]]:
    """Build a radius-one adjacency list for validated, ordered DNA UMIs."""
    codes = [_encode_umi(umi) for umi in ordered]
    code_to_index = {encoded: index for index, encoded in enumerate(codes)}
    neighbors: list[list[int]] = [[] for _ in ordered]
    width = 2 * len(ordered[0]) if ordered else 0
    for index, encoded in enumerate(codes):
        for shift in range(0, width, 2):
            bit = 1 << shift
            for difference in (bit, bit << 1, bit | (bit << 1)):
                candidate = code_to_index.get(encoded ^ difference)
                if candidate is not None and candidate > index:
                    neighbors[index].append(candidate)
                    neighbors[candidate].append(index)
    return neighbors


def indexed_adjacency_edges(umis: Iterable[str]) -> EdgeSet:
    """Generate radius-one edges by enumerating substitutions."""
    ordered = _ordered_equal_length_umis(umis)
    if any(base not in DNA for umi in ordered for base in umi):
        present = set(ordered)
        return {
            (umi, neighbor)
            for umi in ordered
            for neighbor in one_edit_neighbors(umi)
            if neighbor in present and umi < neighbor
        }
    neighbors = _indexed_neighbor_indices(ordered)
    return {
        (umi, ordered[candidate])
        for index, umi in enumerate(ordered)
        for candidate in neighbors[index]
        if candidate > index
    }


def _components(nodes: Iterable[str], edges: Iterable[Edge]) -> list[set[str]]:
    graph: dict[str, set[str]] = {node: set() for node in nodes}
    for left, right in edges:
        graph[left].add(right)
        graph[right].add(left)
    components: list[set[str]] = []
    seen: set[str] = set()
    for start in sorted(graph):
        if start in seen:
            continue
        stack = [start]
        component: set[str] = set()
        seen.add(start)
        while stack:
            current = stack.pop()
            component.add(current)
            new = graph[current].difference(seen)
            seen.update(new)
            stack.extend(new)
        components.append(component)
    return components


def _make_cluster(members: Iterable[str], counts: Mapping[str, int]) -> Cluster:
    ordered = tuple(sorted(members, key=lambda umi: (-counts[umi], umi)))
    return Cluster(ordered[0], sum(counts[umi] for umi in ordered), ordered)


def _adjacency(counts: Mapping[str, int], indexed: bool) -> list[Cluster]:
    if not indexed:
        return [
            _make_cluster(group, counts)
            for group in _components(counts, naive_adjacency_edges(counts))
        ]

    ordered = sorted(counts)
    neighbors = _indexed_neighbor_indices(ordered)
    seen = bytearray(len(ordered))
    clusters: list[Cluster] = []
    for start in range(len(ordered)):
        if seen[start]:
            continue
        seen[start] = 1
        members: list[str] = []
        stack = [start]
        while stack:
            current = stack.pop()
            members.append(ordered[current])
            for candidate in neighbors[current]:
                if not seen[candidate]:
                    seen[candidate] = 1
                    stack.append(candidate)
        clusters.append(_make_cluster(members, counts))
    return clusters


def _directional(counts: Mapping[str, int]) -> list[Cluster]:
    ordered = sorted(counts)
    neighbors = _indexed_neighbor_indices(ordered)
    remaining = bytearray(b"\x01") * len(ordered)
    clusters: list[Cluster] = []
    roots = sorted(range(len(ordered)), key=lambda index: -counts[ordered[index]])
    for root in roots:
        if not remaining[root]:
            continue
        remaining[root] = 0
        members: list[str] = []
        frontier = [root]
        while frontier:
            current = frontier.pop()
            current_umi = ordered[current]
            members.append(current_umi)
            for candidate in neighbors[current]:
                candidate_umi = ordered[candidate]
                if (
                    remaining[candidate]
                    and counts[current_umi] >= 2 * counts[candidate_umi] - 1
                ):
                    remaining[candidate] = 0
                    frontier.append(candidate)
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
