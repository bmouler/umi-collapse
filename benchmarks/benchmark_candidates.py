"""Reproducible candidate and end-to-end public-collapse benchmarks."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import time
from collections.abc import Callable

from umi_collapse import collapse, indexed_adjacency_edges, naive_adjacency_edges

DNA = "ACGT"


def dataset(seed: int, target_umis: int, length: int) -> set[str]:
    """Create seeded parent UMIs and plausible one-base sequencing errors."""
    randomizer = random.Random(seed)
    umis: set[str] = set()
    while len(umis) < target_umis:
        parent = "".join(randomizer.choices(DNA, k=length))
        umis.add(parent)
        position = randomizer.randrange(length)
        replacement = randomizer.choice(DNA.replace(parent[position], ""))
        umis.add(f"{parent[:position]}{replacement}{parent[position + 1 :]}")
    return umis


def count_dataset(
    seed: int, target_umis: int, length: int, errors_per_parent: int = 3
) -> dict[str, int]:
    """Create a seeded, count-skewed table of parents and substitution errors."""
    randomizer = random.Random(seed)
    counts: dict[str, int] = {}
    rank = 0
    family_size = errors_per_parent + 1
    while len(counts) < target_umis:
        parent = "".join(randomizer.choices(DNA, k=length))
        if parent in counts:
            continue
        family = [parent]
        for error_index in range(errors_per_parent):
            position = (rank * 5 + error_index * 7) % length
            original = parent[position]
            replacement = DNA[(DNA.index(original) + error_index + 1) % len(DNA)]
            error = f"{parent[:position]}{replacement}{parent[position + 1 :]}"
            if error in counts or error in family:
                break
            family.append(error)
        if len(family) != family_size:
            continue

        parent_count = max(16, 4_096 // (1 + rank // 25))
        counts[parent] = parent_count
        for error_index, error in enumerate(family[1:]):
            if len(counts) == target_umis:
                break
            counts[error] = max(1, (parent_count + 1) // (2 ** (error_index + 2)))
        rank += 1
    return counts


def _collapse_workload(counts: dict[str, int]) -> tuple[object, object]:
    return (
        collapse(counts, mode="adjacency"),
        collapse(counts, mode="directional"),
    )


def _cluster_checksum(results: tuple[object, object]) -> str:
    serializable = [
        [
            [cluster.representative, cluster.total, list(cluster.members)]
            for cluster in clusters
        ]
        for clusters in results
    ]
    payload = json.dumps(serializable, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def benchmark_collapse(
    seed: int = 2026,
    target_umis: int = 4_800,
    length: int = 12,
    samples: int = 11,
    warmups: int = 3,
    expected_checksum: str | None = None,
) -> dict[str, object]:
    """Time both public collapse modes over one fixed count table."""
    if samples < 1 or warmups < 0:
        raise ValueError("samples must be positive and warmups must be nonnegative")
    counts = count_dataset(seed, target_umis, length)
    result = _collapse_workload(counts)
    checksum = _cluster_checksum(result)
    if expected_checksum is not None and checksum != expected_checksum:
        raise AssertionError(
            f"collapse checksum changed: expected {expected_checksum}, got {checksum}"
        )

    for _ in range(warmups):
        if _cluster_checksum(_collapse_workload(counts)) != checksum:
            raise RuntimeError("collapse result changed during warmup")
    timings: list[float] = []
    for _ in range(samples):
        started = time.perf_counter()
        measured = _collapse_workload(counts)
        timings.append(time.perf_counter() - started)
        if _cluster_checksum(measured) != checksum:
            raise RuntimeError("collapse result changed between timed samples")
    ordered_timings = sorted(timings)
    adjacency, directional = result
    return {
        "seed": seed,
        "umis": len(counts),
        "length": length,
        "total_count": sum(counts.values()),
        "modes": ["adjacency", "directional"],
        "adjacency_clusters": len(adjacency),
        "directional_clusters": len(directional),
        "members_per_mode": len(counts),
        "warmups": warmups,
        "samples": samples,
        "median_seconds": ordered_timings[len(ordered_timings) // 2],
        "min_seconds": ordered_timings[0],
        "max_seconds": ordered_timings[-1],
        "sample_seconds": timings,
        "checksum": checksum,
        "exact_checksum_match": (
            None if expected_checksum is None else checksum == expected_checksum
        ),
    }


def _median_runtime(function: Callable[[set[str]], object], data: set[str]) -> float:
    samples: list[float] = []
    for _ in range(5):
        started = time.perf_counter()
        function(data)
        samples.append(time.perf_counter() - started)
    return sorted(samples)[2]


def benchmark(
    seed: int = 2026, target_umis: int = 1_000, length: int = 12
) -> dict[str, object]:
    """Benchmark both generators and prove their edge sets are identical."""
    data = dataset(seed, target_umis, length)
    naive_edges = naive_adjacency_edges(data)
    indexed_edges = indexed_adjacency_edges(data)
    if naive_edges != indexed_edges:
        raise AssertionError("candidate generators produced different adjacency graphs")
    naive_seconds = _median_runtime(naive_adjacency_edges, data)
    indexed_seconds = _median_runtime(indexed_adjacency_edges, data)
    return {
        "seed": seed,
        "umis": len(data),
        "length": length,
        "identical_edges": True,
        "edges": len(naive_edges),
        "naive_seconds": naive_seconds,
        "indexed_seconds": indexed_seconds,
        "speedup": naive_seconds / indexed_seconds,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--umis", type=int, default=4_800)
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--samples", type=int, default=11)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--expected-checksum")
    args = parser.parse_args()
    print(
        json.dumps(
            benchmark_collapse(
                seed=args.seed,
                target_umis=args.umis,
                length=args.length,
                samples=args.samples,
                warmups=args.warmups,
                expected_checksum=args.expected_checksum,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
