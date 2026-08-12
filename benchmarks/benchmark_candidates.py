"""Reproducible comparison of indexed and all-pairs candidate generation."""

from __future__ import annotations

import argparse
import json
import random
import time
from collections.abc import Callable

from umi_collapse import indexed_adjacency_edges, naive_adjacency_edges

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
    parser.add_argument("--umis", type=int, default=1_000)
    parser.add_argument("--length", type=int, default=12)
    args = parser.parse_args()
    print(json.dumps(benchmark(args.seed, args.umis, args.length), indent=2))


if __name__ == "__main__":
    main()
