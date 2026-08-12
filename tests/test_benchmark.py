from __future__ import annotations

from benchmarks.benchmark_candidates import benchmark, dataset


def test_seeded_dataset_is_reproducible() -> None:
    assert dataset(7, 5, 4) == dataset(7, 5, 4)
    assert len(dataset(7, 5, 4)) >= 5


def test_benchmark_confirms_identical_graphs() -> None:
    result = benchmark(seed=7, target_umis=20, length=6)
    assert result["identical_edges"] is True
    assert result["umis"] >= 20
    assert result["edges"] > 0
    assert result["naive_seconds"] >= 0
    assert result["indexed_seconds"] >= 0
    assert result["speedup"] >= 0
