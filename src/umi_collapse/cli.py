"""Command-line interface for UMI collapse."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from .core import Cluster, collapse


class CliExit(Exception):
    """A requested command-line exit."""

    def __init__(self, status: int, message: str | None = None) -> None:
        super().__init__(message)
        self.status = status
        self.message = message


class UsageError(Exception):
    """A command-line usage or input error."""


class Parser(argparse.ArgumentParser):
    """Argument parser that lets ``main`` honor its integer return contract."""

    def error(self, message: str) -> None:
        raise UsageError(message)

    def exit(self, status: int = 0, message: str | None = None) -> None:
        raise CliExit(status, message)


def _parser() -> Parser:
    parser = Parser(prog="umi-collapse", description=__doc__)
    parser.add_argument("input", type=Path, help="TSV with umi and count columns")
    parser.add_argument(
        "-o", "--output", type=Path, help="output path (default: stdout)"
    )
    parser.add_argument(
        "--mode",
        choices=("adjacency", "directional"),
        default="directional",
    )
    parser.add_argument("--json", action="store_true", help="write JSON instead of TSV")
    parser.add_argument(
        "--naive",
        action="store_true",
        help="use all-pairs candidate generation (adjacency mode only)",
    )
    return parser


def read_counts(stream: TextIO) -> dict[str, int]:
    """Read and validate a two-column ``umi``/``count`` TSV."""
    reader = csv.reader(stream, delimiter="\t", strict=True)
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValueError("input TSV is empty") from exc
    except csv.Error as exc:
        raise ValueError(f"bad TSV: {exc}") from exc
    if header != ["umi", "count"]:
        raise ValueError("TSV header must be exactly: umi<TAB>count")

    counts: dict[str, int] = {}
    try:
        for line_number, row in enumerate(reader, start=2):
            if len(row) != 2 or not row[0] or not row[1]:
                raise ValueError(f"bad TSV row {line_number}: expected two fields")
            umi = row[0].upper()
            if umi in counts:
                raise ValueError(f"duplicate UMI on row {line_number}: {umi}")
            try:
                count = int(row[1])
            except ValueError as exc:
                raise ValueError(f"bad count on row {line_number}: {row[1]!r}") from exc
            counts[umi] = count
    except csv.Error as exc:
        raise ValueError(f"bad TSV: {exc}") from exc
    if not counts:
        raise ValueError("input TSV has no data rows")
    return counts


def _records(clusters: Sequence[Cluster]) -> list[dict[str, object]]:
    return [
        {
            "cluster": number,
            "representative": cluster.representative,
            "total": cluster.total,
            "members": list(cluster.members),
        }
        for number, cluster in enumerate(clusters, start=1)
    ]


def write_clusters(stream: TextIO, clusters: Sequence[Cluster], as_json: bool) -> None:
    """Write cluster records as TSV or JSON."""
    records = _records(clusters)
    if as_json:
        json.dump(records, stream, indent=2)
        stream.write("\n")
        return
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(("cluster", "representative", "total", "members"))
    for record in records:
        writer.writerow(
            (
                record["cluster"],
                record["representative"],
                record["total"],
                ",".join(record["members"]),
            )
        )


def main(argv: list[str] | None = None) -> int:
    """Run the command and return a process status."""
    parser = _parser()
    try:
        args = parser.parse_args(argv)
        with args.input.open(encoding="utf-8", newline="") as input_stream:
            counts = read_counts(input_stream)
        clusters = collapse(counts, mode=args.mode, indexed=not args.naive)
        if args.output is None:
            write_clusters(sys.stdout, clusters, args.json)
        else:
            with args.output.open("w", encoding="utf-8", newline="") as output_stream:
                write_clusters(output_stream, clusters, args.json)
    except CliExit as exc:
        if exc.message:
            print(exc.message, end="", file=sys.stderr)
        return exc.status
    except (UsageError, ValueError, OSError) as exc:
        print(f"umi-collapse: error: {exc}", file=sys.stderr)
        return 2
    return 0
