from __future__ import annotations

import io
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

import umi_collapse.cli as cli
from umi_collapse.cli import main, read_counts, write_clusters
from umi_collapse.core import Cluster


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ('"umi\tcount\n', "bad TSV"),
        ("UMI\tcount\nAAAA\t1\n", "header"),
        ("umi\tcount\n", "no data"),
        ("umi\tcount\nAAAA\n", "expected two"),
        ("umi\tcount\n\t1\n", "expected two"),
        ("umi\tcount\nAAAA\tx\n", "bad count"),
        ("umi\tcount\nAAAA\t1\naaaa\t2\n", "duplicate"),
        ('umi\tcount\n"AAAA\t1\n', "bad TSV"),
    ],
)
def test_read_counts_rejects_bad_tsv(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        read_counts(io.StringIO(text))


def test_read_counts_rejects_empty_stream() -> None:
    with pytest.raises(ValueError, match="empty"):
        read_counts(io.StringIO(""))


def test_read_counts_normalizes_case() -> None:
    assert read_counts(io.StringIO("umi\tcount\naaaa\t2\n")) == {"AAAA": 2}


def test_writers() -> None:
    clusters = [Cluster("AAAA", 3, ("AAAA", "AAAT"))]
    stream = io.StringIO()
    write_clusters(stream, clusters, False)
    assert stream.getvalue() == (
        "cluster\trepresentative\ttotal\tmembers\n1\tAAAA\t3\tAAAA,AAAT\n"
    )
    stream = io.StringIO()
    write_clusters(stream, clusters, True)
    assert json.loads(stream.getvalue()) == [
        {
            "cluster": 1,
            "representative": "AAAA",
            "total": 3,
            "members": ["AAAA", "AAAT"],
        }
    ]


def test_main_stdout_json_help_and_errors(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    source = tmp_path / "counts.tsv"
    source.write_text("umi\tcount\nAAAA\t4\nAAAT\t2\n", encoding="utf-8")
    assert main([str(source), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["total"] == 6
    assert main(["--help"]) == 0
    assert "usage: umi-collapse" in capsys.readouterr().out

    assert main([str(source), "--mode", "bad"]) == 2
    assert "invalid choice" in capsys.readouterr().err
    assert main([str(tmp_path / "missing.tsv")]) == 2
    assert "error:" in capsys.readouterr().err


def test_main_handles_parser_exit_message(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    class ExitingParser:
        def parse_args(self, argv: list[str] | None) -> None:
            raise cli.CliExit(3, "custom exit\n")

    monkeypatch.setattr(cli, "_parser", ExitingParser)
    assert cli.main([]) == 3
    assert capsys.readouterr().err == "custom exit\n"


def test_main_output_file_and_naive(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv"
    output = tmp_path / "clusters.tsv"
    source.write_text("umi\tcount\nAAAA\t4\nAAAT\t2\n", encoding="utf-8")
    assert main([str(source), "--mode", "adjacency", "--naive", "-o", str(output)]) == 0
    assert output.read_text(encoding="utf-8").splitlines()[1] == "1\tAAAA\t6\tAAAA,AAAT"


def test_main_reports_output_error(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    source = tmp_path / "counts.tsv"
    source.write_text("umi\tcount\nAAAA\t1\n", encoding="utf-8")
    assert main([str(source), "-o", str(tmp_path / "missing" / "out.tsv")]) == 2
    assert "error:" in capsys.readouterr().err


def test_module_entrypoint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli, "main", lambda: 7)
    with pytest.raises(SystemExit) as caught:
        runpy.run_module("umi_collapse.__main__", run_name="__main__")
    assert caught.value.code == 7


def test_module_e2e_from_clean_directory(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv"
    output = tmp_path / "clusters.json"
    source.write_text(
        "umi\tcount\nAAAA\t10\nAAAT\t5\nAATT\t3\nCCCC\t7\n",
        encoding="utf-8",
    )
    command = Path(sys.executable).with_name("umi-collapse")
    completed = subprocess.run(
        [
            str(command),
            str(source),
            "--json",
            "-o",
            str(output),
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(output.read_text(encoding="utf-8")) == [
        {
            "cluster": 1,
            "representative": "AAAA",
            "total": 18,
            "members": ["AAAA", "AAAT", "AATT"],
        },
        {
            "cluster": 2,
            "representative": "CCCC",
            "total": 7,
            "members": ["CCCC"],
        },
    ]
