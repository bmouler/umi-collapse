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


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("", "input TSV is empty"),
        ("UMI\tcount\nAAAA\t1\n", "TSV header must be exactly: umi<TAB>count"),
        ("umi\tcount\n", "input TSV has no data rows"),
        ("umi\tcount\nAAAA\n", "bad TSV row 2: expected two fields"),
        ("umi\tcount\nAAAA\t1\nAAAT\n", "bad TSV row 3: expected two fields"),
    ],
)
def test_read_counts_errors_include_exact_context(text: str, message: str) -> None:
    with pytest.raises(ValueError) as caught:
        read_counts(io.StringIO(text))
    assert str(caught.value) == message


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
    assert stream.getvalue() == (
        "[\n"
        "  {\n"
        '    "cluster": 1,\n'
        '    "representative": "AAAA",\n'
        '    "total": 3,\n'
        '    "members": [\n'
        '      "AAAA",\n'
        '      "AAAT"\n'
        "    ]\n"
        "  }\n"
        "]\n"
    )


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


def test_help_text_documents_cli_contract(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--help"]) == 0
    assert capsys.readouterr().out == (
        "usage: umi-collapse [-h] [-o OUTPUT] [--mode {adjacency,directional}] "
        "[--json]\n"
        "                    [--naive]\n"
        "                    input\n"
        "\n"
        "Command-line interface for UMI collapse.\n"
        "\n"
        "positional arguments:\n"
        "  input                 TSV with umi and count columns\n"
        "\n"
        "options:\n"
        "  -h, --help            show this help message and exit\n"
        "  -o OUTPUT, --output OUTPUT\n"
        "                        output path (default: stdout)\n"
        "  --mode {adjacency,directional}\n"
        "  --json                write JSON instead of TSV\n"
        "  --naive               use all-pairs candidate generation "
        "(adjacency mode\n"
        "                        only)\n"
    )


def test_main_handles_parser_exit_message(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    class ExitingParser:
        def parse_args(self, argv: list[str] | None) -> None:
            raise cli.CliExit(3, "custom exit\n")

    monkeypatch.setattr(cli, "_parser", ExitingParser)
    assert cli.main([]) == 3
    assert capsys.readouterr().err == "custom exit\n"


def test_cli_exit_retains_exception_message() -> None:
    error = cli.CliExit(4, "requested exit")
    assert str(error) == "requested exit"
    assert error.status == 4
    assert error.message == "requested exit"


def test_parser_exit_preserves_message() -> None:
    with pytest.raises(cli.CliExit) as caught:
        cli.Parser().exit(5, "stop\n")
    assert caught.value.status == 5
    assert caught.value.message == "stop\n"


def test_main_output_file_and_naive(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv"
    output = tmp_path / "clusters.tsv"
    source.write_text("umi\tcount\nAAAA\t4\nAAAT\t2\n", encoding="utf-8")
    assert main([str(source), "--mode", "adjacency", "--naive", "-o", str(output)]) == 0
    assert output.read_text(encoding="utf-8").splitlines()[1] == "1\tAAAA\t6\tAAAA,AAAT"


def test_main_honors_adjacency_mode_and_json_file(tmp_path: Path) -> None:
    source = tmp_path / "counts.tsv"
    output = tmp_path / "clusters.json"
    source.write_text("umi\tcount\nAAAA\t10\nAAAT\t10\nAATT\t1\n", encoding="utf-8")
    assert main([str(source), "--mode", "adjacency", "--json", "-o", str(output)]) == 0
    assert json.loads(output.read_text(encoding="utf-8")) == [
        {
            "cluster": 1,
            "representative": "AAAA",
            "total": 21,
            "members": ["AAAA", "AAAT", "AATT"],
        }
    ]


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
