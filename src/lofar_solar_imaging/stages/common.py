from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Iterable


def run_command(cmd: list[str], *, cwd: Path | None = None) -> None:
    """Run a command and stream output.

    Raises
    ------
    subprocess.CalledProcessError
        If the command returns a non-zero status.
    """
    printable = " ".join(cmd)
    print(f"$ {printable}")
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def ensure_exists(path: Path, *, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{description} does not exist: {path}")


def strip_all_suffixes(path: Path) -> str:
    stem = path.name
    while "." in stem:
        stem = stem.rsplit(".", 1)[0]
    return stem


def resolve_glob(pattern: str) -> list[Path]:
    return sorted(Path().glob(pattern))


def as_dp3_msin(paths: Iterable[Path]) -> str:
    quoted = ",".join(f'"{p}"' for p in paths)
    return f"[{quoted}]"
