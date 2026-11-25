"""Verify the immune-context inputs against their source checksums."""

import csv
from hashlib import file_digest
from pathlib import Path

CONFIG = Path(__file__).with_name("config")


def records(name):
    with (CONFIG / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def verify_file(path, info):
    expected = f"{info['bytes']} bytes and SHA-256 {info['sha256']}"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {info['path']}; expected {expected}")
    with path.open("rb") as stream:
        actual = file_digest(stream, "sha256").hexdigest()
    if path.stat().st_size != int(info["bytes"]) or actual != info["sha256"]:
        raise ValueError(f"Source fingerprint differs for {info['path']}; expected {expected}, got {actual}")
    return path


def verify_baseline_file(data, name):
    original = next(row for row in records("frozen_baseline_inputs.tsv") if row["historical_alias"] == name)
    canonical = data / original["path"]
    return verify_file(canonical if canonical.is_file() else data / name, original)
