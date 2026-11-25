"""Verify that the release and code carry the historical inputs for directional and phenotype comparisons."""

import csv
import sys
from hashlib import file_digest
from pathlib import Path

RELEASES = {
    "ASD_TREG_CRISPR_processed_data_v1.1.0",
    "ASD_TREG_CRISPR_processed_data_v1.2.0",
}


def verify_release(data: Path) -> None:
    if (data / "VERSION.txt").read_text(encoding="utf-8").strip() not in RELEASES:
        raise ValueError("Expected an extracted ASD_TREG_CRISPR processed-data release v1.1.0 or v1.2.0")


def verify_inputs(data: Path, repo: Path) -> None:
    verify_release(data)
    roots = {"data": data, "repo": repo}
    with Path(__file__).with_name("frozen_inputs.tsv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = roots[row["root"]] / row["path"]
            if row["historical_alias"] and not path.exists():
                path = data / row["historical_alias"]
            with path.open("rb") as source:
                fingerprint = file_digest(source, "sha256").hexdigest()
            if path.stat().st_size != int(row["bytes"]) or fingerprint != row["sha256"]:
                raise ValueError(f"Input differs from the frozen v1.1.0 snapshot: {row['path']}")


if __name__ == "__main__":
    verify_inputs(Path(sys.argv[1]), Path(sys.argv[2]))
