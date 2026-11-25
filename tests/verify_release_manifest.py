"""Validate the current data archive's exact inventory and TSV structure."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path, PurePosixPath
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/release"))
from layout import canonical_path


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def open_data(path):
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def verify(root, report=None):
    root = root.resolve()
    checksums = {}
    for line in (root/"SHA256SUMS.txt").read_text().splitlines():
        digest, rel = line.split(maxsplit=1)
        p = PurePosixPath(rel)
        assert not p.is_absolute() and ".." not in p.parts and rel not in checksums
        assert re.fullmatch("[0-9a-f]{64}",digest)
        checksums[rel] = digest
    files = set()
    for p in root.rglob("*"):
        assert not p.is_symlink(), p
        if p.is_file():
            files.add(p.relative_to(root).as_posix())
    assert files == {*checksums,"SHA256SUMS.txt"}, (files-set(checksums),set(checksums)-files)
    checks = []
    for rel, expected in sorted(checksums.items()):
        path = root/rel
        assert sha(path) == expected, rel
        if not rel.endswith((".tsv",".tsv.gz")):
            continue
        with open_data(path) as f:
            r = csv.reader(f,delimiter="\t")
            header = next(r)
            assert all(header) and len(set(header)) == len(header), rel
            n = 0
            for row in r:
                if not row:
                    continue
                assert len(row) == len(header), (rel, r.line_num, len(header), len(row))
                n += 1
        checks.append((rel,len(header),n,"PASS"))
    paths = (
        "data/processed/whole_blood_pair_program_differences.tsv",
        "06_figure_table_source_data/supplementary/Supplementary_Figure_S4_source_data.tsv")
    for rel in paths:
        with canonical_path(root, rel).open() as stream:
            header = next(csv.reader(stream,delimiter="\t"))
        assert {"anonymous_pair_index","ASD_minus_sibling"} <= set(header)
        assert not {"ASD_score","sibling_score"}.intersection(header)
    if report:
        with report.open("w",newline="") as out:
            w=csv.writer(out,delimiter="\t",lineterminator="\n")
            w.writerow(["path","columns","records","status"])
            w.writerows(checks)
    return len(checksums), len(checks)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("data_root",type=Path)
    cli.add_argument("--report",type=Path)
    args=cli.parse_args()
    print("Integrity verified for %d files and %d TSV schemas." % verify(args.data_root,args.report))
