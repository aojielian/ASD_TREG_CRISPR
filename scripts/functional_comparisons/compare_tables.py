#!/usr/bin/env python3
"""Compare released and replayed TSV values without depending on gzip metadata."""

import argparse
import csv
import gzip
import math
from itertools import zip_longest
from pathlib import Path


def lines(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8", newline="")


def compare(expected, actual, *, rtol, atol, ignore=()):
    if not actual.is_file():
        return ("MISSING", 0, 0, "Replay did not produce this table")
    near = 0
    compared = 0
    with lines(expected) as a, lines(actual) as b:
        frozen, replayed = csv.reader(a, delimiter="\t"), csv.reader(b, delimiter="\t")
        header = next(frozen)
        if header != next(replayed):
            return ("FAIL", 0, 0, "Columns or column order differ")
        for i, (left, right) in enumerate(zip_longest(frozen, replayed), start=2):
            if left is None or right is None:
                return ("FAIL", compared, near, "Row counts differ")
            if len(left) != len(header) or len(right) != len(header):
                return ("FAIL", compared, near, f"Invalid column count at row {i}")
            for column, x, y in zip(header, left, right):
                if column in ignore or x == y:
                    continue
                try:
                    u, v = float(x), float(y)
                except ValueError:
                    return ("FAIL", compared, near, f"Text differs at row {i}, {column}: {x!r} vs {y!r}")
                if not (math.isfinite(u) and math.isfinite(v) and math.isclose(u, v, rel_tol=rtol, abs_tol=atol)):
                    return ("FAIL", compared, near, f"Numeric field differs at row {i}, {column}: {x!r} vs {y!r}")
                near += 1
            compared += 1
    remark = "All compared cells matched" if not ignore else f"Machine-specific provenance columns excluded: {', '.join(ignore)}"
    return ("PASS", compared, near, remark)


def write_report(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("table", "status", "rows_compared", "tolerance_matches", "detail"))
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frozen", type=Path)
    parser.add_argument("replayed", type=Path)
    parser.add_argument("output_report", type=Path)
    parser.add_argument("tables", nargs="*")
    args = parser.parse_args()
    names = args.tables or sorted(p.name for p in args.frozen.iterdir() if p.name.endswith((".tsv", ".tsv.gz")))
    if not names:
        parser.error("No frozen TSV tables selected")
    rows = [(name, *compare(args.frozen / name, args.replayed / name, rtol=1e-7, atol=1e-10)) for name in names]
    write_report(args.output_report, rows)
    for file, status, _, _, remark in rows:
        if status != "PASS":
            print(f"{status}\t{file}\t{remark}")
    if any(row[1] != "PASS" for row in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
