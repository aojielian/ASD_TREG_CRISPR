"""Regenerate final vector figures, PDF-matched TIFFs, and supplemental pages."""
import argparse
import csv
import os
import re
import sys
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS.parent / "release"))
from layout import canonical_path, replay_view


def records(name):
    with (SCRIPTS / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def render(data, out, require_binary_match, submission_legends):
    from figure_common import sha
    from build_figures import build
    from verify_figure_numbers import audit
    from verify_figures import verify
    from write_source_map import write_map

    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("Choose a new or empty output directory.")
    inputs = records("inputs.tsv")
    for item in inputs:
        source = data / item["path"]
        if sha(source) != item["sha256"]:
            raise ValueError(f"Frozen input digest mismatch: {item['path']}")
    qa = out / "07_qa"
    qa.mkdir(parents=True, exist_ok=True)
    for item in inputs:
        if item["staging_path"]:
            destination = out / item["staging_path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            copyfile(data / item["path"], destination)
    with replay_view(data) as compatibility:
        build(out, compatibility)
        verify(out, compatibility / "figure_sources/frozen_r10", qa / "rendered_pages")
        audit(out, compatibility)
    if submission_legends:
        copy_submission_legends(data, out)
    write_map(data, out / "04_figures/FIGURE_SOURCE_MAP.tsv")

    parity = []
    for item in records("expected_figures.tsv"):
        relative = item["path"]
        expected = item["sha256"]
        actual = sha(out / relative)
        same = expected == actual
        parity.append((relative, item["generation"], expected, actual, "BYTE_IDENTICAL" if same else "BYTE_DIFFERENT"))
    write_tsv(qa / "FIGURE_REPLAY_PARITY.tsv", ("path", "generation", "frozen_sha256", "replay_sha256", "comparison"), parity)
    different = [row for row in parity if row[-1] != "BYTE_IDENTICAL"]
    if different and (require_binary_match or any(row[1] == "frozen_vector_copy" for row in different)):
        raise ValueError("Figure byte parity failed for " + ", ".join(row[0] for row in different))
    print(f"{len(parity) - len(different)}/{len(parity)} artifacts byte-identical; source and numeric checks complete.")


def copy_submission_legends(data, out):
    specs = (
        ("Figure_Legends.md", "04_figures", r"^## Figure (\d+)\.", 4),
        ("Supplementary_Figure_Legends.md", "05_supplementary", r"^## Supplementary Figure S(\d+)\.", 23),
    )
    for filename, output_directory, headings, count in specs:
        source = canonical_path(data, f"figure_sources/R11/{filename}")
        text = source.read_text(encoding="utf-8")
        found = [int(number) for number in re.findall(headings, text, flags=re.MULTILINE)]
        if found != list(range(1, count + 1)):
            raise ValueError(f"Incomplete published legend numbering in {filename}")
        copyfile(source, out / output_directory / filename)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="extracted processed-data archive root")
    parser.add_argument("out", type=Path, help="new or empty output directory")
    parser.add_argument("--strict-byte-parity", action="store_true", help="also require exact generated PDF and TIFF hashes")
    parser.add_argument("--submission-legends", action="store_true", help="also export the release's final editorial legends")
    args = parser.parse_args()
    data, out = args.data.resolve(), args.out.resolve()
    if out.is_relative_to(data) or out.is_relative_to(SCRIPTS.parents[1]):
        parser.error("Keep output outside the source and data directories")
    with TemporaryDirectory(prefix="asd_figure_cache_") as cache:
        os.environ["MPLCONFIGDIR"] = cache
        os.environ["XDG_CACHE_HOME"] = cache
        render(data, out, args.strict_byte_parity, args.submission_legends)


if __name__ == "__main__":
    main()
