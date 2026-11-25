"""Verify the v1.2.0 archive, historical member parity and exact release fixtures."""
import argparse
from pathlib import Path
import re
import shutil
import stat
import sys
import tempfile
import zipfile

REPO = Path(__file__).resolve().parents[2]
VERSION = "ASD_TREG_CRISPR_processed_data_v1.2.0"
MANIFEST = "data/manifests/RELEASE_FILE_MANIFEST.tsv"
sys.path.insert(0, str(REPO / "tests"))
from verify_release_manifest import verify  # noqa: E402
from layout import digest, historical_members, rows, safe_path, verify_members


def extract(archive, dest):
    """Check the exact delivered file and reject unsafe ZIPs before writing."""
    sha_file = archive.with_name(archive.name + ".sha256")
    expected, listed_name = sha_file.read_text().split(maxsplit=1)
    assert listed_name.lstrip("*").strip() == archive.name
    assert re.fullmatch("[0-9a-f]{64}", expected) and digest(archive) == expected
    assert not dest.exists() or (dest.is_dir() and not any(dest.iterdir()))
    with zipfile.ZipFile(archive) as source:
        names = set()
        for item in source.infolist():
            rel = safe_path(item.filename)
            assert rel.parts[0] == VERSION and rel not in names
            assert not stat.S_ISLNK(item.external_attr >> 16), rel
            names.add(rel)
        assert source.testzip() is None
        dest.mkdir(parents=True, exist_ok=True)
        for item in source.infolist():
            target = dest / safe_path(item.filename)
            if item.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with source.open(item) as src, target.open("xb") as out:
                shutil.copyfileobj(src, out)
    return dest / VERSION


def validate(root, report=None):
    assert (root / "VERSION.txt").read_text().strip() == VERSION
    files, schemas = verify(root, report)
    manifest = root / MANIFEST
    entries = rows(manifest)
    inventory = {r["path"]: r for r in entries}
    assert len(inventory) == len(entries)
    observed = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    assert set(inventory) == observed - {MANIFEST, "SHA256SUMS.txt"}
    assert manifest.read_bytes() == (
        REPO / "data/manifests/DATA_RELEASE_FILE_MANIFEST.tsv"
    ).read_bytes()
    assert (root / "SHA256SUMS.txt").read_bytes() == (
        REPO / "data/manifests/DATA_RELEASE_SHA256SUMS.txt"
    ).read_bytes()
    for rel, row in inventory.items():
        target = root / safe_path(rel)
        assert target.stat().st_size == int(row["bytes"]), rel
        assert digest(target) == row["sha256"], rel

    history_count = historical_members(root)
    aliases = verify_members(root)
    mirrors = ("COMPONENT_SOURCE_MAP.tsv", "EXTERNAL_REPLAY_INPUTS.tsv")
    for name in mirrors:
        rel = "data/manifests/" + name
        assert (root / rel).read_bytes() == (REPO / rel).read_bytes(), rel

    fixture_count = 0
    for table in rows(REPO / "expected_results/V120_FROZEN_RESULTS_SHA256.tsv"):
        path = root / safe_path(table["path"])
        assert digest(path) == table["sha256"], table["path"]
        fixture_count += 1
    print(f"Release integrity PASS: {files} checksummed files, {schemas} TSV schemas, "
          f"{history_count} historical members, {fixture_count} scientific fixtures, {aliases} replay aliases.", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path)
    p.add_argument("--extract-to", type=Path)
    p.add_argument("--report", type=Path)
    p.add_argument("--legacy-view", type=Path)
    args = p.parse_args()
    with tempfile.TemporaryDirectory(prefix="asd-v120-verify-") as temporary:
        root = args.data.resolve()
        if root.is_file():
            root = extract(root, args.extract_to or Path(temporary))
        elif args.extract_to:
            p.error("--extract-to is only accepted with the ZIP")
        validate(root, args.report)
        if args.legacy_view:
            historical_members(root, args.legacy_view)


if __name__ == "__main__":
    main()
