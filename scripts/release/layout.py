"""Resolve canonical releases and construct temporary views for frozen programs."""
import argparse
from contextlib import contextmanager
import csv
from functools import lru_cache
import hashlib
from pathlib import Path, PurePosixPath
import shutil
import tempfile
import zipfile

COMPATIBILITY = Path("data/manifests/compatibility")
REPLAY_MAP = COMPATIBILITY / "REPLAY_PATHS.tsv"
LEGACY_MAP = COMPATIBILITY / "V110_MEMBERS.tsv"
METADATA_ZIP = COMPATIBILITY / "historical_metadata.zip"


def safe_path(text):
    name = str(text)
    path = PurePosixPath(name)
    if (not path.parts or path.is_absolute() or ".." in path.parts or
            "\\" in name or any(ord(c) < 32 for c in name)):
        raise ValueError(f"Unsafe relative path: {name}")
    return path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@lru_cache(maxsize=8)
def replay_members(data):
    path = data / REPLAY_MAP
    if not path.is_file():
        return {}
    entries = rows(path)
    index = {str(safe_path(item["legacy_path"])): item for item in entries}
    if len(index) != len(entries):
        raise ValueError("Duplicate replay alias")
    return index


def canonical_path(data, old_path):
    """Resolve a canonical or historical relative spelling within the named data root."""
    data = Path(data).resolve()
    relative = str(safe_path(old_path))
    entry = replay_members(data).get(relative)
    return data / safe_path(entry["canonical_path"]) if entry else data / relative


def verify_members(root):
    """Verify alias bytes, sizes and exact eligibility for historical reconstruction."""
    aliases = replay_members(root)
    checked = {}
    for entry in aliases.values():
        rel = safe_path(entry["canonical_path"])
        target = root / rel
        if rel not in checked:
            if not target.is_file() or target.is_symlink():
                raise ValueError(f"Missing canonical replay source: {rel}")
            checked[rel] = (target.stat().st_size, digest(target))
        if checked[rel] != (int(entry["bytes"]), entry["sha256"]):
            raise ValueError(f"Canonical replay source digest mismatch: {rel}")
    return len(aliases)


def materialize_replay(data, view):
    """Write exact alias copies in a new directory outside the canonical input root."""
    data, view = Path(data).resolve(), Path(view).resolve()
    if view == data or data in view.parents or view.exists():
        raise ValueError("Replay view must be a new directory outside the data root.")
    aliases = replay_members(data)
    if not aliases:
        raise ValueError("Canonical release replay map missing.")
    view.mkdir()
    copied = {}
    for legacy, entry in aliases.items():
        source = data / safe_path(entry["canonical_path"])
        if source not in copied:
            if not source.is_file() or source.is_symlink():
                raise ValueError(f"Missing canonical replay source: {entry['canonical_path']}")
            copied[source] = (source.stat().st_size, digest(source))
        if copied[source] != (int(entry["bytes"]), entry["sha256"]):
            raise ValueError(f"Canonical replay source digest mismatch: {entry['canonical_path']}")
        target = view / safe_path(legacy)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


@contextmanager
def replay_view(data):
    """Expose frozen spellings only for the lifetime of a calculation."""
    with tempfile.TemporaryDirectory(prefix="asd-science-inputs-") as scratch:
        view = Path(scratch) / "inputs"
        materialize_replay(Path(data), view)
        yield view


def historical_members(root, output=None):
    """Verify every historical member and optionally reconstruct its exact bytes."""
    entries = rows(root / LEGACY_MAP)
    by_name = {str(safe_path(r["historical_path"])): r for r in entries}
    if len(by_name) != len(entries):
        raise ValueError("Duplicate historical member.")
    if output is not None:
        output = Path(output).resolve()
        if output == root or root in output.parents or output.exists():
            raise ValueError("Legacy view must be a new directory outside the data root.")
    with zipfile.ZipFile(root / METADATA_ZIP) as metadata:
        info = metadata.infolist()
        included = {str(safe_path(i.filename)) for i in info}
        required = {r["metadata_member"] for r in entries if r["kind"] == "metadata_zip"}
        if len(included) != len(info) or included != required or metadata.testzip() is not None:
            raise ValueError("Historical metadata ZIP contents or CRC differ.")
        checksums = dict(reversed(line.split(maxsplit=1)) for line in
                         metadata.read("SHA256SUMS.txt").decode().splitlines())
        historical_inventory = list(csv.DictReader(
            metadata.read("data/manifests/RELEASE_FILE_MANIFEST.tsv").decode().splitlines(),
            delimiter="\t",
        ))
        if set(by_name) != set(checksums) | {"SHA256SUMS.txt"}:
            raise ValueError("Original published checksum membership differs.")
        for key, value in checksums.items():
            if by_name[key]["sha256"] != value:
                raise ValueError(f"Original published checksum differs: {key}")
        for old in historical_inventory:
            current = by_name[old["path"]]
            if old["sha256"] != current["sha256"] or old["bytes"] != current["bytes"]:
                raise ValueError(f"Original published manifest disagrees: {old['path']}")
        if output is not None:
            output.mkdir(parents=True)
        for path, entry in by_name.items():
            relative = safe_path(entry["canonical_path"])
            if entry["kind"] == "canonical" and not entry["metadata_member"]:
                source = root / relative
                if source.is_symlink():
                    raise ValueError(f"Unexpected symlink: {relative}")
                content = source.read_bytes()
            elif entry["kind"] == "metadata_zip" and relative == METADATA_ZIP:
                content = metadata.read(str(safe_path(entry["metadata_member"])))
            else:
                raise ValueError(f"Unknown historical storage mode: {path}")
            if (len(content), hashlib.sha256(content).hexdigest()) != (int(entry["bytes"]), entry["sha256"]):
                raise ValueError(f"Historical member bytes differ: {path}")
            if output is not None:
                target = output / path
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as stream:
                    stream.write(content)
    return len(entries)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    choice = cli.add_subparsers(dest="action", required=True)
    for action in ("replay", "resolve"):
        parser = choice.add_parser(action)
        parser.add_argument("data", type=Path)
        parser.add_argument("target")
    args = cli.parse_args()
    if args.action == "replay":
        materialize_replay(args.data, Path(args.target))
    else:
        print(canonical_path(args.data, args.target))


if __name__ == "__main__":
    main()
